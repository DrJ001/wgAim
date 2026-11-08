# =============================================================================
# demo_checkPanel.R
# Demonstration of checkPanel() -> filterPanel() -> primePanel() pipeline.
#
# Generates a synthetic 150-line, 300-marker panel with known data quality
# issues introduced deliberately so the three functions have something to find.
#
# Filter order in filterPanel() (8 steps):
#   1. Map consistency
#   2. Marker missingness
#   3. Line missingness
#   4. Line heterozygosity   (het.line)
#   5. Marker heterozygosity (het.marker)
#   6. Duplicate lines
#   7. Duplicate markers
#   8. MAF
# =============================================================================

devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim",
                   quiet = TRUE)

# =============================================================================
# 1.  Generate synthetic panel data
# =============================================================================

set.seed(2025)

n_lines   <- 150
n_chr     <- 3
n_mar_chr <- 100        # markers per chromosome
n_markers <- n_chr * n_mar_chr

line_ids <- paste0("Line_", sprintf("%03d", seq_len(n_lines)))
chrs     <- paste0("Chr", 1:n_chr)
mar_ids  <- paste0(rep(chrs, each = n_mar_chr), "_M",
                   rep(sprintf("%03d", seq_len(n_mar_chr)), times = n_chr))

# Allele frequencies from a beta distribution (realistic MAF spread)
p <- rbeta(n_markers, 2, 2)
p <- pmax(pmin(p, 0.95), 0.05)

# Force 10 monomorphic markers
mono_idx <- sample(n_markers, 10)
p[mono_idx] <- 0

# Force 10 additional low-MAF markers (MAF 0.01-0.04)
lowmaf_idx <- sample(setdiff(seq_len(n_markers), mono_idx), 10)
p[lowmaf_idx] <- runif(10, 0.01, 0.04)

# Generate 0/1/2 genotypes under HWE
geno <- matrix(NA_real_, nrow = n_lines, ncol = n_markers,
               dimnames = list(line_ids, mar_ids))
for (j in seq_len(n_markers)) {
    pj <- p[j]
    if (pj == 0) {
        geno[, j] <- 0
    } else {
        pr <- c((1 - pj)^2, 2 * pj * (1 - pj), pj^2)
        geno[, j] <- sample(0:2, n_lines, replace = TRUE, prob = pr)
    }
}

# ---- Introduce known issues ----

# (a) ~3% random missing data
miss_idx <- sample(length(geno), round(0.03 * length(geno)))
geno[miss_idx] <- NA

# (b) 3 high-missingness lines (>35% missing)
for (i in 1:3)
    geno[i, sample(n_markers, round(0.38 * n_markers))] <- NA

# (c) 4 high-missingness markers (first 4 of Chr2, >30% missing)
for (m in paste0("Chr2_M", sprintf("%03d", 1:4)))
    geno[sample(n_lines, round(0.32 * n_lines)), m] <- NA

# (d) 2 duplicate lines (Line_010 and Line_011 are copies of Line_005)
geno["Line_010", ] <- geno["Line_005", ]
geno["Line_011", ] <- geno["Line_005", ]

# (e) 2 duplicate markers (Chr3_M050 and Chr3_M051 are copies of Chr3_M001)
geno[, "Chr3_M050"] <- geno[, "Chr3_M001"]
geno[, "Chr3_M051"] <- geno[, "Chr3_M001"]

# (f) 3 high-heterozygosity lines (Line_020 to Line_022, ~55% het calls)
#     Note: HWE genotypes already carry ~33% natural het (2pq at p=0.5),
#     so seeded values need to be well above that to stand out.
for (i in 20:22)
    geno[i, sample(n_markers, round(0.55 * n_markers))] <- 1L

# (g) 5 high-heterozygosity markers (~55% het calls across all lines)
het_mar_idx <- paste0("Chr1_M", sprintf("%03d", 50:54))
for (m in het_mar_idx)
    geno[sample(n_lines, round(0.55 * n_lines)), m] <- 1L

# ---- Build map (leave 2 markers out to test consistency check) ----
map_df <- data.frame(
    marker = mar_ids,
    chr    = rep(chrs, each = n_mar_chr),
    pos    = rep(seq(0, 200, length.out = n_mar_chr), times = n_chr),
    stringsAsFactors = FALSE
)
map_df <- map_df[!map_df$marker %in% c("Chr1_M098", "Chr1_M099"), ]

cat("\n")
cat("=============================================================\n")
cat("Synthetic panel:  150 lines,  300 markers,  3 chromosomes\n")
cat("Known issues introduced:\n")
cat("  - 3 lines with >35% missing data  (Line_001-003)\n")
cat("  - 4 markers with >30% missing data (Chr2_M001-M004)\n")
cat("  - 3 high-het lines (~55% het)     (Line_020-022)\n")
cat("  - 5 high-het markers (~55% het)   (Chr1_M050-M054)\n")
cat("  - 2 duplicate lines  (Line_010, Line_011 = Line_005)\n")
cat("  - 2 duplicate markers (Chr3_M050, Chr3_M051 = Chr3_M001)\n")
cat("  - 10 monomorphic markers\n")
cat("  - 10 low-MAF markers  (MAF 0.01-0.04)\n")
cat("  - 2 markers in geno not in map (Chr1_M098, Chr1_M099)\n")
cat("=============================================================\n\n")

# =============================================================================
# 2.  checkPanel() — diagnose only
# =============================================================================

cat("──────────────────────────────────────────────────────────────\n")
cat("STEP 1: checkPanel()\n")
cat("──────────────────────────────────────────────────────────────\n\n")

chk <- checkPanel(geno, map_df,
                  map.id   = "marker",
                  map.chr  = "chr",
                  map.pos  = "pos",
                  encoding = "012")

# print() uses a single het.thresh applied to both per-line and per-marker
# het summaries.  Set it to 0.45 to flag what the filterPanel het.line
# threshold will remove — the per-marker filter uses 0.60 (see note below),
# so the flagged marker count here will be conservative.
# NOTE: HWE-generated genotypes naturally carry ~33% heterozygosity
# (2pq ≈ 0.5 at balanced allele frequencies), so the threshold must
# sit well above that to target only the deliberately seeded bad cases.
print(chk,
      miss.line.thresh   = 0.20,
      miss.marker.thresh = 0.20,
      het.thresh         = 0.45)

# =============================================================================
# 3.  filterPanel() — apply thresholds informed by checkPanel output
#
# Filter order (8 steps):
#   1. Map consistency        -- unplaceable markers dropped
#   2. Marker missingness     -- miss.marker = 0.20
#   3. Line missingness       -- miss.line   = 0.20
#   4. Line heterozygosity    -- het.line    = 0.45
#   5. Marker heterozygosity  -- het.marker  = 0.60
#   6. Duplicate lines        -- dup.lines   = TRUE
#   7. Duplicate markers      -- dup.markers = TRUE
#   8. MAF                    -- maf         = 0.05
# NOTE: het.marker is set higher than het.line because HWE-generated markers
# can naturally reach ~50% het at balanced frequencies, so the per-marker
# threshold needs a wider margin to avoid sweeping out real data.
# =============================================================================

cat("──────────────────────────────────────────────────────────────\n")
cat("STEP 2: filterPanel()\n")
cat("──────────────────────────────────────────────────────────────\n\n")

clean <- filterPanel(chk,
                     miss.marker = 0.20,
                     miss.line   = 0.20,
                     het.line    = 0.45,
                     het.marker  = 0.60,
                     dup.lines   = TRUE,
                     dup.markers = TRUE,
                     maf         = 0.05)

print(clean)

# =============================================================================
# 4.  primePanel() — prepare for analysis
# =============================================================================

cat("──────────────────────────────────────────────────────────────\n")
cat("STEP 3: primePanel()\n")
cat("──────────────────────────────────────────────────────────────\n\n")

panel <- primePanel(clean, impute = "knn", knn.k = 5)

cat("\nwgPanel structure:\n")
cat(sprintf("  Class      : %s\n", paste(class(panel), collapse = ", ")))
cat(sprintf("  Lines      : %d\n", nrow(panel$pheno)))
cat(sprintf("  Chromosomes: %d\n", length(panel$geno)))
for (ch in names(panel$geno))
    cat(sprintf("    %-8s : %d markers  [%.1f - %.1f cM]\n",
                ch,
                ncol(panel$geno[[ch]]$imputed.data),
                min(panel$geno[[ch]]$map),
                max(panel$geno[[ch]]$map)))
