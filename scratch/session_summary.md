# wgAim — Session Reference

*Last updated: Session 38 — 2026-11-06*

---

## Session 39 (2026-11-07/08) — rising LRT + AOM stopping-rule research

Colleague's 3-trial oat run: LRT statistic kept *rising* as markers were selected. Diagnosed as the LRT testing "does polygenic background exist?" (always yes with strong background), not "is there another QTL?". Prototyped an **AOM-based stopping rule** (score test of the best marker against the polygenic background, with exact simulated family-wise null) as standalone scratch code — **package code unchanged**. Stage 0 (engine-path verification) and Stage 1 (null calibration, 2,600 replicates: 0.2–4.2% at nominal 5%, conservative) done. **Full details, files, results and next steps: [selection_summary.md](../research/aom_stopping/selection_summary.md).**

**Now tracked in git** on branch `feature/aom-stopping-rule` (from `dev`) under `research/aom_stopping/` (scripts, maths qmd/html, `selection_summary.md`, `results/`); `^research$` added to `.Rbuildignore`. That folder is the canonical copy — the originals in `scratch/` are superseded. Scripts run from the repo root with repo-relative paths.

**`scratch/` is now tracked in git too** (still excluded from the package build via `^scratch$` in `.Rbuildignore`). `.gitignore` keeps out: `scratch/rlib/` (local test install), `*.zip` / `*.pdf` (ASReml installer, ASReml manual, Springer paper — licensed/copyrighted), `scratch/Beata/` and all `*.RData`/`*.Rdata`/`*.rda` (colleague data / workspace images), editor autosaves `*#`. **The repo is public** — check new scratch files for data or licensed material before committing.

---

## ⚠ Current on-disk state (end of Session 38)

- **Trigger:** a colleague's 3-trial `gwasAim()` with `corgh(Trial):id` (or `us`) failed after "Predict step for outlier statistics" with `Error in array(x, ...): 'data' must be of a vector type, was 'NULL'`. Their data is unavailable, so everything was reproduced on simulated 3/4/5-trial data built on the shipped GWAS panel genotypes (`scratch/repro_corgh_3trial.R [struct] [str|NULL] [n.trials]`; 191 genotyped + 9 ungenotyped lines so `.fixLines()`/`Gsave` is exercised; planted Chr2 +0.20 MAIN, Chr5 −0.18 MAIN, Chr8 crossover).
- **Engine untouched for `fa`/`diag` + `vm`** apart from `Perc.Var` (below), so UV results and the FA-based vignette fits are unchanged. Everything else in the MV engine for corh/corgh/us and for the `mbf` path was broken in one way or another and is now fixed — see the Session 38 architecture section.
- **Committed and pushed:** `707ab65` on `dev` ("Fix multivariate engine for corgh/us/corh structures and MV mbf path", 19 files, +664 −166), merged into `main` via `--no-ff` as `47ddf25`. Before committing, local `dev` was created at `main` (fast-forward of `origin/dev`), which **resolved the S37 `deb1f8e` divergence** — `dev` now contains all of `main`. Both pushes were plain fast-forwards. Commit identity passed per command (`git -c user.name=DrJ001 -c user.email=7767208+DrJ001@users.noreply.github.com`) — **no git identity is configured on this machine**; user chose per-command identity, config untouched.
- **Local tooling caveats (this machine):** the R with ASReml is `C:\Users\julia\AppData\Local\Programs\R\R-4.5.3` (ASReml-R 4.2.0.481); it has **no `devtools`, `pkgload`, `testthat` or `roxygen2`**. Package was tested by `R CMD INSTALL` into `scratch/rlib` (git-ignored) and a minimal testthat shim `scratch/run_tests_shim.R` (492 pass / 4 fail — all 4 are shim artefacts: 3 `sink.number()` checks broken by the shim's `capture.output()`, 1 unsupported `label=` arg). `man/gwasAim.Rd`, `man/qtlAim.Rd`, `man/gpAim.Rd` were **hand-edited** to mirror the roxygen changes (pass `tools::checkRd`). `gh` CLI not installed — CI status after the push not checked. The shell is PowerShell (no heredocs; commit messages via `-F file`).
- **Stale artefacts:** MV vignettes (`qtlAim_mv_pipeline`, `gwasAim_mv_pipeline`) use `fa(2)` models, so their `Perc.Var` numbers (summary tables, effects plots) change with the Perc.Var fix — `inst/doc/*.html` and the pkgdown site still show the old values. `scratch/demo_gpAim_mv.R` results also change (gpAim now applies the MV structure upgrade).

---

## Current on-disk state (end of Session 37)

- **GWAS sim tuned for reliable MV detection with 200-line panel** (previous 300-line sim ran into ASReml stability issues with the older reliable ASReml build the user reverted to). Iterative calibration across ~10 rounds converged on parameters that give clean 5/5 QTL detection with correct MAIN/INT classification and zero spurious hits.
- **Panel: 200 lines × ~1691 markers** across 10 chr. `data-raw/simulate_GWASsim.R` now uses `addSnpChip(nSnpPerChr = 175, minSnpFreq = 0.15)` for active MAF-based selection at the AlphaSimR founder level, plus a belt-and-braces F1 MAF > 0.15 filter after `randCross` to catch any drift. Marker renumbering moved to AFTER the filter so marker names remain contiguous. Panel MAF distribution: min 0.152, median 0.320, max 0.500.
- **Response model**: `y = μ_site + Site:Rep + Σ β_jk X_ik + U_ij + L_ij + ε` where U ~ MVN(0, Ga_K ⊗ K) is TRACE K-weighted polygenic (target diag 0.005 — deliberately small so LRT stops once QTL are absorbed) and L ~ MVN(0, Ga_L ⊗ I) is moderate IID-line polygenic (target diag 0.020 — keeps fa(Site,k):id well-conditioned without swamping vm signal). Residual variance 0.142 (SD ~0.377 t/ha, ~8% of ~4.7 t/ha mean). **Ga_L seed changed 1729 → 1730** to fix an Env3 outlier variance (was 2.6× target, causing MAIN QTL to be misclassified as G×E because Env3 estimates wandered).
- **Planted QTL (5, effects in t/ha)** — all with MAF > 0.25 in current shipped data:

  | QTL | Chr | cM | Type | β | Marker (current) | MAF |
  |---|---|---:|---|---:|---|---:|
  | 1 | Chr2 | 45 | MAIN | +0.145 | `SNP_C02_0044` | 0.27 |
  | 2 | Chr5 | 90 | MAIN | −0.145 | `SNP_C05_0094` | 0.43 |
  | 3 | Chr7 | 30 | MAIN-weak | +0.125 | `SNP_C07_0042` | 0.48 |
  | 4 | Chr3 | 110 | GxE-fade | +0.230 → 0 | `SNP_C03_0128` | 0.48 |
  | 5 | Chr9 | 60 | GxE-cross | ±0.180 | `SNP_C09_0068` | 0.48 |

  Marker names shift each regen because AlphaSimR's coalescent picks different segregating sites; positions (cM) are stable.
- **`find_qtl_marker()` MAF floor raised 0.15 → 0.25** with improved highest-MAF fallback (instead of "closest marker regardless of MAF"). Keeps per-QTL variance stable across regenerations.
- **Yield scale: tonnes per hectare (t/ha)**, site means 4.5–4.9, overall mean ~4.7.
- `data/genoGWASraw.rda`, `data/mapGWASsim.rda`, `data/phenoGWASsim.rda` regenerated. `R/genoGWASsim.R` roxygen fully updated (panel size ~1700, MAF > 0.15 mandate, small polygenic components explained, new QTL effect table + shipped marker names).
- **`scratch/smoke_gwas_mv_new_sim.R` header refreshed** to reflect the current sim parameters (200 lines, both polygenic components, current QTL effects).
- **`scratch/smoke_gwas_new_sim.R` (UV smoke test on Env1)** header refreshed with S37 sim details + planted truth table added at end.
- **User confirmed clean detection with `str = NULL`** on the final calibrated sim (previous S36 recommendation of `str = "corh"` no longer necessary — the small Ga_K and healthy MAF distribution avoid the boundary risk).
- **UV GWAS vignette (`vignettes/gwasAim_pipeline.qmd`) regenerated** against the S37 sim: Overview rewritten (200 lines, ~1700 markers, AlphaSimR generation, MAF > 0.15 mandate, t/ha scale, two-component polygenic); filterPanel narrative counts updated (~191 lines, ~1450 markers post-filter, MAF filter now a no-op); gwas_class narrative counts updated; planted-QTL description refreshed with current effects (+0.145/−0.145/+0.125/+0.230→0/±0.180 t/ha); `"Env01"` typo fixed to `"Env1"`; cache cleared and re-rendered; `inst/doc/gwasAim_pipeline.html` refreshed. User confirmed output "spot on".
- **MV GWAS vignette (`vignettes/gwasAim_mv_pipeline.qmd`) regenerated** against the S37 sim with the same treatment: Overview rewritten to match UV (AlphaSimR, MAF, polygenic, t/ha); QTL description updated to new effects; panel_prep narrative counts refreshed; Results narrative updated (Chr7 +0.125, Chr3 fade +0.230→0, Chr9 crossover ±0.180); Summary code-block `TypeI = 0.01 → 0.05` fixed (was inconsistent with the actual `gwasAim()` call in the vignette body); caches cleared and re-rendered; `inst/doc/gwasAim_mv_pipeline.html` refreshed (all 41 chunks executed cleanly, ~4.8MB HTML).
- **README + pkgdown site config updated**: README intro paragraph now says GWAS vignettes are available (only GP+selection remains "in preparation"); Vignettes table has two new rows linking to `articles/gwasAim_pipeline.html` and `articles/gwasAim_mv_pipeline.html`. `_pkgdown.yml` adds `gwasAim_mv_pipeline` to the GWAS articles list (the UV was added earlier). `.github/workflows/pkgdown.yaml` seeds both HTMLs into `docs/articles/` before `build_articles(lazy=TRUE)`. Datasets `genoGWASraw` / `mapGWASsim` / `phenoGWASsim` added to the pkgdown reference datasets section.
- **Package fixes rolled into the same commit** (previously uncommitted from S34–S35 work): `aimTrace.gwasAim()` missing `ncol` arg added + threaded to `.plot_stability()`; `.build_mv_contrast_df` factored `res_lab` to correctly propagate `"Gsave"` when `.fixLines()` triggers; `gwasAim()` exclusion-window message made unit-agnostic; small roxygen typography and Rcpp backtick fixes.
- **Deployed to GitHub end of session**: two commits on `dev` — `6e6aaac` (Rebuild GWAS pipeline: 24 files, +10331 -19) and `c8bd251` (README/pkgdown promotion: 3 files, +9 -4). Merged `dev → main` via `--no-ff` as merge commit `7dde086` "Merge dev: rebuild GWAS pipeline and promote both GWAS vignettes" following the established main-branch pattern. Both branches pushed to `origin`.
- **Divergence flag** (✅ resolved Session 38 — `dev` fast-forwarded to `main`) (surfaced during merge, not addressed): `main` has commit `deb1f8e` "Audit and fix all @return roxygen docs; fix print.qtlAim/gwasAim bug" that never landed on `dev`. It was made directly on main. Not urgent since it's on main and shipped, but `dev` is technically not a strict superset of main history. Worth remembering for a future housekeeping cherry-pick if strict lineage is wanted.

---

## Package Overview

**`wgAim`** — Whole Genome Analyses via Integrated Modelling. Evolves the legacy `wgaim` QTL package into a unified toolkit. Built on ASReml-R. Four analytical functions across three pillars:

| Pillar | Functions | Input object |
|---|---|---|
| QTL mapping | `qtlAim()` | `wgCross` from `primeCross()` |
| GWAS | `gwasAim()` | `wgPanel` from `primePanel()` |
| Genomic prediction & selection | `gpAim()` · `selIndex()` | `wgCross` or `wgPanel` |

**Current branch:** `main` checked out; `dev` = `707ab65`, `main` = `47ddf25` (both pushed). **Test suite:** 997 passing at S37 + 16 new tests in S38 (`test-extractVmGa.R`, `test-engine-helpers.R`) — full `devtools::test()` not yet re-run on a machine with testthat. **Version:** 0.9.0.

---

## R/ File Structure

### Shared Engine (internal, dot-prefixed)
`engine_validate.R` · `engine_genodata.R` · `engine_model.R` · `engine_select.R` · `engine_effects.R` · `engine_results.R` · `engine_lrt.R`

### Data Preparation
`primeCross.R` — `primeCross()`, `fixMap()`, `imputeGen()`, `lambdaf()`
`checkPanel.R` — `checkPanel()`, `print.checkPanel()`
`filterPanel.R` — `filterPanel()`, `.filterPanel_core()`, `print.filteredPanel()`
`primePanel.R` — `primePanel()`, `.primePanel_core()`, `.knn_impute_chr()`

### Analysis
`qtlAim.R` · `gwasAim.R` · `gpAim.R` · `selIndex.R`

### Display — qtlAim
`print.qtlAim.R` · `summary.qtlAim.R` · `plot.qtlAim.R`

### Display — gwasAim
`print.gwasAim.R` · `summary.gwasAim.R` · `plot.gwasAim.R`

### Display — gpAim
`print.gpAim.R` · `summary.gpAim.R` · `plot.gpAim.R`

### Display — selIndex
`print.selIndex.R` · `plot.selIndex.R` (summary in print.selIndex.R)

### linkMap
`linkMap.cross.R` · `linkMap.qtlAim.R` · `linkMap.gwasAim.R` · `linkMap.default.R`

### Other
`fineMap.R` · `print.fineMap.R` · `plot.fineMap.R` · `getQTL.R` · `aimTable.R` · `print.aimTable.R` · `aimTrace.qtlAim.R` · `aimTrace.gwasAim.R` · `theme_scatter.R` · `wgAim-package.R` · `zzz.R`

---

## Key Architecture Decisions

### Genotype encoding
- Biparental (BC/DH/RI): B = −1, A = +1; F2 adds AB = 0
- `primeCross(type = "interval")` → `$interval.data` (lambdaf midpoints) + `$inferred.map`
- `primeCross(type = "marker")` → `$imputed.data` only
- GWAS panel: 0/1/2 → −1/0/+1 internally in `primePanel()`; `$imputed.data`

### QTL internal naming
- State/column keys: `"Chr.<chr>.<idx>"` (e.g. `"Chr.3B.30"`)
- Formula/phenoData columns: `"X.<chr>.<idx>"` (e.g. `"X.3B.30"`)

### Engine paths (vm vs mbf)
- **vm**: markers > lines → `vm(gterm, G)` with genomic relationship matrix
- **mbf**: lines ≥ markers → `mbf('ints')` marker-by-file
- Both paths give identical QTL selection; `force = TRUE` overrides to mbf

### Total genetic value (contrast plot)
Three components summed per line:
1. `vm(gterm, covObj)` BLUP — additive genomic from remaining markers
2. `gterm`/`Gsave` BLUP — residual line-level genetic effect
3. `Σ β_k × x_ik` — detected QTL fixed contributions

Components 1+2 from single `predict(object, classify=gterm, only=c(vm_label, res_label))` call.

### Significance testing
Single LRT per iteration tests the additive variance of the composite genome-wide term — not per-marker. `TypeI = 0.05` is already family-wise. No Bonferroni needed.

### Multivariate (Trait) engine
- `Trait = NULL` → univariate (default, unchanged behaviour)
- `Trait = "Site"` → multi-environment G×E analysis
- Variance structure on the genomic term controlled by `str` argument (`"corh"`, `"fa1"`, etc.)
- LRT uses `pchisq.mixture(stat, ntrait)` after downgrading vm/mbf prefix to `diag(Trait):` for both models
- `waldTest` pruning via internal `.waldTest()` in `engine_results.R` determines MAIN vs INT QTL

---

## Function Signatures (current)

```r
qtlAim(baseModel, genObj, merge.by, fix.lines=TRUE, gen.type=NULL,
       method="fixed", selection="interval", force=FALSE,
       exclusion.window=20, breakout=-1, TypeI=0.05, trace=TRUE,
       verboseLev=0, Trait=NULL, str=NULL, ...)

gwasAim(baseModel, genObj, merge.by, fix.lines=TRUE, force=FALSE,
        exclusion.window=20, breakout=-1, TypeI=0.05, trace=TRUE,
        verboseLev=0, Trait=NULL, str=NULL, ...)
        # method="fixed", selection="interval" hardcoded

gpAim(baseModel, genObj, merge.by, fix.lines=TRUE, gen.type="marker",
      force=FALSE, trace=TRUE, Trait=NULL, str=NULL, ...)

selIndex(object, weights=NULL, desired=NULL,
         type=c("weighted","smith-hazel","desired-gains"),
         standardise=TRUE, prop.select=0.10, threshold=NULL, selected=NULL, ...)

primeCross(object, type="interval", impute="MartinezCurnow",
           consensus.mark=TRUE, id="id", subset=NULL, infer="mid")

checkPanel(geno, map, id="id", map.id="marker", map.chr="chr",
           map.pos="pos", encoding="012")
# Returns "checkPanel" object; pure diagnostic, no filtering

filterPanel(geno, map, id="id", map.id="marker", map.chr="chr",
            map.pos="pos", encoding="012",
            miss.marker=0.20, miss.line=0.20,
            het.line=NULL, het.marker=NULL,
            dup.lines=TRUE, dup.markers=FALSE, maf=0.05)
# geno may be a "checkPanel" object (stored fields extracted automatically)
# Returns "filteredPanel" object with $geno, $map, $removed, $thresholds
# Filter order: (1) map consistency, (2) marker miss, (3) line miss,
#   (4) het.line, (5) het.marker, (6) dup lines, (7) dup markers, (8) MAF

primePanel(geno, map, id="id", map.id="marker", map.chr="chr",
           map.pos="pos", encoding="012", impute="none", knn.k=5L)
# geno may be a "filteredPanel" object (stored fields extracted automatically)
# maf REMOVED — filtering is filterPanel()'s responsibility
# impute = "none" | "mean" | "knn"  (knn is chromosome-wise vectorised KNN)

# Plot methods
plot(x, genObj,
     type = "outlier"|"blups"|"chr"|"effects"|"contrast"|"heatmap",
     data = NULL,    # required for type="contrast"
     iter, chr, chr.lines, sig.col, pt.col, pt.cex, ...)
     # ncol, qtl passed via ... for type="contrast"

fineMap(object, genObj, qtl=NULL, window=50, step=2,
        exclusion.window=20, ...)

linkMap(object, ..., genObj,
        chr, chr.dist, nrow=1, row.chr=NULL, marker.names="none")

# selIndex plot
plot(x, type=c("index","heatmap","weights"),
     prop.select=0.10, threshold=NULL,
     pt.col=c("steelblue","firebrick"), pt.cex=0.9, ...)
```

---

## plot.qtlAim / plot.gwasAim — type options

| `type` | Description | Extra args |
|---|---|---|
| `"outlier"` (default qtlAim) | Genome-wide outlier stats, faceted by iteration | `iter`, `chr`, `chr.lines` |
| `"blups"` | Scaled BLUPs across genome; MV = one line per trial | same |
| `"chr"` | Chromosome-level bar chart (`selection="chromosome"` only) | `iter`, `chr` |
| `"manhattan"` (default gwasAim) | Manhattan plot with chr backgrounds | `iter`, `chr` |
| `"effects"` | Lollipop: effect ± SE, % var. UV = one row/QTL; MV = one row/trial | — |
| `"contrast"` | Allele scatter: total genetic value by allele class | `data` (phenoData) |
| `"heatmap"` | Genome × iteration heatmap: fill = outlier stat | `iter`, `chr` |

**MV dispatch:** activated when `!is.null(x$QTL$Trait)`. Uses `.build_mv_blups_df()`, `.build_mv_effects_df()`, `.plot_mv_effects()`, `.build_mv_contrast_df()`.

**`"effects"` colour:** firebrick = A (+1) favoured; steelblue = B (−1) favoured.

**`"contrast"` allele classes:**
- Biparental: continuous score from `interval.data`/`imputed.data`, x-axis "B (−1)" / "A (+1)"
- GWAS: discrete 0/1/2 mapped back from ±1 column

---

## Shared Internal Helpers — plot.qtlAim.R

| Helper | Used by |
|---|---|
| `.build_cumpos(genObj, gen.type, chr)` | outlier/blups |
| `.build_stat_df(object, slot, iter, chr, cp)` | outlier/blups |
| `.add_sig_labels(gp, object, ...)` | outlier/blups (UV) |
| `.add_mv_blups_flags(gp, object, ...)` | blups (MV) |
| `.plot_chr_bars(object, iter, chr, sig.col)` | chr bar chart |
| `.match_perc_var(summ, object, genObj, effects, gen.type)` | effects + contrast |
| `.build_effects_df(object, genObj)` | effects lollipop (UV) |
| `.plot_effects(edf)` | effects lollipop (UV) |
| `.build_mv_effects_df(object, genObj)` | effects lollipop (MV) |
| `.plot_mv_effects(edf, object)` | effects lollipop (MV) |
| `.build_mv_blups_df(object, iter, chr, cp)` | blups (MV) |
| `.build_contrast_df(object, genObj, data)` | contrast (UV) |
| `.build_mv_contrast_df(object, genObj, data, qtl)` | contrast (MV) |
| `.plot_contrast(cdf_list, ncol)` | contrast (UV + MV) |
| `.build_heatmap_df(object, iter, chr, cp)` | heatmap |
| `.plot_heatmap(object, genObj, iter, chr, sig.col)` | heatmap |

---

## aimTrace

```r
aimTrace(object, iter=seq_along(object$QTL$diag$coef.list),
         lik.out=TRUE, plot=FALSE, sig.col="firebrick", ...)
```

| `plot` value | Returns |
|---|---|
| `FALSE` (default) | `NULL` invisibly |
| `"lrt"` | ggplot — LRT statistic, threshold line, QTL labels |
| `"stability"` | ggplot — QTL effect ± 1 SE; MV = one line per trial |
| `"both"` | `list(lrt=..., stability=...)` |

**MV support:** `.plot_lrt()` uses `qchisq.mixture(1-TypeI, ntrait)` threshold. `.build_stability_df()` reconstructs per-trial effects from main + interaction coefficients.

---

## summary.qtlAim / summary.gwasAim

- Chromosome sort uses leading-integer extraction for alphanumeric names (e.g. "3B", "7A2")
- **MV path:** extracts coefficients from final pruned model (`object$coefficients$fixed`); leading `Trait` column when `$QTL$Trait` non-NULL; per-unique-QTL `perc.var` broadcast to all rows
- `perc.var` denominator: `sum(var.est.qtl) + coef.mark * var.mark.scalar + var.res.scalar` (means over trials for MV)
- **MV (Session 38):** `var.mark` / `var.res` are `diag()` of the fitted genomic / residual-genetic covariance matrices via `.mvPercVarComponents()` (structure-agnostic; corgh/us/fa(1) now agree). Previously the mean of *all* vparameters of each term, incl. correlations/covariances/loadings. Falls back to the old calculation if a structure cannot be identified.

---

## fineMap

- `qtlAim` path: fresh inferred grid anchored at `qtl_pos`; Haldane weights
- `gwasAim` path: scans actual panel markers within ±`window` cM; `step` ignored
- `attr(df, "qtl_pos")` and `attr(df, "clipped")` stored per-QTL
- `plot.fineMap()`: single `geom_text_repel()` for peak + original QTL labels
- **MV support:** fully compatible with multivariate `qtlAim`/`gwasAim` objects
  - Main-effect QTL (`is.interaction = FALSE`): UV z-ratio path unchanged; swaps `X.chr.idx` for candidate
  - Interaction QTL (`is.interaction = TRUE`): swaps `Trait:X.chr.idx` for `Trait:candidate`; joint Wald test via `.waldTest(cc=list(list(coef=whr_nz, type="zero")))` across non-aliased coefficients; LOD = χ²/(2 log 10)

---

## selIndex — index types

| `type` | Formula | Inputs |
|---|---|---|
| `"weighted"` | `I = b'g` | `weights` (named, one per env) |
| `"smith-hazel"` | `b = P⁻¹ Ga w` | `weights` |
| `"desired-gains"` | `b ∝ Ga⁻¹ d`, unit L2 norm | `desired` |

- `standardise=TRUE`: GEBVs scaled to unit variance per env before combining
- `selected=NULL`: top `prop.select` by index; supply character vector for external selection
- `$gain`: `ΔG_j = mean(GEBV_j[selected]) − mean(GEBV_j[all])`, original scale

---

## Multivariate $QTL slots

| Slot | Content |
|---|---|
| `$QTL$Trait` | Character: the Trait factor column name |
| `$QTL$trait.levels` | Character vector: factor levels |
| `$QTL$wald.test` | Data frame: per-QTL Wald zero-equality results |
| `$QTL$final.terms` | Character: post-pruning fixed term names |
| `$QTL$is.interaction` | Logical: TRUE = G×E QTL, FALSE = main effect |

---

## gpAim — GEBV extraction

- **vm path (markers > lines):** `predict(gpModel, classify=merge.by, only=vm_only, sed=TRUE)` → Mrode accuracy `sqrt(1 - PEV/Vg)`; Cullis H² from SED upper triangle (`.cullis_H2()`)
- **mbf path (lines ≥ markers):** GEBVs = `M %*% q.hat`; SE approximate; gen.H2 = NA
- **MV vm path:** `classify = paste(Trait, merge.by, sep=":")`, `only=vm_only`; per-trial Cullis H²; Ga extracted via `.extract_Ga_gp()`
- `str` argument controls variance structure: `NULL` = mirror base model, `"corh"`, `"fa1"`, `"fa2"`, `"diag"`, etc.

---

## engine_results.R — .waldTest() (internal)

Replaces the former exported `waldTest()`. Zero-equality Wald tests only (Contrasts path removed). Returns a plain `data.frame` (not `list(Zero=...)`). Called by `.packResults()` only.

```r
.waldTest(object, cc)
# cc: list of list(coef=integer_indices, group="optional_label")
```

---

## .rhs_terms() — formula parsing helper (engine_model.R)

```r
.rhs_terms <- function(formula.obj) {
    raw   <- paste(deparse(formula.obj[[2L]], width.cutoff=500L), collapse=" ")
    terms <- trimws(unlist(strsplit(raw, "\\s*\\+\\s*")))
    terms[nchar(terms) > 0L]
}
```

Handles leading `+` and long-formula wrapping. Used in `.fixLines()`, `.buildGenomeModel()`, `.rebuildCovObj()`.

---

## CodeCov & CI

- Workflow: `.github/workflows/test-coverage.yaml` (codecov-action@v5), triggers on push to `dev` and `master`
- **Overall coverage (last upload):** ~75% (target: 85%)
- `R CMD check`: 0 errors, 0 warnings, 2 notes (package dependencies + future timestamps — both expected)
- `asreml::` prefix removed from all R/ files; `utils::globalVariables()` in `zzz.R` suppresses runtime-resolved symbols

---

## Test Suite — `tests/testthat/`

| File | Description |
|---|---|
| `helper-fixtures.R` | Shared synthetic fixtures (no ASReml): `make_bc_cross`, `make_wgCross_interval`, `make_wgCross_marker`, `make_wgPanel`, `make_mock_asreml`, `make_mock_qtlAim`, `make_mock_gwasAim`, `make_mock_gpAim`, `make_mock_fineMap` |
| `test-engine.R` | `.validateModel`, `.buildGenoData`, `.lrtTest`, `.packResults`, `.mergeEffect` |
| `test-engine-asreml.R` | Mock-based: `.qtlSelect`, `.buildGenomeModel`, `.rebuildCovObj`, `.vModify`, `.envFix`, `.addEffect`, full loop body, `fineMap` scan |
| `test-primeCross.R` | `primeCross()` (bc/dh/f2, interval/marker, impute, subset, fixMap) + `primePanel()` |
| `test-getQTL-aimTable.R` | `getQTL()`, `aimTable()`, `print.aimTable()` |
| `test-print-summary.R` | `print.*`, `summary.*` for qtlAim/gwasAim/gpAim |
| `test-print-gpAim.R` | UV + MV `print.gpAim` (gen.H2 present/absent, path vm/mbf, Gcor block) |
| `test-plot-qtlAim.R` | All 6 `type=` branches of `plot.qtlAim` + guard clauses |
| `test-plot-qtlAim-extra.R` | MV blups/effects helpers, `plot.gwasAim` MV blups/effects |
| `test-plot-gwasAim.R` | All 6 `type=` branches of `plot.gwasAim` |
| `test-plot-gpAim.R` | All 4 `type=` branches of `plot.gpAim`, `.gp_threshold` |
| `test-plot-fineMap.R` | `plot.fineMap`, `print.fineMap` |
| `test-linkMap.R` | `linkMap.cross`, `linkMap.qtlAim`, `linkMap.gwasAim`, `linkMap.default`, `fineMap` guards |
| `test-qtlAim.R` | `qtlAim()` dispatch, guard clauses, gen.type defaulting, trace file |
| `test-gwasAim.R` | `gwasAim()` dispatch, guard clauses, n.markers capture |
| `test-gwasAim-extra.R` | Trait guards, str guards, hard-coded args, genObj class guards |
| `test-gpAim.R` | `gpAim()` dispatch, guard clauses, vm vs mbf path selection, force=TRUE |
| `test-gpAim-extra.R` | Trait/str guards, `.cullis_H2`, `.extract_Ga_gp`, `summary.gpAim` MV, `.buildGenoData`, `.fixLines` |
| `test-multivariate.R` | 97 tests: LRT mixture helpers, `.compute_vqtilde`, `.lrtTest(ntrait>1)`, `.addEffect(Trait)`, `.packResults(Trait)`, `.waldTest()` internal, guard clauses, summary Env column |
| `test-aimTrace-extra.R` | `.build_stability_df` UV+MV, `.plot_lrt` UV+MV, `.plot_stability` UV+MV, MV aimTrace objects |
| `test-selIndex.R` | 105 tests: all 3 index types, UV path, external `selected=`, threshold variants, `.si_threshold`, `.check_index_vector`, `print/summary/plot.selIndex` |
| `test-waldTest.R` | *(deleted — waldTest is now internal `.waldTest()` in engine_results.R)* |
| `test-fineMap-mv.R` | 32 tests: MV fineMap paths (main-effect z-ratio, interaction Wald, mixed objects, aliased-coefficient guard) |
| `test-checkPanel.R` | 51 tests: return structure, encoding validation, map consistency, duplicates, missingness, MAF, het per line + marker, chr coverage, print, guards |
| `test-filterPanel.R` | 51 tests: all 8 filter steps, correct order, checkPanel dispatch, NULL thresholds, over-filter warnings, print, KNN imputation, filteredPanel dispatch, maf absent from primePanel |
| `test-integration.R` | Full-pipeline smoke tests (ASReml required; skipped if not available) |

**Total: 995 tests, 0 failures** (7 skipped: ASReml licence)

---

## Coverage (Codecov — last upload, branch `dev`)

| File | Coverage |
|---|---|
| `selIndex.R` | ~98% |
| `plot.selIndex.R` | ~100% |
| `print.selIndex.R` | ~100% |
| `engine_results.R` | ~100% |
| `engine_effects.R` | ~100% |
| `engine_genodata.R` | ~80% |
| `engine_select.R` | ~60% |
| `engine_model.R` | ~83% |
| `gpAim.R` | ~65% |
| `plot.qtlAim.R` | ~70% |
| `gwasAim.R` | ~75% |
| `aimTrace.qtlAim.R` | ~85% |
| `summary.gpAim.R` | ~85% |
| `print.gpAim.R` | ~100% |
| `plot.gpAim.R` | ~100% |
| **Overall** | **~75–82%** (target 85%) — new checkPanel/filterPanel/primePanel coverage TBD on next upload |

---

## man/ Help Pages

33 Rd files. All exported functions have dedicated help pages.

**Fixes applied (Session 30):**
- `genoRxK.Rd`, `genoSxT.Rd`: `\usage{data(phenoSxT)}` → correct dataset name
- `qtlAim.Rd` `@examples`: `tr()` → `aimTrace()`, `outStat()` removed (not in wgAim)
- `gpAim.Rd`, `primePanel.Rd`: `\link[wgAim]{primeCross}` → `\link{primeCross}`
- `wgAim-package.Rd`: expanded with selIndex, Trait extensions, supporting functions, datasets

`devtools::check_man()` — **zero warnings**.

---

## inst/CITATION

Three `bibentry` blocks (all with DOIs):
1. Taylor & Verbyla (2011) JSS — software paper
2. Verbyla, Cullis & Thompson (2007) TAG — interval mapping methods
3. Verbyla, Taylor & Verbyla (2012) Genetics Research — RWGAIM random effects

---

## _pkgdown.yml

Bootstrap 5 / flatly. Five reference sections:
- Data Preparation (`primeCross`, `fixMap`, `primePanel`)
- QTL Mapping Pipeline (`qtlAim`, print/summary/plot)
- GWAS Pipeline (`gwasAim`, print/summary/plot)
- Genomic Prediction & Selection Pipeline (`gpAim` + print/summary/plot, `selIndex` + print/summary/plot)
- Supporting Functions (`fineMap`, `aimTable`, `aimTrace`, `linkMap.*`, `getQTL`)
- Datasets (all 6)

Site URL: `https://drj001.github.io/wgAim/`

---

## Demo Scripts (scratch/)

| Script | Content |
|---|---|
| `demo_wgAim.R` | Full pipeline: primeCross → qtlAim → gwasAim → gpAim |
| `demo_effects_contrast.R` | `type="effects"` and `type="contrast"` for qtlAim + gwasAim |
| `demo_fineMap.R` | fineMap for qtlAim (interval/marker) and gwasAim |
| `demo_linkMap.R` | linkMap variants (single/multi-trait, nrow, row.chr) |
| `demo_multivariate.R` | Full MV pipeline: primeCross → qtlAim/gwasAim (Trait="Site") |
| `demo_mv_large.R` | 6-analysis demo: 200 lines, 400 markers, 4 planted QTL, single/two-trial |
| `demo_gpAim_mv.R` | univariate / bivariate (corh) / 8-trial FA(2) gpAim demo |
| `demo_selIndex.R` | All 3 index types, external selection, all 3 plot types |
| `benchmark_QTLaim.R` | Compares qtlAim vs wgaim(gen.type="interval") — ALL CHECKS PASSED |
| `benchmark_gwasAim.R` | Compares gwasAim vs wgaim(gen.type="marker") — ALL CHECKS PASSED |

---

## Session Changelog (compact)

| Session | Date | Key deliverables |
|---|---|---|
| 1–22 | 2025-03 to 2025-05 | Core package built: `qtlAim`, `gwasAim`, `gpAim`, `primeCross`, `primePanel`, `fineMap`, `aimTable`, `aimTrace`, `linkMap`, `getQTL`, `theme_scatter`. Full test suite (526 tests). All man pages. R CMD check clean. CodeCov workflow. |
| 23 | 2025-05-14 | Multivariate extensions: `Trait` + `n.fa` args for `qtlAim`/`gwasAim`; `pchisq.mixture`, `qchisq.mixture`, `.compute_vqtilde` (engine_lrt.R); `waldTest()` exported; `summary.*` Trait-aware. 62 new tests. |
| 24 | 2025-05-15 | MV engine bug fixes: `merge.by` in `.qtlSelect()`, corh/corgh grep, signed BLUPs, vqtilde floor, G.param name, classify term order, atilde ordering, Ga construction. `.constructCM()` ridge. Diag-downgrade LRT. `.packResults()` phenoData scope. `.build_stability_df()` index fix. Summary Trial1/perc.var fixes. Benchmarks. |
| 25 | 2025-05-16 | MV LRT over-detection fix (diag downgrade both sides). `aimTrace` MV support. `plot.qtlAim/gwasAim type="blups"` multi-line. |
| 26 | 2025-05-17 | MV `type="effects"` (grouped lollipop) and `type="contrast"`. `ncol`/`qtl` moved to `...`. blups flat-zero and huge-last-iter fixes. `gpAim` MV (`Trait`, upgrade loop, per-trial GEBVs). `demo_gpAim_mv.R`. |
| 27 | 2025-05-18 | `demo_gpAim_mv.R` Section 2 base model fix (diag→corh). `.extract_Ga_gp()` grep fix (`mbf.*ints`). Design principle: user controls MV structure via base model. |
| 28 | 2025-05-19 | `str` replaces `n.fa`. Critical `oint` bugfix (univariate `qtilde^2/vqtilde`). BLUP accuracy + Cullis H² in `gpAim`. MV `gebv` Trial as factor. MV caterpillar plot fix. `demo_gpAim_mv.R` rewrite. Benchmark scripts. Roxygen audit. R CMD check fixes. Merge feature/multivariate → dev → master. |
| 29 | 2025-05-20 | `selIndex()` + `print/summary/plot.selIndex`. Random formula crash fix (leading `+`, long formulas) via `.rhs_terms()`. `demo_selIndex.R`. |
| 30 | 2025-05-21 | Coverage push (74% → ~82%): 92 new tests across `test-gpAim-extra`, `test-aimTrace-extra`, `test-plot-qtlAim-extra`, `test-gwasAim-extra`. `waldTest` internalised as `.waldTest()` in `engine_results.R` (Contrasts path removed, exports removed). `asreml::` prefix removed; `globalVariables()` in `zzz.R`. Man page fixes: `\usage{}` bugs in genoRxK/genoSxT, `tr()`→`aimTrace()` in qtlAim examples, redundant `[wgAim]` qualifiers. CITATION updated (3 entries with DOIs). `wgAim-package.R` expanded (selIndex, Trait, all helpers, datasets). `_pkgdown.yml` created (Bootstrap 5). README: Active badge, selIndex in table. |
| 31 | 2025-05-22 | **fineMap MV support:** main-effect QTL use UV z-ratio; interaction QTL use `.waldTest()` joint Wald + LOD=χ²/(2log10). `is_mv`/`tgt_is_interaction` dispatch in scan loop. 32 new tests. **Version 0.9.0** set in DESCRIPTION + README (dynamic badge). **README** expanded: three-part pipeline overview (Pre-analysis/Analysis/Post-analysis); complete argument tables for primeCross, primePanel, qtlAim, gwasAim, gpAim. **Three new panel prep functions:** `checkPanel()` (pure diagnostic — encoding, map consistency, duplicates, missingness, MAF table, per-line and per-marker het, chr coverage; `print.checkPanel()`), `filterPanel()` (8-step filter in principled order: map→marker miss→line miss→het.line→het.marker→dup lines→dup markers→MAF; accepts `checkPanel` object directly; `print.filteredPanel()`), updated `primePanel()` (maf removed; accepts `filteredPanel` directly via `inherits()` check; `impute="knn"` added — chromosome-wise vectorised KNN via `tcrossprod()` BLAS + pre-sorted neighbour matrix, 2–4× faster than naive loop; `knn.k` argument; `...` removed). No S3 dispatch — plain functions. 100 new tests; `make_raw_panel()` and `make_updated_model()` moved to `helper-fixtures.R`. Total: 995 tests, 0 failures. |
| 32 | 2025-05-24 | **Vignettes (two new, CRAN-ready).** See full detail below. **`aimTrace`:** `print.out = TRUE` argument suppresses all console output when `FALSE`; 6 new tests. **`plot.fineMap`:** `ncol` argument for facet column count. **`.plot_stability()`:** `ncol` argument; wired through `aimTrace()`. **Datasets:** `genoMVsim` + `phenoMVsim` (120 DH lines, 6 chr × 100 markers, 8 trials, 4 planted QTL). **CRAN vignette infra:** `inst/doc/` pre-built HTMLs + stub `.R` files; `VignetteBuilder: quarto`; `Config/Needs/website: pkgdown`; `vignettes/.gitignore`; `build_vignette.R`. **pkgdown:** `_pkgdown.yml` articles + dataset sections; CI workflow seeds pre-built HTMLs; `--ignore-vignettes` in R CMD check. **README:** restructured (Requirements/Installation/Vignettes moved above function docs); vignette table with live links. |

---

## Vignettes

| File | Title | Status |
|---|---|---|
| `vignettes/qtlAim_pipeline.qmd` | Univariate QTL Mapping with wgAim | ✅ Complete |
| `vignettes/qtlAim_mv_pipeline.qmd` | Multivariate QTL Mapping with wgAim | ✅ Complete |
| `vignettes/gwasAim_pipeline.qmd` | Univariate GWAS with wgAim | ✅ Complete (Session 37) |
| `vignettes/gwasAim_mv_pipeline.qmd` | Multi-Environment GWAS with wgAim | ✅ Complete (Session 37) |
| GP + selection pipeline | — | ⬜ Planned |

### Univariate vignette details (Session 32)
- Expanded univariate math framework (base LMM → additive decomposition → composite model → LRT)
- `linkMap()` 3-row with marker names; QTL map with marker names, height ×1.2
- Outlier/blups height ×1.3; contrast `ncol=2`; stability `ncol=2`; aimTrace `print.out=FALSE`
- `fineMap()` top-4 QTL (`"Chr.3B.30"`, `"Chr.3D.4"`, `"Chr.3B.10"`, `"Chr.4D.2"`) in 2×2 grid
- Removed: print(qtlAim), extracting results section (getQTL/aimTable)
- `#| results: markup` on `qtlAim_class` chunk so engine/QTL count renders

### Multivariate vignette details (Session 32)
- **Math:** MV LMM → FA($k$) structure for $G_a$ → vm/mbf paths → chi-squared mixture LRT → quadratic-form outlier stat → MAIN vs G×E Wald pruning
- **Data:** `genoMVsim` / `phenoMVsim` — 120 DH lines, 6 chr × 100 markers, 8 trials, 2 reps; 2 main-effect QTL (Chr1/Chr3, effects 0.18/−0.15) + 2 G×E QTL (Chr2 fades, Chr5 crossover, max 0.30/0.25); `data-raw/simulate_MVsim.R`
- **Base model development:** `diag(Trial):id` → `fa(Trial,1):id` → `fa(Trial,2):id` with AIC table
- **qtlAim():** `Trait="Trial"`, `str=NULL` mirrors FA(2) base, `TypeI=0.01` (raised from 0.05 to suppress near-threshold spurious QTL with reduced main effects)
- **Plots:** outlier, multi-trial BLUPs, grouped effects lollipop, contrast `ncol=2`, aimTrace LRT + stability `ncol=2`
- **Root cause of cache bug:** stale 4-trial cache passed to 8-trial engine → dimension mismatch; clearing cache resolved it

---

## New/Changed R/ Files (Session 32)

| File | Change |
|---|---|
| `R/aimTrace.qtlAim.R` | `print.out = TRUE` arg; `ncol` arg wired to `.plot_stability()` |
| `R/aimTrace.gwasAim.R` | `print.out = TRUE` arg |
| `R/plot.fineMap.R` | `ncol` arg for facet column count |
| `R/genoMVsim.R` | New — roxygen docs for `genoMVsim` and `phenoMVsim` datasets |

---

## New Files (Session 32)

| File | Description |
|---|---|
| `vignettes/qtlAim_pipeline.qmd` | Univariate QTL vignette source |
| `vignettes/qtlAim_mv_pipeline.qmd` | Multivariate QTL vignette source |
| `vignettes/custom.css` | Custom CSS for both vignettes |
| `vignettes/.gitignore` | Excludes `*.html`, `*_cache/`, `*_files/` from git |
| `inst/doc/qtlAim_pipeline.html` | Pre-built UV vignette (CRAN) |
| `inst/doc/qtlAim_pipeline.R` | Stub R file (CRAN requirement) |
| `inst/doc/qtlAim_mv_pipeline.html` | Pre-built MV vignette (CRAN) |
| `inst/doc/qtlAim_mv_pipeline.R` | Stub R file (CRAN requirement) |
| `data/genoMVsim.rda` | Simulated DH cross object (120 lines, 6 chr × 100 markers) |
| `data/phenoMVsim.rda` | Simulated phenotypic data (8 trials, 2 reps, 1920 obs) |
| `data-raw/simulate_MVsim.R` | Reproducible data generation script |
| `build_vignette.R` | Renders vignettes and copies HTML to `inst/doc/` |
| `.github/workflows/pkgdown.yaml` | pkgdown CI: seeds both pre-built HTMLs, individual build fns |

---

## pkgdown Site

- URL: `https://drj001.github.io/wgAim/`
- Articles: Univariate QTL Mapping, Multivariate QTL Mapping
- Workflow seeds `inst/doc/*.html` → `docs/articles/` before `build_articles(lazy=TRUE)`
- `--ignore-vignettes` + `--no-build-vignettes` in R CMD check workflow

---

## New/Changed R/ Files (Session 33)

| File | Change |
|---|---|
| `R/aimTrace.qtlAim.R` | `.build_stability_df()`: fix ASReml coefficient name lookup — interaction coefficients named `"FactorCol_Level:X.chr.idx"` (was missing `"Trait_"` prefix, causing only reference trial to show). `.plot_stability()` MV: replaced `facet_wrap` + 8 coloured overlapping error-bar lines with `facet_grid(trial ~ qtl_label_fac)` + `geom_ribbon`; colour/legend dropped; `ncol` arg ignored for MV. |
| `R/plot.qtlAim.R` | `.add_sig_labels()`: replaced `geom_point` diamond + `geom_text_repel` nudge with dashed `geom_vline` + `geom_text(y=Inf)` — matches `.add_mv_blups_flags()` style; affects all outlier/blups plots for UV+MV `qtlAim` and `gwasAim`. |
| `R/engine_genodata.R` | `.fixLines()`: `Gomit` fixed-effect term is now `Trait:Gomit` when `Trait` is non-NULL, so the omitted-line correction is applied within each trial level. |
| `R/engine_lrt.R` | New `.buildLRTModels()`: shared diag-downgrade logic extracted from the inline loop code. Uses `.rhs_terms()` for formula robustness; `baseForLRT` now preserves ALL random terms from `baseModel` (was silently dropping spatial/block terms). New `.buildLRTModels()` also added to `engine_lrt.R` (previously in `engine_lrt.R` as blank file). |
| `R/engine_model.R` | `.buildGenomeModel()`: upgrade loop removed — now returns `vm.struct`, `upgrade.to.fa`, `gterm`, `resid.term` in addition to existing fields. New `.upgradeVmStructure()`: extracted upgrade logic, called from `qtlAim`/`gwasAim` only AFTER initial LRT confirms signal. |
| `R/engine_results.R` | `coef.list[[iter]]` indexing corrected (was `iter-1L`, now `iter` consistent with new convention where `iter == n_qtl` at loop exit). |
| `R/engine_select.R` | `.fixLines()` call updated with `Trait = Trait`. |
| `R/gpAim.R` | `.fixLines()` call updated with `Trait = Trait`. |
| `R/gwasAim.R` | Phase 4 fully restructured: pre-loop LRT on `diag(Trait):vm`; if fails → message + return; if passes → `.upgradeVmStructure()`; loop with LRT at bottom; `lrt.stop` local variable for combined break; `cat()` + `message()` for all LRT steps. |
| `R/qtlAim.R` | Same restructuring as `gwasAim.R`; also removed leftover debug `print()` statements. |
| `vignettes/qtlAim_mv_pipeline.qmd` | `aimtrace_stability` chunk: `fig-height` raised from 6.5 → 9.8; caption updated for grid layout. `contrast` plot: `ncol=2` → `ncol=3`. Stability call: `ncol=2` dropped. LaTeX fix: `$\hat{\beta}_j` → `$\hat{\beta}_j$`. `build_vignette.R` updated to include MV vignette. |
| `vignettes/qtlAim_pipeline.qmd` | Prose fixes: heading "and forward selection and" → "and forward selection"; "every every" → "every". |
| `inst/doc/qtlAim_pipeline.html` | Rebuilt with prose fixes + unified annotation style. |
| `inst/doc/qtlAim_mv_pipeline.html` | Rebuilt with all Session 33 changes. |

---

## Key Architecture Changes (Session 33)

### Forward-selection loop restructuring
**Old:** `repeat { qtlSelect → LRT → if fails break → add QTL → iter++ }` — LRT at top, QTL only added when LRT passed.

**New:**
1. Pre-loop: build `diag(Trait):vm` only → initial LRT → if fails: message + return; if passes: upgrade to FA/corgh via `.upgradeVmStructure()`
2. `repeat { qtlSelect → add QTL → LRT at bottom → if (!pass | breakout) break → iter++ }`
3. `iter == n_qtl` at loop exit; `coef.list[[iter]]` is the correct indexing.

This ensures FA/corgh models are only built when there is evidence of genome-wide signal, avoiding expensive model fitting over noise.

### .buildLRTModels() — shared downgrade helper
Extracts the MV `diag(Trait):` downgrade logic used by both the pre-loop LRT and every bottom-of-loop LRT. Key fix: `baseForLRT` now uses `.rhs_terms(baseModel$call$random)` and replaces only the residual genetic term — all other random terms (spatial, block, Rep) are preserved. The previous inline code replaced the entire random formula with a single term.

### ASReml coefficient naming — MV stability plot fix
ASReml names factor:covariate interaction coefficients as `"FactorCol_Level:X.chr.idx"` (column name + underscore + level). `.build_stability_df()` was constructing `paste0(tname, ":", xk)` (missing the `"Trait_"` prefix), so `which(names(cof) == int.nm)` returned `integer(0)` for every non-reference trial, collapsing all per-trial lines to the reference-trial effect only (single colour in plot).

### Stability plot layout
MV stability plot changed from `facet_wrap(~ qtl, ncol=n)` with 8 overlapping coloured error-bar lines per facet to `facet_grid(trial ~ qtl_label_fac)` with one `geom_ribbon` per cell. Rows = trials, columns = QTL with [MAIN]/[INT] suffix. `fig-height` raised to 9.8 in vignette.

### Outlier/blups annotation style unified
`.add_sig_labels()` in `plot.qtlAim.R` now uses dashed `geom_vline` + `geom_text(y=Inf)` (matching `.add_mv_blups_flags()`). Affects: `plot.qtlAim type="outlier"` (UV+MV), `plot.qtlAim type="blups"` (UV), `plot.gwasAim type="outlier"` (UV+MV), `plot.gwasAim type="blups"` (UV).

### User-facing progress messages
`cat()` headers before each LRT model fitting step (go to trace file).
`message()` after each `.lrtTest()` call: `"LRT: statistic = X.XXXX, p-value = Y.YYYY"` — always visible on screen regardless of `trace=` setting. Stopping message appended: `"  -- No further QTL."` or `"  -- Breakout."`.

---

## Session Changelog (compact)

| Session | Date | Key deliverables |
|---|---|---|
| 1–22 | 2025-03 to 2025-05 | Core package built: `qtlAim`, `gwasAim`, `gpAim`, `primeCross`, `primePanel`, `fineMap`, `aimTable`, `aimTrace`, `linkMap`, `getQTL`, `theme_scatter`. Full test suite (526 tests). All man pages. R CMD check clean. CodeCov workflow. |
| 23 | 2025-05-14 | Multivariate extensions: `Trait` + `n.fa` args for `qtlAim`/`gwasAim`; `pchisq.mixture`, `qchisq.mixture`, `.compute_vqtilde` (engine_lrt.R); `waldTest()` exported; `summary.*` Trait-aware. 62 new tests. |
| 24 | 2025-05-15 | MV engine bug fixes: `merge.by` in `.qtlSelect()`, corh/corgh grep, signed BLUPs, vqtilde floor, G.param name, classify term order, atilde ordering, Ga construction. `.constructCM()` ridge. Diag-downgrade LRT. `.packResults()` phenoData scope. `.build_stability_df()` index fix. Summary Trial1/perc.var fixes. Benchmarks. |
| 25 | 2025-05-16 | MV LRT over-detection fix (diag downgrade both sides). `aimTrace` MV support. `plot.qtlAim/gwasAim type="blups"` multi-line. |
| 26 | 2025-05-17 | MV `type="effects"` (grouped lollipop) and `type="contrast"`. `ncol`/`qtl` moved to `...`. blups flat-zero and huge-last-iter fixes. `gpAim` MV (`Trait`, upgrade loop, per-trial GEBVs). `demo_gpAim_mv.R`. |
| 27 | 2025-05-18 | `demo_gpAim_mv.R` Section 2 base model fix (diag→corh). `.extract_Ga_gp()` grep fix (`mbf.*ints`). Design principle: user controls MV structure via base model. |
| 28 | 2025-05-19 | `str` replaces `n.fa`. Critical `oint` bugfix (univariate `qtilde^2/vqtilde`). BLUP accuracy + Cullis H² in `gpAim`. MV `gebv` Trial as factor. MV caterpillar plot fix. `demo_gpAim_mv.R` rewrite. Benchmark scripts. Roxygen audit. R CMD check fixes. Merge feature/multivariate → dev → master. |
| 29 | 2025-05-20 | `selIndex()` + `print/summary/plot.selIndex`. Random formula crash fix (leading `+`, long formulas) via `.rhs_terms()`. `demo_selIndex.R`. |
| 30 | 2025-05-21 | Coverage push (74% → ~82%): 92 new tests across `test-gpAim-extra`, `test-aimTrace-extra`, `test-plot-qtlAim-extra`, `test-gwasAim-extra`. `waldTest` internalised as `.waldTest()` in `engine_results.R` (Contrasts path removed, exports removed). `asreml::` prefix removed; `globalVariables()` in `zzz.R`. Man page fixes: `\usage{}` bugs in genoRxK/genoSxT, `tr()`→`aimTrace()` in qtlAim examples, redundant `[wgAim]` qualifiers. CITATION updated (3 entries with DOIs). `wgAim-package.R` expanded (selIndex, Trait, all helpers, datasets). `_pkgdown.yml` created (Bootstrap 5). README: Active badge, selIndex in table. |
| 31 | 2025-05-22 | **fineMap MV support:** main-effect QTL use UV z-ratio; interaction QTL use `.waldTest()` joint Wald + LOD=χ²/(2log10). `is_mv`/`tgt_is_interaction` dispatch in scan loop. 32 new tests. **Version 0.9.0** set in DESCRIPTION + README (dynamic badge). **README** expanded: three-part pipeline overview (Pre-analysis/Analysis/Post-analysis); complete argument tables for primeCross, primePanel, qtlAim, gwasAim, gpAim. **Three new panel prep functions:** `checkPanel()`, `filterPanel()`, updated `primePanel()`. 100 new tests. Total: 995 tests, 0 failures. |
| 32 | 2025-05-24 | **Vignettes (two new, CRAN-ready).** `aimTrace`: `print.out = TRUE` arg. `plot.fineMap`: `ncol` arg. `.plot_stability()`: `ncol` arg. **Datasets:** `genoMVsim` + `phenoMVsim`. **CRAN vignette infra:** `inst/doc/` pre-built HTMLs, `build_vignette.R`. **pkgdown** articles + CI. README restructured. |
| 33 | 2025-05-24 | **Stability plot fix** (ASReml coef name `Trait_Level:X.chr.idx`); **facet_grid layout** for MV stability; **unified annotation style** (vline+top-label) across all outlier/blups plots; **`Trait:Gomit` fix** in `.fixLines()`; **loop restructuring**: pre-loop LRT → conditional FA/corgh upgrade → LRT at bottom of loop; **`.buildLRTModels()`** shared helper (preserves all random terms); **`.upgradeVmStructure()`** deferred upgrade; **LRT progress messages**; vignette prose + height fixes. 997 tests, 0 failures. |
| 34 | 2025-05-25 | **MV outlier statistic fixed** (oint now ~9+ for true QTL, was ≤2). Two bugs: (1) eigenvalue regularisation of `Ga` in FA branch distorted `Ginv` — removed, now uses `MASS::ginv()` Moore-Penrose pseudoinverse directly; (2) pure-R `.compute_vqtilde()` replaced with compiled **Rcpp/RcppArmadillo** version from old `mvwgaim` code. First Rcpp code in wgAim: `src/compute_vqtilde.cpp`, `src/Makevars[.win]`, `DESCRIPTION`/`NAMESPACE` updated. Also fixed `predict(only=vmterm.raw)` for MV path and removed debug `print()` calls. 997 tests, 0 failures. |
| 35 | 2025-05-26 | **GWAS vignettes (UV + MV):** Univariate `gwasAim_pipeline.qmd` refreshed for new 300-line × 1500-marker × 6-site sim; site names Env1..Env6 (width=1); prose/counts updated; fineMap keys made dynamic (`head(qtl, 2L)`); plot heights reduced 20% (manhattan/outlier/blups 12→9.6, contrast 9→7.2); stability changed to 2×2 (`ncol=2`); math section compressed to reference qtlAim vignette. **New multi-environment vignette** `gwasAim_mv_pipeline.qmd` following qtlAim_mv structure (diag → fa1 → fa2 base progression, `Trait="Site"`, `str=NULL` mirror). **Two package bug fixes:** (1) `aimTrace.gwasAim()` was missing the `ncol` argument in its signature and not threading it to `.plot_stability()` — added and wired; (2) `.build_mv_contrast_df()` hardcoded `classify = "Trait:id"` but `.fixLines()` substitutes `id → Gsave` when phenotype lines lack marker data — factored out `res_lab` to use `"Gsave"` whenever the fallback triggers. **Sim tuning experiments:** GWAS sim `sigma_genetic` reduced 100 → 30 to cut polygenic background (LODs jumped 30–40% cleaner); QTL effects halved (170→85, 180→90, 130→65, 250→125 fade, ±180→±90) — but composite LRT fell to 0.10 (below χ²₁ threshold) so **zero QTL detected**. Diagnosis: composite LRT is a *marginal* test (`vm(id,K)` vs bare `id` random) — halving quarters K's marginal signal beyond `id` and drops it below detection. Currently the halved sim is on disk; user to decide whether to revert. **Knitr cache gotcha:** chunk-source hashing doesn't invalidate downstream chunks when an upstream chunk's parameters change; the fix for parameter changes is `Remove-Item vignettes/<name>_cache` before re-render. 997 tests, 0 failures. |
| 36 | 2025-05-27 | **GWAS sim rewritten to use AlphaSimR** after long diagnostic session on why forward selection kept picking spurious "Chr1 marker 1". Root cause: hand-rolled independent-draw sim gave K essentially = I (off-diag SD ~0.026 from finite-sample noise only), so vm(id, K) was redundant with id and outlier statistics were all zero → R's `which.max(rep(0, p))` returned index 1 → first marker on Chr1 systematically picked. **Intermediate attempts** all failed and are documented for posterity: (a) removed biomAid's IID polygenic and added `U ~ MVN(0, Ga ⊗ K)`, K-shape held LRT high forever (couldn't be exhausted by finite QTL); (b) removed U per user insight (K-signal from QTL alone), but LRT wouldn't fire because K ≈ I; (c) added rank-3 PC latent structure → K became rank-few, all markers projected onto same 3 directions, oint values indistinguishable, first-index tie-break; (d) increased to rank-30 PC — better but still artificial. **Final fix:** AlphaSimR `runMacs2(nInd=50, Ne=100, genLen=1.5)` + 1 gen `randCross` → 300 F1s + `addSnpChip(nSnpPerChr=150)`. K now has realistic full-rank spectrum (off-diag SD 0.085, top-5 eigenvalues 6.4–8.1, ~12% of total). Response model kept simple: `y = μ + Site:Rep + Σ β_jk X_ik + L_ij + ε` with L IID across lines (no U term). **QTL doubled** to 340/−360/260/500-fade/±360 kg/ha (7–11% of mean) after diagnostic showed 5/6 sites had clear vm variance (12k–59k, z-ratios 0.94–2.99) but Env3 hit boundary (0.023 with bound `B`) — MASS::ginv then inflated Ginv[3,3] to ~43 and cascaded into vqtilde ≤ 0 for every marker, poisoning oint. **Discovered `str = "corh"` recommendation** for MV GWAS: pools to 1 variance + 1 correlation, no per-site boundary risk regardless of specific sim draw. **Files changed:** `data-raw/simulate_GWASsim.R` (complete rewrite; AlphaSimR-based), `R/genoGWASsim.R` (roxygen for AlphaSimR generation, doubled QTL table, dropped 3-subpop language, MVN(0, Σ_L)⊗I polygenic wording), `data/*.rda` regenerated, `scratch/smoke_gwas_new_sim.R` and `scratch/smoke_gwas_mv_new_sim.R` (diagnostic smoke tests). **Package engine untouched** — user's guidance "trust qtlAim works" was correct; every problem was in the sim. No test changes needed. **Vignettes unchanged and now stale** (hardcoded old marker names). |
| 37 | 2026-08-03 | **GWAS sim calibration marathon** (~10 iterative rounds with the user) to hit "all 5 planted QTL detected, correct MAIN/INT classification, zero spurious hits" on the final MV smoke test. Scope shrunk to **200 lines** (from 300) after the user reverted to an older reliable ASReml build that couldn't handle the larger fit. Key changes to `data-raw/simulate_GWASsim.R`: (a) **yield scale converted kg/ha → t/ha** so variance components sit near unit scale (better REML gradient behaviour); (b) **two polygenic components** now generated — U ~ MVN(0, Ga_K ⊗ K) as trace K-weighted (target diag 0.005, deliberately small so LRT stops after QTL absorption) + L ~ MVN(0, Ga_L ⊗ I) as moderate IID-line (target diag 0.020, keeps fa(Site,k):id well-conditioned); (c) **QTL effects settled** at Chr2 +0.145, Chr5 −0.145, Chr7 +0.125 (weak), Chr3 +0.230→0 (G×E fade), Chr9 ±0.180 (G×E cross) after progressive tuning informed by |β|/SE ≈ 4-5 rule to keep MAIN QTL classified as MAIN by the Wald pruning; (d) **`find_qtl_marker()` MAF floor 0.15 → 0.25** with improved highest-MAF fallback (previously "closest marker regardless"); (e) **`addSnpChip(minSnpFreq = 0.15)`** for active MAF-based panel selection (was: random from all segregating sites, which oversampled rare variants that cause spurious GWAS hits); (f) **F1 MAF > 0.15 post-filter** for drift; (g) `nSnpPerChr` 150 → 175 to buffer for post-filter drops; (h) **marker renumbering moved to AFTER filter** so `SNP_CXX_NNNN` names stay contiguous; (i) **`Ga_L` seed 1729 → 1730** — the 1729 Wishart draw put Env3 diagonal 2.6× above target (0.052 vs 0.020) causing persistent Chr5 MAIN → INT misclassification because Env3 estimates wandered. `R/genoGWASsim.R` roxygen fully rewritten to match. Final shipped panel: ~1691 markers, all MAF > 0.15 (median 0.32). **Package engine still untouched** — all changes are in the sim + roxygen. `scratch/smoke_gwas_mv_new_sim.R` header refreshed. Vignettes remain stale (next work). |
| 38 | 2026-11-06 | **MV engine overhaul for corh/corgh/us and the mbf path**, triggered by a colleague's 3-trial `corgh(Trial):id` gwasAim crash (`'data' must be of a vector type, was 'NULL'`). Reproduced on simulated 3/4/5-trial data (shipped panel genotypes + 9 ungenotyped lines). Fixes: (1) **`predict(only=)` term name** — raw formula string never matches non-fa terms (ASReml strips `corgh()`/`us()`/`diag()`/`corh()` wrappers: key is `Trait:vm(id, covObj)`); now uses the G.param key → fixes the crash and silent all-zero outlier stats (Chr1.1 picks). (2) **`.extractVmGa()`** — name-based Ga for diag/corh/corgh/us/fa(k); fixes corgh correlations misplaced for ≥ 4 trials (ASReml row-wise vs old column-wise fill — colleague's report, confirmed with 4/5-trial fits vs `us` reference), corh `!Trial!cor` unrecognised (NA Ga), us falling back to wrong diagonal; also replaces gpAim's `.extract_Ga_gp()` internals. (3) **`us` mirrors to `us`** (probe showed us() works on vm()/mbf() in ASReml 4.2 and converges where corgh aborts/stalls for ≥ 4 trials); docs/README recommend `us(Trait):id`. (4) **gpAim never upgraded MV structure** since S33 — now calls `.upgradeVmStructure()`. (5) **gpAim MV mbf branch** (previously fell into UV code). (6) **ASReml NA mbf-key bug** (ungenotyped lines → wrong "converged" fit): `Gmbf` key + zero `.NoGeno` row. (7) **MV mbf outlier stats** always zero (diagonal-only PEV) → exact method via predicted line values + Moore–Penrose `M⁺`. (8) **Formulae held in variables** (`random = rf`) → `.resolveCallFormulae()`. (9) **MV Perc.Var** structure-dependent → diag of fitted covariances. (10) `.buildLRTModels()` skips needless cold-start refit. 16 new tests. Committed `707ab65` (dev) / merged `47ddf25` (main), pushed; `deb1f8e` divergence resolved. |

---

---

## New/Changed Files (Session 34)

| File | Change |
|---|---|
| `src/compute_vqtilde.cpp` | **New** — Rcpp/RcppArmadillo implementation of `compute_vqtilde()` for the MV outlier statistic denominator |
| `src/RcppExports.cpp` | **New** — auto-generated by `Rcpp::compileAttributes()` |
| `src/Makevars` / `src/Makevars.win` | **New** — LAPACK/BLAS linker flags for RcppArmadillo |
| `R/RcppExports.R` | **New** — auto-generated R wrapper calling `_wgAim_compute_vqtilde` |
| `R/engine_select.R` | Removed FA eigenvalue regularisation of `Ga`; replaced `.compute_vqtilde()` call with `compute_vqtilde()` (Rcpp); fixed `only = vmterm.raw` for MV `predict()`; removed debug `print()` calls; comments tightened |
| `R/engine_lrt.R` | Comments tightened; `.compute_vqtilde()` retained for tests |
| `DESCRIPTION` | Added `Rcpp (>= 1.0.0)` to `Imports`; added `Rcpp, RcppArmadillo` to `LinkingTo` |
| `NAMESPACE` | Added `useDynLib(wgAim, .registration=TRUE)` and `importFrom(Rcpp, sourceCpp)` |

---

## Key Architecture Changes (Session 34)

### MV outlier statistic — two bugs fixed

**Bug 1 — Eigenvalue regularisation of Ga (the main culprit)**
The FA branch in `.qtlSelect()` floored small eigenvalues of `Ga` at `0.01 × max_eigenvalue` before computing `Ginv`. For `fa(Trait,k)` models where specific variances are at boundary zero, `Ga = Lam Lam^T` is legitimately rank-deficient (rank `k < ntrait`). The regularisation modified `Ga` by up to O(eps) per element (e.g. 157 units on a real dataset) producing a wrong `Ginv` and deflating the outlier statistics to near-null values (oint ≤ 2 instead of ~9). Fixed by removing the regularisation — `MASS::ginv(Ga)` computes the correct Moore-Penrose pseudoinverse matching the paper's `Ga^-`.

**Bug 2 — Pure-R vqtilde denominator**
The pure-R `.compute_vqtilde()` in `engine_lrt.R` was being used for the MV vm path, producing numerically different results from the Rcpp/Armadillo version in the old working `mvwgaim` code. Fixed by integrating `compute_vqtilde.cpp` from `scratch/vqtilde/` into the package as an Rcpp compiled function.

### Rcpp integration
First compiled code in wgAim. `src/compute_vqtilde.cpp` uses RcppArmadillo to compute `tr(Ginv %*% Ti %*% vatilde %*% t(Ti))` for all markers in a tight C++ loop. `LAPACK_LIBS`, `BLAS_LIBS`, and `FLIBS` linked via `src/Makevars[.win]`.

---

## New/Changed Files (Session 35)

| File | Change |
|---|---|
| `vignettes/gwasAim_pipeline.qmd` | Refreshed for new sim: 300 lines / 1500 markers / 6 sites (was 400/3000/10); `Env01` → `Env1`; fineMap keys made dynamic via `head(gwas_env1$QTL$qtl, 2L)`; plot heights reduced 20% (manhattan/outlier/blups 12 → 9.6, contrast 9 → 7.2); stability set to 2×2 (`ncol=2`); Mathematical Framework section compressed to reference the same section in qtlAim vignette with a short "What's different in GWAS" subsection covering panel-marker source, population structure absorption by K, and vm-path dominance; empty `path` field cat removed (gwasAim doesn't set `QTL$path`, only gpAim does); Results narrative rewritten for actual detection pattern |
| `vignettes/gwasAim_mv_pipeline.qmd` | **New** — Multi-environment GWAS vignette. Structure mirrors `qtlAim_mv_pipeline.qmd`: Overview → Math (referenced) → Data Preparation (checkPanel/filterPanel/primePanel condensed, referring UV vignette) → Base Model Development (diag → fa(Site,1) → fa(Site,2) with AIC table, parameter counts computed dynamically) → `gwasAim(Trait="Site", str=NULL)` mirroring FA(2) → print/summary → 6-plot suite (manhattan default for MV, outlier, blups multi-site, effects, contrast MV, heatmap) + linkMap → aimTrace lrt/stability |
| `R/aimTrace.gwasAim.R` | **Bug fix**: added `ncol = NULL` to signature and threaded through to `.plot_stability()` in both `plot == "stability"` and `plot == "both"` branches. Sibling `aimTrace.qtlAim()` had it; this one didn't — the `ncol` argument was silently ignored |
| `R/plot.qtlAim.R` | **Bug fix** in `.build_mv_contrast_df()`: factored out `res_lab` variable. Was hardcoding `classify_term <- paste(Trait, gterm, sep = ":")` — but when `.fixLines()` triggers (phenotype lines without marker data), the random genetic term becomes `fa(Trait,k):Gsave` not `fa(Trait,k):id`. Now uses whichever label the residual grep found (either `gterm` or `"Gsave"`), and reads the same column from `pv$pvals`. qtlAim MV vignette avoids this because `genoMVsim`/`phenoMVsim` have matched line sets; gwasAim MV hit it because filtering drops some lines |
| `data-raw/simulate_GWASsim.R` | `sigma_genetic` 100 → 30 (rationale: lower polygenic background lets planted QTL stand out cleanly); QTL effects halved (170→85, 180→90, 130→65, 250→125 fade, ±180→±90). **Halving overshot** — see Session 35 architecture note; may need to revert to intermediate values (e.g. divide by 1.3 giving 130/-140/100/±140) |
| `R/genoGWASsim.R` | Roxygen updated: 400→300 lines, 3000→1500 markers, ~0.5→~1 cM spacing, λ=5→3 cM; correct QTL marker names (`SNP_C02_0046`, `SNP_C05_0085`, `SNP_C07_0029`, `SNP_C03_0111`, `SNP_C09_0065`); sigma_genetic and effect tables aligned with new sim |
| `scratch/smoke_gwasAim_mv.R` | QTL effect comments and truth-recovery labels updated to match new sim (85, -90, 65, 125, ±90) |
| `scratch/peek_cache.R` | **New** — small utility to `lazyLoad()` a knitr chunk cache and inspect a fitted object's slots (used to diagnose the "0 QTL after halving" issue) |
| `build_vignette.R` | Added `gwasAim_mv_pipeline.qmd` to the render list |

---

## Key Architecture Changes (Session 36)

### GWAS sim: AlphaSimR replaces hand-rolled construction
`data-raw/simulate_GWASsim.R` now uses `AlphaSimR::runMacs2()` for founder haplotype generation (coalescent-based, realistic LD and MAF spectrum) and `randCross()` for 1 generation of random mating to produce 300 F1 diverse lines from 50 founders. The resulting K matrix has genuine full-rank structure with gradually decaying eigenvalue spectrum — a stark improvement over the previous hand-rolled independent-draw approach (K ≈ I) or rank-few PC hack (K rank-deficient). `AlphaSimR` is a new dependency for `data-raw/` only (not for the package itself). Marker names on the panel drift between sim runs because the coalescent places segregating sites randomly — `find_qtl_marker()` still picks the closest well-behaved marker to each target cM position, so ground-truth QTL positions are stable but marker indices are not.

### Diagnostic: per-site variance boundary poisons MV `oint`
When any single site's `vm(id, K)` variance estimate hits the zero boundary in `str = "diag"` (or any structure with per-site variance parameters), `MASS::ginv(Ga)` produces a Ginv element ~43 for that site (1/0.023) — orders of magnitude larger than other Ginv entries. That inflated element propagates through `vqtilde = tr(Ginv Ti vatilde Ti')`, driving vqtilde ≤ 0 for every marker. The `pmax(vqtilde, 0)` floor + `ifelse(vqtilde > 0, ..., 0)` guard then zeroes the entire `oint` vector. `which.max(rep(0, p))` returns index 1 → "Chr.Chr1.1" is systematically picked regardless of true QTL positions. Recommendation: use `str = "corh"` for MV GWAS which pools to a single variance parameter and cannot boundary.

### Response model: no additive polygenic term
Simulation deliberately omits `U ~ MVN(0, Ga ⊗ K)`. Only sources of between-line variance are (1) planted QTL fixed contributions and (2) `L_i· ~ MVN(0, Σ_L)` drawn independently across lines (compound-symmetric across sites, ρ=0.6). L is IID across lines with zero K-projection by construction, so the K-shaped signal that `vm(id, K)` detects comes ENTIRELY from the planted QTL through their marker columns being in the range of K. Forward selection can then genuinely exhaust the K-signal by moving planted QTL to fixed effects, unlike sims with continuous U where LRT never dropped.

---

## New/Changed Files (Session 37)

| File | Change |
|---|---|
| `data-raw/simulate_GWASsim.R` | Substantial rewrite: `n.lines` 300 → 200; `n.mar.pc` 150 → 175; `n.rows.rep` 15 → 10; `maf.floor <- 0.15` named constant; **yield scale kg/ha → t/ha** (site means 4.5–4.9, residual 0.142, ~0.377 SD/site); **U + L polygenic** (Ga_K trace 0.005 + Ga_L moderate 0.020, with Wishart-generated 6×6 covariances); `addSnpChip(minSnpFreq = 0.15)` for active MAF panel selection; **F1 MAF > 0.15 post-filter** for drift; marker renumbering moved to AFTER filter (contiguous `SNP_CXX_NNNN`); `find_qtl_marker()` MAF floor 0.15 → 0.25 + `maf.floor` arg + highest-MAF fallback; **Ga_L seed 1729 → 1730** (Env3 fix); QTL effects settled at +0.145 / −0.145 / +0.125 / +0.230→0 / ±0.180 t/ha; diagnostic MAF distribution + per-chr count printouts |
| `R/genoGWASsim.R` | Roxygen fully rewritten: panel size (200 lines × ~1700 markers, `SNP_CXX_NNNN` contiguous naming); "All markers MAF > 0.15" statement in both `genoGWASraw` and `mapGWASsim`; U + L polygenic decomposition explained (trace K-weighted + moderate IID-line); LRT-stopping calibration explained (Ga_K deliberately small so LRT drops below threshold after QTL absorption); QTL effect table updated to t/ha scale with note on β/SE ≈ 4-5 rule for MAIN classification; shipped QTL marker names refreshed |
| `data/genoGWASraw.rda` | Regenerated (200 × ~1691 columns, all F1 MAF > 0.15) |
| `data/mapGWASsim.rda` | Regenerated (~1691 rows, ~170 per chr) |
| `data/phenoGWASsim.rda` | Regenerated (2400 rows: 6 sites × 200 lines × 2 reps, yield 2.9–6.4 t/ha, mean ~4.7) |
| `man/genoGWASraw.Rd`, `mapGWASsim.Rd`, `phenoGWASsim.Rd` | Auto-regenerated via `devtools::document()` |
| `scratch/smoke_gwas_mv_new_sim.R` | Header comment refreshed to reflect current sim (200 lines, both U+L polygenic components, current QTL effects in t/ha, updated Yield scale note) |

---

## Key Architecture Changes (Session 37)

### Sim calibration philosophy: separate signal from noise via QTL/polygenic ratio
The core insight from this session: the composite genome-wide LRT stops firing when the residual K-shape variance falls below significance. For clean forward-selection stopping, the sim needs:
- Enough K-shape signal from QTL contributions that the LRT reliably fires at iteration 1 through K (else no detection)
- Small enough non-QTL polygenic (`Ga_K`) that once QTL are absorbed as fixed effects, the residual K-signal drops below LRT threshold (else spurious hits)

Solution: keep Ga_K deliberately trace-level (target diag 0.005). The K-shape signal for the LRT to detect comes primarily from the QTL marker columns being in K's range. Ga_L (IID-line, target diag 0.020) is separate — it doesn't project onto K in expectation but keeps `fa(Site,k):id` well-conditioned so the analysis's FA loadings don't compete with `vm(id, K)`.

### Active MAF selection replaces passive rejection
Low-MAF markers (MAF < 0.15) are the primary source of spurious GWAS hits: their marker score has high variance from few minor-allele carriers, they don't align well with K structure, and the outlier statistic doesn't fully correct for their extreme frequency. Fix applied at two layers:
1. **AlphaSimR `addSnpChip(minSnpFreq = 0.15)`** — samples the SNP chip only from founder segregating sites with MAF > 0.15 (default is 0, i.e. random from all segregating sites)
2. **F1 MAF > 0.15 post-filter** — recomputes MAF on the F1 panel after 1 gen of `randCross` and drops any that drifted below the floor

Result: shipped panel MAF distribution min 0.15, median 0.32, max 0.50 — dominated by informative common variants.

### `find_qtl_marker()`: MAF stability over cM precision
Raised the internal MAF floor from 0.15 → 0.25 for the QTL-position-to-marker selector. Trade-off: QTL marker can now sit a few cM further from the planted target (e.g. Chr7 lands on SNP_C07_0042 which is a bit past 30 cM if the closest are all low-MAF), but the per-QTL variance contribution `2·p·(1−p)·β²` stays stable across regenerations. Before: a "lucky" draw could put Chr7 on a MAF 0.18 marker (var 0.004) or "unlucky" on MAF 0.45 (var 0.010) — 2.5× swing in what the LRT sees. After: all QTL markers reliably land at MAF > 0.25 (var ≥ 0.006) so re-runs are reproducible.

### Ga_L seed selection matters
The 6×6 Ga_L covariance is generated Wishart-style with `n.sites + 2 = 8` degrees of freedom, so per-site diagonal draws have appreciable sampling variability. `set.seed(1729)` drew a matrix where Env3's diagonal came out 2.6× above target (0.052 vs 0.020), which persistently made Env3 observations noisier than the other sites. This did not affect the `vm(id, K)` LRT (L is IID → zero K-projection), but it did make per-site QTL estimates at Env3 noisy — causing MAIN QTL on Chr5 to be **misclassified as G×E** because the Env3 effect estimate wandered from the true constant while other sites stayed near it. `set.seed(1730)` gives a flatter draw (Env3 diag 0.027, closer to target 0.020) and eliminates the Chr5 misclassification.

### Yield scale change: kg/ha → t/ha
All response-scale numbers converted (site means 4700 → 4.7; QTL effects 340 → 0.34 → later tuned to 0.145; residual variance 142000 → 0.142). Reason: variance components sitting near unit scale behave better in ASReml's REML convergence heuristics (per-parameter tolerances are calibrated around unit scale). Analysis is scale-invariant in principle but numerically more robust in practice.

---

## Key Architecture Changes (Session 35)

### `.build_mv_contrast_df` — Gsave label propagation
When some phenotype lines have no marker data, `.fixLines()` substitutes `merge.by → Gsave` in the random genetic term. The affected term becomes e.g. `fa(Site,2):Gsave` instead of `fa(Site,2):id`. Previously `.build_mv_contrast_df` had a fallback for finding the residual G.param entry (grep `Gsave` if `":<gterm>$"` was not found) but continued to use the original `gterm` name for the `predict(classify = paste(Trait, gterm, sep=":"))` call. ASReml then reported "id is not in the model, unable to predict using it". Fix: introduce `res_lab` variable that holds the actually-in-model term name (`"Gsave"` or `gterm`) and use it for both `only` filter and `classify` term. qtlAim MV never hit this because its `genoMVsim` matches `phenoMVsim` exactly.

### `aimTrace.gwasAim()` `ncol` threading
`.plot_stability()` (shared with qtlAim) accepts `ncol`. `aimTrace.qtlAim()` had `ncol = NULL` in its signature and passed it through. `aimTrace.gwasAim()` did not — silently ignored user-supplied `ncol`, always used `min(n_qtl, 3L)`. Simple additive fix (add to signature, thread to both call sites).

### LOD vs LRT distinction (documentation / user-education)
Extensive discussion resulted in a clearer understanding worth recording:
- The reported LOD in `summary.gwasAim()` is a **Wald z-ratio** on the post-pruning coefficient, computed as `LOD = z² / (2 ln 10)` where `z = coef / SE`. The formula is correct. For MV MAIN effects, however, SE reflects the pooled sample size across all sites × lines (effective n ~ 1750 in our sim vs ~960 in qtlAim MV), so LODs come out ~3–5× larger than for the equivalent biparental analysis — this is a **feature of larger n**, not stronger signal.
- The actual **detection** is driven by the **composite genome-wide LRT** — a REML variance-component test comparing `base: y ~ ... + id` against `full: y ~ ... + id + vm(id, K)`. This is a MARGINAL test: it asks "does K-structured between-line variance fit better than unstructured `id`?" It's not "is there between-line variance?" — `id` already models that.
- Consequence: a QTL with a big pooled-n LOD in the summary does not imply a big LRT statistic. Halving QTL effects quarters K's marginal signal against `id`, and (as observed) can drop the LRT below threshold even when individual Wald z-ratios would still be large.
- The `sigma_genetic = 30` + halved-effects combination gave initial LRT stat = 0.10 (p = 0.37) — total detection failure. Diagnosed via a small `scratch/peek_cache.R` script that `lazyLoad()`s the knitr chunk cache and prints `lik.list`.

### Knitr chunk-cache invalidation
Discovered a real gotcha: when a chunk's source code changes (e.g. changing `TypeI = 0.01` to `TypeI = 0.05` in the `gwasAim` chunk), knitr invalidates and re-runs that chunk, but **downstream chunks that consume the changed object keep their stale cache** because their own source code hasn't changed and knitr's default caching is chunk-source-hash-based (no `autodep`). In the MV vignette this manifested as `gwasAim` reporting 4 detected QTL while `print_gwas` / `summary_gwas` still showed the old 3-QTL output from the previous `TypeI=0.01` run. Workaround: `Remove-Item -Recurse -Force vignettes/<name>_cache` whenever an upstream chunk's parameters change. Considered adding `autodep: true` to the vignette execute block as a safety net; deferred.

### GWAS MV recovery pattern
With `sigma_genetic = 30` and original (unhalved) QTL effects, the multi-environment analysis recovered 4/5 planted QTL: three MAINs (Chr2 +170 → LOD 76, Chr5 −180 → LOD 57, Chr7 +130 → LOD 35 — the weakest, only picked up thanks to MV borrowing strength) and Chr9 correctly classified as INT crossover (per-site LODs 12–17). Chr3 fade was not recovered: its trial-averaged effect is +125 but the mixture-χ² joint test doesn't get much signal from `Env6` where the effect is 0. The `.build_mv_contrast_df` fix above was uncovered while trying to render the MV vignette's contrast plot.

---

## New/Changed Files (Session 38)

| File | Change |
|---|---|
| `R/engine_select.R` | `.qtlSelect()`: MV `predict(only=)` uses G.param key `vmterm.only` (+ explicit error if `pv$vcov` is NULL); Ga via new `.extractVmGa()`; **MV mbf branch rewritten** (exact: `predict(classify = "Trait:<mbf key>", only = mbf G.param key, vcov = TRUE)`, drop `.NoGeno`, `trans = MASS::ginv(M)` with dimnames, prior `kronecker(Ga, tcrossprod(M))`, Rcpp `compute_vqtilde()`). New helpers `.extractVmGa(model, Trait, trait.levels, sigma2, term)`, `.mbfEffectRows(rn, Trait, trait.levels, markers)`, `.mvPercVarComponents(object, gterm, scale, sigma2)` |
| `R/engine_model.R` | New `.mbfKey()` / `.mbfCovObj()` (Gmbf key + all-zero `.NoGeno` row); used in `.buildGenomeModel()` and `.rebuildCovObj()` (incl. vm→mbf switch); `attr(intervalObj, "mbf.geno")` set on mbf branches (NULL on vm) for `.qtlSelect()`; **us residual mirrors to `us`** (was corgh); comments rewritten |
| `R/engine_genodata.R` | `.fixLines()` adds `phenoData$Gmbf` (Gsave with NA → `".NoGeno"`) |
| `R/engine_lrt.R` | `.buildLRTModels()` returns `qtlModel` unchanged when the LRT random formula is identical (no cold-start `G.param = NULL` refit) |
| `R/engine_validate.R` | New `.resolveCallFormulae()` (restore symbol/expression call components from `model$formulae`), called first in `.validateModel()` |
| `R/gpAim.R` | `.resolveCallFormulae()` at validation; **`.upgradeVmStructure()` call after `.buildGenomeModel()`** (MV); **new MV mbf branch** (per-trial GEBV/SE/Accuracy, `Ga`, `Gcor`, `gen.H2` NA); `.extract_Ga_gp()` delegates to `.extractVmGa()`; `str` roxygen rewritten |
| `R/fineMap.R` | mbf rescan uses the fitted model's mbf key via `.mbfCovObj()` |
| `R/summary.qtlAim.R`, `R/summary.gwasAim.R` | MV `var.mark`/`var.res` from `.mvPercVarComponents()` |
| `R/gwasAim.R`, `R/qtlAim.R` | `str` roxygen: each structure mirrors to itself; recommend `us` over `corgh` for ≥ 4 groups |
| `man/gwasAim.Rd`, `man/qtlAim.Rd`, `man/gpAim.Rd` | Hand-edited to match roxygen (no roxygen2 on this machine) — run `devtools::document()` to confirm no diff |
| `README.md` | Three `str` table rows: recommend `"us"` / `us(Trait):id` for ≥ 4 groups |
| `tests/testthat/test-extractVmGa.R` | **New** — 8 tests, mock names copied from real ASReml 4.2 fits (diag, corh, corgh 3-trial, corgh 4-trial scrambled order, us, fa(2) with fixed loading, error paths) |
| `tests/testthat/test-engine-helpers.R` | **New** — 8 tests: `.resolveCallFormulae`, `.mbfKey`/`.mbfCovObj`, `.mbfEffectRows` (diag/us, fa + Comp rows, Trait name with `_`, gap error) |
| `tests/testthat/test-gpAim-extra.R` | Mocks use real ASReml names (`!Trial_T1`, corh `!Trial!cor`); corh off-diagonal asserted |
| `tests/testthat/test-multivariate.R` | us-residual test renamed; asserts `vm.struct == "us"` |
| `scratch/` (git-ignored) | `repro_corgh_3trial.R` (main repro, args `[struct] [str\|NULL] [n.trials]`), `check_corgh_order.R` (4/5-trial ordering proof), `probe_us_genomic.R`, `probe_gpAim_upgrade.R`, `probe_fixes_symbol_mbf.R` (A: symbol formulae, B: gpAim vm vs mbf, C: gwasAim vm vs mbf; reduced 255-marker panel), `probe_mbf_predict.R`, `probe_mbf_gsave_fa.R`, `debug_*.R`, `run_tests_shim.R`, `rlib/` (test install), `commit_msg.txt` |

---

## Key Architecture Changes (Session 38)

### ASReml term naming (ASReml-R 4.2.0.481) — the root of most bugs
- **G.param keys / predict `only=`:** non-fa wrappers are stripped: `diag|corh|corgh|us(Trial):vm(id, covObj)` → `"Trial:vm(id, covObj)"`; `fa(Trial,k):vm(...)` keeps `"fa(Trial, k):vm(id, covObj)"`. A raw formula string passed as `only=` silently fails to match (full fitted values returned, or no vcov on some builds). **Always use the G.param key.**
- **vparameter suffixes** (after the term key, i.e. after the first `!`): diag `!Trial_T1`; corh `!Trial!cor` + `!Trial_T1…`; corgh `!Trial!T2:!Trial!T1.cor` … (lower triangle **row-wise**: (2,1)(3,1)(3,2)(4,1)…) + variances; us `!Trial_T1:T1, !Trial_T2:T1, !Trial_T2:T2, …` (row-wise incl. diagonal); fa `!T1!var` (Ψ) + `!T1!fa1`… (Λ, constrained zero loadings included). `.extractVmGa()` places every parameter by name.
- **MV mbf coefficients:** `Trial_T1:mbf("ints")_<marker>` (trial-major); fa: `fa(Trial, 1)_T1:mbf(...)_<marker>` plus `fa(Trial, 1)_Comp1:...` factor-score rows. `.mbfEffectRows()` matches by level + marker name.

### corgh vs us
Same covariance model; identical REML optimum when both converge (verified 3 trials, vm and mbf). With ≥ 4 trials `corgh` (bounded correlations) frequently aborts ("too many exceptions") or stops at a local optimum reported as converged (4-trial example: logL 239.3 vs 265.4 for us), producing negative LRTs. `us` on the genomic term is supported (old code comment was wrong). Recommendation documented: `us(Trait):id` base for unstructured covariance with ≥ 4 groups. Note: two unstructured 4×4 genetic matrices (vm + residual) can still be weakly identified on modest data — `str = "corh"` / `fa(Trait,1)` are the fallbacks.

### mbf path with ungenotyped lines
ASReml mishandles NA mbf keys ("Internal: Unexpected number of random labels, observed 1338, expected 1341" — one per Trait level) and returns a wrong fit flagged converged. Fix: key column `Gmbf` (`.NoGeno` for ungenotyped) + all-zero `.NoGeno` covariate row ⇒ zero genomic contribution, exactly as Gsave = NA on the vm path. After the fix: mbf logL 175.4386 = vm 175.439; gpAim GEBVs vm vs mbf cor = 1.0000 for diag/us/fa(1); gwasAim vm vs mbf identical QTL and LRTs.

### MV mbf outlier statistic (exact)
Needs per-marker full ntrait×ntrait posterior block; `vcoeff` is diagonal only. Now: predict line values `u` of the mbf term (`classify = Trait:Gmbf`; classify on `Gsave` is rejected by ASReml) with full vcov; `q = M⁺ u` (Moore–Penrose; exact for any n/p since q lies in row space of M — verified max diff 7e-6); prior `Ga ⊗ MM'`; then the vm machinery (`compute_vqtilde`). Marker matrix passed via `attr(intervalObj, "mbf.geno")`.

### gpAim MV
Now applies the mirrored/requested structure (was always diag since S33). MV mbf branch: `GEBV_ij = M_i q_j`, SE from diagonal PEV (approximate, lower accuracy than vm — same approximation as UV mbf), `var.genetic_j = Ga_jj × mean(rowSums(M²))`, `gen.H2` NA.

### Base model call hygiene
`model$call$random` etc. hold the *symbol* when the formula is passed as a variable (works only after `update()`). `.resolveCallFormulae()` replaces any non-`~` call component with `model$formulae[[nm]]` (plain formula) at validation — `update()` reproduces the fit exactly.

---

## What Remains Open

| Item | Status |
|---|---|
| Cross-validation for `gpAim` | ⬜ Open |
| GWAS vignette (univariate) | ✅ Done Session 35; regenerated + deployed Session 37 |
| GWAS vignette (multivariate) | ✅ Done Session 35; regenerated + deployed Session 37 |
| GP + selection vignette | ⬜ Planned — likely next session |
| `dev` / `main` history divergence (`deb1f8e` audit fix on main only) | ✅ Resolved Session 38 — local `dev` created at `main` (fast-forward of `origin/dev`) before committing; `dev` now contains all of `main` |
| **Run full `devtools::test()` + `devtools::document()`** on a machine with testthat/roxygen2 | ⬜ Open (S38) — only shim-tested here (492 pass, 4 shim artefacts); `man/*.Rd` were hand-edited, `document()` should produce no diff |
| **Check GitHub CI** (test-coverage, pkgdown) for `707ab65` / `47ddf25` | ⬜ Open (S38) — `gh` not installed on this machine |
| **Re-render MV vignettes** (`qtlAim_mv_pipeline`, `gwasAim_mv_pipeline`) — Perc.Var values change with the S38 fix (fa(2) models) | ⬜ Open (S38) — clear caches first; refresh `inst/doc/*.html` |
| **`demo_gpAim_mv.R`** — results change now that gpAim applies the MV structure upgrade | ⬜ Open (S38) — re-run and check |
| Tell colleague: fixed in `47ddf25`; for ≥ 4 trials use `us(Trial):id` (not corgh) in the base model; fallbacks `str = "corh"` / `fa(Trial,1)` | ⬜ Open (S38) |
| **GWAS sim effect-size tuning** — halved effects gave 0/5 detection (LRT stat 0.10) | ✅ Resolved Session 36 (AlphaSimR + doubled QTL) then further tuned Session 37 (200 lines, trace Ga_K, moderate Ga_L, effects settled at +0.145/−0.145/+0.125/+0.230→0/±0.180 t/ha, MAF > 0.15 panel selection) |
| **Re-render both GWAS vignettes** with new sim + updated QTL effects + current marker names (`SNP_C02_0044`, `SNP_C05_0094`, `SNP_C07_0042`, `SNP_C03_0128`, `SNP_C09_0068`); UV `fineMap` chunks hardcode old `SNP_C03_0111` — make dynamic or update. Also update prose for 200-line panel + ~1700 marker count + t/ha yield scale + Wishart Ga_K + Ga_L polygenic description. | ✅ Done Session 37 — both UV and MV vignettes regenerated against S37 sim; fineMap chunks already dynamic (top-2 from `gwas_env1$QTL$qtl`); README + pkgdown promotions applied; committed and deployed to `main` |
| **Vignette recommendation for MV GWAS `str` argument** — Session 36 diagnosis showed `str = "diag"` risks per-site boundary artifacts; document that `str = "corh"` is the robust default | ⬜ Downgraded — with the small Ga_K and MAF > 0.15 panel from Session 37, `str = NULL` (mirror base model FA(2)) works cleanly on the smoke test. Recommendation can stay as an optional fallback note. |
| **Guard `plot.gwasAim()` against empty QTL** for manhattan/outlier/blups/heatmap types | ⬜ Open — user rejected edit mid-flight; helpful error message drafted but not applied |
| Consider `autodep: true` in vignette YAML | ⬜ Open — safety net for future param changes |
| Last iteration outlier/blups scan (post final QTL) not shown in plots | ⬜ Open — deferred |
| MV outlier statistic scale investigation | ✅ Fixed Session 34 — Ga regularisation removed + Rcpp vqtilde |
| Validate `demo_multivariate.R` end-to-end with ASReml | ⬜ Open |
| Validate bivariate `demo_gpAim_mv.R` with ASReml | ⬜ Open (see S38 row above — results now change) |
| `cran-comments.md` + `inst/WORDLIST` | ⬜ Pre-CRAN |
| `NEWS.md` (README used as interim) | ⬜ Pre-CRAN |
| Push coverage above 85% | ⬜ Ongoing |
