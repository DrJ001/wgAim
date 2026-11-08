# wgAim — AOM-based Stopping Rule for Forward Selection

*Working notes, started Session 39 — 2026-11-07/08.*
*Main project reference: [session_summary.md](../../scratch/session_summary.md) (local working file in the git-ignored `scratch/` folder, so not on GitHub). Maths write-up: [aom_stopping_rule.qmd](aom_stopping_rule.qmd) (rendered: [aom_stopping_rule.html](aom_stopping_rule.html)).*

> **Status:** research prototype on branch `feature/aom-stopping-rule`. **No package code has been changed.** Everything lives in `research/aom_stopping/` (excluded from the package build via `.Rbuildignore`). The diagnostic is attached to `wgAim:::.qtlSelect()` for a single R session by the run scripts and leaves selection unchanged. Deleting this folder removes it completely.
>
> **Running the scripts:** from the repository root, e.g. `Rscript research/aom_stopping/run_stage1.R gwas 1 1 500 2001`. They need ASReml-R and an installed wgAim; `scratch/rlib` is put first on `.libPaths()` if present (the local test install). Results are written to `research/aom_stopping/results/`.

---

## 1. Problem

A colleague's 3-trial oat `gwasAim()` run (11,991 markers, `us(Trial):id`) worked after the S38 fixes, but the **LRT statistic rose** as markers were selected, even though the selected markers looked right:

| Iter | L0 | L1 | Statistic |
|---|---|---|---|
| 1 | 1377.80 | 1452.87 | 150.15 |
| 2 | 1381.19 | 1459.28 | 156.18 |
| 3 | 1374.66 | 1455.11 | 160.91 |
| 4 | 1371.62 | 1453.97 | 164.70 |
| 5 | 1379.86 | 1463.74 | 167.75 |
| 6 | 1377.20 | 1463.93 | 173.47 |

Five markers (2C, 4C, 6A, 6C, 7D) explained about 28% of the genetic variance in total. Effects ran from −0.34 to +0.34, with consistent signs across trials. The map is physical (Mb), so `exclusion.window = 20` means ±20 Mb, and the `dist(cM)` header in `summary()` is misleading for such maps.

## 2. Diagnosis of the LRT

- **Reproduced** on simulated 3-trial `corgh` data with strong K-structured background (`research/aom_stopping/repro_corgh_3trial.R`): the LRT went 48.6 → 19.9 → 21.8 → 23.7 → 28.3.
- **Not a convergence problem.** Both LRT models were refitted to full convergence and gave the same statistics to 3 d.p. (`research/aom_stopping/diag_lrt_rise.R`). The recurring "Log-likelihood not converged" warnings come from `predict(..., maxit = 1)` in `.qtlSelect()`, not from model fits.
- **Not a built-in REML artefact.** Adding null markers to both models moved the LRT both ways (−3.8 to +1.6).
- **What the LRT tests.** `L0` (no genome-wide term) is compared with `L1` (`diag(Trait):vm`), testing H₀: per-trial vm variances = 0, i.e. *does K-structured polygenic variance exist?* With a strong background that's always true, so the LRT never stops. Fitting large QTL (poorly aligned with K, otherwise absorbed by the IID line term) can even *raise* the statistic. With a weak background the LRT can stop too early: one simulation stopped after 1 of 3 planted QTL.
- **How the LRT is run for `us`.** The vm term is fitted as `diag` at entry. After the entry LRT passes it is upgraded to `us`, which is used *only* for marker selection. At each iteration, `.buildLRTModels()` downgrades it back to `diag` and refits from default starting values (`G.param = NULL, R.param = NULL`). `L0` keeps `us(Trial):Gsave` plus the QTL fixed effects and is updated from previous estimates.
- **`pchisq.mixture()` review** (`research/aom_stopping/chibar_weights.R`):
  - The binomial weights assume the per-trial variance estimators are uncorrelated. With correlation 0.3–0.6 the true size is about 5.5–6.3% at nominal 5%.
  - The Stram & Lee citation in the code is inappropriate; the correct reference is Self & Liang.
  - `pchisq(0, df = 0)` returns 0, a harmless edge case.
  - User decision: **no changes for now.**

## 3. The idea: use the AOM outlier statistic to stop

Keep the LRT as the **entry gate** (is there any genetic signal?). At each iteration, test whether the **best marker is exceptional relative to the polygenic background**, and stop when it isn't.

### Theory (details in `aom_stopping_rule.qmd`)

- **AOM null:** marker *j*'s effect is just another polygenic draw, H₀ⱼ: Ψⱼ = 0 in aⱼ ~ N(0, Gₐ/s + Ψⱼ). It is **not** aⱼ = 0, which is vacuous under a polygenic model. That also explains why post-selection Wald p-values are all < 0.001.
- **Use the marginal BLUP variance, not the PEV:** var(q̃ⱼ) = prior − PEV (`vqtilde` in the package). This is the correct reference for testing a variance-component hypothesis.
- **The outlier statistic is the REML score test** for the AOM, evaluated at the null fit.
  - UV: Oⱼ = q̃ⱼ²/vⱼ ~ χ²₁ exactly, given the variance parameters.
  - MV: Oⱼ = q̃ⱼ′Gₐ⁻q̃ⱼ / tr(Gₐ⁻Vⱼ), which is the score for *proportional* inflation Ψⱼ = ψⱼGₐ/s (my derivation; to be checked against Verbyla et al.). Its null is Σλₖχ²₁ with mean 1, not χ²; p-values come from Imhof.
  - Standardised alternative: Sⱼ = q̃ⱼ′Vⱼ⁻q̃ⱼ ~ χ²ᵣ (a different alternative, and possibly different selection).
- **Family-wise null of the maximum, simulated exactly without refitting:** q̃ = (I ⊗ T)g̃ with g̃ ~ N(0, Vg), Vg = Gₐ ⊗ K − PEV. Draw g̃* = L z (pivoted Cholesky), recompute every Oⱼ*, take the maximum, and compare. This handles LD, MV weighting and exclusion windows. Bonferroni is the cheap fallback.
- **Fixed-effect entry is the limiting case** ψⱼ → ∞, so the forward procedure is a coherent step-down.
- **Proposed procedure:** LRT entry gate → each iteration: select the best marker → p_FWE ≤ α? → add it as fixed and continue; otherwise stop. This could become an option such as `stop.rule = c("lrt", "outlier")`.

## 4. Implementation (`research/aom_stopping/`)

| File | Purpose |
|---|---|
| `aom_null_diag.R` | `aomNullDiag()`: vm and mbf paths, UV and MV, marker and chromosome statistics, per-marker p-values (χ²₁ / Imhof), Bonferroni, simulated max-stat null, Sⱼ, M_eff, timings. Helpers `.aom_inputs()`, `.imhof_upper()`, `.sqrt_factor()` (pivoted Cholesky, eigen fallback). |
| `run_aom_diag.R` | Pilot with planted QTL; `TypeI = 0.999`, `breakout = 8`. Args: structure, number of trials, background scale, breakout. |
| `run_aom_null.R` / `summ_aom_null.R` | Pilot null runs: 50 replicates per background, 3 trials, `us`. |
| `run_stage0.R` | Stage 0 verification harness: 11 engine-path cases. |
| `bench_aom_scaling.R` | Scaling benchmark, no ASReml needed. |
| `run_stage1.R` / `summ_stage1.R` | Stage 1 null calibration. Args: `<pop> <ntrait> <kscale> <nrep> <seed0> [struct] [B]`; output `results/stage1/*.csv`. |
| `check_imhof.R` | Imhof vs Monte Carlo / exact checks. |
| `aom_stopping_rule.qmd` / `.html` | Maths write-up. Render with `"C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe" render research/aom_stopping/aom_stopping_rule.qmd` (Quarto isn't on the PATH). |

Practicalities:
- Fit and evaluate in the **global environment**: `gwasAim`/`qtlAim` assign `covObj`/`phenoData` to their caller.
- The UV entry LRT can't be forced open with `TypeI = 0.999`, because `qchisq(1 − 2·TypeI, 1)` is `NaN`. Use 0.49.
- `qtlAim()` returns `lik.mat = NULL` when the entry gate stops; `gwasAim()` still returns it.
- Don't edit a script while Rscript is running it: Rscript reads the file as it goes.

## 5. Results so far

### Pilot with planted QTL (3 trials, `us`, 191 lines, 1,501 markers; planted Chr8 crossover, Chr2/Chr5 main effects)

| Background | Current LRT | AOM rule (simulated threshold) |
|---|---|---|
| Weak (×0.05) | stops after **1** QTL | finds **all 3**, stops immediately after (p = 0.0005, 0.0005, 0.002, then 0.97) |
| Strong (×1) | **never stops** | stops after **2** (misses Chr5: per-marker p = 8e-4, but not exceptional among 1,501 markers); all spurious steps p ≥ 0.37 |

### Stage 0: verification across engine paths (`run_stage0.R`)

- **Exact match** with the package's `oint` (≤ 1e-11): all vm cases (gwas/qtl, UV/MV, `us`/`fa1`, interval/marker), all MV mbf cases, and the chromosome statistic (UV and MV).
- **UV mbf** differs by up to 3–13% per marker (correlation 0.999999, same marker selected). The cause is the diagnostic's extra `predict(maxit = 1)` iteration versus the package's direct use of `coefficients`/`vcoeff`. It doesn't affect validity (the statistic and its null come from the same fit).
- **Forcing mbf on 1,501 markers is impractically slow**, so mbf is tested on thinned panels (lines ≥ markers).
- **Effective tests:** GWAS panel ≈ 90% of markers; DH intervals ≈ 30% (169/594 UV, 323/594 MV). Bonferroni is therefore 2–3.5× conservative for `qtlAim`.
- **Scaling** (reference BLAS):
  - pivoted Cholesky of Vg: 14 s at N = 4,800 (eigen-decomposition: 138 s);
  - simulation: about 10 s per trial at 12,000 markers × 800 lines (B = 2,000);
  - diagnostic ≤ 7 s per step at the test sizes.

### Stage 1: null calibration (2,600 replicates, B = 1,000, `summ_stage1.R`)

| Cell | Replicates | Simulated threshold | Bonferroni | Sⱼ | Chromosome statistic | Full procedure | Effective tests |
|---|---|---|---|---|---|---|---|
| DH UV, weak | 500 | 0.2% | 0.0% | 0.0% | 0.0% | 0.0% | 119 |
| DH UV, strong | 500 | 2.0% | 1.2% | 1.2% | 0.6% | 2.0% | 261 |
| DH 3-trial, weak | 188 | 2.1% | 0.0% | 0.0% | 0.0% | 0.5% | 142 |
| DH 3-trial, strong | 200 | 3.0% | 1.5% | 0.5% | 0.5% | 3.0% | 257 |
| GWAS UV, weak | 500 | 1.6% | 1.4% | 1.4% | — | 0.4% | 1197 |
| GWAS UV, strong | 500 | 3.2% | 2.4% | 2.4% | — | 3.2% | 1181 |
| GWAS 3-trial, weak | 95 | 4.2% | 4.2% | 2.1% | — | 3.2% | 1124 |
| GWAS 3-trial, strong | 100 | 3.0% | 4.0% | 3.0% | — | 3.0% | 1083 |

- **Never above nominal.** It is conservative, and most conservative where there are fewest effective tests (DH).
- **The lower tail is light:** p-values average about 0.5, but P(p ≤ 0.10) is 0.035–0.078. The KS test rejects in the UV cells.
- **Likely cause (hypothesis):** the plug-in REML estimate of σ²ₐ/Gₐ satisfies Σ q̃ⱼ² ≈ Σ var(q̃ⱼ), so the observed statistics are self-normalised, while the simulated null treats Gₐ as known.
- Simulated null mean slightly above 1 in UV (1.007–1.018) and 1.045 for `fa1` (Stage 0), probably from the fitted variance parameters and the `predict(maxit = 1)` PEV not quite matching.
- **Base-model failures:** 17 `us` base-model failures in the weak 3-trial cells ("not positive definite"), unrelated to the rule.

### Earlier pilot null (`run_aom_null.R`, 50 replicates, 3-trial GWAS)
False-positive rate 4.0% (strong) and 4.1% (weak); p-value mean 0.50 / 0.47.

## 6. Open issues and next steps

| Item | Status |
|---|---|
| **Studentised (self-normalised) null:** rescale each simulated draw so Σ num* = Σ den, mimicking REML pinning. Quick test on the DH UV cells (~10 min per 500 replicates). MV needs a per-Gₐ-element version. | ⬜ Proposed |
| Check the eigenvalues of Vg for `fa(k)` / rank-deficient Gₐ (null mean 1.045) | ⬜ |
| **Stage 2: power** with planted QTL (graded main effects, fade, crossover, linked pair). Compare current LRT / simulated / Bonferroni / oracle; QTL recovered, false positives, MAIN/INT classification, step-down behaviour | ⬜ Next after the calibration fix |
| **Stage 3:** current Oⱼ (proportional alternative) vs Sⱼ, especially for G×E QTL | ⬜ |
| **Stage 4:** parametric null on real genotypes (`genoRxK`, `genoSxT`, `genoMVsim`, GWAS panel, colleague's oat data) | ⬜ |
| Other cells: 6–8 trials `fa(2)`, `diag` with a boundary trial, unbalanced trials, ungenotyped lines in DH, spatial residuals | ⬜ |
| Run the colleague's data with the diagnostic attached (`source()` + wrapper) | ⬜ |
| Check the MV score derivation against the AOM literature; add the original AOM reference to the qmd | ⬜ |
| Add an "Inference for random effects" section to the qmd; add the Stage 0/1 results to the qmd | ⬜ |
| If adopted: extend `compute_vqtilde()` to return full t×t blocks Vⱼ; add a `stop.rule` option to `qtlAim`/`gwasAim`; chromosome-selection variant; tests | ⬜ Not started (package untouched) |
| Side items found along the way: `summary()` `dist(cM)` label for physical maps; `pchisq.mixture` weights/citation/`df = 0`; LRT L1 refitted from default starting values every iteration | ⬜ Deferred by user |
