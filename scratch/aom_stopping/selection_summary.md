# wgAim — MOST: an AOM-based Stopping Rule for Forward Selection

**MOST — Maximum Outlier Score Test.** Working name for the new rule. The goal is a **standalone** stopping rule that runs end-to-end with **no LRT anywhere** (it decides both whether to start and when to stop). The composite LRT is kept as the **legacy** alternative, both in the package and for head-to-head comparisons in the planned paper. Verbyla's result stands: the LRT controls type I error for its own hypothesis (no K-structured variance), but that is the wrong hypothesis for stopping when there is a polygenic background.

*Working notes, started Session 39 — 2026-11-07/08; last updated 2026-11-09 (Stages 1b–1d, golden rule).*
*Main project reference: [session_summary.md](../session_summary.md). Maths write-up: [aom_stopping_rule.qmd](aom_stopping_rule.qmd) (rendered: [aom_stopping_rule.html](aom_stopping_rule.html)).*

> **Status:** research prototype on branch `feature/aom-stopping-rule`. **No package code has been changed.** Everything lives in `scratch/aom_stopping/` (tracked in git; `scratch/` is excluded from the package build via `.Rbuildignore`). Run logs are in `logs/`, results in `results/`. The diagnostic is attached to `wgAim:::.qtlSelect()` for a single R session by the run scripts and leaves selection unchanged. Deleting this folder removes it completely.
>
> **Running the scripts:** from the repository root, e.g. `Rscript scratch/aom_stopping/run_stage1.R gwas 1 1 500 2001`. They need ASReml-R and an installed wgAim; `scratch/rlib` is put first on `.libPaths()` if present (the local test install). Results are written to `scratch/aom_stopping/results/`.

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

- **Reproduced** on simulated 3-trial `corgh` data with strong K-structured background (`scratch/repro_corgh_3trial.R`, Session 38): the LRT went 48.6 → 19.9 → 21.8 → 23.7 → 28.3.
- **Not a convergence problem.** Both LRT models were refitted to full convergence and gave the same statistics to 3 d.p. (`scratch/aom_stopping/diag_lrt_rise.R`). The recurring "Log-likelihood not converged" warnings come from `predict(..., maxit = 1)` in `.qtlSelect()`, not from model fits.
- **Not a built-in REML artefact.** Adding null markers to both models moved the LRT both ways (−3.8 to +1.6).
- **What the LRT tests.** `L0` (no genome-wide term) is compared with `L1` (`diag(Trait):vm`), testing H₀: per-trial vm variances = 0, i.e. *does K-structured polygenic variance exist?* With a strong background that's always true, so the LRT never stops. Fitting large QTL (poorly aligned with K, otherwise absorbed by the IID line term) can even *raise* the statistic. With a weak background the LRT can stop too early: one simulation stopped after 1 of 3 planted QTL.
- **How the LRT is run for `us`.** The vm term is fitted as `diag` at entry. After the entry LRT passes it is upgraded to `us`, which is used *only* for marker selection. At each iteration, `.buildLRTModels()` downgrades it back to `diag` and refits from default starting values (`G.param = NULL, R.param = NULL`). `L0` keeps `us(Trial):Gsave` plus the QTL fixed effects and is updated from previous estimates.
- **`pchisq.mixture()` review** (`scratch/aom_stopping/chibar_weights.R`):
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

## 4. Implementation (`scratch/aom_stopping/`)

| File | Purpose |
|---|---|
| `aom_null_diag.R` | `aomNullDiag()`: vm and mbf paths, UV and MV, marker and chromosome statistics, per-marker p-values (χ²₁ / Imhof), Bonferroni, simulated max-stat null, Sⱼ, M_eff, timings. Helpers `.aom_inputs()`, `.imhof_upper()`, `.sqrt_factor()` (pivoted Cholesky, eigen fallback). |
| `run_aom_diag.R` | Pilot with planted QTL; `TypeI = 0.999`, `breakout = 8`. Args: structure, number of trials, background scale, breakout. |
| `run_aom_null.R` / `summ_aom_null.R` | Pilot null runs: 50 replicates per background, 3 trials, `us`. |
| `run_stage0.R` | Stage 0 verification harness: 11 engine-path cases. |
| `bench_aom_scaling.R` | Scaling benchmark, no ASReml needed. |
| `run_stage1.R` / `summ_stage1.R` | Stage 1 null calibration. Args: `<pop> <ntrait> <kscale> <nrep> <seed0> [struct] [B]`; output `results/stage1/*.csv`. |
| `check_imhof.R` | Imhof vs Monte Carlo / exact checks. |
| `diag_lrt_rise.R` | LRT-rise diagnosis: per-step L0/L1, full-convergence refits, null-marker experiment. |
| `debug_stage1_dh3.R` | Debug of the qtlAim gate-stop case (`lik.mat = NULL`). |
| `summ_stud.R` | Stage 1b: original vs studentised null (UV cells), paired seeds. Results `results/stage1_stud/`. |
| `summ_stud_mv.R` / `check_stud_mv.R` | Stage 1c: scalar vs matrix studentisation, MV cells at several α. Results `results/stage1_stud{,_gA,_gB}/`. |
| `check_stud_boundary.R` | Stage 1c: calibration split by whether the REML pin holds. |
| `probe_weak_mv_pin.R` | Stage 1c: why the pin fails in weak MV (non-converged `us` vm fit drifting to Gₐ = 0). Arg: seed. |
| `run_most_null.R` / `run_most_null_v2.R` / `summ_most_null.R` | Stage 1d: **standalone MOST** null grid (no LRT anywhere; legacy LRT computed on the same fits for comparison). Args `<pop> <ntrait> <kscale> <lscale> <nrep> <seed0> [B] [out.sub]`. v2 = MV base falls back `us` → `diag` also on non-convergence. Results `results/most_null{,_v2}/`, summary `results/most_null_summary.txt`. |
| `logs/` | Console logs of all runs. |
| `aom_stopping_rule.qmd` / `.html` | Maths write-up. Render with `"C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe" render scratch/aom_stopping/aom_stopping_rule.qmd` (Quarto isn't on the PATH). |

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

### Stage 1b: studentised null, UV (`summ_stud.R`; same seeds as Stage 1)

**Cause of the Stage 1 conservatism confirmed:** at the REML fit the score equation for the genome-wide variance gives Σⱼ numⱼ = Σⱼ denⱼ exactly (observed ratio 1.0000 in every replicate). The unconditioned simulated null ignores this, so its maximum is too dispersed and its threshold too high. **Fix:** rescale each simulated draw so that Σ num\* = Σ den. This lowers the threshold by 5–13%, most in DH.

| Cell | Original | **Studentised** | KS p (original → studentised) |
|---|---|---|---|
| DH UV, weak (500; 128 gate stops) | 0.2% | **2.4%** (3.2% of ok reps) | <0.001 → 0.12 |
| DH UV, strong (500) | 2.0% | **5.0%** (3.3–7.3) | <0.001 → 0.42 |
| GWAS UV, strong (500) | 3.2% | **4.8%** (3.1–7.1) | 0.04 → 0.99 |

UV is solved.

### Stage 1c: scalar vs matrix studentisation, MV (3 trials, `us` vm; `summ_stud_mv.R`, `check_stud_boundary.R`)

For `us`, REML pins **every element** of Gₐ: Σⱼ q̃ⱼq̃ⱼ′ = Σⱼ Vⱼ (t(t+1)/2 constraints). The **scalar** version imposes only the Gₐ⁻-weighted trace; the **matrix** version transforms each draw so that its cross-product matrix equals Σⱼ Vⱼ.

Pooled over the four MV cells (1,942 replicates):

| α | Original | Scalar | Matrix | Bonferroni |
|---|---|---|---|---|
| 0.01 | 0.006 | 0.009 | **0.010** | 0.004 |
| 0.05 | 0.028 | 0.038 | **0.049** (0.040–0.060) | 0.016 |
| 0.10 | 0.060 | 0.082 | **0.093** | 0.033 |
| 0.20 | 0.137 | 0.166 | **0.196** | 0.065 |

By cell:
- **Strong background:** the pin holds in 91–96% of replicates. **The matrix version is calibrated** (DH 5.0%, GWAS 5.4%; KS p 0.34 / 0.38); the scalar version stays conservative (DH 3.4%).
- **Weak background:** the pin fails in ~90% of replicates. Neither version is exact: DH is conservative (3.2%) and GWAS slightly liberal with the matrix version (6.1%).
- **Why the pin fails in weak MV** (`probe_weak_mv_pin.R`): this is **not** a boundary. The `us(Trial):vm` fit used for selection **has not converged**, and keeps crawling towards Gₐ = 0 even after 400 more iterations (logL 236.79 → 237.90; Σnum/Σden 0.56 → 0.17). The statistic is computed mid-slope, so no null can be calibrated there. This regime only arises because the LRT gate was forced open: with the normal gate at 0.05, the LRT opened in just 13–18% of weak replicates.
- **"Partial" pinning** (pin only the eigen-directions where the observed W = S_target^−½ S_obs S_target^−½ equals 1) and "observed" pinning were both tested: **identical to the matrix version** on DH strong (300 replicates: 6.0% at 0.05, 10.7% at 0.10; 266/300 fits fully pinned). **Dropped.**

### Stage 1d: standalone MOST null grid — no LRT anywhere (`run_most_null*.R`, `summ_most_null.R`)

**Protocol (step 1):**
1. Fit the genome-wide term as **`diag(Trait):vm`** (UV: `vm`) via the engine internals (`.buildGenomeModel(str = "diag")`), continuing the fit if it hasn't converged.
2. Drop trials whose genomic variance is at the boundary (`B`, or ≈ 0). If all trials are bounded, **stop ("no-signal")**.
3. Test the maximum outlier statistic on the remaining trials against the **per-trial (`diag`) studentised** simulated null, which pins each trial's cross-product (the t constraints REML imposes for `diag`).

16 cells: {DH, GWAS} × {UV, 3-trial} × backgrounds {none (κ = λ = 0), IID-only (κ = 0, λ = 1), weak K (κ = 0.05), strong K (κ = 1)}. 500 replicates per UV cell, ~300 per MV cell; **6,389 usable replicates**. Gate stops count as non-rejections (full-procedure rate).

**MOST false-positive rate at α = 0.05:**

| Background | DH UV | DH 3-trial | GWAS UV | GWAS 3-trial |
|---|---|---|---|---|
| None | **1.4%** (285 no-signal) | **3.0%** (57) | **1.8%** (301) | **1.7%** (76) |
| IID only | **2.0%** | **4.4%** | **2.2%** | **5.4%** |
| Weak K | **4.2%** | **3.7%** | **3.4%** | **6.0%** (3.6–9.3) |
| Strong K | **5.2%** | **3.7%** | **5.6%** | **1.3%** (0.4–3.4) |

**Pooled: 216 / 6,389 = 3.4%.** No cell is significantly above 5%.

**Legacy LRT on the same fits:**
- no K-structured variance (its own null): 0.7–4.4% (valid, as Verbyla proved);
- weak K: 10–20%;
- strong K: 99–100%.

As a stopping rule, each of those rejections would start selecting a marker that isn't a QTL. **This is the core head-to-head result for the paper.**

Findings:
- **The no-signal case looks after itself:** 25–60% of no-signal replicates stop at the boundary; the rest are tested at about nominal.
- **Bounded trials are common in weak MV** (1–2 of 3), and dropping them keeps the rule calibrated.
- **Studentisation is essential:** the unstudentised null gives 0–3% and Bonferroni 0–3%, very conservative in DH.
- **The per-trial pin holds** (median deviation < 0.001).
- **The one weak spot is multi-trial with strong background on `diag`:** conservative (GWAS 1.3%, DH 3.7%). The *whole* p-value distribution is shifted up (mean 0.59), and the unstudentised null gives the same, so this looks like **model misspecification**: `diag` ignores the strong cross-trial background correlation (ρ = 0.6). It errs on the safe side; its cost in power is a Stage 2 question.
- **Harness fix (v2):** in the no-genetic-variance MV cells, `us(Trial):id` bases rarely converge. v2 falls back to `diag(Trial):id` on non-convergence as well as on error. This is a base-model choice, not part of MOST.

### Decision: golden rule (one model for UV and MV)

> **MOST:** the genome-wide search term is **`diag(Trait):vm`** (UV: `vm`, the t = 1 case). At each step:
> 1. drop trials whose genomic variance is at the boundary, and stop if none remain;
> 2. select the marker with the maximum outlier statistic **from the same `diag` statistic**;
> 3. test it against the per-trial studentised simulated null at α;
> 4. if significant, fit it as `X + Trait:X` and repeat; otherwise stop.
>
> The user's base structure (`us`, `fa(k)`, …) is kept for the residual genetic term and for final effect estimation and MAIN/INT pruning. No structure-specific studentisation is needed, convergence is reliable, and the boundary handling gives standalone "no-signal" stopping.

Costs to quantify in Stage 2:
- possible lost power in strongly correlated multi-trial backgrounds (the colleague's situation);
- MAIN vs G×E power, since `diag` weights trials independently instead of through Gₐ⁻.

## 6. Open issues and next steps

| Item | Status |
|---|---|
| **Studentised null** (UV scalar; MV matrix for `us`; per-trial for `diag`) | ✅ Stages 1b–1d |
| **Standalone MOST step 1** (no LRT; `diag` vm, bounded trials dropped, no-signal stop) | ✅ Stage 1d: 3.4% pooled at nominal 5% |
| **Golden rule:** MOST always on `diag(Trait):vm` (UV = t = 1) | ✅ Decided |
| **Stage 2: end-to-end MOST (golden rule) vs legacy LRT** with planted QTL: graded main effects, a fade, a crossover, a linked pair; weak/strong background; UV/3-trial; DH/GWAS; plus a colleague-like strong, correlated multi-trial case. Metrics: QTL recovered, false positives, number of steps, MAIN/INT accuracy, step-down error after the true QTL are fitted. Needs an end-to-end MOST loop built from the engine internals (`.mergeEffect`, `.rebuildCovObj`, `.addEffect`, `.packResults`), outside the package | ⬜ **Next** |
| Understand the conservatism of `diag` under a strongly correlated MV background (GWAS 1.3%, DH 3.7%) | ⬜ Investigate alongside Stage 2 |
| Check the eigenvalues of Vg for `fa(k)` / rank-deficient Gₐ (null mean 1.045) | ⬜ Moot under the golden rule (`diag` only) |
| **Stage 3:** current Oⱼ (proportional alternative) vs Sⱼ, especially for G×E QTL | ⬜ |
| **Stage 4:** parametric null on real genotypes (`genoRxK`, `genoSxT`, `genoMVsim`, GWAS panel, colleague's oat data) | ⬜ |
| Other cells: 6–8 trials `fa(2)`, `diag` with a boundary trial, unbalanced trials, ungenotyped lines in DH, spatial residuals | ⬜ |
| Run the colleague's data with the diagnostic attached (`source()` + wrapper) | ⬜ |
| Check the MV score derivation against the AOM literature; add the original AOM reference to the qmd | ⬜ |
| Add an "Inference for random effects" section to the qmd; add the Stage 0/1 results to the qmd | ⬜ |
| If adopted: implement MOST in the package as the new default-capable rule, with `stop.rule = c("most", "lrt")` (LRT = legacy). Includes the `diag` genome-wide term with bounded-trial dropping, per-trial studentised max-stat null (C++), B and seed exposed, a per-step MOST p-value/threshold trace (replacing `lik.mat` / `aimTrace(plot = "lrt")` for MOST), the chromosome-selection variant, docs, vignettes and tests | ⬜ Not started (package untouched) |
| Paper: MOST vs LRT; theory (AOM score test, studentised null as conditioning on the REML-pinned statistic); Stage 1d grid as the main calibration table | ⬜ |
| Side items found along the way: `summary()` `dist(cM)` label for physical maps; `pchisq.mixture` weights/citation/`df = 0`; LRT L1 refitted from default starting values every iteration | ⬜ Deferred by user |
