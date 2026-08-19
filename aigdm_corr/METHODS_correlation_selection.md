# AIGDM with selectable GEE working-correlation structure and GIC-based selection

## Summary

The two original AIGDM implementations — `AIGDM.R` (exchangeable / compound-symmetric
working correlation) and `AIGDM_Toep.R` (Toeplitz working correlation) — have been merged
into a **single** `AIGDM.Cluster()` that accepts a `corstr` argument and supports four GEE
working-correlation structures. When `corstr = "GIC"`, the correlation structure is selected
**per taxon** (and, in the tree/lineage mode, per taxon within each lineage) by a
QIC/CIC-style generalized information criterion.

The merge is behavior-preserving: with `corstr = "exchangeable"` the merged code reproduces
the original `AIGDM.R` score-test p-values, and with `corstr = "toeplitz"` it reproduces the
original `AIGDM_Toep.R` p-values, to machine precision (bit-identical on the validation data).

## Model context

AIGDM tests covariate association with microbiome relative abundances for clustered /
longitudinal designs using an adaptively-inflated (zero-inflated) generalized Dirichlet
multinomial parameterized as a two-part (zero-inflation + Beta) model, fit by an EM /
estimating-equations scheme. The GEE working correlation enters the estimating equations only
through the working covariance of the linearized responses (the mean-model `V_{Z*}` and the
zero-model `V_{Δ}`); under standard GEE theory the point estimates remain consistent even when
the working correlation is misspecified, so the working correlation affects the estimator's
**variance / efficiency**, not the estimand. This is what makes an information criterion the
appropriate tool for choosing among structures.

## Supported working-correlation structures

For subject *i* with observed ordinal study periods `p_i = (p_{i1}, ..., p_{in_i})`, the
working correlation matrix `R_i` is built entry-wise as `R_i[s,t] = f(|p_{is} - p_{it}|)`:

| `corstr`        | correlation params | `R_i[s,t]`, s≠t                                   |
|-----------------|--------------------|---------------------------------------------------|
| `independence`  | 0                  | 0 (identity `R_i`)                                |
| `exchangeable`  | 1 (ρ)              | ρ for all off-diagonal pairs                      |
| `ar1`           | 1 (ρ)              | ρ^{|p_{is} − p_{it}|}                             |
| `toeplitz`      | L (banded)         | ρ_h for lag h = |p_{is} − p_{it}|, h = 1..L; 0 if h>L |

Diagonal entries are 1. The exchangeable and Toeplitz branches use exactly the correlation
constructions and moment estimators of the two original files; independence and AR-1 are new.

### Lag convention (ordinal period)
The lag between two samples is the absolute difference of their **ordinal study-period index**
(`|period_s − period_t|`), *not* real elapsed time. This matches the original `AIGDM_Toep.R`
exactly. Two samples collected in the **same** study period have lag 0.

### Caveat: within-period replication
In the ECAM genus data, most (subject, period) cells contain **more than one sample** (up to
~14; ~92% of cells have ≥2). Under the ordinal-lag convention such same-period replicate pairs
have lag 0 and therefore receive a working correlation of **0** under Toeplitz (they do not
contribute to any estimated lag-h band) and **1** under AR-1 (ρ^0). This is the intended,
documented behavior — it reproduces the original Toeplitz code's treatment of replicates — but
it means the Toeplitz/AR-1 lag bands are estimated only from **across-period** sample pairs.
If one instead wanted same-period replicates to be positively correlated, that would require a
distinct (e.g. nested / block-exchangeable-within-period) structure, which is out of scope
here.

### Data-driven maximum lag (Toeplitz)
The number of Toeplitz bands L is set data-drivenly to the **maximum within-subject ordinal
period distance** actually observed, `L = max_i (max(p_i) − min(p_i))`, after sorting by
(subject, period). On the ECAM genus data (periods 1..5) this gives L = 4. This generalizes to
other designs automatically and avoids estimating empty high-lag bands.

## GIC (QIC/CIC-style) selection

For each taxon and each candidate structure R the criterion is

```
GIC(R) = -2 * Q̂ + 2 * trace( Ω̂_I^{-1} · V̂_R )
```

following Pan's QIC for GEE [5] and the QIC/CIC formulation implemented by Cui [1]:

- **Q̂** — the two-part (zero-inflation Bernoulli + Beta mean/dispersion) **independence
  quasi-log-likelihood** evaluated at the fitted estimates. Because the working correlation
  does not enter this quasi-likelihood, Q̂ is (essentially) shared across the four structures
  for a given taxon; selection is therefore driven by the penalty, which is the intended QIC
  behavior.

- **Ω̂_I** — the model-based ("naive") covariance of the GEE estimator computed under the
  **independence** working correlation, `Ω̂_I = A_I^{-1}`, where `A_I = Σ_i D_i' V_i^{-1} D_i`
  is the GEE sensitivity (bread) with `V_i` built from the identity correlation.

- **V̂_R** — the **robust sandwich** covariance under the candidate structure R,
  `V̂_R = A_R^{-1} B_R A_R^{-1}`, with sensitivity `A_R = Σ_i D_i' V_i^{-1} W_i D_i` and
  variability (meat) `B_R = Σ_i U_i U_i'`, `U_i = D_i' V_i^{-1} W_i r_i`.

`trace(Ω̂_I^{-1} V̂_R)` is the QIC/CIC **generalized number of parameters** (effective degrees
of freedom). The structure minimizing GIC(R) is selected.

### GEE ingredients (mean model)
The penalty is computed on the mean-model GEE block, whose per-subject quantities match the
estimating equations solved inside the fitter:

- `D_i` = ∂μ*/∂α block = `.dma(Xa_i, α̂, μ̂, φ̂)` assembled block-diagonally,
- `V_i` = `diag(sd_{Z*}) · R_i(structure) · diag(sd_{Z*})`, `sd_{Z*}` from `.var_ZStar(μ̂, φ̂)`,
- `W_i` = `diag(1 − E[Δ|·])` (posterior non-zero weight),
- `r_i` = `E[Z*] − μ*` (mean-model working residual).

Because `A = Σ_i D_i' V_i^{-1} W_i D_i` is assembled as a **Gram (sum-of-quadratic-forms)
matrix**, it is positive semi-definite by construction. This guarantees a non-negative
effective-d.f. penalty (an earlier construction based on the posterior-weighted second-
derivative Hessian was not PSD and produced negative penalties; that has been replaced).

### Scope of selection
Selection is applied to the **mean-model** GEE block, which is the primary estimand and the
correlation target; the same selected structure is used for the zero and dispersion blocks in
the merged fit. Restricting the criterion to the mean model keeps the penalty well-defined and
comparable across structures. Penalizing the correlation parameters through the effective d.f.
(rather than the raw parameter count) is consistent with the QIC/CIC literature and tends to
improve selection of the true working correlation [1, 4, 5, 6].

## API

```r
AIGDM.Cluster(ID, OTU, X, X.index,
              period = NULL,                 # ordinal period index, one per sample row
              corstr = "GIC",                # independence | exchangeable | ar1 | toeplitz | GIC
              Tax = NULL, test.type = "Omni",
              min.depth = 0, ZI.pos = "adaptive",
              n.boot = 499, n.perm = NULL,
              n.cores = detectCores() - 1, fdr.alpha = 0.05)
```

- `period` is **required** when `corstr ∈ {ar1, toeplitz, GIC}` (an integer-valued vector, one
  value per sample row, matching `nrow(OTU)`); for `independence` / `exchangeable` it is
  optional and defaults to a single level (all lag structure collapses).
- The maximum Toeplitz lag is computed data-drivenly (see above).
- Returns, in addition to the usual p-value / significant-lineage output:
  - `corr.structure` — the selected structure per taxon (non-tree mode: a named vector; tree
    mode: a named list, per lineage), and
  - `gic.table` (only when `corstr = "GIC"`) — a data frame with, per taxon (and per lineage
    in tree mode), the columns `structure, gic, qhat, penalty, selected`.

## Validation

- **Equivalence.** On deterministic synthetic clustered longitudinal count data, the merged
  code with `corstr = "exchangeable"` reproduces the original `AIGDM.R` asymptotic score-test
  p-values, and with `corstr = "toeplitz"` reproduces the original `AIGDM_Toep.R` p-values,
  with `max|difference| = 0` (bit-identical).
- **New structures.** `independence` yields the identity working correlation; `ar1` yields
  `ρ^{|Δperiod|}` with proper decay.
- **GIC.** Per-taxon selection returns exactly one structure per taxon (= argmin GIC), fits
  all four candidates per taxon, and produces **positive** effective-d.f. penalties that order
  monotonically with structure complexity (independence ≤ exchangeable < ar1/toeplitz).
- **API guards.** Missing / wrong-length / non-integer `period` (when required), and invalid
  `corstr`, raise informative errors; `independence` and `exchangeable` run without `period`.

## References

- [1] Cui J. QIC Program and Model Selection in GEE Analyses. *Stata Journal* 2007.
  doi:10.1177/1536867X0700700205
- [4] Gosho M, Hamada C, Yoshimura I. Modifications of QIC and CIC for Selecting a Working
  Correlation Structure in the GEE. *Japanese Journal of Biometrics* 2011. doi:10.5691/jjb.32.1
- [5] Pan W. Akaike's Information Criterion in Generalized Estimating Equations.
  *Biometrics* 2001. doi:10.1111/j.0006-341X.2001.00120.x
- [6] Nyabwanga RN et al. Consistency Inference Property of QIC in Selecting the True Working
  Correlation Structure. *American Journal of Theoretical and Applied Statistics* 2019.
  doi:10.11648/j.ajtas.20190802.14
