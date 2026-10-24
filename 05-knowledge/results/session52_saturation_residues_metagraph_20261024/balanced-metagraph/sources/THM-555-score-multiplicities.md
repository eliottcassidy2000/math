---
id: THM-555
title: Tiling-model cycle-moment hierarchy and the exact score->OCF boundary
status: MIXED — E[c3] & E[c_k] leading order PROVED; Var[c3], E[c5], E[c7], E[(-1)^c3], max-c3 census VERIFIED exact (adversarially re-derived); E[H] has no polynomial form
source: kind-pasteur-2026-06-20-S20 (application workflow on THM-554)
depends_on:
  - THM-554   # the score partition function Z_n (the engine)
  - THM-553   # two-clock tile address
  - THM-549   # complement-invariance / 2x half-tiling fold
related:
  - THM-027   # BIBD/regular c3-extreme
  - THM-552   # c3-parity dichotomy
  - HYP-2690  # codex address-DP program
---

# THM-555 - Tiling Cycle-Moment Hierarchy

Application of THM-554's score partition function `Z_n` (and the per-subset-linearity template
it exposes). Every value below was adversarially re-derived from scratch against brute
tournament-rebuild enumeration; one workflow self-script had a probability bug, caught and
corrected by brute. Distinguish PROVED (paper proof) / VERIFIED (exact, computationally proven,
confirmed beyond any fit window) / CONJECTURE.

## The cycle-count moments (closed forms)

```text
E[c3]      = (C(n,3) + (n-2)) / 4 = (n^3 - 3n^2 + 8n - 12)/24            [PROVED]
Var[c3]    = (n^3 - 7n^2 + 20n - 16) / 32                                [VERIFIED]
E[c5]      = (n^5 - 10n^4 + 45n^3 - 140n^2 + 294n - 280) / 160           [VERIFIED]
E[c7]      = (deg-7 poly) / 896                                          [VERIFIED]
E[c_k]     leading term = (k-1)!/2^k * C(n,k) ~ n^k / (k * 2^k)          [PROVED, leading order]
E[alpha_2] = 0,0,0,15/16,93/16,173/8,495/8,2395/16  (n=3..10)           [VERIFIED to n=8 full]
```

All by **per-subset linearity** (the template named in CLAUDE.md): an invariant that is a sum
over `k`-subsets of a local indicator has `E = sum_subsets P(local event)`, and the base path
fixes arcs only inside subsets with consecutive integers, so each `P` is a fair-coin computation
over the subset's free arcs. `E[c5]`/`E[c7]` were confirmed at `n=8`/`n=8` via per-subset OUTSIDE
the polynomial-fit window — not fit artifacts. (Full general-`k` lower-order correction series:
CONJECTURE — only the leading order is proved.)

## The c3 distribution structure

```text
min c3 = 0, multiplicity EXACTLY 1 (the unique all-bit-0 transitive tiling)   [PROVED]
max c3 = (n^3 - n)/24  (n odd),  (n^3 - 4n)/24  (n even)                       [VERIFIED n<=10]
max-c3 multiplicity (odd n)  = REGULAR census = 1, 3, 91, 29157 (n=3,5,7,9)    [VERIFIED]
max-c3 multiplicity (even n) = near-regular     = 5, 157, 51949 (n=4,6,8)      [VERIFIED]   (n=10: 204019829)
E[(-1)^c3] = 1 / 2^{floor((n-1)/2)}   (0 at n=3) => Pr[c3 odd] = 1/2 - 1/2^{floor((n-1)/2)+1}  [VERIFIED n<=10]
```

The distribution is a **parity comb** (non-unimodal for `n>=6`), always biased toward EVEN `c3`,
the bias halving every two steps.

## The exact score->OCF boundary

`Z_n` gives no-enumeration access to anything that is a function of the **score vector**, and the
exact **MEAN** of any sum-over-subsets-of-a-local-indicator (all `E[c_k]`, `E[alpha_2]`). The wall:

- **`c3 = C(n,3) - sum_v C(s_v,2)` is the LAST score-determined OCF datum.** Verified score-determined
  at `n<=7`. `c5` and `alpha_2` are NOT score-determined (counterexamples: `n=6` score `(0,2,3,3,3,4)`
  has `c5 in {1,2,3}`; score `(1,1,2,3,4,4)` has `alpha_2 in {0,1}`). Adding `A^2` row-sums does not
  rescue the distribution — fails already at `n=6` (empirical, not a proved impossibility).
- **Means still close** (linearity needs only per-subset probabilities); **distributions beyond `c3`
  need the full `2^F` state.** This is THM-442's "H not cell-affine" as a partition-function statement:
  the linear functionals of `Z_n` are exactly the cut-space (score) observables; `H`/OCF lives in the
  cycle space.
- `E[H]` (OCF, exact via Redei Ham-path count): `2, 4, 79/8, 29, 3175/32` (`n=3..7`). No polynomial
  form (the `4*alpha_2 + ...` disjoint-cycle-packing terms). The `n=5` identity
  `1+2E[c3]+2E[c5]=E[H]=79/8` is a low-`n` coincidence (diverges `n>=6`).

## Repo connections (verified)

- **THM-027 (regular = c3-extreme):** the tiling c3-MAXIMIZER is exactly the regular score, unique at
  odd `n`, multiplicity = the regular census; `max c3` = the Paley/regular 3-cycle count. The ensemble
  recovers regular=extremal from scratch. At `n=7`: mean `10`, `sigma~1.97`, Paley `c3=14 = mean+2sigma`
  — quantifies how extremal the H-maximizer is.
- **THM-552 (c3-parity dichotomy):** the ensemble even-bias `E[(-1)^c3]=1/2^{floor((n-1)/2)}>0` is the
  statistical SHADOW of the pointwise dichotomy.
- **THM-549/THM-280 (complement-halving):** `c3,c5,alpha_2` are complement-invariant (0 violations `n<=6`);
  the `tau<=n` half-tiling address quotient computes every moment at 2x. SC stratum = grid-symmetric =
  `2^{half}` cube exactly. `E_SC[c3] = (C(n,3)+(n-2))/4 + [n odd](n-3)/8` [CONJECTURE, half-cube enum `n<=10`].
- **H-maximizer:** the global H-max IS self-complementary for `n=3..8` (`H = 3,5,15,45,189,661`) — supports
  the half-cube maximizer search (HYP-2687/2688). (`n=8` rests on prior brute; general `n` unproven.)
- Iso-class `V_merged = 2,3,10,34` (`n=3..6`), `A000568=2,4,12,56`, `SC=2,2,8,12` — exact recompute match.

## Open / next

`Var[c3]` and the `E[c_k]` are VERIFIED-exact but the covariance / lower-order correction COUNTS are
confirmed numerically, not written as symbolic proofs (a clean per-run-length gap-signature count would
upgrade them to paper-proof). `Var[c5]`, `E[alpha_2]` full polynomial, and the `E_SC[c3]` half-cube
linearity are open. The regular census `1,3,91,29157` is worth an OEIS submission (its OEIS absence is
NOT independently verified). Scripts: `04-computation/tiling_gf_app_*_kps-Sx-wf.py`.
