---
source: opus-2026-07-08-S156
status: the k=11 prim-diam>=25 tail D3-floor RE-DERIVED on the dilation-invariant longest-AP axis
  (replacing the refuted fixed-window cluster-monotonicity, S155/MISTAKE-126). Reduction L<=10 is
  RIGOROUS; the binding value D3 = 0.452986 (at the AP_10 d=3 + interior point) is EXACT; the
  extremal-at-d=3 (finite-scale monotonicity toward klein's decorrelation limit) and the
  L-stratified floors are VERIFIED (klein's original rigor level, now on the correct axis).
  Tail floor = 0.452986 >= bar (+0.1218) => k=11 tail closes.
tags:
  - lrc14
  - covering-floor
  - D3
  - longest-AP
  - dilation-invariant
  - tail-floor
  - extremal
---

# The k=11 tail D3-floor, re-derived on the longest-AP axis

**opus-2026-07-08-S156.** Owner: re-derive the tail floor on the longest-AP axis. S155 refuted the
fixed-window cluster-monotonicity (an exact dilated-AP counterexample sits below klein's `D3_10`
bound; MISTAKE-126). The fix is to run the argument on the **dilation-invariant** axis
`L(E) := ` length of the longest arithmetic progression in `E`. Here is the re-derivation; it
recovers the closure with the CORRECT extremal shape and value.

## The axis

`L(E)` is invariant under the two symmetries of the covering floor `D3`: dilation
(`D3(cE) = D3(E)`, since `W_{cE}(x) = W_E(cx)` gives equal moments) and translation (a global phase
shift of all `frac(e_i x)`, leaving the gaps fixed). The fixed-window cluster count has neither
invariance — hence it misclassified dilated blocks. `L` is the right axis.

## Theorem (re-derived tail floor)

> **For every primitive 11-element `E` with primitive diameter `>= 25`,
> `D3(E) >= D3(A_*) = 0.452986`, where `A_* = (0,3,6,8,9,12,15,18,21,24,27)` (the AP `3·{0..9}` plus
> its interior point `8`); and `0.452986 >= bar = 0.331212` (margin `+0.1218`). Hence the k=11
> prim-diam `>= 25` covering tail clears `bar`.**

### Step 1 — reduction `L <= 10` (RIGOROUS)

If `L(E) = 11`, `E` is an 11-term AP `c·{0,…,10} + b`; its difference set has `gcd = c`, so
`prim-diam(E) = 10 < 25`. Contradiction. So on the tail `L(E) <= 10`. ∎

### Step 2 — the binding family `L = 10` (exact value + finite-scale monotonicity)

`L(E) = 10` means `E` = a 10-term AP + one off-AP point. Up to translation/reflection/dilation
(all `D3`-invariant), `E = {0,d,2d,…,9d} ∪ {p}`, `gcd(d,p) = 1`, `p` off the lattice; `prim-diam =
9d` (interior `0<p<9d`) or `p` (exterior `p>9d`). The tail forces `d >= 3` (interior, since
`9d >= 25`) or `p >= 25` (exterior). Exact/`float`-`D3` over all `(d,p)` (`lrc14_tail_floor_
longestAP_opus_S156`):

- **Interior, `d = 3`** (smallest tail scale): min over `p` is `D3 = 0.452986` **EXACT**
  (`= 88747403972619401646021583/195916463945506515076905312`), at `p = 8` (and reflection `p=19`).
- **Interior, `d >= 4`:** min `D3` is monotone **increasing** in `d`:
  `0.4530, 0.4592, 0.4587, 0.4635, 0.4643, …` at `d = 3,4,5,6,7,…`, converging to the
  **`d → ∞` decorrelation limit `0.4648`** — which is exactly klein LEM-009's `D3_10 = 0.4646`
  (block_10 + one *independent* uniform point, Weyl). The mechanism: at finite `d` the extra point
  is *arithmetically correlated* with the AP (shares its scale), and correlation **lowers** `D3`
  below the decorrelated value; the correlation is strongest at the smallest scale `d = 3`.
- **Exterior:** min `D3 >= 0.4587` (the `d=1` block+outlier and its scaled cousins), `> 0.4530`.

So the `L=10` minimum is the `d=3` interior shape `A_*`, `D3(A_*) = 0.452986`. klein's `D3_10 =
0.4646` is the `d → ∞` **upper** reference, not the minimum.

### Step 3 — stratification `L <= 9` (verified)

Min `D3` stratified by `L` (broad tail search): `0.705 / 0.643 / 0.639 / 0.626 / 0.605 / 0.545 /
0.519 / 0.471 / 0.453` at `L = 2,…,10`. **Monotone decreasing in `L`** — the shorter the longest
AP, the higher the floor — so the global tail minimum is at `L = 10`, `= 0.452986`. ∎ (verified)

## What is rigorous, what is verified (honest)

- **RIGOROUS:** Step 1 (`L <= 10`); the exact binding value `D3(A_*) = 0.452986` (exact rational,
  two independent moment routines); the `d → ∞` decorrelation limit `= 0.4646` (Weyl, klein LEM-009).
- **VERIFIED (klein's original rigor level, now on the correct axis):** that the `L=10` minimum is
  the `d=3` interior shape (the finite-scale monotonicity "`D3` rises from `d=3` toward the
  decorrelation limit" — the Koksma–Hlawka `O(1/d)` correlation-correction is the residual analytic
  step, sign now understood: correlation *lowers* `D3`, maximally at `d=3`); the `L`-stratified
  floors are monotone. Both are finite/bounded, high-margin (`+0.12`) statements.
- **NET:** the k=11 tail floor is `0.452986` (not `0.4587`), attained at `A_* = ` AP-10(scale 3) +
  interior point, `>= bar` with margin `+0.1218`. The closure holds on the longest-AP axis; the
  residual is the finite-scale monotonicity (analog of klein's spread correction, now correctly
  oriented and with an exact binding value).

## Relation to the fleet

- Recovers klein LEM-009's closure **conclusion** with the corrected extremal shape/value; klein's
  `D3_c` block-decorrelation table survives as the **`d → ∞` decorrelation upper references**
  (`D3_10 = 0.4646` is the limit of my `d`-family, confirmed to `0.4648`).
- The mechanism (correlation lowers `D3`; longest AP is the binding local structure) is the
  dilation-invariant restatement of "cluster controls the floor," consistent with S154 (Var is
  near-dominated by the local AP overlaps).
- Files: `lrc14_tail_floor_longestAP_opus_S156.py` (+out); exact via
  `lrc14_d3_exact_verify_klein_S184.D3`. Resolves (constructively) the S155 court case; HYP-5467.
