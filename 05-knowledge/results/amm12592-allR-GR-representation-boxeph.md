# AMM 12592 — Angle B3: the direct G_R representation (y-box reduction, Lucas diagonal, capacity ledger, multiscale skeleton)

Session: boxeph 2026-08-03 (multifront worktree, post-THM-3329 state).
Conventions: `p = x`, `q = 1-x`, Bernstein cells `B_{d,k} = x^{d-k} q^k`,
epoch problem (*) `sum_{i<R} x^i Delta_i = q^{R-1}`, floor profile
`d_i = floor(gamma*(R+i)) + D0`, `gamma* = log_5(phi^2)`.
Master equation (S2): `Delta_i = (p-q) + c_i` (i <= R-2), `Delta_{R-1} = -1`,
`sum_{i<=R-2} x^i c_i = q^{R-1} - E_{R-1} =: 2 G_R`, `E_m = -1+x+...+x^m`.

Scripts (all exact int arithmetic, no numpy; outputs beside them in
`05-knowledge/results/`):

- `04-computation/amm12592_allR_GR_ybox_reduction_boxeph.py`  (steps A–F)
- `04-computation/amm12592_allR_GR_lucas_diagonal_boxeph.py`  (identities 1–7)
- `04-computation/amm12592_allR_GR_capacity_ledger_boxeph.py` (ledger, R<=256)
- `04-computation/amm12592_allR_GR_witness_anatomy_boxeph.py` (cancellation law)
- `04-computation/amm12592_allR_GR_dyadic_multiscale_boxeph.py` (IDENTITY M)
- `04-computation/amm12592_allR_GR_attractor_bridge_boxeph.py` (sigma/T bridge)

Status of every claim below is one of: PROVED (short proof included/indicated,
machine-verified), VERIFIED (exact machine check on stated range), OBSERVED
(measurement, no proof). No search negative below is used as infeasibility
evidence; the one impossibility statement (T4) is an exact capacity
computation for a precisely delimited scheme class, not a search result.

---

## T1 (PROVED + VERIFIED). Y-box reduction: parity eliminates entirely

S3 recompute (VERIFIED R = 2..1024): `[x^0] 2G_R = 2`;
`[x^j] 2G_R = (-1)^j C(R-1,j) - 1` for `j >= 1`; `2G_R` even coefficientwise
iff `R = 2^t`. (The `j = 0` cell is special: `+1 - (-1) = 2`.)

**Lemma (y-box).** Write `c_i = 2 gamma_i` and let `y_k` be the (unique)
degree-`d` Bernstein cells of `gamma_i`. Then

```
Delta_i = (p-q) + 2 gamma_i  admissible at d
   <=>   -C(d-1,k) <= y_k <= C(d-1,k-1)   for all 0 <= k <= d.
```

*Proof.* Ballot cells `b_k = C(d-1,k) - C(d-1,k-1)` satisfy
`b_k == C(d,k) mod 2` (Pascal), so parity of `b + 2y` is automatic, and
`|b_k + 2y_k| <= C(d,k) = C(d-1,k) + C(d-1,k-1)` unwraps to exactly the
asymmetric box above. QED. (Exhaustive check d <= 5 — 18204 blocks each
direction — plus 4000 random in-box/out-of-box pairs at d <= 40: all pass.)

Consequences (all VERIFIED):
- Box corners are `gamma = -p` (i.e. `Delta = -1`) and `gamma = +q`
  (`Delta = +1`).
- Boundary cells: `y_0 = gamma_i(1) in {-1,0}` (the ballot cut S4 in y-form),
  `y_d = gamma_i(0) in {0,1}` (mirror constraint at x = 0).
- Equivalent form: `gamma_i = q*beta_i - p*alpha_i` with `alpha_i, beta_i`
  subunital at degree `d_i - 1` (cells in `[0, C(d_i-1,k)]`).
- Atom table: `+p^a q^b` alone is in-box iff `b >= 1`; `-p^a q^b` iff
  `a >= 1`; cells of `p^a q^b` at degree d are `C(d-a-b, k-b)`.

**The all-R problem is now parity-free:** find integer cell arrays
`y_i` in the asymmetric boxes with `sum_{i<=R-2} x^i gamma_i = G_R`.
In particular *cross-row carry moves* `(gamma_i, gamma_{i+1}) ->
(gamma_i - x h, gamma_{i+1} + h)` are unconstrained by parity in y-space —
the redistribution program of THM-3026 sec 5 loses its parity lock here.

Witness map (VERIFIED): R = 8, 16, 32 (combined file), R8/R16 ruleB,
R128 ruleB D0=0, R256 ruleA D0=1, R512 ruleA D0=8 — every witness is in
ballot normal form, maps into the y-boxes, and `sum x^i gamma_i = G_R`
exactly.

## T2 (PROVED + VERIFIED). The Lucas diagonal (staircase) expansion of G_R

With `e = pq` and `L_m = p^m + q^m` (Lucas polynomial in e since `p+q = 1`):

- **Step law** (VERIFIED R = 2..399): `2(G_{R+1} - G_R) = -p L_{R-1}`.
- **Doubling law** (VERIFIED sample R <= 256):
  `G_{2R} = 2q G_R^2 + (p-q)(2 G_R - 1) - p^R q^{R-1}`.
  [From `2q G_R = L_R - (p-q)`, `L_2R = L_R^2 - 2e^R`, `(p-q)^2 = 1-4e`.]
- **Staircase expansion** (VERIFIED R <= 1024): since
  `L_R = 1 + sum_{k=1}^{floor(R/2)} a_k e^k` with
  `a_k = (-1)^k [C(R-k,k) + C(R-k-1,k-1)]`, and `(1-(p-q))/q = 2`,
  `e^k/q = p^k q^{k-1}`:

  ```
  2 G_R = 2 + sum_{k=1}^{floor(R/2)} a_k p^k q^{k-1}.
  ```

  The atoms `p^k q^{k-1} = x^k + higher` are triangular, so these
  coefficients are UNIQUE in this basis.
- **Dyadic evenness** (VERIFIED R = 2..1024): all `a_k` (k >= 1) even iff
  `R = 2^t` (equivalently `L_R == 1 mod 2` iff R dyadic, via
  `L_{2m} == L_m^2 mod 2`, `L_1 = 1`, `L_3 = 1 + e`). So at dyadic R,
  `G_R = 1 + sum_k m_k p^k q^{k-1}` with `m_k = a_k/2 in Z`, `m_1 = -R/2`,
  sign `(-1)^k`.
- **x = 1 ledger** (VERIFIED): `G_R(1) = -(R-2)/2`, carried entirely by the
  head `1 + m_1 p` (all k >= 2 atoms vanish at x = 1). Exactly `R/2 - 1`
  rows must have `gamma_i(1) = -1` (marker rows), the rest 0.
- **Peak profile** (VERIFIED table): `argmax_k |m_k| / R -> 0.276...
  = (1 - 1/sqrt5)/2` (Fibonacci-binomial peak), `log2 max_k |m_k| -> R log2
  phi ~ 0.694 R` — the golden ratio enters the CONSTRUCTION side through
  the staircase coefficients themselves.

Closed forms: `G_2 = q`; `G_4 = q^2 - p^3 = 1 - 2p + p^2 q`, giving the first
closed-form epoch witness beyond R = 2 (VERIFIED):
`Delta = ((p-q) + 2q^2, (p-q) - 2p^2, (p-q), -1)`, i.e.
`gamma = (q^2, -p^2, 0)` — at the gamma* floor profile for R = 4.

## T3 (PROVED + VERIFIED). Attractor bridge

With tail targets `T_0 = G_R`, `x T_{i+1} = T_i - gamma_i`, and Delta-space
residuals `sigma_i`:

```
sigma_i = E_{R-2-i} + 2 T_{i+1}      (VERIFIED rowwise on R = 16, 128 witnesses)
```

[Proof: master equation plus `(p-q)(1+x+...+x^i) = E_i + 2x^{i+1}`.]
So Rule A's attractor (`sigma_i` hugging the `E`-shape) is *identically* the
statement that the y-space remaining target is small. The dynamic and the
direct pictures are the same object; B3's "bypass dynamics" = "make T_i
explicit".

## T4 (PROVED, exact capacity computation). Sign-coherent atom-local routing of the staircase expansion is impossible — the R/3 law

Ledger (`amm12592_allR_GR_capacity_ledger_boxeph.out`, R = 8..256, D0 = 0):
for atom `m_k p^k q^{k-1}` placed in row `i <= k` (content
`p^{k-i} q^{k-1}`, cells `mu * C(M, j-(k-1))`, `M = d_i - (2k-1-i)`), the
exact per-row coefficient capacity is

```
cap+(i,k) = min_j floor(C(d_i-1, j-1) / C(M, j-(k-1)))    (mu > 0)
cap-(i,k) = min_j floor(C(d_i-1, j)   / C(M, j-(k-1)))    (mu < 0)
```

**Statement.** In any scheme that splits each staircase atom over rows with
sign-coherent, cellwise non-cancelling contributions (at every cell all
incoming terms share the sign of the total), row shares satisfy
`|mu_{k,i}| <= cap(i,k)`, hence such a scheme requires
`sum_i cap(i,k) >= |m_k|` for every k. The ledger computes both sides
exactly: the condition FAILS for all `2 <= k < k*(R)` with
`k*(R)/R = 0.375, 0.375, 0.359, 0.359` at `R = 16, 32, 64/128, 256` —
converging on `1/3` from above.

**Why 1/3 exactly (entropy level):** private capacity is governed by the
mid-support pinch `~ 2^{a+b} = 2^{2k-1-i} <= 4^k`, while the demand is
`|m_k| ~ C(R-k,k)`; and `C(2R/3, R/3) = 2^{(2R/3) H(1/2)} = 4^{R/3}` — the
transition is the exact solution of `(1-kappa)H(kappa/(1-kappa)) = 2 kappa`,
i.e. `kappa = 1/3`. Maximal deficit `~ 0.26 R` bits at `k ~ R/8` (measured
`2^-66` at R = 256, k = 32).

**Corollary.** `D0 = o(R)` slack cannot repair this (capacity gains
`~ 2^{D0}` against a deficit `~ 2^{0.26 R}`). Any all-R construction must
route mass through *cancelling* contributions — cross-row carries or
opposite-sign atom overlaps (the alternating signs `(-1)^k` of `m_k` make
adjacent-atom cancellation naturally available). Same verdict for the head:
the correction demands `-N_u p^u q` (`N_u ~ R/2`) at `u <= ~log2 R` have
private capacity `<= 4^u` — "NEEDS CARRIES" rows in the ledger.

This *kills route (a) in its naive per-atom form* and delimits exactly what
routes (b)/(c) must supply. (It is consistent with the transportation-LP
feasibility at 21% saturation quoted in the angle brief: that LP has all
cells free, i.e. cancellation enabled.)

## T5 (OBSERVED, exact measurements). How winning witnesses actually route mass

`amm12592_allR_GR_witness_anatomy_boxeph.out`, per verified witness:

- Marker words: exactly `R/2 - 1` rows with `gamma_i(1) = -1`, all values in
  `{0,-1}` (S4 confirmed in y-form). The rule-generated words are strikingly
  periodic early: `0--00 0--00 ...` (density 2/5), then a solid `-...-` band,
  then all-0 tail; R = 512 shows `0--00` * many, then `00---`-type phase
  shifts. (Note 2/5 + solid-band structure; unexplained, worth a follow-up
  against the gamma* ~ 3/5 rate.)
- Cell sparsity: nonzero cells 9.2% (R=128), 7.9% (R=256), 6.7% (R=512);
  max cell saturation = 1 exactly (some cell always rides the box).
- **Cancellation law:** per coefficient t, gross cross-row mass
  `S_t = sum_i |[x^{t-i}] gamma_i|` vs net `|g_t|`:
  `S_t = |g_t|` (bitdiff 0) for `t <~ 7R/8`; deep cancellation is confined
  to the top coefficient band, max `log2(S_t/|g_t|)` at `t = R-1`:
  16 / 36 / 65 / 109 bits at R = 32 / 128 / 256 / 512 (~ 0.21 R,
  slowly declining ratio). So winning constructions are per-coefficient
  LOCAL through the bulk and spend all their cancellation budget on the top
  band, where `|g_t|` is small but many rows' high cells overlap.

## T6 (PROVED + VERIFIED). Dyadic multiscale skeleton (IDENTITY M)

Substitute `u = -x` (`q = 1+u`), `R = 2^t`. Telescoping
`prod_j (1+u)^{2^j} - prod_j (1+u^{2^j})` (the second product =
`1 + u + ... + u^{R-1}`, binary digits) with
`(1+u)^{2^l} - 1 - u^{2^l} = 2 W_l`, `W_l = sum_{0<n<2^l} (C(2^l,n)/2) u^n`
(integer — dyadic Pascal interior is even: the parity miracle in product
form) gives, with `Ghat_R(u) = G_R(-u)` (coefficients
`(C(R-1,n) - (-1)^n)/2 >= 0` for n >= 1, and 1 at n = 0):

```
Ghat_R = 1 + sum_{j odd < R} u^j
       + sum_{l=0}^{t-1} (1+u)^{2^l - 1} W_l * sum_{m < 2^{t-l-1}} u^{m 2^{l+1}}
```

VERIFIED exactly R = 4..1024. Properties: every coefficient receives at most
`log2 R` pieces (one per scale, measured max overlap 8 at R = 512); all
pieces have nonnegative coefficients in u; scale-l piece
`P_l = (1+u)^{2^l-1} W_l` has support width `2^{l+1} - 2` and sup-bits
`~ 2^{l+1}` (1, 5, 12, 28, 59, 123, 250, 506 at R = 512), planted at every
`2^{l+1}`-aligned offset. The top scale is `(1+u)^{R/2-1} W_{t-1}` — an
epoch-R/2 binomial times the half-Pascal row: IDENTITY M is a doubling
recursion skeleton (route (c)) with explicit, positive, O(log R)-overlap
increments.

## Verdict on the three B3 routes

- (a) binomial/per-atom distribution: **dead as stated** (T4, exact); must be
  upgraded with explicit cancellation. The alternating staircase signs and
  the y-space carry moves (parity-free, T1) are the designated mechanisms.
- (b) kernel route: the staircase kernels `K_k = m_k p^k q^{k-1}` exist and
  are unique/explicit (T2) but not row-feasible per-atom (T4).
- (c) recursive route: **most promising**. Two explicit recursions are now
  on the table: the G-doubling law (T2) and IDENTITY M (T6), whose
  increments are positive, dyadically aligned, and O(log R)-overlapping.
  The missing piece is a row-feasible routing of the scale-l pieces; the
  top-scale piece still spans the whole epoch, so the induction must
  consume the R/2-representation (as in the carve-and-carry map) rather
  than re-derive it.

## Next obligations (concrete)

1. Prove or refute: the scale-l pieces of IDENTITY M with `l <= t-2` can be
   routed into rows `[m 2^{l+1}, (m+1) 2^{l+1})` within a FIXED fraction of
   the y-boxes (their sup-bits `~2^{l+1}` sit far under local cap bits
   `~gamma*(R + m 2^{l+1})` for small l; the constraint is per-cell shape,
   not mass — needs the T4-style exact ledger per scale).
2. Explain the `0--00`-periodic marker words (T5) — a 2/5-density opening
   law would pin the head routing for an explicit scheme.
3. Recast the carve-and-carry doubling map in y-space (parity-free carries)
   and re-run the 64 -> 128 failure cases; the parity lock was a plausible
   cause of steering fallbacks (THM-3302 sec 2 lists them as transport
   artifacts).
4. Top-band cancellation budget: prove the `~0.2 R`-bit growth law (T5) or
   bound it; it is the quantitative obstruction any explicit scheme must
   fund in rows `i > 7R/8`.

## Hazards honored

- No search negative is treated as infeasibility; T4 is an exact capacity
  argument whose scheme class is stated precisely (sign-coherent cellwise
  non-cancelling splits of the staircase atoms).
- All identities were recomputed from scratch, including the (S3) j = 0
  special case (`[x^0] 2G_R = 2`); an initial target-array bug at n = 0 in
  the IDENTITY M check was caught by the verifier and fixed (the identity
  itself was never wrong).
- Effective-rate/slack trap: T4's corollary quantifies exactly why o(R)
  slack cannot substitute for cancellation.
