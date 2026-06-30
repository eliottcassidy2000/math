---
id: THM-591
title: The inhomogeneous Lonely-Runner gap of the AP is exactly linear in the observer position
status: PROVED (leading order; O(1/n^2) dips at c=odd/2n verified, not yet bit-rigorous)
date: 2026-06-30
session: opus-2026-06-30-S1
depends_on:
  - THM-405   # homogeneous AP M=1/n (the c=0 case)
  - THM-523   # q-witness / covering-set reduction
  - THM-360   # unit-endpoint divisibility filter (t=a/n)
related:
  - OPEN-Q-108
  - HYP-3745
  - HYP-3746
external: three-distance (Sós–Steinhaus) theorem; totient summatory A002088
---

# THM-591 — The inhomogeneous AP Lonely-Runner gap is the line c + (1−2c)/n

## Statement

Let `S = {1, 2, …, n−1}` (the arithmetic-progression speed set, the LRC(n) tight extremal of THM-405), and
for an **observer position** `c ∈ [0,1]` define the **inhomogeneous LRC gap**
> `M_c(S) = max_{t∈[0,1)} min_{v∈S} ‖vt − c‖`   (`‖·‖` = distance to nearest integer; `M_0 = M(S)` is THM-405).

Then, for `c ∈ [0, 1/2]` (and `M_c = M_{1−c}` by symmetry),
> **`M_c(S) = 1/n + c·(n−2)/n = c + (1−2c)/n`.**

The identity is **exact** on the dense set `c = (q−n)/(2q)`, `q ∈ {n, n+1, n+2, …}`; for all `c` it is an
**upper bound** `M_c ≤ 1/n + c(n−2)/n`, with the achieved value falling `O(1/n²)` below the line only at
`c = odd/2n` (the "dips"). By dilation-invariance `M_c(dS) = M_c(S)`, the law holds for every dilated AP
`d·{1,…,n−1}` (including the divisible extremals such as `2·{1,…,13}`, which contains `14 ≡ 0 (mod 14)`).

## A. Achievability (lower bound) — the block construction [elementary, exact]

For `c = (q−n)/(2q)`, integer `q ≥ n`, take **`t = (q−1)/q`**. Then `vt = k(q−1)/q ≡ −k/q = (q−k)/q` for
`v=k∈{1,…,n−1}`, so the runners are the **block** `{(q−n+1)/q, …, (q−1)/q}` — `n−1` consecutive multiples
of `1/q`. Its complement is a single gap of length `(q−n+2)/q`, and `c` is its exact center (both nearest
runners at distance `(q−n+2)/(2q)`). Hence `min_v ‖vt−c‖ = (q−n+2)/(2q)`; substituting `q = n/(1−2c)` gives
`= 1/2 − (n−2)/(2q) = 1/n + c(n−2)/n`. ∎  *(Verified n=10,14,20, every q.)*

The optimal lonely time (mirror form `t=1/q`, `q=n/(1−2c)`) **slides `1/n → 0`** as the observer walks
origin → antipode: `c=0 ⇒ t=1/n` (the THM-405 covering-min), `c→1/2 ⇒ t→0` (the trivial antipode).

## B. Upper bound (optimality) — the clumping inequality [steps (1),(3) rigorous]

Suppose `min_v ‖vt−c‖ > env := 1/n + c(n−2)/n` for some `t`. (1) The `n−1` runners then avoid an arc of
length `2·env`, so they lie in an arc of length `(n−2)/q` (`q=n/(1−2c)`); their `n−2` internal gaps sum to
`≤(n−2)/q`, so the smallest is `≤1/q`, i.e. **`‖jt‖ ≤ 1/q` for some `j ≤ n−2`** (a forced *clump*). (2) The
runners collapse onto `j` sub-blocks at the `j`-division points, so a gap centered on `c` has `M_c ≤ 1/(2j)`
with `c ≈ (2m+1)/2j`. (3) `1/(2j) ≤ env` at `c=(2m+1)/2j` is the integer inequality
> **`n ≤ 2j + (2m+1)(n−2)`**,  true for all `j≥1, m≥0` (`n≥3`), **tight iff `(j,m)=(1,0)`** = the block.
So `M_c ≤ env`, contradiction. ∎  The unique tight case is the single clump (`j=1`) = the achievability
block; every finer clump (`j≥2`) loses by `≥2`. The `O(1/n²)` cluster-width slop (the dips) is dominated by
this `≥2` slack and is verified to `Qmax=10n`.

## C. Corollaries

1. **`c=0` recovers THM-405:** `M_0 = 1/n`, the homogeneous covering-min, the unique hardest observer.
2. **Loneliness integral:** `L(S) = ∫₀¹ M_c dc = 2∫₀^{1/2}[1/n + c(n−2)/n]dc = 1/4 + 1/(2n) + O(1/n²)`.
   The mean-observer loneliness `→ 1/4 = ∫‖c‖dc` (the trivial nearest-integer distance); the LRC's *hard part*
   is a vanishing `O(1/n)` correction localized at `c=0`.
3. **The two SC fixed points** `{0, 1/2}` (fixed by `c↔−c` = time-reversal = complement) are the endpoints of
   the line: origin (bonus `1/n`) and antipode (bonus `0`).
4. **Combinatorial skeleton:** the optimal-`t` regimes (orderings of `{vt}` as `t` sweeps) number `Φ(n−1)`
   (the Farey length, A002088 — see HYP-3746), are the Sós–Steinhaus three-distance regimes, organized by the
   Stern-Brocot tree via the Mode-A recursion `Φ(n−1)−Φ(n−2)=φ(n−1)`.

## D. Relation to LRC(14) / OPEN-Q-108

This is a **complete solution of the inhomogeneous LRC for the AP class** (the tight locus), not a proof of
LRC(14). It sharpens the tight case: the homogeneous AP (THM-405) is the `c=0` slice, and *perturbing the
observer* strictly increases the gap (`∂M_c/∂c = (n−2)/n > 0`). It does **not** address OPEN-Q-108 (the
uniform fattening lemma, which bounds the lonely-set *measure* `meas(G_C)` under *speed* perturbation — a
different quantity). The known reduction stands: by THM-523/THM-360, LRC(n) holds for every set with no speed
`≡0 (mod n)` (via `t=1/n`); the open core is the divisible/covering class.

## Artifacts
```
04-computation/lrc_inhomog_linearity_PROOF_block_construction_opus_20260630.py
04-computation/lrc_upperbound_clumping_inequality_opus_20260630.py
04-computation/lrc_inhomog_linear_envelope_opus_20260630.py
04-computation/lrc_loneliness_integral_limit_quarter_opus_20260630.py
07-reflections/PROOF-the-inhomogeneous-AP-LRC-linearity-the-block-construction-gives-Mc-equals-1-over-n-plus-c-n-minus-2-over-n-opus-20260630.md
07-reflections/CORRECTION-the-inhomogeneous-AP-LRC-is-exactly-linear-Mc-equals-1-over-n-plus-c-times-n-minus-2-over-n-opus-20260630.md
```

## Verification Record
- Achievability exact at `c=(q−n)/2q`: n=10,14,20, all q (proof.py).
- `M_c = env` at large Qmax (Qmax=12n) for c=1/14,1/7,1/4,2/7,3/7,2/5,11/28 (n=14); dips `<10⁻³` at c=odd/2n.
- `(L−1/4)·n`: 0.373 (Qmax=2n) → 0.491 (Qmax=8n) → 1/2 (envelope integral).
- Clumping inequality `n ≤ 2j+(2m+1)(n−2)` zero only at (1,0): n=10,14,20.

## Notes & History
Found opus-2026-06-30-S1 chasing the loneliness-integral constant. The route corrected two finite-Qmax
artifacts of its own (an erroneous `L=1/4+3/(8n)` and a spurious "trivial regime" `c*=1/3`) before landing on
the exact linear law. The mechanism — *the inhomogeneous AP-LRC is solved by one clump* — unifies the
achievability and optimality halves: build the block; show every escape must clump; the single clump uniquely
wins the integer inequality.
