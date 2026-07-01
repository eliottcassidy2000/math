---
id: THM-591
title: The inhomogeneous Lonely-Runner gap of the AP is exactly linear in the observer position
status: PROVED (both directions ELEMENTARY -- M_c>=env via rational-q achievability + 1-Lipschitz; M_c<=env via the avoided-arc-edge argument)
date: 2026-06-30
session: opus-2026-06-30-S1
depends_on:
  - THM-405   # homogeneous AP M=1/n (the c=0 case)
  - THM-523   # q-witness / covering-set reduction (the upper bound is its dual)
  - THM-360   # unit-endpoint divisibility filter (t=a/n)
  - THM-386   # far_iff_fract: the avoided-arc characterization (Lean-checked)
related:
  - OPEN-Q-108
  - HYP-3747
  - HYP-3748
external: three-distance (Sós–Steinhaus) theorem; totient summatory A002088
---

# THM-591 — The inhomogeneous AP Lonely-Runner gap is the line c + (1−2c)/n

## Statement

Let `S = {1, 2, …, n−1}` (the arithmetic-progression speed set, the LRC(n) tight extremal of THM-405), and
for an **observer position** `c ∈ [0,1]` define the **inhomogeneous LRC gap**
> `M_c(S) = max_{t∈[0,1)} min_{v∈S} ‖vt − c‖`   (`‖·‖` = distance to nearest integer; `M_0 = M(S)` is THM-405).

Then, for `c ∈ [0, 1/2]` (and `M_c = M_{1−c}` by symmetry),
> **`M_c(S) = 1/n + c·(n−2)/n = c + (1−2c)/n`.**

The identity is **exact for every `c ∈ [0,1/2]`** (no error term): achievability holds on the dense set
`c=(q−n)/2q` for all **rational** `q=a/b ≥ n` (§A), and `M_c` is `1`-Lipschitz in `c`, so `M_c ≥ env`
everywhere; with `M_c ≤ env` (§B), `M_c = env`. The earlier "`O(1/n²)` dips at `c=odd/2n`" were **finite-Qmax
artifacts** — the optimal `t` there has denominator `n²` (e.g. `c=1/28, n=14: t=183/196`), so `Qmax<n²`
under-reports. By dilation-invariance `M_c(dS)=M_c(S)` the law holds for every dilated AP `d·{1,…,n−1}`
(including divisible extremals such as `2·{1,…,13}`, which contains `14 ≡ 0 (mod 14)`).

## A. Achievability (lower bound) — the block construction [elementary, exact]

For `c = (q−n)/(2q)`, integer `q ≥ n`, take **`t = (q−1)/q`**. Then `vt = k(q−1)/q ≡ −k/q = (q−k)/q` for
`v=k∈{1,…,n−1}`, so the runners are the **block** `{(q−n+1)/q, …, (q−1)/q}` — `n−1` consecutive multiples
of `1/q`. Its complement is a single gap of length `(q−n+2)/q`, and `c` is its exact center (both nearest
runners at distance `(q−n+2)/(2q)`). Hence `min_v ‖vt−c‖ = (q−n+2)/(2q)`; substituting `q = n/(1−2c)` gives
`= 1/2 − (n−2)/(2q) = 1/n + c(n−2)/n`. ∎  *(Verified n=10,14,20, every q.)*

**Extension to all rational `q` (kills the dips).** For `q = a/b ≥ n` (`gcd(a,b)=1`), take `t = (a−b)/a`.
Then `vt − c ≡ −(vb/a + c)`, so `‖vt−c‖ = ‖(a+wb)/(2a)‖` with `w = 2v−n ∈ {−(n−2),…,n−2}` (step 2); since
`a ≥ nb` the nearest multiple of `2a` to `a+wb` is `0` or `2a` at distance `a−|w|b`, minimized at `|w|=n−2`,
so `min_v ‖vt−c‖ = (a−(n−2)b)/(2a) = (q−n+2)/2q = env`. ∎ The achieving `c=(a−nb)/2a` are **dense** in
`[0,1/2)`. Since `M_c` is `1`-Lipschitz in `c` (a sup of `1`-Lipschitz `min_v ‖vt−c‖`), `M_c ≥ env` on a dense
set ⇒ `M_c ≥ env` **everywhere** — the law is exact, no dips. (See `DIP-UPGRADE-…` reflection.)

The optimal lonely time (mirror form `t=1/q`, `q=n/(1−2c)`) **slides `1/n → 0`** as the observer walks
origin → antipode: `c=0 ⇒ t=1/n` (the THM-405 covering-min), `c→1/2 ⇒ t→0` (the trivial antipode).

## B. Upper bound (optimality) — the avoided-arc-edge argument [ELEMENTARY, fully rigorous]

Set `q := n/(1−2c)` (so `1/q = (1−2c)/n`). The one identity that closes it:
> **`c − env = (2c−1)/n = −1/q`**  and  **`1 − 2·env = (n−2)/q`.**
So the observer's avoided arc `[c−env, c+env]` has **left endpoint exactly `−1/q`**; mod 1 it is
`[0, c+env] ∪ [1−1/q, 1)`, which **contains `(−1/q, 1/q)`** (since `1/q ≤ c+env ⟺ 0 ≤ 2cn`).

**Proof.** Suppose `min_v ‖vt−c‖ > env` for some `t`. Then no runner lies in `[c−env, c+env]`, so the `n−1`
runners lie in the complementary **open** arc of length `1−2·env = (n−2)/q`. Their `n−2` consecutive gaps sum
to *strictly less* than `(n−2)/q`, so the smallest is `< 1/q`; being the arc-distance between two consecutive
runners it equals `‖jt‖` with `j ∈ {1,…,n−2}`. Hence `‖jt‖ < 1/q`, i.e. `jt mod 1 ∈ (−1/q, 1/q) ⊆` the
avoided arc — so the runner `jt` has `‖jt − c‖ ≤ env`, contradicting the hypothesis. ∎

No three-gap theorem, no cluster widths, no `Qmax`. The `c=0` case recovers the classical homogeneous tightness
`M_0 ≤ 1/n` (THM-405): the forced runner with `‖jt‖ < 1/n` lands in `[−1/n,1/n]`. This is the *dual* of the
q-witness (THM-523): the forced near-`0` runner is exactly what the observer's `1/q`-aligned danger arc
catches. (See `UPPER-BOUND-CLOSED-…` reflection. Supersedes the earlier clumping-inequality argument, which had
an `O(1/n²)` structural slop — now removed.)

## C. Corollaries

1. **`c=0` recovers THM-405:** `M_0 = 1/n`, the homogeneous covering-min, the unique hardest observer.
2. **Loneliness integral (exact):** `L(S) = ∫₀¹ M_c dc = 2∫₀^{1/2}[1/n + c(n−2)/n]dc = 1/4 + 1/(2n)`
   **exactly** (the dip upgrade removes the error term). The mean-observer loneliness is `1/4 + 1/(2n)`, with
   `1/4 = ∫‖c‖dc` the trivial nearest-integer floor; the LRC's *hard part* is the `1/(2n)` excess, localized
   at `c=0` (`M_0 = 1/n`).
3. **The two SC fixed points** `{0, 1/2}` (fixed by `c↔−c` = time-reversal = complement) are the endpoints of
   the line: origin (bonus `1/n`) and antipode (bonus `0`).
4. **Combinatorial skeleton:** the optimal-`t` regimes (orderings of `{vt}` as `t` sweeps) number `Φ(n−1)`
   (the Farey length, A002088 — see HYP-3748), are the Sós–Steinhaus three-distance regimes, organized by the
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
04-computation/lrc_dip_upgrade_rational_q_lipschitz_opus_20260630.py
04-computation/lrc_upper_bound_avoided_arc_edge_PROOF_opus_20260630.py
07-reflections/UPPER-BOUND-CLOSED-the-avoided-arc-edge-is-minus-1-over-q-so-the-forced-runner-lands-inside-it-opus-20260630.md
07-reflections/PROOF-the-inhomogeneous-AP-LRC-linearity-the-block-construction-gives-Mc-equals-1-over-n-plus-c-n-minus-2-over-n-opus-20260630.md
07-reflections/CORRECTION-the-inhomogeneous-AP-LRC-is-exactly-linear-Mc-equals-1-over-n-plus-c-times-n-minus-2-over-n-opus-20260630.md
07-reflections/DIP-UPGRADE-the-inhomogeneous-AP-LRC-law-is-EXACT-rational-q-achievability-plus-lipschitz-kills-the-dips-opus-20260630.md
```

## Verification Record
- Achievability exact at `c=(q−n)/2q`: integer q n=10,14,20; rational q=a/b n=8..20 (dipproof.py).
- Dips are Qmax<n² artifacts: `M_{1/28}(n=14)` rises 0.10023(Q=3n) → 0.10198(Q=12n) → 0.10204=5/49=env(Q=n²);
  optimal t=183/196 has denominator n²=196. `M_c` confirmed 1-Lipschitz (max |ΔM/Δc|=1.000).
- `(L−1/4)·n`: 0.373 (Qmax=2n) → 0.491 (Qmax=8n) → 1/2 (= envelope integral, exact).
- Clumping inequality `n ≤ 2j+(2m+1)(n−2)` zero only at (1,0): n=10,14,20.

## Notes & History
Found opus-2026-06-30-S1 chasing the loneliness-integral constant. The route corrected two finite-Qmax
artifacts of its own (an erroneous `L=1/4+3/(8n)` and a spurious "trivial regime" `c*=1/3`) before landing on
the exact linear law. The mechanism — *the inhomogeneous AP-LRC is solved by one clump* — unifies the
achievability and optimality halves: build the block; show every escape must clump; the single clump uniquely
wins the integer inequality.
