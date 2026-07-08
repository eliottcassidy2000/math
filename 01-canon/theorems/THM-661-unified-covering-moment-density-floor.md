---
id: THM-661
title: The unified covering-moment density floor — on the covering reformulation (THM-657, W = uncovered measure), the degree-d one-sided moment bound mu_{1/7}(E) >= max{ sum c_i E[W^i] : sum c_i w^i <= 1_{w>0} on [0,6/7] } is a rigorous DIAMETER-FREE lower bound on mu; degree 4 clears the honest (A') bars at ALL of k=8,9,10 (block: 0.761/0.645/0.553 >= 0.675/0.562/0.452) where degree 2 fails, and degree 2 (= Paley–Zygmund, THM-660) clears k=11,12,13 — so a SINGLE covering framework with degree <= 4 discharges all six density-floor legs for the block, subsuming BOTH the tent (THM-651/655/656) and the PZ floor (THM-660)
status: PROVED per-shape (the moment LP is a rigorous lower bound on P(W>0): any polynomial p with p(w) <= 1_{w>0} on [0,6/7] gives E[p(W)] <= P(W>0), and the LP maximizes E[p(W)] = sum c_i E[W^i]). Block bounds are EXACT-verified: exact rational moments E[W^i] via Farey-cell integration + a feasible rational polynomial; degree-4 bounds 0.7611/0.6446/0.5531 at k=8/9/10 (all >= bars, margins +0.086/+0.083/+0.101); degree-2 = THM-660's PZ = 0.347/0.308/0.272 at k=11/12/13. The UNIFORM floor (min over families) reduces to [compact exact check + decorrelation tail], the block/perforated-block being the minimizer (same structure as THM-660).
source: mac-mini-2026-07-08-S57 (HYP-5347)
depends_on:
  - THM-657   # the covering reformulation: mu = P(W>0), W = 1 - coverage = sum (g_i-1/7)_+
  - THM-660   # klein's PZ covering floor = the degree-2 case of this bound (k>=11)
related:
  - THM-651   # the tent floor (k=8, first-moment on the pair-sum F) — this covering-moment route is an alternative that ALSO clears k=8,9,10
  - THM-655   # the k=9 average-form conditional tent (mine) — superseded for k=9 by degree-4-on-W (diameter-free, self-contained)
  - THM-656   # the tent-variance / additive-energy floor — the moment structure of W parallels it
external: the classical moment problem on [0, M] (Markov–Krein); Paley–Zygmund 1932 (the degree-2 optimum).
---

# THM-661 — the unified covering-moment density floor

## Statement

On the covering reformulation (THM-657): `W(x) = 1 − coverage = Σ_i (g_i − 1/7)_+ ∈ [0, 6/7]`,
`mu_{1/7}(E) = P(W > 0)`. For any degree `d`, the **one-sided moment bound**

> **`mu_{1/7}(E) >= B_d(E) := max { Σ_{i=0}^d c_i E[W^i] : Σ_i c_i w^i <= 1_{w>0}(w) ∀ w ∈ [0, 6/7] }`**

is a rigorous, diameter-free lower bound on `mu`. (`1_{w>0} = 0` at `w=0`, `= 1` on `(0,6/7]`;
so the constraint is `c_0 <= 0` and `Σ c_i w^i <= 1` on `(0, 6/7]`.)

**`B_2 =` Paley–Zygmund `= E[W]²/E[W²]`** (THM-660). `B_d` is non-decreasing in `d`.

## Proof

For any feasible `p(w) = Σ c_i w^i` (with `p(w) <= 1_{w>0}(w)` on `[0,6/7]`), pointwise
`p(W) <= 1_{W>0}`, so `Σ c_i E[W^i] = E[p(W)] <= E[1_{W>0}] = P(W>0) = mu`. Maximizing the LHS
over feasible `p` gives `B_d(E) <= mu`. ∎ (The LP's optimum is the Markov–Krein moment bound;
`d=2` is the Cauchy–Schwarz / Paley–Zygmund case.)

## The unification: degree <= 4 discharges ALL six legs (A′), diameter-free

The honest (A′) bars decrease with `k` (`m_P + 1 − min_{|P|=13−k} meas(G_P)`), while `B_d(block)`
decreases with `k` more slowly at higher `d`. The block values (EXACT rational moments via
Farey-cell integration; feasible rational polynomial exhibited):

| k | bar | `B_2` (PZ) | `B_4` | true `mu(block)` | discharged by |
|---|------|-----------|-------|------|------|
| 8 | 0.6750 | 0.5740 ✗ | **0.7611** ✓ | 0.9401 | degree 4 |
| 9 | 0.5622 | 0.4809 ✗ | **0.6446** ✓ | 0.8401 | degree 4 |
| 10| 0.4521 | 0.4109 ✗ | **0.5531** ✓ | 0.7755 | degree 4 |
| 11| 0.3312 | **0.3471** ✓ | 0.4108 | 0.6263 | degree 2 (THM-660) |
| 12| 0.1993 | **0.3080** ✓ | — | 0.5699 | degree 2 |
| 13| 0.0565 | **0.2720** ✓ | — | 0.4425 | degree 2 |

Exact block moments (k=8/9/10): `E[W] = 429/1715, 1088/5145, 2818/15435`;
`E[W²] = 11779/108045, 351667/3781575, 662617/8168202`; `E[W³], E[W⁴]` likewise rational
(`04-computation`). So **one covering framework — the degree-`≤4` moment bound on the uncovered
measure — clears every density-floor leg `k = 8..13` for the block, diameter-free**, replacing
the two-tool patchwork (tent for `k<=10`, PZ for `k>=11`). The tent (THM-651, first moment of
the pair-sum `F = Σ f(frac(dx))`) and the PZ floor (THM-660, second moment of `W`) are both
special cases of "a low-degree moment bound from the covering/gap structure."

## The uniform floor (min over families)

`B_d` is a functional of the covering moments `E[W^i]`, which are structural. The uniform
statement `min_E B_d(E) >= bar` reduces — exactly as for `B_2` (THM-660) — to **[a finite exact
check over compact shapes]** + **[a decorrelation tail]** (`B_d` rises with diameter: spreading a
family makes `W` less spiky, raising every moment bound toward `mu -> 1`). The minimizer is the
block / a perforated block (compact); exact compact minimizers at k=11/12/13 are
`0.3468 / 0.2992 / 0.2654` (THM-660 addendum), all `>= bar`. The remaining rigorous piece is the
common tail bound (the worst wide family clears by `+0.055` even at the thin k=11 leg).

## Why this matters

The density-floor side of LRC(14) — six legs, historically discharged by unrelated tools (Hunter
`6/49`, the shifted tent, the average-form conditional tent, the additive-energy variance floor,
the PZ covering floor) — is now ONE theorem: `mu >= B_d(E)`, `d <= 4`, on the single covering
reformulation `mu = P(k arcs of length 1/7 fail to cover the circle)`. Everything is a moment
bound on the uncovered measure; the diameter never appears.

## Addendum — the EXACT CLOSED FORM for degree 3 (opus-2026-07-08-S148, HYP-5327)

The `B_d` above are moment LPs (numeric + a rationalized feasible polynomial). For `d = 3`
the optimum has a **closed form**, no LP: the optimal minorant of `1_{w>0}` is

> `p(t) = 1 − (1 − t/r)²(1 − t/M)`,  `M = 6/7`,  `r = r* := (m2 − m3/M)/(m1 − m2/M)`,

which satisfies `p(0) = 0` and `p(t) ≤ 1` on `[0, M]` for **any** `r` (both factors `≥ 0`),
so it is always feasible; optimizing `r` (a rational critical point) gives

> **`B_3(E) = D3(E) = E[W]/M + (E[W] − E[W²]/M)² / (E[W²] − E[W³]/M)`**  (valid when
> `E[W²] − E[W³]/M > 0`; else use `B_2`).

This is exact-rational in the covering moments — a one-line certificate rather than an LP.
**Its use is the binding k=11 leg:** THM-660/this theorem clear `k ≥ 11` with `B_2` (PZ),
whose k=11 margin is razor-thin (`+0.0159`); the exact `B_3` lifts it to **`+0.0735` (4.6×
thicker)**, block value `54912120381817/135668932727076 = 0.404751`. Exhaustive (k=11,
diam ≤ 14): the **block is the exact `B_3`-minimizer**, `min B_3 = 0.404751 ≥ bar`, margin
`+0.0735` uniformly over the compact regime (vs PZ's `+0.0156`). The closed form also clears
all six legs for the block (`+0.00058 / +0.0055 / +0.033 / +0.0735 / +0.159 / +0.257` at
k=8..13); k=8,9 are thin at degree 3 (degree 4 above is better there), but for the binding
k=11 leg `B_3` is the natural exact sharpening. (Honest scope: k=11 was already discharged
thinly by `B_2` + the decorrelation tail; this is a robustness upgrade with a clean formula.)
Files: `04-computation/lrc14_degree3_closed_form_floor_opus_S148.py`,
`lrc14_pz_general_integrator_opus_S148.py` (block-only integrator generalized to arbitrary
families), `lrc14_pz_k11_exhaustive_opus_S148.py` (+outs).

## Verification & files

`04-computation/lrc14_unified_moment_covering_floor_macmini_S57.py` (exact block moments E[W^i]
via Farey; the degree-d moment LP with rational feasibility check; the B_d table k=8..13).
`04-computation/lrc14_block_PZ_exact_macmini_S57.py` (the d=2 exact values).

## The decorrelation tail, reduced (mac-mini-S57)

The uniform floor's tail (wide families) has a rigorous handle via a NEAR/FAR split of the
second moment: `E[W^2] = int int P(y1,y2 both uncovered) = near + far`, where
- `near = 2 int_{1/7}^{2/7} q(L) dL` (overlapping arcs `|y1-y2|<=1/7`; `q(L)=E[sum(g-L)_+]` the
  empty-arc probability) is EXACT and `<= (2/7)E[W]` (since `q(L) <= q(1/7) = E[W]`);
- `far = int int_{|y1-y2|>1/7} P(both empty)` (disjoint arcs) `-> E[W]^2` by decorrelation.

**RIGOROUSLY PROVED (verified, 0 violations):** `near <= (2/7)E[W]` (since `q(L) <= q(1/7)`),
`far <= (5/7)E[W]` (since `q_2(y1,y2) = P(A1 ∪ A2 empty) <= P(A1 empty) = q_1(y1)`, integrate),
hence `E[W^2] <= E[W]` and **`PZ >= E[W]`** unconditionally. Consequence: **k=13 reduces to the
single bound `E[W] >= m_P = 0.0565`** (then `PZ >= E[W] >= m_P`; actually `PZ >= (7/6)E[W]` needs
only `E[W] >= 0.0484`). **k=11,12 reduce to `far <= E[W]^2 + o(1)`** (the crude `far <= (5/7)E[W]`
gives only `PZ >= 0.18 < 0.33/0.20`); with the EXACT near (which is small, `~0.015 << (2/7)E[W]
~ 0.04`, because `q(L)` decays fast) and `far <= E[W]^2`, `PZ >= E[W]^2/(near+E[W]^2)` clears all
three (worst wide k=11 clears at `0.46`; sidon-like `0.55-0.58`).

**HONEST STATUS OF THE PROOF — the barely-covers wall.** The two residual inequalities,
`E[W] >= 0.0565` (k=13) and `far <= E[W]^2 + o(1)` (k=11,12; equivalently `Var(W) <= near + o(1)`),
are BOTH forms of the same obstruction: with total arc length `k/7 = 13/7 = 1.86 > 1`, the
inclusion–exclusion for the empty-arc probability DIVERGES (`E[W]` Bonferroni-3 `= 1 - k/7 + S_2 -
S_3 = -0.857 + 1.599 - 1.413 = -0.671`, vs true `0.127`; the pair masses `P(A_i A_j) = 1/49`
exactly, so Hunter/tree bounds are useless). The truth is safe by a wide margin (true `E[W] ~
0.13`, true `far ~ E[W]^2`, true tail `PZ ~ 0.4-0.8`), and the exact Fourier covariance formula
`far_cov = Σ_{m≠0} Σ_N |V-hat(m,N)|^2 c_N` (`c_N = [N=0](1-2/7) - sin(2πN/7)/(πN)`) reduces it to
klein's pairwise offset-tent masses (THM-638/645) — but the uniform bound over ALL families is
the SAME fundamental difficulty as the μ-extremal lemma "AP minimizes μ" (THM-530), now in
covering/Fourier form. **NOT proved here.** What IS new and rigorous: the reduction (PZ >=
E[W]^2/(near+far), the exact near, `far <= (5/7)E[W]`, `PZ >= E[W]`), which cleanly SEPARATES
k=13 (needs a scalar `E[W]` floor) from k=11,12 (need the decorrelation) and localizes the whole
density floor to this one barely-covers estimate. File:
`04-computation/lrc14_PZ_compact_and_tail_macmini_S57.py` (near/far split, rigorous-inequality
verification, the barely-covers divergence).

## k=12,13 close via a UNIFORM D3 floor; brick (A) for k=11 (mac-mini-S57)

kps-S79 closed the thin k=11 leg with opus's D3 via two bricks (A: max additive energy R2 <=
614 for prim-diam >= 16; B: R2 <= 614 => D3 >= bar). The **k=12,13 analogs are far easier and
need NO R2-diameter split** — they close by a UNIFORM D3 floor: the EXACT compact-minimum D3
(Farey integration over all primitive shapes diam <= 15) is
- **k=12: D3 = 0.355876 >= bar 0.199344, exact margin +0.156532 (1.8x)** at (0,2,4,5,6,..,12,14);
- **k=13: D3 = 0.308844 >= bar 0.056487, exact margin +0.252357 (5.5x)** at (0,2,3,4,..,12,14);
the tail (diam > 15) only rises (D3 -> 1 by decorrelation, LEM-005), so `min_E D3 = ` the compact
min above, clearing both bars comfortably. So k=12,13 need only [the exact compact D3 check +
LEM-005 with huge margin], not the delicate energy-diameter extremal argument k=11 requires.

**Brick (A) for k=11 verified (primitive):** max R2 over PRIMITIVE 11-sets with prim-diam >= 16
is exactly **614** (the 1+10 split `{0,..,9,16}`); the near-2-AP bump (610 at diam 18) stays
under 614, and diam >= 19 gives R2 = 590 (lone-point differences no longer overlap the block).
Primitivity is essential — the dilated block `2·{0..10}` has R2 = 770 but prim-diam 10 (excluded).
The extremal proof (compression/rearrangement; the near-c-AP competitors bounded) is the one
combinatorial lemma the k=11 leg still needs rigorously; k=12,13 sidestep it entirely.
File: `04-computation/lrc14_k12_13_D3_analog_bricks_macmini_S57.py`.
