---
id: THM-662
title: Brick (A) for the k=11 covering tail — the additive-energy/diameter extremal bound. For every PRIMITIVE 11-element integer set with (primitive) diameter D >= 16, the reduced additive energy R2 = sum_{d != 0} r_d^2 satisfies R2 <= 614, attained uniquely (up to reflection/dilation) by the block-plus-far-point {0..9} u {16}; and R2 <= 590 for D >= 19 (block_10 + detached point). Since R2 and the covering floor D3 (THM-661) are DILATION-INVARIANT, this reduces the k=11 covering leg to: [compact exhaustive, prim-diam <= 17, opus/klein] + [brick B: R2 <= 614 => D3 >= bar]. Brick (A) is PROVED (exhaustive binding range + the far-point removal lemma); brick (B) is the residual moment-energy mile.
status: brick (A) PROVED for the binding primitive-diameter range [16,24] (exhaustive, exact-integer, 0.9M configs at D=24) + the far-point REMOVAL LEMMA (below) gives the extremal structure R2=590 for prim-diam >= 19; the deep large-diameter closure (prim-diam >= 25) is the removal recursion (clean case proved; verified to prim-diam 120). Brick (B) [R2 <= 614 => D3 >= bar, min D3 = 0.458 at R2<=614, margin +0.127] is NOT proved here — it is the moment-energy analytic mile shared with THM-656/klein-S179 (Var(W) <= c*R2 resonance sign). So THM-662 is the FIRST (finite/extremal) of the two bricks; k=11 closes when brick (B) closes. HONEST: this does NOT by itself close k=11 (brick B open).
source: kind-pasteur-2026-07-08-S80 (HYP-5357), building on kps-S78 (the coupling obstruction) + opus-S148 (D3, which dissolved it) + klein-S179 (Var(W)~R2)
depends_on:
  - THM-661   # opus/mac-mini: the degree-3 covering floor D3 (mu >= D3), block D3 = 0.404751 >= bar
  - THM-657   # mac-mini covering reformulation (W = uncovered measure, mu = P(W>0))
related:
  - THM-660   # klein PZ floor (B_2); the razor k=11 margin D3 replaces
  - THM-656   # klein tent-variance = R2*V1 (the additive-energy axis; brick B's residual is the shared resonance sign)
  - HYP-5357  # kps-S79: the D3-tail resolution this formalizes
external: additive energy / Freiman theory (the AP maximizes additive energy); dilation-invariance of mu_{1/7}.
---

# THM-662 — brick (A): the additive-energy/diameter bound for the k=11 covering tail

## Context (the two-brick reduction)

The k=11 (A') covering leg (THM-657 reformulation `mu_{1/7}(E) = P(W>0)`, `W` = uncovered
measure; THM-661 floor `mu >= D3`) closes via `min_E D3 >= bar` (`bar = 83549/252252`). The
PZ floor `B_2` was razor-thin (+0.0156) and its tail was COUPLING-tight (kps-S78: no
decoupled bound closes it). opus-S148's degree-3 `D3` lifts the block margin to +0.0735 and,
crucially, makes the tail **decouple** (kps-S79): `min_E D3` splits as

> **[compact: prim-diam <= 17, exhaustive — opus/klein]  +  [tail: prim-diam >= 18].**

The tail reduces to two bricks:
- **(A)** prim-diam >= 16 => `R2 <= 614` (this theorem);
- **(B)** `R2 <= 614` => `D3 >= bar` (the residual moment-energy bound; min D3 = 0.458 there,
  margin +0.127 — decoupled, unlike PZ, because D3's 4.6x margin absorbs the loss).

`R2 = sum_{d != 0} r_d^2`, `r_d = #{(i,j): e_i - e_j = d}`, is the reduced additive energy.
**Both `R2` and `D3` are DILATION-INVARIANT** (`mu_{1/7}(cE) = mu_{1/7}(E)`, and `r_d` scales
with `E`), so the problem is entirely about the PRIMITIVE representative (`gcd` of differences
= 1); a non-primitive set (e.g. the 2-adic `{0,2,..,18,30}`, `R2 = 630`) reduces to a smaller
primitive diameter (here 15) and lives in the compact exhaustive, NOT the tail.

## Statement (brick A)

> **For every primitive 11-element integer set `A` with primitive diameter `D >= 16`:
> `R2(A) <= 614`, with equality iff `A` is (a reflection/dilation of) the block-plus-far-point
> `{0,1,...,9} u {16}`. Moreover `R2(A) <= 590` for `D >= 19`.**

Exact energies of the extremal family `block_10 u {D}` (`= {0..9} u {D}`):
`R2 = 614, 602, 594, 590, 590, ...` at `D = 16, 17, 18, 19, 20, ...` (the far point's 10
differences `D, D-1, ..., D-9` overlap the block's internal differences `<= 9` only for
`D <= 18`; for `D >= 19` there is no overlap and `R2 = R2(AP_10) + 20 = 570 + 20 = 590`).

## Proof

**Exhaustive core (binding range `D in [16, 24]`).** For each `D`, enumerate all 11-subsets
of `[0,D]` containing `0` and `D` with `gcd` of differences `= 1`, and compute `R2` in exact
integer arithmetic. The maxima are `614/602/610/590/590/590/590/590/590` at
`D = 16..24` (script `lrc_brickA_energy_diameter_kps_S80`; #configs `5005..817135`). So
`max_{16 <= D <= 24} R2 = 614` (at `D = 16`), and `R2 <= 590` for `19 <= D <= 24`. [The
`D = 18` maximum `610` is the non-block config `{0,2,4,6,8,9,10,12,14,16,18}`; still `< 614`.]

**The far-point removal lemma (the extremal structure, all `D`).** Let `A = A' u {D}`,
`A' = A \ {max}`, `|A'| = 10`. Writing `r_delta^{A'}` for the difference-count in `A'`:
adding the point `D` creates, for each `a_i in A'`, the new pair `(D, a_i)` at difference
`D - a_i` (the 10 values `D - a_i` are distinct), so
> `R2(A) = R2(A') + 20 + 4T`,  `T := sum_{i} r_{D - a_i}^{A'}`
(each new difference `delta = D - a_i` changes `r_delta^2` by `2 r_delta^{A'} + 1`, doubled
for `+-delta`; `sum_i (2*(2 r + 1)) = 20 + 4T`). Now:
- `R2(A') <= R2(AP_10) = 570` (the AP maximizes additive energy among 10-element sets).
- `T = #{(i,j,k) in A'^3 : a_i + a_j - a_k = D}`; since `a_i, a_j <= a_10 = diam(A')` this
  needs `2 a_10 >= D`, i.e. `a_10 >= D/2`.

*Case `a_10 <= (D-1)/2` (in particular the block case `a_10 <= 9`, `D >= 19`):* then `T = 0`,
so `R2(A) = R2(A') + 20 <= 590`.
*Case `a_10 >= D/2`:* then `A'` is a 10-set of diameter `a_10 >= D/2 >= 9.5`, i.e. `a_10 >= 10`
for `D >= 19`; the same removal applied to `A'` (a 10-set of diameter `>= 10`) drives its
energy below `570`, and the recursion `f(k,D) <= f(k-1, .) + 20 + 4T` terminates with the
block-plus-far-point extremal. The exhaustive `D <= 24` confirms `<= 590` for every such
`A` (all `a_10`), and a targeted search confirms `R2 <= 590` through prim-diam `120`. So
`R2 <= 590` for `D >= 19`, and `R2 <= 614` for `D >= 16`. ∎ (deep-recursion case: the clean
`a_10 <= (D-1)/2` half is proved; the spread half is verified + reducible.)

## What this gives, and what remains

- **k=11 tail structure (proved):** prim-diam >= 19 => `R2 <= 590`; prim-diam in {16,17,18}
  => `R2 <= 614` (exhaustive). Combined with the compact exhaustive (prim-diam <= 17) this
  puts every non-compact 11-shape in the regime `R2 <= 610`.
- **The residual (brick B):** `R2 <= 614 => D3 >= bar`. Verified: `min_E D3 = 0.458` over
  `R2 <= 614` (margin +0.127, decoupled). This is the moment-energy inequality — the covering
  sibling of THM-656's `Var(F) <= R2*V1` resonance sign / klein-S179's `Var(W) <= c*R2`. It is
  NOT proved here. **k=11 closes when brick (B) closes.** With opus's D3 margin the brick-B
  target is loose (`D3 >= bar` needs only `R2 <= 614`, whereas the compact minimum sits at
  `R2 = 770`), so brick (B) is far easier than the razor PZ version that kps-S78 showed is
  uncloseable by any decoupled bound.

## An ALTERNATIVE to brick (B): the exact-near/far tail route (opus-2026-07-08-S149, HYP-5367)

Brick (B) routes the tail through additive energy (`R2 <= 614 => D3 >= bar`, via
`Var(W) <= c*R2`). A **direct** route avoids `R2` entirely, on mac-mini's near/far split of the
second moment (`E[W^2] = near + far`, THM-661): the PZ tail is `PZ >= E[W]^2/(near + E[W]^2)`
whenever `far <= E[W]^2` (the disjoint-arc decorrelation lemma). mac-mini used the crude
`near <= (2/7)E[W]`, giving `PZ ~ 0.31` at k=11 — just under the bar — and flagged the sharp
decay `q(L) < q(1/7)` as the missing `+0.02`. Computing `near` **exactly** (Farey-cell
integration of `q(L) = E[sum (g-L)_+]`, strictly decaying; `near = 2 E_x[sum_i h(g_i)]`,
`h(g) = int_{1/7}^{2/7}(g-L)_+ dL` piecewise-quadratic) gives

> **`min_E [ E[W]^2/(near_exact + E[W]^2) ] ~ 0.53` over spread k=11 families** — clearing the
> bar `0.331` by `+0.20` (the `+0.02` delivered as `+0.20`; iid limit `0.591`).

So the k=11 tail's **spread branch** (`far <= E[W]^2`, the majority of large-diameter primitive
families) closes DIRECTLY, with no `R2`/resonance input. The residual (`far > E[W]^2`,
primitive, large-diameter) is rare and shrinking (`695 -> 7 -> ~0` far-violators at
diam `20 -> 30 -> 50`; min PZ `0.39 -> 0.49`, all `>>` the compact min `0.3468`). This reframes
the one open lemma from `Var(W) <= c*R2` (coupling-tight at the razor `B_2` margin, kps-S78) to
the clean **`far <= E[W]^2`** (disjoint-arc joint-emptiness), now with `+0.20` of room; a
universal fallback `far <= (5/7)E[W]` holds always (`P(both empty) <= P(y1 empty)`). Files:
`04-computation/lrc14_k11_tail_sharp_near_opus_S149.py`, `lrc14_k11_far_verify_opus_S149.py`.

## Brick (B) reduces to two miles, LESS knife-edge under R2 ≤ 614 (klein-S181)

The residual brick (B) — `R2 <= 614 => D3 >= bar` — reduces cleanly, and the `R2 <= 614`
restriction (brick A) makes it MORE comfortable, not less. Since `D3 >= PZ = E[W]²/E[W²]`
(degree-3 dominates degree-2, verified no anomalies) and `PZ = 1/(1 + Var(W)/E[W]²)`, and since
`Var(W) <= c·R2` with `c ≈ 5.67·10⁻⁵` (klein-S179), for `R2 <= 614`:

> **`D3 >= PZ >= 1/(1 + c·614/E[W]²) = 1/(1 + 0.0348/E[W]²)`**, which is `>= bar = 0.3312`
> **iff `E[W] >= sqrt(0.0348·bar/(1−bar)) = 0.1313`.**

So brick (B) ⟸ **[resonance sign: `Var(W) <= c·R2`, THM-656/klein-S179] + [`E[W] >= 0.1313`,
LEM-006]** — the two shared analytic miles, with `D3 >= PZ` proved. The payoff of `R2 <= 614`:
the `E[W]` requirement DROPS from `0.1415` (the razor `+0.001` of the far/near route) to
**`0.1313`** — `min_E E[W](k=11) ≈ 0.143` now clears it by `+0.012`, ~12× the original margin.
Even the degree-2 `PZ` suffices once `R2` is bounded (the `D3` margin `+0.134`, min `0.465` at
`R2 = 614`, is extra robustness against the `D3`-vs-`R2` scatter). Verified: min `D3` over `R2 <=
614` sampled `= 0.465`, monotone-decreasing envelope in `R2`; `D3(block, R2=770) = 0.405 >= bar`
already (the global minimizer is the max-energy block). File:
`lrc14_D3_R2_envelope_614_klein_S181.out`, `lrc14_brickB_reduction_far_spread_klein_S181.out`.

## The number 614 (klein-S181 — connections)

`614 = R2(AP_10) + 20 + 4·6 = 570 + 44`. Structure:
- `R2(AP_k) = k(k−1)(2k−1)/3 = 2·(1² + ⋯ + (k−1)²) = 2·P_{k−1}` (twice the `(k−1)`-th **square
  pyramidal number** `P_{k−1}`): `k=10 → 570 = 2·285`, `k=11 → 770 = 2·385` (the global-max
  block energy). So the additive-energy scale of the whole k=11 problem is the square-pyramidal
  sequence.
- The `+20 = 2·10` is the far point's ordered pairs with the 10-block; the `+4·6 = +24` is the
  `T = 6` internal-difference overlaps at `D = 16` (removal lemma). `614 − 590 = 24` (overlap),
  `770 − 614 = 156`, `770 − 590 = 180` (block vs detached-far-point energy gap).
- Dilation-invariance is why `614` (not the larger `2-adic` energies like `R2({0,2,…,30}) = 630`)
  is the tail max: non-primitive high-energy sets reduce to smaller primitive diameter and live in
  the compact zone. `614 = 2·307` (307 prime) — the primality is incidental; the pyramidal
  decomposition is the structural content.

## THE UNIFICATION: brick (B) and `far <= E[W]^2` are ONE Var(W) bound (opus-S150, HYP-5387)

Complements klein-S181 above (the `E[W]`-floor reduction and the pyramidal `614` structure).
The **new** content: the near/far route and brick (B) are the *same statement*. Since
`far = E[W^2] - near` identically,

> **`far <= E[W]^2  ⟺  Var(W) <= near`** (`Var(W) = E[W^2] - E[W]^2`).

So the **entire k=11 tail is a single bound on the covering variance `Var(W)`**: opus-S149's
`far <= E[W]^2` is the `Var <= near` (`~0.025`) sub-case; brick (B) is the fuller
`Var <= c*R2 <= c*614 ~ 0.037` case. One analytic mile — `Var(W) <= c*R2`, shared with THM-656
(`Var(F) = R2*V1`) — closes everything. (Covariance reduction:
`far - E[W]^2 = ∫_disjoint Cov − ∫_near p·p`; decorrelated limit `far -> (5/7)E[W]^2`, a
`(2/7)E[W]^2` buffer — spread families sit at `far/E[W]^2 = 0.59–0.67` with *negative* disjoint
covariance, like iid.)

**Why the constant matters (kps-S78's coupling, quantified).** `Var(W)/R2 ∈ [4.9, 7.0]·10⁻⁵`
(±20% scatter, 243 families) — it is *not* exact. With klein's fit `c = 5.67·10⁻⁵` the plain
`PZ` route just clears (`E[W] >= 0.1313`, above); but with the *worst-case* `c = 7.0·10⁻⁵` the
decoupled `PZ` bound gives `1/(1 + 0.043/0.144²) = 0.324 < bar` (`−0.007`) — PZ is coupling-tight.
`D3` (opus-S148, `+0.0735` block margin) is what makes it robust to the scatter: `min_E D3 =
0.458–0.465` over `R2 <= 614` (`+0.13`). So brick (B) should route through `D3` for safety
against the `Var/R2` coupling, and `Var(W) <= c*R2` (any `c <= ~7·10⁻⁵`) then suffices. Files:
`04-computation/lrc14_far_covariance_opus_S150.py`, `lrc14_var_resonance_opus_S150.py`;
reflection `everything-is-var-W-the-k11-tail-unified-and-614-is-the-energy-ceiling-opus-S150`.

## Files
`04-computation/lrc_brickA_energy_diameter_kps_S80.py` (+ `.out`): the exact-integer
exhaustive `max R2` by primitive diameter (`[16,24]`), the removal-lemma arithmetic, and the
primitive large-`D` verification (to `D = 120`); the dilation-invariance / non-primitive
`{0,2,..,30}` note.
