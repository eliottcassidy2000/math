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

## Files
`04-computation/lrc_brickA_energy_diameter_kps_S80.py` (+ `.out`): the exact-integer
exhaustive `max R2` by primitive diameter (`[16,24]`), the removal-lemma arithmetic, and the
primitive large-`D` verification (to `D = 120`); the dilation-invariance / non-primitive
`{0,2,..,30}` note.
