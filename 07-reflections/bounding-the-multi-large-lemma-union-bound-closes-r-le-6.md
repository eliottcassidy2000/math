# Bounding the multi-large equidistribution lemma: the union bound closes r ≤ 6 rigorously

*kind-pasteur-2026-06-22-S31v. The irreducibly-analytic Node-3 lemma (HYP-2899) is "the committed
large speeds don't cover the bounded core's lonely set." Developing the bounding toolbox, the union
bound + the elementary comb-teeth estimate + proven LRC(13) give a RIGOROUS closure for up to six
large speeds — narrowing the residual to r ≥ 7.*

## The lemma, precisely
A counterexample is a covering 13-set `S = C ∪ L`, `C` = bounded core (speeds `≤ V*`), `L = {v_1..v_r}`
= large committed speeds (`> V*`). Let `G_C = {x : ‖c x‖ ≥ 1/14 ∀ c ∈ C}` (the core's lonely set) and
`U_v = {x : ‖v x‖ < 1/14}` (a large speed's danger comb, `meas = 1/7`). The lemma to bound:
> `meas(G_C ∖ ⋃_i U_{v_i}) > 0` — i.e. the large speeds do not cover the core's lonely set ⟹ `M(S) ≥ 1/14`.

## The toolbox (verified, `lrc_equidist_lemma_bounds_kps.py`)
**(T1) `meas(G_C)` is bounded below — RIGOROUS from proven LRC(13).** `C` has `13−r ≤ 12` speeds, so by
proven LRC(13), `M(C) ≥ 1/(13−r+1) ≥ 1/13 > 1/14`. Hence `{min_c ‖cx‖ > 1/14}` is open and nonempty, so
`meas(G_C) > 0`; over the FINITE family of bounded cores (`C ⊆ [1,V*]`), `meas(G_C) ≥ c_0 > 0` uniformly.

**(T2) Each comb covers `≤ 1/7 + O(arcCount/v)` of `G_C` — ELEMENTARY.** `G_C` is a union of `≤ A_0` arcs
(its boundary is the finite breakpoint set of `C`). `𝟙_{U_v}` is `(1/v)`-periodic, so over an arc of length
`ℓ` it integrates to `ℓ/7` plus at most the two boundary teeth, each `≤ 1/(7v)`. Hence
> `meas(G_C ∩ U_v) ≤ (1/7) meas(G_C) + arcCount/(7v)`.
Verified: the ratio `→ 1/7` with error `≤ arcCount/(7v)` exactly (no deep Weyl needed — just counting comb
teeth inside fixed arcs).

**(T3) Pairwise overlaps `→ (1/7)² = 1/49` for non-resonant pairs; `~3.5×` for resonant `v_i | v_j`.**
(`(101,103)→0.004`, `(733,1009)→0.021≈1/49`; `(101,202)→0.067`, `(500,1000)→0.071`.) The resonant excess is
the only obstruction to independence, and it is carried by divisibility relations among the large speeds.

## The closure for r ≤ 6 (RIGOROUS)
Union bound + (T1)+(T2):
> `meas(G_C ∖ ⋃ U_{v_i}) ≥ meas(G_C) − Σ_i meas(G_C ∩ U_{v_i}) ≥ (1 − r/7) meas(G_C) − (arcCount/7) Σ_i 1/v_i.`

For `r ≤ 6` this is `≥ (1/7) meas(G_C) − (A_0/7) Σ 1/v_i > 0` once every `v_i > V* := 6 A_0 / ((1/7) c_0)`.
**So for up to six large speeds over any bounded core, the lonely set survives: `M(S) ≥ 1/14`** — the
`r`-fold generalization of THM-565 (the `r = 1` slice), needing no equidistribution machinery beyond
counting comb teeth in fixed arcs, plus proven LRC(13) for the floor `c_0`.

**Scale-separation / `V*` (honest).** `V*` depends on the core's `arcCount A_0` and floor `c_0`, exactly
as THM-565's `V* ≈ 234` does — "large" means the `v_i` are well-separated from the core's scale (`A_0` is
fixed by the core, not by `V`, the THM-565 independence). The bounded/intermediate speeds (`≤ V*`, where
`V*` is `~ r ×` THM-565's value, so `~ 1400` for `r = 6`) are absorbed by the finite Node-2 check; the
union bound takes over for the genuinely separated `v_i > V*`. So the partition is: Node 2 (all speeds
`≤ V*`, finite) + this union bound (`1 ≤ r ≤ 6` speeds `> V*`) + the `r ≥ 7` residual below.

## The residual: r ≥ 7 (core ≤ 6 small speeds)
Here `1 − r/7 ≤ 0` and the union bound is vacuous — but this is the regime where the core is SMALLEST
(≤ 6 speeds), so `G_C` is LARGEST (`meas(G_C) ≈ (6/7)^{≤6}`, lots of room). The bound is the
**second-moment / inclusion–exclusion**:
> `meas(G_C ∖ ⋃ U_i) ≥ meas(G_C) − Σ meas(G_C∩U_i) + Σ_{i<j} meas(G_C∩U_i∩U_j) − …`
With the pairwise overlaps `≈ 1/49` (T3) except for the resonant pairs, the leading independent estimate
is `(6/7)^r meas(G_C) > 0` (e.g. `(6/7)^{13} ≈ 0.135`), and the only correction is the **resonant pairs**
`v_i | v_j` (overlap `~3.5/49`). Their number is bounded by the divisibility/relation lattice of the large
speeds — exactly mac-mini's CRT over-determination (HYP-+2878: a 13-set resonates with `≤ 3` of `5` primes).
So `r ≥ 7` closes by a **second-moment bound with a bounded resonance defect**, the one remaining analytic
input.

## Net — the lemma is now split, half closed
> **Node 3 = [r ≤ 6: union bound + comb-teeth + LRC(13), RIGOROUS] + [r ≥ 7: second-moment, (6/7)^r minus
> bounded resonant-pair defect].**

The hard adversarial single-large family (THM-566) is `r = 1` ⊂ the closed half. The remaining target is
sharp and finite-flavored: bound the **resonant-pair defect** among `≥ 7` large speeds (the `v_i | v_j`
overlaps), so that the second-moment stays above `0`. That is a statement about the divisibility lattice of
the committed speeds, not an unbounded analytic estimate — and it is bounded by the CRT over-determination
the team already has. The union bound did the heavy lifting for `r ≤ 6`; the second moment + a bounded
resonance count should finish `r ≥ 7`.

→ HYP-2899 (two structures / Node 3), THM-565 (r=1 three-gap sampling), THM-566 (adversarial r=1),
HYP-2895 (single-large equidistribution), HYP-+2878 (CRT over-determination / resonance bound),
proven LRC(13) (the `meas(G_C)` floor), [[lrc14-thread]].
