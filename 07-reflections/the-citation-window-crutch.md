# The citation-window crutch (an assumption audit that closed four strata)

**opus-2026-07-05-S91** (HYP-4196). The owner asked: reframe your assumptions; where do
they point at missing structure? This is the record of the audit and what fell out.

## The assumption

Every window argument in the l >= 7 program (S76 onward) had the shape: cite LRC(|B|+1)
on the unlifted base for margin beta = 1/(|B|+1) at some t*, transfer by Lipschitz to a
window of half-width delta = (beta - 1/13)/12, then dodge everything else inside the
window. The window width delta was load-bearing in every co-incidence table -- and it
was the reason the B <= 4 strata saturated (42+ zones of irreducible width 2*delta
cover the circle) and the degenerate types looked hard.

## The reframe

The final point of the whole proof needs every runner at distance >= 1/13. The base
needs margin 1/13 + eps AT THAT POINT -- not 1/6, not anywhere else. The citation and
its window were only a SEED: a guarantee that a good point exists. But for <= 5 base
runners the union bound is VALID at band 1/13 (measure of good >= 1 - 2|B|/13 >= 3/13):
existence is free, fat, and structured. The margin-1/6 citation was a crutch carried
from the peel lemmas into a setting that never needed it.

## What closed

Working directly with the components of good(B, 1/13 + eps) minus the exact desert
zones (no window term): every stratum |B| = 2..5 keeps clear components of length
>= 0.0111 in the worst case OVER-GRANTING all degeneracy types simultaneously; the
descent entry bar is 2/component <= 179 everywhere; the bounded box [14, 179] covers
below it. THE CO-INCIDENCE VARIETY IS CLOSED AT EVERY STRATUM -- the B <= 4 saturation
was an artifact of the crutch, not a feature of the problem.

## The pattern worth keeping

Three of the last four corrections in this program have the same shape:
* the fee criterion (S81): a location-oblivious budget where adaptivity was available;
* the Lipschitz window (S91): a worst-case transfer where the true good set was fat;
* the naive incidence table (S87): worst-case radii where per-cluster refinement existed.

Each time, the assumption imported a UNIFORM bound (fees, windows, radii) into a place
where the STRUCTURED object (leftover sets, good components, per-type zones) was
directly available and strictly better. The missing insight was never a new theorem --
it was that the crude bound had quietly become an architectural assumption. The
remaining fleet frontier (the loose branch's second-value induction) should be audited
for the same pattern: klein's drag windows and kps's peeling windows both carry
Lipschitz seeds; the band-level good sets there are fat too (12 runners at 2/25:
union bound invalid -- but at intermediate bands and for sub-families it is not).

## Convergence note

kps-S15 (HYP-4177) independently reached the multi-leg domination law (deep well
14/169 = the global fourteen-fold minimum) hours after S90 -- two proofs, one law;
cross-referenced. The fleet's convergence rate on this endgame is itself a signal:
the remaining surface is small.
