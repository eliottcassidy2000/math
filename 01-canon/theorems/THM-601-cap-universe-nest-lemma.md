---
id: THM-601
title: The Cap-Universe Nest Lemma — every d ≥ 3 danger intersection within {1..13} at r = 1/14 is gcd-nest-exact (|∩_Q D_v| = 2r/max(Q/gcd Q)); hence Λ_P(1/14) is closed-form rational (singles − THM-594(B) pairs + nest chain) for every P, and the cap ladder, THM-530's m_P, and the pentagon census minimum are one object
status: PROVED (complete exact-rational verification over the finite universe: all 8100 subsets, zero violations; mechanism lemma — the mediant criterion — proved for the boundary and verified on all 79 violations outside the universe)
source: klein-2026-07-01-S92 (canonizing HYP-3837/HYP-3838, klein-S91)
scripts: 04-computation/exact_cap_ladder_decomposition_klein.py (+ .out)
related: THM-594(B), THM-576, THM-530, THM-596, HYP-3837, HYP-3838, HYP-3848, kps HYP-3954 (the c-averaged/torus complement)
---

# THM-601: the Cap-Universe Nest Lemma

## Statement

Let r = 1/14 and, for a finite speed set Q ⊂ Z⁺, let D_v = {t ∈ R/Z : ||vt|| < r} and
g = gcd(Q). Then for every Q ⊆ {1, ..., 13} with |Q| ≥ 3:

    | ∩_{v ∈ Q} D_v |  =  2r / max(Q/g)      (the gcd-reduced origin nest).

Consequently, for every P ⊆ {1..13}, by inclusion–exclusion truncating exactly at depth 2:

    Λ_P(1/14) = 1 − 2r·|P| + Σ_{pairs} |D_u ∩ D_v| − NestChain(P),

where the pair terms are THM-594(B)'s exact two-branch values and NestChain(P) =
Σ_{|T| ≥ 3} (−1)^{|T|} · 2r/max(T/gcd T) is elementary (for gcd-free P it telescopes; see
the NestTelescope Lean lemma). Every cap constant is therefore derived rational
arithmetic; in particular the ladder rungs (HYP-3837) cap_10 = 14249/252252 = THM-530's
m_P and cap_11 = 313/9702 = the pentagon 11-core census minimum.

## Proof

**(1) Reduction and lower bound.** t ↦ gt preserves Lebesgue measure on R/Z and maps
∩_{v∈Q} D_v onto ∩_{v∈Q/g} D_v, so assume gcd(Q) = 1. The interval of the fastest speed
w = max Q around 0 has half-width r/w, minimal among the 0-intervals of all v ∈ Q, so the
origin nest contributes exactly 2r/w to the intersection.

**(2) Any excess is a mediant event.** A point x of the intersection outside the origin
nest satisfies ||vx|| < r for all v ∈ Q with x bounded away from 0 (mod 1). Fix any two
speeds v < w of Q: x lies in D_v ∩ D_w away from 0, i.e., in a wrapped coincidence
component of the pair — an interval centred at a point x_c of the pair's coincidence
family, where fractions a/v and b/w with d = |bv − aw| ≥ 1 approach within r/v + r/w;
these x_c are the mediant-family points ((a+b)/(v+w)-type for the first layer). For the
excess to survive in the full intersection, x must also lie in D_u for the remaining
speeds u — i.e., **a pair-coincidence point must fall inside a third speed's danger**
(the mediant criterion).

**(3) The universe is below every threshold.** Coincidences at pair (v, w) exist away
from 0 only when r(v + w) is large enough for the d ≥ 1 approach, and their positions
x_c must additionally fall within r/u-width of a u-fraction. Onset examples (all 79
violations among triples ≤ 20 verified to carry the mechanism): the family (1, v, v+1)
wraps iff its first mediant 2/(2v+1) < r ⟺ v ≥ 14 — onset (1, 14, 15) via 2/29, the
same mediant governing the final window (THM-596, MISTAKE-093); (2, 13, 15) wraps via
the (13,15) coincidence at ≈ 0.4948 inside D_2's arc at 1/2. Within {1..13}, all pair
sums are ≤ 25, below the observed onset thresholds.

**(4) Completion by exhaustive exact verification.** The claim is a finite statement:
for each of the 8100 subsets Q ⊆ {1..13}, |Q| ≥ 3, the intersection measure is an exact
rational computed by interval arithmetic over Q(⊂ ℚ)-endpoints. The verification
(script above, exact rationals throughout) reports **zero violations**. ∎

*Remark (scope).* The lemma is a PINNED (common-start) statement: the origin nest is the
common zero's artifact. Under free phases (shifted LRC) the nest disappears; the phased
analogue is THM-598's territory and kps's torus-band theorem gives the phase-averaged
d-folds. Above the universe (speeds > 13 or other radii) the wrapped corrections are the
Bernoulli/subtorus slices (mac-mini S96 §1; kps HYP-3954).

## Consequences

1. **hp0cap's cap side is closed-form.** No Vitali estimate, no measure theory: the
   sector/cap arithmetic of THM-534/576 reduces to THM-594(B) pairs + the nest chain.
2. **One ladder.** The census minimum (kps-S27), the pigeonhole floor m_P (THM-530), and
   the caps (THM-576) are rungs of the single function j ↦ min_{|P|=j} Λ_P(1/14),
   computed exactly for all j = 1..13 in HYP-3837.
3. **The L_y extremality comparison** ("consecutive maximizes", THM-534's open half)
   becomes a finite rational comparison with all values in closed form.
