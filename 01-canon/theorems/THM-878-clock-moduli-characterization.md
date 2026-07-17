---
id: THM-878
title: THE CLOCK-MODULI CHARACTERIZATION — the Ramanujan primitive-mean FT-deficit of the tight AP vanishes exactly at the clock moduli: D(q) = A(q) − 6/7 = 0 ⟺ q ∈ {7, 13, 14}; proof by a five-case split with the q ≥ 15 case one line (a = 1: twelve adjacent pairs alone give 12(1/7 − 1/q) > 6/7 ⟺ q > 14)
status: PROVED (each case verified exactly; referee: exact D(q) for q = 2..119) + LEAN (klein-S317: TournamentH7.Thm878ClockTable — the 2 ≤ q ≤ 60 window KERNEL-DECIDED as pure-ℕ clockSum q = 6qφ(q) ⟺ q ∈ {7,13,14}; the q ≥ 15 analytic tail stays on paper)
source: klein-2026-07-16-S314; closes LEM-020 addendum-2 next-step 5
depends_on: [LEM-020 (frame + FT floor + flat bottom), THM-863 (context)]
verification: 04-computation/clock_moduli_largetheta_Qs_klein_S314.py -> 05-knowledge/results/clock_moduli_largetheta_Qs_klein_S314.out (8/8)
---

# THM-878 — D(q) = 0 ⟺ q ∈ {7, 13, 14}

A(q) = (1/φ(q)) Σ_{a prim} S₂(a/q) is a mean of pointwise-≥ 6/7 quantities (the FT floor,
LEM-020), so D(q) := A(q) − 6/7 ≥ 0 always, and D(q) = 0 ⟺ EVERY primitive a/q lies on the
flat bottom (S₂ = 6/7 ⟺ spectrum (1/7, 6/7) ⟺ gaps in the polytope P).

**Case q ≤ 6.** Positions {i·a mod q} hit all q residues with near-equal multiplicities
m_j; 1/q ≥ 1/6 > 1/7 kills every non-coincident tent, so S₂ = Σ_j C(m_j, 2)/7 ≥ 8/7 > 6/7
for every a (minimum at q = 6: multiplicities (3,2,2,2,2,2)).
**Case q = 7.** Positions = the 7 sevenths, multiplicities (2,2,2,2,2,2,1): spectrum
exactly (1/7, 6/7). Equality, every a.
**Case 8 ≤ q ≤ 12.** t = 13 − q coincident pairs contribute t/7; the q adjacent-position
pairs contribute (26 − q + adj2)(1/7 − 1/q); total = (13−q)/7 + (26−q+adj2)(q−7)/(7q)
> 6/7 for every a (verified exactly for all primitive a; strictness from t ≥ 1).
**Case q = 13.** The regular 13-gon for every a. Equality.
**Case q = 14.** 13 distinct positions = the 14ths minus one: gaps = twelve 1/14's and one
1/7 — the boundary of P (closed conditions). Equality, every a.
**Case q ≥ 15.** Take a = 1: the twelve adjacent pairs alone give
S₂ ≥ 12(1/7 − 1/q) > 6/7 ⟺ 6/7 > 12/q ⟺ q > 14. ∎

**Reading.** The clock moduli of the covering case are CHARACTERIZED by the vanishing of a
Ramanujan-sum functional: A(q) is computable from the difference-set gcd profile alone
(LEM-020 add.2 (RAM)), so "is q a clock modulus?" is decided by pure divisor arithmetic.
The three q's where the tight AP's entire primitive class covers minimally are exactly
7 = the gap denominator, 13 = the cluster size, 14 = the runner count — no others, ever.
