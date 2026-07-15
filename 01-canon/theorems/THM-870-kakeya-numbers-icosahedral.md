---
id: THM-870
title: KAKEYA NUMBERS OF THE ICOSAHEDRAL GROUPS — needle directions = the 31 maximal cyclic subgroups (= the icosahedral axes 15+10+6); K(A₅) = 15 EXACTLY (lower bound from the six 5-fold needles alone: 6·5 − C(6,2) = 15, matched by an explicit 15-element witness hitting all 31 axis-cosets); every maximal cyclic of 2I contains −1, so the 2I problem is the exact ±-pullback: K(2I) = 2·K(A₅) = 30; the ODD (tournament-hostable) shadow already costs the full 15 — the 15 involution axes are FREE (Feit–Thompson face)
status: PROVED (lower bounds: trivial-intersection Bonferroni; upper bounds: explicit verified witnesses; the ± descent: every maximal cyclic of 2I contains the unique involution −1)
source: klein-2026-07-15-S313 (cont.3); owner directive "combine icosian ring with Kakeya needle ideas"
depends_on: [THM-868 (icosian/E8 face), classical: 2I = SL(2,5), icosian ring]
related: [THM-856 (LRC Kakeya-rank formalism — the direction-incidence analogy), the-solvability-tower s703]
verification: 04-computation/icosian_kakeya_milgram_sedenion_klein_S313.py -> 05-knowledge/results/icosian_kakeya_milgram_sedenion_klein_S313.out (24/24)
---

# THM-870 — Kakeya numbers of A₅ and 2I

**Definition.** In a finite group G, a *needle direction* is a maximal cyclic
subgroup H; a *needle in direction H* is a left coset gH; a **Kakeya set** is an
S ⊆ G containing a needle in every direction; K(G) = min |S|. (The group model of
"a segment in every direction"; for A₅ the directions are literally the 31
icosahedral rotation axes: 15 two-fold, 10 three-fold, 6 five-fold.)

## 1. K(A₅) = 15

Lower bound: distinct maximal cyclics intersect trivially, so two needles in
different directions share ≤ 1 point. The six 5-fold needles alone force
|S| ≥ 6·5 − C(6,2)·1 = 15 (Bonferroni). Upper bound: an explicit 15-element set
(stored in the .out) contains a full coset of **all 31** directions. ∎

The optimum is forced by the 5-fold needles; the 3-fold and 2-fold needles ride
along free. **Odd shadow:** restricting to odd-order directions (6 Z₅ + 10 Z₃ — the
only ones a tournament automorphism group can host, by Feit–Thompson) the minimum
is still 15: the even needles cost nothing. The obstruction content of the Kakeya
problem in A₅ is entirely odd/solvable-visible, even though A₅ itself is not.

## 2. K(2I) = 30 (exact ±-descent)

Every maximal cyclic of 2I (15 Z₄, 10 Z₆, 6 Z₁₀) contains the UNIQUE involution −1;
hence every needle is ±-closed and equals the full preimage of an A₅ needle under
2I → 2I/{±1} = A₅. So Kakeya sets of 2I are exactly preimages of Kakeya sets of A₅
augmented arbitrarily, and K(2I) = 2·K(A₅) = 30; the lifted witness verifies. ∎
(The annealer found 32 before the descent was noticed — the theorem beat the search.)

## 3. Readings

- 15 = the number of 2-fold axes = the edge-pairs of the icosahedron: the Kakeya
  budget of the icosahedral group equals its involution count — but the involutions
  are not what costs; the 5-fold needles are. (Numerology guard: 15 also = C(6,2)
  from the bound; the coincidence is the bound's arithmetic, 30 − 15.)
- In the icosian ring the 6 five-fold directions are the 6 golden axes of the
  600-cell; a minimal Kakeya set is a 30-icosian configuration hitting every
  rotation axis — a finite model of "needle turning through all directions" inside
  the same 2I whose McKay correspondent is THM-868's E8.
- LRC analogy (THM-856): "uncovered mass = deficit + rank of the needle incidence"
  — here the incidence matrix of (chosen cosets × directions) has the trivial-
  intersection property that makes the first-moment (Bonferroni) bound TIGHT; the
  icosahedral case is the clean solvable-visible model of a tight needle system.
