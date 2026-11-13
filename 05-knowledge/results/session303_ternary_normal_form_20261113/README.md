# The reversible part of the ternary rays

**Proved in note; independently controlled by subset DP.** This does not prove that the 7-ray is the only delayed ray. It identifies the exact reversible triangle mechanism, its finite capacity within each core, and the difference between a local normal form and the missing global descent.

## 1. Disassemble the two holes

The [exclusion proof](../exclusion_7_21_20261103/PROOF.md) has three roles. Insertion injects Hamilton paths into a larger tournament, but does not control the excess. Strong deletion lets the seven-vertex lower bound 25 control every larger strong tournament; the strong spectra at orders 3–6 exclude 7 and 21. Finally, strong components multiply their counts, and any factorization of 7 or 21 needs a factor 7 or 21. Only this last step is freely reversible: an ordinal sum realizes any sequence of available factors. The first attained member of the 7-ray is 63, and every H=63 realization has exactly one nontrivial strong component, since splitting off a triangle or a second nontrivial factor would require 21 or 7.

Write e(d) for the first exponent such that d·3^e is attained, with d odd and 3-free; use infinity if the ray is empty. Appending a directed triangle makes the realized portion a tail. The conjecture is exactly e(7)=2 and e(d)=0 for all other d. The [multiplicative reduction](../session299_three_and_two_holes_20261113/README.md) shrinks the open part to the odd prime and seven-times-prime lanes. Starting new rays is different from traveling along them.

## 2. A genuine confluent exact /3 rewrite

A module is a vertex set seen uniformly by each exterior vertex. Call a directed triangle module *peelable* if the exterior cannot visit it in three runs. Contract it to a root v. In the bipartite graph of arcs from out-neighbors of v to in-neighbors of v, peelability is exactly absence of a matching of size two. Two disjoint excursions each contain such a cross arc; conversely use a 2-matching as two length-two excursions, placing all unused in-neighbors before the first run and all unused out-neighbors after the last. Each leftover set has a Hamilton path. As in [Session 302](../session302_ternary_ray_birth_20261113/README.md), the ordered cover profiles of TT2 and C3 are (1,2,0) and (3,6,6). With exterior skeleton counts a_k, the exact identity is H(C3)=3H(TT2)+6a_3 in that context. Deleting any one vertex of a peelable C3 therefore divides the full count by exactly three. A triangle which is an entire SCC is a special case; its leftover TT2 is ordinal padding.

**Normal-form theorem.** Every tournament has a unique reduced isomorphism type under these induced deletions. In fact a canonical deletion choice within each named triple makes the result literally independent of order.

* Distinct directed triangle modules cannot overlap. An exterior vertex of one triangle cannot see two of its own cyclic neighbors uniformly. If two triangles share just z, the two other vertices of the first see z oppositely, and therefore see the nonshared vertices of the second oppositely; those vertices cannot then see the first triangle uniformly.
* After a deletion, the retained pair is a module. A new triangle module cannot contain both vertices of that pair (the third would see both alike), or just one (the other would see both triangle mates as the retained vertex does, i.e. oppositely). A disjoint triangle module is also a module in the original graph, since the deleted vertex sees its exterior as the two surviving twins do.
* For every other triangle module the matching test remains the same: the shrunk triple is uniform to the tested module and sits on one side of the matching graph, so a matching of size two can use at most two members of it. Swap in the two survivors without changing any cross-incidence. Matching feasibility can only decrease under deletion, so this proves equivalence. Thus no new eligible triangle appears, and no other eligible triangle ceases to be eligible.
* The deletion choices inside a triple are isomorphic relative to its exterior (rotate the original triangle). One may peel all and only the initially peelable triples, in any order. Their simultaneous removal decreases H by exactly 3^t, where t is their number.

Both first-born 63 tournaments are module-reduced. Indeed *every* H=63 tournament is reduced, since a peel would yield the forbidden H=21. Conversely reduced does not mean 3-free: the order-five bit rows (16,1,3,7,14) give a strong module-free H=9 tournament. Likewise (24,17,3,6,12) gives a strong module-free H=15 tournament. A direct labelled head shows no smaller reduced 3-divisible tournament. The distinction is the entire fibre, not any single reduced representative.

## 3. Inverse grammar and an exact resource

The inverse of peeling is precisely completing a uniform TT2 module into a C3 with the matching guard above. At a fixed reduced named core, all inverse expansions are obtained by selecting disjoint eligible pairs; there are no further options created by expansion. Indeed a TT2 module cannot cross a C3 module: a vertex outside C3 cannot be a twin of one of its vertices because the other two see that vertex oppositely. A disjoint TT2 module persists in the reduced graph, and the same 2-matching argument preserves eligibility. These statements prove the claimed inverse grammar by reversing any reduction history.

The eligible-pair graph is a disjoint union of paths. Two overlapping TT2 modules have transitive union with a forced linear order, and extend to a maximal transitive twin run. For detail, two neighbors of a vertex in the pair graph must lie on opposite sides of it in the tournament order, else each twin condition demands the opposite arc between the two neighbors. Thus degrees are at most two. A path of pairs is consistently oriented, and each pair of nonadjacent vertices has the forward arc by repeated use of the module condition. A cycle is impossible: deleting one edge would give an ordered path whose endpoints cannot be twins with an interior vertex present. Every outside vertex sees every member of a run uniformly by transitivity of equality along its pairs.

Eligibility is constant in each run. In the matching graph for a pair contracted to v, earlier twins are in-neighbors of v, later twins are out-neighbors. None can take part in an out-to-in arc: all exterior in-neighbors beat the whole run, all exterior out-neighbors lose to it, and earlier run vertices beat later ones. The matching therefore depends only on the exterior of the run.

A run of length L supplies exactly binom(L−j,j) selections of j disjoint pairs (shift the selected edge indices by 0,1,…,j−1). Thus the polynomial for the marked choices is the product of these path matching polynomials over eligible runs, and the exact capacity is the sum of floor(L/2). These are named expansion choices, not necessarily distinct unlabelled tournaments. At the order-eight H=63 core with apex word 1001100, there are three eligible runs of length two, giving (1+z)^3 and the strong counts 63,189,567,1701. The saturation is structural, not just a failed extension search. Infinite repeats would require changing the reduced core or supplying fresh runs rather than using this rewrite alone. More generally if a fixed module-slot architecture preserves F_2/H, the split-path injection bounds each block order by 1+F_2/H even under adaptive internal replacements.

## 4. What really makes the chain special?

For a realised h, let D(h) be the maximum number of peels possible over all its realisers. Then
```
D(d·3^j)=j−e(d).
```
Every peel divides by 3, so more steps would realise a point below the first. Equality comes by adjoining j−e(d) separate triangles to a first witness. For the 7-ray there are always exactly two factors of 3 unavailable even after optimizing over the fibre. On an ordinary ray starting at d, all factors are removable in at least one representative, but not necessarily the one in hand. This identity is a reformulation, not a shortcut: a free replacement by any member of the same count fibre already presupposes the smaller witness. The nontrivial problem is a constructive, bounded description of count-preserving replacements that reach a peelable core.

An exact hostile test shows that no purely 3-adic substitute for that replacement works. In an apex over a transitive spine, a cyclic triangle must contain the apex and two spine vertices. For it to be a module, they must be adjacent, all earlier vertices must beat the apex, and all later ones must lose to it. Thus the only active word with a C3 module is 10. The two strong, C3-module-free families
```
1 0^a 1100 : H=13·2^a+11    (a≥1)
1^b 0       : H=1+2^b        (b≥2)
```
follow from the inversion formula of Session 302. Modulo every 3^s the first family hits exactly the residues not 2 mod 3, the second exactly those not 1 mod 3, since 2 generates the unit group. Every such residue has arbitrarily large witnesses. Their union hits every 3-adic residue, including unbounded valuation, with only O(log X) exact targets up to X; each exact ternary ray meets it finitely by Session 302. These families are triangle-module-free, not necessarily prime tournaments. This supplies a second kernel for the missing residue class without pretending to fill its prime targets. It is the same distinction that arises in relaxed Collatz: a complete modular address system need not give a positive integer route with the required descent.

The promising next target is a **local neutral replacement preserving a named interface that creates an eligible TT2 or C3 where none existed**, with exact guard and resource accounting. Proving that such replacements cover all prime-indexed rays (or finding a different inverse) is still open.

## Replay

Run `python3 test.py` and `python3 -O test.py`. The checks use subset endpoint DP, literal pair additions, all labelled tournaments of orders 2–5, all eight marked 63 expansions, and the order-five reduced witnesses. Finite checks are not used to establish all-order confluence.
