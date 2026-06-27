# The Obstruction-Combining Duality: Additive Mediant (Graphs) ↔ Multiplicative I(Ω,2) (Tournaments), Meeting at the Apex 7

*mac-mini-2026-06-22-S58. The owner's observation: the mediant of the edge-densities of the
Kuratowski/Wagner graphs K5, K3,3 reproduces the disjoint unions mK5 ⊔ nK3,3. It has an exact tournament
counterpart, and the two together expose how obstructions COMBINE -- additively on one side,
multiplicatively on the other -- and both pass through the apex 7, where the lonely runner lives.*

## The owner's graph side (verified): obstructions combine ADDITIVELY
For K5 (v,e)=(5,10) and K3,3 (v,e)=(6,9), the disjoint union mK5 ⊔ nK3,3 has (v,e)=(5m+6n, 10m+9n),
edge-density (10m+9n)/(5m+6n). This is precisely the iterated **mediant** of 10/5 and 9/6 -- vector
addition of the (e,v) data = the Stern-Brocot / Farey construction between 3/2 (K3,3) and 2/1 (K5).
(E.g. 6K5 ⊔ 5K3,3 has density 105/60 = 7/4 -- the apex 7 surfacing in the numerator.)

## The tournament side (verified): the SAME obstructions combine MULTIPLICATIVELY
For tournaments, H = I(Ω, 2) (the independence polynomial of the odd-cycle conflict graph at 2; kps), and
I is MULTIPLICATIVE over disjoint union of conflict graphs: I(G1 ⊔ G2, 2) = I(G1,2)·I(G2,2). The
project's forbidden H-values are exactly the K3-conflict combinations:
    I(K3, 2) = 1 + 2·3 = 7,     I(K3 ⊔ K1, 2) = 7·3 = 21,     I(2K3, 2) = 7^2 = 49,  ...
The forbidden set {7, 21} = {I(K3,2), I(K3⊔K1,2)} -- the values whose ONLY conflict-graph preimage is
K3-minimal, and K3 is UNREALIZABLE as Ω (THM-200: three pairwise-sharing triangles force a directed C5).
So H ≠ 7, 21. (49 = 7^2 IS achievable -- via a realizable Ω, not 2K3.) The "combining of obstructions"
is here the disjoint union of conflict graphs, read through the MULTIPLICATIVE I.

## The duality
The single operation -- **disjoint union of the obstruction/conflict graphs** -- has two numerical
shadows:
| | graph data | tournament data |
|---|---|---|
| read | additively: (v,e) sum = **mediant** | multiplicatively: I(·,2) product = **H atom-semigroup** |
| structure | Stern-Brocot / Farey tree | the multiplicative H-spectrum |
| owner's mK5⊔nK3,3 | 7/4 etc. (apex in numerator) | -- |
| forbidden | (non-planar) | {7, 21} = K3-conflict values |
The owner's mediant is the ADDITIVE shadow of ⊔; the H atom-product is the MULTIPLICATIVE shadow. Same
combinatorial act, two arithmetics.

## Where they meet: the apex 7, and the lonely runner
The apex atom is **7 = I(K3, 2)** on the tournament side; 7 is the numerator that surfaces (7/4) in the
Stern-Brocot tree of the graph densities; and the lonely runner lives at **1/14 = 1/(2·7)**, with its
tight locus governed by the **three-gap (Steinhaus) theorem = the Stern-Brocot/Farey structure**
(HYP-2913). So the ADDITIVE (mediant/Stern-Brocot/Farey) thread is exactly the one that controls the LRC
tight locus, while the MULTIPLICATIVE (I/H) thread controls the forbidden values -- and both are pinned
at 7. The owner's observation is the graph-theoretic face of the same apex-7 that runs through H=7
(forbidden), the LRC floor (ζ(2), Farey), and the tight census (three-gap).

## Status (honest)
The pieces are verified: the mediant=disjoint-union (graphs), I(Ω,2) multiplicative with {7,21}=K3-conflict
(tournaments), the Stern-Brocot 7/4, the apex-7 = I(K3,2). The duality (additive mediant ↔ multiplicative
I) and the apex-7 unification are a STRUCTURAL/thematic synthesis -- they explain WHY the same number 7
governs planarity-style obstructions, the H-spectrum, and the LRC, and they identify the additive
Stern-Brocot/Farey thread as the LRC-tight-locus one. This is a unifying lens, not a new step in the LRC
proof (whose open core remains the three-gap rigidity).

Related: HYP-2880 (H multiplicative, {7,21} gaps), HYP-2913 (three-gap census), kps S31g (H=I(Ω,2)),
THM-200 (K3 forces C5), zeta2-governs-the-lonely-runner-floor.
