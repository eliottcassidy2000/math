# The RM duality square: level permutations and the d⁺ fork (claudebox-2026-06-11-S2)

User dispatch: RM(r,m)⊥ = RM(m−r−1,m); repetition↔SPC, biorthogonal↔extended Hamming,
k = n/2 self-dual — "see if you can find equations or recurrences ... keep encircling the
truth of the underlying structure," with pointers at the tiling metagraph, blue line pairs,
and the most symmetric tournaments.

## What the dispatch unfolded into

**1. The switching world is an RM-duality square with an odd twist (THM-479).** The graph
side is 𝔽₂-linear: cut space ⊥ cycle space inside the edge space is K_n's avatar of
RM(r)⊥ = RM(m−r−1) (the repo already had this as base-path/wiggly arcs). Mallows–Sloane
(#two-graphs = #Euler graphs) is its Burnside shadow. Tournaments are the AFFINE TORSOR
over that square: the cocycle t_π is the odd twist (LEM-004's "tournaments are odd
functions" at the group-action level), and π fixes a switching class iff a 2-adic
obstruction vanishes. Deriving the obstruction from scratch: the holonomy of t_π is
supported exactly on the DIAMETER pair-orbits {i, π^{ℓ/2}i} of even cycles — antipodal /
half-turn / twin pairs, the σ-locus the repo keeps meeting (THM-430, THM-447 twins, the
blue-line hypotenuse) — and cut-compensation works iff all cycle lengths share one
2-adic valuation. That law turned out to be Babai–Cameron's "level" Lemma 7.1 (proved
1981-82, lost, published 2000) — an honest independent rediscovery with a more elementary
proof, verified exhaustively (n ≤ 10, all cycle types) and against all 16 OEIS terms of
A049313. NEW: the branch split N_odd + N_lev (separately integral for n ≥ 3, N_lev = 0 at
odd n, neither in OEIS), which sharpens OPEN-Q-060: the "odd Mallows–Sloane partner" must
itself split along the 2-adic seam.

**2. The two doublings generate the two ends of the RM chain (THM-480).** Sylvester
doubling = Plotkin (u, u+v) = the biorthogonal end: row code RM(1,m), dim m+1. The skew
tower (THM-447) in the TOURNAMENT GAUGE (M = I+S, no renormalization) = pair-doubling +
glue = the SELF-DUAL MIDDLE: k = n/2 Type II at every level — ê₈ = RM(1,3) (the r = 1
self-dual RM point) at order 8, then d₁₆⁺ at 16, d₃₂⁺-enumerator at 32. The user's
"k = n/2 self-dual ↔ most symmetric tournaments" is literally what skewness does to the
code: the diagonal observer +1 and the twin corrections pump the row space from the bottom
of the RM filtration to its self-dual fixed point. At 16, RM has no self-dual member, the
tower must leave the family, and the exit is THE fork: the two Type II [16,8] codes
e₈⊕e₈ / d₁₆⁺ share one weight enumerator (Gleason/Milnor isospectral pair; E₈×E₈ vs SO(32))
— and the tower picks d₁₆⁺, which EXPLAINS kps1's Hall-class pin (tower-16 ≅ had.16.3,
excluded had.16.4 only by a 10⁷-node witness search; the weight-4 support graph's
connectivity separates them in milliseconds). THM-451's Hadamard split at 16 is, at code
level, the isospectral fork.

**3. Convergence.** mac-mini-2026-06-11-S1 took the same dispatch into TILING space
(THM-477 blue code: ker(1+σ), glue 𝔽₂^f on the hypotenuse; THM-478, HYP-2406..2408).
This session took it into AMBIENT and GROUP-ACTION space. The two self-dual defects sit on
the same locus (hypotenuse = twin anti-diagonal = diameter orbits): HYP-2409(3) asks for
the exact glue map. Three sessions, one dispatch, three spaces — the blue line pairs were
the right pointer: they are the σ-fixed/antipodal locus in every incarnation.

## Honest scope

- THM-479 parts 1-2 are Babai–Cameron 1981/2000 rediscovered (cited); the holonomy proof,
  the branch split, branch integrality (unexplained), and the OPEN-Q-060 sharpening are ours.
- THM-480's order-32 identification is enumerator-level (Type II + d₃₂⁺ enumerator), not a
  certified equivalence; order 16 is rigorous (two-code classification + indecomposability).
- Equations delivered: the level law; A049313 = Σ_odd 2^{o₂−k+1}/z + Σ_lev 2^{o₂−k}/z;
  o₂ = Σ⌊ℓ/2⌋ + Σgcd; dim(C∩Im) = n−k−[no odd cycle]; the d⁺ pair-doubling mechanism;
  C(H₈) ≅ RM(1,3) with explicit witness (0,1,2,4,6,3,7,5).
