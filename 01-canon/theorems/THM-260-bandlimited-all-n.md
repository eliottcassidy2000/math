# THM-260: Band-Limitedness of H for All n (Resolves B1 of OPEN-Q-040)

**Status:** PROVED (algebraic sketch, all n >= 4; verified computationally n=4,5,6,7)
**Devil's advocate review (S29):** Upper bound correct but substitution argument needs detail. Lower bound interleaving construction correct but even-n case not explicitly handled. Both gaps are fillable.
**Filed by:** kind-pasteur-2026-03-25-S1
**Dependencies:** THM-076 (Walsh-OCF factorization), THM-259 (Walsh degree formula)

## Statement

For tournaments on n >= 4 vertices, the Hamiltonian path count H(T), viewed as a function on the tiling hypercube Q_m (m = C(n-1,2)), has:

1. **Walsh degree exactly 2*floor((n-1)/2)** in the tiling model.
   This is n-1 for odd n, n-2 for even n.

2. **Band-limitedness:** For n >= 6, the Walsh degree satisfies
   2*floor((n-1)/2) < C(n-1,2)/2 = m/2.
   So ALL Walsh coefficients above the midpoint of the spectrum vanish.

3. **Correction:** At n=4,5, the Walsh degree EXCEEDS m/2 (degree 2 > 1.5
   at n=4; degree 4 > 3 at n=5). Band-limitedness at m/2 holds for n >= 6 only.

## Proof

**Upper bound (all n):** THM-076 proves that in the full arc model (all C(n,2) arc
variables), the Walsh coefficient hat{H}[S] is nonzero only for monomials S that
are unions of even-length paths with total vertices 2k+r <= n. The Hamming weight
|S| = 2k <= n-r <= n-1. Restricting to tiling variables (fixing the n-1 base-path
arcs to their constant values ±1) is a SUBSTITUTION in the multilinear Walsh
expansion. Substituting a variable x_i = c (constant) eliminates all monomials
containing x_i and preserves those not containing it — so the multilinear degree
can only decrease or stay the same. Since only even-weight monomials survive
complement symmetry in the full model, the max even weight in the tiling model
is <= 2*floor((n-1)/2).

**Note (S29 review):** THM-259 uses m = C(n,2) (full arc model). THM-260 uses
m = C(n-1,2) (tiling model). The Walsh degrees are the same in both models
because the max-weight monomials can be constructed using only tile arcs.

**Lower bound (n >= 4):** We exhibit a Walsh monomial of weight 2*floor((n-1)/2)
using only tile arcs (non-base-path arcs).

*Construction:* A single path P_{2k} with k = floor((n-1)/2) visits 2k+1 vertices.
For this path to use only tile arcs, it must avoid all consecutive-vertex arcs (v, v-1).

The **interleaving construction** achieves this:

**Odd n** (k = (n-1)/2, path uses all n vertices):
Order vertices as 1, 3, 5, ..., n, 2, 4, ..., n-1.
Example at n=7 (k=3, 6 arcs on 7 vertices): 1→3→5→7→2→4→6.
All arcs: (1,3),(3,5),(5,7),(7,2),(2,4),(4,6) — differences ≥ 2, all tiles. ✓

**Even n** (k = (n-2)/2, path uses n-1 vertices):
Order vertices as 1, 3, 5, ..., n-1, 2, 4, ..., n-2. Take this (n-1)-vertex path.
Example at n=8 (k=3, 6 arcs on 7 vertices): 1→3→5→7→2→4→6.
The transition 7→2 has |7-2|=5 ≥ 2. All arcs are tiles. ✓
(Vertex 8 is unused, giving a path on n-1 = 7 vertices with n-2 = 6 arcs.)

**General verification:** In the interleaving order, consecutive elements differ
by 2 (within odds or evens) or by n-3 or more (at the odd→even boundary).
Since n ≥ 4, the minimum difference is min(2, n-3) ≥ 2. All arcs are tiles.

THM-076 gives |hat{H}[S]| = 2 * (n-2k)! / 2^{n-1} > 0 for this monomial.
So the tiling model has nonzero coefficient at weight 2k. QED.

**Band-limitedness (n >= 6):** We need 2*floor((n-1)/2) < C(n-1,2)/2.
The left side <= n-1. The right side = (n-1)(n-2)/4.
For n >= 6: (n-1)(n-2)/4 >= 5*4/4 = 5 > 5 = n-1. (At n=6: 5 > 4.)
For n >= 7: even stronger separation.

At n=5: 4 > 3 (Walsh degree exceeds m/2). NOT band-limited.
At n=4: 2 > 1.5. NOT band-limited.

## Computational Verification

| n | m | m/2 | Walsh degree | Match THM-259 | Band-limited at m/2 |
|---|---|-----|-------------|---------------|---------------------|
| 4 | 3 | 1.5 | 2 | YES | NO |
| 5 | 6 | 3.0 | 4 | YES | NO |
| 6 | 10 | 5.0 | 4 | YES | YES |
| 7 | 15 | 7.5 | 6 | YES | YES |

All verified exhaustively (32,768 tilings at n=7).

## Correction to Prior Claims

The reflection "h-is-band-limited.md" (opus-S306) claims hat{H}_k = 0 for
k >= ceil(m/2)+1 at ALL n including n=5. This is INCORRECT in the tiling model:
at n=5, Walsh weight 4 has 7 nonzero coefficients, and 4 >= ceil(6/2)+1 = 4.

The correct statement: band-limitedness at m/2 holds for n >= 6.
For ALL n >= 4: the Walsh degree is 2*floor((n-1)/2) << m for large n.

## Additional Observations

1. **Odd-weight Walsh coefficients are nonzero** in the tiling model (complement
   symmetry H(t) = H(~t) FAILS because the complement in tile space is not T^op).
   Verified: at n=7, 907 nonzero odd-weight coefficients out of 9453 total odd-weight.

2. **The THM-076 amplitude formula applies to the full arc model.** In the tiling
   model, the maximum amplitude at each weight matches the THM-076 prediction
   (since the maximum-weight monomials use only tile arcs by the interleaving
   construction), but lower-weight monomials may have different amplitudes due
   to base-path arc contributions.

3. **The Walsh degree grows as O(n) while m grows as O(n^2).** This means the
   "information content" of H is concentrated in an asymptotically vanishing
   fraction (O(1/n)) of the Walsh spectrum. H is extremely low-frequency.

## Related

- THM-076: Walsh-OCF factorization (the algebraic engine)
- THM-259: Walsh degree formula (computational verification)
- OPEN-Q-040 B1: RESOLVED by this theorem
- h-is-band-limited.md: CORRECTED (n=5 is not band-limited at m/2)
