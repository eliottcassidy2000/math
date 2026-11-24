# The PG(2,4) branch needs dualities and cannot have cyclic signs

**PROVED conditional filters.** The primitive degree-56 catalogue leaves the alternating and symmetric actions and the Gewirtz actions after the elementary dimension tests. In the Gewirtz branch, use the realization of sheets as one PSL(3,4) orbit of 56 hyperovals in PG(2,4); its automorphism group acts by collineations or correlations. A hyperoval is a six-point set meeting every line in zero or two points, and its six exterior lines form a dual hyperoval. See the [hyperoval model in the Gewirtz atlas](https://aeb.win.tue.nl/graphs/Sims-Gewirtz.html). These restrictions apply to transport, not global graph automorphisms.

## Incidence denominator and its signed refinement

Let M be the 21 by 56 incidence of projective points in the hyperovals, and N the square 21 by 21 point-line incidence. Each point is on 16 hyperovals and each pair on 4, by double transitivity of PSL(3,4) on projective points and 56 double counts. Thus MM^t=12I+4J, so M has full row rank. Also NN^t=4I+J with all row and column sums five; exterior-line incidence is E=J−N^t M/2.

In the hypothetical mapping sum of 56 sharply transitive permutations on sheets (e.g. direct and two-hop maps multiplied by the inverse of the direct match), write c,d for the numbers of collineations and correlations, respectively, and C,D for the corresponding sums of point-to-point and exterior-line-to-point permutation matrices. Sum the column images under M. The defining identity on the 56 sheets is sum(P)=J_56, so

    C M + D E = M J_56 = 16 J,
    6 C + d J − 3 D N^t = 16 J,

where the second identity uses J M=6J and full row rank. Multiplying by (N^t)^(-1)=N/4−J/20 and using C J=c J gives

    D = C N/2 + (16−c)J/6.                                           (4)

In particular c≡d≡1 mod 3 and C N≡c J mod 2. An all-collineation family is impossible since 56 is not 1 mod 3; directly it would require C=(8/3)J. This excludes the PSL(3,4) action and every extension acting only by collineations.

Neither c nor d can be 1 or 4. For c=1 the parity relation requires the 0/1 row weights of a permuted N to be all ones, but they are five. For c=4, in each row the points of odd multiplicity in C form an even set (meeting every line evenly). A nonempty even set must have at least six points: five lines through a point of the set need one other point each. Therefore the four images of every point have even multiplicity, forcing at least two pairwise coincidences per point, or at least 42 over all six pairs of collineations. Two distinct collineations agree at at most seven projective points (six if their relative map is linear). To see this, four fixed points forming a projective frame conjugate the map to a field automorphism, fixing seven points if it is nontrivial. Without a frame, a set containing a noncollinear triple is confined to one side of that triangle plus the opposite vertex: adding nonvertices on two different sides gives a frame. This bounds its size by six. Among four collineations at least two of the pairs have the same field part, so the upper bound on coincidences is at most 2·6+4·7=40, a contradiction. To exclude d=1,4, if correlations occur, compose every member of the sheet set with one fixed correlation; this preserves sharp transitivity and interchanges c,d. Thus d is one of 7,10,…,49 in every sharp family.

## All endpoint pairs at once: a ternary sign test

Let S be the symmetric 57 by 57 matrix with diagonal zero and S_ij=+1 for collineation and −1 for correlation. The sharp family for pair i,j, made relative to the direct matching, has identity and 55 triangle holonomies. Its number of correlations is exactly

    (55 − S_ij (S²)_ij) / 2.

By d≡1 mod3 we get (S−I)²=0 mod3 (diagonal uses 56 mod3). In the normalized gauge where the first row is all +1 off-diagonal, this gives a binary edge constraint on the remaining 56 branch labels: encode the minus edges by their 0/1 matrix K, then each row has sum ≡1 modulo 3, and K²≡K−J−I modulo 3. Directly: for edge pair not involving the anchor, the number of relative correlations is either the number of other bits differing between the two rows, or its complement in 55; expanding row sums gives the stated equation. This is an exact numerical filter on an intrinsic binary relation on branches, not a second graph symmetry assumed for free.

In particular the correlation sign matrix cannot have a cyclic chart of length 57 in *any* gauge where the seams lie in this collineation/correlation model. If S is circulant, let f in F_3[x]/(x^57−1) represent its first row minus I. The relation says f²=0. Now x^57−1=(x^19−1)^3 and x^19−1 is square-free, so (x^19−1)² divides f. This square equals 1+x^19+x^38 in F_3. Thus the coefficient sequence of the degree<57 representative repeats every 19, and rows separated by 19 of S−I coincide. Such i,j have S_ij=−1 and S_ik=S_jk for all other k, so every one of their 55 relative triangle holonomies is a correlation. This contradicts d≤49. Circulant signs are a construction ansatz, not a necessary symmetry of the missing graph.

The noncyclic duality-mixed Gewirtz models and A56/S56 are not settled. The split and the congruences are intended for a simultaneous matching search rather than a claim that a single mixed sharp set has been found.
