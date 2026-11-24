# Which tournament representations can preserve the native symmetries?

A one-tournament representation on the original vertex set is much stronger than a marked observer tournament, a family of fixed-spine tilings, or a representation with a choice of order. These notions have different obstructions.

## A two-value Seidel Gram is impossible at degree 3, 7 or 57

**PROVED conditional for 57.** Let A be the adjacency of a diameter-two Moore graph of degree d, and let n=d²+1. Let S be a tournament's real skew Seidel matrix on those very same n labels, i.e. S_ii=0 and S_ij=−S_ji∈{−1,+1}. Suppose the full row Gram distinguishes adjacency by two values, equivalently

    SS^T = aI+bA+cJ,      b≠0.

Then d cannot be 3, 7 or 57.

For these odd degrees n≡2 mod 4. Let p_i be the product of the n−1 signs in row i. For i≠j, the product of the n−2 terms S_ik S_jk, k outside i,j, is −p_i p_j. Therefore the off-diagonal Gram entry is 2 modulo 4 when p_i=p_j and 0 otherwise. In particular both off-diagonal values c and b+c are even. If b is 2 modulo 4, this formula would make A either union of two cliques or complete bipartite across the parity split. Neither is a Moore graph of degree at least three and girth five. Hence b≡0 modulo 4. The off-diagonal residue must then be 2 (it cannot be 0, since at least two of the n products have the same sign), so c≡2 modulo 4.

On the other hand S is invertible for even n: modulo 2 it is J−I, invertible in even dimension. Each eigenvalue of SS^T has even multiplicity, since a nonsingular real skew matrix is an orthogonal direct sum of two-dimensional rotation blocks. The Moore eigenvalues and multiplicities are

| d | n | spectrum |
|---:|---:|---|
| 3 | 10 | 3 once, 1 five times, −2 four times |
| 7 | 50 | 7 once, 2 twenty-eight times, −3 twenty-one times |
| 57 | 3250 | 57 once, 7 seventeen hundred twenty-nine times, −8 fifteen hundred twenty times |

The table follows from A²+A=(d−1)I+J and trace A=0. As b≠0 the two nonprincipal Gram eigenvalues are different. The principal Gram value must coincide with the nonprincipal value of odd multiplicity, forcing b(d−τ)+cn=0 for τ=1,−3,7 respectively. Thus c=−b/5,−b/5,−b/65 respectively. These odd denominators force c≡0 modulo 4, a contradiction.

The pentagon genuinely passes this stronger representation: orient vertex i to i+1 and i+2 modulo 5. Its Gram offdiagonal is 1 for cyclic neighbors and −3 for the other pairs. The same cyclic-half orientation on an odd cycle of order n yields Gram value n−4h at cyclic distance h≥1, so the largest off-diagonal value recovers the cycle even at larger diameter. For h between 1 and (n−1)/2 the two rows have 2h−1 disagreements away from the two zero sites (count the two interval boundaries of the cyclic sign word), proving the formula. Complete graphs are the trivial Moore diameter-one case. The focus here is the nontrivial diameter-two list.

The [observer construction](moore-tournament-representations.md) evades this obstruction: an off-diagonal block of a larger tournament is allowed to be symmetric in its paired labels. Its rectangular Gram preserves both geodesic counts exactly. That is not a claim about the full skew Gram.

## Faithful symmetry has an exact criterion

A noncomplete diameter-two Moore graph of degree d≥2 has an ordinary uncolored tournament encoding with its *entire* automorphism group realized faithfully, and no extra automorphisms, iff its automorphism group has odd order. In any tournament, an involution exchanging a vertex pair would reverse their arc; hence the tournament group is odd by Cauchy, proving necessity.

For sufficiency choose a tournament T on the graph vertices whose arcs are invariant under Aut(G): every unordered pair orbit can be oriented consistently, as an odd-order group has no pair stabilizer exchanging endpoints. Take copies L,R of T. Orient L_u→R_v exactly if u=v or u∼v. Add two markers a→b with L→{a,b}→R.

The two markers are the unique twins (vertices with identical arc directions to every other vertex). Distinct vertices in one layer have different closed-neighborhood incidence to the other layer; vertices in opposite layers are distinguished by the markers. A marker and a layer vertex are distinguished by a nonneighbor in the opposite layer. Hence the markers and both layers are intrinsically recovered, with the markers individually fixed.

In the bipartite support of cross arcs L→R, every diagonal support entry lies in d rectangles of four support entries and every other support entry lies in exactly one. Indeed any two row neighborhoods intersect in their two labels if adjacent, and in their unique common neighbor otherwise. Rectangles therefore recover the pairing L_v↔R_v, and then every graph edge is recovered. All graph automorphisms preserve both copies of T and extend uniquely; all tournament automorphisms preserve the recovered graph and act faithfully on its labels.

For the known pentagon, Petersen and Hoffman–Singleton graphs there are involutions, so one cannot require this symmetry of ANY ordinary encoding. Nevertheless all admit the rigid observer constructions. For degree 57, the symmetry-preserving construction applies conditional on existence together with [Ishida's no-involutions result](https://arxiv.org/abs/2606.29183); no such premise is needed for the marked observer formula. Allowing anti-automorphisms changes the native action, as detailed in the [observer note](moore-tournament-representations.md).

## A directed-triangle version keeps the actual middle label

For verification, take three distinguished copies L,M,R of the graph labels and any internal tournament orientations. Orient L_u→M_w if u=w or u∼w, orient M_w→R_v if w∼v, and orient every R_v→L_u. Reverse unspecified arcs between their respective layers. This is an ordinary tournament with three marked roles. The directed triangles of the form L_u→M_w→R_v→L_u record exactly the candidate middle labels w. With w=u allowed, they include the direct route. Their count is ((I+A)A)_(u,v), i.e., d when u=v and one when u≠v in a Moore graph. Thus the same family of directed-triangle tests certifies every endpoint pair jointly, without counting all cycles of the tournament or claiming a Hamiltonian-height equivalence. In a partial graph the observer set of each pair says precisely which intermediate routes exist and which are missing.
