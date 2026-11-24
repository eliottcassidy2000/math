# The second neighborhood forces one of the two proposed Moore(57,2) critical groups

Conditional throughout on the existence of a Moore graph of valency 57
and diameter 2. No existence or nonexistence conclusion follows.

Let L be its Laplacian and e0 its rank modulo 5. Then the 5-primary
critical group is necessarily

    (Z/5)^(1521−e0) ⊕ (Z/25)^(1730−e0) ⊕ (Z/125)^(e0−2).

This selects the second of the alternatives in [Ducey, Theorem 3.2](https://arxiv.org/abs/1509.00327) and [Ducey et al., Example 3.5](https://arxiv.org/abs/1910.07686). No independent literature-priority
claim is made. The argument also gives 196≤e0≤1521, an explicit smaller
matrix for e0, and constraints on the local eigenlattices.

## The local split

Fix a vertex o. Its 57 neighbors are independent. Its 3192 distance-two
vertices split into 57 fibers F_i, 56 vertices each, according to their
unique neighbor in the first neighborhood. There are no edges within a
fiber, and two different fibers are joined by a perfect matching.
Indeed two neighbors in the other fiber would make a 4-cycle, and every
second-layer vertex has 56 second-layer neighbors to place across 56
other fibers. Let B be the second-layer adjacency. If P_ij is the
matching block with rows in F_i and columns in F_j, then

    P_ij + Σ_(k≠i,j) P_ik P_kj = J_56      (i≠j).

The 56 summands partition all positions, by the unique common-neighbor
condition on nonadjacent vertices and triangle-freeness on adjacent
vertices. Conversely these identities and P_ji=P_ij^t are sufficient:
the diagonal blocks of B² are 56I; adding the prescribed root and
first layer gives a regular graph with A²+A=56I+J, hence the stated
Moore parameters. One fiber's incident matchings may all be set to I
by independent relabeling of the other fibers.

Over Z_5, the full vertex lattice is a direct sum of a local fluctuation
part W and its complement Q. The first is zero on the root and the
neighbor layer and has sum zero in each fiber. The second is arbitrary
on the root and neighbor layer and constant in each fiber. This is an
integral local splitting because 56 is a unit. Both parts are preserved
by adjacency. Over Q the adjacency on W is B|W and satisfies
(B−7)(B+8)=0, since the global J acts trivially. Its trace is zero:
B has zero diagonal, and its action on the complementary fiber constants
is J_57−I, also trace zero. Since dim W=3135, the multiplicities
are 1672 for 7 and 1463 for −8.

## 5-primary Smith calculation

On W put D=L|W. Its eigenvalues are 50 and 65, so

    D²−115D+3250I=0,
    125 D⁻¹=(115I−D)/26.

The right side is integral over Z_5. Hence every elementary-divisor
valuation of D is at most three. The valuations of 125D⁻¹ are those
of D reflected around 3/2 by taking its Smith form. But 125D⁻¹ and D
differ modulo 5 by multiplication by a unit. Therefore if
r=rank_F5(D), both valuation zero and valuation three occur exactly
r times. Dimension and determinant valuation 2·1672+1463 give
1463−r factors of valuation one and 1672−r of valuation two.

On Q, with variables (root t, first-neighbor coordinates x, and
fiber constants y), adjacency acts by

    (t,x,y) ↦ (sum x, t·1+56y, x+(J_57−I)y).

Since 57 is a unit at 5, split x,y into their mean and zero-sum
parts. On the 56 copies of the zero-mean pair, L has matrix

    [[57,−56], [−1,58]].

Its determinant is 3250=26·125 and it has a unit entry, so it contributes
one unit and one valuation-three factor per copy. On the three
constant variables, L has matrix

    [[57,−57,0], [−1,57,−56], [0,−1,1]].

It has rank two and the minor in its last two rows and first two
columns is 1. Thus it contributes two units and one zero. Altogether,
e0=58+r and the nonunit finite factors are 1463−r at 5, 1672−r at 25,
56+r at 125. Substitution proves the result and e0≤1521.

## Three-fiber rank witnesses at 5 and 3

On three chosen fibers choose transport labels so both matchings from
the first are identity. Then the last matching P is a derangement
(otherwise a triangle occurs). On the direct sum of the three
55-dimensional augmentation spaces consider the principal compression
of N=B−7. Eliminate the first diagonal block aI, a=−7. The remaining
diagonals are (a−a⁻¹)I and the off-diagonals P−a⁻¹I, P⁻¹−a⁻¹I.

At 5 this gives diagonal I and an additional Schur complement
I−(P⁻¹−2I)(P−2I), equivalent after multiplication by P to
2(P−I)^2. On the full permutation space each cycle of length l
contributes one kernel dimension if 5∤l and two otherwise (the
cyclic module is F_5[X]/(X^l−1)). This is at most floor(l/2) since
l≥2, and the augmentation space removes one constant dimension
because 56 is a unit. Hence the last Schur kernel has dimension
at most 27. The three-fiber block has rank at least 165−27=138,
so r≥138 and e0≥196.

At 3 the remaining diagonals vanish, and the two off-diagonals are
P+I, P⁻¹+I. On the full permutation space their nullity is the number
of even cycles (at most 28), and the constant dimension is not in the
kernel. Thus the analogous rank on W at 3 is at least
55+2(55−28)=109. These are necessary filters, not claimed sharp.

## Discriminant groups of the fluctuation eigenlattices

Give W its inherited integral Euclidean form, so W is the orthogonal
sum of 57 copies of the augmentation root lattice A_55 and its
discriminant group is (Z/56)^57. Every permutation of one fiber acts
trivially on its cyclic discriminant; hence the action of B on the
whole discriminant is J_57−I. Let E_+, E_- be the integral eigenlattices
in W for 7 and −8 respectively, and let r_p=rank_Fp(B−7|W) at p=3,5.

At p∤15 the eigenprojectors split W over Z_p, orthogonally; its
discriminant decomposes accordingly. At 2 the 56-dimensional
augmentation of the fiber-discriminant group has B eigenvalue −1,
which agrees with 7 mod 8, and its constant line has eigenvalue
56, which agrees with −8 mod 8. At 7 these roles reverse. Thus E_+
has 56 copies of Z/8 and one Z/7 in its discriminant; E_- has one
Z/8 and 56 copies of Z/7.

For p=3,5, W is unimodular over Z_p. The lattice quotient
W/(E_+⊕E_-) is then naturally isomorphic to the discriminant
group of either eigenlattice: orthogonal projection into its dual
is onto, since its primitive inclusion allows all integral
linear forms to extend to W, which is itself unimodular.
The same lattice quotient identifies with
(B−7)W/15E_- by the kernel and its action on E_-. It is killed
by p, and its dimension is r_p because E_- is primitive in W and
15 is a unit times p. In particular

    detGram E_+ = 2^168 · 7    · 3^r3 · 5^r5,
    detGram E_- = 2^3   · 7^56 · 3^r3 · 5^r5.

These conditional restrictions and the matching identities reduce
any finite-field certificate to a structured second-neighborhood
calculation. They do not determine r_3 or r_5 or solve existence.
