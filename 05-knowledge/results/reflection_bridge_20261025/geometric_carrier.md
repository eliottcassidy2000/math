# A projective carrier for the divisor action, and where its seven lives

2026-10-25. **PROVED** by the explicit construction below. The pointed-monoid projective model and its point conventions are **CITED**, not new. Finite coordinate controls are **FINITE-EXACT**. This is not a divisor-to-tournament equivalence.

## 1. An actual geometric action

Fix distinct primes p,q,r and k>=2. Let L_k be all divisors p^a q^b r^c of p^kqr, and D_k=L_k minus its unit and maximum. Set

    X_k = P^k x P^1 x P^1,
    i(a,b,c) = ([e_a],[e_b],[e_c]),

where e_a is the a-th coordinate vector, indexing from zero. Work either over Z or in the standard pointed-monoid projective model over A_1={0,1}. The coordinate points are the rank-zero points of the monoid model and the fixed points of the standard split torus after base extension. Thus i identifies L_k with precisely these 4(k+1) points.

The following coordinate permutations are defined over the pointed monoid, hence also over Z:

    c([x_0:...:x_k],[y_0:y_1],[z_0:z_1])
      = ([x_k:...:x_0],[y_1:y_0],[z_1:z_0]),
    s([x],[y],[z]) = ([x],[z],[y]).

They commute, square to the identity, and on i(L_k) act exactly as

    C(a,b,c)=(k-a,1-b,1-c),  Q(a,b,c)=(a,c,b).

Delete the two coordinate points i(0,0,0) and i(k,1,1) in the monoid model. Over Z, delete their two coordinate sections, each isomorphic to Spec Z; over every field fibre these are two rational points. They form a closed, c- and s-invariant union. Therefore the open model

    U_k = X_k \ {i(0,0,0), i(k,1,1)}

has honest commuting geometric automorphisms c,s, and its rank-zero locus is exactly D_k. There is no split fibre, no identification of distinct divisors, and no ambiguity about which complement is used. The two axis reflections are s and cs; their product c is divisor complement. This is the positive geometric carrier, as distinct from the obstructed graph carrier.

For k=2, C,Q,CQ have respectively 0,4,2 fixed points **on this ten-point rank-zero locus**. The entire geometric fixed loci can be larger. Even the literal pointed-monoid Hom set is larger, as section 3 shows.

## 2. Recovering 3-4-3 and 4-2-4 with a flag sidecar

Fix the coordinate order on each projective factor. The standard cell with last nonzero coordinate a is A^a and contains exactly the coordinate point e_a. The product cells consequently have dimensions a+b+c. Their enumerator is

    R_k(t)=(1+t+...+t^k)(1+t)^2.

At k=2 this is 1+3t+4t^2+3t^3+t^4. Removing the two **coordinate-point labels** at the extreme dimensions gives the rank-level counts 3,4,3 on D_2. These are the same divisibility ranks as in the divisor lattice. They are also its rank enumerator, without any geometry assumption.

Separately, the faithful plane coordinates

    f_k=2a-k(b+c),  g=b-c

are weights read from the marked coordinate labels. At k=2 the signs of g have sizes 3,4,3, and the signs of f_2 have sizes 4,2,4. Rank and sign(g) are two different 3-4-3 decompositions; neither should be substituted for the other. On these coordinates s sends (f,g) to (f,-g), and cs sends it to (-f,g).

What is geometric and what requires the sidecar? Coordinate permutations are actual automorphisms. Divisibility and the rank on the closed coordinate points require the chosen coordinate flags: all these points are closed, so their divisibility order is **not** their specialization order. The reversal c takes the standard flags to the opposite flags; it does not take each whole standard cell to the standard cell with complementary index. Already [x_0:x_1] with x_1 nonzero can have x_0 nonzero or zero, so its coordinate reversal has no fixed last-nonzero index.

In particular U_k is obtained by deleting two points, not by deleting the two extreme cells. Its finite-field count is

    #U_k(F_Q)=(1+Q+...+Q^k)(1+Q)^2-2.

For k=2 it is -1+3Q+4Q^2+3Q^3+Q^4, **not** 3Q+4Q^2+3Q^3. These happen to have the same formal value 10 at Q=1, but differ by Q^4-1. This hostile control prevents a rank enumerator from silently becoming the wrong geometric count. Subtracting two rational points commutes with every finite-field extension, so this is a uniform all-prime counting polynomial.

## 3. The same P^2 factor contains the honest seven-support object

In the pointed-monoid convention, a literal A_1-valued point of P^j is a nonzero {0,1}-vector, with trivial scalar quotient. It is therefore a nonempty subset of its j+1 coordinate labels. The seven such points of P^2 are the supports

    {p},{q},{r};  {p,q},{p,r},{q,r};  {p,q,r}.

This is also the underlying seven-point monoid spectrum. Its specialization direction is: if S is contained in T, the point x_S lies in the closure of {x_T}; larger supports are more generic. The closed/rank-zero points are only its three singleton supports. Its usual finite-field counting polynomial is 1+Q+Q^2, whose formal value at Q=1 is 3. The literal pointed Hom count is 7. These are different, explicitly specified constructions.

Consequently X_k has (2^(k+1)-1)*3*3 literal A_1-valued points, while its rank-zero count is 4(k+1). U_k has 9(2^(k+1)-1)-2 literal points and 4k+2 rank-zero points: at k=2 these are **61 and 10**, respectively. Only the maps whose images are one of the deleted closed points are removed; every source Spec A_1 has one point, so every remaining map factors through the open complement.

The seven-support object is thus present inside the *same projective factor* that supplies three coordinate labels to the ten-divisor product. This is a typed connection, not the equation 7=3. The extra structure on the seven supports is inclusion/specialization; the extra structure on the ten rank-zero points is the chosen coordinate-flag indexing.

The involution on seven obtained by collapsing the boundary pair of a Boolean three-cube (see [seven.md](seven.md)) is **not** the restriction of c or s. In this P^2 factor, c exchanges coordinate labels 0 and 2. On nonempty supports it fixes three supports: {1}, {0,2}, and {0,1,2}. The collapsed-complement involution fixes just the full support. Different fixed counts rule out a conjugacy of these two seven-point actions. Further, a coordinate permutation preserves support size, whereas collapsed complement exchanges sizes one and two.

## 4. Connection contract

| Source and target | Actual map | Preserved | Lost or separately retained |
|---|---|---|---|
| L_k to rank-zero X_k | exponent tuple to coordinate-point tuple | faithful set, C,Q,CQ and all their fixed points | divisibility and rank need coordinate flags |
| D_k to rank-zero U_k | restriction of the same map | all ten original states when k=2, both reflections | no tournament arcs, Hamiltonian paths, or graph complement |
| Three-coordinate supports to P^2(A_1) | indicator vector modulo its trivial scalar group | supports, coordinate permutations, specialization | field addition and Fano lines absent |
| Boolean boundary quotient to seven supports | collapsed boundary represented by full support | the specified one-fixed-point involution | inclusion and four of seven Fano lines fail |
| Rank enumerator to geometric point count | no equality after deleting only coordinate points | both give ten at formal Q=1 for k=2 | their full polynomials differ by Q^4-1 |

Cheapest decisive tests: the ten fixed counts (0,4,2); the seven fixed-count mismatch (3 versus 1); and #U_2(F_2)=61 rather than 46, while the rank-zero count remains ten. The proof of the scheme automorphisms is the coordinate formula and its chart compatibility, not these finite tests.

## 5. Sources and scope

- Inherited arithmetic: [THM-4733, divisor plane versus marked-path tiling cube](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md), sections 0 and 3, and [reflection resolution](../divisor_reflection_resolution_20261021/REFLECTION_RESOLUTION.md). Their graph obstruction is retained.
- Inherited model conventions and chart proof: [Session 59 cyclotomic support, sections 1–2](../session59_f1_geometry_and_information_20261025/cyclotomic-support/CYCLOTOMIC_SUPPORT_DYNAMICS.md).
- **CITED**, primary text checked 2026-10-25: Koen Thas, [Projective Spaces over F1^ell](https://backoffice.biblio.ugent.be/download/8589084/8603796), sections 1–4: distinguish closed points from morphism points and use pointed projective frames. Anton Deitmar, [Schemes over F1](https://arxiv.org/html/math/0404185v7), sections 1–3: monoid schemes and base extension; his original unpointed convention is not silently substituted for ours.
- **CITED** only for the counting-polynomial convention: Deitmar–Koyama–Kurokawa, [Counting and zeta functions over F1](https://arxiv.org/html/1306.3681v2), section 1. The product, removal, and action statements here are proved directly.

No historical-priority claim. This explicit monoid model is not a claim that an ordinary one-element field exists, that every tournament becomes transitive over F1, or that cardinality alone recovers graph/tournament geometry.
