# Fano cycle supports, integral circulation decorations, and rational SCC counting

Session 23b, 2026-10-21. Status: proved in the stated scopes, with exact finite controls. The new contribution is the oriented-integral lift and the correctly typed decoration/quotient formulas, not a claim that cycle spaces, Fano planes, or rational cone enumeration are new. No metagraph adjacency or Euler-tour enumeration is asserted.

## 1. The canonical tetrahedral Fano plane and the orientation information it loses

Let E be the six edges of K4, and let

    Z = {S subset E : every vertex has even S-degree},

with addition symmetric difference. This is the binary cycle space, of dimension 6-4+1=3. Its seven nonzero elements are precisely the four triangular faces and three Hamiltonian quadrilaterals. Projectivizing Z gives a Fano plane: its lines are the triples {a,b,a+b}. This construction needs no chosen basis. Identifying it with a particular labelled PG(2,2) does require one.

The tetrahedron supplies additional structure: a four-point triangle set and a three-point quadrilateral set. Its vertex permutations induce S4 of order 24, exactly the subgroup of GL(3,2) of order 168 preserving that 4+3 partition. To see exactness without counting all matrices, the four face points span Z and have a unique dependence, their sum zero; therefore every permutation extends uniquely to a linear map, and every map preserving the four points is determined by its permutation. All 24 permutations are vertex-induced.

For an orientation T of K4, retain only the edge supports of its directed simple cycles. The complete census is:

- 24 transitive tournaments: no directed-cycle support.
- 8 tournaments consisting of a source followed by a directed triangle, and 8 with a triangle followed by a sink: one triangle point.
- 24 strong tournaments: two triangle points and one quadrilateral point, and these three points are NOT collinear.

Here is an exact finite proof of completeness, independent of a search cutoff. The four representatives have degree multisets (0,1,2,3), (1,1,1,3), (0,2,2,2), and (1,1,2,2), so are nonisomorphic. Their automorphism group orders are 1,3,3,1. Their labelled orbits therefore have sizes 24,8,8,24, totaling all 2^6 orientations. The strong representative is treated explicitly below, including its three cycles and noncollinearity.

Take the strong tournament with ordered arc list

    e = (01,12,20,13,30,23),

where 01 denotes 0 -> 1. Its only directed simple cycles are

    A = 0 -> 1 -> 2 -> 0,
    B = 0 -> 1 -> 3 -> 0,
    Q = 0 -> 1 -> 2 -> 3 -> 0.

A XOR B is the undirected quadrilateral with edge set {12,20,13,30}, not Q. Its inherited arc directions are not a directed cycle. Thus A,B,Q are linearly independent in Z and form a Fano noncollinear triple. Mod-2 cancellation of a shared directed arc is not cancellation in the oriented nonnegative circulation cone.

There are exactly 12 such support triples in the tetrahedral model: choose the two faces in 6 ways, and choose either of the two quadrilaterals other than their XOR. Each triple is produced by exactly two tournaments, T and T^op. The S4 stabilizer of a strong support triple has order 2: identity and the unique anti-isomorphism of T. The actual tournament automorphism group is trivial. Consequently a support-only encoding cannot identify these two orientations or identify its stabilizer with Aut(T).

## 2. A free integral cone on the three directed cycles

A nonnegative integer circulation on a digraph D is an arc vector f in Z_{≥0}^E satisfying inflow=outflow at every vertex. Its weight is

    |f| = sum_e f_e.

This is total arc use with multiplicity, not the number of nonzero arcs, the number of cycles in a decomposition, or the number of Euler tours.

For the displayed T4, write a=f_20, b=f_13, c=f_23. The balance equations force

    f = (a+b+c, a+c, a, b, b+c, c).                 (1)

Conversely every a,b,c>=0 gives a circulation. The same equations hold over Z without positivity, and a,b,c are actual edge coordinates, so the three primitive cycle vectors are a lattice basis of ker(B) intersect Z^6. In particular the cone is simplicial and unimodular in its circulation lattice, and the integral semigroup is free. This is stronger than merely knowing that every integral circulation decomposes into cycles; in larger tournaments such decompositions are usually nonunique.

Since |f|=3a+3b+4c, and Aut(T4)=1, its ordinary flow-orbit generating function is

    f_T4(q) = 1 / ((1-q^3)^2 (1-q^4)).             (2)

The unique anti-isomorphism s=(01)(23) acts on arc vectors by (u,v) -> (s(v),s(u)). In the coordinates (1), it swaps a and b and fixes c. Therefore

    k_T4(q) = 1 / ((1-q^6)(1-q^4))                (3)

counts flow-decorated T4 types isomorphic to their converse. It is not the generating function for all flows over the underlying self-converse T4. For example A alone has a=1,b=c=0 and is not fixed.

The original tournament has H=5 Hamiltonian paths, h=H/|Aut|=5, and b=3 blue Hamilton-word orbits; these are underlying tournament quantities. No correspondence identifying a circulation with a Hamiltonian path, Hamilton word, or metagraph edge is used.

## 3. A general rational fixed-cone and anti-coset lemma

Let C be a finite tournament and G=Aut(C). Let K_C be the cone

    K_C = {x in R_{≥0}^{E(C)} : B_C x=0}.

For any permutation a of the arcs induced either by an automorphism or by an anti-isomorphism, put

    R_(C,a)(q) = sum_{x in K_C intersect Z^E, ax=x} q^{sum x_e}.

This is a rational power series. Here is a finite-cone proof sufficient for this setting. The constraints x>=0, B_C x=0, (a-I)x=0 are homogeneous rational linear constraints. The resulting cone is rational and pointed, being contained in the nonnegative orthant. Its span is rational; use the lattice obtained by intersecting that span with Z^E. A rational triangulation divides it into finitely many simplicial cones. Half-open facets, or inclusion-exclusion on the finite intersections, make the subdivision exact. For a simplicial cone generated by independent integral lattice vectors v_1,...,v_r, its lattice-point series is

    (sum_{p in Pi intersect lattice} q^{sum p_e}) /
       product_i (1-q^{sum (v_i)_e}),              (4)

where Pi is an appropriate half-open fundamental parallelepiped. It has finitely many lattice points. Every nonzero generator has positive total coordinate sum, so every denominator factor has positive exponent and specialization to this one-variable series is well-defined. The zero cone contributes 1. Thus a finite sum of such rational functions proves the claim, including lower-dimensional fixed cones. Unimodularity is NOT needed in this lemma.

This is the standard rational-cone argument, not a new polyhedral theorem. A primary exposition with the parallelepiped formula and the triangulation argument is Beck, Haase, and Sottile, *Theorems of Brion, Lawrence, and Varchenko on rational generating functions for cones*, Mathematics Intelligencer 31(1), 2009, 9-17, proof section, equations preceding and including (5):
https://arxiv.org/abs/math/0506466v4
https://arxiv.org/html/math/0506466v4

Because each fixed-weight slice is finite, ordinary Burnside gives the rational series

    f_C(q) = (1/|G|) sum_{g in G} R_(C,g)(q).       (5)

If C is self-converse, tag each ordinary map with a plus sign and each anti-isomorphism with a minus sign. The anti-tagged maps form a coset A of the ordinary-tagged group G, of size |G|. The tags keep these copies disjoint even for K1 (and the empty object), where their underlying vertex permutations coincide. Converse induces an involution on the G-orbits of flows, even when a chosen representative anti-isomorphism is not itself an involution. The series counting converse-fixed G-orbits is

    k_C(q) = (1/|G|) sum_{a in A} R_(C,a)(q).       (6)

To prove (6), apply Burnside to the disjoint tagged extended group G union A and compare with the two-element action on the G-orbit set. An orbit fixed by converse contributes one to the right side; a nonfixed pair contributes zero. This avoids the incorrect assumption that a single chosen anti-isomorphism fixes every self-converse decorated orbit. If C is not self-converse, put k_C=0.

These series have nonnegative integral coefficients despite the displayed averages. Also f_C(0)=1 for every C, and k_C(0)=1 precisely for self-converse C. Setting q=0 forgets the decorations by retaining only the zero flow. Setting q=1 generally diverges and is NOT the forgetful specialization.

## 4. Why circulation decorations factor over ordered SCCs

Every nonnegative circulation has zero weight on every arc between different SCCs. For a tournament, list its SCCs in their unique forward order. Sum the balance equations over a proper initial union of SCCs. No arc enters this union, so the sum of the nonnegative weights leaving it is zero; each such weight vanishes. Apply this to every SCC cut.

Consequently the flow semigroup is a Cartesian product over the SCCs. Every tournament isomorphism preserves the unique SCC order, so its automorphism group is the direct product of the SCC automorphism groups, without a permutation factor for repeated isomorphic SCCs. Flow-orbits therefore also factor componentwise. Converse reverses the order and sends each decorated atom to its opposite decorated atom.

For the alphabet {K1,C3,T4}, the atomic series are

    z,   z^3/(1-q^3),   z^4/((1-q^3)^2(1-q^4)).

The C3 automorphism group rotates all three equal edge weights and therefore acts trivially on its flow parameter. The ordinary decorated-word series is exactly

    W(z,q) = 1 /
       (1-z-z^3/(1-q^3)-z^4/((1-q^3)^2(1-q^4))).   (7)

A converse-fixed decorated word is a free first half, an optional fixed center, and the forced opposite second half. Hence

    J(z,q) = W(z^2,q^2) *
       (1+z+z^3/(1-q^3)+z^4/((1-q^6)(1-q^4))).     (8)

Both vertex count and total edge use double in a mirrored pair. The ordinary converse-merged series is (W+J)/2. Equations (7)-(8) enumerate actual decorated isomorphism types, not labelled edge assignments or sets of cycles.

## 5. Extension at fixed UNDERLYING h

Inherit the audited finite-strong-alphabet and center laws from

- `fixed_degree_metagraph_rationality_20261021/session22-fixed-degree/PROOF.md`;
- `effective_black_defect_order_bound_20261020/PROOF.txt`;
- `metagraph_classification_20261021a/README.md`;
- `THM-2454-pure-blue-alphabet-and-center-law.md` and its current supersession.

Here h(C)=H(C)/|Aut(C)| refers to the underlying undecorated tournament; h multiplies over ordered SCCs. The only h=1 strong atoms are K1,C3. For every fixed odd j there are finitely many strong atoms with h=j. The established effective order bound is sufficient; no new classification is assumed.

Set

    A(z,q)=1/(1-z-z^3/(1-q^3)),
    P_j(z,q)=sum_{strong C: h(C)=j} z^|C| f_C(q).

Each P_j is a finite sum of rational series. For j odd, let F_j count all ordinary flow-decorated tournament types with underlying h=j, including empty at j=1. Unique ordered factorization gives

    F_j = sum_{j_1...j_r=j, j_i>1}
                 A^(r+1) product_i P_(j_i).        (9)

There are finitely many tuples, r<=Omega(j). Thus F_j is rational. This is an effective rationality theorem using finite alphabets, not an assertion that those alphabets have all been computed.

Let C_a^fix(z,q) be the sum of z^|C| k_C(q) over strong self-converse underlying C with h=a, together with the empty center of weight 1 when a=1. Then the actually self-converse decorated series at underlying h=j is

    J_j = sum_{x^2|j} F_x(z^2,q^2) C_(j/x^2)^fix.  (10)

Therefore the converse-merged flow-decorated count at fixed underlying h=j is rational, with series (F_j+J_j)/2.

It would be wrong to replace h(C) with H(C)/|Aut(C,f)| without changing the entire stratification. Decorations can shrink stabilizers. A directed triangle flow in R5 has stabilizer 1 rather than 5, so H(R5)/|Aut(R5,f)|=15 although its underlying h is 3. The present theorem does not count by such a decorated h or assign a new metagraph degree to decorations.

## 6. Fixed UNDERLYING fundamental black degree: the third carrier

This distinction is essential: a decoration over a self-converse underlying tournament need not be self-converse. To count by the degree of the UNDERLYING converse-merged node, one must count independently decorated opposite halves as well as mirrored halves.

For a half-atom C, independent decoration of C and C^op has weight

    z^(2|C|) f_C(q)^2,

not z^(2|C|) f_C(q^2). Define

    A_ind = 1/(1-z^2-z^6/(1-q^3)^2),
    P_j^ind = sum_{strong C: h(C)=j} z^(2|C|) f_C(q)^2,
    G_j = sum_{j_1...j_r=j, j_i>1}
                  A_ind^(r+1) product_i P_(j_i)^ind.          (11)

The finite tuple rule is the same as (9). Equivalently G_j sums z^(2|U|) f_U(q)^2 over underlying first-half types U with h(U)=j. It is not the square of F_j, which would independently choose two possibly different underlying halves.

Let C_a^all sum z^|C| f_C(q) over the same underlying self-converse centers used for C_a^fix, with the same optional empty center. All ordinary decorations over underlying self-converse h=j have series

    U_j = sum_{x^2|j} G_x C_(j/x^2)^all.             (12)

Next recall the underlying center law: if T=U -> C -> U^op, with x=h(U), a=h(C), beta=b(C), then its underlying deficit is

    d=h(T)-b(T)=x(xa-beta).

For d>0, define E_(x,d)^all and E_(x,d)^fix by summing z^|C| f_C(q) or z^|C| k_C(q), respectively, over actual empty/strong self-converse centers satisfying this equation. The center alphabets are finite: for x=1 use the inherited positive-deficit strong-order bound; for x>1 use

    x divides d,   x(x-1)<=d,   a<=d/[x(x-1)].

Define

    L_d = sum_x G_x E_(x,d)^all,
    M_d = sum_x F_x(z^2,q^2) E_(x,d)^fix.            (13)

Thus L_d counts all ordinary decorations over underlying self-converse deficit d, whereas M_d counts the actually self-converse decorated ones. Both are rational.

For k>=1 let B_k^flow count decorated converse-orbits whose UNDERLYING merged tournament node has fundamental paired-black degree k. Then

    B_k^flow = ( L_(2k) + M_(2k)
                       + 1_(k odd)(F_k-U_k) )/2.  (14)

Proof: over an underlying self-converse type of deficit 2k, quotienting decoration-orbits by converse gives (all+fixed)/2. Over underlying non-self-converse types with h=k, converse has no fixed decorated types and pairs the two underlying types, giving (F_k-U_k)/2. The strata are disjoint. A factor 1/2 on L alone or the substitution M for L would miss valid decorations.

The inherited degree-zero classification has rigid halves and empty/K1/C3/R5 centers. Hence, writing f_R5 and k_R5 for the rational series (5)-(6),

    B_0^flow = 1/2 * [
       A_ind*(1+z+z^3/(1-q^3)+z^5 f_R5(q))
       + A(z^2,q^2)*(1+z+z^3/(1-q^3)+z^5 k_R5(q)) ].          (15)

This proves rationality at every fixed nonnegative underlying degree. It does not assert an unchanged denominator from the undecorated theory: atom-cone denominators appear, and independent pairing introduces the genuinely different carrier A_ind. All formulas recover the inherited undecorated formulas at q=0.

Minimal hostile: the underlying tournament C3 -> C3 is self-converse with degree zero. Give one cycle weight 0 and the other weight 1. The decoration is not self-converse; its two ordered assignments form one converse pair. At z^6 q^3 the neutral-alphabet all-over-underlying-SC series has coefficient 2, while its actually self-converse series has coefficient 0. Thus the two denominator carriers cannot be interchanged. In the larger K1,C3,T4 test alphabet the same coefficient is 4 versus 0, because K1 -> T4 -> K1 adds two further nonfixed decorations.

## 7. R5 and the fractional stable-set cone of the pentagon

For cyclic R5 on Z/5 with arcs i->i+1 and i->i+2, let x_i,y_i be the corresponding flow values. The five balance equations are equivalent to

    x_i = t-y_i-y_(i-1),
    y_i>=0,   t>=y_i+y_(i-1),                       (16)

with one common integer t. The equivalence follows by subtracting successive expressions x_i+y_i+y_(i-1), and the total arc use is 5t-sum y_i. Thus its circulation cone is integrally isomorphic to the cone over the fractional stable-set polytope of the pentagon, not the ordinary integral stable-set polytope.

The twelve extreme cycle rays comprise five triangles, five quadrilaterals, and two Hamiltonian cycles. In (t,y) coordinates they are

    (1,0), (1,e_i), (1,e_i+e_(i+2)), (2,1,1,1,1,1),

where the latter notation means t=2 and all five y coordinates are 1. The first ray has weight 5, the singles weight 4, the independent pairs weight 3, and the last ray weight 5. The half-integral pentagon vertex matters: dropping it removes one Hamiltonian cycle ray.

A nontrivial rotation fixes exactly the flows with x_i=x and y_i=y, giving weight series 1/(1-q^5)^2. If R(q) denotes the raw integral-flow series, then

    f_R5(q) = (R(q)+4/(1-q^5)^2)/5.                (17)

All five anti-isomorphisms are rotation-conjugate. Under i->-i, the y coordinates have form (a,b,b,a,c). Therefore

    k_R5(q) = sum_{a,b,c>=0, t>=max(a+b,2b,a+c)}
                      q^(5t-2a-2b-c).             (18)

This is an exact rational-cone expression, not yet a reduced rational numerator formula. The fixed-cone lemma proves it rational. The accompanying test verifies the coefficients independently through total use 14. This connection gives an explicit polyhedral target for further formula work, without pretending that independent primitive-cycle coordinates exist for R5 as they do for T4.

## 8. Scope and provenance guardrails

The local prior search found inherited binary-cycle-space and basis-choice results in THM-4069 and `tournament_even_transfer_20261011b.md`; the K4 section XOR code in the results index already distinguishes an orientation section from the actual cycle space. General integral-circulation and decomposition ideas occur in `fivefold_boundary_20261018g.md`, `finite_extraction_and_inverse_clocks_20261017i.md`, and THM-3044. These are acknowledged background, not rediscovered claims. The targeted search did not find this exact oriented T4 semigroup plus correctly stratified SCC flow-orbit extension. This limited local search is not a literature-priority claim.

A circulation is not an Euler trail. It records multiplicities, forgets traversal order, and may have disconnected support. On C3 -> C3 with positive flow in both triangles, the support is disconnected because every inter-SCC edge has zero circulation, so there is no single Euler tour of its nonzero support. Even when support is connected and Eulerian, different tours can yield the same circulation. Neither cycle-cone rationality nor Fano incidence creates metagraph ports, adjacency, bridges, connectivity, Collatz trajectories, or a solution to another open conjecture.

## 9. Reproducible controls

`check.py` uses explicit failure checks that survive Python -O. Normal and optimized JSON outputs agree byte-for-byte. It checks all 64 K4 orientations, all Fano incidences and 168 projective automorphisms, the actual S4 subgroup, all 134,596 nonnegative T4 edge vectors of total at most 18, 5,054 decorated SCC words at (vertices<=10,total use<=12), their literal converse orbits, both independent and mirrored half formulas, and the corrected fixed-underlying-degree formula. It additionally checks 611 R5 raw flows through total use 14, their 127 automorphism orbits, anti-coset Burnside, and a separate direct edge-vector enumeration through total use 8. These are controls of the proofs, not extrapolated proofs or a census outside the stated bounds.
