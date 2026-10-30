# Session110 checkpoint 2: a closed transitive-action carrier and exact symmetry boundaries

2026-10-30. **Verified extension; the overnight session remains active.**
The portable JavaScript checks below were replayed directly after the shell service became unavailable. They use exact integers, explicit failures, and BigInt where necessary. They do not depend on pending filesystem writes.

## 1. Prime counts on transitive odd-group actions

Let G be a finite group of odd order and K a subgroup. The vertices are the **named left cosets** G/K. Let A(G,K) count G-invariant tournaments on that set, and P(G,K) count those with no nontrivial proper module. Set A(G,G)=1 and P(G,G)=0.

If d(G,K) is the number of double cosets K\\G/K, then

    A(G,K)=2^((d(G,K)-1)/2).

Ordered-pair orbits correspond to double cosets; reversing an ordered pair corresponds to inversion. A non-diagonal orbit cannot equal its reverse: that would produce an element interchanging two vertices, impossible for an odd-order permutation. Thus each paired orbit supplies one orientation bit.

The exact recursion is

    A(G,K)=P(G,K)+sum_(K<H<G) P(G,H) A(H,K).

This is incidence-algebra inversion on the **embedded subgroup interval**, not a recurrence on abstract subgroup orders alone. With incidence convolution,

    A=I+P*A,   P=I-A^(-1).

### Why the gluing is exact

Every invariant tournament is vertex-transitive and therefore strongly connected: its strong components would otherwise form a rigid linear order. In a strongly connected tournament, distinct maximal proper modules are disjoint. Overlapping modules have modular union; if two maximal ones covered all vertices, their intersection and differences would form a linearly ordered three-block decomposition, contradicting strong connectivity. Hence the maximal modules give a canonical partition with prime quotient.

Let B be its block containing K, and H its setwise stabilizer in G. Then K<=H and B=H/K. Indeed, if xK belongs to B, the blocks xB and B intersect, so coincide. The outer tournament is G-invariant on G/H; the inner tournament is H-invariant on H/K.

Conversely transport the inner tournament to every block gH/K using g. Changing g to gh, h in H, leaves the transported arcs unchanged by H-invariance. Combining these with a prime outer tournament gives the unique inflation with that canonical first block system. There is no section multiplicity and no assumption that H or K is normal.

G/H and H/K denote coset spaces, not quotient groups unless normal. The closure operation is modular inflation, not arbitrary attachment. Every vertex-transitive tournament fits some odd group action by taking its full automorphism group, but presentations overlap. This is not an unrestricted isomorphism census or an A000568 formula.

### Exact order21 controls

The independent [direct graph checker](action_carrier/parent_graph_certificate.js) enumerates all1024 inverse-paired connection sets and tests minimal module closures, without using subgroup counts.

| Regular group | Prime labelled banks | Composite labelled banks |
|---|---:|---:|
| C21 |992|32|
| C7 semidirect C3 |980|44|

For the nonabelian group, its seven subgroups of order3 each contribute2 prime seven-coset orientations times2 inner orientations, giving28. Its single subgroup of order7 contributes2 times8, giving16. These seven-coset orientations are the two Paley7 orientations. The cyclic group instead has only one subgroup of each order and contributes16+16.

This is an exact role for7 and21 through group actions. It does not identify these numbers with the separate Hamilton-path exclusions.

Classical modular decomposition, double-coset orbit counting, and incidence inversion are credited. The nearby paper [Chudnovsky et al., On prime Cayley graphs](https://arxiv.org/html/2401.06062v1) is background; its stated directed theorem excludes tournament connection sets, so it is not used as a direct theorem covering this argument.

## 2. Full graph classification at125 vertices inside the ramified bank

The source-labelled golden bank at uniformizer depth3 has16384 states:
16320 prime and64 composite. Its clock normalizer supplies408 prime orbits.

The [portable full graph verifier](portable_checks/graphs125.js) does not assume that normalizer or affine equivalence equals graph isomorphism. For each representative it builds the actual125-vertex adjacency matrix and performs graph-intrinsic ordered refinement with vertex0 marked. The resulting partition has five singletons and24 cells of size5. A canonically selected nontrivial cell is verified to be an orbit of the known order5 clock automorphism. Individualizing one vertex in it makes the refinement discrete.

All five branches give the same ordered adjacency certificate because that clock acts transitively on the selected cell. The discrete refinement shows that a rooted automorphism fixing the individualized vertex is trivial; hence the rooted stabilizer has exactly5 elements. Translation transitivity then gives full automorphism order625.

The checker compares **full adjacency byte strings**, not merely hashes, and finds408 distinct certificates. Any unrooted isomorphism can be followed by a target translation to fix0. Therefore these are exactly408 prime graph types. All composite states have the established type T5[T5[T5]], giving one further graph type, with automorphism order5^31.

Thus this specified family has **409 full graph types at order125**. It is not the full universe of125-vertex tournaments.

The same canonical composition theorem gives exactly817 composite graph types in the625-vertex bank:408 of T5[prime125],408 of prime125[T5], and one fully coherent T5 tower. Their full automorphism orders are5^21,5^129,5^156 respectively. Prime625 graph isomorphism remains a separate problem.

## 3. Nilpotent automorphisms change the625 normalizer

In R=Z[delta]/(25,delta^2-5), the tempting assumption that only delta->plus-or-minus delta survives is false. The map

    delta -> 6 delta

is an additional automorphism. In fact all ring automorphisms have delta->b delta with b=plus-or-minus1 modulo5, giving10 choices.

To see completeness, write an image as a+b delta. Invertibility and reduction modulo5 force b to be a unit and5 to divide a. Squaring and comparing coefficients then forces a=0 modulo25 and b^2=1 modulo5. The regular multiplication algebra is generated by the named clock, so a clock normalizer is a unit multiplication followed by one of these algebra automorphisms. All preserve its norm-one principal-unit subgroup H of order25. There are500 units, so the normalizer has order5000.

The [exact signed-action checker](portable_checks/normalizers.js) computes its200 effective actions on24 orientation bits. Their fixed-state histogram is

    175 actions fix0,
    16 actions fix256,
    8 actions fix65536,
    the identity fixes16777216.

It also checks all eight ideal-cut intersections for every action. Only the identity effective action fixes prime states. Burnside therefore gives86528 normalizer orbits in the full bank and83232 in its prime part.

### Why83232 also counts full affine prime classes

A prime state's normalizer stabilizer is exactly H. Its full additive stabilizer K lies in GL2(Z/25), whose odd part has order3*5^5. Tournament automorphism groups have odd order. A Sylow5 subgroup containing H cannot be larger: in a finite p-group every proper subgroup is properly contained in its normalizer, which would contradict the computed normalizer stabilizer. (The normalizer condition follows by the subgroup action on cosets: the number of fixed cosets is positive and divisible by p.)

Thus the Sylow5 subgroup is H. Its number of conjugates divides3 and is1 modulo5, hence is1. Therefore H is normal in K, so K lies in the clock normalizer and K=H. Any additive equivalence between prime states must normalize this same H.

Consequently there are **83232 full affine prime classes**, each with affine automorphism order625*25=15625. This does not assert83232 full graph types.

At order125 the composite states already show the distinction:64 labels give8 normalizer classes,4 full additive-affine classes, and one graph class. The same checker verifies all2000 additive maps and both orbit partitions.

## 4. A hostile at the next depth prevents extrapolation

At uniformizer depth5, on3125 vertices, scalar multiplication by26 has order5 and lies outside H: its norm is51 modulo125, not1. Nevertheless it fixes prime states.

The [depth-five certificate](portable_checks/depth5_hostile.js) checks the74-bit bank, every fixed-cut intersection, and an explicit prime mask. It finds

    all scalar-fixed states =17179869184,
    prime scalar-fixed states =17045913600,
    explicit prime mask =0x80000801000.

The displayed mask fails all four proper ideal-cut tests. Thus small-depth affine rigidity does not continue automatically. The missing coordinate was an actual nilpotent automorphism/clock lift, not a numerical fluctuation.

## Reproduction and boundaries

Run each linked file with Node.js. No packages are needed. Every mathematical check throws explicitly on failure. All floating-point-number calculations in these four programs are within their exact integer range; the larger counts and masks use BigInt.

The first prime-factor checkpoint and the independent audits supply the all-depth ideal-module criterion used by the signed-cut computations. The full graph classification has its own adjacency certificate and does not rely on an unproved identification of affine and graph isomorphism.

Further enumeration, spectral and arithmetic work continues. In particular, no claim here completes unrestricted A000568 enumeration or proves a major number-theory conjecture.
