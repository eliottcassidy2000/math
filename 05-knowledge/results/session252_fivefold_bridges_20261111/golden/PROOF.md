# Five cores, an exact associahedral interface, and a finite/affine ADE split

2026-11-11. New elementary identities below are **PROVED IN THIS NOTE**. Small controls are **FINITE-EXACT**. The five-core all-height covering assertion and Chan–Pak remain **OPEN**. This note does not identify five arbitrary cores with five Platonic solids.

## Results first

1. The Session251 nine-coordinate signature is the pair of evaluations at ±1 of a **five-coordinate upper-triangular path-interface algebra**. Ordinal gluing is exactly matrix multiplication. Its reassociations give the genuine five-vertex Catalan/A2 associahedral pentagon.
2. Directed-triangle **arc-overlap graphs** of the five actual cores are respectively empty, A1, A2, A3, and affine A2. In particular the two H9 cores are separated by a finite-versus-affine Dynkin boundary although their odd-cycle conflict graphs and H agree.
3. The proposed D5 birth fork is only a selected spanning skeleton: its actual padded-parent birth graph has eight edges, not four.
4. The binary g transfer has an exact digit sum and Mahler equation; its triangular dyadic matrices do not secretly carry golden-ratio eigenvalues. Golden geometry enters the reassociation pentagon honestly, rather than by relabelling the five seeds.

Inheritance: Session251 `ternary/PROOF.md` and `compiler.py`; Session250 rooted gluing; THM-002 OCF. Canonical hostile: scalar H loses marked path endpoints, and both H9 cores have equal uncoloured path evaluations. Closest corrected near miss: historical Platonic correspondence equating a six-dimensional tournament cube's polar to a three-dimensional octahedron. Underused sidecar: actual triangle overlap multiplicities and polynomial path-interface matrices.

Live board: core birth grammar / endpoint ownership / ordinal interface / Catalan reassociation / triangle overlap / digit registers.

## 1. Polynomial transfer matrix and its pentagon

For a coloured tournament C, use precisely Session251's d=H(C) and polynomials K,I,O,J, with path weight H(C minus path vertices), exponent the number of path vertices, O marking starts and I marking ends. Define

    T_C(z) = [[d, O(z), J(z)],
              [0, d+K(z), I(z)],
              [0, 0, d]].

For disjoint coloured tournaments C,D with every arc C→D,

    T_(C→D)(z) = T_C(z) T_D(z).                         (T)

Proof: H multiplies under ordinal join. Every directed path meets each part in at most one segment. A path may start in the first part, the second part, or meet both. Thus

    K_CD = K_C d_D + d_C K_D + K_C K_D,
    O_CD = O_C(d_D+K_D) + d_C O_D,
    I_CD = d_D I_C + (d_C+K_C)I_D,
    J_CD = J_C d_D + d_C J_D + O_C I_D.

Unused vertices in each part contribute the product of residual H values. These four equations and d_CD=d_C d_D are exactly (T). They are polynomial identities at all z, not just ±1.

For four ordered blocks A,B,C,D there are exactly five binary parenthesizations:

    (((AB)C)D), ((A(BC))D), (A((BC)D)), (A(B(CD))), ((AB)(CD)).

Elementary rotations connect them in that cyclic order, giving the associahedron K4 (the A2 cluster pentagon). Each vertex evaluates to the identical marked interface T_A T_B T_C T_D. This is a real map from decomposition trees to interfaces, not an identification of the five cores with those trees. The map preserves H and the path polynomial endpoints, and forgets decomposition history. Actual ordered coloured blocks are the restoration sidecar.

**Golden geometry with a specified source.** The flip graph is C5. Its adjacency eigenvalues are 2, (sqrt(5)-1)/2 twice, and -(sqrt(5)+1)/2 twice. Its planar first Fourier mode gives the regular pentagon, whose diagonal/side ratio is phi. Thus the golden triangle is a faithful geometric representation of the *reassociation relation*, not of the seed bank or H. Polygon geometry forgets every coefficient and core marking; matrix labels restore them. An exact interface implementation has zero discrepancy around this pentagon. Coherence alone is insufficient: dropping O_C I_D can produce a different still-associative algebra. Therefore each chart must also agree with direct marked-path evaluation; the pentagon tests consistency, not correctness by itself.

Possible proof use (**OPEN**): build divide-and-conquer inverse certificates with tree-rotation normalization; prove a finite family of marked interface rewrite rules preserves the target image. Associativity supplies coherence automatically, but does not supply the arithmetic covering invariant. This is the right locus to seek a five-chart proof atlas.

## 2. A genuinely non-cosmetic ADE fingerprint of the core bank

For a tournament C, form Γ_arc(C): vertices are its directed triangles, and two are adjacent iff they share a directed arc. This is different from the OCF conflict graph Ω(C), whose vertices include *all* directed odd cycles and whose adjacency means sharing any vertex.

The exact core data are:

- K1: no directed triangles, Γ_arc empty.
- C3=(2,4,1): triangle 012, Γ_arc=A1.
- S4=(2,4,9,3): triangles 012,123, Γ_arc=A2.
- S5a=(2,4,9,19,7): triangles 012,123,234 and pentagon 01234. Γ_arc=P3=A3. The outer triangles meet at vertex2, but share no arc.
- S5b=(2,4,25,3,11): triangles 012,123,124 and pentagon 01243. All three triangles share arc1→2; Γ_arc=K3, the underlying graph of affine A2.

Every two odd cycles in a tournament of order≤5 intersect. Therefore Ω(C) is complete on respectively 0,1,2,4,4 vertices, and

    I(Ω(C),x) = 1 + (0,1,2,4,4) x.

At x=2 these give the five core H values 1,3,5,9,9. At x=1 they give 1,2,3,5,5; this small Fibonacci prefix is an evaluation fact, not a Fibonacci recurrence or a bijection with canonical Fibonacci representations. The actual mechanism is low-order pairwise cycle intersection.

The two H9 cores therefore have identical complete odd-cycle conflict graphs K4. Yet their triangle arc-overlap adjacency spectra are

    S5a: sqrt(2),0,-sqrt(2),
    S5b: 2,-1,-1.

Equivalently, 2I-A(Γ_arc) is positive definite for S5a (determinant4), while for S5b it is positive semidefinite with null vector (1,1,1) and determinant0. This is exactly the elementary finite/affine ADE distinction. It is a natural, explicitly constructed invariant rather than attaching Dynkin names to a five-element list.

The local distinction has an endpoint meaning: S5b has a common directed-arc bottleneck for its three triangle excursions, whereas S5a has a chain with only a single-vertex overlap between its outer excursions. This is a promising explanation of why coloured attachment languages differ while scalar H does not. It is **not yet a proof** that this invariant alone determines the nine-coordinate cut language.

Cheap decisive next test: among more cores, compare full cut-signature languages after grouping by triangle arc-overlap graph and triangle-to-pentagon incidence. A collision with unequal languages reveals the next missing coordinate. A useful stronger retained structure is the cycle-by-arc incidence matrix, with rows labelled by complete directed cycles and columns by directed arcs; its row Gram records overlap multiplicity and does not collapse the common-arc versus chain configurations.

## 3. The D5 fork is a skeleton, not the induced birth relation

Deleting each vertex of the literal cores gives H lists:

    C3 : (1,1,1)
    S4 : (3,1,1,3)
    S5a: (5,3,1,3,5)
    S5b: (5,1,1,5,5).

At these sizes, H1 means transitive, H3 means padded C3, H5 means S4. Therefore the directed padded-parent relation contains

    K1→C3,S4,S5a,S5b;
    C3→S4,S5a;
    S4→S5a,S5b.

These eight edges include the aesthetically appealing D5 tree K1–C3–S4 with fork S5a,S5b. Selecting only that tree discards four actual birth operations. The full graph is K5 minus edges C3–S5b and S5a–S5b. There is no induced D5 classification here. This hostile probe leaves the arc-overlap ADE fingerprint intact.

## 4. The binary transfer's real arithmetic

Write n=sum b_i 2^i. Iterating Session251's recurrence gives exactly

    g(n)=1+sum_(i:b_i=0) [floor(n/2^(i+1))+b_(i+1)]
        =1+2 sum_(i:b_i=0) [floor(n/2^(i+2))+b_(i+1)],

where only positions below n's bit length need appear, and omitted high bits are0. In the second line the first line's floor is split into twice the shorter prefix plus b_(i+1). This is a nonnegative digit-pair statistic, immediately showing g odd.

For G(z)=sum_(n≥0) g(n)z^n, the exact formal Mahler equation is

    G(z)=(1+z)G(z²)+z²/(1-z²)²+z²/(1-z⁴).

The coefficient pullbacks from Session251 act on (a,b,c,e) by

    D0(a,b,c,e)=(a,a+2b,a,e),
    D1(a,b,c,e)=(a,2b,0,b+c+e).

Both are lower-triangular linear maps in this ordering, with diagonal (1,2,0,1). A word of length k has diagonal (1,2^k,0,1). In particular no golden eigenvalue or faithful order-five rotation is hiding in these *unaugmented* linear digit maps. A literal linear conjugacy to Fibonacci transfer is ruled out by spectra. Fibonacci/golden mechanisms require a different state space or constrained digit-language count, not an unsupported renaming.

Possible new route (**OPEN**): use the digit-pair expression to organize targets by weighted zero/one incidences, and retain the core's arc-overlap incidence as the missing combinatorial sidecar. Seek a controlled exchange of zero-blocks between A and B that preserves the bilinear AB term while moving the boundary coefficients between the two H9 charts. This is a specific arithmetic surgery problem, not an assumption that a finite-state automaton already exists.

## 5. Historical recovery and use boundaries

- `07-reflections/the-five-platonic-tournaments.md`: useful early ambition, with an explicit current correction at its top. The six arc coordinates of K4 give a six-dimensional tournament cube, not a three-dimensional cube; L(K4) is the valid octahedral graph. Equal f-vectors do not identify the irregular five-vertex tournament class graph with the icosahedron.
- `05-knowledge/results/platonic_interchange_s204.out`: older concrete tetrahedral/cubic/interchange computations. Retain the specific maps, not its concluding analogy as a proof.
- `05-knowledge/results/fivefold_boundary_20261018g.md`, section3: genuine A4→two complex Fourier coordinates, with fivefold action and golden conjugate window. Planar projection loses the bounded-conjugate window that certifies discreteness. It could host a five-chart count vector, but cyclic chart order has not been constructed, so no core-equivariance is claimed.
- `07-reflections/zeckendorf-non-consecutivity-pairing.md`: old OCF/independence-fugacity program. Its useful primitive is I(Ω,2), contrasted with golden path independence at fugacity1. The five-core computation above makes the precise low-order survivor explicit; broad historical realizability assertions are not imported.
- `05-knowledge/results/zeckendorf_guard_automaton_20261003.md`: retained two-register Fibonacci arithmetic, a potentially usable target representation for future inverse search, with actual guards and markers required.
- `05-knowledge/results/three_sixties_mechanism_boxeph_S156.out`: older Fibonacci-mod5/A5 discussion, but equal order60 is not a map from the five-core grammar. No new A5 symmetry is claimed here.

Reproduction: `python3 golden/probe.py` from this packet (relative Session251 compiler import). Output `golden/results.txt` records all triangle/pentagon lists and deletion counts, 200 exact ±1 ordinal matrix composition checks, and 4096 digit-formula checks. No large inherited census was rerun. The proof of (T) is independent of these checks.


## 6. The long stem: a genuine golden/Fermat carrier, recovered from the tower

Define T_n on 0,...,n−1 by i→i+1 and j→i for j≥i+2. Then T3=C3, T4=S4 and T5=S5a literally in the supplied labels. T1=K1; T2 is transitive padding.

Every directed cycle is a consecutive interval, traversed increasingly and then closed from its maximum to its minimum. Proof: start at the cycle's minimum a. Its only out-neighbour larger than a is a+1. At each subsequent vertex b, all smaller vertices ≥a already lie on the path, so the cycle either closes to a or continues to b+1. No skipped vertex or differently ordered cycle is possible. Thus the triangles are exactly consecutive triples, and

    Γ_arc(T_n)=P_(n−2)=A_(n−2), n≥3.

The two triangles at distance2 share one vertex but not an arc. Other distinct non-neighbour triangles share no arc. Consequently

    ρ(Γ_arc(T_n))=2 cos(π/(n−1)).

At n=6 this is φ. The positive eigenvector of P4 is proportional to (1,φ,φ,1), directly encoding a golden ratio between the central and end triangle coordinates. This is a genuine golden observable on the next member of the five-core stem, not an externally imposed pentagon.

For the cycle-packing polynomial Z_n(x)=I(Ω(T_n),x), decomposition at the first vertex yields

    Z_n=Z_(n−1)+x sum_(odd ℓ≥3,ℓ≤n) Z_(n−ℓ),
    sum_(n≥0) Z_n(x)z^n=(1−z²)/(1−z−z²−(x−1)z³).

Hence Z_n=Z_(n−1)+Z_(n−2)+(x−1)Z_(n−3) from n≥3, with Z0=Z1=Z2=1. At x=1, Z_n(1)=F_n for n≥1; at x=2, H(T_n) is tribonacci. These are two evaluations of the *same actual cycle incidence structure*. Counting unweighted cycle packings has golden growth; weighting each chosen cycle twice gives tournament Hamiltonian counts and tribonacci growth. This is more useful and honest than pretending H itself is Fibonacci.

**Recovered, not new:** precisely this polynomial recurrence is THM-4952(C), `01-canon/theorems/THM-4952-ribbon-odd-cycle-packing-polynomial-is-finite-state-base-2-is-the-unique-positive-rank-drop-x-bonacci-tower-base-3-ribbon-words-miss-exactly-10-and-31-and-64-has-no-strong-realiser.md`. The rooted generalization is THM-4865, `01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci-h-k-equals-h-k-1-plus-h-k-2-plus-h-k-3-for-every-rooted-base-the-lineage-c2-law-at-any-vertex-and-the-three-root-response-numerations-attain-every-odd-number-in-69-702189.md`. That theorem traces the special tower back to March5 and THM-4731(6). We use these identities as inherited mechanisms, with the elementary interval argument supplied to expose the map.

The new synthesis is:

    K1,C3,S4,S5a = the initial non-padding stem;
    S5b = the affine-A2 arc-bottleneck alternative;
    T6 = next finite A4 stem, golden spectral radius φ, H=17.

Indeed direct H counts for n=1,...,13 give 1,1,3,5,9,17,31,57,105,193,355,653,1201. Thus the Fermat prime17 and golden ratio occur as distinct rigorously specified observables of one canonical tournament T6. It does not follow that Fermat primes classify the stem; its later values quickly leave that family. The natural Chan–Pak question is whether the affine branch supplies the cuts missing from the finite-type stem's arithmetic images. The two H9 holes21901 and86641 make that question testable without extending any census.


## 7. A direct icosahedral incidence marker, not one solid per core

The older `07-reflections/the-icosahedral-fifteen-s6.md` supplies a canonical
five-label dictionary. Two 3-subsets sharing exactly one label v partition
the other four labels into two unordered pairs {a,b},{c,d}. Send that pair
of triples to the double transposition (ab)(cd), fixing v. There are
5·3=15 such objects, exactly the fifteen involutions in A5. In its
icosahedral rotation action these are the fifteen half-turn axes. The
incidence map is exact and relabelling-equivariant; a spatial realization
of A5 is needed only to draw the axes.

Apply this dictionary only to pairs of actual directed triangles of a
five-vertex core. S5a has one such pair, 012 and234, so its marked axis is
(01)(34), fixing2. S5b has none: every triangle pair shares the full
arc1→2. Thus the finite/affine overlap split also has a genuine Platonic
incidence shadow, one icosahedral axis versus no such axis. The map forgets
edge directions and all endpoint weights, so it is not a complete
attachment invariant or an icosahedral symmetry action on the core.

There is also a direct golden-triangle realization of the T6 Perron
vector: its values 1,φ,φ,1 assign short weight1 to the two end triangles
and long weightφ to the two interior triangles. The side lengths
(φ,φ,1) form the 36°–72°–72° golden triangle. Folding the A4 path by
reversal gives the two-weight equations λa=b and λb=a+b, hence
λ²=λ+1. This pins the source, quotient, and preserved ratio; the plane
triangle by itself forgets the four labelled triangle positions.
