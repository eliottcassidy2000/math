# Short planar interfaces: exact Hamiltonian probes and the six-port parity boundary

**Session 190, 4 November 2026. Status: PROVED (elementary derivations below), FINITE-EXACT (reproducibility packet), CITED (external existence theorems). Not a literature-priority claim.** These results specialize standard matching-connectivity / meander methods to the small boundaries relevant to planar coloring reductions. They do not transfer coloring reducibility to Hamiltonicity, prove a near-linear Hamiltonian algorithm, or settle an open all-order square-sum problem.

## 1. Current literature and inheritance

The October 1, 2026 revision of Inoue–Kawarabayashi–Miyashita–Mohar–Thomassen–Thorup is stronger than the March abstract formulation: Theorem 3.3 produces linearly many non-touching separating 3- or 4-cycles or reducible configurations; Corollary 3.4 instead makes the configurations induced and permits obstructing cycles. Their configuration rings can have 18 vertices. Obstructing cycles mean separating triangles or quadrilaterals, or pentagons having at least two vertices on each side. These are 4-coloring reductions, not Hamiltonicity reductions. [Primary v3](https://arxiv.org/html/2603.24880v3), Sections 2–3. This packet uses the geometric boundaries, not the coloring extension claim.

The established algorithmic object is the matching-connectivity matrix: its entry records whether two perfect matchings form one cycle. Cygan–Kratsch–Nederlof use its modulo-two rank and permutation submatrices. [Primary paper](https://arxiv.org/abs/1211.1506). Noncrossing pairings with component weights are the meander / Temperley–Lieb framework of Di Francesco–Golinelli–Guitter. [Primary paper](https://arxiv.org/abs/hep-th/9602025). A newly checked October 8, 2026 paper by Caroppo–Da Lozzo–Di Battista–Vihrovs explicitly links sphere-cut dynamic programming to meanders; it claims improved subexponential, not near-linear, algorithms for planar Hamiltonian cycle. [Primary paper](https://arxiv.org/abs/2610.11944). Our local formulas are independently proved here, rather than claimed as a new general framework.

Targeted repository inheritance:

- `session189-publication/S189_final.md`: exact named response versus lost marginals, physical probes, and first dyadic obstruction are useful methodological comparisons, not dependencies on tournament theorems.
- `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md` and `THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md`: square-sum interval graphs and corrected endpoint discipline.
- `session32-square-sum/README.md` and `verify_square_sum.py`: exact endpoint-marked chain kernels; the hostile Q_23 path ending at its degree-two vertex 22.
- `07-reflections/hypohamiltonian-criticality-and-the-lonely-runner.md`: historical motivation only. No lonely-runner conclusion is inherited.

Closest proved mechanism: path connectivity across a separator. Canonical hostile: two pairings of four ports with the same degree vector. Corrected near miss: graph degree two does not force both edges in a Hamiltonian path unless the vertex is certified not to be an endpoint. Least-used sidecar: the entire boundary pairing and its deletion-indexed profiles.

## 2. Precise setup

Let A and B be simple graphs drawn in opposite closed disks on the sphere. They intersect exactly in k distinct, cyclically ordered boundary vertices R=(r_0,...,r_(k-1)), have disjoint edge sets, and each has at least one vertex outside R. Edges along a separating cycle are assigned to exactly one side. The graphs need not themselves contain the boundary cycle. There are no crossings on either side.

A *partial Hamiltonian forest* F in A is a spanning edge subset such that:

1. every vertex outside R has degree two;
2. each boundary vertex has degree zero, one, or two;
3. F has no cycle.

Isolated boundary vertices are allowed. Every positive-edge component is a path with two distinct degree-one boundary endpoints. Write d(F)=(d_0,...,d_(k-1)); let M(F) pair the two endpoints of each path. M is a noncrossing perfect matching of the nonempty, even set D={i:d_i=1}. Noncrossing follows because vertex-disjoint paths in a disk cannot join alternating boundary endpoints. Conversely all formal states (d,M) are realizable, as proved in Section 5.

Let f_A(d,M) count such edge subsets. Counts are unoriented and paths are not ordered. Since A has a nonempty interior, a forest with D empty is impossible: its positive-edge components would be cycles. This excludes the 2^k formal degree-only states with no ones. These exclusions are essential; a cycle already completed inside one side cannot later cover the other side's interior.

## 3. Exact gluing theorem

**Theorem 1 (all k).** The number h(A union B) of unoriented Hamiltonian cycles is

    h(A union B) = sum_(d,M),(e,N) f_A(d,M) f_B(e,N) K((d,M),(e,N)),

where K=1 exactly when d_i+e_i=2 at every boundary position and the multigraph M union N on D is connected. Matching edges from opposite sides remain distinguishable: two equal matching edges on two endpoints form one alternating two-edge cycle.

**Proof.** A Hamiltonian cycle restricts on each side to a partial forest. A closed component in one restriction would already be the entire global cycle, impossible because the other side has an interior vertex. Its degrees add to two and contracting each positive-edge path gives M union N. The original union is connected iff this contracted multigraph is connected. Conversely compatible forests glue to a connected spanning 2-regular graph, hence a Hamiltonian cycle. Edge disjointness makes restriction and union inverse maps, with no orientation or factorial factor. ∎

The usual noncrossing matching count gives

    S_k = sum_(m>=0) binom(k,2m) Catalan(m) 2^(k-2m) = Catalan(k+1).
    L_k = S_k - 2^k.

For completeness, the weighted-Motzkin generating function F=1+2xF+x^2F^2 equals (C(x)-1)/x, where C=1+xC^2 is the Catalan generating function. Thus S_k=C_(k+1). L_k is the number of live states, not the count of nonzero entries of every particular patch.

At k=2,3,4,5,6, the live counts are respectively 1,6,26,100,365; ambient counts are 5,14,42,132,429.

## 4. Unique exterior completion through five ports

**Theorem 2.** For 2<=k<=5, K is a permutation matrix on the L_k live states, and is symmetric. Define the involution tau as follows:

- replace d by 2-d;
- if D has two members, keep its sole matching;
- if D has four members, exchange its two noncrossing matchings.

Then

    h(A union B) = sum_s f_A(s) f_B(tau(s)).

**Proof.** D has either two or four members. At two, the two path components glue to one cycle. At four, write the cyclic order as a,b,c,d. The possible noncrossing matchings are ab|cd and ad|bc. Equal matchings form two components; unequal matchings form one. The degree complement is forced. ∎

**Consequences.** For separating triangles only six interface coordinates can be live; for quadrilaterals there are 26; for pentagons there are 100. These are exact arbitrary-patch interface dimensions, not a claim that evaluating one large patch is easy. The v3 3/4-cycle outcome therefore exposes a 6/26-coordinate Hamiltonian interface, but says nothing about how expensive its interior is to evaluate, and the alternative reducible-configuration rings can be much larger.

**Explicit hostile to forgetting the matching.** Four ports with d=(1,1,1,1) and matching 01|23 versus 03|12 have identical degrees and identical total partial-forest count. Realize each paired connection as a length-two path through its own private interior vertex. Glue against another copy of 01|23. The first gives two disjoint 4-cycles and no Hamiltonian cycle; the second gives one 8-cycle. The verifier checks this physically. Degrees, scalar path counts, and coloring extension data cannot replace the pairing coordinate.

## 5. Realizable unit probes and exact sufficiency

**Theorem 3.** Every live formal state s=(d,M), for any k, is the entire profile of a simple planar disk graph P_s: f_(P_s)=e_s. It can be constructed with at most 2k-1 total vertices and with no edge joining two boundary vertices. For k<=5, gluing P_(tau(s)) to A returns exactly f_A(s). For k<=5, two patches with the same ordered boundary have equal Hamiltonian-cycle counts in every planar disk exterior with nonempty interior if and only if their full profiles agree. The live-state connection matrix has rank L_k over every field for k<=5. Any factorization of the formal response operator through a d-dimensional vector space has d>=L_k; in particular any linear/bilinear gluing representation needs at least L_k coordinates.

**Construction and proof.** List the degree-one boundary positions cyclically. To each such endpoint u assign all degree-two boundary positions encountered after u and before the next endpoint, skipping degree-zero positions. Join u and its assigned vertices in a path in a thin boundary collar. For each pair of M, join the tips of its two collar paths by a noncrossing interior arc. The tips retain the endpoints' cyclic order, so these arcs can be drawn disjointly. Subdivide every edge just constructed once, using different private vertices. Each private vertex has graph degree two and is required to have selected degree two, so every edge is forced. The resulting profile is exactly the desired single state. If |D|=2m, there are (#degree-two positions)+m unsubdivided edges, at most k-m<=k-1; thus the total vertex bound is 2k-1. Degree-zero ports remain isolated on this side. At least one private vertex exists. No boundary-boundary edge survives subdivision, so probes never duplicate a boundary edge of A.

Apply Theorem 2 to obtain coefficient recovery. Equality of profiles clearly suffices by Theorem 1, and the probes prove necessity. The unit profiles and the permutation K also give the rank / formal-response-factorization lower bound. ∎

The universal exterior class here is sparse disk graphs. The theorem does NOT assert that each coordinate has an isolating exterior that is itself a triangulation. It also does not give a compression lower bound for arbitrary encodings of integer vectors, including a real linear encoding followed by an unrestricted nonlinear decoder. Irrationally independent real weights can encode an integer vector in one real coordinate; the rank claim instead concerns factorization of the formal linear response operator.

## 6. The first parity-blind direction is at six ports

All degree patterns with at most four ones retain their unique mate. At k=6 the all-ones pattern has the following five matchings, in this order:

    M0 = 01|23|45
    M1 = 01|25|34
    M2 = 03|12|45
    M3 = 05|12|34
    M4 = 05|14|23.

Its K block is

    0 0 0 1 0
    0 0 1 0 1
    0 1 0 0 1
    1 0 0 0 0
    0 1 1 0 0.

It is the adjacency matrix of K_2 on {M0,M3} plus K_3 on {M1,M2,M4}. Hence its determinant is -2, its rational rank is five, and its F_2 rank is four. The full 365-state K has rational rank 365 and F_2 rank 364; its determinant has absolute value two. There are 362 rows with one mate and exactly three rows with two mates.

The vector e_(M1)+e_(M2)+e_(M4) lies in the modulo-two kernel. Every exterior therefore sees an even count from this *formal sum* of three unit probes. This is a linear-algebraic blindness statement, not an assertion that a single connected graph realizes that sum as its profile. Integer counts still see it. It is also not the false claim that parity alone detects Hamiltonicity through five ports: an individual nonzero integer entry can already be even. What changes at six is injectivity of the entire exterior-response operator acting on profile vectors modulo two.

This supplies a precise comparison with S189: the object, map, and lost coordinate are explicit. Both use full named interfaces and distinguish integral information from a modular observer. It proves no theorem about the S189 tournament response.

## 7. Exact deletion-criticality interface

Suppose both A and B have at least two interior vertices, so deleting one vertex never empties either interior. For an interior vertex v of A, compute f_(A-v) on the same boundary; for an interior vertex of B use f_(B-v). For a boundary vertex r, delete r from *both* sides and use the shorter cyclic order on R minus {r}. Theorem 1 applies to every deletion.

**Corollary 4.** For k<=5, A union B is hypohamiltonian iff

    sum_s f_A(s) f_B(tau(s)) = 0,

and each corresponding sum after every one-vertex deletion is positive. Because entries are nonnegative integers, the first condition says the two supports are disjoint after applying tau, and every deletion restores a compatible pair. The changed profile must be recomputed; it is not inferred from the undeleted one.

If a side has just one interior vertex, deletion can empty it and an already closed cycle may become valid. The present cycle-free state definition must then be augmented by a completed-cycle state or the exceptional deletion tested directly. Omitting this case would invalidate the criterion.

For hypotraceability and arbitrary graphs G of order at least two, one can adjoin a universal vertex z. Unoriented Hamiltonian paths of G correspond bijectively to Hamiltonian cycles of G+z. If G is separated by R, put z in both sides, assign each z-r edge once, and the separator becomes R union {z}. This graph need not remain planar. Therefore unrestricted matchings, not Catalan noncrossing matchings, are required: at original |R|<=5 the largest augmented alphabet has

    sum_m binom(6,2m) (2m-1)!! 2^(6-2m) = 499

ambient states, of which 435 are live under the same nonempty-interior hypotheses. This is a correct finite interface for path/deletion experiments, not a 100-state planar claim for square-sum graphs.

## 8. Square-sum graphs: a useful negative bridge

Let Q_n have vertices 1,...,n and edges xy when x!=y and x+y is a square. Every hypohamiltonian graph has minimum degree at least three: deleting a neighbor of a vertex of degree<=2 leaves that vertex with degree<=1, incompatible with a spanning cycle. Every hypotraceable graph of order at least three also has minimum degree at least three. Degree zero or one immediately violates a deletion card. If v has degree two with neighbors a,b, any Hamiltonian path of G-a must end at v; appending a at v would give a path of G, a contradiction.

For all 2<=n<=31, Q_n has a vertex of degree at most two. Explicit witnesses are vertex 1 for n<=3, vertex 4 for 4<=n<=7, vertex 8 for 8<=n<=15, and vertex 16 for 16<=n<=31. In the last band vertex 16 can only use square sums 25 and 36. The other bands are equally immediate. Thus none of these graphs is hypohamiltonian, and none with n>=3 is hypotraceable.

Gerbicz's own January 21, 2018 entry records Hamiltonian paths for every n>=25 and cycles for every n>=32. [Author's statement at OEIS A090461](https://oeis.org/A090461). We cite, but do not independently verify, the infinite construction and its finite base/glue tables. The linked Mersenneforum page returned HTTP 403 in this session. **Using that cited theorem**, no square-sum interval graph Q_n, n>=3, is hypohamiltonian or hypotraceable. There is a genuine convention-sensitive exception: Q_2 is two isolated vertices, so if a singleton is traceable and disconnected two-vertex graphs are allowed, Q_2 satisfies the raw hypotraceable definition. It is excluded only by an explicit nontrivial/connected convention, not by calling it disconnected. This rules out treating the interval failures as hypo-critical examples.

A new bounded deletion experiment gives finer information, independently of the cited infinite theorem. The nontraceable interval graphs below have exactly these vertex deletions that restore traceability:

    Q18, Q19, Q20, Q21: {18}
    Q22: {2,18}
    Q24: {11,18,22,24}.

Thus Q24 is not even close to deletion-universal: 20 of its 24 deletion cards still fail. Full positive path witnesses and all negative decisions are in RESULTS.json. The underlying useful research target is now a prescribed-endpoint / deletion-profile extension of the existing square-sum construction, not rediscovery of ordinary path existence.

## 9. An infinite hostile: many reductions with robust Hamiltonian failure

**Theorem 5 (complete face-stellation).** Let T be any n-vertex simple plane triangulation, n>=4. Insert a new vertex in each face and join it to that face's three corners. Call the resulting triangulation K(T). It has N=3n-4 vertices and 2n-4 pairwise non-touching degree-three configurations. Nevertheless:

- for n>=5, K(T) is non-Hamiltonian;
- for n>=6, K(T) is nontraceable and every one-vertex-deleted card is non-Hamiltonian;
- for n>=7, every one-vertex-deleted card is nontraceable as well.

**Proof.** Euler's formula gives 2n-4 faces. The new vertices form an independent set I of size 2n-4; deleting the n original vertices leaves exactly that many singleton components. A Hamiltonian cycle needs |I|<=n, and a Hamiltonian path needs |I|<=n+1. These fail at n>=5 and n>=6, respectively. After deleting an old vertex there remain 2n-4 independent new vertices but only n-1 old vertices. After deleting a new vertex there remain 2n-5 new vertices and n old vertices. Apply the same cycle/path capacity bound to obtain the stated thresholds. Every new vertex has degree three, distinct new vertices are nonadjacent, and each neighbor triple is a triangle. Such a vertex is a 4-color-reducible configuration: any proper coloring of its neighbors uses three different colors and extends with the fourth. ∎

The exact abundance is (2N-4)/3 non-touching reducible configurations. Hence even a positive linear density of these reductions does not imply traceability, Hamiltonicity, or proximity to a deletion-critical graph. This hostile changes the proposed bridge from an existence transfer to an *interface-preservation obligation*.

At n=5 take the triangular bipyramid with base triangle 0,1,2 and apices 3,4. Its six faces in order are

    (3,0,1), (3,1,2), (3,2,0), (4,0,1), (4,1,2), (4,2,0).

Insert vertices 5,...,10 in those faces. The resulting 11-vertex graph is non-Hamiltonian by the six-versus-five independent-set bound. But after deleting the degree-three reducible vertex 10 it has Hamiltonian cycle

    0,5,3,6,1,8,4,9,2,7,0.

Thus even one legitimate coloring reduction can change the yes/no answer to Hamiltonicity. The script `stellation.py` independently enumerates cycles, verifies this witness, and records all deletion counts. Its six new-vertex deletions restore Hamiltonicity, while its five old-vertex deletions do not. This is a concrete non-hypohamiltonian partial-deletion-critical control, not a newly claimed graph construction.

## 10. Reproducibility, scope, and next decisive experiment

Run:

    python3 verify.py > RESULTS.json
    python3 -O verify.py > RESULTS.optimized.json
    cmp RESULTS.json RESULTS.optimized.json
    python3 independent_kernel.py > KERNEL_REPLAY.json
    python3 stellation.py > STELLATION.json

No installs, repository mutation, floating-point calculation, or external library is used. Explicit exceptions, not Python assertions, guard all tests.

Finite universes:

- all 498 live formal states for k=2,...,6, with a physically constructed unit probe independently enumerated by edge subsets;
- all compatibility entries for those state alphabets;
- 189 physical gluing tests: k=3,4,5; seven fixed full/edge-deleted two-interior-vertex wheel patches per k; original graph and every vertex deletion, checked by both forest edge-subsets and global vertex-order cycle enumeration;
- the four-port degree-only hostile;
- all degree witnesses through n=31;
- all single-vertex deletion cards of Q18,Q19,Q20,Q21,Q22,Q24, with direct path search and retained witnesses.

The six-port matrix is derived in full above, so the computation is not its sole proof. The square-sum negative deletion decisions are FINITE-EXACT, not human cut certificates; an independent endpoint-kernel replay is supplied separately.

**Most useful next test.** Choose one concrete reducible configuration from the 2026 list, retain its actual ordered ring and all endpoint pairings, and test a proposed smaller replacement for equality of the full Hamiltonian profile. A profile mismatch gives a physical exterior counterexample using the probes when the ring has at most five vertices. On larger rings Theorem 1 still applies, but the uniqueness and parity-unimodularity theorem no longer does. A successful local replacement must be established separately before attempting simultaneous replacements. Many non-touching configurations alone do not make their exterior Hamiltonian requirements independent.
