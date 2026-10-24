# G27: exact birth response, not a new Hamiltonicity threshold

Status: **FINITE-EXACT** for the G24..G28 census below; **PROVED (elementary; no priority claim)** for the graph-birth identities and endpoint-pair transfer. The degree-two path-deletion bijection itself is **INHERITED**, already explicit in Session32. No Collatz convergence result or cross-domain equivalence is asserted.

## 1. The mechanism at27

Let G_N be the simple undirected graph with vertices1..N and edges between distinct vertices whose sum is a square. Paths are unoriented, counted up to reversal. A path cover is an unordered collection of disjoint paths partitioning all vertices; singletons are allowed.

The new vertex27 has precisely the neighbors9 and22: 27+9=36 and27+22=49. Its arrival degree2 is not special: every arrival19..32 has degree2. In G26, vertex18 is a leaf, hence one endpoint of every Hamiltonian path. There are12 old Hamiltonian paths, with other-endpoint multiset

    2,4,8,9,11,13,13,15,22,22,26,26.

Exactly three extend to27: the one ending at9, and the two ending at22. The remaining32 new Hamiltonian paths place27 internally, necessarily in the segment9–27–22. Since9+22=31 is not square, these are not subdivisions of an existing G26 edge: suppressing27 creates a marked non-square edge9–22. Deleting27 yields32 different spanning two-path covers of G26, in each of which9 and22 are exposed endpoints of different paths. Conversely each cover yields exactly one new unoriented Hamiltonian path. The fibres have size1 because there is only one possible neighbor-pair. This is not a universal fibre-size claim for higher-degree births.

The independent edge-subset census finds1116 total spanning two-path covers of G26, of which exactly32 have this endpoint compatibility. These are near-minimum covers, not minimum covers: pc(G26)=1. Retaining only old Hamiltonian paths would miss32 of35 new paths.

The three extendible G26 paths, each written from its useful port, are:

    9,16,20,5,11,25,24,12,4,21,15,1,8,17,19,6,10,26,23,13,3,22,14,2,7,18
    22,3,13,12,24,25,11,14,2,23,26,10,6,19,17,8,1,15,21,4,5,20,16,9,7,18
    22,14,2,23,26,10,15,21,4,12,13,3,6,19,17,8,1,24,25,11,5,20,16,9,7,18

A genuine internal-join witness is the cover

    (1,8,17,19,6,3,13,12,24,25,11,14,22)
    (9,16,20,5,4,21,15,10,26,23,2,7,18).

Joining its displayed22 and9 through27 gives a Hamiltonian path. All35 paths, all three extensions and all32 compatible covers are retained in results.json.

## 2. Adjacent births expose the response, not just total count

| New N | New neighbors | Old HPs | Endpoint extensions | Internal joins | New HPs | All old two-covers |
|---|---|---:|---:|---:|---:|---:|
|25|11,24|0|0|10|10|192|
|26|10,23|10|2|10|12|668|
|27|9,22|12|3|32|35|1116|
|28|8,21|35|5|47|52|2546|

Every internal deletion fibre in these four rows has multiplicity1. At25 the compatible covers are minimum covers of G24, which has pc=2. At26,27,28 the compatible covers have one more path than the old minimum. This distinction is why the minimum-cover lowering criterion alone is not a path-count recurrence.

The35 count was already known in the September wave-six census and Session52. The new finite contribution here is its independently checked birth/deletion decomposition, not discovery of35 or of the first square-sum Hamiltonian path. The first path is at15; at27 traceability was already present at25,26. Session52's Fortunate comparison fails at27, but that is a failure of a proposed cross-domain bit match, not a graph phase transition.

## 3. All-size path-cover birth polynomial

For any finite simple graph H, let LF(H) be its spanning linear forests (equivalently path covers, identified by selected edges). Write k(F) for the number of components. An exposed vertex has selected degree0 or1; a singleton is one exposed vertex, not two endpoint slots.

Adjoin a new vertex v with arbitrary old-neighbor set A. Define

    e_A(F) = number of exposed vertices of F in A;
    j_A(F) = number of unordered pairs {a,b} of exposed vertices in A
             lying in distinct components of F;
    L_H(t) = sum over F in LF(H) of t^k(F).

Then, as a polynomial identity,

    L_(H+v)(t) = sum_F [t^(k(F)+1) + e_A(F)t^k(F)
                                      + j_A(F)t^(k(F)-1)].

Proof: in any selected linear forest the new vertex has degree0,1,or2. Deletion gives respectively the old cover plus a singleton choice, an exposed endpoint attachment, or a pair of exposed endpoints in distinct old components. Degree2 into the same component creates a cycle and is prohibited. These cases are disjoint and deletion recovers the old forest and attachment choice uniquely. Conversely each permitted choice creates a linear forest. This also proves the identity for an empty H, with L_empty=1. No negative exponent contributes because j_A(F)=0 when k(F)<2.

Writing c_k(H) for the number of k-path covers, E_k(H,A) for the sum of e_A over those covers and J_k(H,A) for the sum of j_A, the coefficient recurrence is

    c_k(H+v) = c_(k-1)(H) + E_k(H,A) + J_(k+1)(H,A).

For nonempty H the Hamiltonian-path specialization is

    h(H+v) = E_1(H,A) + J_2(H,A).

This is an exact response hierarchy, not a closed recurrence in h(H) alone. At27 the two terms are3 and32.

### Exact minimum-cover trichotomy

For nonempty H put m=pc(H). Then

    pc(H+v)=m-1 iff J_m>0;
    pc(H+v)=m   iff J_m=0 and (E_m>0 or J_(m+1)>0);
    pc(H+v)=m+1 otherwise.

This follows by reading the first nonzero coefficient; the isolated-v construction always gives an (m+1)-cover. Session52 already proved the first line and the bound m-1<=pc(H+v)<=m+1. The displayed middle line makes explicit the missing near-minimum-cover channel. In particular, absence of usable endpoints in old minimum covers does not by itself force an increase: a compatible near-minimum cover may join instead. A smallest transparent example is H=1–2–3–4 with A={2,3}. Its unique minimum cover has no usable endpoint, but deleting edge2–3 yields a compatible two-cover; insertingv produces1–2–v–3–4 and preservespc=1.

## 4. A sufficient state for arbitrary future vertex births

A cover F can be compressed to the unordered collection of unordered pairs of its path endpoints. A singleton contributes(a,a). Retain the integer multiplicity of covers realizing each such state. This loses internal order but preserves the exact number of extensions under every specified future vertex-neighborhood sequence.

The transition from a state S on adjoining v with neighborhood A is the sum of:

1. Add(v,v).
2. For each component-pair containing a distinct exposed vertex a in A, replace that endpoint a by v. For(a,a), there is one choice, producing(a,v).
3. For each two distinct component-pairs and each choice of one distinct exposed vertex a,b from them in A, remove the two pairs and insert their remaining endpoints as one pair. For singleton(a,a), the remaining endpoint is a. The two new selected edges are a–v andv–b.

Aggregate equal resulting states with multiplicity. Starting from the empty state of weight1, this transition computes exact endpoint-marked cover counts for every finite graph presented by vertex births. For G_N the neighborhood at birth is {s²-N: N<s²<=2N-1}. The number of pairs is the path-cover count; coefficients of one-pair states count Hamiltonian paths by endpoints.

Why the state suffices: the validity of an attachment depends only on selected degree<=1 and whether two exposed vertices lie in the same selected path. The endpoint pairs record precisely those facts, including the singleton case. Different old covers with the same state admit exactly the same state transitions; summing their multiplicities therefore remains exact. This is a transfer identity, not an efficiency or bounded-state theorem.

Why less information fails: covers {(1,2),(3,4)} and{(1,3),(2,4)} have the same exposed vertex set but differ for A={1,2}. The first has j_A=0, the second j_A=1. Forgetting endpoint pairing admits a cycle or loses a legitimate join. Forgetting singleton identity doubles some transitions. For the cover{(1,2),(3)} and A={1,2,3}, there are two internal joins, so unmarked deletion fibres can have size2 when the new degree is3.

## 5. Provenance, connection contract, and boundaries

Closest mechanism: Session32 square-sum README, section3, the exact marked two-path-cover correspondence. Session52 theorem sections4–5 adds the componentwise-traceability distinction and minimum-cover birth law. The canonical hostile remains a graph-degree2 vertex that is an endpoint: degree2 does not force both edges. Useful formerly omitted sidecar here: the near-minimum-cover response and pairing of exposed endpoints.

The square-sum growth map sends an old endpoint-marked cover plus an allowed attachment choice to a new cover; it preserves acyclicity, spanning, and selected degree<=2. Deletion is its exact inverse with the attachment marked. Discarding the marks loses multiplicity; retaining component endpoint pairing repairs that loss.

The older summand-closure note distinguishes three founders{1,2,4} from closure holes{1,4,6}; its dyadic27 coefficient belongs to the distinct summand closure, not G27. Its additive comparison with Collatz also has a concrete obstruction: the shortcut target(3n+1)/2 for odd n cannot be a square, because 2k²=1 mod3 is impossible. Thus the natural shortcut summand hyperedge does not belong to the square-target filter. Berggren transport concerns odd coprime root pairs with multiplier/exponent guards; these are not the arbitrary square-sum graph vertices or its endpoint-cover states. No map preserving the relevant dynamics is supplied by equality of the label27, square vocabulary, or a matching finite bit word.

These identities are elementary graph enumeration, with no literature-priority claim. No claim is made about statistical significance of27, Collatz convergence, LRC, a scalar closed recurrence, or an efficient all-N Hamilton-path algorithm.

## 6. Reproduction and independent tests

Run from this directory:

    python3 check.py > results.json
    python3 -O check.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 test_rejections.py > test.out
    python3 -O test_rejections.py > test.optimized.out
    cmp test.out test.optimized.out

The positive census uses two different search spaces: vertex-order Hamiltonian-path DFS, and selected-edge acyclic degree<=2 spanning-forest DFS with rollback union-find. Equality compares the full compatible-cover sets, not just their counts. The general transfer, polynomial identity and minimum-cover trichotomy are independently checked on every labeled simple graph of order1 through5. Hostiles exercise endpoint pairing, same-component cycles, singleton overcount, higher-degree deletion multiplicity, invalid vertices, repeated vertices, nonedges, and reused birth labels. Explicit checks remain active under Python -O. The finite tests support the elementary proof; they are not its all-size justification.
