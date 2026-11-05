# Every strong tournament with three Hamiltonian cycles has a strong bridge

Session197. **All-order structural theorem, with an exact finite low-cycle census through order12.** No literature-priority claim is made. The proof uses the independently audited C1 and C2 classifications in `C1_COMPLETE.md` and `C2_COMPLETE.md`, the one-seam endpoint lemma in `ROOTED.md`, and the separating-vertex product in `C2_CUT_REDUCTION.md`. It does not use a published quantitative order bound or infer stabilization from a census.

We use standard Hamiltonian-cycle existence for strong tournaments (Camion). Let C(T) count directed Hamiltonian cycles modulo rotation, and H(T) count directed Hamilton paths as ordered vertex words. A strong bridge is an arc whose deletion destroys strong connectivity. The conclusion is

    C(T)=3 and T strong  =>  T has a strong bridge.

Equivalently within this class, T has an outdegree-one vertex or an indegree-one vertex. The two orientations need not both occur. No closed all-order class-count generating function is asserted here.

## 1. Strong induced extension and an exact two-vertex ear

**Lemma1.** If S is a strong induced subtournament of a strong tournament T, then C(S)<=C(T). Furthermore, when |S|>=4 and C(T)<=3, there is a chain from S to T adding one vertex at a time, with every intermediate induced tournament strong and having at most three cycles.

**Proof.** If an outside vertex v has both an in-neighbor and an out-neighbor in S, insert it into a cyclic in-to-out gap of every Hamilton cycle of S. Choosing one gap deterministically gives an injection: deleting the same named v recovers the old cycle. The enlarged induced graph is strong.

Otherwise every outside vertex uniformly dominates S or is uniformly dominated by S. Both types must occur. Let P dominate S and let S dominate Q. Some q in Q beats some p in P: if not, P would dominate both S and Q, contradicting strongness of T. The induced graph on S,q,p has

    S -> q -> p -> S.

Every Hamilton cycle of this graph consists of one Hamilton path of S followed by q,p. Hence its exact cycle count is H(S). Cutting each old Hamilton cycle at each of its |S| edges gives distinct Hamilton paths, so

    C(S+q+p)=H(S)>=|S| C(S).                         (1)

This is another strong induced extension with nondecreasing cycle count. Repeating the one- or two-vertex construction reaches T and proves monotonicity.

If |S|>=4 and C(T)<=3, no two-vertex step can occur in any such growth: equation(1) would give at least four cycles, and subsequent steps cannot decrease the count. Thus a mixed-incidence outside vertex always exists until T is reached. Every intermediate count is at most C(T). QED.

**Lemma2.** Every strong tournament of order at least four contains a strong induced four-vertex tournament.

**Proof.** A strong tournament contains a directed triangle, by taking a shortest directed cycle and using its chords. Apply the preceding extension choice to that triangle. A mixed-incidence vertex immediately gives a strong four-vertex graph. Otherwise use q,p as above. Choose any two triangle vertices a,b with a->b. The induced four vertices have the Hamilton cycle a,b,q,p,a and are strong. QED.

There is exactly one strong four-vertex isomorphism class, with rows

    [2,4,9,3].

Indeed all degrees lie between one and two, so there are two vertices of each outdegree. If the degree-one vertices are a->b, then b beats a third vertex c. The remaining vertex d beats a,b, and strongness forces c->d. All arcs are determined. The displayed graph has C=1.

**Consequence.** Every strong tournament with C<=3 and order at least four is generated from this one seed by strong singleton extensions, retaining only outputs with C<=3. This is the all-order completeness justification for the finite census below.

## 2. Exact finite census with complete isomorphism handling

`verify_c3_complete.cpp` performs that generation through order12. For every retained representative A, it considers every nonempty proper subset O of old vertices as the outneighbors of a new vertex. Because A is strong and both incidence signs occur, the extension is strong. No strong singleton extension is omitted.

For an old Hamilton cycle, every cyclic transition from an old vertex outside O to one inside O permits an insertion of the newborn. Different gaps and different old cycles give distinct new cycles. The sum of these insertion counts is therefore a valid lower bound. The verifier discards a mask only when this lower bound exceeds three, then counts the remaining extensions exactly with a Hamilton-cycle recursion capped at four. Capping is sound because only the distinction between1,2,3 and at least4 is needed.

For each retained graph it lists all its Hamilton cycles with vertex0 fixed first. To canonicalize, it places each of these cycles in turn in the natural order, with every cyclic rotation allowed, and chooses the lexicographically smallest complete adjacency encoding. This is an exact isomorphism canonicalization: every isomorphism maps the set of Hamilton cycles and their rotations bijectively, while equality of full adjacency encodings supplies an isomorphism. No degree-only heuristic or arbitrary tie-break is used.

The complete unlabelled inventory is:

    order:  4   5   6   7    8    9    10    11     12
    C=1:    1   3   8  21   55  144   377   987   2584
    C=2:    0   2   5  18   60  191   589  1775   5257
    C=3:    0   1  11  38  124  390  1193  3574  10537

Every retained C3 graph has an indegree-one or outdegree-one vertex. Such a vertex's sole incoming or outgoing arc is a strong bridge. Thus no bridge-free C3 graph exists through order12. At orders below four there is none either: a singleton or a two-vertex tournament has no Hamilton cycle, and a triangle has one.

The C1 and C2 rows agree at every order with their independently proved all-order enumerations. The certificate additionally records per-order and aggregate SHA256 hashes of every canonical adjacency record, exact extension counters, and explicit tournament/degree validation. The finite C3 row alone is not claimed as an all-order enumeration.

## 3. A bridge at cycle count three forces one degree-one orientation

**Lemma3.** A strong C3 tournament with a strong bridge has an indegree-one or outdegree-one vertex.

**Proof.** Split a strong bridge b->a into nonempty sides A,G with every cross arc A->G except b->a. The exact cycle law gives

    C(T)=s_a(A)e_b(G)=3.

One factor is one. If e_b(G)=1, the endpoint-one lemma gives an outdegree-one vertex in T. To include weak G explicitly, all SCCs before its terminal component containing b are singletons. If that terminal component is the singleton b, then b has only the outgoing arc b->a. Otherwise its strong endpoint-one structure has an outdegree-one vertex different from b; it loses to every earlier G vertex and every A vertex, so the degree remains one globally. If s_a(A)=1, reverse all arcs and apply the same argument to obtain an indegree-one vertex. QED.

The converse holds without the C3 restriction: a degree-one vertex's sole incoming or outgoing arc is a strong bridge. Thus, for strong C3 graphs, bridge existence is equivalent to the disjunction of the two degree-one conditions.

## 4. The smallest hypothetical bridge-free atom has order at most12

**Lemma4.** If any bridge-free strong C3 tournament exists, a smallest one has order at most12.

**Proof.** Let T be a smallest counterexample. If its order is already at most12 there is nothing to prove, so suppose n>=13. The separating-vertex cycle product in `C2_CUT_REDUCTION.md` shows that a bridge-free strong tournament with a separating vertex has at least four cycles. Hence T is 2-strong.

For every vertex v, T-v is strong. Single-vertex cycle insertion implies1<=C(T-v)<=3. If its count is one, the complete C1 theorem gives a degree-one vertex. If its count is two, the complete C2 theorem gives one unless T-v is the regular five-tournament; its order n-1>=12 excludes that exception. If its count is three, minimality gives a bridge, and Lemma3 gives a degree-one vertex of at least one orientation.

Let L be the set of outdegree-two vertices of T and R the set of indegree-two vertices. Because T is 2-strong, its minimum indegree and outdegree are at least two. An outdegree-one vertex of T-v must therefore be some u in L with u->v; an indegree-one vertex of T-v must be some w in R with v->w. Thus every vertex v belongs to

    N^+(L) union N^-(R).                              (2)

Write l=|L|. The 2l outgoing arcs from L contain binom(l,2) internal arcs. Therefore l<=5, and

    |L union N^+(L)| <= l+2l-binom(l,2)
                       = l(7-l)/2 <=6.

This includes l=0, when the set is empty. Reversing arcs gives the same bound for R union N^-(R). Equation(2) now implies n<=6+6=12, a contradiction. QED.

**Theorem.** Every strong tournament with exactly three Hamiltonian cycles has a strong bridge and hence an indegree-one or outdegree-one vertex.

**Proof.** If a counterexample existed, choose a smallest one. Lemma4 puts it at order at most12, where the complete census in Section2 excludes it. QED.

This is a finite-reduction proof with a justified all-order search cutoff, not a claim that an observed pattern stabilizes.

## 5. Scope controls and reproduction

Unlike the C2 theorem, the two degree-one orientations are not equivalent. The six-vertex tournament

    [14,52,56,50,33,1]

has exactly three Hamilton cycles and one outdegree-one vertex, but no indegree-one vertex. It is obtained from r->C3->v->r by adding b that beats only r. Its converse gives the reverse asymmetry. Thus neither the start-two indegree lemma nor the C2 rigidity/enumeration conclusions are imported here.

Run `python3 verify_c3_complete.py` and `python3 -O verify_c3_complete.py`. The wrapper uses only the Python standard library and an installed C++17 compiler, compiles the adjacent source into a temporary directory, and rebuilds the entire finite universe. It writes `c3_complete_certificate.json`; normal and optimized Python receipts agree. No precomputed graph bank, external package, repository history, network access or dirty checkout is used.

The structural theorem is all-order. The numerical C3 class counts are proved only through order12. The unrestricted Hamiltonian-path spectrum `{7,21}` remains open.
