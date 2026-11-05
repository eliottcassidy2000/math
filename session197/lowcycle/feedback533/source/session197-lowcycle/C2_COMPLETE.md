# The regular five-tournament is the only bridge-free C=2 tournament

Session 197. **Complete all-order classification and enumeration of strong tournaments with exactly two Hamiltonian cycles.** The proof uses the bridge-class theorem in `BRIDGE_TWO.md`, an elementary degree bound, and an explicit structural argument on at most six vertices. A complete independently audited census is retained only as a verification control. It does not use Busch's quantitative order bound. No literature-priority claim is made.

Let C(T) count directed Hamiltonian cycles modulo rotation. The regular five-tournament R_5 has rows

    [6,12,24,17,3],

equivalently i->i+1,i+2 modulo five.

## 1. Two elementary lemmas

**Lemma 1 (cycle insertion is injective).** If T and T-v are strong, then C(T-v)<=C(T).

**Proof.** In a Hamilton cycle of T-v, vertex v has both an in-neighbor and an out-neighbor, since T is strong. There is an adjacent cyclic pair x,y with x->v->y. Insert v there. Choose one such insertion by any deterministic convention. Distinct old cycles cannot produce the same new cycle, because deleting v recovers the old cyclic word. This gives the injection. This is the counting form of the insertion argument in `C1_COMPLETE.md`. QED.

**Lemma 2 (degree-two cover).** Suppose T is 2-strong and every T-v has an outdegree-one vertex. Then |T|<=6.

**Proof.** A 2-strong tournament has minimum outdegree at least two. If u has outdegree one after v is deleted, then u has outdegree exactly two in T and u->v. Let L be the set of outdegree-two vertices of T, and put l=|L|. Every vertex is beaten by at least one vertex of L.

There are exactly 2l outgoing arcs from vertices of L. Exactly binom(l,2) of them remain inside L, so

    arcs(L,V\L)=2l-binom(l,2)=l(5-l)/2.

Nonnegativity implies l<=5, and l>=1 because L covers all vertices. Every vertex outside L receives at least one of these arcs. Hence

    n-l <= l(5-l)/2,
    n <= l(7-l)/2 <= 6,       for l=1,...,5.

This is an integer bound. It does not count distinct outside vertices more than once, and it does not presume the outgoing arcs have distinct endpoints. QED.

## 2. The elementary small-order argument

A 2-strong tournament has order at least five. At order five, minimum outdegree two and the total of ten arcs force every outdegree to be two. Fix vertex0, label its outneighbors1,2 with1->2, and its inneighbors3,4 with3->4. Degree two forces2->3,4, then1->3 and4->1. Thus the graph is exactly R_5. Its two Hamilton cycles are the step-one and step-two cycles modulo five; direct inspection lists no others.

Suppose now that T is 2-strong, has order six, and C(T)=2. By Lemma1, every vertex deletion has C=1 or2. By the complete C1 theorem and the bridge-C2 theorem, each deletion either has an outdegree-one vertex or is R_5. Here the assertion for order-five C2 is already established: a bridge-free one is 2-strong by the separating-vertex product, hence R_5.

If T-v is R_5, vertex v has two or three outneighbors by minimum semidegree two in T. The step-one and step-two Hamilton cycles of R_5 partition all its arcs. For a two-element subset, the two vertices are consecutive on exactly one of these cycles. The numbers of cyclic transitions entering the subset are therefore one and two. Taking complements leaves these numbers unchanged. Hence the two old cycles admit altogether three different insertions of v, producing three distinct Hamilton cycles of T. This contradicts C(T)=2.

Otherwise every deletion has an outdegree-one vertex. Use the set L from Lemma2. The inequality at n=6 forces l=3 or4.

- If l=4, exactly two arcs leave L for the two outside vertices x,y. Suppose x->y. Minimum indegree two requires at least two incoming arcs from L at x and at least one at y, a contradiction.
- If l=3, every vertex of L receives an arc from L, so L is a directed triangle. Every vertex of L therefore has exactly one outgoing arc to the outside triple. All three outside vertices have indegree at least two, and their total indegree is exactly six: three internal arcs plus three arcs from L. Thus each has indegree two. Each also receives an arc from L by the cover property, so its internal indegree is at most one. The three internal indegrees sum to three, so all equal one: the outside tournament is a directed triangle, and the arcs from L form a perfect matching onto it.

Label L as0->1->2->0 and the matching i->i+3. All other crossing arcs point into L. There are only two orientations of the outside triangle. Both have the alternating Hamilton cycles

    (0,3,1,4,2,5), (0,3,2,5,1,4).

If3->4->5->3, a third cycle is (0,1,2,5,3,4). If the outside triangle is reversed, a third is (0,1,4,3,2,5). Every closing arc is included. Again C(T)>=3, a contradiction.

This proves the order-six exclusion without a census. The complete labelled verifier independently finds 24 two-strong order-five graphs, all R_5, and the order-six histogram C=4:720, C=5:1920. Those finite checks test the elementary proof rather than supply an extrapolation premise.

## 3. The all-order classification

**Theorem 3.** Every strong tournament with C=2 either is R_5 or has a strong bridge.

**Proof.** Suppose otherwise, and choose a counterexample T of minimum order. By `C2_CUT_REDUCTION.md`'s separating-vertex product, a bridge-free C=2 tournament is 2-strong. The elementary small-order argument above implies n>=7.

For every vertex v, T-v is strong. Lemma 1 and Hamiltonian-cycle existence give C(T-v)=1 or 2. If it is one, `C1_COMPLETE.md` gives an outdegree-one vertex in T-v. If it is two, minimality says that T-v either has a strong bridge or is R_5. Its order n-1>=6 excludes R_5. The bridge-C2 equivalence in `BRIDGE_TWO.md` therefore again gives an outdegree-one vertex.

Thus every vertex deletion has an outdegree-one vertex. Lemma 2 forces n<=6, contradicting n>=7. QED.

Conversely, R_5 is 2-strong and has no strong bridge, while every graph in the complete bridge-C2 grammar has C=2. The classification is therefore exact in both directions.

In particular every strong C=2 tournament other than R_5 has an outdegree-one vertex, an indegree-one vertex, a start-two root and an end-two root. The start-two grammar followed by the outdegree-one-root construction in `BRIDGE_TWO.md` supplies all these graphs. The exceptional R_5 has three Hamilton paths from every vertex, so it is genuinely outside that rooted grammar.

## 4. Exact unlabelled and labelled enumeration at every order

Write Delta(z)=1-3z+z^2. The bridge class generating function from `BRIDGE_TWO.md` is

    z^5(1-z)(1-z^2)/Delta(z)^2.

Theorem 3 leaves exactly one extra unlabelled graph, at order five. Thus the complete strong-C=2 generating function is

    sum_n a_n z^n
      = z^5 + z^5(1-z)(1-z^2)/(1-3z+z^2)^2.          (1)

The counts at orders 5 through 11 are

    2,5,18,60,191,589,1775.

Every bridge-C2 graph is rigid. R_5 has automorphism group of order five: cyclic rotations act, and fixing one vertex fixes its two ordered out-neighbors and then the remaining vertices. Consequently the complete labelled count is

    n! * b_(n-1,2) + 24*[n=5],                      (2)

where b_(m,2) is the rooted start-two count from `START_TWO.md`. In particular there are 144 labelled C=2 tournaments at order five, and every C=2 tournament of order at least six is rigid.

The partial polynomial of degree at most73 left undetermined in `BRIDGE_TWO.md` is now proved to be exactly z^5. That older packet remains a valid separately frozen proof; this theorem strengthens it rather than modifying its source.

## 5. Verification and proof boundary

`verify_c2_complete.py` independently exhausts the labelled tournaments through order six, verifies the complete 2-strong cycle histogram, checks the R_5 isomorphism and automorphism statements, checks the insertion inequality on the relevant physical graphs, and verifies the exact degree arithmetic. It also checks the initial enumeration against the complete small universe.

The classification proof is all-order and elementary, including the small-order cases. The finite census, insertion controls and literal cycle witnesses are independent checks of its steps. No stabilization or census extrapolation is used. The classification of the C=1 and C=2 structural classes is now complete within these packets. The unrestricted Hamiltonian-path spectrum `{7,21}` remains open; these cycle-classification results do not claim that every allowed path count belongs to either low-cycle class.
