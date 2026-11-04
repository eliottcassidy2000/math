# A physical cut port detects every permutation of arbitrary ordered blocks

**PROVED:** A nonidentity permutation of an ordered join of arbitrary nonempty tournament blocks is visible to the full named two-port response in every fixed ambient quotient, even when each block is replaced by another block with the same complete ordered-path-cover profile. The blocks need not be transitive or isomorphic. Consequently every nonempty laminar simultaneous original-module family having a not-strongly-connected maximal selected module is visible. The still-unresolved fully covered laminar case can therefore be restricted to strongly connected maximal selected modules.

Here “strongly connected” always refers to directed reachability, not to the distinct modular-decomposition term “strong module” (a module that overlaps no other module). Names and ports remain physical throughout. This is an exact integer separator; it does not discard cubic or quartic terms or claim a universal reconstruction theorem.

## 1. Profiles and a slotwise injection

For a nonempty tournament A, let q_h(A) count its ordered spanning covers by h nonempty directed paths. An ordered cover has named slots 1,...,h. Two elementary profile facts are needed.

**Substitution profile lemma.** In a fixed quotient Q with distinguished root z, for every h there are nonnegative integer coefficients c_(h,r), independent of A, such that

    q_h(Q[A]) = sum_r c_(h,r) q_r(A).

Compress all maximal A-runs in each of the h path slots. Keep the ordered slots and all exterior names. Reading the slots in order supplies a canonical ordering of the r A-runs. Conversely, an exterior skeleton and an ordered r-cover of A reconstruct the original cover. This gives the claimed disjoint counting decomposition. For h=1 its one-run coefficient is c_(1,1)=H(Q)>0.

In particular, substituting blocks with equal complete q-profiles preserves every ambient q-profile. Several blocks may be replaced one at a time. A module converse preserves each q_h exactly: reverse every maximal run of that module in every path slot. This includes the full-module converse and all legal sequences of current module reversals. The invariance is hereditary: on any fixed physical subset S, the intersection of a currently valid module with S remains a module of the induced tournament. Reverse that intersection (with sets of size zero or one doing nothing), and apply the same run-reversal argument. Thus every such legal sequence preserves H and every ordered-cover profile on every induced subset S.

**Ordered-quotient domination.** For fixed nonempty blocks A_1,...,A_k, any quotient K satisfies

    q_h(K[A_1,...,A_k]) >= q_h(A_1 -> ... -> A_k)

for every h. In each source slot, a block occurs at most once, as a contiguous segment. Reorder those segments in a fixed chosen Hamilton order of K restricted to the occupied block indices. Such an order exists by the elementary tournament insertion argument. Internal segment orders and the physical slot remain unchanged. This is injective: the occupied names identify the segments, which can be sorted back into the original block order. Distinct slots may use different subsets of a block; no independence assumption on endpoint marginals is made.

Also

    H(K[A_1,...,A_k]) >= H(K) product_i H(A_i).

Concatenate one Hamilton path from each block along each Hamilton ordering of K. The physical blocks identify that ordering, so these counted paths are distinct. A nontransitive tournament K has H(K)>=3: start with a directed triangle's three Hamilton paths and insert the remaining vertices by any fixed legal insertion rule; deletion of each inserted name makes the maps injective.

## 2. The named cut-port theorem

Let M be partitioned into fixed nonempty physical blocks A_1,...,A_k, with k>=2. In T[M] their quotient order is

    A_1 -> A_2 -> ... -> A_k.

In U[M] the same blocks occur in a different linear order pi, with possibly different internal tournaments B_i on the same names, satisfying

    q_h(A_i)=q_h(B_i) for every i,h.

Outside M, T and U have the same rooted ambient quotient R and the same incidences to M.

Choose t in {1,...,k-1} such that the original prefix {1,...,t} is not a prefix set of pi. Such a t exists: if all prefix sets agree, the two linear orders agree. Then some j>t occurs before some i<=t in pi.

Use these identical physical port attachments on both sides:

* x->y, and x dominates every old name of M;
* every A_i with i<=t dominates y, and y dominates every A_j with j>t;
* both x and y have M's incidences to every old exterior name.

Thus M together with x,y is a common module in both actual extensions. The source extended quotient is transitive, with x first and y inserted at cut t. In the target, the three blocks y,A_j,A_i form a directed triangle:

    y -> A_j -> A_i -> y.

First replace B_i by A_i in the target extension using the substitution profile lemma. This changes no q_h. The ordered-quotient injection then gives coordinatewise domination of the source extended module's complete q-profile. Moreover its q_1 strictly increases:

    H(target extended module) - H(source extended module)
        >= 2 product_i H(A_i) > 0.

The ambient substitution coefficients are nonnegative and their one-run coefficient is H(R)>0. Therefore the actual, identical physical query satisfies

    F_U(a,c)-F_T(a,c) >= 2 H(R) product_i H(A_i) > 0.

This proves visibility in every ambient quotient, including strongly connected ones, with no bound on the number of module runs. It uses complete profile domination, not an invalid exact scalar multiplication rule for a general ambient contrast. The lower bound arises from the one-run contribution only.

The result is stronger than the earlier transitive-source-module rank theorem in this respect: each A_i here may be an arbitrary nontransitive tournament. It does not assert that a port cut works for a prime quotient converse.

## 3. Consequence for fully covered laminar families

Let F be a nonempty laminar family of distinct nontrivial modules of the original tournament T, with simultaneous edge-XOR semantics. Suppose a maximal selected module B is not-strongly-connected. No extra scalar-equality restriction is needed: a laminar family is always a legal sequence of current module reversals. Reversing a containing set reverses a contained module's uniform boundary, reversing a contained set changes no boundary, and reversing a disjoint set changes no boundary. All selected modules therefore remain modules throughout. By Section 1, this already preserves H and every q_h on every physical induced subset.

If B has an exposed name, the separate `simultaneous_response/ANCHORED.md` theorem already supplies a separator. Otherwise its maximal proper selected submodules A_1,...,A_k are pairwise disjoint and cover B, with k>=2. They form a modular partition. Let K be their quotient.

K is not-strongly-connected. Indeed, substituting arbitrary nonempty tournaments into a strongly connected quotient of at least three vertices gives a strong tournament: paths between different blocks follow quotient paths, and a directed quotient cycle through a block connects any two names in that same block via other blocks. Thus strongly connected K would imply strongly connected B.

The strongly connected components of K form a linearly ordered quotient. Group the A_i belonging to each component into R_1,...,R_s, with s>=2 and

    R_1 -> R_2 -> ... -> R_s.

Every R_j is a module of B and hence of T, since each quotient SCC is a module of K. Each selected set is contained in one R_j, contains all R_j (the selected B), or is disjoint from B. Therefore every R_j remains a module after every reversal, by the same uniform-boundary argument. Reversing selected B reverses this component order. All other selected sets inside B are contained in one A_i and hence in one R_j, so they do not alter the reversed inter-R_j arcs. Within each R_j the induced operation is a whole converse, followed by the selected reversals contained there. These are a legal sequence of module reversals: laminarity preserves all relevant module conditions. Therefore the complete q-profile of each R_j is preserved.

Apply the cut-port theorem to this reversed ordered join. Since no selected set properly contains maximal B, no reversal changes any incidence from B to its exterior. Both x and y keep B's original exterior incidences on both sides, as well as the fixed x->y arc. The extended B is consequently a common module in the same ambient quotient. Selected modules outside B are disjoint from B, since B is maximal in a laminar family. Each lifted port sees any such module uniformly. Their reversals preserve every actual query count exactly and therefore cannot cancel the strict contrast from B.

Hence this family is separated. In particular, a hypothetical full-response collision from a fully covered laminar family must have every maximal selected module strong. No claim here settles the remaining strong-quotient converses or arbitrary legal inversion sequences.

## 4. Verification boundary

`verify_ordered.py` independently computes complete q-profiles by subset Hamilton-path DP and ordered subset-cover convolution. It checks all specified block orders, actual cut queries, profile domination, the quantitative bound, and fixed strongly connected ambient lifts. It also exhaustively checks the covered-laminar not-strongly-connected-maximal criterion through five old names against the retained full response banks. These finite controls support the implementation; the all-order claims are the explicit injections and module arguments above.

No finite reconstruction census is used in this theorem. The fully covered case with strongly connected maximal selected modules remains open here.
