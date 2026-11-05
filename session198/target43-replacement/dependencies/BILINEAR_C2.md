# The full C2 scalar class as two bilinear images of physical start-one sources

Session197. This packet replaces an arbitrary post-pair gadget word by two actual rooted sources. It gives an exact all-order characterization and a finite four-coordinate inverse algorithm for the **full strong C2 class**. It does not prove that this class, together with C1, covers every odd integer outside7,21,67.

Here H counts directed Hamilton paths and C counts directed Hamilton cycles modulo rotation. Empty paths are permitted only in the explicitly optional cover slot, and H(empty)=1.

## 1. The physical source inventory

Let S be the set of quadruples

    (D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c))             (1)

of actual nonempty rooted tournaments(V,c) with exactly one Hamilton path starting at c. Here R_c counts ordered covers(P,Q) with Q nonempty and starting at c, and P possibly empty. The singleton has quadruple

    o=(1,1,1,1).

No free choice of d or formal inverse-closure assumption is made. Let S* exclude this singleton quadruple. It excludes precisely the singleton graph: a fixed root-start Hamilton path gives |V| different reserved covers by putting its initial segment in Q and its remaining suffix in P, so U>=|V|.

The reviewed ROOTED grammar and the SCC decomposition give an exact physical generator of S. Start at o and allow:

    O(D,E,U,d)=(D,0,2U,d),                          (2)

which appends a universal sink while retaining c, and

    M_L(D,E,U,d)=(g(L)D-(L mod2)E+LU,
                  D-[L=1]E, Q(L)D+B(L)U, D),       (3)

for every L>=1, which attaches the ordinary rooted one-apex gadget and uses its new root. Its deletion is transitive, its ending count is one, and its root-ending two-cover count is L. The exact g,Q,B recurrences are inherited from ROOTED.md and implemented in core.py.

Every move preserves root-start count one. Conversely, every strong start-one source is supplied by the complete ordinary grammar. Any weak start-one source consists of a strong start-one initial SCC, or a singleton, followed by singleton SCCs: its starting-path count is the initial root-start count times the H counts of later SCCs, and a nontrivial strong SCC has H>=3. Thus its weak tail is a sequence of O moves. This proves that(2)-(3) generate exactly S.

For a nontrivial source, inserting c immediately before its first out-neighbor in any path of V-c gives a non-ending path. Such an out-neighbor exists because the unique c-start path visits another vertex. Hence also `D-E>=d` for every source in S*.

All source states satisfy

    1<=d<=D<=U,       0<=E<=D.

The inequality D<=U follows by cutting a Hamilton path immediately before c; d<=D follows by deterministic insertion of c into a path of V-c. Both moves strictly increase U. In(3), Q,B>=1 and gD-(L mod2)E>=0, so U_new>U and D_new>=LU.

Therefore S_B={s in S: U<=B} has an exact finite enumeration: at each visited state try O and all `1<=L<=floor(B/U)`, retaining only successors with U_new<=B. Every route to an accepted state stays below B, and every possible last L is tested. At most B^4 quadruples are visited. This is a forward physical enumeration with recorded predecessor recipes, not a numerical cone asserted to be reachable.

## 2. A general two-entry gluing identity

For any tournament T and distinct a,b, put

    h=H(T), E=e_a(T)+e_b(T), R=R_a(T)+R_b(T),
    K=K_ab(T)+K_ba(T),                              (4)

where K_ab counts ordered covers by two nonempty paths, the first ending at a and the second starting at b.

For any nonempty rooted tournament(G,r), orient every cross edge T->G except `r->a,b`. Put

    g=H(G), s=s_r(G), L=L_r(G), d=H(G-r),

where L_r counts covers(X,Y) with X nonempty and ending at r and Y possibly empty. The resulting actual graph F has

    H(F)=h g-E s+R L-K d.                            (5)

**Proof.** A path using no reversed seam is a T path followed by a G path; the forbidden joins subtract Es from hg. A path using a reversed seam has T pieces P,Q, with Q starting at one selected entry, and G pieces X,Y, with X ending at r. The optional P and Y may be empty. These choices contribute RL. Its only remaining invalid join is P->X when X is the singleton r and P ends at the other selected entry. The remaining G-r vertices form one directed Hamilton path in d ways, giving the correction dK. No path can use two reversed seams because both leave r. This proves(5) for arbitrary G, not just a transitive deletion.

If |G|>=2, the cycle count is

    C(F)=(s_a(T)+s_b(T))*e_r(G).                    (6)

Every cycle uses one G-to-T seam, so its T and G portions are single paths. The G path ends at r and starts elsewhere, making its incoming cross join automatically legal. Conversely every such pair of paths closes uniquely. When G is a singleton this argument fails: the preceding T endpoint can be the other selected entry. That case is deliberately excluded from the main bilinear theorem.

For later bounds, two elementary inequalities on(4) are useful:

    E<=h,        0<=K<=R,        R-K>=h.             (7)

The endpoint counts are disjoint, and the K covers are a subset of the R covers. To inject Hamilton paths into the remaining covers, cut before a unless its predecessor is b; in that case cut before b. The second cut cannot have predecessor a. Thus the chosen cover is not counted by K, and concatenating its two pieces recovers the original word.

If |G|>=2 and e_r(G)=1, then L>=d+1: the d covers with X=singleton r are disjoint from the one cover with Y empty. Since s<=g, (5) and(7) imply

    H(F)>=h L,             H(F)>=R.                 (8)

Indeed `hg-Es>=0`, and `LR-dK=L(R-K)+(L-d)K`, with L-d>=1. These bounds are valid without cancelling signed terms numerically.

## 3. Every marked C1 pair has one of two actual source forms

Use J_c(V) for the birth of a vertex whose only old out-neighbor is c.

### One-step pair N

Take (V,c) with quadruple v=(D,E,U,d) in S*. Form T=J_c(V), call the newborn a, and select the pair{a,c}. Its response(4) is

    A_N(v)=(D-E+U, D-E+d, D+3U-d, U+E).             (9)

T has exactly one cycle, and both selected vertices have one starting path. The nontriviality of V is required here. The adjacent pair response follows from the exact J and retained-root identities:

    e_a=D-E, R_a=D+U, e_c=d, R_c=2U-d,
    K_ac+K_ca=R_a-e_a=U+E.

The latter is the reviewed adjacent-pair cover split, also checked directly in the verifier.

### Two-step pair M

Take v in S, including the singleton. Form S=J_c(V), call its newborn b, then T=J_b(S), with newborn u. Select{u,c}. Its response is

    A_M(v)=(D+2U, U+d, 2D-E+6U-d, U+E).             (10)

This graph is also obtained by appending a universal sink u to V and then appending b beating only c. The retained-root O/Z laws give

    e_c=d, R_c=4U-d,

while two J calculations give

    e_u=U, R_u=2D-E+2U.

For K_uc, b cannot lie in the path ending at u: its only out-neighbor c belongs to the other path, and b does not beat u. Hence b is last in the c-start path, giving U minus the d covers whose old reserved path is singleton c. Thus K_uc=U-d. For K_cu, the u-start path is either singleton u, leaving d possibilities, or exactly u,b, leaving E possibilities. Thus K_cu=d+E. This proves(10).

T is strong C1. Its u-start path is unique, and its c-start path is the unique c-start path of V followed by u,b. For singleton V, T is the triangle; this boundary is a two-step pair, not a genuinely nonadjacent pair.

These two forms exhaust every strong C1 graph with an unordered pair of distinct start-one vertices. The proved START_TWO theorem says that these vertices form a cyclic-consecutive block of two or three, with all but its last vertex having outdegree one; the triangle is handled separately. For an adjacent selected pair a->c, a has sole out-neighbor c. Deleting a gives exactly the physical source V in(9), and its unique c-start path is obtained by cutting the old unique cycle. For the nonadjacent first/last pair u,c, the intervening b has `u->b->c`, with u,b each of outdegree one. Deleting u,b gives exactly V in(10). Cutting the unique cycle again proves s_c(V)=1. Every edge is fixed by the degree-one constraints. No arbitrary graph rebasing is inserted into this reduction.

## 4. Two bilinear forms characterize the full C2 scalar class

For a second actual source w=(g,s,L,d) in S*, let G be its converse, marked at the same vertex r. Then e_r(G)=1 and |G|>=2. Define

    N(v,w)=(D-E+U)g-(D-E+d_v)s
            +(D+3U-d_v)L-(U+E)d,
                 v in S*, w in S*,                 (11)

    M(v,w)=(D+2U)g-(U+d_v)s
            +(2D-E+6U-d_v)L-(U+E)d,
                 v in S*, w in S*.                 (12)

Here d_v is the fourth coordinate of v; d is the fourth coordinate of w. Both are actual deletion counts. These forms are exactly(5) applied to(9)-(10).

**Theorem.** The set of Hamilton-path counts of all strong tournaments with exactly two Hamilton cycles is

    {15} union N(S*,S*) union M(S*,S*).               (13)

Every value in the two bilinear images has an explicit physical C2 graph, by(6). The exceptional value15 is realized by the regular five-tournament.

For the converse, use the completed C2 classification: every C2 graph except that regular five-tournament has an outdegree-one vertex v. Let a be its sole out-neighbor. Deleting v leaves a rooted graph with exactly two paths starting at a, by cutting the two cycles at v. Its initial SCC U is strong with start count two; all later SCCs are singletons, because their H product divides two and every nontrivial strong component has H>=3.

Apply the canonical START_TWO decomposition to (U,a), following its one-entry case backwards until its first two-entry case. This leaves a strong C1 graph T with a selected start-one pair. All vertices subsequently added, the tail SCCs, and the final closing vertex v form one tournament G. There is exactly one vertex r of G with reverse edges into T, namely the first pair-entry root, and its reverse edges go precisely to the selected pair. Every other cross edge is T->G. Thus the original graph is literally a graft of Section2.

G contains both r and v, so it is nontrivial. The cycle identity(6) and C=2 force e_r(G)=1. Its converse is therefore an actual source w in S*. This argument avoids any assumption that arbitrary endpoint-one graphs have a particular chosen suffix decomposition. Finally Section3 places the marked T pair in one of the two source forms, giving(11) or(12). This proves(13), including actual graph reconstruction, not just equality of scalar recurrences. The triangle pair can always be put in the one-step N form, so both sources may be required to be nontrivial. If the singleton left M boundary is retained for convenience, it contributes nothing new: `A_M(o)=A_N(1,0,2,1)=(3,2,6,2)`.

`factor_graph` implements this proof on a promised C2 graph: remove an outdegree-one vertex, take SCCs, follow the unique-entry reductions until the pair, split the two sources, and reconstruct every edge under an explicit map. It recognizes the regular-five exception separately.

### Ordered graph presentations are bijective

The preceding proof is a surjection from ordered pairs of nontrivial rooted start-one graph isomorphism classes, with a choice of N or M, onto all nonregular C2 graph classes. It is in fact a bijection at every order.

The proved rooted-source generating function, including the singleton and all weak tails, is

    W(z)=z(1-2z)/(1-3z+z^2).

Indeed FIBONACCI_CLASS gives the strong/singleton rooted count `B(z)=z+z^3/(1-3z+z^2)`, and the unique weak-tail decomposition gives `W=B/(1-z)`. N adds one vertex to its two sources and M adds two. Therefore the ordered presentation generating function is

    (z+z^2)(W(z)-z)^2
      = z^5(1-z)(1-z^2)/(1-3z+z^2)^2,

exactly the independently proved nonregular C2 generating function in C2_COMPLETE. Every coefficient counts finite sets. A surjection between finite sets of equal cardinality is injective, proving the all-order graph normal form. The construction cannot produce R5: its only order-five presentation has two transitive two-vertex sources and gives H13, whereas R5 has H15.

Ordered sources are essential for this graph statement. The scalar expressions happen to be symmetric, so the scalar inverse below may use unordered numerical source pairs. It does not assert that equal source quadruples identify rooted graphs.

## 5. Positive symmetric forms and an exact bounded inverse

For a nontrivial physical source define

    x=(D-E-d, U-D, E, d).

These four coordinates are nonnegative, the last is at least one, and their sum is U. This change of coordinates is invertible: `D=x0+x2+x3`, `E=x2`, `U=sum(x)`, `d=x3`. It does not turn the entire nonnegative cone into a physical source inventory.

Direct expansion of(11)-(12) gives `N(v,w)=x(v)^T A_N x(w)` and likewise for M, with the symmetric positive matrices

    A_N = [ [6,4,5,5], [4,3,4,3],
            [5,4,5,3], [5,3,3,4] ],

    A_M = [ [11,8,10,10], [8,6,7,7],
            [10,7,9,8],  [10,7,8,9] ].             (14)

In particular both scalar forms are symmetric in their two physical sources. Subtracting three times the all-ones matrix from A_N leaves a nonnegative matrix with last diagonal entry1. Subtracting six times the all-ones matrix, then the last row and last column of ones, then the last diagonal unit, from A_M also leaves a nonnegative matrix. Since each last coordinate d is at least one, this proves

    N(v,w)>=3U1 U2+1,
    M(v,w)>=6U1 U2+U1+U2+1.                        (15)

Every nontrivial source has U>=2. Thus if either form equals target n, both source cover counts are at most `floor((n-1)/6)`. By scalar symmetry the sources may be ordered so U1<=U2; the smaller one then satisfies

    U1<=floor(sqrt((n-1)/3))                       (16)

for N, with the stronger M product bound available separately.

For comparison, the direct graph inequalities(8) remain useful: `L<=n/h` and `R_pair<=n`. In the N left source `R_pair>=3U_left`; in the M left source `R_pair>=6U_left`. The matrix bounds strengthen the common source cutoff without any signed numerical pruning.

The exact inverse is now finite. Enumerate S_B by Section1 with `B=floor((n-1)/6)`, remove the singleton, and test the two bilinear forms on unordered pairs, retaining only products allowed by(15). Also allow the exceptional15. Every omitted source or pair is excluded by a proved bound, and every retained state has a physical predecessor recipe. The test is therefore complete for the full C2 scalar class and constructs an actual graph for every accepted target.

The crude inventory bound is B^4. Its generation and at most twice B^8 bilinear comparisons are polynomial in the **numerical** target under integer arithmetic; intermediate values have O(log n) bits. This is not polynomial in input bit length. Predecessor pointers suffice to reconstruct sources. No graph-order cutoff, Busch estimate, numerical saturation or unrestricted formal inverse is needed.

This gives a precise target-independent arithmetic reformulation of the remaining low-cycle coverage problem: after adding the complete C1 scalar image, does the union in(13) contain every odd number outside7,21,67? The theorem does not answer that representability question. It removes arbitrary post-pair word depth, not physical source membership.

## 6. Scope guards and finite validation

The previous generalized N/M families are recovered by restricting w to the converse of a single one-apex gadget. Here w is any actual start-one source, including weak ones, with its genuine deletion sidecar. This is why the two-source reduction has a full-class completeness proof rather than only a subclass construction.

The graft entry is not a returning start-two interface in general. Its start count is `2d`, because a path beginning there first enters one of the two left roots, traverses the left graph, then traverses G-r in d ways. For example the actual graph

    [58,124,57,16,96,8,45], marked root6,

has H41 and C2 but six paths starting at the marked graft entry. Its left source is `(1,0,2,1)` and right source `(5,2,6,3)`. The theorem constructs the graph and can factor it again; it does not silently retain the old three-coordinate start-two state.

The verifier checks the general arbitrary-G scalar formula, both left pair response formulas, actual source recipes, the complete bounded image and physical reverse factorization. The finite benchmark at target bound1861 has2,048 source quadruples,10,678 tested bilinear pairs and884 distinct C2 values including15. Every stored value has an actual graph counted directly with C=2. These bounded results do not establish cofinite coverage.

The structural dependencies are the independently reviewed ROOTED.md, START_TWO.md, C2_COMPLETE.md and FIBONACCI_CLASS.md snapshots copied under dependencies/. ROOTED supplies the ordinary source transfer and complete strong start-one grammar; START_TWO supplies the unique-start block and canonical start-two decomposition; C2_COMPLETE supplies the regular-five/degree-one dichotomy and the unrooted graph count; FIBONACCI_CLASS supplies the rooted strong-source count. All new gluing, pair response, positive bound and bilinear-reduction arguments are given above. The unchanged core.py is the reviewed graph/counter helper, not a new external dependency.
