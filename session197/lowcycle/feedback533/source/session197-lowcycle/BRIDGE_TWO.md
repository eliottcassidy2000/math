# Complete classification and enumeration of the bridge C=2 class

Session 197. **All-order classification of strong C=2 tournaments with a strong bridge, and exact enumeration of all strong C=2 tournaments above order73.** Bridge-free atoms of orders at most73 are not classified here.

This packet uses the completed C=1 theorem in `C1_COMPLETE.md`, the start-two grammar in `START_TWO.md`, and the published residual order bound in `C2_CUT_REDUCTION.md`. No literature-priority claim is made.

## 1. An indegree-one vertex outside any two selected start-one roots

Let T be strong with C(T)=1, and let S be its set of start-one vertices.

**Lemma 1.** If n>=4, some indegree-one vertex of T lies outside S. If n=3, any pair of start-one vertices leaves an indegree-one vertex unselected. Consequently, at every order, a pair of distinct start-one vertices cannot contain all indegree-one vertices.

**Proof.** The converse of the completed C=1 degree-one theorem ensures that T has an indegree-one vertex. If none is in S, the claim follows. Otherwise choose r in S with unique in-neighbor v.

Suppose T-r were strong. Since r beats every vertex except v, the root-start identity would give

    s_r(T)=sum_(a in T-r, a!=v) s_a(T-r) >= n-2 >= 2,

contradicting s_r=1. Thus T-r is not strong. In its ordered SCC decomposition, v must lie in the terminal component: any later component would have no outgoing route to r or to an earlier component, contrary to strongness of T. Hence r beats every vertex of the initial component, and the same root-start identity gives

    s_r(T)=H(T-r)=1.

Therefore T-r is transitive, with v last. The whole graph has the natural order

    r,x_1,...,x_(n-2),v

with all arcs forward except v->r. Vertex x_1 has indegree one. It has at least two Hamilton paths starting there:

    x_1,x_2,...,x_(n-2),v,r;
    x_1,v,r,x_2,...,x_(n-2).

These are distinct and valid for n>=4. Thus x_1 is outside S, proving the claim. The triangle case is immediate because all three vertices have indegree one. QED.

## 2. Every start-two graph has another indegree-one vertex

**Theorem 2.** If T is strong and s_r(T)=2, then T has an indegree-one vertex different from r.

**Proof.** Induct along the complete canonical start-two grammar from `START_TWO.md`. Its initial strong component U has one entry vertex a with s_a(U)=2, or two entry vertices a,b each with start count one. The new root is r and the remaining new vertices form a transitive tail dominated by U.

In the one-entry case, induction supplies an indegree-one vertex u of U different from a. The only new incoming arc into U is r->a, so u keeps indegree one.

In the two-entry case U has C=1. Lemma 1 supplies an indegree-one vertex u outside {a,b}. Again its indegree is unchanged, because the only new incoming arcs into U are r->a and r->b. The triangle input supplies the base case, and every induction step strictly reduces the order. In both cases u is an old vertex, hence differs from r. QED.

A useful independent corollary is that every vertex of a 2-strong tournament starts at least three Hamilton paths. Start count one would give C=1 and an outdegree-one vertex; start count two is excluded by Theorem 2 because an indegree-one vertex prevents 2-strongness. The same holds for ending paths by reversal. Therefore every 2-strong order-n tournament has H>=3n. The regular five-tournament attains equality.

## 3. Equivalent descriptions of the bridge C=2 class

**Theorem 3.** For a strong tournament T with C(T)=2, the following are equivalent:

1. T has a strong bridge.
2. T has an outdegree-one vertex.
3. T has an indegree-one vertex.
4. Some vertex starts exactly two Hamilton paths.
5. Some vertex ends exactly two Hamilton paths.

**Proof.** An outdegree-one vertex starts exactly C(T)=2 Hamilton paths, by cutting and closing cycles. Thus (2) implies (4), and Theorem 2 gives (4) implies (3). Reversal gives (3) implies (5) implies (2). A degree-one vertex exposes a strong bridge through its sole outgoing or incoming arc.

For the remaining implication, split a strong bridge b->a into nonempty sides A,G with all cross arcs A->G except b->a. The exact cycle product is

    C(T)=s_a(A)e_b(G)=2.

One factor is one. The endpoint-one lemma proved in `ROOTED.md` gives an outdegree-one vertex of T when e_b(G)=1, and its dual gives an indegree-one vertex when s_a(A)=1.

For completeness, the first assertion also covers weak G. Its SCCs before the terminal component containing b must be singletons, since e_b(G)=1. If the terminal component is the singleton b, then b's only outgoing arc in T is b->a. Otherwise the strong endpoint-one lemma gives an outdegree-one vertex different from b in that terminal component; it loses to all earlier G vertices and to all A vertices, so its degree stays one in T. The dual argument treats A. This proves (1) implies (2) or (3), which are already equivalent. QED.

Thus the start-two grammar, restricted to its C=2 outputs and then forgetting the root, is a complete recursive classification of the bridge C=2 class. Its C=1 outputs remain separate. No assertion that every C=2 tournament has a bridge is used.

## 4. Rigidity and the unrooted count

A strong tournament of order at least four has at most two outdegree-one vertices, as proved in `FIBONACCI_CLASS.md`. Every bridge C=2 graph has one or two of them.

**Lemma 4.** Every bridge C=2 tournament is rigid.

**Proof.** An automorphism preserves the one- or two-vertex outdegree-one set. In the two-vertex case its mutual arc prevents an exchange, so in either case an outdegree-one vertex v is fixed. The automorphism permutes the two Hamilton paths starting at v. If it fixes one path it fixes every vertex. If it exchanges the two, its square fixes both and is the identity. But a nontrivial involution cannot be a tournament automorphism: any exchanged pair of vertices would have its mutual arc reversed. Thus only the identity is possible. QED.

Write b_(n,2) for the exact rooted start-two counts of `START_TWO.md`. Let c_n count a C=2 graph marked at an outdegree-one vertex. In the canonical start-two grammar, an outdegree-one new root requires the one-entry case and all tail bits zero. Therefore

    c_n=sum_(m<n) b_(m,2).                              (1)

The two-outdegree-one bijection from `FIBONACCI_CLASS.md` preserves **any** fixed positive cycle count, not only one. Add u that beats only a marked outdegree-one vertex v of a strong graph S. The result has exactly two outdegree-one vertices u,v, and

    C(T)=s_u(T)=s_v(S)=C(S).

Conversely, if u,v are the two outdegree-one vertices, orient their names so u->v. Deleting u gives a strong graph S: the path argument in that packet uses only the degree-one constraints and strongness. Vertex v still has outdegree one, and cutting a cycle of T at u gives precisely a Hamilton path of S starting at v. Hence C(S)=C(T). The two operations are inverse on isomorphism classes.

Thus the number of unrooted C=2 graphs with two outdegree-one vertices is c_(n-1). If a_n denotes the number of all unrooted bridge C=2 graphs, rigidity and the one-or-two alternative give

    c_n=a_n+c_(n-1),       a_n=b_(n-1,2).              (2)

**Theorem 5.** The unlabelled bridge C=2 generating function is

    sum_n a_n z^n
      = z^5(1-z)(1-z^2)/(1-3z+z^2)^2.                (3)

At orders 5 through 11 the counts are

    1,5,18,60,191,589,1775.

All these graphs are rigid, so the labelled count at each order is n! a_n.

## 5. Exact total enumeration at every order above73

`C2_CUT_REDUCTION.md` proves that a bridge-free C=2 tournament is 2-strong and, by Busch's published Corollary 2, has order at most73.

Consequently **for every n>=74**, equation (2) counts every strong C=2 tournament, and every such tournament is rigid. At smaller orders the total unlabelled count is

    a_n + e_n,

where e_n counts the residual bridge-free atoms and is zero outside orders5 through73. The regular five-tournament gives e_5=1. The independently verified complete inventories through order8 give e_6=e_7=e_8=0. The remaining coefficients are not determined here.

Equivalently, the full C=2 class generating function differs from (3) by a polynomial of degree at most73. This is an exact all-order tail result, not an assertion that this polynomial is just z^5.

## 6. Verification

`verify_bridge_two.py` uses the frozen standalone graph constructors from `verify_start_two.py`. It checks the indegree lemma on all generated rooted start-two graphs, the degree-one equivalences on their C=2 outputs, rigidity, the two-vertex inverse, marked/unmarked multiplicities, and complete labelled inventories through order six. The earlier independently audited order-seven/eight census supplies the stated finite atom coefficients there.

The finite checks support the constructions and conventions. Lemmas 1–4, the canonical grammar, and the published order bound prove the all-order conclusions. The unrestricted `{7,21}` spectrum and the remaining bridge-free C=2 atoms are still open here.
