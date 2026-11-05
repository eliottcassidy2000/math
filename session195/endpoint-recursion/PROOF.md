# Exact recursive classification of endpoint-one tournaments

Session 195, 2026-11-05. **PROVED elementary classification and transfer identities; FINITE-EXACT complete scalar census.** Independent review PASS in `../endpoint-recursion-audit/AUDIT.md`; the audit also checks all labelled tournaments through order six. The unrestricted Hamiltonian-path spectrum with exceptions `{7,21}` remains **OPEN**.

## 1. Inheritance and scope

The closest inherited mechanism is the rooted one-arc join law in the S193 `session193-joins/PROOF.md`, routed from the S194 guide. The decisive hostile is H=13: it is attainable, but not in the recursive class below. The corrected near miss is believing that an infinite physical returning class automatically covers the prime atoms. The underused coordinate is the endpoint count exactly one, together with the opposite endpoint and opposite two-path-cover response.

S194 `session194/ternary-fans/PROOF.md` requires an endpoint-one input for its coefficient-matched branches; its endpoint-floor theorem explains why some outputs cannot return. Here we classify **all** endpoint-one inputs, rather than selecting a bank. This yields actual recursive closure with a complete inverse decomposition and explicit arithmetic reachability failures. There is no literature-priority claim.

The live concepts are unique endpoints, source/sink strong components, one reversed cross arc, binary chain words, exact integer reachability, and atom coverage. The anchor is a returning physical classification; the niche is a complete all-order obstruction without a cited strong-order bound; the wildcard is a finite binary transfer implementation of arbitrary attachment words.

## 2. Coordinates and elementary facts

A directed Hamiltonian vertex word is counted once in its displayed order. Write H(empty)=1. For a tournament C with distinguished u, let

    h = H(C),
    e_u = number of Hamilton paths ending at u,
    s = number of Hamilton paths starting at u,
    L = sum_{S subset V(C), u in S} e_u(C[S]) H(C minus S).

Thus L counts ordered spanning two-path covers (P,Q), with P nonempty and ending at u; Q may be empty. A singleton has (h,e_u,s,L)=(1,1,1,1). Always 0<=s<=h and L>=h: split every Hamilton word immediately after u. This split is injective and lands in the covers counted by L.

Strong components of a tournament have a unique linear order, with all arcs pointing from earlier to later components. Hamilton words concatenate component Hamilton words. The only vertices with positive terminal counts are those in the sink component. In a strong tournament every vertex is an endpoint of a Hamilton path: take a longest path ending at the prescribed vertex. Any omitted vertex that cannot be inserted before its first vertex or between consecutive vertices is dominated by every vertex of the path. If every omitted vertex has that property, the path's vertices form a nonempty proper set with no incoming arc, contradicting strongness. Hence the longest path spans.

A tournament has H=1 exactly when it is transitive. One direction is immediate. A nontransitive tournament contains a directed triangle. Its three Hamilton words inject into those of the full tournament by inserting missing named vertices at the first available legal slot, in a fixed order. Therefore its H is at least three. The same insertion argument establishes induced monotonicity whenever needed.

## 3. Complete unique inverse decomposition

**Theorem 1.** Let (T,v) be a strong tournament with e_v(T)=1. Either T is a singleton, or the following decomposition exists and is unique on the labelled graph:

1. T-v is Z -> C, where Z is a transitive chain, possibly empty, and C is its sink strong component.
2. There is exactly one vertex u of C with u->v. It satisfies e_u(C)=1. Every other vertex of C is dominated by v.
3. The incidences of v with Z are arbitrary except that, when Z is nonempty, v must beat its first vertex.
4. If Z is empty, C is not a singleton.

Conversely, every construction satisfying these conditions from a singleton or a strong endpoint-one (C,u) produces a strong endpoint-one (T,v).

**Proof.** Let C be the sink strong component of T-v, and let P be the union of preceding components. Terminal-removal and SCC factorization give

    e_v(T) = H(P) sum_{u in C, u->v} e_u(C).

Every endpoint count in C is a positive integer, and the displayed expression equals one. Thus H(P)=1, so P is the transitive chain Z; there is precisely one summand, and it equals one. This proves 1 and 2 and makes C,Z,u unique. The chain order is unique because Z is transitive.

If Z is nonempty and its first vertex beats v, that vertex is a source of T, contradicting strongness. If v beats that first vertex, then v reaches every vertex of Z and C; every vertex of C reaches u and then v; every vertex of Z reaches u and then v. Thus the output is strong. If Z is empty and |C|>1, v reaches C-u directly, C is strong, and C reaches u->v, again proving strongness. If C is a singleton and Z is empty, the output is a two-vertex transitive tournament and is disallowed. Finally the same displayed endpoint formula proves e_v=1 for every allowed output. QED.

This is an exact recursive structural classification: the child core C is smaller, and recursive decomposition terminates at a singleton. Each forward step preserves every old core vertex and arc. Different histories may yield the same numerical state; no uniqueness of a tournament from its three numbers is claimed. At graph level the inverse core, root and word are unique.

For a general, possibly nonstrong tournament with an endpoint-one mark, all components before its sink component are transitive singletons and the marked sink component has endpoint count one and the same H. Thus the scalar support of all endpoint-one tournaments equals the scalar support of the singleton plus the strong class in Theorem 1.

## 4. Exact closed three-coordinate transfer

Write Z=(z_1->...->z_m), and encode w=b_1...b_m by b_i=1 iff v->z_i. For nonempty w define

    A(w) = sum_i b_i 2^(m-i),
    C(w) = sum_i (1-b_i) 2^(i-1),
    q(w) = 2^m,
    epsilon(w) = 1-b_m.

Let g(w) be the Hamilton count of the chain Z plus its binary apex v, ignoring C. It is computed exactly by starting (g,A,p)=(1,0,0) and appending a bit b:

    b=1:  (g,A,p) -> (g,2A+1,1),
    b=0:  (g,A,p) -> (g+A+p,2A,0).

For the empty word use g=1,A=C=0,q=1,epsilon=1. For completeness, appending bit 1 adds a universal sink, preserving g. Appending bit 0 reverses the new sink-to-apex seam: paths lost from the ordered join are precisely the old paths ending at the apex, counted by 1-p, while the new paths correspond to ordered two-path covers whose second path starts there, counted by A+1. Thus the increment is (A+1)-(1-p)=A+p. This also proves the inherited recurrence directly.

**Theorem 2.** The graph in Theorem 1 has e_v=1 and the following new state:

    w nonempty, b_1=1:
    h' = g(w) h - epsilon(w) s + (A(w)+1)L,
    s' = h,
    L' = (1+C(w))h + q(w)L.                       (1)

For the empty word, allowed only when h>1,

    (h',s',L') = (h+L-s, h-s, h+L).               (2)

In particular this state closes exactly under arbitrarily many returns, without losing the endpoint-one condition.

**Proof of h'.** Put D=Z+v. All D-to-C arcs point to C except u->v. Paths avoiding u->v are a D Hamilton word followed by a C Hamilton word, except those attempting the removed seam v->u. Their count is g(w)h-e_v(D)s. Here e_v(D)=epsilon(w), including the empty-word convention.

Paths using u->v split C into an ordered cover with first path ending at u, counted by L. They split D into an ordered cover with second path starting at v. Remove v from the second path: choose any subset of Z for that second increasing run, subject to its first vertex being beaten by v, or choose the empty subset. The number of such subsets is 1+sum_i b_i 2^(m-i)=A+1. This gives the remaining term in h'. The split is unique, so no path is omitted or counted twice.

**Proof of s'.** If Z is nonempty, a path starting at v must then traverse every Z vertex in order and then a Hamilton word of C. It cannot enter C before all Z are used, since C has no outgoing arc to Z; and it cannot return to v. The initial edge is legal precisely because b_1=1. Thus s'=h. If Z is empty it is any C Hamilton word except those starting at u, giving h-s.

**Proof of L'.** Delete v from the first path of a two-path cover ending at v. The singleton first path contributes H(Z->C)=h. A nonempty predecessor must be either u or a zero-bit z_i. The predecessor-u contribution is 2^m L: each chain vertex independently chooses one of the two path slots, and is placed in increasing order before that slot's C part. For predecessor z_i, its first path uses an arbitrary subset of z_1,...,z_(i-1), giving 2^(i-1) possibilities; all C vertices belong to the second path and contribute h. Summing zero-bit positions gives C(w)h. This proves (1) and also (2). QED.

### A concrete infinite growth subclass

Starting from C3 with (h,s,L)=(3,1,3), repeatedly use the empty-word step. The successive H values are

    3,5,9,17,31,57,105,193,... .

The common matrix is

    M = [[1,-1,1],[1,-1,0],[1,0,1]],

whose characteristic polynomial is t^3-t^2-t-1. Thus h_(j+3)=h_(j+2)+h_(j+1)+h_j for every j>=0. Every graph is strong, every new root has endpoint count one, and every old arc is retained. This is a proved nonconstant infinite physical return, not merely a numerically fitted sequence.

### Finite alphabet implementation, not finitely many numerical states

There are infinitely many integer states (h,s,L), as there must be for unbounded H. The arbitrary word macros have a finite-dimensional linear implementation over a fixed binary alphabet. Retain

    (h,s,L,U,X,Y,Q_h,Q_L,V,P_h,P_s)
    =(h,s,L,gh,Ah,AL,qh,qL,Ch,ph,ps).

Initialize U=h, X=Y=V=P_h=P_s=0, Q_h=h, Q_L=L. For b in {0,1}, set

    U'=U+(1-b)(X+P_h),
    X'=2X+b h,       Y'=2Y+b L,
    Q_h'=2Q_h,       Q_L'=2Q_L,
    V'=V+(1-b)Q_h,   P_h'=b h,   P_s'=b s,

and retain h,s,L. These are two fixed integer-linear maps. After a nonempty word beginning in 1, the return map is

    (h_new,s_new,L_new)=(U-s+P_s+Y+L, h, h+V+Q_L).

The empty-word map is (2). At the graph level, after creating v, read a bit by inserting the next z at the end of the current Z chain, before all C vertices, with the indicated v incidence. This uses only births, keeps the old core and root marks while the word is read, and forgets the old mark upon return. The transient empty-word graph from a singleton is weak, but the first allowed bit 1 makes the completed operation strong. No claim of a finite set of bounded numerical states is intended.

## 5. Exact all-order bounded reachability and an atom obstruction

For every allowed nonempty word the transfer obeys

    h' >= (A+1)L,                                  (3)

because gh-epsilon s>=0. Therefore every child with h'<=X has

    1<=A<=floor(X/L)-1.                            (4)

Canonical positive binary expansions enumerate exactly the allowed words beginning in 1. For the empty-word step on a nontrivial strong core,

    h'-h=L-s >= h-s >= |C|-1 >=2.                  (5)

The penultimate inequality uses positivity of all start counts in a strong tournament. Nonempty steps satisfy h'>=2L>=2h>h. Thus every recursive child strictly increases h, and (4) gives a finite exhaustive branch list at each parent.

**Theorem 3.** Starting from (1,1,1), applying (1)-(2), merging identical triples, and retaining only h<=X computes exactly the complete state set of all strong endpoint-one rooted tournaments with H<=X. This is an all-order certificate, not a vertex-order cutoff or a bounded-depth search.

Indeed the forward direction constructs an actual graph at every state. The inverse decomposition theorem and strict h decrease prove by induction that every possible state is reached. No large-order tournament or long hidden ancestry can escape this argument.

For X=13 the full state set is exactly

    (1,1,1), (3,1,3),
    (5,1,5), (5,1,7), (5,2,6),
    (9,1,9), (9,1,11), (9,1,15),
    (9,3,9), (9,3,11), (9,4,10),
    (11,1,13), (11,4,12).

Consequently **no tournament with H=13 has an endpoint-one vertex**, at any order. H=13 itself is attained by the strong tournament with outgoing masks

    [2,12,9,17,7],

whose terminal counts are [4,3,2,2,2]. This is the smallest admissible prime atom absent from the returning class. By taking converses, no H=13 tournament has a start-one vertex either. The classification is closed and nontrivially infinite, but cannot alone cover the target atom set, and cannot supply an endpoint-one presentation for every attainable scalar parent.

The reproducible X=1000 census gives 9,293 states and 475 distinct odd values. The exact 25 absent odds are

    7,13,21,41,43,55,67,77,83,163,173,231,233,
    247,251,277,299,421,449,503,533,569,583,637,821.

The independently checked census through 10,000 has 284,844 states and 4,971 distinct values. Its exact absent set is the preceding 25 values together with

    1067,1339,1861,2647.

These are absent from the endpoint-one class at every order. They are not asserted absent from the unrestricted Hamilton-path spectrum. Their all-order status follows from Theorem 3, not from an external order bound. In particular this supplies a precise limitation on any scheme requiring the current graph to keep an endpoint-one mark after every recursive return.

## 6. Verification and reproduction

Run `python3 verify.py` and `python3 -O verify.py`; compare stdout and `receipt.json` bytes. The verifier uses exceptions rather than removable assertions. It checks every labelled tournament of orders one through five, every endpoint-one strong root, unique SCC decomposition, graph reconstruction, and exact state transfer. Small counts and cover responses are checked by literal permutation/cut enumeration. It separately constructs all short binary macros from physical inputs, validates the finite-alphabet implementation, supplies every X=13 state with graph rows, and certifies the exact X=1000 and X=10000 closures using the proved exhaustive parameter bound.

A useful complementary witness is the one-apex word `1110101110`, which has H=511. It is obtained from a singleton by the macro w=`111010111`, with state (511,1,553). The receipt gives its eleven adjacency rows and verifies its count directly.

The independent `crosscheck_join.py` checks the other S195 worker's changing-root one-arc transfer on every labelled pair of tournaments of orders one through three and every root pair, 841 physical outputs, using terminal-removal subset DP. This cross-check is supplementary; its general formula is not a dependency of Theorems 1-3 here.
