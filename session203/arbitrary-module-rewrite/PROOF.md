# A shared-module rewrite with an exact variable correction

**PROVED, independently audited in the stated scope.** This is an all-order physical identity for arbitrary tournament modules. It is not a spectrum occupancy theorem, and no eventual threshold is claimed.

## 1. Statement and inherited physical cores

Let U be any finite tournament, including the empty tournament. Put h=H(U), with H(empty)=1. For a nonnegative integer q, let B_U(q) count ordered partitions into q labelled bins, where each occupied bin carries a directed Hamilton path on its assigned vertices and empty bins are allowed. Thus B_U(1)=h and B_U(0) is 1 exactly when U is empty.

Use the three independently audited core rewrites in the frozen dependency `dependencies/TRANSITIVE_REPAIR.md`. For a core Q with a mask M, Q[M;U] retains U literally as a module: every vertex of U beats the core vertices selected by M, loses to all other core vertices, and keeps all its original internal edges. The parent and intermediate child data are:

| Family | Parent core P | Intermediate child core T | Parent mask | Child mask |
|---|---|---|---:|---:|
| A | (2,24,3,37,13,23) | (2,8,67,101,15,87,19) | 15 | 79 |
| B | (32,5,1,39,47,70,31) | (0,5,129,39,47,199,159,27) | 5 | 5 |
| C | (2,4,41,3,15,91,31) | (0,5,41,131,15,219,159,23) | 27 | 155 |

All rows are outgoing bitmasks. Core vertices come first, followed by the vertices of U in any fixed order.

In A, reverse old core arc 1->4 and add vertex 6, beating old core vertices 0,1,4 and losing to U. Take the intermediate child as the final child R_A(U). This keeps every old vertex and changes exactly one old arc.

In B, reverse old core arc 0->5 and add vertex 7, beating old core vertices 0,1,3,4 and every vertex of U. In C, reverse 0->1 and add vertex 7, beating old core vertices 0,1,2,4 and losing to U. In both B and C, old vertex 0 is now a universal sink. Delete it to obtain the final child R_F(U). These are same-order replacements of one vertex: no edge between retained old vertices, including any edge within U, is changed.

**Theorem.** For each of A, B, C, and every finite tournament U, the displayed physical construction satisfies

    H(R_F(U)) = 2 H(P_F(U)) + H(U).                  (1)

The correction is the actual Hamilton-path count of the entire shared module, not an independently prescribed integer. Every parent and every final child is strong. Families B and C have the same order |U|+7 before and after replacement and satisfy

    C(P_F(U))=H(U),     C(R_F(U))=3H(U),             (2)

where C counts Hamilton cycles modulo rotation. Family A increases order from |U|+6 to |U|+7 and satisfies

    C(R_A(U))=2 C(P_A(U))-H(U).                      (3)

The construction, retained vertices, and exact scalar correction are all physical. Formula (1) is not permission to subtract counts to manufacture another tournament.

## 2. Why a transitive run calculation extends to arbitrary modules

Let A_r(U) count ordered partitions of all vertices of U into r nonempty directed paths. The paths are ordered and together use every vertex exactly once. Set A_0(empty)=1 and A_0(U)=0 for nonempty U. In particular A_1(empty)=0, while A_1(U)=H(U) for nonempty U.

For a fixed core and mask, delete U from a Hamilton word. This leaves one core order and a specified set of occupied gaps. Within every occupied gap, the removed vertices form a directed path in U. Conversely, choose a permitted core order, a permitted set of r occupied gaps, and an ordered r-path partition of U; insert the r paths into those gaps. Since all U vertices have the same core incidences, every joining edge is legal. These two operations are inverse.

Consequently, if a_r is the inherited core skeleton coefficient for r occupied gaps, then for arbitrary U

    H(Q[M;U]) = sum_r a_r A_r(U).                   (4)

No transitivity assumption on U is used here. The same argument works for cycles by fixing one core vertex as the cyclic starting point. Every spanning cycle contains that core vertex, so this removes precisely rotational multiplicity, even if U has its own internal cycles.

The relation to empty-allowed bin counts is

    B_U(q)=sum_(r=0)^q binom(q,r) A_r(U),
    A_r(U)=sum_(j=0)^r (-1)^(r-j) binom(r,j) B_U(j). (5)

The first identity chooses which q labelled bins are occupied; the second is its finite binomial inverse. These are identities on actual path partitions, not formal realizability assertions.

## 3. Exact correction, including the empty boundary

The independently audited path skeleton coefficient vectors are:

- A parent: (25,59,41,8); child: (51,119,82,16).
- B parent: (25,45,20); child before or after sink deletion: (51,91,40).
- C parent: (25,53,27); child before or after sink deletion: (51,107,54).

In every pair the child vector minus twice the parent vector is (1,1,0,...). Equations (4) and (5) give

    H(child)-2H(parent)=A_0(U)+A_1(U)=B_U(1)=H(U).

This proves (1) uniformly. When U is empty, A_0=1 and A_1=0. When U is nonempty, A_0=0 and A_1=h. Neither boundary is silently extrapolated from the other.

For explicit formulas write b_j=B_U(j). The same coefficients give:

| Family | Parent H | Final child H |
|---|---|---|
| A | h+17b_2+8b_3-b_0 | 3h+34b_2+16b_3-2b_0 |
| B | 5h+20b_2 | 11h+40b_2 |
| C | 27b_2-h-b_0 | 54b_2-h-2b_0 |

Thus the all-module statement retains the actual coupled coordinates h,b_2,b_3. They may not be selected independently. For example U=C3 has (h,b_2,b_3)=(3,12,33), giving parent/child path counts (471,945), (255,513), and (321,645), respectively. These scalar values are illustrative and carry no novelty claim.

## 4. Cycle counts and strongness

The inherited cycle skeleton vectors are A: (2,3,1) and (3,5,2); B,C: (1,1) and (3,3) for the parent and final child. Therefore

    C(P_A(U))=h+b_2,       C(R_A(U))=h+2b_2,
    C(P_B(U))=C(P_C(U))=h,
    C(R_B(U))=C(R_C(U))=3h.

They prove (2) and (3), including the empty module. Every tournament U has at least one Hamilton path, and H(empty)=1, so every displayed count is positive. A spanning directed cycle witnesses strongness of the corresponding whole graph. The intermediate children B,C are weak because of their universal sink; deletion preserves H by the bijection removing the forced last vertex.

## 5. A returned count ray and necessary coupling

For a nonempty m-vertex module U, splitting any Hamilton path at each of its m+1 cuts injects into the ordered two-bin covers. The original path and cut are recovered by concatenating the two resulting paths. Hence

    b_2 >= (m+1)h.

Also b_2 is even: swapping the two bins is a fixed-point-free involution when U is nonempty. In particular the parent values obey

    H(P_B(U)) >= (20m+25)h,
    H(P_C(U)) >= (27m+26)h.

These inequalities illustrate why the correction h and the parent count are not freely interchangeable.

Let U^+=U->singleton be the actual ordered join obtained by adding a universal sink. Then H(U^+)=h and B_(U^+)(2)=2b_2: in each two-bin cover the new sink can be placed at the end of either bin. More generally, ordered joining a transitive k-chain gives B_(U->Tr_k)(q)=q^k B_U(q), since each new chain vertex independently chooses its labelled bin and appears there in the unique increasing order.

For nonempty U the C formulas consequently give the exact count equality

    H(R_C(U)) = H(P_C(U^+)).                        (6)

This is a physical count-preserving rebase between two explicitly constructed graphs, not an isomorphism or induced extension. Their orders are |U|+7 and |U|+8, and their Hamilton-cycle counts are 3h and h respectively. The empty-module exception is genuine: R_C(empty) has count51, while P_C(singleton) has count53.

Thus for every fixed nonempty U, define U_k=U->Tr_k. The actual strong parents have

    N_k=H(P_C(U_k))=27*2^k B_U(2)-h,
    N_(k+1)=2N_k+h.

The same-order replacement R_C(U_k) realizes N_(k+1), and (6) supplies the next parent representation. This is an explicit returning count ray for every physical module U. The module's fixed count h remains attached throughout; it is not an induction covering arbitrary odd integers.

## 6. Verification, inheritance, and scope

The core graphs, their exact local edits, and their finite skeleton coefficients are inherited from the separately audited transitive-repair packet. Frozen proof, coefficient receipt, audit, and provenance are under dependencies/. The new deduction is that those same physical core rewrites supply correction H(U), with the displayed path and cycle formulas, for every tournament module U.

`verify.cpp` independently builds literal parent/intermediate/final graphs for every labelled tournament U of orders 0 through 4. It checks every edge and named edit, preservation of the retained graph, universal-sink deletion, all path and cycle formulas, the correction identity, the padding rebase, its empty-module hostile, and the elementary coupling bounds. Hamilton paths/cycles are counted by vertex-subset DP. The b_j coordinates are computed by direct subset convolution of induced-subgraph path counts, without using the core skeleton formulas. Normal and debug builds must agree.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug > DEBUG.json
    cmp NORMAL.json DEBUG.json

This is a shared-module correction operation, not a universal operation on arbitrary parents with a given count. It neither realizes arbitrary tuples (h,b_2,b_3) nor proves that every large odd target lies in one of the displayed images. In particular the free forward sector u of E(u) is not H(E(u)); substituting E(u) as U does not remove its correlated deletion contribution. Extracting that sector into a total count and proving an eventual threshold remain open.
