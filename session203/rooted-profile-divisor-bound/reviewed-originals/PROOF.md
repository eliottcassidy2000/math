# A sharper rooted-profile order cap from the published strong minimum

**PROOF CANDIDATE, pending independent audit.** This is an elementary corollary of an already cited theorem, applied to the physical rooted profile. No novelty or priority claim is made for the inherited strong-order bounds. It does not prove F(8,137) or an eventual spectrum threshold.

## 1. Published input and exact integer functions

The frozen `dependencies/BUSCH_CORE_BOUND.md` is the Session197 unrestricted-target audit. It records the verified primary source Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics13 (2006), N3, Theorem1 and Corollary1:

https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf

For a nontrivial strong tournament of order n>=3 the sharp minimum Hamilton-path count is

    f(3k)=3*5^(k-1),
    f(3k+1)=5^k,
    f(3k+2)=9*5^(k-1).

There is no strong order2. Put S(t)=max{n>=3:f(n)<=t}, taking S(t)=0 when the set is empty. This is computed with exact integers.

The same inherited packet proves the exact minimum g(N) among tournaments with no singleton strong components:

    g(0)=1,
    g(3k)=3^k                         (k>=1),
    g(3k+1)=5*3^(k-1)                (k>=1),
    g(3k+2)=25*3^(k-2)               (k>=2),
    g(5)=9.

Positive admissible orders are N>=3. Define G(t)=max({0} union {N>=3:g(N)<=t}). In particular G(1)=0. This is only an order upper bound when a specific count t is absent; it does not infer that every count below a maximum occurs.

## 2. Rooted quantities and reduction

For a nonempty rooted tournament (U,v), define h=H(U), s=start_v(U), d=H(U-v), and let L count ordered spanning two-path covers (P,Q) with P nonempty ending at v and Q allowed empty. H(empty)=1. Put B=h+L-s. The physical edge extension E(U) adds a vertex a with v->a and a->U-v. The inherited path-split identity is

    H(E(U))=B.

The elementary injections give L>=d>=s: designated singleton P=(v) gives d covers, and deleting v injects root-starting Hamilton paths into paths of U-v. Thus B>=h and B>=s.

Suppose s>0. The root belongs to the initial strong component C, since a Hamilton path must visit strong components in dominance order. Delete every singleton component after C, and call the resulting induced rooted tournament U'. This literal deletion preserves B and s. More precisely, if the remaining tail components have Hamilton-count product q, then

    B(U')=q B(C),  s(U')=q s(C),  H(U')=q H(C).  (1)

The cover identity follows because a path ending at v stays inside C, while the other path traverses all later components in order. In particular q divides gcd(B,s). The retained tail has no singleton components and has order at most G(q).

The case C singleton is exactly B=s. In one direction B(C)=s(C)=1. Conversely, if C is nontrivial, strong connectivity implies h(C)>s(C), since every vertex other than v also starts a Hamilton path; L(C)>=s(C) then gives B(C)>s(C), and (1) gives B>s.

## 3. The edge extension sharpens the strong-core bound

If C is nontrivial strong, E(C) is strong too: a reaches C-v directly, and reaches v through a path in C; every vertex of C reaches v and then a. It has order |C|+1 and Hamilton count B(C). Therefore

    |C| <= S(B(C))-1.                         (2)

This uses the actual edge extension, not merely B(C)>=H(C), and saves one vertex.

**Theorem.** Every realized positive-start pair (B,s) has an induced representative U' obtained by the stated singleton deletion, satisfying the following bounds.

If B=s, then

    |U'| <= 1+G(B).                           (3)

If B>s, define

    Q(B,s)={q positive : q divides gcd(B,s), S(B/q)>=4},
    K(B,s)=max_(q in Q(B,s)) [S(B/q)-1+G(q)].  (4)

Then Q(B,s) is nonempty and

    |U'| <= K(B,s) <= G(B)-1.                 (5)

If Q(B,s) is empty, the pair (B,s) is not realized. These conclusions are necessary bounds, not existence claims.

**Proof.** Formula(3) follows from the singleton initial component and the reduced tail. In the other case C is nontrivial. Its order is at least3, so E(C) has order at least4 and the actual tail product q lies in Q(B,s). Equations(1)-(2) and the tail bound G(q) give the first inequality in(5).

For the second inequality, fix any q in Q(B,s). By the sharp definitions, choose a strong tournament of order S(B/q) and count f(S(B/q))<=B/q, and a reduced tournament of order G(q) and count g(G(q))<=q; when G(q)=0 use the empty tail of count1. Their ordered join has no singleton components, has order S(B/q)+G(q), and has count at most B. By the definition of G(B), S(B/q)+G(q)<=G(B), proving the second inequality after subtraction of1. This step uses sharp minimum realizations only to compare order bounds, not to realize count B. QED.

There is also a direct explanation of the uniform bound: E(U') has nontrivial initial component E(C) and the same nonsingleton tail, so it is a reduced tournament with count B and order |U'|+1. Thus |U'|<=G(B)-1.

Combining with the elementary Session203 bound gives, when B>s,

    |U'| <= min(K(B,s), s+1+floor(log_2 B)).    (6)

Both bounds hold for the same literal singleton-deleted representative. No selection of incompatible alternative representatives is involved.

## 4. The first staggered base has a ten-vertex complete search cap

At base X=137 the F(8,137) candidates are

    B=137-2j,  s=r+3j,  r in {1,2}, 0<=j<=8.

All have B>s. The exact values of K(B,s), listed in increasing j, are

    r=1: 9,9,9,9,9,9,9,8,8;
    r=2: 9,10,9,9,9,9,10,8,8.                 (7)

Thus, if any of these eighteen pairs is realizable at any order, that pair has an induced representative of order at most10. The largest entries correspond to (B,s)=(135,5) and(125,20): the upper-bound divisor q=5 allows S(27)-1+G(5)=S(25)-1+G(5)=6+4=10. The table is a completeness cap, not a claim that the bound is attained for each pair.

This improves the direct elementary cap for this base, whose maximum is33. The already published physical certificates at base137 are unchanged. A smaller complete universe does not prove all-base occupancy; F(8,137) remains OPEN.

## 5. Verification and scope

`verify.cpp` checks literal profiles and SCC deletion on every labeled root-zero tournament through order6. It validates B>=h, B>=s, the B=s boundary, both new order bounds, preservation of (B,s), the edge-extension count, and strong connectivity of the extended initial component. It computes the eighteen-entry cap table using exact integers. Normal and unoptimized debug receipts must agree.

The proof depends on Busch's published strong minimum and the inherited component minimum; the finite controls do not replace those results. The contribution here is only their direct application to the coupled physical (B,s) certificate and its actual edge extension. No general occupancy, spectrum exclusion beyond an explicitly empty cap set, or eventual threshold is inferred.
