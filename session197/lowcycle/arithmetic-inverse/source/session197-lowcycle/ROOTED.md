# A rooted low-cycle gluing calculus

Session 197. **All-order physical identities and a recursive classification of the rooted start-one class.** This supplement does not classify all tournaments with one or two Hamiltonian cycles.

## 1. Inheritance and notation

The scalar one-seam formula and the returning three-coordinate state below are inherited from the explicitly linked [S195 rooted-interface proof](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/rooted-interfaces/PROOF.md). This supplement adds the cycle and root-start responses, an all-order SCC characterization of the start-one class, its complete arbitrary-gadget compiler, and the resulting bounded image. No literature-priority claim is made.

For a rooted tournament (A,a), let h=H(A), e=e_a(A) be its number of paths ending at a, and s=s_a(A) its number of paths starting at a. Let R_a(A) count ordered spanning path covers (P,Q), where P may be empty and Q is nonempty and starts at a. Define L_a(A) dually: P is nonempty and ends at a, while Q may be empty. A directed path of one vertex is allowed, and H(empty)=1.

Given disjoint nonempty rooted tournaments (A,a),(G,b), orient all cross arcs A->G except the single reverse arc b->a, and move the root to b. Denote the actual result K(A,a;G,b).

Write the gadget constants

    g=H(G), t=s_b(G), E=e_b(G), L=L_b(G), Q=R_b(G),
    d=H(G-b), B=P_(G-b)(2), delta=1 if |G|=1 and 0 otherwise,

where P(2) counts two ordered path slots with empty slots allowed.

The inherited transfer is

    h' = g h - t e + L R,
    e' = E h - delta e,
    R' = Q h + B R.                                  (1)

The scalar identity follows from the disjoint cases avoiding or using b->a. Avoiding it gives an A-path followed by a G-path, with exactly e*t forbidden joins. Using it gives the unique run pattern A_before,G_before,A_after,G_after, with A_after starting at a and G_before ending at b, hence R*L choices. This is recalled for clarity, not claimed as a new law.

## 2. Two additional exact responses

**Proposition 1.** The new root-start count is

    s_b(K)=d*s_a(A).                                 (2)

For total order at least three,

    C(K)=s_a(A)*e_b(G).                              (3)

**Proof of (2).** A Hamiltonian path beginning at b must immediately take b->a. Otherwise it enters G-b and can never reach A, because b->a is the only reverse cross arc and b has already been used. After a it traverses all A in a path starting at a, then all G-b in an arbitrary path. It cannot leave A and subsequently return. Every such pair is valid, giving (2).

**Proof of (3).** Every spanning cycle uses b->a, the only way back from G to A. Starting just after this arc, it consists of an A-path starting at a and a G-path ending at b. The other cross boundary is legal. It could be the missing a->b arc only if both pieces were singletons, which is precisely the excluded order-two case. This gives a bijection and (3). QED.

For order at least three, K is strong exactly when the product in (3) is positive. Positivity directly supplies a Hamiltonian cycle. Conversely, a strong tournament has a Hamiltonian cycle, so the product is positive. In that event b->a is a strong bridge: deleting it leaves no directed route from G to A.

The cycle count is not determined by the old cycle count. The extra response s, rather than C(A), is the correct gluing interface.

## 3. Exact recursive classification of a unique starting path

**Theorem 2.** Every nontrivial strong rooted tournament (T,r) with s_r(T)=1 is obtainable from the singleton by repeated one-seam gluings whose new gadget G satisfies:

- G-b is transitive;
- e_b(G)=1.

Conversely, every output of this grammar of order at least three is strong and has C=1 and s_root=1. The only smaller outputs are rooted transitive tournaments of orders one and two.

**Proof of completeness.** Decompose T-r into its uniquely ordered strong components. Let U be the initial component, and let q be the product of the Hamiltonian-path counts of all later components. Every Hamiltonian path starting at r first enters an out-neighbor of r in U. Thus

    s_r(T) = q * sum_(a in U, r->a) s_a(U).           (4)

A vertex of a strong tournament starts at least one Hamiltonian path, by cutting a Hamiltonian cycle. This also holds for a singleton component. Every summand is therefore positive. A nontrivial strong component has at least three Hamiltonian paths, since it contains a directed triangle and vertex insertion injects that triangle's three paths into its paths.

Equation (4) equal to one forces q=1, so all components after U are singletons and form a transitive chain W. It also forces exactly one out-neighbor a of r in U, with s_a(U)=1. Consequently T is literally K(U,a;G,r), where G=T[W+r]. If W is nonempty, its last vertex must point to r: it loses to U and all earlier W vertices, so otherwise it would be a sink in the strong tournament T. A Hamiltonian path of G ending at r is uniquely the increasing W-path followed by r. If W is empty, G is a singleton and e_r(G)=1 as well. This proves the gadget conditions.

The old piece U is strong or singleton and is smaller, so repeat until reaching the singleton.

**Proof of soundness.** Every gadget has d=1. Starting with s=1, (2) preserves s=1. Every gadget also has E=1, so (3) gives C=1 whenever total order is at least three. This directly proves strongness. QED.

This is a complete classification of the explicitly rooted class, not an inference from a finite list. Whether every C=1 tournament admits such a root is a separate question not settled here.

### A useful graph-theoretic equivalence

For a strong tournament with C=1, the following are equivalent:

- it has a strong bridge;
- it has a vertex of outdegree one;
- it has a root with s_root=1;
- it belongs to Theorem 2's rooted grammar after choosing a suitable root.

Here is a proof of the only additional implication needed. A nontrivial rooted strong tournament (B,b) with e_b(B)=1 has an outdegree-one vertex different from b. To see this, decompose B-b into its ordered strong components. The terminal component U has exactly one vertex a pointing to b with e_a(U)=1, and all earlier components are singletons. If U is a singleton, that vertex points only to b and has outdegree one. If U is nontrivial, apply induction to (U,a): its outdegree-one vertex v is different from a, loses to b, and loses to every earlier component, so it remains outdegree one in B.

Now split a strong bridge b->a into nonempty sides A,G with all cross arcs A->G except that bridge. Such a split exists by taking the vertices reachable from b after deleting the bridge as G and their complement as A. Equation (3) and C=1 give e_b(G)=1. In G's ordered SCC decomposition, b lies in the last component and all earlier components are singletons. If that last component is a singleton, b has no internal outgoing arc and gains only b->a in the full tournament, so it has outdegree one. Otherwise the preceding lemma supplies an outdegree-one vertex different from b in the last component; its outdegree stays one in the full tournament.

A vertex of outdegree one has its unique outgoing arc as a strong bridge. If its unique out-neighbor is u, every Hamiltonian path from that vertex closes back to it, since every other vertex points to it; therefore its start count equals C. Finally Theorem 2's decomposition exposes a strong bridge, giving the other directions.

## 4. A binary integer compiler for every allowed gadget

Write a one-apex gadget G as a root b and a transitive base x_1,...,x_m. Let bit w_i=1 mean b->x_i. The condition E=1 means the final bit is zero when m>0; m=0 is the singleton gadget.

The root-ending two-cover count is

    L = 1 + sum_(i:w_i=0) 2^(i-1).                  (5)

Indeed choose the base vertices in the first path ending at b. If nonempty, their largest index i must point to b, and the earlier chosen vertices are arbitrary. The remaining base vertices form the second increasing path. The empty base selection contributes one.

Because the final bit is zero, (5) is a bijection between all allowed gadgets and the positive integers L: L=1 is the singleton; for L>=2, write L-1 in binary from least significant to most significant bit and complement each bit to obtain w. In particular, the gadget has base size m=bit_length(L-1).

Let g(L)=H(G), Q(L)=R_b(G), and B(L)=2^m. The starting values are g(1)=Q(1)=B(1)=1. For ell>=1,

    g(2ell)=g(ell), Q(2ell)=Q(ell), B(2ell)=2B(ell).

For ell>=2,

    g(2ell-1)=g(ell)+ell-(ell mod 2),
    Q(2ell-1)=Q(ell)+B(ell), B(2ell-1)=2B(ell).        (6)

Also t=L mod 2 and delta=1 exactly at L=1. Thus the complete physical transfer is

    h' = g(L)h - (L mod 2)e + L R,
    e' = h - [L=1]e,
    R' = Q(L)h + B(L)R.                              (7)

**Proof of the recurrences.** Prefixing a zero base bit adds a universal source, so g and Q do not change; L doubles and B doubles. Prefixing one gives L'=2L-1 and adds B to Q. For H, split by whether the new first base vertex lies before or after the root in the Hamiltonian path's two increasing base runs. In the first case the old left run must be nonempty, giving g-t choices. In the second case the new vertex satisfies the right-run first-edge condition automatically, giving exactly L choices. Hence g'=g-t+L. For nonempty words, t is their first bit, equivalently L mod 2. This proves (6) by direct path counting. QED.

The integer compiler produces actual gadgets, not merely formal matrices or arbitrary integer profiles.

## 5. Exact bounded image with no gadget-size or depth cutoff

The inherited inequalities R>=h>=1 and e<=h give, for every allowed gadget,

    h' >= L R >= R,      R' >= h+R.                   (8)

The first follows because g h-t e>=0. The second follows because Q,B>=1. Also h' cannot be smaller than h, either by the inherited transfer inequalities or by induced-path insertion.

For a bound X, explore all states (h,e,R) with h<=R<=X, starting at (1,1,1). For each state try exactly the integers 1<=L<=floor(X/R), compute (7), and record h' when h'<=X. Queue the new state only if R'<=X and it has not been visited.

There are finitely many states. Every discarded larger L has h'>X by (8). Every valid preterminal state of a construction ending at H<=X has R<=X, because the next h is at least R; h never decreases. A terminal state with R'>X is recorded before being discarded. Hence this computes **exactly** all bounded H-values of the entire arbitrary-gadget grammar. Theorem 2 identifies that grammar, at orders at least three, with all strong tournaments having a root-start count one.

The independently implemented `unique_start.py` gives:

- X=501: 2,727 retained states, 3,900 admissible transitions, maximum discovery depth 9.
- X=5001: 81,585 retained states, 117,980 admissible transitions, maximum discovery depth 13.

The missing odd values from 3 through 5001 are exactly

    7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,
    277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647.

These are all-order exclusions from this rooted class, not global Hamiltonian-path gaps. In particular H=67 is excluded here for an independent structural reason; the stronger unrestricted low-cycle exclusion in `PROOF.md` uses Busch plus a complete order-eight census.

The small coordinate s=1 is physically maintained, but does not enter the three transfer rows (7). It therefore does not repair their content obstruction. For example using L=2 twice gives (h,e,R)=(9,3,9) while s=1. All future h-values from that marked state under (7) remain divisible by three.

## 6. One Hamiltonian cycle with arbitrarily large feedback vertex number

Take L=3, the directed triangle gadget, at every step. The matrix is the inherited

    J = [[3,-1,3],[1,0,0],[3,0,4]].

From the singleton, let T_t be the result after t steps. For t>=1:

- T_t is strong, has order 1+3t, and has exactly one Hamiltonian cycle;
- its feedback vertex number is exactly t;
- its first H-values, starting at t=0, are

      1,5,35,229,1483,9605,62219,403045,...;

- for t>=3,

      h_t = 7h_(t-1)-4h_(t-2)+4h_(t-3).             (9)

The cycle assertion follows from (2)-(3), since d=E=1. Each added triangle is an induced block disjoint from the other added triangles, forcing every feedback vertex set to meet all t blocks. Conversely delete the new marked vertex from each added triangle. This removes every backward seam and leaves two vertices from each block in forward transitive order, together with the original singleton. Thus t deletions suffice. The recurrence is Cayley-Hamilton for `det(zI-J)=z^3-7z^2+4z-4`.

The graph family and three-state matrix were already constructed in S195. The precise additional conclusion here is that this entire unbounded-feedback family is **uniquely Hamiltonian**. It shows why a fixed two-apex description cannot classify all low-cycle tournaments.

## 7. Finite structural evidence, not a theorem

A separate complete generation through order eight, followed by all 1,761,280 one-vertex extensions of its 6,880 representatives to order nine, gave the following C=1 class counts at orders 3 through 9:

    1,1,3,8,21,55,144.

The corresponding C=2 counts are

    0,0,2,5,18,60,191.

Every C=1 example had a strong bridge. The only bridge-free C=2 example was the regular order-five tournament, with H=15. The apparent Fibonacci pattern and global bridge classification remain **finite evidence only**. Neither is used to prove the all-order results above or the H=67 classification.

The publication verifier for the main H=67 packet reconstructs the complete orders through eight. The order-nine exploratory census is not needed by any asserted theorem in this supplement.

## 8. Reproducible physical checks

Run from this packet:

    python3 verify_rooted.py > rooted.normal.stdout
    python3 -O verify_rooted.py > rooted.optimized.stdout
    cmp rooted.normal.stdout rooted.optimized.stdout
    python3 unique_start.py 5001

The standard-library checker imports only this packet's `verify.py` and `unique_start.py`. It performs 255 integer-compiled gadget versus adjacency-DP checks, 1,609 actual one-seam gluing checks (including gadgets with nontransitive root-deletion and hence nonunit start multiplication), 36 complete recursive decompositions across all rooted strong classes through order six having start count one, four literal triangle-family checks through order 13, and the two closed bounded-image searches. The family checks include actual disjoint induced triangle blocks and an explicit feedback set leaving a transitive graph.

Normal and optimized runs agree byte-for-byte. The proofs establish the all-order statements; the finite tests are independent physical checks of their implementation. The primary H=67 census is a separate verifier, so its result does not depend on the start-one classification or its arithmetic compiler.
