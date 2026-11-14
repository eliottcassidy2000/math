# A Fibonacci grammar for the unit-defect inverse branch

The [ordered-overlap analysis](BUSCH_INVERSE.md) reduces its defect-one branch to strong tournaments A obtained by adding a vertex v with only exit v→y to a strong tournament U in which exactly one Hamiltonian path ends at y. Here is a complete recursive classification of the pointed sources (U,y); it is not merely a finite census.

## 1. Terminal peeling with a marked predecessor

For a nontrivial strong U with one path ending at y, consider the SCC chain of U−y. A Hamilton path ending at y is obtained by appending y to a path on U−y. Only endpoints in the terminal SCC C are possible, and every vertex of C is a possible endpoint. Each available endpoint contributes the product of the path counts of the earlier SCCs times its endpoint count inside C.

Because U is strong, some vertex of C beats y. The end count equals one if and only if:
1. every earlier component is a singleton, so U−y=TT_b→C;
2. exactly one p in C beats y;
3. p is itself the end of exactly one Hamiltonian path of C.

On the prefix TT_b, the arcs to y are arbitrary except that, when b≥1, y must beat the first vertex. This condition is necessary (otherwise the first vertex is a source of U) and sufficient: y reaches the whole prefix and hence C, and every vertex can reach p and then y. If b=0 and C is strong of order at least three, y beats the rest of C and loses only to p, so U is strong as well. C a singleton with b=0 would make order two and is excluded.

Conversely all these constructions yield a strong rooted (U,y) with exactly one path ending at y. Peeling y recovers the ordered singleton prefix, the pointed terminal (C,p), and the bit word of the new arcs. The recursive code therefore identifies rooted isomorphism classes exactly. The initial rooted singleton is allowed as the terminal base, although it cannot itself be used as U for a strong defect-one extension A.

## 2. Enumeration, and a discrete exact solver

Let u_n count the strong rooted order-n classes with unique Hamiltonian end, using u_1=1, u_2=0. At n≥3, b=0 contributes u_(n−1); when b=n−1−c≥1 there are 2^(b−1) words on the prefix for each rooted C of order c. Thus

    u_n = u_(n−1) + Σ_(c=1)^(n−2) 2^(n−c−2) u_c,
    u_n = 3u_(n−1) − u_(n−2)   (n≥4),
    u_n = F_(2n−4)              (n≥3; F_0=0,F_1=1).

The order series is z+z³/(1−3z+z²). There is no isomorphism screening in the decoder. A rooted class is rigid as a rooted object because a root-fixing automorphism fixes its unique ending path pointwise.

For exact prescribed source H in the Busch D=1 branch, enumerate this grammar only through orders |U|≤3 log_5 H (since |A|=|U|+1), attach v, and retain exactly the H matches. The finite solver is complete for that branch; the path count can be computed by subset dynamic programming. This replaces enumeration of all tournaments by an exact Fibonacci-sized candidate bank, not an assertion that the bank hits every H.

## 3. Transfer coordinates and why the endpoint is real

Let W(C,p) count two-slot covers whose second path starts at p; W(C*,p) is its dual. For the b=0 step, the unique-end mark y has the actual profile

    H(U) = H(C) − S(C,p) + W(C*,p),
    W(U,y) = Q_C(2) − W(C,p).

The start count S(C,p) is not necessarily one. The first formula uses the one-incoming version of the one-exit count. The second deletes y from the beginning of the second path; only a nonempty second path starting p is forbidden. Therefore after adding the one-exit v to U, the resulting unit-defect source satisfies

    H(A) = H(C) + Q_C(2) + W(C*,p) − W(C,p) − S(C,p) − 1.

This is a count on a literal two-step core, with the signed start/end difference retained. In particular even the exact Fibonacci grammar does not let us substitute only the integer H(C). The rooted inverse analogue of an isothermal net is the unique terminal path and the bits at its changing prefix, not the total path count alone.

## 4. When can the same root have both unit endpoints?

The peeling formula also gives the number of Hamiltonian paths starting at y. If b≥1, y beats the unique initial singleton of U−y and every path from that singleton finishes through C, so the start count is exactly H(C). If b=0, the possible first vertices after y are all vertices except p, hence the start count is H(C)−S(C,p). The latter is at least |C|−1≥2 by the spanning-cycle cuts of C.

It follows that a nontrivial strong rooted tournament has both endpoints equal to one iff C is a singleton and b≥1 in the peeling grammar; equivalently U−y is transitive, with y beating its first vertex and losing to its last. Thus the two-sided-unit class has exactly 2^(n−3) rooted isomorphism classes at order n≥3, compared with F_(2n−4) for the one-sided class. The second endpoint requirement changes the growth base from (3+√5)/2 to 2. In particular it cannot be silently added to a one-sided recursion as an equivalent constraint.

## 5. Full prefix transfer, not just the empty-prefix case

For the arbitrary prefix P=TT_b let L_P,H_P,E_P be the one-way overlap coordinates of the small tournament on P and y. If the bits are u_i=[y→i], then L_P=1+Σ_i u_i 2^(b−i), and E_P=[u_b=0] for b≥1. At b=0 put L_P=H_P=E_P=1. H_P is the ordinary one-apex-over-chain count; it can be computed directly by summing weighted inversions of x_i=1−u_i as in [Session 323](../session323_computability_and_native_descent_20261113/PAIRED_END_DESCENT.md) (the cycle-packing formula there works without requiring an active cut).

Write h=H(C), S=S(C,p), W=W(C,p), W*=W(C*,p), Q=Q_C(2). Then the exact transfers are

    H(U)=L_P W* + H_P h − E_P S,
    W(U,y)=L_P Q − W,
    H(A)=L_P (W*+Q) + H_P h − E_P S − W − 1.

For H(U), use the four terms of the seam formula: on the C side there are W* covers with nonempty legal first part, all also legal on the second side since p cannot be used twice; the number with first part empty and legal second is h−S. On the prefix side h_P=1. For W(U,y), deleting y from the start of the second path leaves a cover of P→C: if the second path's P part is nonempty and lawful, its initial vertex is already fixed by P, so all Q options on C work. If that part is empty, Q−W options remain. Finally A has H(U)−1+W(U,y) because y is its unique exit neighbor and has unit end count.

Thus each retained core gives a target-side equation in one binary prefix word; its coefficients and parity are known before the prefix is chosen. The coefficients of L_P and H_P are positive, but H_P is not in general monotone in the binary value L_P. No binary-search shortcut is assumed.
