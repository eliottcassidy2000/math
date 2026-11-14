# All new port positions in a Fermat seam

**Status:** exact four-state transfer with orientation boundary, literal strong atom examples and an all-order narrow-gap family. The shortest-letter scarcity theorem applies to this grammar, not to all tournaments.

Start with a tournament T rooted at v. Remember H=H(T), h=H(T−v), L equal to the relaxed two-piece response of the root, and E the number of full Hamilton paths ending at v. As in S328, adjoin x₁→…→x_j→y, oriented transitively within the chain, all arcs from T to the chain except y→v. This is attachment of R_j along v. A strong T stays strong. Its output count is F=H+2^jL, independently of which root is selected next. Write O for retaining v, T for selecting y, P for selecting x_j, and I_(j,k) for arbitrary new x_k. The empty deleted core has count one.

| next root | h' | L' | E' |
|---|---|---|---|
| old v | h | 2^(j+1)L−h | h |
| terminal y | H | F | H |
| interior x_k, t=j−k | H+2^(j−1)L when j≥2; H+L−E when j=1 | 2^(t+1)H+2^k3^t L−[k=1]E | 2^(k−1)L |

A term depending on the *first* interior position is essential, even for a long chain; the adjacent x₁–y boundary when j=1 makes that term especially easy to miss. For j≥2 the penultimate position k=j avoids it.

## Exact derivation of the interior row

Removing x_k leaves the same seam one length shorter if j≥2. If j=1 only y remains, so non-returning paths ending at y number H−E (the edge v→y does not exist), and paths using y→v number L. The full paths ending at x_k must use y→v, and exactly the k−1 earlier x vertices can be allocated freely to either side of the return; the later x vertices must precede y. This gives 2^(k−1)L. For the old root, a path ending v has to use the unique return y→v with all other vertices preceding y; this gives E'=h. Terminal end-paths number H, and the deleted terminal is ordinal with H(T−y)=H and no defect at y.

For L' at interior x_k let s=k−1, t=j−k. Start with the genuine paths through x_k counted by F. A relaxed two-cover that is not a path through x_k either has empty after-path and a before-path ending at a vertex beaten by x_k, or has both pieces nonempty, the before piece ending and the after piece starting among the t+1 later vertices. In the latter case the two pieces cannot both use the unique return y→v. Before the one explicit adjacency correction, the four defective-cover counts are:

| class | count |
|---|---:|
| empty after piece | H + 2^s(2^t−1)L |
| both pieces and no return | 2(2^t−1)H |
| both pieces and return in before piece | 2^s(3^t−2^(t+1)+1)L |
| both pieces and return in after piece | 2^s(3^t−2^t)L |

For the first line use the possible later terminal x_i or y, summing 2^(i−2)L for each x_i after deleting x_k. For the no-return case, the t+1 later vertices are split into two nonempty sets and the old core precedes the chain. If the before piece has the return, choose a nonempty subset of the t later interiors for the after piece; of those remaining, at least one must follow the return in the first piece. The early vertices have two choices; summing over nonempty subsets gives the third line. If the after piece has the return, each later interior goes to the first piece or to either side of the return in the other piece, with at least one in the first piece; the early vertices go to the first piece or after the return only. This gives the fourth line. The return-path core split is exactly the L counted at v.

The correction subtracts E exactly once if k=1: when the before piece ends in y without using a return and contains no other chain vertex, it cannot come directly from v. This case is possible precisely with no early vertices and all late x interiors in the other piece (empty if t=0). Summing the table, subtracting the correction, and adding F gives the stated L'.

## The shortest-block grammar is still sparse

Allow all three root choices O,P,T but use j=1 alone and a finite seed bank. There is an s<1 such that at most O(X^s) counts below X occur, including nonprime counts. Initially L≥H. After any step, writing r=L/H, one has 1≤r≤2, h/H≤2/3, E/H≤1/2 directly from the rows. A next O leaves r≥10/9; a next P leaves r=1+(1−E/H)/(1+2r)≥11/10. T resets r=1. Thus each two-step block increases H by at least 9 if the first root selector is T (three possibilities) and at least 48/5 otherwise (six possibilities). The inverse-scale sum at exponent 1 is 3/9+6/(48/5)=23/24<1; choose s<1 continuously, sum the geometric word weights, and count the finite prefix/suffix separately. This bounds distinct integers by word count. It does not purport to bound all genuine lengths at once.

## Ready states and exact odd-gcd gates

Just after selecting y the four-state is (H,e,H,e), where e is the count before that terminal move and e/H≤1/3. Call a ready pair primitive when gcd(H,e)=1. Integer transfer ensures any common divisor of a state divides all subsequent counts.

1. O_j then T_ℓ, with d=2^j and b=2^ℓ, returns ready

       ((1+d+2db)H−b e, (1+d)H).

   For primitive input its exact gcd is gcd(1+d,2H+e): since b is a power of two and (1+d)H is odd, take gcd with 2dH−e, which is coprime to H. For the genuine j=1,2,4,… the moduli 1+d are the pairwise relatively prime Fermat numbers. Hence at least one and indeed infinitely many j are safe for any primitive input. Repeating O_j,T_j with one safe fixed j remains primitive, because the new e is divisible by 1+d and the new H is coprime to it.

2. For j≥2 one can insert any number r≥0 of P_j's between O_j and T_j, and retain the same exact gate. On the two working coordinates (H,L) each P_j has matrix [[1,d],[2,d]] with determinant −d; the final T_j has determinant −d too. These do not change the odd gcd. Modulo F=1+d the first O projects to (0,−(2H+e)); P_j reduces to [[1,−1],[2,−1]], whose square is −I. After final T_j, even and odd r give (unit,0) and (0,unit), respectively, modulo F. The next return passes the same gate. This is an arbitrary-length physical primitive return with a fixed genuine length, not a theorem that its H is prime.

3. P_1,O_j,T_ℓ is a genuinely different gate. Just before the terminal it has

       (H*,L*)=((3+4d)H−d e, (8d−2)H+(1−2d)e).

   Its determinant before T_ℓ is 3−4d. For primitive input the new gcd can only divide 4d−3, and modulo that number the two pre coordinates are (3/4)(8H−e), (1/2)(8H−e). The final multiplier b is a unit at every odd divisor. Thus the exact gcd is gcd(4d−3,8H−e). For j=ℓ=1 the output H is 39H−8e. Starting ready (3,1) gives prime 109 with e'=31; starting ready (17,5) gives 623=7·89 with e'=177. The latter seed is obtained from the four-vertex source (2,4,9,3) with root 0, by T₁; both are checked on literal full graphs.

For either type of gate, an unsafe outcome from a primitive ready H≥3 cannot reach a prime or 7p later through this four-state grammar. In the fixed-j return of (2), the new common divisor is ≤1+d and 7 does not divide it (d modulo 7 is 2 or 4). If it were the prime p in 7p, the new H is already ≥(1+d)²H>7(1+d). In (3), similarly 7 does not divide 4d−3, and H_new≥3[(3+4d)H−d e]>7(4d−3) since e≤H/3. Future counts only grow. These are necessary gates, not sufficient prime tests.

## A same-order, gap-four primitive triple at unbounded heights

Take three order-six sources with root 2 and rows

    (2,8,3,21,39,15), (2,16,3,7,45,15), (2,4,41,19,7,27).

Their (H,L) are (15,29),(17,28),(19,27); T₁ produces three strong order-eight ready states (73,e) with e=15,17,19. For every genuine j≥4, apply O_j,T₁. The three outputs are

    365·2^j+35, 365·2^j+39, 365·2^j+43,

at the same strong order j+11, with common terminal deletion count 73(1+2^j), all primitive. The three odd-gcd gates involve 161=7·23, 163, and 165=3·5·11. None of 7,23,163,11 divides a Fermat number, since the order of 2 at any such divisor would have to be a power of two while p−1 has 2-adic part just 2, forcing order 2 (impossible for p>3). The factors 3,5 occur only in j=1,2 respectively by pairwise coprimality. All j≥4 therefore pass.

The replay compares literal graph/root profiles for 494 old and fresh placements and independently counts these small third-port and gap-four graphs. Longer lengths follow by the transfer, not a sampling extrapolation.
