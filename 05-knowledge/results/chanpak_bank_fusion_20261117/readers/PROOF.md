# Fuse a reflected gap with an exact sevenfold reader

The growing reflected bank from the earlier two-flank work has a genuine fixed-order obstruction to a +2/equal-two-cover fusion. But changing the marked strong flank lifts the horizontal gap to the least compatible integer, six, and a finite rational-base carry loop gives exact equality of the two-cover profile at arbitrarily large orders. Inserting this pair into the prior sevenfold scaffold produces a new four-corner family; a second independent ordinal padding parameter supplies an exact prime-to-7p example.

## What a two-cover reader sees

Write P_k(T) for ordered covers of a tournament by k paths, allowing empty paths, so P_1=H. The odd-cycle pack identity represents it by disjoint directed odd cycles, with weight 2 for each nontrivial cycle and a factor k for each cycle or singleton. In particular every power of k has the parity of the tournament order n, giving P_2(T) == (-1)^n H(T) modulo 3. Consequently no pair of tournaments **of the same order** with H differing by two can have equal P_2, no matter how the reflected architecture is changed. At fixed order and equal P_2 any positive odd-count gap is a multiple of six.

To see and cancel the extra interface, for any marked block V and t in {+1,-1}, let L_t(V), R_t(V), Z_t(V) count nonempty paths beginning marked, ending unmarked, or both, respectively. Each path is weighted by t to its number of vertices times P_2 of the remaining induced block. Let Q_t(V) be the same unrestricted sum allowing the empty path; in particular Q_1=P_3(V). Write p(V)=P_2(V).

In X -> D -> Y*, where * reverses arcs and complements marks, append the apex beating precisely the marked vertices. Its P_2 is twice p(X)p(D)p(Y) plus twice the sum of the two spanning-triple-register Z_t values for t=+1,-1. Indeed a directed odd cycle containing the apex contributes 4 times P_2 of the leftover graph, and traverses an even-length path from a marked vertex to an unmarked vertex through these ordinal blocks; the parity projection is half the sum of the signed and unsigned reads. For an ordinal union A -> B, the full transfer for (p,Q,L,R,Z) is

    (p_A p_B, Q_A Q_B, L_A Q_B+p_A L_B,
         R_A p_B+Q_A R_B, Z_A p_B+p_A Z_B+L_A R_B).

Each formula separates a path within a block from a path crossing the forward boundary. Reflection * swaps L and R and preserves p,Q,Z. So on exchanging X,Y, the difference of P_2 (reflected minus original) is exactly

    2 sum_(t=+1,-1) [L_t(D)-R_t(D)] [L_t(X)p(Y)-L_t(Y)p(X)].

Use the original marked center D=(54,32,2,119,6,20,55), apex mark26 (outgoing masks, zero-based). Its two-cover-weighted L,R are 1397,1641 at +1 and -109,-109 at -1 (verified by direct weighted subset DP). Thus the difference simplifies to -488*(L_1(X)p(Y)-L_1(Y)p(X)). The same reflection changes H by 2*(h(Y)e(X)-h(X)e(Y)), where h=H and e counts full paths ending marked. This follows either by using complement weights P_1 in the same derivation, or from the earlier signed flank identity (the center's H-weighted unsigned ports are equal and the signed difference is -2). Both signs have been checked on literal graphs to fix the orientation.

## Rational-base carry loop with a physical gap-six marker

Take X to be the transitive address u followed by the four-vertex extreme-back-edge block E: transitive except 3 -> 0. Mark vertices 2,3 of E (mask12). Let Y be just the transitive address w. Directly, E has (p,Q_1,L_1,e)=(32,117,36,3). The replay checks every number, in particular that the end counts of vertices 2 and 3 sum to 3. With w even, e(Y)=0 and h(Y)=1; always H(reflected)=H(original)+6.

For a binary address with n bits b_i in their ordinary left-to-right order, put D_n(b)=sum_(i=0..n-1) b_i*2^i*3^(n-1-i). This is the marked start count L_1 on its transitive spine: after selecting starting vertex i, every later vertex independently belongs to the first path or to one of two leftover path colors, while every earlier vertex must be leftover. It is odd for every positive address (the leading digit is 1). By the ordinal transfer,

    L_1(X)/p(X) = (117 D(u)/2^n_u + 36)/32,
    L_1(Y)/p(Y) = D(w)/2^n_w.

Therefore equality becomes an exact rational-base equation, with u_j,w_j now read least significant digit first and padded with zeros:

    sum_j (32*w_j - 117*u_j)*(3/2)^j = 72.

The carry recursion starts at C_0=72 and reads one digit column by C_(j+1)=2*(C_j-32*w_j+117*u_j)/3, finishing at zero. Finite binary solutions force the intermediate carries to be integers: inductively the next expression has denominator at most three but, as a remaining rational-base sum, has denominator a power of two. They also stay bounded because the numerical carry map contracts by 2/3 with bounded translations. Thus the equation has a finite state inverse for this specific coefficient pair, and a literal productive cycle:

    u, least-significant first: 1 (1001)^k 1
    w, least-significant first: 0 (0000)^k 000001

The initial column sends 72 to 126; the four-column block returns 126 -> 162 -> 108 -> 72 -> 126; the final columns (padding the final u digit by five zeros) send 126 to 0. This proves feasibility for every k>=0. The lengths are 2+4k and 7+4k as required by the valuation: any equality for this marker forces n_w=n_u+v_2(32), since 117 and D(u) are odd. The length disparity means this construction does not solve the near-balanced source selector; it uses a physical flank away from that regime. Both reflected graphs have the same total order and are strong, since the apex beats the first left address vertex and loses to the last right address vertex, with the intervening vertices in ordinal direction.

The first addresses are u=3,w=64. Their ordinary counts are 1815977 and 1815983; both have P_2=4027766500. For comparison, the former +2 bank with mark0 has L_1(E_m)=3^(m+1)-2^m rather than 3^(m+1): the selected two-vertex set {0,last} is the exceptional reversed edge. For even w, the left side of the putative rational-base equality is then 0 modulo 3 while the right side is not. This is the particular carry version of the general same-order modulo-3 obstruction.

## Exact sevenfold fusion and a prime example

Use the earlier monotone ten-vertex quotient reader (750,16,18,6,41,78,158,318,127,510), slot9, with exterior reversals (1,8),(2,6),(3,6),(5,7),(3,7),(0,7),(4,7). Its eight profiles are (292,1),(791,117),(1150,147),(1307,195),(1315,273),(1315,335),(2017,19),(2044,7), reading alpha*P_2+beta*H. The replay independently checks the exact coefficients by TT1 and TT2: the slot's sole exterior entry is vertex 0, so any spanning path visits the slot in at most two runs. Each quotient is strong and so is every nonempty literal substitution. Since P_2>=2H, the stage increments are strictly positive for every actual input. Since beta is positive at all stages, the reflected input's horizontal +6 turns into +6*beta, always positive; the vertical exterior edits never touch any inner module arc.

Insert either reflected host, optionally with the same transitive ordinal prefix of size j, giving P_2 scaled by 2^j and H unchanged. This yields exact four-corner counts

    N, N+6, 7N, 7N+42,
    N=292*2^j*P_2(T_k)+H(T_k),   k,j>=0.

The j-padding does not preserve strength of the input module itself, but strength of the final quotient substitution is preserved, which is the claimed physical construction. For k=0,j=21, N=2466476862736151977 is prime and 7N=17265338039153063839 is physically attained. The proof of primality is Lucas's criterion with base 5, using the complete factorization N-1=2^3*3*8779*829097*14119373 and trial division for these factors. No claim of infinite prime quotients is made. The alternative earlier *induced* factor-seven reader 34*P_2+5*H yields +30,+210 horizontal gaps on the same bank, hence also gives a literal retained induced source for each upper corner.

Run `python replay.py` for independent small raw graph and weighted local signatures, exact carry equality, twelve pumped levels by the proved transfer, exterior reader replay and the elementary prime certificate. The result constructs new exact families while identifying the precise obstruction to carrying a native +2 unchanged through an equal-P_2 reader.
