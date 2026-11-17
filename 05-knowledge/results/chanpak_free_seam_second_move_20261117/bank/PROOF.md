# Three physical translations with two independent flanks

Take the nine-vertex center R=(52,365,176,69,266,280,437,315,13) in outgoing-mask notation. This is the earlier single-exit three-state rotor. Its stronger property is equality of the **full signed two-flank profile**: for apex tags 113,43,39 respectively,

    (h,P,L,R,Lminus,Rminus,kappa)
      = (1521,26412,10422,12894,-750,-890,7285+2j), j=0,1,2.

Here h and P are the center's ordinary Hamilton-path and ordered two-cover counts, L sums the ordered two-covers whose second path starts tagged, and R sums those whose first path ends untagged. The signed versions weight the distinguished path by (-1) to its number of vertices. Kappa is the Hamilton-path count of the center together with the tagged apex. The companion standard-library replay computes these numbers by a subset DP with endpoint incidences; it also checks 24 literal assembled tournaments covering all four address parity regimes.

With the two-flank rule (left binary address as apex out-tags, right reversed complemented address as apex out-tags), the inherited signed formula at any positive addresses u,w is

    H_j(u,w) = kappa_j + h(F(u)+F(w)-2) + R*u + L*w + P*u*w
                         - Rminus*[u odd] - Lminus*[w odd].

In this formula F(0)=1, F(2x+1)=F(x), and F(2x)=F(x)+x+(x mod 2) for positive x. In particular all F values are odd, and F(x)<=x+1 by induction (for odd x one can use the stronger bound (x+1)/2). Thus the two exact +2 moves are repeated by retagging the same center, with no change in either of the addresses or any continuation coordinate other than kappa. Unlike the prior signed reflection switch, these are not swaps of the two addresses: the unsigned ports are unequal. The moves produce three consecutive odd counts jointly, on equal-order strong tournaments. For strength the apex beats the first left vertex and loses to the last right vertex; both flanks and every vertex of the middle lie between these in the ordinal order and hence in the apex's SCC.

For every fixed modulus M and every large prescribed **odd** N, the bank has a three-count progression with a term H congruent to N modulo M and with absolute error O_M(N^(2/5)). The term need not be the lower term; j is selected to repair the mod-3 phase when that is necessary. This is stronger than either an absolute approximation with no phase or an unrestricted congruence with unbounded error. It does not show H=N, and the implicit constant is not uniform in M.

## Complementary suffix steering

Let x be odd and y be even. Write Q=x+y+1, Z=x+1, and
f=F(x)+F(y)-Q+Z. Append a word to x and the bitwise complementary word to y. Use the following three choices of word on x, all ending in 1. Directly from the F recursion their effects are

    A = 1   : (Q,Z) -> (2Q,2Z),     f -> f
    B = 01  : (Q,Z) -> (4Q,4Z-2),   f -> f + 2Z-Q
    T = 001 : (Q,Z) -> (8Q,8Z-6),   f -> f + 6Z-3Q-4.

They maintain x odd, y even. Modulo the odd part m of M they are invertible. An inverse is implemented as a positive power of the same word on this finite state space. Let d be the order of 2 modulo m (d=1 when m=1). In right-to-left composition order, C=A^(-2)B, D=A^(-1)CA, E=A^(-3)T all fix Q and move Z by -1/2,-1/4,-3/4 respectively. Their f increments are 2Z-Q, 4Z-2Q, 6Z-3Q-4. Hence V=D^(-3)E fixes Z as well and changes f by -6Z+3Q-1. Move first to Z=Q/2 with C, control f there by V (step -1), and translate by C to the desired Z, compensating its known f drift. All these regulator words have lengths divisible by d. For inverse C,D,E, their whole-map inverses as positive powers retain this length divisibility. This steers the two coordinates freely with Q fixed modulo m.

Suppose x+y+1 is initially 2 modulo 4. The two first choices AA and B have equal length, but their f values differ by 2 modulo 4 since Z is even. Choose this bit, then append enough A's to put Q in 0 modulo 4; this does not change f. Every subsequent A,B,T keeps f fixed mod 4 since Q is 0 mod 4 and Z even. If 2^s is the dyadic part of M, the choice provides the required mod-4 class for any admissible even final f (for s<=1 there is no choice to make). Achieve higher dyadic bits with W=1^r001, r>=max(s,2), d dividing r+3. Its first run of ones resets Q,Z modulo 2^s, so W subtracts 4 from f modulo 2^s; modulo m it fixes Q, translates Z by -6, and adds an affine function of Z to f. Therefore W^m is the identity on all three odd-residue coordinates (sum the arithmetic progression, using m odd) and subtracts 4m from f modulo 2^s. This generates every increment divisible by 4. Between odd steering and any prescribed final suffix we then append a long A-guard of length divisible by d: Q and Z are reset modulo 2^s, but f and the odd state are preserved.

Fix a final s-bit suffix ending in 1 on x (an empty suffix if s=0); its complement is appended to y. It fixes x mod 2^s and its affine action on (Q,Z,f) is invertible modulo m, because each bit has diagonal (2,2,1) and previous parity is known. Its dyadic drift in f from reset Q,Z is also known. Thus any desired (x mod m,F(x)+F(y) mod M) with the sum even at the dyadic part can be imposed, while the word length is held to a **constant L depending only on M**. Indeed odd steering routes have lengths 0 mod d; the initial toggle has fixed length 2; odd-state neutral guards, W^m and padding A^d have lengths 0 mod d. There are only finitely many starting residues and requested final residues. After the reset guard and before the final suffix, pad with A^d as needed to give all routes the same total length. One chooses the target odd-state Q before the final suffix according to the total length, but this is not circular: the total length is fixed modulo d before padding, and Q only doubles at each step, so its final odd residue is already determined. Require L>=s. The odd steering and the high dyadic adjustment do not disturb the chosen mod-4 class. For the fixed final suffix, its drift mod 4 after the reset is independent of the toggle choice.

In ordinary integers this keeps the sum *exact*: with D=2^L and fixed initial prefixes X,Y the new numbers lie in

    x = D X + e,   y = D Y + (D-1-e),    x+y = D(X+Y)+(D-1).

That is the key difference from independently changing both addresses. A fixed-sized adjustment of x and y in the same direction can cause a much larger product error.

## Choose a congruent quadratic line from the target

Use a power of two V comparable to N^(1/10), with V>=max(4,2^(s+1)), and set u=Vx+1, w=Vy. The fixed suffix expansion of F and the parity x odd, y even give

    H_j = A*x*y + B*x + C*y + K_j + h(F(x)+F(y)),
    A = P*V^2, B=R*V+h(V/2-1), C=(L+P)*V+h(V-1),
    K_j=kappa_j-h+R-Rminus.

In the definitions of B,C, the letters R,L still denote the **signed-profile unsigned ports**, not binary lengths. The use of x odd, y even accounts for the +h inside K_j relative to expansion with no prefix parity contributions.

First select modular (x,y,Fsum) so this read is congruent to N. For an odd prime dividing M other than 3 or 13, h is a unit, so pick x,y at will and use Fsum. At p=13 take w=0 and use u, whose coefficient R is a unit there. At p=3 also take w=0; R has valuation one, and every variable coefficient of the original formula is divisible by 3. Select j from the three kappa values to match N modulo 3, then vary u to solve every higher digit. Translating the read between (u,w,F(u)+F(w)) and (x,y,Fsum) is invertible modulo odd primes, by the suffix formula and V being a power of two. At the dyadic part x+y≡-1 because L>=s, so the displayed read reduces to K_j+h(1+Fsum) modulo 2^s, and every odd residue is possible with even Fsum. Combine the odd and dyadic requirements on Fsum by CRT. The dyadic address suffix is arbitrary with x odd; choose any fixed one.

Let S denote the exact line total x+y. It must match the chosen x+y mod m, and be of form D S0+(D-1) with S0≡1 mod 4, so that compatible prefix parities X odd, Y even can be chosen and Q starts 2 mod 4. Since D is a unit modulo m, these conditions are one progression of spacing 4Dm. Define p=C/A, q=B/A, and E0=K_j-BC/A. Choose S as the **least** such line at or above 2sqrt((N-E0)/A)-p-q. Then on x+y=S the uncorrected quadratic A(x+p)(y+q)+E0 equals N at a real point with x,y comparable to S. The square root of the quadratic discriminant, in address units, is O_M(sqrt(S)) because S exceeds the threshold by at most 4Dm. Choose an odd prefix X within O(1) of x/D there, and Y=S0-X is necessarily even. They are positive for all sufficiently large targets. Use the complementary suffix controller with these prefixes and the selected residue of x and Fsum. The final shift of x from the real point is O_M(1), and that of y is its negative. The simultaneous constraints on y follow from the exact line total.

Rounding on this line changes the quadratic by O_M(A sqrt(S)), by the derivative at the real root and the exact quadratic remainder. The last term h(F(x)+F(y)) is O(S). As A=Theta(N^(1/5)) and S=Theta(N^(2/5)), both terms cost O_M(N^(2/5)). The entire choice can use exact integer/rational comparisons rather than floating square roots. This establishes the claimed absolute-and-congruence bound. To land exactly, one still needs a way to choose growing M or remove this actual residual without sacrificing the physical profile; constants from the fixed-modulus controller do not justify taking M larger than the error.

## Scope and replay

The same complementary controller also upgrades the original *balanced* reflected seam (h,P,L,R)=(5,256,142,142), signed ports (-2,0), and kappa=91. Its only exceptional odd prime for Fsum is 5, where R is a unit and w can be set to zero modulo the prime power. Thus its fixed designated **lower** member can be congruent to N within the same O_M(N^(2/5)) error, retaining the independent reflected address move. The rotor result instead realizes two consecutive +2 moves by retagging and does not need reflection.

A coarse quantitative route-length bound is L=O(M^3): the odd controller uses O(m^2 d) bits using repeated regulator words (D^m is the identity on the odd state); the higher dyadic corrections use O(2^s m d) bits. The line spacing is O(2^L m); the absolute rounding constant is O(2^(3L/2) sqrt(m)) when that spacing is small relative to the line height. Consequently M=o((log N)^(1/3)) can grow while giving error N^(2/5+o(1)). This does not put the error below the chosen modulus.

Run `python replay.py` for the literal center signatures and raw graph counts, and `python steer.py` for complementary suffixes and exact rational/integer scale checks. The latter does not enumerate the giant tournaments: its constructed words are checked by the exact readout theorem. Both scripts were replayed in normal and optimized Python, with identical outputs.
