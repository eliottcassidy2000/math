# Every compatible joint boundary profile modulo every odd modulus

For a tournament U write H=Q_U(1) and q₂=Q_U(2), where Q is the optional-slot/ordered-path-cover polynomial (a slot may be empty). For every odd M≥3 and every prescribed parity b of the tournament order, a pair (A,B) occurs as (H,q₂) modulo M in a **strong tournament** iff, when 3 divides M,
$$
B\equiv (-1)^b A\pmod3.
$$
There is no restriction when 3 is a unit. This includes zero and nonunit coordinates, both parities, and mixed prime powers. It does not give absolute targets, nor a general even-modulus joint characterization. The p≥5 **ordinal** packet subroutine is from Session 326's joint residue theorem; the 3-adic exceptional derivative, common scaling wrapper, and neutralizers below allow everything to coexist physically.

## 1. Coordinates and the non-generic derivative at three

Q has only powers of the same parity as the order N, by the odd-cycle packing formula. Hence Q(2)=Q(−1)=(−1)^N H modulo 3. Put (h=H,\ q=2^{-N}Q(2),\ D=(q-h)/3) in 3-adic arithmetic. Products of ordinal components act by
$$
(h,D)(u,E)=(hu, hE+uD+3DE).
$$
Singletons are identity in these normalized coordinates.

Use the inherited one-apex patterns (P_d(n)) with transitive chain of length dn, the apex beats the positions starting each length-d block, and d=2,3. Define
$$
\Phi_n(r)=\sum_{0\le k<n}(n-k)r^k
 =\frac n{1-r}-\frac{r(1-r^n)}{(1-r)^2}.
$$
The two profiles in normalized coordinates follow by summing the path-gap correlation by residue class within the d-block:
$$
h=1+n+(2^{d-1}-1)\Phi_n(2^d),\qquad
q=1+\tfrac12 A_d\Phi_n((3/2)^d)+\tfrac12 B_d\Phi_n(2^{-d}),
$$
with ((A_2,B_2)=(1/2,1/2),(A_3,B_3)=(5/4,3/4)). Reproduction evaluates these as a division-free polynomial and does not invert 3.

For shifts (Delta=3^t k), t≥1, hold n parity fixed when the ratio is a unit −1 modulo 3. If (r\equiv−1\pmod3), the normalized difference (Delta\Phi/\Delta) modulo 3 is (1/(1-r)) in the sense (Delta\Phi\equiv \Delta/(1-r)\pmod{3^{t+1}}), since (r^\Delta=1\pmod{3^{t+1}}). If r is a nonunit, the same assertion holds whenever n is large enough that (r^n=0\pmod{3^{t+1}}). The constructions can start with all variable n≥2e+20 for precision 3^e.

For (r=1+s, 3\mid s) the derivative is instead
$$
\Delta\Phi\equiv \Delta(n+1/2+s/3)\pmod{3^{t+1}}.\tag{1}
$$
Indeed (Phi_n(1+s)=\sum_{j\ge0}\binom{n+1}{j+2}s^j). Vandermonde with (v_3\binom{\Delta}{i}\ge t-v_3(i)) leaves the usual j=0 terms i=1,2 and the **extra term** j=1,i=3: (s\binom{\Delta}{3}\equiv s\Delta/3\pmod{3^{t+1}}). Every j≥2 vanishes because (j-v_3(i)\ge1) for i≤j+2. These increments use nonnegative shifts only; if k is divisible by 3 the congruence also follows by raising t.

For P₂, the three ratios are 4, 9/4 and 1/4. Thus h has slope n+1, q has slope n+2, modulo 3; the split of slopes is due to the s/3 term in (1). Therefore shifting n by (3^{t+1}k) changes D by (3^t k) modulo (3^{t+1}) while fixing h at that precision. For P₃, the h slope is 4/7 and the q slope 22/133, both 1 modulo 3, holding n parity. Thus a shift (3^t k) gives h increment (3^t k) and D increment divisible by (3^t). At first digits, direct evaluation of the same formulas gives:
* for P₂ with even n divisible by 3, h=1 mod 3 and D takes 0,1,2 at n=0,3,6 mod 9, respectively;
* for P₃, (h=1+n\pmod3), and D=0 for even n and D=2 for odd n modulo 3. The latter holds for positive odd n (the nonunit ratio powers already vanish at the required precision). For any fixed n mod 3 and parity, h mod 9 ranges through all three lifts under n→n+6, by the h slope;
* P₃(8) has h=3 mod 9 and D=0 mod 3. The congruences h=3 mod 9 and D=0 mod 3 persist at large n=8 mod 18. At higher precisions work from a large representative, since the nonunit exponential can survive at n=8.

## 2. Ordinal synthesis at 3

Prescribe h,D modulo (3^e). If h is a unit, take (P_2(n)\to P_3(m)) with n even and 0 mod 3, m even. Match the first h digit by m and then the first D digit by n mod 9. At digit t, t≥1, vary n by (6\cdot 3^ti) and m by (2\cdot3^t j) with i,j=0,1,2. The h derivative in i is zero, and its derivative in j is a unit. The D derivative in i is nonzero because the P₃ h is a unit. Thus h chooses j and D chooses i. No derivative formula for D in m is needed.

If h is a nonunit and D a unit, take (P_3(n)\to P_3(m)) with n odd and n=2 mod 3, m even with h unit. Their normalized total D mod 3 is (2h(P_3(m))), so choose m mod 3 to match it and n mod 9 to match h/3 mod 3. At digit t, vary n and m as before but use (h/3,D). The h/3 derivative in i is (2h(P_3(m))≠0\pmod3); the D derivative in i is zero, and its derivative in j is (2D(P_3(n))\ne0\pmod3). Consequently choose j first, then i. Since n is already fixed modulo 3, h/3 has no divisibility problem. The h/3 comparison may be made one digit beyond the final h requirement.

If h and D are both nonunits use a scaling letter (F=P_3(n_0)) with even (n_0=8\pmod{18}), large enough for the current precision. Thus (h_F=3u) and (D_F=3v) with u and (u+3v) units. To realize (h,D), recurse at precision (3^{e-1}) for
$$
h_0=(h/3)/u,\qquad D_0=[D/3-vh_0]/(u+3v),
$$
and ordinally prepend F. Depth e ends at the empty (algebraic) identity. This proves the 3-local ordinal result. Padding by singleton SCCs allows (2^N) to be either sign modulo (3^e), compatible with either desired parity. In view of the nonunit exponential, **all** 3-packet letters used in the mixed synthesis, including F, must start at large n at the original precision, not the smaller recursion precision.

## 3. Simultaneous integer lengths at every odd prime

Write (3^e\Vert M) and process 3 first (omit when e=0). Make the initial 3-packet have the desired source parity, and pad by at least e singleton SCCs so that (Q_U(3)=0\pmod{3^e}). The raw required q₂ divided by the chosen sign (2^N) is precisely a normalized q compatible with h. The following details ensure later primes are not inadvertently made singular and later packets keep the 3-profile.

**3-packet unit avoidance.** Each length can be shifted by multiples of (2\cdot3^{e+2}) without disturbing either normalized 3-coordinate or the parity. At each later prime p in increasing order, freeze its phase modulo p−1 and choose n mod p to preserve both units. The p≥5 nonzero-polynomial claims for d=2,3 are those in the inherited packet proof: h and q each have degree at most two in n on a phase, with nonzero slope unless there is a nonconstant quadratic or the q₂ slope exception d=3,p=11. At 11, directly evaluate q normalized by (2^N) for phases n=0,...,9 to get
$$
1,2,0,1,5,5,4,6,5,1.
$$
The sole bad phase n=2 mod 10 is avoided by forbidding n=2 mod 5 at the digit-5 step, performing that preparatory step even if 5 is not requested. For d=3 both h and q are affine nonconstant at 5 on a fixed phase, so this extra forbidden digit leaves a choice. Later changes preserve all earlier residues; p−1 has no factor p, and no previous prime step is divisible by p. Only after adjusting all lengths choose the singleton padding/both parity and (2^N\equiv\pm1\) at 3.

**Later packets and neutralizing 3.** At p≥5 use the inherited S326 ordinal packet, dividing by all preceding unit values at the current p. Require every new pattern length to be divisible by (2\cdot3^{e+2}) when e>0, in addition to the usual retained periods for processed primes. This is coprime to p and so does not restrict its local interpolation. At 3, each such positive-length d-letter has fixed normalized pair ((1,\kappa_d)), where
$$
\kappa_d=1+\tfrac12 A_d\frac{r_d}{(1-r_d)^2},\qquad r_d=(3/2)^d,
$$
reduced modulo (3^e). In the nonunit term (r_d^n=0\pmod{3^e}); all linear-n terms vanish at 3 and the other unit terms vanish by the retained progression (or repeated increment from n=0). Thus (kappa_d\equiv1\pmod3) is always a unit. To cancel this constant, add enough neutral d-letters that the total multiplicity of each d is divisible by the order of (kappa_d) modulo (3^e). The neutralizing lengths are positive multiples of all processed (r^{f+2}(r-1)), **including the current prime** r=p with f its requested exponent. At p≥5 these new letters read (H,q₂)=(1,2), so the normalized packet product there is unchanged. At future primes select good digits by the same avoidance; d=3 is now used only at p=5 and a neutralizing d=3 letter has n=0 mod 5. Singleton padding makes total letter-plus-padding count a multiple of all processed Euler periods, including the 3-Euler period. Since the additional dn lengths are 0 modulo the 3-Euler period, the raw Q₂ effect of this packet at 3 is also identity. Each later packet has even order (pattern lengths even, count with padding even). Therefore the source parity is retained. When e=0 an initial zero or one singleton sets its parity, and no 3 neutralizers are needed. The initial ≥e singletons make (Q_U(3)) vanish through all later packets at 3.

## 4. One wrapper with scalar actions at every prime

The literal strong gluing is (C_3[U,TT_m,1]); it adds m+1 vertices. For every m≥1 its full polynomial is
$$
\begin{split}
Q_\mathrm{out}(z)={}&z^{m+1}Q(z)\\
 &+z((z+1)^m-z^m)(Q(z+1)-Q(z))\\
 &+z((z-1)^m-z^m)(Q(z-1)-Q(z)).\tag{2}
\end{split}
$$
To prove (2), use the odd-cycle packing identity. Every cross-module cycle must use the singleton, so there is at most one such cycle. It contains a nonempty directed path on a subset S of U, a nonempty increasing path on a subset of TT_m, and then the singleton. The cross-cycle is odd iff the sizes of those two subsets have the same parity. The odd-cycle term contributes twice a factor z, while unused chain vertices contribute their singleton factors z. Averaging the two parity signs gives exactly the two terms. To justify the signed subset sum over U, the set convolution on disjoint vertex subsets satisfies (Q_{U}(z+w)=\sum_S Q_{U[S]}(w)Q_{U-S}(z)) by recording which vertices go to the first w slots, and then by polynomial continuation. At w=1 this counts a marked path; at w=−1 the parity of each induced polynomial gives the signed version. Remove the empty S to subtract Q(z). When no cross-cycle exists the contribution is (z^{m+1}Q(z)). This proves (2), including all polynomial evaluations near z=0.

Choose m positive, at least e and divisible by all Euler periods (p^{f-1}(p-1)) for prime powers (p^f\Vert M). It is even. At z=1, (2) gives (H'=(3-2^m)h+(2^m-1)q_2\equiv2h\pmod M). At z=2: for p≥5 all brackets vanish and (Q_2'\equiv2q_2); at 3, (2^m=1,3^m=0,Q_U(3)=0\pmod{3^e}), so (Q'_2\equiv4q_2). All scale factors are units. Prescribe the source profile by dividing by these scale factors, and its order parity as 1−b. At 3 the output compatibility transforms exactly to source compatibility because 4/2=−1 modulo 3. This completes the theorem.

## Scope

This is an exact local-global statement **only for congruences**: the full odd-modulus response space of an order-parity-conditioned strong source has one residue constraint, already at 3. It removes modular impossibility as an explanation for failure of these two-coordinate interfaces. It is not a positive-integer selection theorem. The constructions use large lengths for CRT and repetition, a cycle wrapper with an m chosen by the modulus, and no absolute height control sufficient for Chan–Pak. What remains is an Archimedean/physical budget mechanism, not a further test modulo an odd integer for this pair alone.

Two consequences separate modular reachability from structural capacity.

* **Common scaling on every unit slot coordinate.** Under the same modulus-dependent m and at least max f singleton padding, (2) gives, whenever z is a unit mod an odd prime power p^f,
  $$
  Q_\mathrm{out}(z)\equiv
  \begin{cases}
  2zQ_U(z)&z\equiv\pm1\pmod p,\\
  zQ_U(z)&z\not\equiv\pm1\pmod p.
  \end{cases}
  $$
  Precisely one of z±1 is a nonunit in the first case, and its source profile is killed by padding; in the second both brackets vanish. At a nonunit z the wrapper instead reads (z[Q_U(z+1)+Q_U(z-1)]\), so extension to simultaneous multiple slots must preserve this unit/nonunit boundary.

* **The literal bank is sparse in absolute height.** All the constructed sources are ordinal joins of P₂, P₃ and singleton letters, with an outer (C_3[U,TT_m,1]). The output coordinates depend only on the inner component multiset and m, by (2) and ordinal multiplicativity. Since (Q_\mathrm{out}(2)\ge 2^{|V_\mathrm{out}|}), for height at most X the order budget is O(log X). There are at most two component types of each size, so the number of possible unwrapped multisets up to that order is exp(O(sqrt(log X))) by the elementary colored partition bound. Multiplying by the O(log X) choices of m does not change the exponent. Thus even surjectivity onto every permitted residue class at every odd modulus does not make this bank an absolute joint selector.
