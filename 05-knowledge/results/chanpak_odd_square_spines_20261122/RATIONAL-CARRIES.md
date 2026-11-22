# Simultaneous rational approximation of the two binary reads

**PROVED.** For any two prescribed suffixes of width \(\ell\), the chart equation reduces to an exact family of integer divisor problems. A common denominator chosen from the two suffixes guarantees at most \(1+O(Z/B^{1/3})\) possible carries when \(X+Y\le Z\), where \(B=2^\ell\) and \(X,Y\) are the high words. The constant is absolute. The [integer-slope inverse](../chanpak_farey_two_three_clocks_20261122/CARRY-FACTOR.md) gave small carries at nearly constant tails. The estimate here applies to every pair of tails, including aperiodic words. It bounds the work **per prescribed pair**. The number of pairs at width \(\ell\) is \(B^2\).

Write the physical lower count as
\[
 H=k+qxy+tx+cy+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y)).
\]
Use \(x=BX+i,\ y=BY+j\) with padded low words \(0\le i,j<B\). Put \(z_i=B-1-\operatorname{rev}_\ell(i)\) and \(g_i(X)=f(X)+[i<B/2][X\text{ odd}]\). The exact binary split is \(f(x)=f(i)+z_iX+g_i(X)\) with \(0\le g_i(X)\le X\), also at zero. Here \(q,t,c,h\) are positive, \(h\) is odd, and all other chart coefficients are nonnegative. [Binary split and absolute reader](../chanpak_decimal_kernel_20261122/REVERSED-PHASE.md)

Choose a positive denominator \(d\), and numerators \(0\le u,v\le d\). Set \(c_i=dz_i-Bu, c_j=dz_j-Bv\), and
\[
\begin{split}
 E&=c_iX+c_jY+d(g_i(X)+g_j(Y)),\\
 K&=H(i,j),\quad P=dqB,\\
 A&=d(t+qj)+hu,\quad C=d(c+qi)+hv .
\end{split}
\]
Expanding gives the integer identity
\[
 d(H-K)=B(PXY+AX+CY)+hE.
\]
Thus \(hE\equiv d(N-K)\pmod{B\gcd(P,A,C)}\) is necessary. The simpler modulus \(B\) already selects a single class of \(E\). Whenever \(X+Y\le Z\),
\[
 Z\min(0,c_i,c_j)\le E\le Z\max(0,c_i+d,c_j+d).
\]
For a surviving integer carry, put \(V=(d(N-K)-hE)/B\). Solve the divisor equation
\[
 (PX+C)(PY+A)=PV+AC .
\]
Keep the ordered nonnegative integer solutions with the required residues and check \(H(x,y)=N\). These operations are necessary and sufficient. The factors and the modulus record the absolute scale discarded by the rational phase.

## Why the common denominator exists

Take \(L=\lceil B^{1/3}\rceil\). The \(L^2+1\) points \((\{jz_i/B\},\{jz_j/B\})\), \(0\le j\le L^2\), occupy \(L^2\) half-open squares of side \(1/L\). Two share a square. Their positive index difference \(d\le L^2\) admits integers \(u,v\in[0,d]\) with \(|c_i|,|c_j|<B/L\). The displayed carry interval consequently has length at most \((2B/L+L^2)Z\), proving the claimed bound after division by \(B\). This is simultaneous Dirichlet approximation in two coordinates of the actual reader. With one free phase the same proof in dimension one gives a denominator and errors of order \(\sqrt B\).

The associated integer lattice has basis \((1,z_i,z_j),(0,B,0),(0,0,B)\) and determinant \(B^2\). Every vector directly supplies the candidate denominator and both residuals for the same physical pair. With LLL parameter \(3/4\) the first vector has norm at most \(\sqrt 2 B^{2/3}\). For \(B\ge4\) this is less than \(B\), so its first coordinate is nonzero. Negate if necessary to obtain \(d>0\). The other coordinates are signed residuals of magnitude less than \(B\), hence their numerators lie in \([0,d]\). This supplies the uniform carry bound in polynomial time in the bit width of the suffixes before the absolute search. The implementation checks the norm bound cubed and both integral reconstructions exactly. This uses the same small-vector operation that the [Leech projection replay](../../../04-computation/experiments/deep_holes_theta_20261022/audit_S53/audit_common.py) uses with a different lattice and norm. [Lenstra, Lenstra and Lovász, Proposition 1.6](https://ir.cwi.nl/pub/9304) For a particular suffix bank the replay below uses known fractions, even if their denominators differ, by clearing their least common multiple.

## The finite exponent-61 search

**FINITE-EXACT.** Set \(N=2^{60}(2^{61}-1)-1\). For each reduced \(u/d\in[0,1]\), \(d\le L_0\), form the nearest integer to \(Bu/d\), rounding ties upward. Allow \(z_i\) to differ from it by at most the tabulated radius, staying in \([0,B-1]\). Choose the two low words independently. Both high words range over **all** nonnegative integers, including zero on either coordinate, at every aspect ratio.

| Width | \(L_0\) | Radius | Distinct low-word pairs per spine | Positive-high carries, spine 4 / spine 8 | One-high-word reads, spine 4 / spine 8 | Hits on either |
|---:|---:|---:|---:|---:|---:|---:|
|34|11|0|1,764|17,604 / 258|0 / 0|0|
|37|21|2|490,000|70,210 / 1,108|520,038 / 132|0|

Neither spine had a positive-high divisor pair surviving the prescribed classes. For \(X,Y\ge 1\), the proof uses \(xy\le D=\lfloor(N-k)/q\rfloor\) and starts at \(x_0=B+i,y_0=B+j\). Convexity bounds \(x+y\) by the larger of \(x_0+\lfloor D/x_0\rfloor\) and \(y_0+\lfloor D/y_0\rfloor\). To enumerate divisors without a factorization package, the smaller of \(X,Y\) is at most \(\lfloor\sqrt D/B\rfloor\). Fixing it determines the other by a single exact division. When \(X=0\), the equation reads \(N-K=sY+h g_j(Y)\), with \(s=B(c+qi)+hz_j\). Its whole inverse lies in \([\lceil(N-K)/(s+h)\rceil,\lfloor(N-K)/s\rfloor]\), with \(sY\equiv N-K\pmod h\). The reflected face is handled the same way, counting the common corner once.

Unlike doubling on odd residue rings, multiplication by 2 modulo \(B\) has a kernel. The carry congruence above inverts the odd coefficient \(h\) modulo \(B\), and retains the discarded high words for the absolute divisor equation. It also works when \(d\) is even or shares primes with \(h\).

From this directory run the following. The first command compares complete small fibers to independent exhaustive loops, including deliberate nonoptimal denominators.

~~~sh
python3 -B rational_carries.py --check
python3 -B rational_carries.py --width 34 --limit 11
python3 -B rational_carries.py --width 37 --limit 21 --radius 2
~~~

## Repairing the near witnesses

**FINITE-EXACT.** For either [order-118 witness](NEAR-WITNESSES.md), concatenate its two low words of width 40. Every possible spine-four repair of the exponent-61 midpoint has Hamming distance at least **three** from that 80-bit vector. This holds even after replacing both high words by arbitrary nonnegative integers. The two fixed-width balls each contain \(1+80+\binom{80}{2}=3{,}241\) distinct suffix pairs. The LLL-selected common denominators reduce their absolute search to 180,213 carries at the negative-offset witness and 136,904 at the positive-offset witness. Neither ball contains a compatible integer divisor pair or a matching one-high-word face. The test concerns digit edits to an address, a separate operation from the distinguished single *tournament arc* whose reversal always changes the count by two.

Run the following with SymPy installed. The smaller example illustrates the selector on two new arbitrary suffixes. Each returned common denominator and both residuals are checked in integers. The two-bit repair test uses a fresh rational approximation for each changed pair.

~~~sh
python3 -B rational_carries.py --width 40 --suffix 512345678901 109876543210
python3 -B repair_suffixes.py --width 40 --radius 2
~~~
