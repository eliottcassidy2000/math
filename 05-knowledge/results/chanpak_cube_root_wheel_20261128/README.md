# An arbitrary finite wheel inside the literal +2 cube-root window

2026-11-28. **Theorem.** Fix any integer \(M\ge1\) and a residue \(R\bmod M\) attained by an odd integer. For every sufficiently large integer \(N\) the same [six-site Farey core](../chanpak_literal_square_shells_20261128) admits a strong labelled pair differing by one arc, with counts
\[
N<L<L+2<N+C_MN^{1/3},\qquad L\equiv R\pmod M.
\]
The hosts have logarithmic order. This extends the [mod-29 window](../chanpak_complementary_suffix_window_20261128) to arbitrary fixed mixed prime powers, in particular \(2^a3^b11^c\). Choosing the wheel to exclude prime factors simultaneously is a consequence, not a prime certificate. The construction has a conservative threshold and does not choose an exact prescribed integer.

## 1. A switch invisible to both addresses

Use physical flank addresses \(2a,2b\) in the published core. The exact lower count is
\[
L(a,b)=111+1280ab+287a+353b+29(f(a)+f(b)+[a\text{ odd}]+[b\text{ odd}]),
\]
with \(f(0)=0,\ f(2u+d)=f(u)+(1-d)(u+[u\text{ odd}])\). Positive addresses retain strong connectivity and the single arc in the core.

For a padded \(K\)-bit word \(v\) and its bitwise complement \(\bar v\), put \(g(v)=f(v)+f(\bar v)\). The inversion law gives
\[
g(v)=\sum_{i<j,\, v_i\ne v_j}(2^{j-i-1}+[j=i+1]).
\]
Equivalently, sum the incident weights of the one-sites and subtract twice the weights between pairs of ones. At position \(i\) (numbering from zero), the incident sum is \((2^i-1)+(2^{K-1-i}-1)\) plus its number of immediate neighbors.

Write \(M=2^e m_{\mathrm{odd}}\). Choose any even \(d\ge2\) with \(2^d\equiv1\bmod m_{\mathrm{odd}}\), for instance the least positive even such exponent (take \(d=2\) if the odd part is 1). A block of \(d+4\) bits has ones at positions \(1,2\) in its first setting and \(1,2+d\) in its second. Pad it by \(e+2\) zero bits on both ends. The two settings agree modulo \(M\) in both forward and reversed binary value: the displaced bit has the same power mod \(m_{\mathrm{odd}}\) and both changed powers are zero mod \(2^e\). In any larger word the two moved-bit incident sums also agree mod \(M\), and every weight from the moved bit to an outside one-site agrees mod \(M\), by the same padding and displacement. The only internal pair weight changes from \(2\) to \(2^d\). Hence one such switch changes \(g(v)\) by \(4-2^{d+1}\) modulo \(M\). The corresponding height step is
\[
I=29(4-2^{d+1}),\qquad G=\gcd(M,I)=\gcd(M,116),
\]
since the parenthesis has valuation two at 2 and is 2 modulo the odd part.

## 2. Realize the correct class modulo \(G\), then step within it

Let \(J=M/G-1\). Put \(J\) padded blocks in the low word \(v\), prepend two zero bits, and append the five-bit tail \(00000\) or \(11101\). Let the resulting width be \(K=7+J(d+4+2(e+2))\), and \(P=2^K\). Use \(a=Pu+v,\ b=Pw+\bar v\), with the same higher addresses \(u,w\) across all choices. This is possible whenever the sum \(s_0=a+b\) is \(-1\bmod D\) for a dyadic \(D\ge P\). The [exact suffix identity](../chanpak_complementary_suffix_window_20261128) shows that the part of \(f(a)+f(b)\) depending on \(v\) is
\[
(u-w)(P-1-\operatorname{rev}_K(v))+g(v)
 +[u\text{ odd}][v_0=0]+[w\text{ odd}][v_0=1].
\]
Thus every internal block switch changes \(L\) by \(I\bmod M\), including when the other tail is used.

If \(29\mid M\), choose \(s_0\) with \(111+353s_0\equiv R\bmod29\), and choose \(u\) so \(Pu+v_{\rm base}\equiv0\bmod29\) for the all-adjacent word with tail \(00000\). The other tail adds 29 to \(a\), so both tails give the same residue mod 29. Internal switches preserve it as well. This works because the nonlinear term vanishes mod 29, and with \(a=0\bmod29\) the polynomial reduces to \(111+353s_0\), with unit slope.

For mod 4 only distances 1 and 2 in \(g(v)\) matter, each with weight 2. The parities of their unequal-pair counts are the XORs of the boundary bits of the respective distance graphs. Adding them cancels the first and last bits, so \(g(v)\equiv2(v_1\mathbin{\rm xor}v_{K-2})\bmod4\). Both tails have penultimate bit zero and the leading two bits stay fixed, including the reverse value mod 4. Hence \(f(a)+f(b)\) is unchanged mod 4. The sum \(s_0\) is odd, so the parity terms add to one. Meanwhile \(287a+353b\equiv s_0+2a\bmod4\) changes by 2 on replacing the tail. Select the tail as needed when \(4\mid M\); oddness handles the smaller two-part.

The resulting baseline is \(R\bmod G\). Switch the first \(j\) blocks, \(0\le j\le J\), to move by \(jI\). This exhausts that class modulo \(M\), and the correct \(j\) can be obtained by inverting \(I/G\bmod M/G\).

## 3. Absolute placement with exact comparisons

Take \(N\ge P^6,\ q=N^{1/3},\ D=2^{\lfloor(\log_2N)/6\rfloor}\) and \(T=N+100000P\lceil q\rceil\). Then \(q\ge P^2\), \(P\ge128\), \(\sqrt q/2<D\le\sqrt q\), and \(D\ge P\). Use the previous proxy
\[
Q=1280(a+p)(b+t)+E,\quad p=735/2560,\quad t=603/2560,\quad E=140-1280pt.
\]
Write \(r=\sqrt{(T-E)/1280}\). Choose the first \(s_0\) in the required CRT class with \(s=s_0+p+t\ge2r\), so \(s-2r<29D\). On this line the smaller real solution is \(a_0=s/2-p-\Delta\) with \(\Delta=\sqrt{s^2/4-r^2}\). Choose \(Pu+v_{\rm base}\) nearest \(a_0\) on its lattice of spacing at most \(29P\), and use the tail and the \(j\) switches within that cell. The total movement from \(a_0\) is less than \(16P\). Keep \(b=s_0-a\); its low word is the actual complement of \(v\) and its higher address \(w\) does not change within the cell.

Indeed \(\sqrt N/36<r<\sqrt N/25\) from \(0<E<140\) and \(T<2N\), and \(\Delta\le\sqrt{29rD}+29D/2<2q\). As \(r>(P/36)q\) and \(16P\le q/8\), both adjusted addresses are positive. The proxy changes by at most \(1280(2\Delta\cdot16P+(16P)^2)<86000Pq\). Since complementarity holds modulo all of \(D\), the published centered cancellation estimate applies on those low bits. Its error in \(L\) is at most \(29(s_0/(2D)+D+\Delta+16P+1)<100q\). Thus \(|L-T|<90000Pq\). In particular \(C_M=200000P\) suffices for the stated one-sided window. Addresses are \(O(\sqrt N)\) for the fixed wheel.

The helper uses precisely the previous integer-ceiling square roots for the minimum sum and the nearest-root step, then selects a lattice cell by integer division and solves the one-variable congruence. No floating-point calculation certifies a case. Run node replay.js: it checks the switch identity for \(M=1,\ldots,100\) on varied prefixes, both tails and multiple substitutions, and runs 44 full integer decodes, including \(792=2^3\cdot3^2\cdot11\) and \(841=29^2\).

## 4. An unbounded wheel still fits a near-cube-root window

For a fixed finite set of primes, taking \(R=-1\) modulo their product gives \(L\equiv-1,\ L+2\equiv1\) and avoids all of them. One can even let the wheel grow with \(N\). In the above construction \(d\le2M\), \(J<M\) (or \(J=0\)), and \(e\le\log_2 M\), so \(K=O(M^2)\) with an absolute constant. For \(X=\log_2 N\), all moduli \(M\le X^{1/4}\) have \(K=O(\sqrt X)\), satisfy \(N\ge P^6\) for all sufficiently large \(N\), and give a uniform window \(N^{1/3}2^{O(\sqrt X)}=N^{1/3+o(1)}\). In particular one can use the product of primes up to
\[
B=\left\lfloor\frac{\log_2 X}{8}\right\rfloor .
\]
The prime product through \(n\) is at most \(4^n\) by induction: for \(n=2k\), all primes \(k<p\le2k\) divide \(\binom{2k}k\le4^k\); for \(n=2k+1\), use \(k+1<p\le2k+1\) and \(\binom{2k+1}k\le4^k\), since the two equal middle coefficients together are at most \(2^{2k+1}\). Thus the wheel product is at most \(4^B\le X^{1/4}\). Therefore both counts can avoid every prime at most \(B\) in a single physical pair. This uses no prime-distribution estimate and does not imply either count is prime. The near-cube-root power still gives strict Mills-style branching at every exponent \(c>3/2\), even with this growing wheel.
