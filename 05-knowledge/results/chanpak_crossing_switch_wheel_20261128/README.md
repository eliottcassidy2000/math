# Binary crossing switches compress the literal +2 cube-root wheel

2026-11-28. The [arbitrary-modulus window](../chanpak_cube_root_wheel_20261128) used up to \(M/\gcd(M,116)-1\) repeated blocks. The replacement below uses **\(O(\log M)\) blocks**, each of length \(O(d+\log M)\), where \(d\) is a suitable even multiple of the binary order modulo the odd part of \(M\). The fixed Farey core, one-arc property, and \(N+C_MN^{1/3}\) theorem are unchanged. No primeness or exact target claim is made.

For \(M=792=2^3\cdot3^2\cdot11\), the padded suffix drops from 14,585 bits in the prior replay to **910 bits**, and the sufficient target length drops from 87,511 to **5,461 bits**. These are proof thresholds, not the first target sizes where the construction works. For \(29^2=841\), the new suffix has 12,292 bits and the sufficient target length is 73,753 bits.

## 1. Move one bit across another one

Write \(M=2^e m\) with \(m\) odd. Choose an even multiple \(d\) of \(\operatorname{ord}_m(2)\), with order 1 when \(m=1\). In a block of \(d+3\) bits put a movable one at position \(1\) or \(1+d\) and a fixed one at \(1+\ell\), numbered from zero. Pad by \(e+2\) zeros on both ends. The displacement preserves forward and reverse address contributions modulo \(M\): its powers differ by a period modulo \(m\), and the padding kills both endpoint powers modulo \(2^e\). By the identical argument it preserves the moved bit's incident sum and all its interactions with the other blocks. In the identity for \(f(v)+f(\bar v)\), sum incident weights at one-sites and subtract twice the one-one pair weights. When \(2\le\ell\le d-2\), the internal distance changes from \(\ell\) to \(d-\ell\), so the change in the tournament count is
\[
29C_\ell,\qquad C_\ell=2^\ell-2^{d-\ell}.
\]
This statement is modulo \(M\) and holds irrespective of the other settings. If \(\ell=1\), the adjacency bonus makes the step \(29(4-2^{d-1})\). It is a multiple of four and is 1 modulo 3. For \(\ell\ge2\) the crossing is zero modulo 3 whenever \(3\mid m\), since powers of two are self-inverse mod 3. For \(2\le\ell<e\) and \(d-\ell\ge e\), its change mod \(2^e\) is \(29\cdot2^\ell\).

All blocks are placed in the complementary low word, with two leading zeros and one of the two five-bit tails from the [prior wheel](../chanpak_cube_root_wheel_20261128). The construction of the mod-29 sum and starting-address lattice is unchanged, using the actual baseline word (including any reversed initial blocks). A block substitution keeps both addresses modulo \(M\); the bit-reversal identity from that proof therefore sends the stated change in \(f(v)+f(\bar v)\) to 29 times that change in the true height, including for unequal high prefixes. The same tail still controls mod 4.

## 2. Three crossings write one odd-modulus binary digit

Put \(g=\gcd(m,87),\ z=m/g,\ n=\lceil\log_2 z\rceil\), taking \(n=0\) for \(z=1\). Let \(d\) be the first positive multiple of \(\operatorname{lcm}(2,\operatorname{ord}_m(2))\) at least \(2e+n+6\). First set mod 29 and mod 4 as above. If \(3\mid m\), install two adjacent crossings and use zero, one, or two to match the residue mod 3. Then use the crossings at \(\ell=2,\ldots,e-1\) in increasing order to match each remaining dyadic bit: at step \(\ell\) the increment has valuation \(\ell\) modulo \(2^e\), and so does not unset a lower bit. These dyadic switches preserve the residue mod 3, and every step so far preserves mod 29 when needed. Oddness gives the initial bit when \(e\ge1\).

The remaining groups use \(\ell=e+2+i\) for \(i=0,\ldots,n-1\). Each group consists of two initially unmoved crossings at \(\ell+1\), and one initially moved crossing at \(\ell\). Switching all three gives
\[
29(2C_{\ell+1}-C_\ell)\equiv87\cdot2^\ell\pmod m
\]
since the inverse-power terms cancel. All their participating distances are at least \(e\), so this preserves mod \(2^e\). If \(L_0\) is the current lower count, it already equals \(R\bmod g\). Solve
\[
X\equiv \frac{R-L_0}{g}
 \left(\frac{87\cdot2^{e+2}}{g}\right)^{-1}\pmod z.
\]
The inverse exists because \(87=3\cdot29\) is squarefree and \(g=\gcd(m,87)\). The bits of \(0\le X<z\) choose the groups. Thus the final count equals \(R\bmod M\).

The number of reserved suffix bits is exactly
\[
K=7+\left(2[3\mid m]+\max(e-2,0)+3n\right)\left(d+3+2(e+2)\right).
\]
All choices stay in the same cell of width \(P=2^K\). Consequently the old absolute placement proof applies without modification: for \(N\ge P^6\), center at \(N+100000P\lceil N^{1/3}\rceil\), use the same square-root decoder, and obtain \(N<L<L+2<N+200000PN^{1/3}\). Run node replay.js to reproduce 285 exact decodes, including 792, 841, and the simultaneous avoidance of primes through 11 with \(M=2310\).

## 3. Growing moduli and a larger small-prime cutoff

Already \(d=O(M+\log M)\) gives \(K=O(M\log(2M))\). In particular the near-cube-root window \(N^{1/3+o(1)}\) and the threshold \(N\ge P^6\) hold uniformly for \(M\le\log N/(\log\log N)^2\). The logarithms in the rest of this section are natural. A primorial wheel does substantially better because its binary period depends on the shifted-prime lcm
\[
H(B)=\operatorname{lcm}\{p-1:p\le B,\ p \text{ prime}\}.
\]
Indeed the order modulo the odd primorial divides \(H(B)\), and the elementary primorial bound gives \(K=O(B(H(B)+B))\).

We use [Kevin Ford, *Integers divisible by shifted primes*, Lemma 4.4](https://www.ford126.web.illinois.edu/wwwpapers/divshift.pdf) to bound \(\log H(B)=O(B\log\log B/\log B)\). Set \(y=B/\log B\). Prime bases at most \(y\) contribute \(O(y)\) for their first powers by the elementary Chebyshev bound, with at most \(O(\sqrt B\log B)\) extra for higher powers (each higher-power base is \(\le\sqrt B\) and its entire contribution is \(\le\log B\)). A prime \(q>y\) dividing \(H(B)\) has \(p=kq+1\le B\) prime with \(k\le\log B\). For large \(B\), \(k\) must be even. The first estimate of Ford's lemma with \(x=B/k,\ z=\sqrt x\) bounds these prime \(q\) by \(O(B/(\phi(k)\log^2 B))\): \(q>y>\sqrt x\) is among the sifted integers, and \(\log x\sim\log B\) uniformly. Multiply by \(\log B\) and sum over \(k\). The estimate \(\sum_{k\le R}1/\phi(k)=O(\log(2R))\) follows by expanding \(k/\phi(k)=\sum_{d\mid k}\mu^2(d)/\phi(d)\) and using convergence of \(\sum_d\mu^2(d)/(d\phi(d))=\prod_p(1+1/(p(p-1)))\). This proves the claimed lcm bound.

It follows that \(\log K=O(B\log\log B/\log B)\). With \(t=\log\log N\), take \(B=\lfloor\eta t\log t/\log\log t\rfloor\) for a sufficiently small positive absolute \(\eta\). Then \(K\le(\log N)^{1/2}\) eventually, making \(N\ge P^6\) and the window \(N^{1/3+o(1)}\). Prescribing \(R=-1\) modulo the primorial makes the two counts \(-1,+1\) on the entire wheel. Both therefore avoid every prime up to \(B\), a cutoff beyond any fixed multiple of \(\log\log N\). The sieve is only used to shorten the primorial period; it does not assert that the resulting counts are prime.
