# Complementary suffixes give a cube-root window for literal +2 pairs

2026-11-28. For every integer \(N\ge10^{12}\), the [fixed six-site Farey core](../chanpak_literal_square_shells_20261128) gives a deterministic pair of strong tournaments on \(O(\log N)\) vertices, with one differing arc and counts \(L,L+2\) satisfying
\[
N<L<L+2<N+160N^{1/3}.
\]
There is also a version that prescribes \(L\bmod29\) simultaneously, with window \(25000N^{1/3}\). As a consequence, every sufficiently large strict power shell with fixed exponent \(c>3/2\) admits a literal sibling pair. This includes \(c=5/3\), which the previous \(N^{2/5}\) window left open. This does not settle \(c=3/2\) or exact prescribed counts.

## 1. Complement the low bits before solving the size constraint

In the existing chart take even physical addresses \(x=2a,y=2b\). Its lower height is exactly
\[
L=111+1280ab+287a+353b+29\big(f(a)+f(b)+[a\text{ odd}]+[b\text{ odd}]\big),
\]
where \(f(0)=0,\ f(2u+d)=f(u)+(1-d)(u+[u\text{ odd}])\). The common apex and the flanks are unchanged across the arc reversal. The earlier paper supplies the labelled construction and exact profile. In particular it proves strong connectivity for positive even flanks.

Let \(D=2^k\), and in \(z=Du+v\) write exactly \(k\) bits for \(v\), including leading zeros. Put \(Y_v=D-1-\operatorname{rev}_k(v)\). The prefix recurrence gives
\[
f(z)= f(u)+uY_v+f(v)+[u\text{ odd}][v<D/2].
\]
Indeed the large prefix contributes \(u2^i\) at the \(i\)th zero, and it affects the parity term only for a zero in the first position. This is the suffix estimate used in [the square-phase calculation](../chanpak_square_one_apex_flip_density_20261128), now imposed as a construction constraint rather than averaged over. It involves weighted reversed bits, not just Thue–Morse parity.

Since \(f(u)+[u\text{ odd}]\le u\) and \(0\le f(v)\le v\), the difference \(f(z)-zY_v/D\) lies between \(-v\) and \(u+v\). Require
\[
m=a+b\equiv D-1\pmod D.
\]
The low remainders \(v,w\) of \(a,b\) then add to \(D-1\) and their padded strings are bitwise complements. Thus \(Y_v+Y_w=D-1\), and the sum of the two errors just bounded lies in \([-(D-1),\,m/D+D-1]\). Also
\[
(aY_v+bY_w)/D=m/2-m/(2D)+(a-b)(Y_v-Y_w)/(2D).
\]
Together these give the centered cancellation bound
\[
\left|f(a)+f(b)-m/2\right|\le m/(2D)+D+|a-b|/2.
\]
The factor \(1/2\) in the first term uses the two one-sided suffix estimates, not two independent triangle inequalities.

## 2. Decode in the correct sum class

Write \(q=N^{1/3}\) and take \(D=2^{\lfloor(\log_2 N)/6\rfloor-4}\), so \(\tau=D/\sqrt q\) belongs to \((1/32,1/16]\). The sum \(m\) is odd, hence the parity terms add to one. Replacing \(f(a)+f(b)\) by \(m/2\) gives the quadratic proxy
\[
Q=140+1280ab+(603/2)a+(735/2)b
 = 1280(a+p)(b+t)+E,
\quad p=735/2560,\quad t=603/2560,\quad E=140-1280pt.
\]
Set \(T=N+80\lceil q\rceil\), \(r=\sqrt{(T-E)/1280}\), and take the least permitted \(m\) with \(s=m+p+t\ge2r\). On \(a+b=m\), the smaller real solution of \(Q=T\) is \(a_0=s/2-p-\Delta\), where \(\Delta=\sqrt{s^2/4-r^2}\). Round \(a_0\) to the nearest integer and take \(b=m-a\).

At our threshold \(q\ge10^4\), \(0<E<140\), \(T<1.001N\), and \(\sqrt N/36<r<0.029\sqrt N\). Since \(s-2r<D\),
\[
\Delta\le\sqrt{rD}+D/2\le \sqrt{0.029\tau}\,q+D/2 .
\]
These estimates give \(a,b>0\) even after rounding (in particular \(r/q>\sqrt q/36\), while \(\Delta/q<0.044\)). Also \(m<2r+D\), \(|a-b|/2\le\Delta+1\), and \(|Q-T|\le1280(\Delta+1/4)\). Consequently
\[
|L-T|
 \le 1309\Delta+29r/D+29D+364
 \le \left(223\sqrt\tau+0.841/\tau\right)q+43\sqrt q+364 .
\]
The coefficient is convex as a function of \(\log\tau\), and is below 70 at both endpoints. The right side is less than \(71q\). This puts \(L,L+2\) strictly above \(N\) and below \(N+160q\). The mechanism balances \(O(\sqrt{rD})\) rounding against \(O(r/D)\) leftover binary dependence; \(D\) of order \(N^{1/6}\) makes both errors of cube-root size.

### Exact integer comparisons

With \(A=2560,\ B=603,\ C=735,\ M=A(2T-280)+BC\), compute
\[
m_0=\left\lceil\frac{\lceil\sqrt{4M}\rceil-B-C}{A}\right\rceil,\qquad
m=m_0+((D-1-m_0)\bmod D),\qquad S=Am+B+C .
\]
The residue is nonnegative. Compute \(a=\lfloor(S-2C+A-\lceil\sqrt{S^2-4M}\rceil)/(2A)\rfloor\) and \(b=m-a\). Replacing a square root by its integer ceiling before the outer ceiling or the indicated subtracted floor does not affect the result, since all remaining numerators and the divisors are integers. This takes only integer operations. Both addresses are \(O(\sqrt N)\), so the word lengths are logarithmic. At \(N=10^{12}\) the decoder gives \((x,y)=(55454,56352)\) and \(L=1{,}000{,}000{,}755{,}937\).

## 3. A fixed congruence and an absolute window at once

Modulo 29 the nonlinear term vanishes. To prescribe any residue \(R\), enforce \(a\equiv0\bmod29\) and \(111+353m\equiv R\bmod29\); \(353\equiv5\bmod29\) is a unit. Combine the latter condition with \(m\equiv D-1\bmod D\) by CRT. With \(T=N+12000\lceil q\rceil\), take the first allowed \(m\) above \(2r-p-t\) and round the same smaller root to the nearest multiple of 29. The lattice spacings are \(29D\) for \(m\) and 29 for \(a\), so \(m<2r+29D\), \(s-2r<29D\) and the rounding error in \(a\) is at most \(29/2\). Now \(r<0.04\sqrt N\) and \(\Delta\le\sqrt{29rD}+29D/2<0.30q\). The lower bound \(r>\sqrt N/36\) still ensures positivity. The proxy error is at most \(1280(29\Delta+29^2/4)<11164q\). For the binary error, use \(m/(2D)<1.28q+14.5\) and \(|a-b|/2<\Delta+15\); multiplying the centered suffix bound by 29 gives less than \(47q\). Thus the total is \(<12000q\). Centering and allowing for the ceiling and the extra two gives the stated \(25000q\) window. This argument addresses a fixed useful modulus, not arbitrary moduli or an exact inverse.

## 4. Strict power shells and the \(5/3\) witness

For fixed real \(c>3/2\), take \(N=\lfloor h^c\rfloor+1\). The power-shell width is at least \(c h^{c-1}\), and \(N^{1/3}=O(h^{c/3})\). Since \(c-1>c/3\), eventually both child intervals and their upper endpoints are strictly inside the parent's interval. As in the [literal square tree](../chanpak_literal_square_shells_20261128), the rooted intervals strictly nest and shrink, so their common point gives the actual floors. Distinct binary choices give distinct real seeds; no primeness is implied.

For \(c=5/3\) the comparison works for every \(h\ge2^{63}\): \(h^{1/9}\ge128\) and \(160N^{1/3}+2<202h^{5/9}+2<(5/3)h^{2/3}\). Start at the strong one-apex count \(2^{63}+1\), from a 1 followed by 63 zeros. The first two children, on 105 vertices, are
\[
40564819207303340857117708104011,\qquad
40564819207303340857117708104013 .
\]
The replay verifies their strict inclusion by cubing the lower count and the upper endpoint and comparing with the corresponding fifth powers.

Run node replay.js. It validates 2,040 suffix identities, 17,505 exact cube-root-window cases, 116 simultaneous residue checks, and the first \(5/3\) shell. These tests use integer inequalities only and are not used as evidence for the infinite bounds.
