# A linear-length nonlinear bypass

The [separation theorem](SEPARATION.md) is about triangular **one-control-per-digit** decoding,
not about the existence of more strongly correlated clocks. There is an
explicit \(O(k)\)-bit guard with a nonlinear deterministic decoder,
which can replace the \(O(k^2)\) reset pulses in the prior bracket.

Work in the complementary pre-address rail of the
[both-one clock](../chanpak_collatz_norm_clock_20261115/DYADIC_CLOCK.md).
Let \(U,V\) be the pre-addresses, \(S=U+V\) and
\(K=p(U)+p(V)-U+[U\text{ odd}]\). Append complementary bits \(e,1-e\).
On the mixed-parity sheet the exact increment is
\[
 \Delta K=e(S+1)-2e(U+[U\text{ odd}]).
\]
Assume \(S=-1\bmod 2^k\), reset \(U=0\bmod 2^k\) with \(k\)
zero bits on its address, and then read a finite binary mask from high
to low. If its bits from low to high are \(b_i\), the total increment
modulo \(2^k\) is \(-4F(b)\), where \(F\) is the integer quadratic
polynomial
\[
 F(b)=\sum_{i<j}c(j-i)b_i b_j,\qquad
 c(1)=1,\quad c(d)=2^{d-2}\ (d\ge2).
\]
Indeed each current 1 meets a previous 1 at distance \(d\), and its
contribution before dividing by \(-4\) is \(2^d\) for \(d\ge2\)
and \(4\) for neighbors (the parity guard doubles the neighbor term).

## Each packet has a complete Fourier pairing

Put \(n=k-2\), \(q=2^n\) and suppose \(k\ge3\). A packet consists
of two freely chosen \(n\)-bit words separated by exactly one zero.
If \(u_i\), \(v_j\) are bits of the low and high words
(\(0\le i,j<n\)), the cross term has coefficient
\(2^{n-1-i+j}\). Thus the packet value is exactly
\[
  F_{0}(u)+F_{0}(v)+X Y,\quad
  X=\sum_{i=0}^{n-1}2^{n-1-i}u_i,\quad Y=\sum_{j=0}^{n-1}2^jv_j.
\]
Here \(F_0\) accounts for the two independent within-word quadratic
terms. Let \(u,v\) be independent uniform words. For a primitive
additive character of modulus \(D=2^r\), \(1\le r\le n\), first
condition on \(X,Y\bmod D\); both are uniform and independent.
The two inner conditional phase averages \(a_s,b_t\) have
absolute value at most one. The Fourier coefficient is
\[
   \widehat\mu(a/D)
    =D^{-2}\sum_{s,t\bmod D}a_s b_t \exp(2\pi i a s t/D),
       \qquad (a,D)=1.
\]
The \(D\times D\) matrix in this bilinear form has operator norm
\(\sqrt D\) by character orthogonality. By Cauchy–Schwarz,
\[
  |\widehat\mu(a/D)|\le D^{-1/2}.
\]
This does not require any cancellation in either unary phase.
For \(D=2\), that bound improves to \(1/2\): the conditional averages
are real in \([-1,1]\), and the maximum occurs at a sign-vector
corner. The dot products of \((1,1),(1,-1)\) against a sign
vector are \(0,\pm2\) in some order, giving the bound.

Use **three independent packets**, with at least \(k\) zero bits
between packets, so cross-packet coefficients vanish modulo \(q\).
Their total value has law \(\mu^{*3}\). At any \(t\bmod q\), the
sum of the absolute values of its nonconstant Fourier contributions
in the numerator of Fourier inversion is bounded above by
\[
 \frac18+\sum_{r\ge2}2^{r-1}2^{-3r/2}
    =\frac18+\frac{1}{4(1-2^{-1/2})}<1.
\]
Every residue of \(F\bmod q\) therefore occurs. The same argument
gives a uniform positive constant times \(1/q\) for the probability
of any prescribed residue. This is an existence theorem on a
finite *literal* packet space: no heuristic randomness in tournament
counts and no analytic NT hypothesis is needed.

## Exact decoder and size

Use an initial \(k\)-zero reset, then the three packets, two
separators of \(k\) zeros, and a final reset of \(k\) zeros; finish
with the prior two-bit selector \(00\) or \(10\) and append 1 to both
addresses. Every realization has the same total length and the same
actual complementary sum before closure. Given the required residue,
evaluate \(F\bmod q\) for all \(q^2\) packets and keep one witness
per value. Enumerate all ordered pairs of supported values, storing
one witness for each sum modulo \(q\); look for the complementary
one-packet value. Fourier inversion proves this exact decoder
succeeds. It costs \(O(q^2\operatorname{poly}(k))\) bit operations
and \(O(q\operatorname{poly}(k))\) storage without FFT or floating point.

The bit count is
\[
 B=4k+3(2(k-2)+1)+2=10k-7 .
\]
Selecting the prior mod-3 center index and the mod-4 terminal choice,
then using the packets for \(K\bmod 2^k\), gives every odd target
class modulo \(3\cdot 2^k\), keeping both exact \(+2\) steps.

In the prior square-crossing bracket use \(q_{\rm grid}=2^B\),
separate from \(q=2^{k-2}\). That argument only needed fixed guard
length, \(q_{\rm grid}\le N^{1/24}\) and \(q_{\rm grid}\cdot
2^k\) dividing the complementary dyadic width \(g\asymp
N^{1/6}\). With \(k=\lfloor(\log_2 N)/240\rfloor\), all conditions
hold for large \(N\), including the near-midpoint constraint
\(q_{\rm grid}^2=o(\sqrt N)\). Consequently one physical three-count
window on each side of \(N\) has a common designated member
congruent to \(N\bmod 3\cdot2^k\), with absolute error \(O(N^{3/8})\).
The precision modulus is now polynomial in \(N\), not
\(2^{O(\sqrt{\log N})}\), within the same budget. This still
does not make the modulus larger than the residual.

**Reusable mechanism:** factor the mixed interaction as a perfect
finite Fourier pairing to control arbitrary internal phase, then
compose a bounded number of orthogonally isolated packets. This is
different from transporting separately selectable high-bit controls:
it works precisely because the new decoder is not triangular in the
individual toggles, escaping the quadratic separation theorem.

## Optional coexistence with a fixed odd modulus

This linear cost also holds in \(k\) when a fixed odd factor \(m\)
(including the mod-3 center requirement) must be preserved. Use the
[previous odd selector](../chanpak_collatz_norm_clock_20261115/FIXED_CRT.md)
before the packets, with an extra mod-5 condition on the sum if
needed. Replace every packet by the packet followed by at least \(k\)
zeros, padding with more zeros until the total length \(\ell\) is a
multiple of \(\mathrm{ord}_m(2)\). Repeat that **same** padded packet
\(m\) times. The dyadic increments add and multiply the packet value by
\(m\), a unit modulo \(q\). At the same time the repeated operation on
\((U,V,W,1)\bmod m\) is the identity. Indeed the even-\(U\) sheet
is both input and output, the address multipliers are
\(2^\ell=1\bmod m\), and \(W=p(U)+p(V)\) gains only a linear form
and a constant. Thus the matrix is \(I+J\) with \(J^3=0\).
Since \(m\) is odd, \(m\mid {m\choose2}\), hence \((I+J)^m=I\)
modulo \(m\). The same padding length works for arbitrary packets.

Three such repeated packets, starting and ending with zeros where
needed, have total length \(O_m(k)\). Choose the required packet sum
by dividing the dyadic target correction by \(m\) modulo \(q\).
The odd selector, the mod-4 terminal choice, and the feasibility of
the bracketing interval follow by the same order of operations as in
the fixed-CRT construction. This proves simultaneous matching modulo
\(m2^k\) with \(k\ge c_m\log_2 N\) and the original
\(O_m(N^{3/8})\) absolute budget (interpreting \(m\) as enlarged to
a multiple of 3). It is not uniform in \(m\).

Replay: run the adjacent `three_packets_check.py` proof control.

## Abstract transfer test

For any prime p, uniform independent X,Y on a prime-power cyclic
ring, and arbitrary unary phases arising from integer-valued
F(X),G(Y), the same pairing estimate applies to
F(X)+G(Y)+XY. For odd p its triple sum has full support: the
sum of nontrivial Fourier moduli after cubing is bounded by

    (1−1/p) Σ_(r≥1) p^(−r/2) = (sqrt(p)+1)/p < 1.

For p=2 the generic constant exceeds one, but the sharper
two-state bound above repairs precisely the problematic mode.
This transfer requires a full bilinear pairing and independent
packets. The binary physical proof additionally needs the packet
separators to kill cross terms; one cannot use binary zero padding
to infer the same physical isolation at odd prime powers.
