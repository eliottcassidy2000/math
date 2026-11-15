# When the dyadic pattern cannot be packed: a separation theorem

This is a structural obstruction to a precise architecture, not a lower bound
for all tournament constructions. It applies to the checksum
\(p(2u)=p(u)+u+[u\text{ odd}],\ p(2u+1)=p(u)\) used by the physical
three-center bank. The [isolated-column construction](../chanpak_collatz_norm_clock_20261115/WALSH.md)
gets \(O(k^2)\) guard length for \(k\) bits of precision. Can one shorten
it to linear length simply by allowing each bit to depend on the lower
bits rather than strictly isolating the pairs? **No**, if each digit
still has one independently flippable column or disjoint interval
of flippable columns. Interleaved controls and target-dependent
selection are not covered by this no-go.

## A mixed carry computes a valuation of distance

Fix at each column the sum of the two address bits. At a sum-one column
write \(b_i\) for its bit on address 1; the other address gets \(1-b_i\).
The part of \(W=p(x)+p(y)\) depending on two active columns \(i<j\) is
\[
  w(j-i)(b_i\mathbin{\mathrm{xor}} b_j),\qquad
  w(d)=2^{d-1}+[d=1].
\]
Terms crossing fixed columns are independent of every active bit.
A mixed difference in the two binary choices has absolute value
\(2w(d)\) and valuation \(\max(2,d)\). Thus passing from pairwise
coincidence to composition needs the mixed difference, not
just a one-coordinate readout.

More generally, let controls \(s,t\) act in two disjoint intervals,
all changed bits of \(s\) below all changed bits of \(t\). In each
interval, the two settings of the control specify arbitrary patterns
with the given column sums. Write \(\delta_i\in\{-1,0,1\}\) for
the change of the first control's bit at column \(i\), and \(\eta_j\)
for the second. Let \(u=\max\{i:\delta_i\ne0\}\),
\(v=\min\{j:\eta_j\ne0\}\), \(d=v-u\). For \(d\ge2\), the mixed
checksum difference up to sign equals
\[
 \sum_{i,j}2^{j-i}\delta_i\eta_j
  =2^d\left(\sum_i2^{u-i}\delta_i\right)
        \left(\sum_j2^{j-v}\eta_j\right).
\]
Both factors in parentheses are odd, so the valuation is **exactly d**
no matter how complicated each switch is within its interval. At
\(d=1\), the nearest interaction acquires another factor of two;
interactions with \(j-i=2\) can cancel it. For example \(\delta_{u-1}=-1,
\delta_u=+1,\eta_{u+1}=+1\) gives exact zero mixed difference.
There is no universal valuation claim for this special gap.

## The quadratic obstruction

For a triangle of physical outputs on the prior exact bank, fix both
addresses modulo \(2^k\), fix the center index \(j\) by the mod-3
condition, and vary only a guard above the fixed low bits. Then modulo
\(2^k\), the difference between counts is exactly 45 times the checksum
difference. Since 45 is odd, mixed valuations are unchanged. (This is
also a direct application to a stand-alone checksum clock.)

Say a control of level \(\ell\), \(1\le\ell<k\), is **a triangular
2-adic digit** if flipping it always changes the output by
\(2^\ell\bmod 2^{\ell+1}\), whatever the other controls are set to.
This is the weak condition sufficient for greedy programming from
the lowest bit upwards. Taking another control's difference of this
congruence shows the mixed difference is divisible by \(2^{\ell+1}\).
For two levels \(a,b\), it must therefore be divisible by
\(2^{1+\max(a,b)}\). In particular any two singleton switches at levels
\(a,b\), with \(\max(a,b)\ge2\), require distance at least
\(1+\max(a,b)\). This includes neighboring controls in their spatial
order.

For interval switches assume facing changed bits of any two neighboring
controls have gap at least 2. The exact factorization then gives the
same bound on every facing gap. Sum over neighbors. The sum of
\(\max(a,b)\) over consecutive pairs of any permutation of \(1,\ldots,k-1\)
is at least \(\sum_{\ell=1}^{k-1}\ell-1\): root the path at the control
of level 1, then charge each other level to its unique edge toward the
root. Each edge maximum covers its charge. Thus the span of the guard,
even ignoring intervals' internal widths, is at least
\[
 (k-2)+\frac{k(k-1)}2-1=\frac{k(k+1)}2-3 \qquad (k\ge3).
\]
For singleton controls the same bound follows directly. Even arbitrary
interval controls **do not evade quadratic width by touching**. Divide
the spatially ordered intervals into its even and odd subsequences.
Any two consecutive members of either subsequence have at least one
other control between them, so their facing columns are at distance
at least two and the exact valuation applies. The sum of the facing
gaps in each subsequence is at most the guard's span. Applying the
same root-charge bound to both subsequences, with minima \(m_0,m_1\),
gives
\[
  2\,\mathrm{span}\ \ge\ \frac{k(k-1)}{2} - (m_0+m_1) + (k-3)
  \ \ge\ \frac{k(k-1)}{2}-3,\qquad k\ge 3,
\]
where one minimum is 1 and the other is at most \(k-1\).
This proves that the prior quadratic exponent in the guard length is
necessary for this triangular disjoint architecture, not an artifact
of requiring complete independence of output bits.

## Where the analogy works and where it stops

The distance appears as an exact mixed valuation and so survives
repeated composition of disjoint carries. The \(d=1\) exceptional
coefficient is a concrete place where naive periodic transport fails;
it is not safe to extend the generic valuation to that case. A scheme
using overlapping or interleaved controls could cancel multiple
leading terms allowing shorter guards; to demonstrate it, exhibit
the cancelling mixed differences and prove constructive decoding.
Another route is a nonlinear decoder not requiring each elementary
flip to be a triangular digit, even when all switches commute as
physical actions. Neither route is excluded by this theorem.
