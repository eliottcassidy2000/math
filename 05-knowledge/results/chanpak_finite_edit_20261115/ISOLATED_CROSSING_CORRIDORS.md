# Infinite corridors of adjacent crossing pairs still have sparse support

This is the unrestricted-length, **isolated zero-gap** refinement of
the descending anchor chain. Its sparsity has a different cause from
the fixed-budget theorem: unboundedly many edits are allowed here.
A **positive gap** has at least one unmarked vertex; a zero gap is
between two consecutive anchors. We require that zero gaps not be
consecutive in the anchor chain, but do not bound their number.
The prior general-gap lane treats a larger chain model; the sharper
exponent here only applies to this subclass.

## Exact geometry and transfer

Start with a transitive order and anchors \(a_0<\cdots<a_q\), with
\(x_i\ge0\) ordinary vertices in the \(i\)-th gap. Reverse precisely
the consecutive-anchor arcs \(a_i\to a_{i-1}\). Strip any exterior
padding. Every directed cycle has an inclusive anchor interval
\([a_i,a_j]\). Indeed, across any cut between successive anchors
in this interval it must cross backward at least once; there is
only one backward edge crossing that cut. Therefore the cycle uses
the entire descending chain \(a_j\to\cdots\to a_i\), followed
by one increasing path through a subset of the gap vertices
inside the interval. No other anchor returns along the ascent.
If its anchor interval spans \(s\) edges and there are \(K\)
gap vertices inside, the parity condition on the increasing
path's subset size is \(s\bmod2\). The weighted odd-cycle
contribution from this interval is \(2^K\) for \(K>0\), and
\(1+(-1)^s\) for \(K=0\), using the
[odd-cycle expansion](https://arxiv.org/html/2307.05569#S1.SS8.SSS3).
Two cycles are vertex-disjoint exactly when their inclusive
anchor intervals are disjoint. So a subset of selected anchor
edges splits uniquely into its runs, each run assigned the
weight just given, and \(H\) is the sum of the run products
over all such subsets.

Now assume the zero gaps are isolated. The only zero-only run that
could occur is one zero gap by itself. Its weight is zero. Thus,
first choose the positive gaps, with fugacity \(a_i=2^{x_i}\)
at each chosen site \(i\) and spin \(\sigma_i\in\{0,1\}\) recording
this choice. Each zero gap can be selected exactly when a neighboring
positive gap is chosen; whether selected or not introduces no other
weight. An interior zero between positive sites consequently has
matrix \(B_{\sigma,\tau}=1\) for \(\sigma=\tau=0\), and \(2\) otherwise:
\[
B=\begin{pmatrix}1&2\\2&2\end{pmatrix}.
\]
An exterior zero acts as if adjacent to a fixed spin zero beyond
the end. Consecutive positive gaps have no connecting factor; the
all-one matrix is rank one, so the count factors into independent
maximal alternating blocks. Every block has positive sites
\(a_1,\ldots,a_r\), one isolated zero between each adjacent pair,
and optionally an exterior zero at either end. This is an exact
multiplicative decomposition, not a relaxation. The two all-zero
boundary cases (no gaps, or one gap) give \(H=1\).

There is a literal geometric interpretation: for every zero gap,
swap its two adjacent anchors in the reference transitive order.
These are disjoint adjacent swaps. The zero backarc becomes forward;
the physical tournament is unchanged, and all other original
backarcs stay backward. An internal isolated zero now turns its
adjacent positive backarcs into a crossing pair with central gap
of size zero. For a block with exactly two positive gaps and no end
zeros this recovers the \(y=0\) crossing formula
\(1+2a_1+2a_2+2a_1a_2\). It illustrates why crossing pairs can compose
and still fail to supply absolute coverage.

## A \(13/20\) support exponent despite unlimited length

No endpoint zero can lower a block's count. With neither endpoint
zero, initialize the transfer row by \(U=1,\ V=a_1\), where the
coordinates describe the current spin. After appending another
positive site \(a\ge2\) with a zero before it, the update is
\[
U'=U+2V,\qquad V'=2a(U+V).
\]
Let \(p=U/(U+V)\), initially at most \(1/3\). Then the multiplier of
the total \(U+V\) is \(2(a+1)-p\), and
\(p'=(2-p)/(2(a+1)-p)\le1/(a+1)\le1/3\).
Therefore its Hamilton count is at least
\[
\frac12\prod_{i=1}^r (2a_i+5/3).
\tag{1}
\]
This includes \(r=1\) because \(1+a_1\ge(2a_1+5/3)/2\).
It is an ordered growth bound, not just a count of selectable vertices.

For \(\sigma=13/20\), define \(w(x)=2^{x+1}+5/3\) and
\(R=\sum_{x\ge1}w(x)^{-\sigma}\). Then \(R<1\).
For an exact certificate, the first eight terms are bounded above
by \((3239,2289,1547,1017,659,424,271,174)/10000\),
as checked by raising positive numbers to the twentieth power.
For the tail, \(w(x)>2^{x+1}\) and \(2^{-\sigma}<16/25\), so its
sum is less than \((16/25)^{10}/(1-16/25)\).
The sum of these two rational bounds is less than \(0.995\).

Consequently the sum of \(H^{-\sigma}\) over all literal connected
blocks is bounded by \(4\cdot 2^\sigma \sum_{r\ge1}R^r<\infty\),
where four accounts for endpoint choices. In particular the
series converges over the distinct factor values \(A>1\) of blocks.
The unordered multiplicative products of block counts have
\[
\sum_{\text{distinct products }N} N^{-\sigma}
\ \le\ \prod_{A\text{ distinct block value}>1}(1-A^{-\sigma})^{-1}
\ <\ \infty .
\]
The final convergence follows from summability of \(A^{-\sigma}\);
different factorizations of one integer can only increase the right
side. Thus the whole isolated-zero corridor realizes only
\(O(X^{13/20})\) distinct values at most \(X\), in fact
\(o(X^{13/20})\) by using the summable series tail. This is a support
statement even with arbitrary numbers and lengths of blocks.
The decoder remains exact on a prescribed word, but complete
Chan–Pak coverage must add qualitatively different composition
than these synchronized crossings.
