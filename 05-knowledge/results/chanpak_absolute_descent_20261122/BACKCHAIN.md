# A backedge chain compiles a positive three-register count

Take boundary vertices \(b_0,\ldots,b_k\) in increasing order, with \(m_i\ge0\) neutral spine vertices strictly between \(b_{i-1}\) and \(b_i\). Start with all arcs forward and reverse precisely \(b_i\to b_{i-1}\) for \(1\le i\le k\). The reversed chain remains fixed while gaps are chosen.

## Exact cycle partition and the transfer

Every directed cycle has extremes \(b_a,b_b\), where \(a<b\). It crosses each cut between them backwards, and the only backwards edge through a cut in gap \(i\) is \(b_i\to b_{i-1}\). It therefore traverses the whole backwards chain from \(b_b\) to \(b_a\), then an increasing path through some of the \(\sum_{i=a+1}^b m_i\) neutral vertices. Its vertex count is \(b-a+1\) plus the neutral subset size. A singleton backedge with no chosen neutral vertex cannot be a cycle, and its length would be even, so it contributes nothing below. All other candidate cycles of odd length are valid. (A forward path through an intervening boundary cannot occur because it has already been visited on the backwards chain.)

By the [Grinberg–Stanley odd-cycle formula](https://arxiv.org/abs/2307.05569), a family of vertex-disjoint odd directed cycles has weight \(2\) per cycle. Thus \(H\) is the partition function of disjoint closed integer intervals \([a,b]\) on the boundary indices. An interval with neutral mass \(m>0\) has weight \(2^m\) by parity balance. One with \(m=0\) has weight \(2\) when \(b-a\) is even and zero otherwise.

Let \(F_i\) be the partition function up to \(b_i\), with \(F_{-1}=F_0=1\). Let \(A_i\) be the terminal-interval contribution if **every** interval has hypothetical weight \(2^m\) but earlier intervals use the true weights in \(F_{a-1}\). Let \(S_i\) be its correction for the intervals in the final zero-gap run. Then
\[
F_i=F_{i-1}+A_i+S_i,\quad
A_i=2^{m_i}(F_{i-1}-S_{i-1}),\qquad
S_i=\begin{cases}0&m_i>0,\\-F_{i-2}-S_{i-1}&m_i=0.\end{cases}
\]
Indeed in \(A_i\) the new shortest interval contributes \(2^{m_i}F_{i-2}\), and all other terms extend the preceding \(A_{i-1}\). In a trailing zero run the correction to the hypothetical weight 1 is \(+1\) for even-length intervals and \(-1\) for odd-length intervals, giving the signed relation. The initial values are \(A_0=S_0=0\).

The signed correction disappears in the coordinates \(D_i=F_i-F_{i-1},\ U_i=-S_i,\ V_i=F_{i-1}+S_i\). Start at \((D,U,V)=(0,0,1)\); the output is \(D+U+V\). The complete update is nonnegative:
\[
\begin{array}{ll}
m_i=0:&(D,U,V)\mapsto(D+2U,\ V,\ D+U),\\
m_i>0:&(D,U,V)\mapsto(2^{m_i}(D+2U+V),\ 0,\ D+U+V).
\end{array}
\]
These are actual count constructions, not only formal residue updates.

## Products, one-zero coupling, and exact backward tests

When every gap is positive, the state resets \(U=0\) at each step and the count is the product \(\prod_i(1+2^{m_i})\). All vertices lie in a strong component: \(b_0\) can reach all boundaries via forward two-edge paths through positive gaps, the backwards chain returns all boundaries to \(b_0\), and each neutral vertex lies between its two boundaries. This uses \(1+k+\sum m_i\) vertices, \(k-1\) fewer than taking the ordinal product of the usual near-transitive factors separately. The uniform width-two specialization reaches Busch's strong minimum \(5^k\) on \(3k+1\) vertices, so this part should also be viewed alongside Moon's extremal construction.

If exactly one zero gap is placed between two positive gaps with factors \(a,b\), these two factors become \(2ab-1\), while earlier and later positive factors multiply as before. More generally, take words that begin and end with positive gaps and have only isolated zero gaps. At each positive boundary keep \(F=D+V,\ G=V\), so \(F\) is the ordinary count and \(G\) the preceding boundary prefix count. Appending positive factor \(1+w\), \(w=2^m\), either directly or after one zero gives:
\[
(F',G')=((1+w)F,F)
\quad\text{or}\quad
(F',G')=(2(1+w)F-G,\ 2F-G).
\]
In the second case, \(F'-G'=2wF\). Because \(F\) is odd, the **ordinary integer** difference decodes \(m=v_2(F'-G')-1\), \(F\) as the odd part of the difference, and \(G=2F-G'\). An inadmissible recovered state rejects that branch; no modular search is needed. In the first case \(F=G'\) and the ratio \((F'-G')/G'\) must be exactly a positive power of two other than one. The direct step forgets the old auxiliary register, so this is a decoder for the coupled step and a necessary test for the reset step, not a scalar inverse for all odd targets.

Every new positive gap increases the count by at least a factor 3. At such a gap the multiplier is at least \(2^m\), so the sum of positive gap widths is bounded by \(\log_2 H\), and the total number of gaps (isolated zeros included) by \(2\log_3 H\). Thus a successful absolute certificate from this grammar has logarithmic order. What is missing for cofiniteness is a proof that each sufficiently large odd target admits a compatible terminal register and a descending sequence; the forward construction and coupled inverse do not give that coverage.

## Exact finite census, and a non-apex witness

At height at most \(100001\), this entire backchain grammar (now
allowing arbitrary zero runs) gives 13,263 distinct counts. Of
these, 6,437 are absent from the one-apex/transitive-spine chart;
the latter gives 25,117 distinct counts in the same range. These
are grammar comparisons, not new gaps or a census of all tournaments.
For example \(99991\) has backchain word
\((2,0,0,2,0,0,0,2,0,1,0,1)\) on 21 vertices but no one-apex word.
The one-apex check is complete because trailing ones can be stripped
without changing its read, and an even address \(2u\) has
\(f(2u)\ge u\). For the backchain enumeration the output never
decreases; the first zero step preserves height 1, but followed
by a zero it gives exactly the one-positive-width-one state,
and followed by positive width \(m\) it gives exactly the
one-positive-width-\((m+1)\) state. Every other step increases
the height strictly. Thus one can discard the initial zero,
prune at the height cutoff and enumerate all remaining words
and all possible positive widths.

With isolated zeros strictly between positives the grammar is
still strong: one can cross such a gap forwards by skipping its
right boundary and then returning to it on the reversed arc from
the next boundary, which is reachable through the adjacent
positive gap. More explicitly, \(b_{i-1}\to b_{i+1}\to b_i\) at
zero gap \(i\), and each positive gap allows forward passage through
one interior vertex. Backwards all boundaries reach \(b_0\).
