# A shared-apex seam with an exact two-reservoir inverse

This construction does not use disjoint matched return arcs. Its many
return arcs share a vertex, which makes the formula a sum of **all**
out-to-in intervals rather than a single color profile of the reservoir.
The main concrete payoffs are a universal local transfer, exact seam
composition, and a prescribed-count arithmetic decoder. Nothing here
establishes eventual coverage of primes or seven times primes.

## 1. Uniform in each block, nonuniform on the whole reservoir

Take any nonempty tournaments \(B_1,\ldots,B_k\), with all arcs between
them forward in this order. Give a single apex \(a\) arbitrary contact
bits \(b_i\) to the blocks: \(b_i=1\) means \(a\to B_i\), \(b_i=0\) that
\(B_i\to a\). Write \(F_m(B)\) for the color observable (the sum of the
products of the Hamiltonian-path counts on the \(m\) color classes;
empty classes have count 1), including \(F_0(B_i)=0\). By the
odd-cycle expansion, \(F_w(B)\) is the polynomial
\[
 \sum_{\mathcal C}2^{|\mathcal C|}
             w^{|B|+|\mathcal C|-|\bigcup\mathcal C|},
\]
where \(\mathcal C\) runs over disjoint odd directed cycle packings.
Use \(F_1=H\). For every \(m\ge1\), write \(h_i=F_m(B_i)\),
\(p_i=F_{m+1}(B_i)\), \(r_i=F_{m-1}(B_i)\). The augmented graph \(T\)
has
\[
\frac{F_m(T)}{m}=
\prod_i h_i+
\sum_{\substack{i<j\\ b_i=1,\ b_j=0}}
 \left[
 (p_i-h_i)(p_j-h_j)\prod_{i<t<j}p_t+
 (h_i-r_i)(h_j-r_j)\prod_{i<t<j}r_t
 \right]
 \prod_{t<i\text{ or }t>j}h_t.                                      \tag{1}
\]
Every term is nonnegative. No assumption on directed cycles within
the blocks is used.

**Proof.** Fix the apex's color (there are \(m\) choices). A packing
not using the apex is the base term, since no directed cycle can cross
forward ordinal blocks. In a packing using the apex, the unique cycle
containing it visits a nonempty ordered path of vertices starting in
block \(i\) of type 1 and ending in block \(j>i\) of type 0. Its selected
vertices have even total cardinality. In each block of the closed
interval we split the apex color into two distinguished slots, selected
and unused. Remaining color classes are unchanged. A selected-slot
subset gets its ordinary Hamiltonian count: its entire path is used
in the apex cycle. The unsigned evaluation is \(F_{m+1}\). The signed
evaluation, by \((-1)^{|S|}\) for the selected subset \(S\), is
\(F_{m-1}\): in the odd-cycle expansion for the selected-slot Hamilton
count, each directed odd cycle of length \(\ell\) contributes
\((-1)^\ell=-1\), and each free vertex has the same color-sign rule,
so the effective color number is \((m-1)+1-1=m-1\).
At the endpoints one must subtract the unselected contribution \(F_m\)
from either evaluation. The two resulting negatives of the signed
endpoint differences cancel. The even-cardinality projector contributes
one half, canceled by the factor 2 for the apex cycle. This gives (1).

For \(m=1\) the middle \(F_0\) term survives only for adjacent blocks.
Write \(h=H(B_i)\), \(\ell=F_2(B_i)\), \(d=\ell-h\) for the new block
and \(b=b_i\). Maintain prefix product \(P\), weighted past start \(E\),
whole prefix count \(K\) (including the apex), and previous bit \(a\).
Then appending the new block updates
\[
 P'=hP,\qquad E'=\ell E+b\,dP,\qquad
 K'=hK+(1-b)(dE+a hP),\qquad a'=b.                                  \tag{2}
\]
Start at \((1,0,1,0)\). Indeed \(E\) sums \(d_i\) over prior type-1
blocks, times \(H\)'s before \(i\) and \(F_2\)'s after \(i\);
the next type-0 block collects precisely these cross terms, with a
separate nearest-neighbor correction. Singleton-specialization
recovers the inherited one-apex inversion law. More generally,
if two prefixes ever agree on \((P,E,a)\) with \(K\) different by 2,
**any** suffix of transitive blocks with any contact bits retains that
exact difference. This sufficient repair criterion does not assert
such a pair exists.

## 2. Compose around a nonuniform seam

For two lists of blocks, let their base products and apex counts be
\((P_U,K_U),(P_V,K_V)\). Let \(E_U\) be the forward left weighted
start count of (2), and let \(G_V\) be the corresponding weighted
finish count, using the reversed list with the bits complemented.
Then, with \(b_{\rm last}(U),b_{\rm first}(V)\) the two exposed bits,
\[
K_{UV}=
 K_U P_V+P_U K_V-P_U P_V+
 E_U G_V+
 b_{\rm last}(U)(1-b_{\rm first}(V))P_U P_V.                         \tag{3}
\]
Indeed, the first three terms collect intervals within each side
once; the separated left and right interval factors multiply for
an interval crossing the boundary. The last term is the only short
interval that crosses it. This is a target-selected baseline
condition involving actual starters and finishers, not merely scalar
output counts \((K_U,K_V)\).

## 3. A prescribed-count factorization with two independent reservoirs

Take the four-block pattern \(1,0,1,0\). The first and third blocks
are singletons; the second is \(B\), and the fourth is \(C\). Write
\(h=H(B),\ell=F_2(B), q=H(C),L=F_2(C)\). Substitution in (2) gives
\[
 H(T)=2\ell L+hL-\ell q+hq,\qquad
 \boxed{2H(T)-3hq=(2\ell+h)(2L-q)}.                                 \tag{4}
\]
Both reservoirs are nonempty and may be arbitrary. The graph is
strong: the first vertex is reached from the apex, every vertex
is reached from that first vertex, and every vertex can reach a
last-block vertex, which points to the apex.

For any two fixed nonempty cores \((h,\ell_0),(q,L_0)\), independently
pad them with forward transitive vertices within their contact classes.
Padding preserves \(h,q\) and replaces \(\ell_0,L_0\) by
\(\ell_0 2^r,L_0 2^s\) respectively. Given an exact target \(N\),
factor \(2N-3hq\). For each positive ordered divisor pair \((f,g)\)
test
\[
    (f-h)/(2\ell_0)=2^r,\qquad (g+q)/(2L_0)=2^s,\qquad r,s\ge0 .
\]
This is a complete exact inverse of the chosen two-core chart, and
may reject targets even when local congruences permit them. If the
second core itself is a variable transitive chain of positive size,
use \(q=L_0=1\) in the divisor test and require \(s\ge1\).

A concrete primitive source is the four-vertex tournament with
out-masks \((6,12,8,1)\), counted from bit 0, giving \((h,\ell_0)=(5,32)\).
For the second reservoir a transitive \(s\)-chain, formula (4) becomes
\[
    N=2^{r+s+6}+5\,2^s-2^{r+5}+5,\qquad r\ge0,\quad s\ge1.
\]
It realizes primes (e.g. \(2069\) at \((1,4)\)) and seven times primes
(e.g. \(32837=7\cdot4691\) at \((3,6)\)). Arithmetic of the two-width
chart further shows that any \(7p\) with \(p\ne3,7\) must be in
\((r,s)\equiv(3,0)\) or \((5,4)\pmod6\). These two residues are
necessary, not a prime-producing guarantee. This family is
not being proposed as dense; a fixed pair of profiles with dyadic
widths still has only logarithmically many width choices in each
dimension subject to any absolute height ceiling. The usable
advance is the factorization with exactly specified offsets and
two *independently* variable reservoir profiles.

**Provenance.** Formula (1) uses the Grinberg–Stanley directed
odd-cycle packing identity and signed coloring already used in the
paired-port packet. Session317/322's one-apex formulas specialize
to transitive singleton blocks; the two-reservoir cyclic triangle in
Session329 has contact pattern \(1,0\), not the \(1,0,1,0\) readout
in (4). The new point is the arbitrary-block interval transfer and
the resulting alternating two-reservoir inverse.
