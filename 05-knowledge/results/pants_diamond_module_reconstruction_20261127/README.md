# A four-vertex certificate for every proper module reversal

2026-11-27. This closes the open regular case of the
[local reversal calculation](../pants_module_reversal_visibility_20261127/README.md).
A module is a vertex set with uniform arrows to each outside vertex.
For a nontrivial proper module, the exact thresholds are:

| Internal tournament | First switched path count that detects reversal |
|---|---|
| Irregular | Three vertices |
| Regular of order at least seven | Six vertices |
| Regular of order three or five | None, even allowing every path and cycle length |

Counts are for directed paths on a specified number of distinct vertices,
summed over all subsets of that size. The data keep the cut labels.

## 1. What the square detects

Write the signed adjacency as
\[
 A=\begin{pmatrix} C&\mathbf1 b^\top\\-b\mathbf1^\top&E\end{pmatrix}
\]
relative to the module and its exterior. Replacing \(C\) by \(-C\)
changes the upper-right block of \(G=A^2\) by \(-2(C\mathbf1)b^\top\).
It changes no other independent block. Hence it preserves \(G\)
iff the module is regular, since the exterior and \(b\) are nonempty.
The [short-path square formulas](../pants_three_path_square_inverse_20261127/README.md)
give a three-path coefficient of \(G_{ij}/2\) at each pair and show
that regular reversal preserves every path response through five
vertices. Counts on one or two vertices are independent of arrows.
For orders three and five, regularity gives the carousel, whose
cutwise reflection was established in the preceding note.

## 2. The diamond bit

Fix one exterior contact \(d\) and put \(\tau=A_{ad}\) for internal
\(a\). Write \(p=2m+1\) for the order of the regular module.
For a triple \(T=\{a,b,c\}\), let \(\Delta_T\) be the change when
reversing the module in the six-path Fourier coefficient on
\(T\cup\{d\}\). The preceding quartet calculation gives
\(\Delta_T=-\tau(3H_T-J_T)/2\), where
\[
 H_T=\sum_x C_{ax}C_{bx}C_{cx},\qquad
 J_T=(C^2)_{ab}(C_{ca}+C_{cb})+
 (C^2)_{ac}(C_{ba}+C_{bc})+(C^2)_{bc}(C_{ab}+C_{ac}).
\]
Let \(t_T=1\) for a transitive triple and \(0\) for a cyclic one.
Let \(r_{ab}\) be the number of common predecessors of \(a,b\).
Equal outdegrees make it also their common successors, so
\((C^2)_{ab}=p-2-4r_{ab}\).
On a cyclic triple \(J_T=0\); on \(a\to b\to c,\ a\to c\),
\(J_T=2((C^2)_{bc}-(C^2)_{ab})\), a multiple of eight.

To compute \(H_T\) modulo eight, expand the product of the three
signs outside \(T\) as \(\prod (1-2u_i)\). There are \(3m-3\)
negative single incidences, and \(r_{ab}+r_{ac}+r_{bc}-t_T\)
negative pair incidences, because the removed internal source
contributes one pair precisely in a transitive triple. Hence
\[
 H_T\equiv4(1-m+r_{ab}+r_{ac}+r_{bc}-t_T)\pmod8 .
\]
In particular \(\Delta_T\) is even and
\(b_T=\Delta_T/2\bmod2=m-1+r_{ab}+r_{ac}+r_{bc}+t_T\bmod2\).

A **diamond** is a directed triangle with a source or sink attached.
Let \(D_U\) be its indicator on a four-set \(U\). Then
\[
 D_U=\sum_{\substack{T\subset U\\|T|=3}} b_T,\qquad
 b_T=\sum_{x\notin T}D_{T\cup\{x\}}\quad\text{in }\mathbb F_2 .
\]
For the first, each pair term occurs twice and the constant
cancels. The remaining sum is the parity of the number of
cyclic triples: of the four tournament types of order four,
only the two diamonds have an odd number. For the converse,
a pair \(a,b\) extends to \(m-r_{ab}\) cyclic triples (use the
outneighbors of the loser not dominated by the winner), and
to \(m-1+r_{ab}\) transitive triples. In the sum over extensions
of \(T\), its own face appears \(p-3\) times, an even number.
The three other pair sums add to \(m-1+\sum r_{ab}+t_T\)
modulo two, where each pair has omitted the third vertex of
\(T\). This yields the second identity.

A diamond therefore certifies a nonzero change, equal to two
modulo four, on one of its four triples and **any** outside
contact. Conversely a change equal to two modulo four certifies
a diamond extension. In particular the number of flagged triples
is at least the diamond count divided by \(p-3\) for \(p>3\).

## 3. The diamond-free branch and its exact spectrum

A diamond-free tournament has transitive in- and outneighborhoods
(the [locally transitive characterization](https://arxiv.org/abs/1503.04124)).
In the regular case this forces a carousel: choose vertex 0 and
order the two neighborhoods as \(a_1\to\cdots\to a_m\) and
\(b_1\to\cdots\to b_m\), respectively out and in. If
\(b_j\to a_i\to b_k\) for \(j<k\), then \(b_j\) dominates the
triangle \(0,a_i,b_k\), a contradiction. So the \(b\)'s beaten by
\(a_i\) form an initial segment, of length \(i\) by regularity.
These are precisely the carousel arrows in the circular order
\(0,a_1,\ldots,a_m,b_1,\ldots,b_m\). Regular orders three and five
are automatically diamond-free since both neighborhoods have
size at most two.

For a carousel \(H_T=0\) for every \(T\). Use its increasing
orientation switched with signs \(s_i=(-1)^i\), equivalent
to the carousel by \(i\mapsto i/2\) modulo \(p\). For ordered
endpoints \(a<b<c\), \(H_T\) is \(s_as_bs_c\) times the sum of
\(s_x\) on the four complementary intervals with signs
\(-,+,-,+\). In terms of the prefixes \(F(j)=\sum_{i<j}s_i\)
this is
\[
 F(p)-[F(a)+F(a+1)]+[F(b)+F(b+1)]-[F(c)+F(c+1)]=0.
\]
The cyclic triples have \(\Delta_T=0\). A transitive one has
circular successive gaps \(u,v>0\) from source through the
middle to the sink, with \(u+v\le m\). By the [carousel square
formula](../pants_module_reversal_visibility_20261127/README.md)
its change has magnitude \(4|u-v|\). For \(1\le k\le m-2\)
there are exactly \(2p\lfloor(m-k)/2\rfloor\) triples with
magnitude \(4k\), by the two gap orders and the \(p\) choices
of source. Thus the total detector count is
\(p[\binom m2-\lfloor m/2\rfloor]\), and the total squared
change is \(32p\binom{m+1}{4}\), since the sum of \((u-v)^2\)
over \(u+v=s\) is \(2\binom s3\). For \(p\ge7\) the
detector count is positive, completing the classification.

## 4. Simultaneous recovery of arrows between known modules

Consider a fixed labelled partition into at least two regular
modules, each of odd size at least three. Suppose each internal
matrix is known up to reversal, but the whole tournament on
the modules (the quotient) may be unknown. Call a module large
if it has size at least seven. Assume there is at least one large
module. Then the labelled three- and six-path responses together
recover every arrow touching a large module, up to a common
global converse. They recover the entire quotient if at most
three modules are small, or more generally under the rank
criterion below. Small internal directions remain invisible.

First take a known detecting triple in a large module \(i\)
and an exterior vertex \(d\) in module \(j\). On the quartet,
the absolute Fourier coefficient can be written
\[
 F(G,C_i,T,d)+s_i A_{ij}(3H_T-J_T)/4,
\]
where \(A_{ij}\) here denotes the quotient sign, and \(F\) is
known from \(G\), the reference internal matrix \(C_i\) and
the labels. This follows from the preceding
[two-run expansion](../pants_module_reversal_visibility_20261127/README.md):
the pair terms within the triple use only sign-independent
internal products and the square of the quotient sign. In
the quartic column product, the exterior-center contribution
is exactly \(-G_{ad}\) for \(a\in T\), because the internal
row sum at \(a\) vanishes. Thus there is no additional unknown
quotient dependence. As the multiplier is nonzero, one
quartet reads \(s_i A_{ij}\). Read it at both large ends of
their mutual arrows to determine all large signs relative to
one root. Fixing that root sign as a global converse convention
determines the large signs and all their incident quotient arcs.

For the residual small set, write \(Z\) for the matrix of known
large-to-small signs, \(W\) for the diagonal matrix of small
module sizes, and \(D\) for the difference of two proposed
small quotient skew matrices. Computing the large-to-small
block of \(G\) gives \(ZWD=0\): intermediate vertices in a
third module contribute their size times the corresponding
quotient product, while same-module contributions vanish by
regularity. With \(r\) small modules, if
\(\mathrm{rank}(Z)\ge r-1\), the skew matrix \(D\) has rank at
most one, hence is zero. This covers \(r\le2\) (the empty
case is trivial) and \(r=3\) with two independent rows of \(Z\).

The remaining \(r=3\) case has rank one. Its sign rows are
proportional. If \(D\ne0\), its one-dimensional kernel contains
\(Wz^\top\) for a row \(z\) of \(Z\). All entries of this
vector are nonzero. The kernel of a nonzero three-by-three
skew matrix is spanned by its three signed off-diagonal
entries in complementary order. As the entries here are
\(0,\pm2\), all three must have magnitude two, and all
small sizes are equal. Switch each entire small block by its
bit in \(z\), in both putative tournaments. This just translates
both cut-response functions. The three small blocks now form a
proper module, because each large row to them is uniform.
Their quotient arcs all reverse, and the kernel condition
makes their quotient a directed triangle in this gauge.
Align their small internal signs by cutwise invisible
reversals if necessary. Reversing that directed triangle
of blocks is then indistinguishable from reversing the
entire composite module and undoing the three invisible
internal reversals. But the composite is regular of order
nine or fifteen, so the six-path response sees it by
Sections 2–3. Contradiction. This completes the three-small
case even when the linear square test alone is insufficient.

## Replay

Run `node replay.js`. It verifies the incidence identities and
carousel spectra through order 35, counts actual six-paths over all
cuts in the eight-vertex Paley and carousel hosts, and gives a
six-path discrepancy for a concrete three-small-block residual
whose labelled squares agree. The residual consists of three
triangles forming a directed triangle of blocks, joined to a
seven-vertex carousel. The separate three-path control detects
an irregular reversal.
