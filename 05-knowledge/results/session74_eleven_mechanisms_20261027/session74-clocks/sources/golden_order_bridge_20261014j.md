# Golden cyclic words, marked poset cuts, and exact comparison weights

**Status: PROVED scoped; FINITE-EXACT controls; CITED primary context.**
The new instrument transports the eleven golden states into conditioned order
ideals and conditioned linear extensions, retaining the precise closing-edge
defect. It also restores genuine uniform-linear-extension probabilities by
source-bound prefix/suffix counts. It does not settle the 1/3–2/3 conjecture,
Gold Partition, Borsuk(4), or Collatz.

## 1. Inheritance and the actual eleven-state object

The corrected source is
[THM-4592, eleven-squares Fibonacci torus / period-five cycle set](../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md),
with the MISTAKE-584 repairs retained. Its map is the **standard** rational
Collatz map

\[
C(x)=\begin{cases}3x+1&x\text{ odd},\\x/2&x\text{ even},\end{cases}
\]

on rationals of odd denominator. It is not the shortcut map. At period
dividing five, its parity words are exactly

\[
\mathcal C_5=\{w\in\{0,1\}^5:w_iw_{i+1}=0\quad(i\bmod5)\}.
\tag{1}
\]

There are eleven. The golden label is
\(\kappa(w)=\sum_{i=0}^4w_i4^{4-i}\pmod {11}\), and left shift multiplies
the label by four. The three standard-map cycles are

\[
\{0\},\quad
1/13\to16/13\to8/13\to4/13\to2/13\to1/13,\quad
-5\to-14\to-7\to-20\to-10\to-5.
\]

The script independently authenticates all affine fixed points, odd-denominator
parities, labels and shifts. This is a finite rational-cycle dictionary, not
coverage of positive integer orbits. General odd-length bijectivity of the
golden reader is not a dependency here.

The closest earlier transfer warning is
[THM-4503, Collatz–poset height selection](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md):
conditioning an extension set by an arithmetic restriction need not leave the
complete extension set of any poset. The least-used coordinate below is the
**closing edge**, followed by the **correct sampling weight**.

## 2. A lossless marked-cut compiler

For any \(n\ge3\), cut the cycle between positions \(c-1,c\), and write its
local word \(x_i=w_{c+i\bmod n}\). First retain only the path constraints
\(x_ix_{i+1}=0\), \(0\le i<n-1\). Set

\[
b_i=x_i\mathbin\oplus(i\bmod2).
\tag{2}
\]

Let \(F_n\) be the height-two fence with each odd position below its adjacent
even positions. Its ideal indicators obey \(b_{\rm even}\le b_{\rm odd}\).
Thus (2) is a bijection between path-admissible words and ideals of \(F_n\).
The inverse is the same coordinate complement. The retained seam is
\(\eta=x_0x_{n-1}\); cyclic admissibility is exactly \(\eta=0\).

The API stores `(word, cut, ideal, seam)` and reauthenticates all four fields.
Under a left rotation by \(s\), the cut becomes \(c-s\pmod n\); local ideal
and seam remain unchanged. Omitting the cut would change the coordinate map.

**Odd-cycle obstruction, precisely scoped.** Turning every forbidden pair
`11` into an implication by coordinate complements requires opposite
complement choices at the two ends of each edge. This is an edgewise XOR-one
system, solvable iff the cycle length is even. For odd length, one retained
seam is sufficient and zero is impossible **within this coordinate-complement
representation**. This does not prohibit other abstract bijections or posets.

For \(n=5\), the fence has 13 ideals. The closing condition excludes exactly

\[
I=\{0,1,3,4\},\qquad I=\{0,1,2,3,4\}.
\tag{3}
\]

Hence the uniform golden law becomes the uniform ideal law **conditioned** on
excluding (3). Its coordinate probabilities are
\((3,8,3,8,3)/11\). None lies in \([1/3,2/3]\). The unconditioned fence-ideal
probabilities are \((5,10,4,10,5)/13\). These are membership probabilities;
section 4 supplies an actual comparison-event map.

## 3. One Fibonacci partition function, with its observable retained

Give occupation at position \(i\) positive rational activity \(z_i\). The
cyclic partition function is exactly

\[
Z_n(z_0,\ldots,z_{n-1})=
\operatorname{Tr}\prod_{i=0}^{n-1}
\begin{pmatrix}1&z_i\\1&0\end{pmatrix}.
\tag{4}
\]

The matrix index retains the preceding bit; a closed walk contributes its
occupation weight once. Alternatively, subtract from the path partition
function the weight of words occupying both endpoints. This defect weight is
\(z_0z_{n-1}\) times the path partition function on positions \(2,\ldots,n-3\),
with an empty interval contributing one. This also covers \(n=3,4\).

Through (2), the correct ideal weight is
\(\prod_{i\text{ even}}z_i^{b_i}\prod_{i\text{ odd}}z_i^{1-b_i}\), not
\(z^{|I|}\). At constant activity two, the eleven-state partition function
is 31, whereas the latter incorrect weight sums to 73.

For constant activity \(z>0\), put
\(\lambda_\pm(z)=(1\pm\sqrt{1+4z})/2\). Then

\[
Z_n(z)=\lambda_+(z)^n+\lambda_-(z)^n,\qquad
p_n(z)=\frac z n\frac{d}{dz}\log Z_n(z).
\tag{5}
\]

Rotation invariance makes \(p_n\) the occupation probability at any specified
position. At \(z=1\), with Fibonacci convention \(f_0=0,f_1=1\),

\[
p_n(1)=\frac{f_{n-1}}{L_n},\qquad
\lim_{n\to\infty}p_n(1)
=\frac1{\varphi\sqrt5}=\frac{5-\sqrt5}{10}.
\tag{6}
\]

The numerator counts the remaining path after deleting an occupied vertex
and its neighbours. Also \(L_n=3f_{n-1}+f_{n-2}\), so \(p_n(1)<1/3\)
for every \(n\ge3\). In particular, \(p_5=3/11\).

These constants do have a concrete poset connection. The primary source
[Brightwell–Felsner–Trotter, *Balancing Pairs and the Cross Product Conjecture*,
§§1,7–8](https://trotter.math.gatech.edu/papers/97.pdf)
uses the thin order on integers defined by \(i<j\) iff \(j\ge i+2\), whose
bulk comparison limit is \((5-\sqrt5)/10\). Its universal bound additionally
uses correlation inequalities; that proof is not replaced by (4). The next
section gives the elementary finite map behind the shared Fibonacci law.

## 4. An actual comparison map, and an exact closure obstruction

Let \(Q_N\) have labels \(0,\ldots,N-1\), with \(i<j\) iff \(j\ge i+2\).
Every linear extension is the identity order with some vertex-disjoint
adjacent pairs swapped. Indeed, only consecutive labels can be inverted;
two such inversions cannot share a label because that would invert a
comparable pair. Conversely, disjoint adjacent swaps respect all relations.

Thus uniform extensions of \(Q_N\) are exactly uniform independent sets of
the path with \(N-1\) possible swap positions. Occupation of edge \(i\) is
the comparison event that \(i+1\) precedes \(i\). Consequently,

\[
e(Q_N)=f_{N+1},\qquad
\Pr(i+1\prec i)=\frac{f_{i+1}f_{N-i-1}}{f_{N+1}}.
\tag{7}
\]

If both distances to the boundary tend to infinity, (7) tends to (6).
At the first edge it instead tends to \(\varphi^{-2}=(3-\sqrt5)/2>1/3\).
The transfer matrix and the event are shared; the boundary conditions differ.

**Closing-edge theorem.** For \(n\ge3\), select from the extensions of
\(Q_{n+1}\) just those whose swap words also obey the closing constraint
\(w_0w_{n-1}=0\). This is a proper subset: the word with just the first and
last swaps is excluded. Nevertheless, its intersection of order relations
is exactly \(Q_{n+1}\). The identity extension survives, and every singleton
swap survives, so neither orientation of any consecutive pair becomes forced.
All nonconsecutive comparisons remain those of \(Q_{n+1}\).

It follows that this selected family is **not** the full extension set of
any poset on these same labelled elements: a finite poset is the intersection
of its linear extensions. For every \(n\ge3\), all adjacent inversions have
probability \(p_n(1)<1/3\), while nonadjacent comparisons are deterministic.
Thus the selected uniform law has no balanced pair despite preserving every
common order relation. For \(n=5\), 13 extensions of \(Q_6\) become 11, and
the nonconstant probabilities are \(3/11\) and \(8/11\). This is a conditioned
law, not a counterexample to a conjecture about the full uniform extension law.

The two faithful interpretations of the same 13 path words are now explicit:
ideals of the five-element fence, and extensions of the six-element thin
order. They are different posets and different observables unless (2) and
the adjacent-swap map are retained.

**An all-length Gold certificate.** For every \(N\ge3\), compare \(0,1\)
first and \(1,2\) second, with the natural orientations designated true.
The exact four-cell count table for \(Q_N\) is

\[
\begin{pmatrix}f_{N-1}&f_{N-2}\\f_{N-1}&0\end{pmatrix}.
\tag{7a}
\]

The zero cell would require two adjacent swap edges. The other cells follow
by deleting the forced initial edges and counting the remaining path.
The heavy first branch has size \(H=f_N\), while the largest final cell
has size \(M=f_{N-1}\). Therefore the Gold budget holds with equality:
\(H+M=f_{N+1}=e(Q_N)\). In particular the first comparison has minority
probability \(f_{N-1}/f_{N+1}\in[1/3,1/2]\). These are complete-extension
certificates for this family, not a statement for arbitrary posets.

For \(Q_6\), the table is \(\left(\begin{smallmatrix}5&3\\5&0\end{smallmatrix}\right)\).
The cyclic closing condition deletes exactly two orders from its lower-left
cell, giving \(\left(\begin{smallmatrix}5&3\\3&0\end{smallmatrix}\right)\).
Now \(H+M=8+5=13>11\). The marked seam precisely identifies the deletion
which destroys this certificate.

## 5. Restoring genuine extension probabilities by exact flow

For an arbitrary finite supplied poset \(P\), let \(e(P)\) count its linear
extensions. For each ideal \(I\), set
\(f(I)=e(P[I])\), \(g(I)=e(P\setminus I)\). Concatenation gives

\[
\Pr(\text{first }k\text{ labels form }I)
=\frac{f(I)g(I)}{e(P)},\quad |I|=k.
\tag{8}
\]

The directed ideal graph has edge \(I\to I\cup\{x\}\) whenever \(x\) is
available. Give it the exact nonnegative flow

\[
\mathsf F(I,x)=\frac{f(I)g(I\cup\{x\})}{e(P)}.
\tag{9}
\]

Forward/backward recurrences prove conservation, unit mass across every rank,
and the comparison formula

\[
\Pr(x\prec y)=\sum_{I:\ x,y\notin I,\ x\text{ available}}\mathsf F(I,x).
\tag{10}
\]

Every extension crosses exactly one adding-\(x\) edge; (10) selects those
crossings before \(y\). This is a reusable exact compiler, with exponential
finite-state cost; no new efficient counting bound is claimed.

For any retained set \(A\) of ideals, `rank_transport` returns **both** its
mass \(q_k\) under (8) and the normalized conditional law. The full
unconditioned prefix law transfers unchanged iff \(q_k=1\). If \(q_k=0\), no conditional
law is fabricated. If edges are discarded, their lost flow must likewise be
retained as an obligation; renormalization is not conservation.

For the fence \(F_5\), there are 16 extensions. Keeping only the eleven
golden ideals gives retained rank masses

\[
(q_0,\ldots,q_5)=(1,1,1,1,5/8,0).
\tag{11}
\]

All predicates of the prefix ideal through rank three therefore transport
exactly after reweighting. At rank four, the missing ideal in (3) has mass
\(6/16=3/8\); at rank five the whole terminal mass is absent. Adding back
the two marked defect states restores the full source law. For example,
the actual fence comparison \(\Pr(0\prec4)=1/2\) is computed by (10).
Uniform ideals, conditional ideals, and uniform extension prefixes are never
silently substituted for one another.

## 6. The affine group is a separate symmetry, not an order

In the eleven-state chart,
\(H=\langle4\rangle=\langle3\rangle=\{1,3,4,5,9\}\), and
\(G=\{x\mapsto ax+b:a\in H,b\in\mathbb F_{11}\}\) has order 55.
This is the finite affine quotient of \(BS(1,3)\) inherited from THM-4592.
Translations do not commute with multiplication by four: the full group is
a cell-label action, not a group of autonomous Collatz conjugacies.

**Obstruction.** A finite poset with a transitive automorphism group must be
an antichain. All principal ideals would have equal size, whereas \(x<y\)
forces strict inclusion of principal ideals. Thus no nontrivial poset on
these eleven labels preserves the full affine action. In contrast, that
action preserves the Paley tournament, which contains
\(0\to1\to2\to0\). A directed cyclic comparison cannot be read as a
partial-order relation. Multiplication by four alone is not transitive on
all eleven labels, so the obstruction is deliberately for the full group.

The Golay theme supplies no additional map here: the identities involving
11 coordinates and 55 pairs do not specify a code, its Hamming metric,
or an extension measure. A code-theoretic transfer would have to retain
those objects explicitly. The eleven-squares torus also retains its physical
marks and clearances; the accepted paper is not a statement about all square
packings. See the inherited
[routing bridge](eleven_squares_routing_bridges_20261007.md).

## 7. Executable contract and finite universe

The owned [script](../../04-computation/experiments/golden_order_bridge_20261014j.py)
uses exact integers and rational weights. It rejects Boolean/float aliases,
invalid posets, unauthenticated count tables, wrong source masks, internal
`11`, inconsistent cut data, and zero activity. `compile_cut` deliberately retains
both closed and defective path states; only seam zero is a cyclic state.

Controls cover every bit word on cycles of lengths 3–13; two exact activity
laws at each length; all 408 distinct naturally labelled posets of orders
0–5 generated by forward-edge subsets; literal permutations as an independent
extension/flow/comparison check; the thin-order swap codec through order 11
(literal permutation equality through order 7); and all eleven rational
sources and 55 finite affine maps. Exhaustive coordinate gauges verify the
odd/even obstruction in this finite range; its all-length proof is section 2.
No general odd-length golden-reader bijection or universal poset conjecture
is inferred from these controls.

Reproduce:

```text
python -B 04-computation/experiments/golden_order_bridge_20261014j.py
python -B -O 04-computation/experiments/golden_order_bridge_20261014j.py
```

The [saved JSON](golden_order_bridge_20261014j.json) records exact counts,
weights, scope and a missing selected-family extension. Normal and optimized
runs agree on **63,095 checks**, LF SHA256
`252a5f13fd96758e5cf2db1196ea02615e1379394e8d50f5f7367d0e055a237b`.
Independent proof/API review and optimized-vs-saved replay **PASS**, including
the final all-length Gold-table delta and its exact source/type guards.
