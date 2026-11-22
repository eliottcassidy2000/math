# Reverse shells: exactly what the path-count probes can see

**2026-11-22 · PROVED (structural and observation results); Chan–Pak OPEN.**
A two-vertex reversal construction gives an exactly observable quotient of
the Rédei–Berge function. Adding vertex-deletion traces gives a strictly richer
two-variable quotient, and marking directed cycles recovers the full function.
None of the decoding operations makes an arbitrary decoded integer a path count.

## 1. The shell operator

Write \(F_T\) for the [Grinberg–Stanley Rédei–Berge symmetric function](https://arxiv.org/abs/2307.05569),
and \(Z_T(r)=F_T(1^r)\), with \(H(\varnothing)=1\). In coloring
coordinates, the weight of a color class is the Hamiltonian-path count of
its induced tournament. In power sums, the tournament formula is
\[
F_T=\sum_{\mathcal C}2^{|\mathcal C|}\,p_1^{u(\mathcal C)}
                    \prod_{C\in\mathcal C}p_{|C|},
\]
where \(\mathcal C\) runs over disjoint collections of directed odd cycles
of length at least three, and \(u\) counts uncovered vertices. This also
follows by expanding the monochromatic packing formula for each color class.

Let \(S(T)\) adjoin \(x,y\) with \(x\to T\to y\to x\). For a nonempty host
it is strong, and any supplied anti-involution extends by swapping \(x,y\).
Splitting or joining their assigned color gives
\[
 Z_{S(T)}(r)=(LZ_T)(r),\qquad
 Lf(r)=r^2 f(r)+r[f(r+1)-f(r-1)].
\]
The scalar transfer is **inherited** from the even binary family proof
in [Signed fold, path-slot transfer](SIGNED-FOLD.md#4-total-count-in-the-even-binary-family-a-central-difference);
the polynomial statement applies to every host. To retain the color
weights, if \(x,y\) have different colors they are a source and sink in
their respective classes. If their color agrees, add the path orders
using the arc \(y\to x\) by replacing their shared old color with two
distinguishable copies (old prefix and suffix). Exclude the case in
which both copies are empty, which is already counted in the baseline.
For the common color of weight \(x_i\), this is adding one copy of
\(x_i\) minus removing its sole copy. The factor for the new vertices
is \(x_i^2\). The resulting universal operator is
\[
\mathcal S F
=p_1^2 F+
 2\sum_{\substack{\alpha\ge0\\|\alpha|\ {\rm odd}}}
 p_{2+\sum_{j\ {\rm odd}} j\alpha_j}\,
 \frac{\partial^\alpha F}{\alpha!}.
\]
Indeed adding or removing color \(i\) translates \(p_j\) by
\(\pm x_i^j\); Taylor's formula and the sum over \(i\) give the identity.
Only finitely many shifts occur for a particular polynomial. In packing
coordinates the correction chooses an odd number of cycles or singletons
and merges them into one odd part, adding two vertices. This is a
differential operator with explicit positive terms, before specializing.

For separated interior blocks, \(F\) and \(Z\) multiply under ordinal
join, since any directed cycle lies in one block. In particular
\[
 H(S(T_1\Rightarrow\cdots\Rightarrow T_s))
    =\prod_i H(T_i)+\prod_i Z_{T_i}(2)
\]
for nonempty blocks. A transitive \(m\)-block has \(Z(2)=2^m\);
reversing its extreme arc for \(m\ge3\) gives \(Z(2)=2^m+2\cdot3^{m-2}-2\).
The unmarked outer arc flip alone changes \(H\) by \(Z_T(2)\ge 2^{|T|}\),
since each of the \(2^{|T|}\) splits carries at least one path pair.
Thus it is not a general unit-step seam.

## 2. Exactly observable by unmarked shells and joins

For an \(n>0\) vertex host, the first \(\lceil n/2\rceil\) numbers
\(H(S^k(T))\), starting at \(k=0\), determine \(Z_T\). In the operator
\(L^k\) evaluated at 1, the largest evaluation point is \(k+1\), attained
only by \(k\) upward shifts, with coefficient \(k!\); any evaluation at
zero is harmless (\(Z_T(0)=0\)). Therefore these counts first recover
\(Z_T(1),\ldots,Z_T(\lceil n/2\rceil)\). Its exponents are
\(n,n-2,\ldots\), stopping at 1 or 2, by the packing formula. Divide by
\(r\) or \(r^2\) according to parity and interpolate in \(r^2\).
Conversely \(Z_T\) predicts every count obtained by combining \(S\) and
ordinal joins with fixed graphs, even if the host is used more than once.
Thus equality of \(Z\) is exactly equality under all those observations.

This invisibility extends from ordinal joins to arbitrary fixed-quotient
lexicographic substitutions, provided the observed host is used as a
**factor**, not as the quotient. Indeed \(Z_T(r)=\sum_d \binom{r}{d}P_d(T)\),
where \(P_d\) is the number of ordered nonempty \(d\)-path covers. To
compute the \(q\)-slot profile of a substituted tournament, choose a
quotient walk for each of the \(q\) slots, allowing empty walks. The
total visits to module \(i\), ordered across all slots, are filled by
\(P_{d_i}\) choices. Thus the new profile depends only on the module
profiles and the fixed quotient. This is the multi-slot version of
[the module-cover law](EULER-GLUING.md#4-the-general-module-law-euler-flows-paired-against-ordered-path-covers).

The dimension bound is real, not an artefact of interpolation. Over
\(\mathbb Q\), the order-\(n\) profiles span the entire indicated parity
space: use a transitive tournament and the tournaments with one reversed
outer arc inside a transitive block of order \(2j+1\), padding by an ordinal
join of singletons as necessary. All odd cycles in the block use the
reversed arc; there is one full cycle of length \(2j+1\), no longer one,
and no disjoint cycles. The least power is consequently \(r^{n-2j}\)
with coefficient 2. These profiles form a triangular basis. The stated
number of evaluations is optimal for linear reconstruction on this span;
no information-theoretic lower bound for nonlinear decoding on the finite
set of actual graphs is claimed.

## 3. The vertex-deletion trace refines this quotient exactly

Put \(B_T(x,y)=F_T(p_1=x,\ p_3=p_5=p_7=\cdots=y)\). In the packing
formula, \(\partial_{p_1} F_T=\sum_v F_{T-v}\): differentiate by marking
one uncovered singleton and deleting it. Higher derivatives similarly
delete ordered tuples of distinct vertices. Thus \(B_T\) is closed under
vertex-deletion traces, under ordinal product, and under the shell:
\[
 B_{S(T)}(x,y)
 = x^2B_T(x,y)+y[B_T(x+1,y+1)-B_T(x-1,y-1)].
\]
In the multivariate operator every correction has merged index at least
three and therefore specializes to \(y\). Summing the divided derivatives
of total odd order along the all-ones direction proves the displayed law.

Conversely, for every \(d\), take the trace of the shell observations over
all ordered \(d\)-vertex deletions of \(T\). The same triangular
interpolation recovers \((\partial_x^d B_T)(r,r)\) (the empty host is
known). Taylor expansion in \(x\) about \(x=r\) recovers \(B_T(x,r)\) for
all \(r\), then \(B_T\). So the joint data of all such traced probes
are *exactly* the two-variable quotient; they remember the number of
singletons and of nontrivial cycles, but not the individual longer lengths.

More generally, for odd \(\ell\ge3\),
\[
 \partial_{p_\ell}F_T=2\sum_{C:\,|C|=\ell} F_{T-V(C)}.
\]
The derivative marks one particular directed cycle in a collection.
Iterating these identities counts ordered choices of marked disjoint
cycles, with optional marked singletons. Evaluating all derivatives at
\(p_j=1\) reads a finite Taylor expansion and recovers the full \(F_T\).
This is a way to restore the lost grading using **traces of marked
deletions**, not a claim that a single unmarked tournament count does so.

## 4. Small self-converse witness and physical probes

Using bit \(j\) to indicate an arc toward \(j\), take the two six-vertex
tournaments
\[
 U=(58,4,1,6,14,30),\quad V=(38,8,2,5,15,30).
\]
Anti-involutions are \((2,5,0,4,3,1)\) and \((3,5,4,0,2,1)\), as image
arrays in labels \(0,\ldots,5\). Both are strong. Each has four directed
triangles. The counts of 5-cycles and of disjoint triangle pairs are
respectively \((4,0)\) and \((2,1)\). Consequently
\[
 Z_U(r)=Z_V(r)=r^6+8r^4+8r^2,\qquad
 F_U-F_V=4(p_5p_1-p_3^2).
\]
The common first four unmarked shell counts are \(17,241,4001,80929\).
But their summed vertex-deleted counts differ by 4, since
\((\partial_{p_1}F_U-\partial_{p_1}F_V)(1)=4\).

One can also separate the pair **without deleting anything**. Adjoin a
new fixed center, choosing one vertex of each antipodal pair to be beaten
by the center and having the other vertex beat the center. This preserves
both self-converseness and strength. For the eight possible sign choices
the total counts are the following disjoint sets:
\[
 U:\ \{33,49,51,67,93,105,111,123\},\qquad
 V:\ \{27,41,65,69,75,87,115,129\}.
\]
The smallest-labelled member of each mate pair is one common rule
giving 49 and 41, respectively. Disjointness of the whole sets avoids
any dependence on that convention.

There is also a mark-free operation that changes the role of each
host. [THM-3235, §§6.3–6.4](../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md) already separates same-profile quotient twins by uniform inflation; here the witnesses remain strong and self-converse, and the separation is strict at every subsequent shell depth. Specifically, use it as the **quotient** and replace each vertex by a transitive
pair. The result is again strong and self-converse, reversing each pair
while applying the old anti-involution. Write \(U',V'\) for the two
inflations. Exact evaluation gives \(H(U')=27513\), \(H(V')=23361\), and
\[
 Z_{U'}(r)-Z_{V'}(r)=336r^6+2760r^4+1056r^2.
\]
The difference has nonnegative coefficients and is strictly positive on
positive inputs. Since \(L\) is coefficientwise nonnegative on monomials
(expand the central difference by the binomial theorem), the separation
persists after every number of further shells. The same two hosts are
therefore interchangeable as modules but are distinguishable as routing
graphs even with identical two-vertex factors.

A second named probe uses the unique vertex \(v\) of outdegree two.
In the shell reverse the two arcs \(x\to v\) and \(\sigma(v)\to y\),
where \(\sigma\) is the host's anti-involution. The changes are exchanged
by the involution, and the host remains strong: \(x\) still reaches the
old core, the core still reaches \(y\), and \(y\to x\). The two ordinary
counts are 333 and 325. Thus a small canonical asymmetric contact also
breaks the equality.

The [one-center endpoint law](CENTER-PROBE.md) proves the exact
all-order quadratic response and a graded orientation-neutral trace.

This is a concrete separation between scalar shells, marked vertex
traces, and a reflection-preserving change of the interface. The next
arithmetic question is whether changing marked contact can select
inverse branches with disjoint **absolute** images, not merely separate
two hosts or recover a grading by linear combinations.

Finite checks on the displayed new instances: run the accompanying
`check_shell_information.py` from this directory.
