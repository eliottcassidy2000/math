# The 67/69 reversal is an exposed domino contrast

Let \(A=(2,12,17,37,75,23,47)\) and
\(B=(2,12,17,37,43,71,31)\) be outgoing-row masks on vertices 0 through 6.
The reversed triangle is 4,5,6; in \(A\) it is \(4\to6\to5\to4\).
Simultaneously swapping \(2,3\) and \(4,5\) preserves every contact other
than reversing the triangle and the edge \(2,3\). Thus \(B\), after this
relabeling, differs from \(A\) only in the latter edge; this remains
true for all extensions below.

## Domino calculation

Paths not using the changed edge cancel. Contract the used oriented pair
\(3,2\) in the old arm or \(2,3\) in the new arm to a domino \(D\).
Write \(a=4,b=5,c=6\); their edges are fixed. For old \(z=a\) and new
\(z=b\), the only contacts from the three to \(D\) are \(c\to D,z\to D\).
The only contacts from \(D\) to the three are \(D\to z\). The resulting
two opposite incidences with \(z\) are intentional: they represent the
separate entry and exit tests at the two domino endpoints.

Allow an arbitrary nonempty painted gate tournament \(G\). All of
\(a,b,c\) beat every gate; each gate either beats both 2 and 3 or loses
to both, independently. The old gate group is \(0\to1\) with contact bits
zero, one, respectively. Hold \(G\) and its contacts fixed between arms;
gate vertices need not themselves form a homogeneous module. Since no
gate beats \(a,b,c\), a nonterminal gate segment must immediately precede
\(D\). Such a segment can enter \(D\) precisely when its final gate has
bit one; it is uniformly accessible from any of \(a,b,c\). A terminal
gate segment after \(a,b,c\) is unrestricted. A terminal gate segment
after \(D\) must start at a bit-zero gate. Counts for any fixed partition
of the gate vertex set between these segments therefore depend on the
type of slots, not on which triangle vertex precedes them.

The exterior comparison needs four numbers:

* With \(D\) last, the three triangle Hamilton paths offer three possible
  orders if the edge entering \(D\) is ignored, and two if it is tested.
  Both counts are equal in the two arms.
* With \(D\) nonterminal and no test of its incoming edge, there is one
  order for each of its three possible positions. Its immediate successor
  must be \(z\); the remaining vertices have their unique required order
  in the directed triangle. This count is three in either arm.
* With \(D\) nonterminal and its incoming edge tested, the \(D\)-first
  order exists in both arms. For either nonfirst position, the predecessor
  of \(D\) is the triangle predecessor of \(z\). Since \(z\) itself is
  already the successor of \(D\), this incoming arc exists if and only if
  that predecessor is \(c\). This is true precisely in the new arm.
  Therefore the corresponding counts are one and three, respectively.

Every partition involving use of a pre-\(D\) gate slot cancels, as does
every terminal-\(D\) partition. What remains is two additional orders,
both of which put all gates in a single unrestricted terminal path.
Therefore
\[
H(B_G)-H(A_G)=2H(G).
\]
This holds for arbitrary internal contacts of \(G\); in particular every
binary word on a transitive gate group preserves the gap exactly two.
If vertices 0 and 1 are inflated independently by homogeneous tournaments
\(X,Y\), then \(G=X\to Y\), so the gap is \(2H(X)H(Y)\). This explains
directly why all multiple-visit coefficients cancel in those two slots.
For arbitrary nonempty homogeneous occupants it is exactly two precisely
when both are transitive.

## A new independent clock that preserves the exact contrast

Adjoin a nonempty painted source bank \(W\) that beats \(G,a,b,c\).
Each vertex of \(W\) independently either beats both 2 and 3 or loses to both.
It may have arbitrary internal structure. In a domino path it may occur
initially and/or directly after \(D\); \(W\) and \(G\) contacts never allow it
in any other interior slot. Its visits are independent cover choices inside \(W\).

If \(W\) occurs initially but not after \(D\), the comparison above is
unchanged: it can change only the legality of the \(D\)-first order,
depending on the last source-bank vertex or a pre-\(D\) gate segment.
That exterior order has equal counts in the two arms.
If \(W\) occurs immediately after \(D\), the outgoing exterior edge
from \(D\) need no longer be tested. For either internal position of
\(D\), fixing its predecessor determines exactly one order (the other
two triangle vertices have their unique order). There are therefore three
orders before testing the incoming edge and two after testing it, in
both arms. The first and last positions also contribute equally.
This remains true whether or not there is an initial \(W\) visit.
Thus every post-\(D\) \(W\) term cancels, leaving two additional orders
with \(W\) wholly at the beginning and \(G\) wholly at the end:
\[
H(B_{G,W})-H(A_{G,W})=2H(G)H(W).
\]

For the original two singleton gates and \(W=TT_n\), the lower count is
\[
277+351 C_2(TT_n)+124 C_3(TT_n)=124\,3^n-21\,2^n-53
\]
and the upper count is two higher. The baseline formula uses the checked
local path-cover tensor, independently replayed from literal adjacency
below. The discrepancy formula does not depend on this finite tensor
fitting. This clock creates \(277/279,979/981,3127/3129,\ldots\) but does
not locate a clock at every absolute target.

## Why inflating the original c=6 clock is different

Inflate the two original gates homogeneously by \(X,Y\) and replace
\(c=6\) by an arbitrary nonempty tournament \(Z\). The module \(Z\) may
be visited at most twice: a second visit requires the single outside
vertex \(a\) between the visits. The single-visit lane contributes the
original \(2 H(X)H(Y)H(Z)\).

Write \(Z_1,Z_2\) for the two runs in order. When there are two, there is
one common exterior order ending at \(D\), namely \(Z_1,b,a,Z_2,D\).
For nonterminal \(D\), the one normal old order is
\(Z_1,D,a,Z_2,b\); the one normal new order is \(Z_1,D,b,a,Z_2\).
When a pre-\(D\) gate segment makes the test of the entry arc unnecessary,
there is one extra old exterior order \(Z_1,b,D,a,Z_2\), and no extra new
one. These lists exhaust the possibilities because \(a\) must be
immediately before \(Z_2\), \(D\) can only leave to its active neighbor,
and \(Z_i\) can only leave to \(b\) or \(D\) among exterior types.

The pre-\(D\) gate weight is \(y(1+x)\) for visits to \(X,Y\); the final
gate weight after nonterminal \(D\) is \((1+x)(1+y)\).
Requiring both gate types to be present, the two-\(Z\) signed visit
polynomial is therefore \(-z^2(2x+x^2)(y+y^2)\).
Writing \(C_2(U)\) for the ordered two-path cover count, this gives
\[
\Delta
=2H(X)H(Y)H(Z)-C_2(Z)(2H(X)+C_2(X))(H(Y)+C_2(Y)).
\]
As a consequence, exact gap magnitude two is possible only when
\(Z\) is singleton and both \(X,Y\) are transitive (sign \(+2\)), or when
\(X,Y\) are singletons and \(Z=TT_2\) (sign \(-2\)).
Indeed \(C_2(U)\ge (|U|-1)H(U)\) by splitting full Hamilton paths at
every index; the maps at a fixed split index are injective and indices
are disjoint. For \(|Z|\ge2\), \(C_2(Z)-H(Z)\ge1\), with equality only
at size two (size three separately gives \(C_2=6,H=1\) or 3).
The other terms in the negative summand enforce strictness if either
gate has more than one vertex.

## Replay

Run "python triangle_6769_controls.py". It computes literal tournament
counts by subset DP, checks arbitrary small painted gate configurations,
the new clock, and the exact third-old-module obstruction. It also
derives the coefficients \(277,351,124\) via exterior permutations
and ordered run counts, without using the predicted terms.

## Two programmable tapes still do not cover the odds

Take both painted banks transitive, and delete terminal zero bits of the
sink-gate tape \(G\), since those vertices are universal sinks; delete
initial one bits of the source tape \(W\), universal sources. Write \(n,m\)
for the normalized lengths. A single tester odd cycle using the final
one of \(G\) or the first zero of \(W\), and an arbitrary even subset of
the remaining tape, gives \(H\ge2^{n-1},2^{m-1}\) for a nonempty tape.
For the \(G\) cycle use \(2\to4\to G\cdots\to2\), and for the \(W\) cycle
use \(2\to W\cdots\to5\to2\). Thus \(n,m=O(\log X)\) if \(H\le X\).

When both banks have at least two active contacts, fix the first two
zero endpoints of \(W\) and last two one endpoints of \(G\). Use two
disjoint cycles of form
\[
\text{tester}\to W\text{-start}(0)\to\cdots\to
 \text{triangle vertex}\to\cdots\to G\text{-end}(1)\to\text{tester},
\]
using testers 2,3 and triangle vertices 4,5 respectively. The four
mandatory vertices per cycle are distinct. Ignore the fixed middle
runs between the two selected endpoints in each bank. Write \(S\) for
the number of tape positions after the second zero in \(W\) or before
the penultimate one in \(G\). Each of those \(S\) positions is assigned
independently to either cycle or neither. Requiring an odd number of
extras in both cycles gives exactly
\((3^S-2+(-1)^S)/4\ge3^S/9\) choices for \(S\ge2\).
Every resulting pair of odd cycles has weight four in the odd-cycle
packing formula. Hence \(S\le \log_3(9X/4)\).
With total \(S\) there are at most \((S+1)2^S\) words on the two variable
parts, and at most \(O((\log X)^2)\) middle-run length choices.
If one bank is empty or has one active contact, there is at most one
such bank per length; two disjoint odd cycles on the other bank give
the same \(3^S\) bound with even instead of odd extras. Consequently
the total number of distinct output pairs with lower height at most
\(X\) is \(O(X^{\log_3 2}(\log X)^3)\). This bound concerns only the
fixed exterior triangle, not arbitrary changed exterior cores.

## An exterior criterion for changing the baseline

Let \(C\) be any tournament of exterior vertices, split into disjoint
sets \(I,O,\{a,b\}\). Use two hubs with common incoming neighbors \(I\)
and common outgoing neighbors \(O\), and opposite contacts on \(a,b\)
so that the directed domino in one arm has both an incoming and
outgoing contact at \(a\), and in the other arm at \(b\).
This is implemented by reversing the mutual hub edge.
Write \(s(u),e(u)\) for whole-\(C\) start and end path counts;
\(P(u)\) for ordered complementary path pairs whose second starts at
\(u\); \(Q(u)\) for pairs whose first ends at \(u\), allowing the other
path empty in either case; and \(M(u,v)\) for pairs with both endpoints
specified, first ending \(u\) and second starting \(v\).

Suppose these four local metrics match between \(a,b\):
\[
(s(a),e(a),P(a),Q(a))=(s(b),e(b),P(b),Q(b)).
\]
The same two painted banks now have gap \(H(G)H(W)\) times
\[
\sum_{u\in I\cup\{b\},v\in O\cup\{b\}}M(u,v)
-\sum_{u\in I\cup\{a\},v\in O\cup\{a\}}M(u,v).
\]
Indeed a pre-domino bank segment erases the entry test and leaves
\(P\); a post-domino segment erases the exit test and leaves \(Q\).
Terminal orders are controlled by \(e\), initial ones by \(s\).
With both tests erased there is no difference at all. Hence only
the fully tested internal domino placements remain.
That is the displayed signed two-path-cover difference.

For example, the cyclic order-five tournament \(C=(6,12,24,17,3)\)
has the four metrics constant by translation. With \(a=0,b=1\),
\(I=\varnothing\), \(O=\{2,3,4\}\), its signed sum is two.
The resulting bare arms are \((38,76,24,17,3,94,29)\) and
\((38,76,24,17,3,30,61)\), with literal counts 89 and 91.
The equality of the \(Q\) rows gives the signed sum alternatively
as the skew local entry \(M(a,b)-M(b,a)=2\).

Add arbitrary ordinal tails \(U\to C\to V\), and independently paint
each of their vertices as belonging to \(I\) or \(O\). The cyclic
automorphism on \(C\) fixes the tails, so the contributions from their
single endpoints cancel by symmetry. Complementary paths with both
specified endpoints in \(C\) place all \(U\) on the first path and
all \(V\) on the second, so their signed sum is multiplied by
\(H(U)H(V)\). Consequently this gives four independently specified transitive tape segments with exact gap two.
The adjacent ordinal segments literally concatenate with the source
and sink tapes, so this is not a larger graph family. See
[PAINTED_DOMINO](PAINTED_DOMINO.md) for the general two-bank formula
and explicit exponent upper bound.
