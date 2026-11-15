# Bounded repairs over bounded counts, and a sharp crossing normal form

**Scope.** This is an all-order obstruction and a constructive exact
inverse for a two-edit class, not an exclusion from the full tournament
spectrum. Within this strategy for Chan–Pak, one must admit unboundedly
many edits or bases of unbounded initial count.

## The bounded-edit theorem

Fix integers \(B\ge1\) and \(q\ge1\). Allow an arbitrary initial
tournament with \(H\le B\), on any number of vertices, then reverse at
most \(q\) arcs. The number of **distinct final values** at most \(X\) is
\[
O_{B,q}\big((1+\log X)^{\lfloor\log_3 B\rfloor+2q-1}\big),\quad X\ge1.
\]
There is no bound on order before or after editing. Width here is
actual edit budget, not distance between endpoints.

**Proof.** Use the [Grinberg–Stanley odd-cycle expansion](https://arxiv.org/html/2307.05569#S1.SS8.SSS3):
each collection of \(r\) vertex-disjoint directed odd cycles of length
at least three contributes \(2^r\) to \(H\), including the empty collection
once. Each vertex of a nonsingleton strong tournament lies on some
directed triangle: for \(v\), both neighbor sets are nonempty, and if no
outneighbor beats any inneighbor then the outneighbors have no arc to
their complement, contrary to strongness. Thus a strong component of
order \(m\ge3\) has at least \(m/3\) distinct triangles, and
\(m\le 3(H-1)/2\). Across the initial ordinal chain of SCCs the Hamilton
count is the product of the component counts, because paths exhaust
each component in order. Their nonsingleton counts \(h_i\) satisfy
\(h_i\ge3\), so there are at most \(k=\lfloor\log_3 B\rfloor\) nonsingleton
SCCs and at most \(3\sum_i(h_i-1)/2\le 3(B-1)/2\) vertices in them.

Keep those SCCs and mark every distinct singleton SCC incident to an
edited arc. Treat each whole old nonsingleton SCC and each marked
singleton as one ordered anchor; there are at most \(k+2q\) anchors.
All other vertices form transitive consecutive intervals before,
between and after the anchors. Each interval is a module even after
the edit: each outside vertex sees it uniformly. The fixed internal
kernels, edits between or within them, order of anchors, and the set of
nonempty intervals have finitely many possibilities depending only
on \(B,q\). The two exterior intervals contribute no cycle and do not
affect \(H\), since they are pure forward prefixes/suffixes.

Shrink every positive interval to one representative. If an interval's
representative belongs to no directed cycle, that interval can be deleted
from the actual tournament without changing \(H\). Indeed, any cycle
touching that interval must go outside because it is transitive.
The trip from leaving to returning projects to a closed directed walk
through the representative, and shortening repetitions produces a
cycle through it. Status therefore depends only on the finite kernel
and the set of positive intervals.

Conversely, for a representative lying on a directed cycle, choose a
simple such cycle, lifting fixed representatives elsewhere. Replace
the chosen representative by any nonempty subsequence of its original
interval. All arcs still exist, in forward internal order. If its
length is \(t\), at least \(2^{t-1}-1\) subsets make the new cycle odd
(half the subsets have each parity, with at most the empty one removed).
The empty packing plus these one-cycle contributions gives
\(H\ge2^t-1\). Thus every active interval has \(t\le\log_2(X+1)\).
Only the \(k+2q-1\) internal intervals can be active, giving the bound.
If there are no anchors then the base is transitive, no arcs change,
and \(H=1\).

### Exact finite kernel evaluation

For a fixed edited kernel with \(a\) retained individual vertices,
the dependence on interval lengths \(t_i\) is a finite integer linear
combination of functions \(\prod_i j_i^{t_i}\), \(0\le j_i\le a+1\),
with \(0^0=1\). Choose the order of the retained vertices on a Hamilton
path. There are \(a+1\) slots, including before and after. Within
a slot, interval types occur in original order, and vertices of
each interval in ascending order. Choose which interval occupies
which slots, reject if the compressed path uses a missing arc, and
count actual assignments into the chosen slots. If interval \(i\)
occupies \(s_i\) slots, its factor is the number of surjections of its
\(t_i\) vertices to those \(s_i\) labeled slots:
\(\sum_{j=0}^{s_i}(-1)^{s_i-j}\binom{s_i}{j}j^{t_i}\).
The path recovers its retained order and slot occupancy. This is an
exact inverse for fixed budgets: enumerate the logarithmic active
ranges, and set each positive inactive interval to length one.

## Two reversals from a transitive order: the five shapes

Put the anchors in original order. If there are four distinct endpoints,
write \(a<b<c<d\), with gaps \(x,y,z\ge0\). For three distinct endpoints
use \(a<b<c\) and \(x,y\ge0\). Let \(f(t)=2^t\) for \(t\ge1\) and \(f(0)=0\);
write \([P]\) for the indicator of \(P\). Every possible pair of reversed
arcs is in the table, up to converse. Exterior padding is irrelevant;
a single reversal across gap \(t\) has \(H=1+f(t)\).

| Shape | Reversed arcs | Exact \(H\) |
| --- | --- | --- |
| Disjoint | \(b\to a,d\to c\) | \((1+f(x))(1+f(z))\) |
| Nested | \(d\to a,c\to b\) | \(1+f(y)+2^{x+z+1}(2^{y+1}+3^y-1)\) |
| Crossing | \(c\to a,d\to b\) | \(1+2^{x+y+1}+2^{y+z+1}+2^{x+z+1}3^y-[x=0]2^z-[z=0]2^x\) |
| Same head | \(b\to a,c\to a\) | \(1+f(x)+2^{x+y+1}-2^y-[y=0]\) |
| Chain | \(b\to a,c\to b\) | \(1+f(x)+f(y)+2^{x+y}+[x+y=0]\) |

**Derivation.** All cycles use a reversed arc, and between reversals
they are increasing. A single reversed arc requires an odd number of
internal vertices. Its weighted contribution is \(f(t)\) if all \(t\)
interior vertices are free.

* Disjoint intervals lie in separate SCC regions, so the counts multiply.
* With common head (or, by converse, common tail) no pack can use both
  arcs. For the outer arc remove increasing routes using the missing
  \(a\to b\), i.e. containing \(b\) and no vertices in the first gap.
  Their remaining \(y\)-subset is even, contributing \(2^y+[y=0]\)
  after weighting.
* With a shared chain middle, one cycle using both arcs must contain
  the descending run \(c,b,a\); its return path can take any even
  subset of the two gaps, giving \(2^{x+y}+[x+y=0]\).
* For crossing and nesting, a cycle using both reversed arcs has two
  increasing runs between its four mandatory endpoints. In the central
  gap every vertex is unused or in either run; in each outer gap it is
  unused or in its only available run. The total assigned is odd.
  Its weighted contribution is
  \(2^{x+z}3^y-[x=z=0](-1)^y\) by the even/odd character.
  If two disjoint odd cycles use one reversal each, the central gap
  again has three choices; the \(x,z\) gaps belong to distinct cycles
  in the crossing shape and both to the outer cycle in the nested
  shape. The two odd parity filters, after weighting by four, give
  respectively
  \[
  2^{x+z}3^y-[x=0]2^z-[z=0]2^x+[x=z=0](-1)^y
  \]
  and
  \[
  2^{x+z}3^y-[x=z=0]-2^{x+z}+[x=z=0](-1)^y.
  \]
  These ascending routes are allowed: no forbidden forward arc has
  both endpoints in the same chosen run, except the unused forward
  copy of the cycle's own reversed arc in a possible 2-cycle.
  The odd filters remove any such 2-cycles.
  The individual cycles in the crossing shape have free gap lengths
  \(x+y+1,y+z+1\). In the nested shape the inner has length \(y\),
  and the outer has \(x+y+z+2\) free interior vertices except that
  it must not traverse the missing \(b\to c\).
  This removes subsets including \(b,c\) and no central vertex, of odd
  cardinality overall. Their weighted contribution is
  \(2^{x+z}-[x=z=0]\). Summing proves the table.

## Prescribed-target inverse, not a residue test

In the crossing case \(x,z\ge1\), put \(s=3^y\), \(t=2^y\).
For a prescribed odd \(N\), the following is equivalent to a physical
witness in this chart:
\[
(s2^x+t)(s2^z+t)=t^2+s(N-1)/2,\qquad y\ge0,\ x,z\ge1.
\]
Each factor must have the displayed form: subtract \(t\), divide
by \(s\), and require a power of two at least two. No factorization
oracle is needed. Try powers for a first factor no larger than the
integer square root of the right side; then compute the other
factor by divisibility and check its address. Only \(O(\log N)\)
values of \(y\) are possible by the \(8\cdot3^y\) term, and only
\(O(\log N)\) powers per \(y\) are possible. For \(y=0\) the same
formula reads \(H=2(2^x+1)(2^z+1)-1\); \(y>0\) need not be a
Fermat product.

There is a direct inverse in the chamber \(1\le x<z<y\). The
smallest dyadic valuation occurs once, so \(\nu_2(H-1)=1+x+z\) and
\[
(H-1)/2^{1+x+z}=3^y+2^{y-z}+2^{y-x}.
\]
The right side is in \([3^y,3^{y+1})\). Recover \(y\) by exact powers
of three, subtract \(3^y\), and recover \(x,z\) from the ordered two
set-bit positions. Check the valuation and chamber inequalities.
This family is injective. Choosing all positive \(x<z<y\le c\log X\)
with small absolute \(c>0\) gives \(\Omega((\log X)^3)\) distinct counts
at most \(X\). With the bounded-edit theorem at \(B=1,q=2\), the full
at-most-two-edit spectrum thus has \(\Theta((\log X)^3)\) elements.
The four other shapes contribute \(O((\log X)^2)\) together by the
table and the active-length bound (nesting uses \(x+z\) only).
Every graph in the chamber has a Hamilton frame in the original
order, and its adjacency differs from that frame's transitive
completion at precisely two pairs.

**Mechanism.** The scale \(3^y\) is the number of choices per free
middle vertex between two simultaneously open runs in the odd-cycle
expansion. Nesting has the same ternary choice but makes both outer
gaps equivalent, losing one exponent. The converse preserves \(y\)
and exchanges \(x,z\). The repair question is which changing kernels
can carry arbitrarily many active gaps while keeping an exact decoder.
Adding only neutral padding and bounded edits to a bounded-count
source provably cannot suffice.

## Joint inflation bound: the necessary number of edits grows

The active intervals satisfy a stronger bound than the separate
logarithmic cutoffs:
\[
H\ \ge\ \prod_{\text{active }i}(2^{t_i}-1),
\qquad
\sum_{\text{active }i}(t_i-1)\ \le\ \log_2 H .
\tag{1}
\]
To prove it, pick a directed Hamilton cycle in each nonsingleton SCC
of the quotient. For completeness, in a strong tournament take a
longest directed cycle \(C\). A vertex with both kinds of adjacency
to \(C\) inserts between consecutive vertices at an incoming-to-outgoing
transition. If there is no mixed outside vertex, partition the outside
into \(W\) beating \(C\) and \(L\) losing to \(C\). Both are nonempty by
strongness, and some \(l\to w\) must cross from \(L\) to \(W\) since
otherwise \(L\) has no outgoing arc to its complement. Inserting
\(l,w\) together into a single cycle edge contradicts maximality.
This proves the cycle claim.

Each active representative on its SCC's cycle can be replaced
independently by any nonempty subsequence of its interval. For
a \(t\)-interval the number of even and odd nonempty subsets
differ by exactly one. Consequently the parity imbalance in the
product of these choice sets has absolute value one. The number of
odd lifted cycles is at least \((\prod_i(2^{t_i}-1)-1)/2\) for that
SCC. The empty pack plus these one-cycle terms gives the required
factor; the factors multiply over disjoint SCCs. When no interval
occurs in a SCC, keep the factor one by the empty packing. Finally
\(2^t-1\ge2^{t-1}\) for \(t\ge1\). This proves (1).

Now keep \(B\) fixed but allow the edit budget \(q\) to grow with \(X\).
Write \(L=\lfloor\log_2 X\rfloor\) and \(k=\lfloor\log_3 B\rfloor\).
There are at most
\[
O_B((q+1)^{k+1})\,[C_B(q+1)]^q\,
2^{2q+k}\binom{L+2q+k}{2q+k}
\tag{2}
\]
count presentations, for some fixed \(C_B\) (including a possible
constant count from the no-anchor case). To see this, the small old
nonsingleton SCCs have finitely many ordered internal structures.
Interleave those at most \(k\) anchors with at most \(2q\) marked
singletons, using at most \(O_B((q+1)^{k+1})\) orders including the
choice of number of singletons. The distinguished number of vertices
is at most \(2q+3(B-1)/2\); choosing up to \(q\) reversed pairs on
them has at most \([C_B(q+1)]^q\) possibilities by the binomial
bound (absorbing a factor \(q+1\) into \(C_B^q\)). Choose the positive
internal gaps, at most \(2^{2q+k}\) masks. Each mask fixes its list
of active intervals. If that list has size \(s\), the nonnegative
excesses \(t_i-1\) fit in the simplex of sum at most \(L\) by (1),
giving at most \(\binom{L+s}{s}\) choices. Set all positive inactive
gaps to one. This yields (2) without needing to bound the size of
the original tournament.

In particular if \(q\le c\log X/\log\log X\) for fixed \(0<c<1\)
(natural logs), (2) is \(O_B(X^{c+o(1)})\): the edit set contributes
\(q\log(q+1)+O_B(q)\), and the length simplex contributes only
\(O_B((\log X)(\log\log\log X)/\log\log X)\) to the logarithm. Thus all
but \(o(X)\) odd targets at most \(X\) have no witness in this edit
budget, even with arbitrary-sized baselines of bounded count.
For a single transitive baseline count, this translates to ordinary
minimum feedback-arc distance for any potential witness. This is
an exact support obstruction, not a proof that all remaining
targets have some tournament witness.

If \(B\) varies, the [Busch bound](https://doi.org/10.37236/1141)
gives \(a_0\le 3\log_3 B\) total vertices in old nonsingleton SCCs
(the triangle case is immediate). This imports the reduced-order
cutoff from [fixed-count tightness](../chanpak_crossing_scales_20261115/FIXED_COUNT_TIGHTNESS.md).
There are \(\exp(O((1+\log B)^2))\) ordered reduced kernels:
encode each on its ordered vertex set and mark SCC boundaries. Their
interleavings with the marked singletons contribute at most
\((2q+1)\,2^{k+2q}\). The edit-selection term can be bounded by
\(\exp(q\log(e(a_0+2q)^2/q)+O(\log(q+1)))\) when the possible-pair
count is at least \(q\), and by \(\exp(O(q))\) otherwise. The same
positive-gap mask and simplex bound apply. In particular the
\(X^{c+o(1)}\) conclusion still holds for
\(\log B=o(\sqrt{\log X})\) and
\(q\le c\log X/\log\log X\): substitute the maximal \(q\) in the
monotone budget bound and note that \(a_0=o(q)\) at that maximum.
So even subexponentially many baseline counts, in this precise
range, do not compensate for a too-small edit budget.

## Rooted consequence: full bounded-baseline \(H+2E\) is only dyadic rays

This is the rooted analogue of the earlier finite dyadic-ray description of unrooted \\(Q_2\\) at fixed count; it now gives a direct test for the outlet inverse. Let
\(E_u(W)\) and \(O_u(W)\) count, respectively, even- and odd-size
nonempty paths starting at \(u\), each weighted by the Hamilton
count of its complement. Adding an outlet \(v\) beating only \(u\)
gives \(N=H(W)+2E_u(W)\): the new cycles use \(v\to u\), and
an odd one corresponds to precisely such an even path in \(W\).
Other cycles are disjoint and counted on the complement.

Write \(W=P\to C\to S\), where \(C\) is the SCC containing \(u\).
Let \(Q(S)=\sum_{A\subseteq S}H(S[A])H(S[S\setminus A])\),
the number of ordered two-path covers, permitting empties.
Then
\[
N=
\begin{cases}
H(W)+2H(P)E_u(C), & S=\varnothing,\\
H(W)+H(P)\,[E_u(C)+O_u(C)]\,Q(S), & S\ne\varnothing .
\end{cases}
\tag{3}
\]
A rooted path in \(W\) cannot use \(P\), and splits into a nonempty
rooted path in \(C\) and a possibly empty path in \(S\). The suffix's
selected subsets have each size parity with equal weighted sum
\(Q(S)/2\) when \(S\ne\varnothing\). To see this, expand their two
Hamilton counts by odd cycle packs. Any packing of cycles in \(S\)
splits the vertices into blocks, one per cycle and one per unused
vertex, each of odd size. All assignments of those blocks to the
two parts carry the same weight. Changing the part of one fixed
block reverses the parity of the selected part, proving equality.
For an empty suffix, there is only an even empty part. This proves (3).

The quantity \(Q\) multiplies across an ordinal sum, because the
two slots separately concatenate along its components. Thus
\(Q(S)=2^t\prod_{D\text{ nonsingleton SCC in }S}Q(D)\) with \(t\)
the number of suffix singletons. If \(H(W)\le B\), the proof above
bounds the nontrivial components up to isomorphism and quantity,
and \(C\) is one of these or a marked singleton. Formula (3)
therefore describes the **entire** bounded-baseline rooted image
as finitely many constants and finitely many rays
\[
h+d\,2^t,\qquad h,d\text{ fixed positive integers},\quad t\ge t_0,
\]
where \(t_0=0\) for a nonempty nontrivial suffix and \(t_0=1\) for
a suffix constrained to be nonempty and all singleton.
All literal ray points exist by appending singleton SCCs to the
same suffix, and all such choices have unchanged baseline count.
The coefficients and finite exceptional constants can be obtained
from the finite SCC kernels without guessing the padded order.
Consequently the full bounded-baseline rooted image has at most
\(O_B(\log X)\) distinct values at most \(X\). Even changing the
remainder's size arbitrarily, but keeping its Hamilton count
bounded, cannot close the absolute-target gap using one outlet.
