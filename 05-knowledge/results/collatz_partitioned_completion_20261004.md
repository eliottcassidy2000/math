# Different constructions for different Collatz proof obligations

2026-10-04. **PROVED:** a composite construction supplies new guarded
smaller-child dependencies; a whole arithmetic cell escapes every fixed
pair of search-depth bounds; each expanding fixed-word pair has a finite
child-cut obligation. **FINITE-EXACT:** the component experiments.
**OPEN:** guaranteed completion for every positive integer.

The productive partition is by applicable proof mechanism, with the
original source and its final smaller dependency retained. The new result
is an infinite collection of actual common-future constructions, including
one arithmetic progression with a proved positive portion outside the
previous sixteen-debt and ternary-sibling bank. This does not assert home
certificates for all members of that progression: their children remain
separate obligations.

## The composition that improves coverage

The [complement-routing construction](collatz_complement_routing_20261004.md)
uses three operations in order: follow a growing source prefix, select a
ternary sibling connection, then apply a binary reset rule to the resulting
intermediate child. Test the final child against the original source.
Requiring the intermediate child itself to be smaller would discard valid
compositions.

A concrete family is

\[
\begin{aligned}
 n_t&=27727075633746555+79062194724345216t,\\
 h_t&=25270565367447551+72057594037927936t,
 \qquad t\in\mathbb Z_{\ge0}.
\end{aligned}
\]

For every parameter, the exact valuation words give

\[
 n_t\xrightarrow{(1,2,1)}J_t
 \xleftarrow{(1^{32},2,19)}h_t,
 \qquad 0<h_t<n_t.
\]

The child slope is \(2^{49}/3^{31}<1\). Here \(1^{32}\) means 32
successive valuation-one steps. Source and child are positive odd
integers, and all indicated valuations are exact. A supplied first-hit
certificate for \(h_t\) transports to the original \(n_t\).

The source's first four valuations are \((1,2,1,2)\); each corresponding
iterate exceeds the original source. Every source is divisible by nine,
so it has no odd predecessor of its own. The common-future construction
works by moving forward first.

**PROVED relative coverage.** At least \(77/78\) of this progression,
in relative natural density, lies outside the specified previous binary
and ternary bank. The binary words exclude its sixteen rows directly.
The first ten ternary addresses are disjoint, and the later ones have
total relative occupancy at most \(1/78\). The proof includes the
source-size cutoff needed to control the infinite tail. This fraction is
relative to this one progression, not to all integers or all currently
unresolved inputs; other repository rules have not been exhausted.

The construction extends to an arbitrary initial run of \(r\ge1\)
valuation-one steps ending at a valuation-two reset. Its sibling index
satisfies \(k\equiv1\pmod{3^{r+1}}\). The required word depth grows with
these parameters. A height bound makes testing applicability to any one
supplied source finite. These are unbounded word templates, rather than
a larger fixed-depth lookup table.

## The exact role of modulo 19

In the general construction, let \(\ell\) be the chosen inverse-run
length and \(h\) its final child. Then

\[
 h+1\equiv (2/3)^{\ell-r-1}(n+1)\pmod{19}.
\]

This follows from \(4^k\equiv4\pmod{19}\) for the allowed indices and
\((2/3)^3\equiv1\pmod{19}\). It is an exact three-phase law. Increasing
\(\ell\) by at most two makes \(\ell-r-1\) a multiple of three. The
corresponding ternary guard becomes finer by a factor at most nine. The
child does not increase, decreases strictly when depth increases, and
satisfies \(h\equiv n\pmod{19}\).

For the displayed family this selects the refinement

\[
\begin{aligned}
 n_t&=502100243979817851+711559752519106944t,\\
 h_t&=203384946486673407+288230376151711744t.
\end{aligned}
\]

All 19 source phases occur, and each pair preserves its phase. The
integer inequality proves descent; the congruence records information
that survives it. Membership in the same residue class does not prove
that the child satisfies the same complete construction guard.

## Which methods take responsibility for which sets

A dispatcher can give overlapping rules a priority and thereby partition
their union into disjoint domains. It must retain all valid alternatives
when actual child certificates are available: the numerically smallest
child need not be the child whose proof is already known.

| Mechanism | Applicable source condition | Resulting obligation |
|---|---|---|
| Ordinary halving and immediate odd descent | Even input; or odd input with first valuation at least two | A smaller actual iterate, stopping separately at ROOT1 |
| Small inverse words | Odd \(n\equiv2\pmod3\), or \(n\equiv4\pmod9\) | Respectively \((2n-1)/3\) or \((8n-5)/9\), each smaller with a checked route to \(n\) |
| Unbounded binary reset rule | Initial run of ones followed by a reset at least three | The half-child with a checked common future |
| Existing debt and sibling rules | Their exact binary or ternary source guards | A specified smaller child; sixteen debt rows and an unbounded sibling-address family are retained |
| New composed rule | Compatible reset-two prefix and ternary address, with the final size bound | A smaller final child even when the intermediate child is larger |
| Reused grounded suffix | A retained actual path reaches an already certified component | A completed first-hit certificate for the original request |

The residual domain consists of sources for which none of the selected
rules has supplied the required evidence. A direct descent rule on a
class does not itself prove home for every member: its child may lie in
another class. A complete global cover by strict smaller-child rules
would close this through induction. We have not proved that cover.

As a concrete boundary split, the new composite schema misses
\(n=2^{2j+1}-1\), \(j\ge1\), because its inverse-depth size requirement
is too large at these least binary representatives. The inherited
\(4\pmod9\) rule nevertheless handles every member \(2^p-1\) with
\(p\equiv5\pmod6\), via the smaller child \((8n-5)/9\). Other exponent
classes require other rules. A failure of one construction therefore
remains attached to that construction.

## Why a fixed-depth partition cannot be the whole proof

The [bounded-depth obstruction](collatz_partition_cover_20261004.md)
proves the following without assuming any source reaches the root.
For any nonnegative forward and inverse depth bounds \(R,S\), every
positive odd \(n>1\) satisfying

\[
 n\equiv-1\pmod{2^{R+1}},\qquad n\equiv0\pmod{3^S}
\]

has no smaller common-future source within those bounds:

\[
 U^r(n)\ne U^s(m)
 \quad(0<m<n,\ 0\le r\le R,\ 0\le s\le S).
\]

Here \(m\) ranges over positive odd integers.

Inverse valuations are unrestricted. Thus the obstruction does not come
from an artificial cap on individual exponents. This entire arithmetic
cell has relative density \(1/(2^R3^S)\) among odd integers. Finer residue
filters or more rules of the same bounded depths cannot cover it.

The proof keeps the formal source-zero intercept of a candidate inverse
word. Its required sign conflicts with integrality after the dyadic
denominator clears. Explicit unbounded reset templates enter subfamilies
of the obstructed cell after exceeding the forward bound. The result
therefore identifies a necessary unbounded coordinate, not a class of
nonconvergent integers.

## A second partition makes each fixed word pair finite

The [child-cut obligation analysis](child_cut_obligation_partition_20261004.md)
normalizes compatible source and child words as
\(n=n_0+Pt\), \(m=m_0+Qt\), with \(0<n_0<P\), \(0<m_0<Q\).
If the lift expands, \(Q>P\), the original child is smaller only for

\[
 0\le t\le
 \left\lfloor\frac{n_0-m_0-1}{Q-P}\right\rfloor.
\]

All child cuts whose slopes can repair the family collectively cover a
possibly empty upper parameter ray. Their complement in the displayed finite interval
is therefore an exact finite obligation for that particular word pair.
This replaces another indefinite search by a specified finite check.
There are still infinitely many word pairs.

A rational example with a completed route defeats a proposed argument
based only on positive values and record order. Its exact ternary
integrality guard excludes that example from the integer problem. The
integer guard is essential; the universal integer child-cut assertion
remains open.

## Inheritance and the remaining target

Closest mechanisms: exact common-future lifts and retained grounded
suffixes. Canonical hostiles: unpaid \(27\to41\to31\), expanding
\(233/231\), and the new whole-cell depth obstruction. Corrected near
miss: local descent was once incorrectly described as impossible on
every residue class; \(5\pmod8\) is an immediate counterexample, and the
historical row and mistakes ledger are repaired. Least-used sidecars:
intermediate versus final child, two independent word depths, and the
integer endpoint lattice.

There is also a useful literature distinction. Every complete arithmetic
progression is a sufficient set in the common-future sense: proving home
for every member of one such progression would settle Collatz. This
does not forbid descent rules on arithmetic progressions; it shows why
their dependency closure matters. [Monks et al., introduction and
Theorem 4.1](https://arxiv.org/pdf/1204.3904).

The live board is original source / composite dependency / uncovered
domain / unbounded depth / integer guard / modular phase. Anchor:
partitioned global completion. Niche: finite obligations for an expanding
word pair. Wildcard: a three-phase refinement that preserves modulo 19.
The research moves used are to separate local support from coverage,
compose operations before testing their final constraint, and preserve
the coordinate lost by a quotient.

The next proof target is **entry into the union of paid constructions
with no infinite residual continuation**. Finite applicability tests,
smaller dependencies, first-hit transport and fair search are available.
Universal coverage or a well-founded rule for the remaining domain is
the missing statement. A larger census, denser residue coverage or the
modulo-19 phase law alone does not establish it.
