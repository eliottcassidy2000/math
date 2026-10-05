# Quadratic escape components and the rank coordinate that Collatz cannot discard

2026-10-04. **INHERITED PROVED:** the three small quadratic graphs,
Chebyshev trace identity, rational denominator growth, and the finite
rational-chart conjugacy obstruction. **PROVED:** the complete escaping
integer-component codec and the finite monotone-rank-atlas obstruction,
including bounded forward macros. **FINITE-EXACT:** the declared controls.
Universal Collatz completion remains **OPEN**. No priority claim.

## 1. Recover the objects before transplanting their pattern

The closest sources are [arithmetic-braids geometry](arithmetic_braids_20260917_geometry.md),
the [row/braid sign and PCF audit, section5](collatz_mod6_20260917_row_braid_typing.md),
and [crossroads arithmetic, sections3--5](crossroads_crossing_20260926_arithmetic.md).
They already identify the user's three maps as
\(x^2,\ x^2-1,\ x^2-2\), prove their finite rational preperiodic sets,
and keep the sign, denominator and map degree separate.

The relevant incoming mechanism is the [branch-toll graph rank](collatz_branch_toll_rank_20261004.md).
Its proper input-defined rank allows some larger integers as lower-ranked
proof obligations, but retains an explicit infinite critical set. Its
anchor controller uses unbounded valuation precision as fuel. Neither
result claims that every critical state has been grounded.

Canonical hostile: \(6/5\) has bounded real motion under \(x^2-2\) and
unbounded rational denominator. Corrected near miss: finite critical orbit
does not mean every integer is attracted to that orbit. Least-used
sidecars: the root of an escaping integer component, and the number of
distinct monotone charts needed along a growing run.

The live board is **integer component / analytic critical point / arithmetic
height / graph rank / chart label / valuation fuel**. Anchor: identify a
usable boundary for Collatz ranking. Niche: classify all integer components
of the three quadratics. Wildcard: make the loss in a signed fold explicit.

Put \(f_c(x)=x^2-c\), \(c\in\{0,1,2\}\). The recovered finite parts are:

| Map | Integer preperiodic set | Complete edges in that set |
|---|---|---|
| \(x^2\) | \(\{-1,0,1\}\) | \(0\to0,\ -1\to1\to1\) |
| \(x^2-1\) | \(\{-1,0,1\}\) | \(1\to0\to-1\to0\) |
| \(x^2-2\) | \(\{-2,-1,0,1,2\}\) | \(1\to-1\to-1,\ 0\to-2\to2\to2\) |

The derivative-critical point is zero, and its forward orbit is finite
in each case. That is a statement about one marked starting point.
It does not assert global attraction, connectedness of the integer
functional graph, or vanishing of a proof rank's local minima.

## 2. A complete codec for every escaping integer component

Set \(B_0=B_1=2\), \(B_2=3\). For every integer \(|x|\ge B_c\),
\[
 f_c(x)>|x|,\qquad f_c(x)\ge B_c.
\]
Thus all integers outside the table escape to positive infinity; after
the first iterate the values increase strictly.

**PROVED component classification.** Each escaping component of the
undirected functional graph is uniquely indexed by an integer
\[
 b\ge B_c,\qquad b+c\ \hbox{is not a square}.                         \tag{1}
\]
Its vertices are exactly
\[
 \{\ \pm f_c^j(b):j\ge0\ \}.                                        \tag{2}
\]
Every escaping integer has a unique address \((b,j,\varepsilon)\),
\(\varepsilon\in\{-1,1\}\), meaning \(x=\varepsilon f_c^j(b)\).
The forward map on addresses is
\[
 (b,j,\varepsilon)\longmapsto(b,j+1,+1).                            \tag{3}
\]

**Constructive proof.** Begin with \(|x|\). Whenever its value \(y\)
has \(y+c\) an integer square, replace it by the positive integer
\(\sqrt{y+c}\). This is an actual inverse edge and is smaller than \(y\)
throughout the escaping region. It stays in that region, since a
preperiodic inverse would make \(y\) preperiodic. Magnitude therefore
decreases to a unique \(b\) satisfying (1). Count the inverse steps
and retain the original sign to recover \(j,\varepsilon\).
The only integer preimages of a positive value are the two signed
square roots. This proves exhaustiveness and uniqueness of (2).
Distinct labels \(b\) cannot have a common future.

For every integer \(X\ge B_c\), the number of component labels through
\(X\) is exactly
\[
 \#\{b\le X:\text{(1)}\}=X-\lfloor\sqrt{X+c}\rfloor.                 \tag{4}
\]
Subtract the squares in \([B_c+c,X+c]\) from the \(X-B_c+1\) candidates.
In particular, there are infinitely many escaping components, despite
the finite derivative-critical orbit.

The least examples are
\[
\begin{array}{ll}
c=0:&2\to4\to16\to256,\\
c=1:&2\to3\to8\to63,\\
c=2:&3\to7\to47\to2207.
\end{array}
\]
The familiar squaring forest is inherited; the same inverse-strip
description covers the two translated maps as well.

The codec preserves integer identity and common-future component;
discarding \(j\) loses height, and discarding the sign loses the original
input. The equality \(f_c(x)=f_c(-x)\) makes the sign fold a genuine
common-future move. After that fold, inverse stripping is a proper
decreasing reduction, but its terminal label usually represents an
escaping component rather than a periodic basin. A terminating graph
reduction is therefore not itself a theorem of attraction to a chosen root.

## 3. Chebyshev preserves a trace, not arithmetic boundedness

For \(J(z)=z+z^{-1}\), direct multiplication gives
\[
 J(z^2)=J(z)^2-2.
\]
This is the inherited semiconjugacy for the third map. It identifies
\(z\) and \(z^{-1}\); it is not a faithful phase coordinate without a
sheet choice.

Use the signed Gaussian point \(z=(3+4i)/5\). Under repeated squaring,
the integer lift
\[
 (a,b,C)\longmapsto(a^2-b^2,2ab,C^2)
\]
retains \(a^2+b^2=C^2\). Its trace \(x=2a/C\) stays in \([-2,2]\),
but the reduced denominator after \(j\) steps is \(5^{2^j}\).
The orbit beginning at \(6/5\) is therefore not preperiodic. The
underlying general argument is already in the recovered notes: a
nonintegral reduced \(a/b\) under any integer monic quadratic has
successive reduced denominators \(b^{2^j}\).

The [existing rational-chart obstruction](crossroads_crossing_20260926_arithmetic.md#5-a-finite-rational-chart-obstruction-to-a-collatz-quadratic-model)
also already rules out an infinite-range encoding by finitely many rational
charts of shortcut Collatz into any fixed rational map of degree at least two:
along a recurrent chart transition the degrees would multiply by that
degree. More precisely, any such finite-chart semiconjugacy has finite
range. We do not claim this older result as a new obstruction, nor
exclude nonrational coordinates, infinitely many charts or unbounded
return times.

## 4. A new obstruction for rank atlases, with no conjugacy assumption

The following argument uses only finite actual increasing runs. It does
not assume Collatz convergence or divergence.

**PROVED finite monotone-atlas obstruction.** Fix a finite number \(q\)
of rank charts \(H_1,\ldots,H_q\) taking values in a common ordered set.
Suppose each is eventually nondecreasing in its positive integer input,
on the domain where it is used. Chart selection may depend arbitrarily
on the integer or on finite controller history. These charts cannot
supply a rank that strictly decreases on every positive odd forward
Collatz edge.

More generally, fix any \(L\ge1\). They cannot supply a total dispatcher
that always chooses an actual forward macro of between one and \(L\)
odd steps and makes the rank strictly decrease across that macro.
This includes a finite library of fixed forward words.

**Proof.** For any \(N\ge1\) and positive odd \(u\), the actual values
\[
 n_j=3^j2^{N+1-j}u-1,\qquad 0\le j\le N,                            \tag{5}
\]
give \(N\) consecutive valuation-one edges, and increase strictly.
Choose \(N=qL\), and choose \(u\) large enough that every value is
above every chart's monotonicity threshold and any finite exceptional
domain. Follow \(q\) purported dispatcher moves. Each uses at most
\(L\) edges, so they remain within (5). Their \(q+1\) endpoints
increase strictly and must repeat a chart. At that repeated chart
the later rank is at least the earlier one, contrary to strict
decrease across the intervening macros. For single edges take \(L=1\).

The proof does not need regular, periodic or computable chart selection.
The number of charts needed to rank a particular increasing block is
at least the number of its ranked endpoints.

**Corollaries.** No finite atlas of proper nonnegative polynomial or
rational height functions can rank every odd forward edge or every
chosen bounded-length forward macro. A proper rational function tends
to positive infinity and is eventually increasing. Finitely many
quasipolynomial residue pieces are covered too. The same obstruction
applies to finite polynomial-tuple atlases in lexicographic order:
the first nonconstant coordinate is eventually increasing if the
coordinates are eventually nonnegative; an entirely constant chart
is nondecreasing as well.

This is different from the older degree obstruction: it does not
assume a quadratic model or any semiconjugacy. It strengthens the
incoming same-residue affine-rank hostile to arbitrary finite
monotone chart selection and bounded forward macros.

## 5. What escapes the obstruction, and what it still owes

The incoming rank's unbounded valuation coordinate is essential to
the scope. For a simple example write
\[
 n+1=2^Kt,\qquad t>0\text{ odd},\qquad
 R_-(n)=(3^Kt,K).
\]
This is an input-defined proper lexicographic rank, since
\(3^Kt\ge n+1\). On every valuation-one edge,
\[
 (t,K)\longmapsto(3t,K-1),
\]
so the first coordinate stays fixed and the second decreases. This
one-mode valuation chart pays for arbitrarily long growing one-runs.
It is not eventually nondecreasing in \(n\), and cannot be replaced
by finitely many eventually nondecreasing charts while retaining
those rank comparisons. This is the negative-anchor special case
of the incoming anchor mechanism, not a newly claimed universal rank.

The reset debt remains:
\[
 9\to7:\quad 3^Kt\text{ rises from }15\text{ to }27,
\qquad
 27\to41\to31:\quad 63,63,243.
\]
Even a numerical decrease can increase this particular energy.
A correct controller still needs covered reset guards, paid exits,
and grounded terminal dependencies.

Nor does the theorem prevent inverse or common-future moves, variable
unbounded return lengths, or a nonmonotone valuation rank. The incoming
root-centered rank \(R_+(n)\) operates on the undirected proof graph and
is explicitly outside the all-forward hypothesis.

One sign distinction is decisive. For a quadratic, \(x\) and \(-x\)
have the same next value. For the incoming Collatz rank,
\(R_+(n)=R_+(2-n)\) is only rank equality. Already
\[
 R_+(3)=R_+(-1),\qquad 3\to5\to1,\quad -1\to-1
\]
has two distinct signed basins. Folding this equality into an actual
proof edge would lose the basin being certified. The incoming work
correctly retains the signed state and oriented carry.

Thus the transferable mechanism is a proper graph reduction with a
terminal-set obligation and an explicit chart sidecar. The quadratic
component labels show why terminal classification matters; the finite
atlas theorem shows why the valuation fuel cannot be replaced by a
finite collection of ordinary polynomial height pieces. Neither
statement grounds every positive Collatz critical state.

## 6. Exact controls

Program:
[quadratic_escape_rank_atlas_20261004.py](../../04-computation/experiments/quadratic_escape_rank_atlas_20261004.py).
Output:
[quadratic_escape_rank_atlas_20261004.out](quadratic_escape_rank_atlas_20261004.out).

    python 04-computation/experiments/quadratic_escape_rank_atlas_20261004.py
    python -O 04-computation/experiments/quadratic_escape_rank_atlas_20261004.py

Controls cover all escaping integers in \([-5000,5000]\) for the three
maps: 29992 lossless codec and transition checks, with 4930 component
labels through 5000 for each map. The Gaussian lift checks nine exact
stages, retaining bounded real trace and the explicit denominator.

The finite rank check enumerates every assignment of \(q\) polynomial
charts to \(q+1\) increasing states for \(q=1,\ldots,4\): 1114 cases.
Finite unrolled charts before repetition are positive controls.
The bounded-macro check enumerates all macro lengths in \(\{1,2,3\}\)
and all chart assignments for \(q=1,2,3\): 2262 cases.
The valuation rank is checked on 6240 actual growing edges, with
reset and signed-fold hostiles retained. These finite tests support,
but do not replace, the general proofs. All checks survive optimization.
