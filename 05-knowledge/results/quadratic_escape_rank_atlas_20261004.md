# Quadratic escape components and the rank coordinate that Collatz cannot discard

2026-10-04. **INHERITED PROVED:** the three small quadratic graphs,
Chebyshev trace identity, rational denominator growth, and the finite
rational-chart conjugacy obstruction. **PROVED:** the complete escaping
integer-component codec and the finite monotone-rank-atlas obstruction,
including bounded forward macros; the word-doubling application of the
inherited Chebyshev identity. **FINITE-EXACT:** the declared controls.
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

## 6. The first and third quadratics act on repeated word presentations

The later incoming [branch-gap audit, section7](collatz_branch_toll_rank_20261004.md#the-branch-gap-and-the-negative-cycle-denominator-are-the-same-invariant)
identifies an exact cancellation: doubling an ordered Collatz word multiplies
its carry and its raw fixed-point denominator by the same factor. Combining
this with the inherited Chebyshev identity gives a precise connection to two
of the user's maps. The object being iterated is a **word presentation**,
not a Collatz starting integer.

For a nonempty positive valuation word \(w=(a_1,\ldots,a_r)\), put
\[
 A=\sum_i a_i,\qquad u=3^r,\qquad v=2^A,\qquad
 M_w=\begin{pmatrix}u&B_w\\0&v\end{pmatrix},\qquad
 F_w(n)=\frac{un+B_w}{v}.
\]
These are formal affine data; using \(w\) on an integer requires its actual
valuation guard. Concatenating \(w\) with itself gives
\[
 M_{ww}=M_w^2
 =\begin{pmatrix}u^2&B_w(u+v)\\0&v^2\end{pmatrix}.                 \tag{6}
\]
Let \(\lambda=u/v\), and use the exact rational trace coordinate
\[
 J(\lambda)=\lambda+\lambda^{-1}
 =\frac{u^2+v^2}{uv}>2.
\]
The strict inequality follows from \(u\ne v\), by unique factorization.
Then word doubling simultaneously realizes
\[
 \boxed{\lambda\longmapsto\lambda^2,\qquad
 J\longmapsto J^2-2.}                                           \tag{7}
\]
Indeed \(J(\lambda^2)=(\lambda+\lambda^{-1})^2-2\).
Equivalently the determinant-one normalization has trace
\(X=(u+v)/\sqrt{uv}\), with \(X^2=J+2\), and
\(X(M_w^2)=X(M_w)^2-2\). No irrational arithmetic is needed for the
rational version (7). The middle map \(x^2-1\) remains a separate member
of the quadratic classification; this construction gives no corresponding
Collatz-word operation for it.

For the one-letter word \((1)\),
\(\lambda=3/2\) and \(J=13/6\), so the first doubling gives \(97/36\).
After \(j\) doublings the reduced denominator of \(J\) is exactly
\[
 (uv)^{2^j}=(3^r2^A)^{2^j},                                    \tag{8}
\]
because its numerator \(u^{2^{j+1}}+v^{2^{j+1}}\) is coprime to both
2 and 3. Thus these rational coordinates escape in real size as well as
arithmetic denominator. This describes increasing presentation length.
It does not imply escape of any positive Collatz source: even when
\(\lambda<1\), its trace \(J>2\) grows under (7).

The affine anchor is unchanged:
\[
 \rho_w=-\frac{B_w}{u-v}
 =-\frac{B_w(u+v)}{u^2-v^2}=\rho_{ww}.                           \tag{9}
\]
Every prime factor supplied by the common multiplier \(u+v\) cancels in
this reduced quotient. The incoming example \(w=(1,2)\) has
\((u,B_w,v)=(9,5,8)\); doubling gives \((81,85,64)\). The raw gap
acquires the factor17 but the anchor stays \(-5\). This is a cancellation
in the repeated presentation, not evidence of a new signed basin or
new primitive prime phase.

There are two separate information issues. On the ambient positive
rationals, \(J\) identifies \(\lambda\) with \(\lambda^{-1}\).
The positive-word type repairs this sheet ambiguity: the reciprocal
of \(3^r/2^A\) cannot be another permitted slope. Indeed the reduced
denominator \(3^r2^A\) of \(J\) recovers both counts. A large trace
by itself still does not imply an expanding slope.
Even retaining the typed slope loses the carry:
the words \((1,2)\) and \((2,1)\) have the same
\(J=145/72\), but anchors \(-5\) and \(-7\). Retaining the ordered word,
or at least its full affine data together with exact valuation guards,
repairs the relevant losses. An anchor alone also discards length and
the binary precision consumed by an application.

In particular, replacing a checked word by its formal square does not
prove that square is legal at the same source. The actual valuations
starting from27 are \(1,2,1,1\); the word \((1,2)\) is legal, but
\((1,2,1,2)\) is not. The doubled word uses \(2A\) halvings, and needs
its own source guard. The [guarded-pumping result](collatz_guarded_pumping_memory_20261004.md)
separately measures the finite repeated-word fuel. Neither a fixed anchor
nor the trace identity supplies that fuel or a root certificate.

**Connection contract.** Source: ordered affine word presentations under
concatenation. Target: rational squaring and Chebyshev dynamics.
Map: \(w\mapsto\lambda_w\mapsto J(\lambda_w)\). Preserved predicate:
the doubling law and, with the carry retained, the affine fixed point.
Lost data: carry, order, source legality and first-hit status; the ambient
sheet is also lost if the positive-word arithmetic type is forgotten.
Needed sidecar: that type, ordered word and its guarded source or supplied proof.
The smallest tests are one-letter squaring, the two carries \(12/21\),
and the failed doubled word at27. This connects the user's quadratic
maps exactly while staying outside an integer-orbit conjugacy claim.

## 7. Exact controls

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
The word-matrix controls cover every word of length1--4 with letters1--4
at doubling depths0--3: 1360 exact checks of coefficients, both trace
identities, denominator growth and anchor cancellation. They retain the
same-trace/different-carry and illegal-repeat hostiles.

The finite rank check enumerates every assignment of \(q\) polynomial
charts to \(q+1\) increasing states for \(q=1,\ldots,4\): 1114 cases.
Finite unrolled charts before repetition are positive controls.
The bounded-macro check enumerates all macro lengths in \(\{1,2,3\}\)
and all chart assignments for \(q=1,2,3\): 2262 cases.
The valuation rank is checked on 6240 actual growing edges, with
reset and signed-fold hostiles retained. These finite tests support,
but do not replace, the general proofs. All checks survive optimization.
