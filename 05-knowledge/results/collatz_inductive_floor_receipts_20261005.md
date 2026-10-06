# Positive measurement bounds from a well-founded descent certificate

2026-10-05. **PROVED:** guarded first-descent induction, a canonical receipt
compiler, proper infinite binary and mixed four-child families, and a positive
localized measurement bound on its designated leaves. **FINITE-EXACT:** the
stated controls. **OPEN:** covering every positive odd source by grounded
descent obligations. Production algorithms take neither a target ROOT word
nor measured moment values as inputs.

[Script](../../04-computation/experiments/collatz_inductive_floor_receipts_20261005.py)
and [exact output](collatz_inductive_floor_receipts_20261005.out).

## 1. Inheritance and the missing premise

The closest mechanism is the [floor transport
theorem](collatz_floor_transport_deadlines_20261005.md): a truthful positive
floor at an actual endpoint propagates backward through a checked finite
word. The [source-refinement note](collatz_source_refinement_floor_20261005.md)
adds signed common-future bills, cancellation before scalarization, and a
proper binary tree whose individual forward edges decrease.

The new [shadow note](collatz_refinement_floor_shadows_20261005.md) emphasizes
smaller **ancestors**. That differs from our forward path to a smaller
**descendant**. For example, 7 has the exact first descent

\[
 7\longrightarrow11\longrightarrow17\longrightarrow13\longrightarrow5,
 \qquad w=(1,1,2,3).
\]

A floor at 5 therefore gives a floor at 7. A classification of smaller
ancestors cannot exclude this move. We inherit the guarded affine and
counter identities, not empirical densities as all-height statements.

The canonical hostile is an ungrounded circular implication:
3 and 13 have the same future 5, but replacing both unknown positivity claims
by one another supplies no positive premise. The corrected near miss is to
treat a numerical assumed floor as already established. The least-used
sidecar is a strict ordinary-source rank at every proof obligation.

The board is **actual source / first-descent cut / smaller obligation /
retained counters / analytic family floor / localized readout**.

## 2. Canonical induction and explicit obligations

Let $U(n)=\operatorname{oddpart}(3n+1)$, stopping at ROOT 1. A
first-descent receipt for odd $n>1$ consists of its actual nonempty valuation
word $w=(a_1,\ldots,a_r)$ and endpoint $h$, with

\[
 U^i(n)\ge n\quad(1\le i<r),\qquad U^r(n)=h<n.                 \tag{1}
\]

All valuations are checked exactly. Determinism gives at most one such
receipt per source. It exists exactly when that source eventually visits a
smaller odd integer. This is uniqueness of an actual first-descent word,
not uniqueness of arbitrary common-future proofs or source programs.

For a rooted nonroot source write

\[
 L=\tau-1,\quad K=\sum_i\lfloor(a_i-1)/2\rfloor,\qquad
 W(n)=w(L,K)=\frac{2K!(L+1)!}{(L+K+2)!}.
\]

Set $W(1)=1$ and $W(n)=0$ when $n$ is not rooted. For a receipt ending at
$h>1$, put $k=\sum_{i=1}^r\lfloor(a_i-1)/2\rfloor$. If an independently
justified floor is $W(h)\ge\epsilon>0$, let

\[
 B=\max\{b\ge0:\epsilon(b+1)(b+2)\le2\}.
\]

The inherited counter ceiling and factorial ratio give

\[
 W(n)\ge
 \epsilon\,\frac{(k+1)!(r+1)!}{(B+3)_{r+k}}>0.                \tag{2}
\]

Here $(x)_j$ is the rising factorial. Indeed $K(h)\ge1$,
$L(h)+K(h)\le B$, and the exact ratio is
$(K(h)+1)_k(L(h)+2)_r/(L(h)+K(h)+3)_{r+k}$.
If $h=1$, the receipt itself gives its exact weight; (2) is not used at ROOT.

**Well-founded induction theorem.** If every odd $n>1$ is assigned a checked
receipt (1), then every positive odd integer is rooted and obtains a positive
floor recursively from $W(1)=1$. A least unrooted source would have a receipt
to a smaller rooted source, a contradiction. Conversely universal rootedness
supplies a first descent for every $n>1$. Universal existence of these receipts
is therefore an equivalent remaining coverage obligation, not a proved premise.

The executable compiler consumes a finite dictionary of supplied receipts.
Every dependency strictly decreases, so checking cannot loop. It returns:

* **grounded:** the chain ends at 1; its accumulated word gives an exact weight;
* **conditional:** it ends at an explicitly supplied numerical floor assumption,
  which remains an exact point obligation;
* **pending:** it ends at a source with neither a receipt nor an assumed floor.

No missing suffix is explored by the compiler. It adds lengths and sibling
costs before applying (2), preserving the original anchor. A ROOT certificate
can subsequently be exported by concatenating the checked segments. It is an
output, not an input recycled to justify its own floor.

A finite rule collection is not a proof that its guards cover all integers.
Finitely many symbolic schemas can also have infinitely many instances;
each instance still needs its guard and decreasing dependency.
Actual common-future relations $3\leftrightarrow13$ cannot be installed as a
cycle here: the increasing dependency violates the rank condition.
Homogeneous conditional inequalities alone remain compatible with both
unknown weights being zero.

## 3. An infinite domain that rises before it descends

This construction extends the kind of induction domain, not the known ROOT
basin. It uses the standard inverse-word mechanism with a finite palette.

At a unit parent $y>1$, choose the unique $a_0\in\{3,4,5,6,7,8\}$ satisfying
$2^{a_0+1}y\equiv5\pmod9$. Uniqueness follows from the order six of 2 modulo 9.
For $a\in\{a_0,a_0+6,a_0+12\}$ set

\[
 n_a=\frac{2^{a+1}y-5}{9}.                                  \tag{3}
\]

All three are positive odd integers. Exactly one is divisible by three;
retain the other two. Indeed

\[
 n_{a+6}-n_a=7\,2^{a+1}y\equiv2\pmod3,                       \tag{4}
\]

using $2^{a+1}y\equiv5\pmod9$. The three residues modulo 3 are distinct.

Every retained child has $n_a\equiv3\pmod8$, and

\[
 U(n_a)=\frac{2^a y-1}{3}>n_a,\qquad U^2(n_a)=y<n_a.         \tag{5}
\]

The first valuation is exactly one, the intermediate value is odd, and its
next numerator is exactly $2^a y$. Thus $(1,a)$ is the canonical first-descent
word, with no early ROOT occurrence.

Start at 5 and repeatedly retain these two children. Unique forward parents
and strict size increase give a genuine infinite binary tree $\mathcal T_2$.
Its first children are

\[
 35\xrightarrow{(1,5)}5,\qquad 2275\xrightarrow{(1,11)}5.
\]

Every nonroot vertex starts with valuation one. Hence, apart from 5, this tree
is disjoint from the source-refinement note's SF5 tree, whose child-to-parent
edges have valuation at least two. It is proper: the already rooted sources
7 and 27 are rejected. No new Collatz basin is claimed.

Membership is decidable from the source alone. Check its two actual steps,
verify that it belongs to the parent's palette, and replace it by that smaller
parent. Stop successfully at 5; reject at the first failure. Every successful
replacement strictly decreases a positive integer, so recognition terminates
on nonmembers as well as members.

## 4. An analytic floor forces a positive localized measurement

At tree depth $t$, the rank and sibling cost satisfy

\[
 \tau=2t+1,\qquad L=2t,\qquad K\le9t+1,                     \tag{6}
\]

because $a\le20$ and source 5 has counters $(0,1)$. Furthermore

\[
 n_a-1\ge\frac{16}{9}(y-1),\qquad n-1\ge4(16/9)^t.          \tag{7}
\]

Since $(16/9)^2>2$, $b=\operatorname{bitlength}(n)$ gives $t\le2b$.
Monotonicity in both counters proves the source-size floor

\[
 W(n)\ge w(4b,18b+1)>0,\qquad n\in\mathcal T_2.              \tag{8}
\]

Choose the least positive inverse exponent $e\in\{1,\ldots,6\}$ with
$2^e n\equiv1\pmod9$, and put

\[
 z=\rho(n)=\frac{2^e n-1}{3},\qquad m=(z-3)/6.
\]

Then $z$ is a positive odd multiple of three, $U(z)=n$ with exact exponent
$e$, and the added sibling cost is at most two. The inherited injection
identity gives $\lambda(z)=W(z)$. Therefore

\[
 p_m:=\lambda(6m+3)\ge\epsilon_n:=w(4b+1,18b+3)>0.          \tag{9}
\]

This all-height theorem uses induction and counter bounds. Neither individual
target ROOT records nor numerical measurements are premises.

Apply the inherited [localized
selector](collatz_localized_resolvent_floor_20261005.md):

\[
 h_m(j)=\frac{4\,2^{m-j}}{(1+2^{m-j})^2},\quad
 H_{m,d}=\sum_jp_jh_m(j)^d,\quad A_{m,d}=9H_{m,d+1}-8H_{m,d}.
\]

Its universal error theorem is

\[
 p_m-8(16/25)^d\le A_{m,d}\le p_m.
\]

Choose the least integer $d\ge0$ satisfying

\[
 8(16/25)^d\le\epsilon_n/2.                                 \tag{10}
\]

It follows without evaluating either moment that

\[
 A_{m,d}\ge\epsilon_n/2>0.                                   \tag{11}
\]

For $\epsilon_n=u/v$, the exact integer loop starts with left $=16v$,
right $=u$, and repeatedly multiplies them by 16 and 25 until left is at
most right. It computes the least sufficient degree without logarithmic
rounding or evaluating any $H$.

The direction is deliberate: prove a floor on a constructive family, then
force a localized readout to be positive. The readout is not an unproved
premise. Evaluating actual moments is not free, and a source outside this
family obtains no floor from this theorem. Membership validation itself
constitutes a finite structural proof; we do not claim positivity is easier
to prove than finding some ROOT certificate.

Truthful information refinement at this same source preserves (9) and (11).
Changing the source, or using an arbitrary probability in place of the actual
$\lambda$, is not licensed by this result.

### Combining the independent palettes

At every parent also allow the two old SF5 unit children
$(2^a y-1)/3$: use $a\in\{2,4,6\}$ when $y\equiv1\bmod3$ and
$a\in\{3,5,7\}$ when $y\equiv2\bmod3$, omitting the one child divisible
by three. Increasing $a$ by two changes the child modulo three by one,
so exactly two survive. Their actual forward valuations are at least two,
and their parent is strictly smaller.

Together with (3), these give **four** children at every parent. A source's
first valuation distinguishes the palettes: valuation one uses the new
two-step rule, while a larger valuation uses the old one-step rule.
Each accepted block ends at its first strict descent. Thus the mixed tree
$\mathcal T_{1,2}$ has a canonical parent, four distinct children, disjoint
generations, and a terminating source recognizer.

It strictly exceeds the union of the two pure trees. The exact witness is

\[
 739\xrightarrow{(1,8)}13\xrightarrow{(3)}5\xrightarrow{(4)}1.
\]

Its first block excludes the pure one-step tree, while its intermediate
parent 13 excludes the pure two-step tree. No new seed was added: the root
remains 5 with the single checked edge to 1.

Every child satisfies $n-1\ge(4/3)(y-1)$. Since $(4/3)^3>2$, macro depth
$t$ is at most $3b$, where $b=\operatorname{bitlength}(n)$. Each macro
contributes at most two odd steps and nine sibling-cost units. Consequently

\[
 L(n)\le6b,\quad K(n)\le27b+1,\quad
 \lambda(\rho(n))\ge\epsilon_n^{\rm mix}:=w(6b+1,27b+3)>0.
 \tag{12}
\]

The same degree loop gives $A_{m,d}\ge\epsilon_n^{\rm mix}/2$.
For 739 its designated leaf is 15765 and the displayed bound gives $d=371$.
These bounds allow any interleaving of the two guarded constructions,
rather than only an initial choice of one pure tree. They remain a proper
domain: for instance 7 fails both permitted first-descent palettes.

## 5. Reproduction, controls, and limits

Run from the repository root:

    python 04-computation/experiments/collatz_inductive_floor_receipts_20261005.py
    python -O 04-computation/experiments/collatz_inductive_floor_receipts_20261005.py

The finite induction universe is the 256 odd sources from 1 through 511.
Every source is observed only to its first strict descent or to the declared
per-source cap. No completed nonroot routes are supplied as seeds.

| Local cap | Checked descent receipts | ROOT-grounded requests | Literal odd-step queries |
|---:|---:|---:|---:|
| 1 | 127 | 21 | 255 |
| 2 | 159 | 29 | 383 |
| 4 | 203 | 126 | 543 |
| 8 | 229 | 150 | 696 |
| 16 | 241 | 165 | 848 |

These are separate runs, not incremental-query totals. Grounded words are
compared with fresh literal replay only after the proof output has been
produced. This audit does not feed the compiler or close pending obligations.

The binary-tree control contains every vertex through depth eight: 511
distinct sources, up to 114 bits. Exact guards, parent identities, first
descents, depths and analytic floors are checked. The first four levels give
15 independent leaf and measurement-degree controls. Source 35 selects leaf
93; (9) and the exact loop give $d=160$. These conservative constants come
from the same all-height formula, not an individual stored ROOT word.

The mixed-tree control contains every vertex through macro depth five:
1,365 sources, with level sizes $1,4,16,64,256,1024$. It independently
checks source parsing, distinct generations, counters and the floor, and
replays the mixed witness's designated leaf after the analytic output.

Sixteen hostile controls reject malformed types, mismatched sources,
incorrect endpoints or valuations, an increasing dependency, a cut later
than the first descent, invalid floors and invalid palettes. An assumed
floor remains conditional until grounded receipts discharge it.

There are 14,225 exact checks. The result gives reusable positive localized
measurement guarantees on a proper infinite induction domain. Global entry
or descent coverage remains open. A finite benchmark, uniquely selected
certificate, or conditional positive number does not discharge that task.
