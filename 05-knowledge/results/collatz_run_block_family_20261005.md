# Unbounded run blocks with a source-bounded induction expense

2026-10-05. **PROVED:** exact inverse phases, a canonical unbounded-run
grammar, a terminating membership test, source-only counter and measurement
bounds, and strict enlargement of the previous mixed tree.
**FINITE-EXACT:** the specified finite controls. **OPEN:** entry for every
source. Neither a finite ROOT bank nor measured moment values are inputs
to the production floor routine.

[Script](../../04-computation/experiments/collatz_run_block_family_20261005.py)
and [exact output](collatz_run_block_family_20261005.out).

## 1. Inheritance and the stronger question

The [inductive-floor companion](collatz_inductive_floor_receipts_20261005.md)
combines two guarded inverse palettes: one actual odd step, or a valuation-one
step followed by a terminal reset. Their mixture has four children per
parent and a canonical first-descent parser. The older
[complement-routing construction](collatz_complement_routing_20261004.md)
and [early-reroute compiler](collatz_early_reroute_20261004.md) already keep
the actual initial-one run and ternary integrality when constructing
common-future dependencies. Here the new parent is an **actual descendant**;
no common-future suffix substitution is needed.

The anchor is to allow unbounded initial-one runs without losing source
rank or a quantitative floor. The niche is an arithmetic expense attached
to each terminal reset. The wildcard is whether the same expense survives
a weaker, exact-size descent rule. The board is **maximal one-run /
ternary exponent phase / first descent / multiplicative height /
source expense / positive localized measurement**.

The old binary and four-child trees are sufficient domains, not universal
coverage. The canonical hostile remains 7: its first two valuation-one
steps end at 17, whose valuation-two reset reaches 13, still above 7.
The corrected near miss below is to keep a useful size-decreasing block
while assuming it obeys the stronger expense inequality.

## 2. An exact all-run inverse block

Fix an odd unit parent $y>1$, a run $r\ge0$, and a terminal valuation $a\ge2$.
Solving the word $(1^r,a)$ backward gives

\[
 n+1=\frac{2^{r+a}y+2^{r+1}}{3^{r+1}}
     =\frac{2^{r+1}(2^{a-1}y+1)}{3^{r+1}}.             \tag{1}
\]

Thus integrality is equivalent to

\[
 2^{a-1}y\equiv-1\pmod{3^{r+1}}.                      \tag{2}
\]

When (2) holds and $n$ is positive, the quotient
$(2^{a-1}y+1)/3^{r+1}$ is odd. Consequently

\[
 v_2(n+1)=r+1.                                        \tag{3}
\]

The first $r$ valuations are exactly one. After these steps the checkpoint is
$x=(3^r/2^r)(n+1)-1$, and $3x+1=2^a y$ gives the exact last valuation.
Every initial-one step from a positive source strictly grows.

We impose the simple stronger guard

\[
 a\ge r+2,\qquad 3\nmid n.                             \tag{4}
\]

Let $\rho=2^{r+a}/3^{r+1}$ and $\beta=(2/3)^{r+1}$.
Then $n=\rho y+\beta-1$. The key estimates are

\[
 \begin{split}
 \rho+\beta
 &\ge(4/3)^{r+1}+(2/3)^{r+1}\ge2,\\
 n-1&\ge\rho(y-1),\\
 \rho
 &=2^{a-1}(2/3)^{r+1}\ge(4/3)^{a-1}>1.
 \end{split}                                         \tag{5}
\]

For the first line, the sum is two at exponent one and is increasing:
its consecutive difference is
$[(4/3)^j-(2/3)^j]/3>0$.
The final line uses $r+1\le a-1$. In particular $n>y$.
Thus the terminal reset is the **first strict descent** below $n$.
All earlier states exceed $n$, and no ROOT padding occurs.

The guard is sufficient; it is not claimed necessary for every paid block.

## 3. Exponent phases and an infinite native palette

The order of 2 modulo $3^{r+1}$ is $2\cdot3^r$.
An elementary proof starts with $v_3(4-1)=1$ and notes that cubing a number
$1+3^j u$, $3\nmid u$, increases the valuation of its difference from one
by exactly one. The order is therefore the full number of units.

For each $(y,r)$, (2) selects one exponent class modulo

\[
 P_r=2\cdot3^r.
\]

Let $a_0$ be its least representative at least $r+2$.
Every $a=a_0+P_r t$, $t\ge0$, is integral and satisfies the strong guard.
Among each three consecutive values of $t$, exactly one child is divisible
by three and two are units. Indeed

\[
 n_{a+P_r}-n_a
 =2^{r+a}y\,\frac{2^{P_r}-1}{3^{r+1}},
\]

and the quotient is a unit modulo three by the same valuation calculation.
Successive differences are the same nonzero residue modulo three, since
$2^{P_r}\equiv1\bmod3$.

The program obtains the phase by one ternary-digit lift at each level.
It need not search an interval containing $3^r$ exponents.
Literal construction of the resulting integer can still require an enormous
number of bits; phase compression does not make those source bits disappear.

For $r=0$, the first three phase exponents are exactly the old SF5 palette
$\{2,4,6\}$ or $\{3,5,7\}$ before deleting the three-divisible child.
For $r=1$, they are exactly the old two-step palette.
Keeping **all** permitted exponents and **all** $r\ge0$ strictly generalizes
both palettes, including their entire mixed tree.

## 4. The canonical tree and its source-only expense

Start at 5. At each unit parent retain every child satisfying (1)--(4).
This is an infinite tree, with countably many children per parent.
Different runs are distinguished by (3); within a run, the actual terminal
valuation and actual endpoint determine the parent uniquely.
Strict source growth prevents generation cycles or merging across depths.

For an input source, recover $r=v_2(n+1)-1$, compute its checkpoint and
terminal valuation, and check (4), the unit parent, and the exact source
identity. Replace an accepted source by its smaller parent. Stop at 5;
reject at a failed guard. This parser always terminates and does not consult
a stored ROOT suffix. The final source 5 uses its one checked word $(4)$.

The selected root boundary is deliberate. Other direct ROOT rays, such as
21 or the unit source 85, are not thereby members. ROOT 1 is not given a
self-loop. Adding a different terminal rule would be a separate extension.

For a member built from blocks $(r_i,a_i)$, define the total expense

\[
 E=\sum_i(a_i-1).
\]

Multiplying (5) along its decreasing parent chain gives the all-height
source inequality

\[
 4(4/3)^E\le n-1.                                     \tag{6}
\]

Define the computable integer ceiling

\[
 D(n)=\max\{d\ge0:4\cdot4^d\le(n-1)3^d\}.              \tag{7}
\]

This is the exact integer ceiling associated with (6), not a claim that
some valid family route attains it for every source.
It is found by a terminating integer loop. In particular
$D(n)<3(\operatorname{bitlength}(n)-2)$ because $(4/3)^3>2$.

Every block has length $r_i+1\le a_i-1$. Including the terminal word $(4)$,

\[
 L=\tau-1=\sum_i(r_i+1)\le E\le D(n),\qquad
 K=1+\sum_i\lfloor(a_i-1)/2\rfloor
   \le1+\lfloor D(n)/2\rfloor.                         \tag{8}
\]

Thus the larger domain has sharper counter bounds than the earlier
independent maxima on run length and terminal valuation.
The odd ROOT deadline is $D(n)+1$.

The standard weight is
$w(L,K)=2K!(L+1)!/(L+K+2)!$. Its decrease in each counter proves

\[
 W(n)\ge
 w\bigl(D(n),\,1+\lfloor D(n)/2\rfloor\bigr)>0.          \tag{9}
\]

These are source-size formulas justified by the family theorem, not evaluations
of a supplied target ROOT word. Recognition is itself a structural proof;
there is no claim that all integers satisfy its guards.

A small new member is

\[
 23\longrightarrow35\longrightarrow53\longrightarrow5,
 \qquad (r,a)=(2,5).
\]

Its first two steps grow, so the old mixed tree rejects it. The new tree accepts
it. A complete check of the odd integers below 23 confirms that it is the least
new member relative to that old tree.

## 5. A positive localized measurement from the expense bound

Use the established least three-divisible predecessor
$z=\rho(n)=(2^e n-1)/3$, $1\le e\le6$, and put $m=(z-3)/6$.
Its extra length is one and its extra sibling cost is at most two.
The inherited identity $\lambda(z)=W(z)$ and (8) imply

\[
 p_m\ge\epsilon_n:=
 w\bigl(D(n)+1,\,3+\lfloor D(n)/2\rfloor\bigr)>0.        \tag{10}
\]

The [localized selector](collatz_localized_resolvent_floor_20261005.md) has
$A_{m,d}=9H_{m,d+1}-8H_{m,d}$ and universal error

\[
 p_m-8(16/25)^d\le A_{m,d}\le p_m.
\]

Choosing the least integer $d$ with $8(16/25)^d\le\epsilon_n/2$ proves

\[
 A_{m,d}\ge\epsilon_n/2>0.                              \tag{11}
\]

No moment values are computed or assumed by this routine. The bound is for the
actual injection measure on a specified constructive basin; arbitrary
probability data cannot replace it.

| Source | Expense ceiling | Designated leaf | Guaranteed degree |
|---:|---:|---:|---:|
| 5 | 0 | 3 | 14 |
| 23 | 5 | 15 | 26 |
| 35 | 7 | 93 | 30 |
| 739 | 18 | 15765 | 56 |

For the same leaves, the previous separate-counter formulas required degrees
160 at source 35 and 371 at source 739. This is a comparison of sufficient
selector degrees, not of measured moment values or computational access costs.
The source-dependent deadline compiler may provide another valid bound;
neither route claims a best possible floor.

## 6. A weaker paid block and the precise missing budget

The exact-size condition $n>y$ can hold when (4) fails.
For instance

\[
 55\xrightarrow{(1,1,3)}47<55.
\]

Here the block length is three but its proposed expense is only two.
Also

\[
 \frac{55-1}{47-1}=\frac{27}{23}<(4/3)^2.
\]

Thus both the length bill and the multiplicative expense bound fail if this
block is admitted without a different account. The failure is in the proposed
bound, not in the validity of the actual descent.

A grammar allowing every exact first-reset descent remains decidable by
strict size decrease, but (6)--(11) do not automatically extend to it.
The next productive target is a different source-bounded expense that pays
such weaker resets. This package does not claim that the displayed block's
parent 47 lies in its present domain.

## 7. Reproduction and exact controls

Run from the repository root:

    python 04-computation/experiments/collatz_run_block_family_20261005.py
    python -O 04-computation/experiments/collatz_run_block_family_20261005.py

The experiment checks:

* 168 native child blocks: parents 5,13,35,739; runs 0 through 6; three
  successive exponent groups; two unit children per group;
* every vertex through depth three using runs 0 through 3 and the first
  phase group: 585 sources, level sizes 1,8,64,512, up to 364 bits;
* all 341 old mixed-tree vertices through depth four, all retained;
* five exact measurement-budget controls, including one run-five source;
* the weaker paid block, excluded direct ROOT boundaries, 7 and 27, and ten
  malformed input/guard controls.

Each generated source's word is compared with fresh literal replay only after
the structural receipt is built. The independent replay is a test, not a
production seed. There are 6,390 exact checks, unchanged under optimized Python.

The constructive domain is strictly larger and permits unbounded run lengths
and exponents. Universal source entry remains open.
