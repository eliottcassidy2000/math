# Gold partition as a two-comparison defect: a failed greedy proof and its repair target

2026-10-14. **PROVED:** the exact heavy-child and signed-correlation
reformulations. **FINITE-EXACT / REFUTED:** the proposed maximally-balanced
first-comparison rule, and its conditional best-response repair, on the
explicit source below. **OPEN:** the Gold Partition and universal
1/3--2/3 conjectures. No novelty claim about the elementary reformulations
or the already published source poset is made.

## 1. Inheritance and proof plan

The closest repository mechanism is [THM-4503, the two-chain poset bridge
with source-height sidecar](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md).
It preserves a uniform extension law, but selecting integer sources can
change its weights and destroy the poset support. The hostile there is the
six-letter cell whose height cutoff selects only its first and last orders.
The corrected near miss is the 2026-09-26 crossings entry in
[MISTAKES](../../01-canon/MISTAKES.md): a balanced comparison gives a
log_(3/2) sorting bound, not a log_2 bound. The least-used coordinate in
the present session is the **conditional second-comparison distribution**.

Our board is: extension counts / conditional correlation / insertion fibres /
golden transfer matrices / geometric realizability / strict slack.
The anchor is a genuine gold-partition proof attempt; Borsuk4 is the
orthogonal lane and the eleven-state golden torus is the wildcard.

We use Peczarski's weak inequality from
[Comments on the Golden Partition Conjecture (2017)](https://cdm.ucalgary.ca/article/download/62349/46785/178291):
after two comparisons, the extension counts satisfy t0>=t1+t2 on every
possible branch. A chain after the first comparison has t2=1. The
[2010 workshop formulation](https://dwest.web.illinois.edu/regs/UL2010.pdf)
also states this inequality. A strict inequality would already fail at
the two-element antichain.

We accept the supplied/input literature, including
[Gupta's order-14 result](https://arxiv.org/html/2607.23926v1), without
repeating its census. Its low-balance examples are inputs to a new question:
does minimizing the first imbalance also certify the two-step budget?
The established order-ideal forward/backward counting method is credited
to De Loof, De Meyer and De Baets, as in Gupta Section3; our checker is
an independent implementation aimed at the new conditional statistic.

## 2. Only the heavier branch has an unpaid obligation

Let E=e(P)>1. Choose an incomparable first pair and orient its heavier
outcome as Q, with m=e(Q), n=E-m, and m>=n. Let

\[
 B(Q)=\max_{u\ne v}\min\{e(Q+uv),e(Q+vu)\}.
\]

For a chain set B(Q)=0. Define the integer defect

\[
 D_P(x,y)=m-n-B(Q).                                      \tag{1}
\]

**Proposition.** This first comparison can begin a gold certificate iff
D_P(x,y)<=0. Consequently P satisfies gold iff some incomparable pair
has nonpositive defect.

*Proof.* The light branch contains n extensions, and any second outcome
contains at most n. Thus t1+t2<=2n<=E there. In the heavy branch the smallest
achievable worst second outcome is m-B(Q), so the desired inequality is
m+(m-B(Q))<=m+n, which is exactly(1). If m=n=1, the stipulated terminal
value gives equality; the same formula holds. This proves both directions.

It also reconciles fixed and adaptive second comparisons: choose in advance
the pair needed on the heavier branch. The light branch accepts every
second pair, including an already decided comparison. No stronger adaptive
quantifier is needed for this particular two-step inequality.

The 1/3--2/3 consequence is immediate: B(Q)<=m/2, so a successful pair has
m-n<=m/2, or m<=2E/3. The converse needs the second-comparison condition.

If delta=min(m,n)/E and b(Q)=B(Q)/m, the exact required amplification is

\[
 b(Q)\ge g(\delta):=\frac{1-2\delta}{1-\delta}.          \tag{2}
\]

The map g has projective matrix [[-2,1],[-1,1]], determinant-1, and
eigenvalues 1/phi and -phi. Its fixed point in[0,1/2] is phi^(-2).
This is a genuine golden algebra arising from the budget inequality;
it does not prove that conditional poset balances satisfy that inequality.

## 3. An exact mixed-correlation target

Let X,Y be the signed outcomes (+1/-1) of two comparisons in a uniformly
chosen linear extension. Orient X so a=E[X]>=0, and put b=E[Y], r=E[XY].
Their four joint masses are

\[
 p_{st}=(1+sa+tb+st r)/4,\qquad s,t\in\{-1,1\}.
\]

The two comparisons give a gold certificate precisely when

\[
 \boxed{3a+|b+r|\le1.}                                \tag{3}
\]

Indeed the heavy branch has mass(1+a)/2 and its larger leaf has
mass(1+a+|b+r|)/4. The light branch is automatic by Section2. Reversing
the second comparison negates both b and r and leaves(3) invariant.
All these quantities are exact ratios of integers; no statistical estimate
is needed. Merely making a small leaves b+r uncontrolled.

The [Hexagon six-spin paper](https://hexagonmath.org/2610.00241v1) uses a
mixed-covariance mechanism to distinguish a triangle from biclique
corrections. We take that theme as inspiration for retaining r; none of its
spin-model theorems is asserted to apply to a poset's comparison variables.
Here equation(3) follows directly from the four cells.

This is the concrete remaining target for a proof: select two actual
comparison events whose means and mixed moment lie in(3). Arbitrary events
or arbitrary binary partitions are not permitted substitutes.

## 4. A genuine failed proof route, with an explicit escape

The first proposed strengthening was: *some maximally balanced first
comparison always starts a gold certificate*. A stronger version said
*every* maximizing comparison does. Both are false on the same14-element
poset. The source is the retained144/403 low-balance witness in
[Gupta's census-n14 data](https://github.com/agupta/gold-partition-conjecture/blob/main/data/census-n14.txt).
We used this published poset to measure a new conditional defect, not to
recheck the source's global census.

Its exact predecessor masks, four hex digits per labelled vertex, are

```text
0000 0000 0003 0001 000f 0009 003f 002f 00ff 00af 03ff 02ff 0fff 0aff
```

There are806 extensions. Its two best first pairs are(4,7) and(6,9), each
splitting them into518 and288. Every second comparison leaves at least322
extensions in some outcome of the heavier branch, because B(Q)=196.
Hence t1+t2>=518+322=840>806.

| First comparison, heavy orientation | Heavy m | Light n | B(Q) | Defect |
|---|---:|---:|---:|---:|
|4<7|518|288|196|34|
|6<9|518|288|196|34|
|0<1|525|281|244|0|

The last row succeeds: compare0 with1, then1 with3. Its four cell counts,
with the heavier first outcome in the first row, are

\[
 \begin{pmatrix}281&244\\281&0\end{pmatrix}.
\]

Thus525+281=806 and525+244<806 in the heavy branch; the light branch is
automatic. The less balanced first split buys a better second split.
There is no counterexample to either original conjecture here.

The corresponding signed data are
(a,b,r)=(244,318,-244)/806, saturating(3). For the failed first pair4<7
and its best second pair6<7 they are(230,-414,162)/806, so
3a+|b+r|=942/806>1. The four cells are[[196,322],[0,288]].

A natural local repair also fails. Map a first comparison to any best
comparison inside its heavier child, then try that comparison first in P.
The three comparison vertices form a closed best-response component:

```text
{4,7} -> {6,7} <- {6,9}
{6,7} -> {4,7} or {6,9}
```

The middle first comparison, oriented7<6, has defect126. All arrows remain
inside positive defects; none reaches the successful pair{0,1}. This is an
actual directed-graph obstruction to that recursive selector. A global
selection rule, an additional retained witness, or a justified escape move
is needed.

Deleting any single vertex from this example repairs the proposed greedy
property. This is deletion-minimality for this witness, not a claim that
fourteen is the smallest possible counterexample order.

## 5. Positive continuation and reproduction

The [context-transport packet](gold_partition_context_transport_20261014j.md)
proves a constructive way to preserve a *strict* certificate while adding
an arbitrary upper context and a sufficiently long buffer chain. The
insertion fibres have explicit bounded weight distortion; their slack pays
for it. Its zero-slack example fails for every buffer length, despite
convergence of normalized measures. This is a theorem for a constructed
unbounded family, not a proof that the original unbuffered source is gold.

The [golden order bridge](golden_order_bridge_20261014j.md) separately
identifies the eleven cyclic no-11 states with a boundary-conditioned
fence-ideal family and supplies the exact rank/edge-flow reweighting.
Its thin posets, with i<j exactly when j>=i+2, have an all-order Gold
certificate with cells [[F_(N-1),F_(N-2)],[F_(N-1),0]]. This is the same
equality-certificate shape as the successful [[281,244],[281,0]] above.
For the six-element thin poset, closing the five swap positions into a
cycle deletes exactly two orders from [[5,3],[5,0]], leaving
[[5,3],[3,0]] and breaking the Gold budget. The seam records precisely
which missing extensions destroyed the proof. The
[Borsuk interface](borsuk_order_interfaces_20261014j.md) retains the
geometric realization and near-diameter margin that a bare order forgets.

The finite exploratory universe is972 seeded random DAG draws on3..11
vertices,42 selected low-balance14-vertex inputs, and the specified
cover-deletion/comparison-addition perturbations of those inputs. It has
2,153 distinct labelled posets,2,144 nonchains. This is deliberately not
an unlabeled-poset census. Both greedy claims fail on the displayed source.
Every sampled poset nevertheless has a successful two-step certificate.

The [input JSON](gold_low_balance_inputs_20261014j.json) records the source
URL, full-source SHA256 and all42 exact masks, so reruns need no network.
The [output](gold_heavy_child_20261014j.json) preserves the source relations,
conditional counts, four cells and all selected hard-input results. Small
controls enumerate every permutation directly. The hostile's806 actual
linear extensions are independently generated by minimal-element recursion,
then used to count every comparison and every relevant heavy-child pair.

```text
python -B -X utf8 04-computation/experiments/gold_heavy_child_20261014j.py
python -B -O -X utf8 04-computation/experiments/gold_heavy_child_20261014j.py
```

The final normal and optimized runs perform **5,655 checks**. The
LF-normalized JSON SHA256 is
`11885a903a644eef4392deb4b48fa27ce45a65f8bdf1d8024c1c4a5882a984b2`.
Independent peer enumeration confirms the 806-order witness, its closed
best-response component, and the mathematical reformulations. The final
API validates exact source types before cache lookup; warmed-cache controls
reject Boolean/float aliases and malformed comparison endpoints.

The next proof obligation is **global certificate selection with conditional
weights retained**. The positive buffer theorem and negative best-response
component specify what such an induction must do: transport strict slack
when available, and handle exact equality without silently replacing it
by an asymptotic approximation.
