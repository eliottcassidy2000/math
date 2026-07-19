---
id: THM-1135
title: The r=6 unbounded horn has a harmonic-discrepancy tail and a finite ratio box
status: PROVED as a uniform finite reduction.  The exact one-interval discrepancy gives an order-free harmonic tail; a dependency-free Fraction referee checks all 792 core measures and component counts; settled LRC for proper prefixes gives the recursive ratio horn.  Every r=6 killer tuple outside the displayed finite box is certified.  This is NOT an r=6 closure: the finite box has not been exhaustively verified, although THM-1121 closes its subbranch 92<=k_i<=332.
source: codex-2026-07-18-S67 (r6 tail salvage)
depends_on:
  - THM-933    # exact singleton primitive discrepancy / interval block-gluing interface
  - THM-955    # prior exact 6/(49k) danger-comb discrepancy and cluster-gap application
  - THM-1121   # universal weighted atlas for the bounded 92..332 subbranch
  - THM-1134   # refutation of the order-sensitive constant-max-T extrapolation
  - LRC(<=13)  # settled proper-prefix clearance used in the ratio horn
script: 04-computation/r6_harmonic_tail_finite_box_codex_S67.py
output: 05-knowledge/results/r6_harmonic_tail_finite_box_codex_S67.out
---

# THM-1135 — the harmonic tail and finite ratio box

## Statement

Let `P` be a seven-element subset of `{1,...,12}`, put `M=max(P)`, and let

\[
 13M<k_1<k_2<\cdots<k_6
\]

be the ordered killers.  Then the thirteen-speed family `P union {k_1,...,k_6}`
has a time at which every speed is at distance at least `1/14` from an integer
whenever either of the following holds.

1. `k_1>=514` or `k_2>=951`.
2. For some `s in {2,3,4,5}`,
   \[
     \frac{k_{s+1}}{k_s}>\frac{6(8+s)}{s+1}.
   \]

Consequently, every tuple not certified by these horns lies in the finite box

```text
k1 <=       513
k2 <=       950
k3 <=    19,000
k4 <=   313,500
k5 <= 4,514,400
k6 <=58,687,200.
```

This theorem converts the unbounded `r=6` tail into a finite problem.  It does
**not** assert that this box has been enumerated, and therefore does not close
the `r=6` branch.

## 1. The exact one-interval loss

For a speed `k`, write

\[
 D_k=\{t\in\mathbb R/\mathbb Z:\|kt\|<1/14\}.
\]

The function `1_{D_k}` has mean `1/7`.  On one period of length `1/k`, the
primitive of `1_{D_k}-1/7` rises with slope `6/7` through a danger tooth of
length `1/(7k)` and falls with slope `-1/7` through the complementary interval
of length `6/(7k)`.  Its oscillation is therefore exactly

\[
 q_k=\frac6{49k}.
\]

Thus every circular interval `I` satisfies

\[
 |I\cap D_k|\le \frac{|I|}{7}+\frac6{49k}.                 \tag{1}
\]

This is the singleton specialization of THM-933's primitive-discrepancy
interface, but the two-slope proof above is self-contained.

If `E` is a disjoint union of `C` circular intervals and has measure `mu`,
summing (1) gives

\[
 |E\cap D_k|\le \frac\mu7+\frac{6C}{49k}.                 \tag{2}
\]

Apply the union bound to six danger sets.  A point of `E` survives provided

\[
 \begin{aligned}
 |E\setminus\textstyle\bigcup_{i=1}^6D_{k_i}|
 &\ge \mu-\sum_{i=1}^6|E\cap D_{k_i}|\\
 &\ge \frac\mu7-\frac{6C}{49}\sum_{i=1}^6\frac1{k_i}>0.
 \end{aligned}
\]

Equivalently, the **harmonic tail criterion** is

\[
 \boxed{\displaystyle
   \sum_{i=1}^6\frac1{k_i}<B(P):=\frac{7\mu(P)}{6C(P)}.}   \tag{3}
\]

Here `E=S(P)` is the exact core-safe set, `mu(P)=|S(P)|`, and `C(P)` is its
number of components.  Unlike the failed max-`T` extrapolation, (3) is
order-free and global in scale.  Its tail variable is the reciprocal-loss
ledger, not a constant upper bound on an interval threshold.

This also reconciles THM-1134's infinite translated obstruction with a valid
tail theorem.  Along that family the sequential threshold `T` tends to
infinity, but any six-killer completion with smallest killer at least `514` is
immediately dispatched by (3).  Unbounded `T` and a uniform LRC tail are
therefore compatible; the scalar `T` was the wrong global coordinate.

## 2. Exact 792-core constants

The dependency-free verifier constructs `S(P)` from its rational endpoints
`j/p +/- 1/(14p)` for every one of the `binom(12,7)=792` cores.  It does this
twice, once by breakpoint/midpoint classification and independently by
unioning all clipped danger teeth and taking the complement; the exact
component lists agree on every core.  The smallest budget is attained at

\[
 P_*=\{1,2,7,9,10,11,12\},\qquad
 C(P_*)=26,\qquad
 \mu(P_*)=\frac{12559}{48510},
\]

and equals

\[
 B(P_*)=\frac{12559}{1081080}.
\]

The killers are distinct ordered integers, so `k_1>=K` gives the sharper
estimate

\[
 \sum_{i=1}^6\frac1{k_i}\le\sum_{j=K}^{K+5}\frac1j.
\]

At the exact worst core the endpoint comparison is

\[
 B(P_*)-\sum_{j=514}^{519}\frac1j
 =\frac{114421278347}{370205035514594280}>0.             \tag{4}
\]

The same estimate at `K=513` is nonpositive, so 514 is the sharp integer
threshold delivered by this uniform comparison.  Thus (3)-(4) prove the first
large-tail claim.

There is a useful second exact consequence.  The killer convention gives
`k_1>=13M+1`.  If `k_2>=951`, distinctness gives

\[
 \sum_{i=1}^6\frac1{k_i}
 \le \frac1{13M+1}+\sum_{j=951}^{955}\frac1j.
\]

The verifier compares this expression with all 792 exact core budgets.  Its
smallest margin is again at `P_*` (where `13M+1=157`) and equals

\[
B(P_*)-\frac1{157}-\sum_{j=951}^{955}\frac1j
=\frac{4670546220583}{4412023437154312980}>0.            \tag{5}
\]

The analogous margin starting at 950 is nonpositive.  Hence (3) applies at
`k_2>=951`, and every harmonic-uncertified tuple has

\[
 k_1\le513,\qquad k_2\le950.                              \tag{6}
\]

No covering or divisibility assumption is used in (1)-(6).

## 3. Proper prefixes absorb every later killer at once

Fix `s in {1,...,5}` and let

\[
 A_s=P\cup\{k_1,\ldots,k_s\}.
\]

This prefix has `7+s<=12` speeds.  Settled LRC for proper prefixes supplies a
time `t_0` with

\[
 \min_{v\in A_s}\|vt_0\|\ge\frac1{8+s}.
\]

The minimum-clearance function is `k_s`-Lipschitz.  Hence the prefix remains
safe at level `1/14` on a circular interval of length

\[
 L_s\ge \frac2{k_s}\left(\frac1{8+s}-\frac1{14}\right)
       =\frac{6-s}{7(8+s)k_s}.                             \tag{7}
\]

There are `m=6-s` remaining killers.  Applying (1) on this single interval,
using `k_j>=k_{s+1}` for every `j>s`, leaves positive measure whenever

\[
 L_s-\frac{mL_s}{7}-\frac{6m}{49k_{s+1}}>0.
\]

Substituting (7) and cancelling `m=6-s` yields the clean ratio condition

\[
 \boxed{\displaystyle
   \frac{k_{s+1}}{k_s}>\frac{6(8+s)}{s+1}.}                \tag{8}
\]

The five thresholds are

```text
s=1: 27,   s=2: 20,   s=3: 33/2,   s=4: 72/5,   s=5: 13.
```

The `s=1` bound is weaker than the exact harmonic cap (6).  If none of the
remaining `s=2,...,5` conditions applies, recursively using (6) and the
negations of (8) gives

\[
\begin{aligned}
k_3&\le20(950)=19000,\\
k_4&\le\frac{33}{2}(19000)=313500,\\
k_5&\le\frac{72}{5}(313500)=4514400,\\
k_6&\le13(4514400)=58687200,
\end{aligned}
\]

which proves the finite box.

## 4. The explicit covering max-T failure is easy by another chart

The failure requested by the tail audit is the genuinely covering family

\[
 V=\{1,2,4,7,9,11,12,290,292,294,296,298,338\}.
\]

It covers every modulus `2,...,14`.  After removing the first five killers,
the exact residual data are

\[
 N=370,\qquad
 \mu=\frac{214450238441}{1981564300485},\qquad
 L=\frac1{1043},\qquad
 T=\frac{1043}{3}>338.
\]

Thus THM-1102's last-killer sufficient condition really does fail in the
covering branch; replacing `333` by the next `13`-multiple `338` is not merely
the non-covering issue isolated in THM-1134.

Nevertheless the particularly simple exact rational time

\[
 t=\frac5{16}
\]

has

\[
 \min_{v\in V}\|vt\|=\frac18>\frac1{14},
\]

with equality owners `v=294,298`.  The originally supplied witness is also
checked exactly:

\[
 t=\frac{106}{303},\qquad
 \min_{v\in V}\|vt\|=\frac{15}{101},
\]

with equality owners `v=9,294`.  The denominator-16 chart makes especially
clear that the displayed family is a warning about **certificate selection**,
not a hard LRC instance.  Its scale belongs inside the finite box above, where
arithmetic charts must complement the analytic horn.

## 5. What changed in the underlying object

The failed route asked for the largest residual interval after five ordered
removals and scalarized that complex to `T`.  This loses both permutation
symmetry and the fact that six danger sets each have a tiny, exactly known
primitive discrepancy on every core component.  The replacement object is a
two-layer loss ledger:

- core components carry baseline measure;
- a danger-set edge has bulk cost `mu/7` and boundary cost `6/(49k)` per
  component;
- scale gaps are handled by restarting on a proper-prefix safe interval.

The harmonic ledger preserves all six killers symmetrically.  The prefix
restart then uses only where symmetry is genuinely broken: at a large scale
jump.  This explains why the two bounds compose cleanly.

## Tournament Analysis and assumption challenge

On the six danger-set vertices, orient an edge by the harmonic observable
`1/k_i-1/k_j`, with integer label as the tie gauge.  The tournament is
transitive: score histogram `0,1,2,3,4,5`, no directed cycles, six singleton
SCCs, and one Hamiltonian path (`290 -> 292 -> 294 -> 296 -> 298 -> 338`) for
the explicit family.  It preserves the scalar loss ordering used in (3), but
destroys component owners, phases, and every higher danger-set intersection.

Thus the tournament is only a diagnostic quotient.  The proof-facing
vertices are better taken to be **core components and danger-set obligations**;
their weighted incidence ledger, not the runner tournament, proves the tail.
This explicitly challenges the assumption that either runners or arcs must be
the vertices of the useful combinatorial object.

## Reproducibility

The referee uses `fractions.Fraction` throughout, checks all 792 cores, freezes
the exact extrema (4)-(5), verifies the ratio-box arithmetic, recomputes the
covering predicate for the displayed family, reproduces its exact max-`T`
failure, and verifies all thirteen distances at both `5/16` and `106/303`.
Every proof-bearing check uses an explicit exception rather than Python
`assert`, so verification remains active under `python -O`; normal and `-O`
outputs are byte-identical to the frozen artifact.  SHA-256:

```text
script  192e6b3eb3ecac895345611bea315d39512ebe952771db4e4370f1cfb1dd581b
output  dd75bc8fd8dad1f3a2feac9a2cfeb9580edd9bce12ba86ed97e67459b6dad8c7
```

## Correct frontier

THM-1121 closes the universal subbox in which all six killers lie in
`[92,332]`.  This theorem removes infinity but leaves a large mixed-scale
finite box.  The next target is now exact and honest:

1. extend the weighted obligation atlas to the two bounded leading coordinates
   `k1<=513`, `k2<=950`, allowing later charts to depend on their residues;
2. use (8) to stratify only the bounded adjacent-ratio regions, rather than
   scanning arbitrary unbounded quintuples;
3. treat high harmonic load by arithmetic charts such as the `5/16`
   witness, and low harmonic load by (3);
4. formalize (1), (3), and (8); they are short reusable analytic lemmas and do
   not require reflecting the 792-core census into the kernel immediately.
