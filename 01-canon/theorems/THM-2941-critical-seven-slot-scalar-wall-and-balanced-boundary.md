---
id: THM-2941
title: "Critical seven-slot scalar wall, empty balanced boundary, and the A6 Gram/address state"
status: >
  PROOF CANDIDATE / SOURCE AND CANON AUDIT PENDING.  On every one of the
  3,003 literal six-body carriers, the seven-slot pair-Hunter scalar has
  first crossing exactly h/7: q7>h/7 and B2>2h/7 uniformly, so it supplies
  no discrepancy-finite first-centre core.  Nevertheless a literal
  pointwise seven-tail cover cannot have zero multiplicity excess.  The
  longest carrier component would have to belong to one open tooth, forcing
  an external owner w<=6 although w>=15.  The hypothetical balanced boundary
  is the rank-six regular-simplex Gram (h/49)(7I-J).  Arbitrary
  positive-excess covers and LRC(14) remain open.
source: root-lrc14-j7-critical-wall-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2923-complete-seven-body-six-slot-recursive-pair-hunter-closure
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
related:
  - THM-856-hunter-tree-seven-comb-crossing
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-2184-two-scale-tail-continuation-profile
  - THM-2893-complement-cap-finite-core-flag-lemma
verification:
  - 04-computation/lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.py
  - 05-knowledge/results/lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.out
  - 04-computation/lrc14_j7_critical_scalar_wall_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_critical_scalar_wall_independent_thm2941.out
  - 04-computation/lrc14_j7_top7_overlap_graph_scout_thm2941.py
  - 05-knowledge/results/lrc14_j7_top7_overlap_graph_scout_thm2941.out
---

# THM-2941 -- critical seven-slot scalar wall and empty balanced boundary

**PROOF CANDIDATE / SOURCE AND CANON AUDIT PENDING.**

This theorem identifies what the successful THM-2923 recursion loses one rung
earlier.  The complete scalar state `(q_1,...,q_7,B_2)` stops at its exact
critical density on every six-body root.  The critical boundary is not a
counterexample: pointwise topology makes it empty.  What can still escape is a
sequence of positive-excess covers approaching that empty boundary.  Its
natural state is a `7 x 7` restricted-overlap Gram deformation of the
rank-six boundary, joined to the endpoint-grid owner/transition word.

Nothing here proves the six-body/seven-tail rung, the sector with fewer body
speeds, or unrestricted LRC(14).

## 1. Literal six-body carrier

On `T=R/Z`, put

```text
D_w={t in T: ||wt||<1/14}.
```

Fix

```text
E in C({1,...,14},6),
C_E=T minus union_(e in E) D_e,
h=mu(C_E)>0.                                             (1)
```

The allowed external labels are the distinct integers `w>=15`.  Write

```text
c(w)=mu(C_E intersect D_w)
```

and let

```text
q_1>=...>=q_7
```

be the exact global top seven allowed singleton coverages.  Let

```text
B_2=max_(15<=u<v) mu(C_E intersect (D_u union D_v)).       (2)
```

The computations below prove that these are maxima over the infinite label
set, not maxima inside an arbitrary search box.

## 2. General critical scalar lemma

For any mass `h`, nonincreasing ranks `q_1,...,q_7`, and pair cap `B_2`,
define

```text
G_7(a)=a+sum_(j=2)^7 min(a,q_j,B_2-a),
0<=a<=min(q_1,B_2).                                      (3)
```

Suppose

```text
q_7>=h/7,                     B_2>=2h/7.                  (4)
```

Then

```text
G_7(h/7)=h.                                               (5)
```

Indeed, all seven terms in `(3)` equal `h/7`.  If `a<h/7`, every term is at
most `a`, so

```text
G_7(a)<=7a<h.                                             (6)
```

Consequently the first hostile level is exactly

```text
lambda_7=min{a:G_7(a)>=h}=h/7.                            (7)
```

This is a no-go for the scalar certificate, not evidence for a cover.
THM-735 only bounds

```text
c(w)<h/7+gamma/w.
```

At a threshold strictly above `h/7` this gives a finite core; at `(7)` it
does not.

## 3. Infinite-tail seals

If `C_E` has `r` interval components, THM-735 gives the explicit strict tail

```text
c(w)<h/7+gamma/w,                   gamma=99r/490.         (8)
```

For a finite head through `M`, put

```text
tau=h/7+gamma/(M+1).                                      (9)
```

If its seventh rank satisfies `q_7>=tau`, every omitted singleton lies
strictly below the retained top seven.

For the pair cap, let `H_2` be the largest pair union in the finite head and
`q_1` its largest singleton.  The verifier accepts a pair seal only when

```text
H_2>=q_1+tau.                                             (10)
```

A head-tail pair is then at most `q_1+tau`.  Subadditivity also gives

```text
H_2<=2q_1,
```

so `(10)` implies `tau<=q_1`.  Therefore a pair with two omitted endpoints is
at most

```text
2tau<=q_1+tau<=H_2.                                      (11)
```

Thus `(10)` seals the actual global maximum `(2)`, including the
two-omitted-label case.

## 4. Exact all-root scalar census

The primary rational-interval verifier and an independent integer-ruler
implementation both reconstruct all

```text
binom(14,6)=3,003
```

carriers.  Every root satisfies the strict form of `(4)` and hence `(7)`.
There are no scalar closures and no positive hostile-level gaps.

The exact extrema are

```text
min(q_7-h/7) = 11/2548
  at E=(1,2,3,4,6,13),
max(q_7-h/7) = 26077/1284192
  at E=(4,6,9,12,13,14);

min(B_2-2h/7) = 1756/77175
  at E=(1,3,5,7,11,14),
max(B_2-2h/7) = 203624/2207205
  at E=(2,4,8,10,13,14);

min(max_a G_7(a)-h) = 6847/120120
  at E=(1,2,3,4,6,13),
max(max_a G_7(a)-h) = 200138/945945
  at E=(2,4,6,10,12,14).                                (12)
```

The primary verifier's smallest **finite-horizon** pair seal margin

```text
H_2-(q_1+tau)
```

is `1/21351330`.  The independent verifier separately records the smallest
**asymptotic** margin

```text
H_2-q_1-h/7=1613/3783780.
```

These are different quantities.  The maximum pair horizon is `6,635`.  The
primary implementation starts every top-seven scan at `2,000`; the
independent dynamic implementation proves that `717` is the largest actually
required top-seven horizon.

## 5. Multiplicity excess and the exact Gram deformation

Now assume, only in this section, that seven allowed labels give a literal
pointwise cover of `C_E`.  Put

```text
A_i=C_E intersect D_(w_i),
c_i=mu(A_i),
delta_i=c_i-h/7,
m(t)=sum_i 1_(A_i)(t),
Delta=sum_i delta_i.                                     (13)
```

Because the cover is pointwise, `m>=1` on `C_E`, and

```text
Delta=integral_(C_E)(m-1)>=0.                            (14)
```

If every `c_i<=h/7+epsilon`, then

```text
0<=Delta<=7epsilon,
-6epsilon<=delta_i<=epsilon,
mu{m>=2}<=Delta,
sum_(i<j)mu(A_i intersect A_j)
  =integral binom(m,2)<=(7/2)Delta,
integral_(C_E)(m-1)^2<=6Delta.                           (15)
```

Thus a near-critical cover is quantitatively an approximate measurable
partition.

There is an exact covariance form.  In `L^2(C_E)`, put

```text
f_i=1_(D_(w_i))-1/7,
p_ij=mu(C_E intersect D_(w_i) intersect D_(w_j)).
```

Direct expansion gives

```text
<f_i,f_i> = 6h/49+5delta_i/7,
<f_i,f_j> =-h/49+p_ij-(delta_i+delta_j)/7,     i!=j.     (16)
```

For positive excess this `7 x 7` Gram matrix may have rank seven:
`sum_i f_i=m-1` need not vanish.  Rank six belongs to the balanced partition
boundary below, or to an explicitly projected `1^perp` quotient; it is not a
property of every near state.

If a hypothetical cover were both exact and balanced, then

```text
delta_i=0,             p_ij=0,
Gram(f_1,...,f_7)=(h/49)(7I-J).                          (17)
```

This is the Laplacian of `K_7` scaled by `h/49`: it has rank six, kernel the
all-ones line, and realizes the regular `A_6` simplex.  Equation `(17)` is a
description of a hypothetical boundary.  It does not assert that balance
alone gives a cover, and the next section proves that no literal allowed
cover can attain this boundary.

## 6. The zero-excess boundary is empty

The danger combs are open, so `C_E` is a finite union of closed intervals.
All body danger combs contain a neighborhood of zero; hence every
positive-length component has a nonwrapping closed representative.

Suppose the pointwise cover in Section 5 had `Delta=0`.  Equation `(14)` and
the integer-valued inequality `m>=1` imply

```text
m=1 almost everywhere.                                  (18)
```

On a positive-length carrier component `I`, every pair

```text
I intersect D_(w_i) intersect D_(w_j)
```

is relatively open.  By `(18)` it has measure zero, so it is empty.  The
seven relatively open owner sets are therefore a disjoint pointwise cover of
the connected set `I`.  One owner must contain all of `I`.

If `I` has length `ell` and lies in `D_w`, connectedness puts it inside one
open tooth of length `1/(7w)`.  Since `I` is closed,

```text
ell<1/(7w),                     w<1/(7ell).              (19)
```

The exact all-root component census gives

```text
min_E longest_component(C_E)=23/1092
  at E=(1,6,7,8,10,13).                                  (20)
```

Thus the owner of a longest component would obey

```text
w<1092/161<7,
```

so `w<=6`, contradicting the allowed-label condition `w>=15`.  The
per-root integral owner-bound histogram is

```text
bound   1     2      3     4    5   6
roots  39  1,601  1,049  283   28   3.                  (21)
```

Therefore every hypothetical literal six-body/seven-tail cover has

```text
Delta>0,                         max_i c_i>h/7.           (22)
```

If `w_*` realizes the largest singleton excess, then

```text
delta_*>=Delta/7.
```

Combining this with the strict tail `(8)` gives the packetwise bound

```text
w_*<7gamma/Delta.                                        (22a)
```

There is a complementary transition-width form of the same obstruction.
On a longest component `I`, no one allowed owner contains all of `I`, by
`(19)`--`(20)`.  The relatively open owner sets cover the connected interval
`I`; if all pair intersections were empty, they would disconnect it.
Therefore some pair intersection contains a nonempty relatively open
subinterval `J`.  Put `kappa=mu(J)>0`.  Since `m>=2` on `J`,

```text
Delta>=kappa,             max_i delta_i>=kappa/7,
w_*<7gamma/kappa.                                      (22b)
```

Thus every counterexample packet carries a positive owner-transition width
and a correspondingly bounded excess owner.  This is still not a uniform
bound: `kappa` depends on the packet and may tend to zero along an escaping
sequence.  The only possible escape is therefore an endpoint-transition
collapse, rather than an exact disjoint partition.

The strictness in `(22)` is packetwise, not yet uniform.  A sequence of
putative covers could still have `Delta` tend to zero while its labels
escape.  Quantifying or classifying that escape is the open problem.

## 7. What pair compatibility recovers

As a scoped control, the second verifier retains the complete labelled
restricted-overlap graph on each root's seven actual top-ranked singleton
labels.  The Hunter lower bound is

```text
h-sum_i c_i+max_T sum_(ij in T)p_ij.                     (23)
```

Among the `3,003` controls:

```text
positive three-edge matching margin          47,
positive best-star margin                   256,
positive maximum-spanning-tree margin     2,200,
nonpositive maximum-tree controls            803.        (24)
```

Every one of the `3,003` literal top-seven unions nevertheless has positive
complement.  Hence compatible pair incidence is real information—the tree
recovers `1,944` controls missed by every star—but maximum-tree Hunter remains
incomplete on `803` of these actual-top-seven controls.  This is only a test
on the actual top-seven labels; it is not a reduction of arbitrary
seven-label packets to those labels, nor a proof that every possible
pair-graph functional fails on literal LRC packets.

There is also an abstract hostile control for any attempt to make first and
second moments equivalent to coverage.  On the even- and odd-parity laws on
`{0,1}^3`, the three coordinate events have the same singleton masses `1/2`
and pair masses `1/4`.  Their union masses are respectively `3/4` and `1`.
Padding by four empty labelled events gives the same seven-event Gram data
with different union predicates.  This is an information-loss example, not
a claim that either law is realized by an LRC carrier.

## 8. The multi-slope Gram/address frontier and a uniform finite sector

Let `L` be the endpoint-grid denominator of the body carrier and write a tail
speed as

```text
w=La+b,                       0<=b<L.                    (25)
```

THM-2928 proves that a literal six-body carrier cannot be covered when all
seven residues `b` vanish, or when exactly one is nonzero.  Combining it with
`(22)` yields:

> Every remaining hypothetical literal six-body/seven-tail cover has
> positive multiplicity excess and at least two nonzero endpoint-grid drift
> residues.

For exactly two drifts, THM-2928 proves that every chart with fixed body
carrier and aligned multiplier set reduces to a finite exact pair-clock
quotient.  Its pointwise argument bounds the smaller slope by `585/154` and
the larger by `13 max A`; its load-bearing extra coordinate is the
carrier-local pair endpoint current, not the global full-circle overlap.
This does not prove those charts empty.

There is, however, a uniform **finite reduction** that does not require
bounding the aligned multiplier set.  Suppose there are `k` aligned tails
and `r=7-k` drifts, with `2<=k<=6`.  The component-free tariff
THM-2928(25h) supplies a drift

```text
z<C_k L,                    C_k=(7-k)/b_k,              (25a)
```

where `L` is the fixed canonical denominator of the body carrier and `b_k`
depends only on `k`.  Hence only finitely many integral `z` occur for each
of the `3,003` body roots.  Delete such a `D_z`.  Six tail labels remain,
with no alignment assumption needed.  At every proper node of their literal
residual tree, the body speeds and already chosen tails total at most twelve
speeds.  Settled `LRC(<=13)` gives a point at clearance at least `1/13`;
the strict margin `1/13-1/14` gives a positive interval in that residual.
It therefore satisfies the positive-mass interval hypothesis of the
cap-free `p<=6` first-apex recursion THM-2893(7a)--(7b).  Recursing to depth
at most six produces a finite exact decision tree.

With

```text
L_E=lcm(14e:e in E)=14 lcm(E),
max_E L_E=5,045,040,
```

the strict integral first-drift caps can even be chosen globally as

```text
k                         2           3          4          5         6
z <=              281,080,799  74,954,879 42,988,828 19,164,599 7,665,839.
                                                               (25b)
```

These large numbers are reduction bounds, not a completed finite census.
Together with THM-2928's empty `k=7` branch, this proves:

> Every six-body/seven-tail branch with at least two aligned tails is either
> already empty or uniformly reducible to a finite exact decision tree.
> Consequently any sector not yet known to admit such a finite reduction
> has at most one aligned tail, hence at least six drifts.

The next faithful object must therefore join

```text
node data:
  delta_i, residue b_i, normalized slope a_i+b_i/L;
edge data:
  restricted overlap p_ij, gcd/reduced ratio, tooth indices;
component data:
  ordered owner word, transition positions and widths;
hyperedge data:
  the endpoint relation producing the width kappa in (22b).              (26)
```

If tooth endpoints are

```text
(14k+sigma)/(14u),       (14l+tau)/(14v),
sigma,tau in {-1,+1},
```

their signed separation is exactly

```text
((14k+sigma)v-(14l+tau)u)/(14uv).                        (27)
```

Thus vanishing transition widths are not anonymous analytic errors: they
carry an integer endpoint relation.  Equations `(16)`, `(25)`, and `(27)`
define the multi-slope Gram/address transition ledger to be classified.

A tournament on the seven labels is not an equivalent quotient.  It can
orient a chosen owner precedence, but it loses symmetric overlap magnitudes,
component order, repeated transitions, ties, endpoint signs, and the
multi-residue clock.  The lawful tournament-related object is a transition
path/graph decorated by the Gram and address sidecars.

Likewise, permuting equal-mass disjoint owners among carrier components
preserves all singleton and pair masses while changing the component owner
word.  This second abstract control isolates the address loss from the
higher-moment loss in the parity example.  Neither control asserts literal
LRC realizability.

## 9. Scope and audit state

The exact scalar census and component bound have two independent
implementations.  The top-seven overlap graph is a scoped structural scout.
Promotion still requires frozen-source ordinary/optimized replay and final
canon audit.

This theorem does not give a uniform lower bound for `Delta` or `kappa`,
turn the `803` nonpositive actual-top-seven tree margins into certificates,
handle arbitrary packets, run the new finite decision trees, classify the
zero/one-aligned multi-drift address hypergraph, close the
six-body/seven-tail rung, or prove LRC(14).
