---
id: THM-2941
title: "Critical seven-slot scalar wall, projected aligned-sector closure, and the A6 boundary"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.  On all 3,003
  literal six-body carriers the seven-slot pair-Hunter scalar first crosses
  exactly at h/7, while a literal zero-excess cover is impossible and its
  hypothetical Gram is (h/49)(7I-J).  Aligned safe surplus bounds the first
  drift absolutely for k>=2.  The lossless projected residual closes the
  five-aligned/two-drift face independently of THM-2928, and exact suffix
  filters make k=2,3,4 uniformly finite-reducible; two frontier addenda
  improve the k=3 first-drift cap from 380 to 378.  THM-2928's later
  divisor-status/local-needle chain empties k=4, so only k=2,3 remain
  finite-but-uncensused; k=5,6,7 are also empty.  The zero/one-aligned
  sector, the remaining finite censuses, the full
  six-body/seven-tail rung, and LRC(14) remain open.  Verification is
  internal exact computation and proof audit; there is no Lean or external
  peer-review claim.
source: root-lrc14-j7-critical-wall-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-1094-exact-two-comb-component-theorem
  - THM-1234-sharp-five-comb-compatibility-floor
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2923-complete-seven-body-six-slot-recursive-pair-hunter-closure
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - LRC(<=13)
related:
  - THM-856-hunter-tree-seven-comb-crossing
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-2184-two-scale-tail-continuation-profile
verification:
  - 04-computation/lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.py
  - 05-knowledge/results/lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.out
  - 04-computation/lrc14_j7_critical_scalar_wall_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_critical_scalar_wall_independent_thm2941.out
  - 04-computation/lrc14_j7_top7_overlap_graph_scout_thm2941.py
  - 05-knowledge/results/lrc14_j7_top7_overlap_graph_scout_thm2941.out
  - 04-computation/lrc14_j7_aligned_projected_arc_suffix_thm2941.py
  - 05-knowledge/results/lrc14_j7_aligned_projected_arc_suffix_thm2941.out
  - 04-computation/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_frontier_fibre_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_frontier_fibre_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_next_frontier_scalar_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_next_frontier_scalar_closure_thm2941.out
  - 04-computation/lrc14_j7_k3_next_frontier_scalar_independent_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_next_frontier_scalar_independent_thm2941.out
  - 04-computation/lrc14_j7_k3_uniform_ray_status_closure_thm2941.py
  - 05-knowledge/results/lrc14_j7_k3_uniform_ray_status_closure_thm2941.out
---

# THM-2941 -- critical scalar wall, projected aligned closure, and A6 boundary

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY AUDITED.**

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

The word **pointwise** is essential.  An almost-everywhere partition of
strict-open teeth may have an uncovered exit/entry seam, but an LRC
counterexample may not: such a seam belongs to neither open tooth and is a
literal safe point unless a further tooth covers it.  Thus the zero-excess
case has no endpoint handoff or gcd-capacity switching graph.  Its connected
carrier components are chamber-locked to single owners.

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

Thus every counterexample packet carries a positive owner-overlap width and a
correspondingly bounded excess owner.  This is still not a uniform bound:
`kappa` depends on the packet and may tend to zero along an escaping sequence.
That limiting degeneration is an overlap-width collapse, not a lawful
strict-open handoff.  At every actual packet an endpoint seam must remain
protected by another open tooth; an unprotected limiting seam is itself a
safe point.

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
By itself this does not prove those charts empty; the lossless projected
residual below supplies the missing terminal.

The multiplicity deformation in fact gives an **absolute** first apex
throughout the `k>=2` aligned sector.  Suppose there are `k` aligned tails,
`d=7-k` drifts, and let `m_A` and `m` count the active aligned and total tail
combs.  Pointwise on a cover,

```text
m-1 >= (m_A-1)_+.                                      (25a)
```

If `u_A` is the normalized safe mass of the `k` aligned multipliers, Boolean
grid tensorization gives

```text
Delta=integral_(C_E)(m-1)
 >=h integral_T(m_A-1)_+
 =h[k/7-mu(union_(a in A)D_a)]
 =h[u_A-d/7].                                          (25b)
```

The safe floors in THM-2928 therefore give

```text
k                 2       3         4          5          6
d                 5       4         3          2          1
u_k             66/91   55/91    558/1183   478/1365   61/273
eta_k:=u_k-d/7   1/91    3/91     51/1183    88/1365   22/273. (25c)
```

Here `u_A>=u_k`, so `(25b)` gives `Delta>=h eta_k`.  The positivity begins
exactly at `k=2`: for `k=1`, `u_A=6/7=d/7`, so this excess mechanism
vanishes.  This explains, rather than merely observes, why the zero/one-
aligned sector is the remaining infinite frontier.

Each aligned singleton has restricted mass `h/7`.  Ordering the drifts
`z_1<...<z_d` and writing `delta(z)=c(z)-h/7`,

```text
Delta=sum_(q=1)^d delta(z_q)
 <=(6r_E/49)sum_(q=1)^d 1/z_q
 <=6d r_E/(49z_1).                                    (25d)
```

Here `r_E` counts positive-length carrier components, and THM-1094 supplies
the componentwise discrepancy inequality.  The two independent carrier
reconstructions give

```text
max_E r_E/h_E=3993990/32029
```

at `E=(1,10,11,12,13,14)`, where
`h_E=32029/105105` and `r_E=38`.  Combining `(25c)`--`(25d)` gives the
body-uniform integral caps

```text
k                         2      3      4      5     6
smallest drift z_1 <=   6,947  1,852  1,062    473   189.       (25e)
```

There is a complementary lower bound on the largest drift scale.  Put

```text
M=max(E union {z_1,...,z_d}).
```

The `6+d=13-k` nonaligned/body speeds have, by settled `LRC(<=13)`, a point
with clearance at least `1/(14-k)`.  Their distance functions are
`M`-Lipschitz, so the closed arc `I` of radius

```text
R=k/[14(14-k)M]
```

around that point is safe at level `1/14`.  Let

```text
phi_L(t)=Lt mod 1,                 P=phi_L(I).
```

The body and drift teeth miss `I`, while each aligned tooth pulls back from
its normalized danger set.  A cover therefore forces

```text
P subset U_A=union_(a in A)D_a.
```

The interval-image formula gives `mu(P)=min(1,2LR)`.  Moreover `P` is
compact and `U_A` is open and proper, so the inclusion has strict measure:

```text
min(1,2LR)<mu(U_A)=1-u_A<=1-u_k.
```

Consequently every such cover must satisfy the **projected-safe-arc wall**

```text
M > alpha_k L,              alpha_k=k/[7(14-k)(1-u_k)],

k                         2        3          4           5          6
alpha_k               13/150   13/132   2366/21875   2275/18627   819/5936.
                                                               (25f)
```

The earlier whole-cell test is the weaker consequence
`M>kL/[14(14-k)]`.  Thus the aligned sector is squeezed between a bounded
first drift and a quantitatively forced second scale; it is not merely
placed in an undifferentiated finite box.

The arc is only a lower bound for a more faithful quotient.  For any fixed
drift packet `Z`, put

```text
S_(E,Z)=C_E minus union_(z in Z)D_z,
P_(E,Z)=phi_L(S_(E,Z)).
```

For a fixed aligned multiplier set `A`, the original tail-cover predicate is
equivalent to

```text
{D_z:z in Z} union {D_(La):a in A} covers C_E
 iff P_(E,Z) subset U_A.                                      (25g)
```

Indeed `D_(La)=phi_L^(-1)(D_a)`, and the points outside `S_(E,Z)` are
already covered by a drift.  Since `S_(E,Z)` and its image are compact while
`U_A` is open and proper, a completion in particular forces

```text
mu(P_(E,Z))<1-u_A<=1-u_k.                                 (25h)
```

Thus `P_(E,Z)` is a multiplier-free, lossless reduction of the aligned
completion problem, not merely a necessary statistic.  Equation `(25f)`
uses only one guaranteed arc inside `S_(E,Z)`, whereas a finite drift census
can retain the whole projected set, its components, and its endpoint
address.

There is also a finite clause representation.  Let `J_E` be the body-safe
`1/L` cells and, in the normalized coordinate `u` on cell `j`, put

```text
E_z(j)={u in [0,1]: ||z(j+u)/L||<1/14}.
```

De Morgan's law gives, modulo finitely many endpoints and therefore exactly
in Lebesgue mass, the interval identity

```text
T minus P_(E,Z)
 = intersection_(j in J_E) union_(z in Z) E_z(j).        (25i)
```

Closures may be used when evaluating this displayed Lebesgue mass, but the
pointwise predicate always retains the open `D_z`: every drift-tooth seam is
uncovered by that drift and therefore remains in `S_(E,Z)` and in its
projection.  Any isolated body-safe grid point omitted by the positive-cell
ledger projects to zero, which belongs to every aligned `D_a`; hence the
pointwise equivalence `(25g)` is unchanged.

Thus projection changes the drift problem into a finite
cell-by-drift intersection-of-unions with rational interval literals.  The
choice of a drift owner in each cell is precisely the address sidecar that
singleton excess and the Gram forget.

Unlike the safe-surplus bound, the projected statement still has content at
`k=1`.  Here `u_1=6/7`, so the same proof gives

```text
max(E union Z)>L/13,
P_(E,Z) subset D_a
```

for the one aligned multiplier `a`.  This locates the scale and forces a
single-comb shape on the projected residual, but it does not bound the first
of the six drifts.  It is therefore a genuine constraint on the remaining
infinite sector, not a finite reduction of it.

There is an exact all-scale functional form behind the discrepancy tail.
Write the carrier components as

```text
C_E=union_s [a_s/L,b_s/L]
```

with integral endpoints, and put

```text
P(x)=integral_0^x 1_(D_1)(u)du,       delta_E(z)=c(z)-h/7.
```

The integrand in `P` is one-periodic with mean `1/7`, so
`P(x+m)=P(x)+m/7` for every integer `m`.  Direct substitution into

```text
z c(z)=sum_s[P(zb_s/L)-P(za_s/L)]
```

therefore gives the exact recurrence

```text
(z+L)delta_E(z+L)=z delta_E(z).                         (25l)
```

Equivalently, there is a residue amplitude

```text
A_E(b)=z delta_E(z),              z=b mod L,
delta_E(La+b)=A_E(b)/(La+b).                            (25m)
```

Here `A_E(0)=0`.  Extending the even coverage function `c(z)=c(-z)` to
negative labels also gives

```text
A_E(-b)=-A_E(b).                                        (25n)
```

The amplitude itself is a boundary object.  Center the primitive by

```text
g(x)=P(x)-x/7.
```

Then `g` is one-periodic, odd, and piecewise linear, and cancellation of
the linear part gives the finite endpoint formula

```text
A_E(b)=sum_s[g(b*b_s/L)-g(b*a_s/L)].                     (25o)
```

Thus the apparently analytic drift tail is the pairing of a fixed
toothpick wave with the signed endpoint boundary of `C_E` on `Z/LZ`.
Passing from `b` to its denominator remembers the cyclic quotient on which
this pairing lives but forgets the unit direction that evaluates it.  On
one period `g` ranges exactly from `-3/49` to `3/49`; hence each carrier
component contributes at most `6/49` to `A_E(b)`.  This recovers the
component constant in `(25d)` directly from the toothpick amplitude and
identifies the sharp oscillation scale behind that bound.

Thus every drift residue is an exact hyperbolic ray, not merely an
asymptotic `O(1/z)` tail.  Opposite residue directions have opposite excess;
the self-opposite denominator-two direction has zero excess.  If

```text
d=L/gcd(L,b),             b=(L/d)u,
```

then `u` is a unit modulo `d`, and moving along the ray preserves both `d`
and `u`.  For any fixed denominator multiplicity, its all-label excess
maximum is consequently obtained by merging finitely many decreasing
nonnegative rays.  In fact this maximum is attained: every negative unit
direction has a positive opposite, every zero opposite is already an
attained zero ray, and the denominator-two ray is identically zero.  The
reversal-paired unit directions form the oriented sidecar that a bare
denominator multiset forgets.  This exact ray law will replace finite-horizon
padding in the quotient-first addenda below; it also explains the equal
positive/negative ray counts in their audits.

The exact suffix verifier composes `(25c)`--`(25f)` without treating a
search horizon as exhaustive.  For every root and proposed `z_1`, it
integrates all allowed suffix labels through `H=7,000`, retains the largest
`d-1` distinct exact excesses, and pads omitted labels by
`6r_E/[49(H+1)]`.  When `(25f)` forces a later high label, the corresponding
tail starts at the larger of `H+1` and the strict projected-wall floor.  The
result is the necessary-filter census

```text
k                              2          3        4      5     6
suffix-only max z_1         2,340        432      260    130    44
projected-wall max z_1      2,142        380      182     66  EMPTY
projected surviving rows 2,239,853    376,020   87,975  4,702     0. (25j)
```

A row here is a pair `(E,z_1)` surviving rigorous upper envelopes, not a
realized cover.  The empty `k=6` column independently reproduces the known
one-drift closure; the other columns are finite banks, not emptiness
censuses.

The unique projected `k=3` row at the maximum `z_1=380` can now be removed
exactly.  It has body

```text
E=(1,4,8,10,12,14)
```

and exact suffix dynamic programming through `H=7,000`, with the inherited
omitted-label tail, leaves only

```text
(z_1,z_2,z_3)=(380,410,492),
z_4 in {1164,1220,1358,1500,1836}.
```

All five packets pass scalar capacity, but each has a body-support fibre
larger than the sum of the four sharp divisor-fibre caps: four have margin
`3` at `q=420`, and the `z_4=1358` packet has margin `1` at `q=140`.
Thus, before testing the next integer, the projected `k=3` first-drift cap
was `379`.  This is one closed frontier row, not a census of the remaining
`376,019` necessary rows.  Ordinary and optimized replays are byte-identical; script/output
SHA-256 are

```text
64f98439f677668c82045e7f9107cbfdff467afd8f16975c7e37d8ae5c5c9f26
a1c77b24488240f1ee0295e427ee4583b7d8215caf6615f424bf325350fb56b6.
```

The next integer row is empty already at the scalar suffix level.  At
`z_1=379`, the rigorous all-tail bound leaves `2,579` of the `3,003` body
rows for exact evaluation.  Exact singleton integration through `H=7,000`,
with three distinct later drifts, the projected high-label slot, and
THM-1094 padding for every omitted label, leaves zero rows.  The closest
strict rejection is again

```text
E=(1,4,8,10,12,14),
upper-lower=-4741191283/1316479619000.
```

The live positive control reconstructs the preceding `z_1=380` body and its
positive scalar margin `437649/1736780500`.  An independent rational-interval
carrier and guarded vector primitive reproduce the complete `2,579`-row
semantic digest

```text
48ab29334a93fd0087d9645513be14f884a30bd014c2f05329c1f7d0c295d4ee.
```

Therefore the current projected `k=3` first-drift cap is `378`.  The source
and output hashes for the integer-ruler referee and independent audit are

```text
88c563a247d59b2d9feb552935d91a2bbc5018beeed56df74c84a37a1174894b
12c1d60a6f1caf7f3a36a9bc890c388b4e44833a474e233038b5f79599715ae3
bd22ce0f86d9f5e359c2a940e0f8133849616e9aa9fc67eb823632ec9371f16d
ea7f6f2c9b189ffa4940fc25c58c74b13af905aed0fc7a6dc02266869775de77.
```

This removes no additional member of the old `376,020` necessary-row ledger:
there was no `z_1=379` row in it.  After the separate `z_1=380` packet
closure, `376,019` necessary rows remain.

The ray law also makes an all-label quotient-first closure possible far
inside that bank.  Fix

```text
E=(1,4,8,10,12,14),        z_1=250,        L=11760.
```

The first denominator is `L/gcd(L,250)=1176`.  For each multiset containing
this denominator and three arbitrary nontrivial divisors of `L`, merge the
first three eligible points on every nonnegative unit ray.  This gives the
exact attained scalar maximum over all later labels, while relaxing their
order and the projected high wall.  Of the `35,990` denominator multisets,
`1,965` survive this exact ray maximum.

The remaining quotient test strengthens THM-2928's common status table
`(37tg+)` by retaining forced pair overlap inside each status.  At an outer
divisor `q`, put `M=D/q`.  If two inner needles have denominators `d,f` and
lengths `e,z`, write

```text
g=gcd(d,f),       e=Ag+r,       z=Bg+s,       0<=r,s<g.
```

Distribution among the `g` common CRT classes and the Fréchet intersection
floor give

```text
I_(d,e;f,z)>=M/lcm(d,f) *
 (gAB+As+Br+max(0,r+s-g)).                              (25p)
```

For any number of needles, Hunter's tree inequality therefore gives the
status-wise union upper bound

```text
|union_i N_i|
 <=sum_i |N_i|-max_(spanning trees T)sum_(ij in T) I_ij.  (25q)
```

One nonnegative `2^p`-cell status table must have the exact activity
marginals and must dominate every tail of the target-load histogram.
Allowing a different table at each load would be unsound; infeasibility of
the common table is instead certified by an exact rational Farkas vector.
The crude all-divisor capacities remove `699` of the `1,965` scalar states,
and the common four-needle Hunter-status test removes all remaining
`1,266`.  Thus this entire `(E,z_1=250)` row is empty uniformly, with no
finite label horizon.  The source/output and semantic SHA-256 values are

```text
dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36
5abccb7ef700cec83b9989e8abcd83bc24f51c0a35f7f9054522da0dd62109fe
bcfa48e8b59080ced069a794d02cc04f62db8137f94b163b2fe4c98c3b3f77fa.
```

Ordinary and optimized transcripts are byte-identical.  Consequently this
separate uniform closure leaves `376,018` rows in the old necessary ledger;
it does not by itself change the current first-drift cap.

For `k=5`, there is a second, Gram-facing derivation.  Pointwise

```text
(m_A-1)_+ >=(2/5)binom(m_A,2),
```

while tensorization and THM-1234 give

```text
integral_(C_E)binom(m_A,2)
 =h sum_(i<j)rho(a_i,a_j)>=44h/273.
```

This reproduces `eta_5=88/1365`; the safe-surplus and pair-Gram views are
the same pressure in two coordinate systems.

The `k=5` finite bank is in fact empty.  Its `4,702` rows split exactly
according to the first-drift excess.

In the `4,084` high-excess rows

```text
delta(z_1)>=88h/1365,
```

apply THM-2893's six-tail first-apex gate to
`R=C_E minus D_(z_1)`.  If `h_R` and `r_R` are its mass and component
count, one of the five aligned labels or `z_2` is at most

```text
floor(36r_R/(7h_R)).
```

The aligned labels are at least `L`, while `(25f)` forces
`z_2>2275L/18627`.  These typed lower bounds immediately close `3,827`
rows.  The other `257` rows leave exactly `42,912` integral `z_2`
candidates.  On every candidate, the cell formula `(25i)` gives an exact
rational prefix lower bound

```text
mu(P_(E,{z_1,z_2}))>=887/1365=1-u_5,                    (25k)
```

contradicting the strict containment inequality `(25h)`.

On the subcritical side, `2,290` first rows have a nonempty finite analytic
interval for `z_2`; exactly `618` also survive the projected-suffix predicate.
For such a row put

```text
g=88h/1365-delta(z_1)>0.
```

The component discrepancy bound forces

```text
z_2<=floor((6r_E/49)/g).
```

Across the `2,290` analytic rows there are `7,218,110` row-labelled `z_2`
candidates.  Exact singleton integration together with the suffix predicate
leaves `194,073` admissible drift pairs, supported on `590` `(E,z_1)` rows
inside the `618`-row suffix bank.  Every one again satisfies `(25k)`.  The
smallest certified prefix margin over both banks is `1/378105`; equality at
the five-comb union cap would already be impossible because `P` is compact
and the aligned union is open.

A separate typed recursion independently checks the high-excess bank:
`39,913` of its `42,912` second-drift rows close at the first aligned gate;
the remaining `2,999` close before a one-label terminal, with no multiplier
above four.  Therefore

> No literal six-body/seven-tail cover has five aligned tails and two
> drifts.

THM-2928 now supplies a genuinely independent closure of the same face by
body/divisor projection and relaxed arithmetic-progression address masks.
The two proofs reverse the quantifier order.  THM-2928 retains the aligned
safe set and asks whether two drift masks cover every selected body address;
`(25g)` first removes the drifts and asks whether the aligned union covers
their existential projected residual.  They are dual cell-address
projections, and both succeed precisely because they retain the address
coordinate erased by singleton and Gram statistics.

More explicitly, let `I(j,u)=1` when the two drifts cover phase `u` in body
cell `j`.  THM-2928 fixes `u` in the aligned safe set and obstructs a full
column `I(j,u)=1` for every `j`.  Here

```text
P={u:there exists j with I(j,u)=0}
```

is the set of non-full columns, and its measure is too large for the aligned
union.  The two proofs are arithmetic column obstruction versus measure of
the De Morgan-dual column set on the same body-by-phase incidence object.

The caps in `(25e)` give a uniform **finite reduction** without bounding
the aligned multiplier set.  Delete the bounded `D_(z_1)`.  Six tail labels
remain, with no alignment assumption needed.  At every proper node of their
literal residual tree, the body speeds and already chosen tails total at
most twelve speeds.  Settled `LRC(<=13)` gives a point at clearance at least
`1/13`; the strict margin `1/13-1/14` gives a positive interval in that
residual.  It therefore satisfies the positive-mass interval hypothesis of
the cap-free `p<=6` first-apex recursion THM-2893(7a)--(7b).  Recursing to
depth at most six produces a finite exact decision tree.  This is a
finiteness theorem, not an emptiness census.  At each node the component
estimate inherited from THM-1094 is non-strict with coefficient `6r/49`;
replacing that coefficient by the explicit larger rational `6r/49+1`
supplies the strict form required by THM-2893 without changing finiteness.
The `1/13-1/14` margin is what keeps every proper-node residual positive.

Together with THM-2928's fully aligned and one/two-drift closures, this
proves:

> Every six-body/seven-tail branch with at least five aligned tails is empty.
> Each branch with `k=2,3,4` aligned tails is uniformly reducible to a finite
> exact decision tree.
> Consequently any sector not yet known to admit such a finite reduction
> has at most one aligned tail, hence at least six drifts.

THM-2928's later divisor-status transport and local `98/99` needle terminal
run the `k=4` census and prove it empty.  Thus the current composition
strengthens the first two conclusions to: every `k>=4` branch is empty, and
only `k=2,3` remain finite-but-uncensused.

The next faithful object must therefore join

```text
node data:
  delta_i, residue b_i, normalized slope a_i+b_i/L;
edge data:
  restricted overlap p_ij, gcd/reduced ratio, tooth indices;
component data:
  ordered owner word, transition positions and widths;
quotient data:
  projected residual P_(E,Z), its component word and measure;
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

Thus vanishing overlap widths are not anonymous analytic errors: they carry
an integer endpoint relation.  Because every `D_w` is open, an exit/entry
coincidence is an uncovered seam, not a two-owner handoff.  Equation `(27)`
records collision address only; it is not a persistence or capacity law.
Any event mesh must use the reduced winding `w/gcd(w,L)`, not raw `w`.
Equations `(16)`, `(25)`, and `(27)` define the multi-slope Gram/address
transition ledger to be classified.

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
implementations.  The aligned suffix census and five-aligned closure have
independent exact reconstructions and hostile proof audits.  The top-seven
overlap graph is a scoped structural scout.

The consolidated five-aligned closure verifier has LF-normalized
source/output SHA-256
`76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662`
and
`9aecfd75893a537278dcc4e50af7bd45fa2b7925d017748781a18c7163bb716d`;
ordinary and optimized replays are byte-identical.

This theorem does not give a uniform lower bound for `Delta` or `kappa`,
turn the `803` nonpositive actual-top-seven tree margins into certificates,
handle arbitrary packets, run its `k=2,3,4` finite decision trees, classify
the zero/one-aligned multi-drift address hypergraph, close the
six-body/seven-tail rung, or prove LRC(14).  The independent THM-2928
divisor-status route now closes `k=4`; `k=2,3` remain open.
