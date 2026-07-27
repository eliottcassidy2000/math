---
id: THM-2613
title: "Canonical root-diagonal opposite-shift section"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  positive k,a with 13 not dividing k and 13 dividing a, for either an
  ordinary or guard target, the explicit physical inverse-root section
  q(r)=-kr mod 13 has a target-danger/blocker-safe slice of mass at least
  1/(182ak) on every root r.  In the THM-2604 range this is equivalently the
  statement that the canonical invariant phase v=kr+q=0 is absent from
  THM-2605's complete null-cell graph.  The section is exactly
  equivariant under (r,q)->(r+Delta,q-kDelta), realizes every finite root
  word on one base-13 ancestry, and with Delta=k^(-1)alpha gives the target
  correction q_j=q_0-j alpha.  This canonically chooses one physical
  root-to-local-shift identification, but does not identify that shift with
  THM-2585's next-target section or THM-2334's old relation-residue action,
  provide an adjacent chart clock or semantic endpoint, exclude a row, or
  prove LRC(14).
source: root-long-frontiers-2026-07-28-zero-phase-corollary
depends_on:
  - THM-2605-inverse-root-dipole-connection-and-mixed-square-invoice
related:
  - THM-2599-rootwise-opposite-shift-paired-slice-law
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
---

# THM-2613 -- the unshifted target tooth canonically selects the root section

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Two bare `C_13` torsors with fixed generators admit thirteen equivariant
identifications: choose the image of one point and equivariance determines the
rest.  The physical paired gate has more structure than two bare torsors: it
contains an unshifted target tooth.  Relative to the fixed inverse-root origin,
the label `q=0`, the chosen role `k`, and THM-2599's shift orientation, that
tooth selects phase zero.  A direct arithmetic argument proves that this
phase is always positive, with an explicit uniform rational floor; THM-2605's
complete bad-graph classification gives an independent check in the LRC
all-root range.

This closes the choice problem for the **physical root/local-shift gate**.
It does not transport the old relation action through positive time, nor does
it type the local shift as the next THM-2585 target state.

## 1. The canonical section

Put

```text
d_L(x)=1_(||x||<L/14),            u_1=1-d_1,

iota_r(z)=(z+r)/13,

a=13b,                            13 does not divide k,   (1)
```

where `k,b` are positive integers and `L` is either one or two.  The LRC
application uses the narrower THM-2604 all-root ranges `L=1,k>=12`
(ordinary) and `L=2,k>=10` (guard), but the local theorem below has no such
lower threshold.

The physical paired gate is

```text
P_(r,q)(z)
 =d_L((kz+kr+q)/13)u_1(bz-q/13).                         (2)
```

Its invariant phase is

```text
v=kr+q.                                                   (3)
```

Set

```text
q_*(r)=-kr mod 13.                                        (4)
```

Then `v=0` and the target factor in (2) is the same unshifted tooth on every
root:

```text
P_(r,q_*(r))(z)
 =d_L(kz/13)u_1(bz+kr/13),                               (5)
```

where changing the representative of `kr mod 13` changes the second argument
only by an integer.

## 2. Direct all-range positivity and its floor

Write `y=(r+z)/13`, with `0<z<1`.  On the canonical section (4), the two
factors are

```text
d_L(kz/13),                    u_1(bz+kr/13).             (6)
```

Suppose first that `r!=0`.  Since `k` is a unit modulo thirteen, the blocker
phase at `z=0` is a nonzero multiple of `1/13`.  Its distance from the open
danger arc is at least

```text
1/13-1/14=1/182.                                        (7)
```

Consequently (6) is one on an initial `z`-interval of length at least

```text
min(13L/(14k),1/(182b)).                                (8)
```

After the Jacobian `dy=dz/13`, both entries in (8) are at least
`1/(182ak)`.  This proves the claimed floor for every nonzero root.

It remains to treat `r=0`.  If `k<La`, choose `z` strictly between

```text
1/(14b)=13/(14a)             and             13L/(14k). (9)
```

This is a positive target-danger/blocker-safe collar.  Equality `k=La` is
impossible because `13|La` but `13` does not divide `k`; its `y`-mass is at
least `1/(14ak)`, stronger than the advertised floor.

Now suppose `k>La`, and put `alpha=a/k`.  There is an integer

```text
1<=j<k/13                 with                 ||j alpha||>=1/14. (10)
```

Indeed, if `alpha` lies in `[1/14,13/14]`, take `j=1`.  If
`alpha<1/14`, take `j=ceil(k/(14a))`; then
`j alpha in [1/14,1/7)`.  If `alpha>13/14`, put `d=k-a` and take
`j=ceil(k/(14d))`.  Here `a>13d`, so `j<=floor((k-1)/13)`, while
`||j alpha||=jd/k` lies in `[1/14,1/7)`.  These choices prove (10).

At

```text
z_j=13j/k,                                               (11)
```

the target phase is the integer `j`, while the blocker phase is `aj/k` and
is safe by (10).  If it is strictly safe, its distance from the nearest
`1/14` seam is a nonzero integer multiple of `1/(14k)`; if it is on a seam,
use the safe one-sided interval.  Dividing this phase length by `b`, and
intersecting with the target tooth, gives a `z`-interval of length at least
`1/(14bk)`.  The target radius is `13L/(14k)`, and each root-end distance is
at least `1/k`; both exceed the required `1/(14bk)`.  Thus the inequalities
`1<=j<k/13` keep the required one-sided interval inside `(0,1)`.  After
division by thirteen its `y`-mass is at least
`1/(14ak)`.  Thus, in all cases,

```text
mu(I_r intersection {P_(r,q_*(r))=1}) >= 1/(182ak)>0.   (12)
```

This direct argument is independent of any finite speed census.

## 3. Agreement with the complete bad graph

THM-2605 proves that the guard bad graph is empty.  In the ordinary case its
complete list is

```text
k=15: v=5,6;
k=16: v=2,3,7,8;
k=17: v=1,4,5,8;
k=18: v=2,3,5,6;
k=19: v=1,3,4,6;
k=20: v=1,2,4,5;
k=21: v=2,3;
k=22: v=1,3;
k=23: v=1,2,                                     (13)
```

with no bad pair at every other legal speed or blocker depth.  In particular

```text
0 notin projection_v B_(L,k,b)                           (14)
```

for every legal role.  By the definition of the bad graph, (14) is exactly

```text
mu{z:P_(r,q_*(r))(z)=1}>0             for every r.       (15)
```

This independently recovers the positivity part of (12) in every legal
all-root role, with no upper bound on `k` or `b`.  The finite list (13) is
complete because
THM-2605 first reduces the infinite speed/blocker family analytically by its
opposing full-tooth and full-safe-gap inequalities.

As in THM-2605, this is cellwise positive measure, not pointwise repair of
every `z`.  A positive-measure set of base points may fail for all thirteen
shifts while the canonical section still contains a positive chamber in every
root.

## 4. Equivariance and every finite root word

For every `Delta in F_13`,

```text
q_*(r+Delta)=q_*(r)-kDelta,                              (16)
```

so the graph of `q_*` is invariant under THM-2605's physical connection

```text
(r,q)->(r+Delta,q-kDelta).                              (17)
```

Choose one positive chamber from (12) in each root and a common-depth
base-13 cylinder inside every chamber.  THM-2605's digit-block specification
then gives, for every finite root word `(r_0,...,r_m)`, one common ancestry of
positive measure carrying

```text
(r_j,q_j)=(r_j,-kr_j)                    at every stop.  (18)
```

For a prescribed target correction `alpha!=0`, choose

```text
r_j=r_0+j k^(-1)alpha.                                  (19)
```

Equations (4) and (19) give

```text
q_j=q_0-j alpha,                  q_7=q_0-7alpha.        (20)
```

For the LRC all-root roles, the digit-block specification from THM-2605
therefore makes the canonical section pass THM-2602's positive twisted-return
support test without choosing a phase or initial shift.

## 5. Which torsor is fixed, and which one is not

Normalize the local shift by

```text
qbar=-k^(-1)q.                                           (21)
```

On the graph (4), `qbar=r`.  Hence the physical target tooth canonically
chooses one fixed-generator equivariant identification between

```text
physical inverse-root digit r
     and
local THM-2599 paired-shift state qbar.                   (22)
```

This is genuine progress beyond an abstract thirteen-choice torsor in the
fixed physical gauge.  It is not invariant under an independent rebasing of
the root and shift torsors, and it supplies no THM-2542 chart gauge.  The
following arrows remain unproved:

```text
local paired-shift q       -/-> THM-2585 next section q',

physical root r            -/-> adjacent THM-2542 chart root,

future local character     -/-> old THM-2334 relation residue.  (23)
```

THM-2608 requires the first arrow together with the adjacent clock; THM-2610
shows positive chronology erases the old `C_13` deck before the third arrow
can be read off.  THM-2611 proves the abstract principal-bibundle lift and
holonomy invoice, but no physical instance is used as a dependency here.
The section (4) supplies none of those missing types merely by sharing the
alphabet `F_13`.

Independent audits checked the all-range arithmetic lemma including strict
seams and the floor in (12), the legal-role phase-zero table, the signs in
(16)--(21), and the precise fixed-gauge meaning of "canonical".  No row is
excluded, the ledger remains `165`, and LRC(14) remains open.  QED.
