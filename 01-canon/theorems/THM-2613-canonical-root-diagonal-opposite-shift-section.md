---
id: THM-2613
title: "Canonical root-diagonal opposite-shift section"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  THM-2604 all-root role k and its
  linked 13-divisible blocker 13b, the explicit physical inverse-root
  section q(r)=-kr mod 13 has a positive target-danger/blocker-safe slice
  on every root r.  Equivalently, the canonical invariant phase v=kr+q=0
  is absent from THM-2605's complete null-cell graph.  The section is exactly
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
tooth selects phase zero.  THM-2605's complete bad-graph classification proves
that phase zero is always positive.

This closes the choice problem for the **physical root/local-shift gate**.
It does not transport the old relation action through positive time, nor does
it type the local shift as the next THM-2585 target state.

## 1. The canonical section

Use THM-2605's notation

```text
d_L(x)=1_(||x||<L/14),            u_1=1-d_1,

iota_r(z)=(z+r)/13,

a=13b,                            13 does not divide k.   (1)
```

The legal all-root roles are `L=1,k>=12` (ordinary) and `L=2,k>=10`
(guard), exactly as in THM-2604/2605.

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

## 2. Phase zero is never bad

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
k=23: v=1,2,                                      (6)
```

with no bad pair at every other legal speed or blocker depth.  In particular

```text
0 notin projection_v B_(L,k,b)                            (7)
```

for every legal role.  By the definition of the bad graph, (7) is exactly

```text
mu{z:P_(r,q_*(r))(z)=1}>0             for every r.        (8)
```

This holds for every legal all-root role in (1), with no upper bound on `k`
or `b`.  The finite list (6) is complete because
THM-2605 first reduces the infinite speed/blocker family analytically by its
opposing full-tooth and full-safe-gap inequalities.

## 3. Equivariance and every finite root word

For every `Delta in F_13`,

```text
q_*(r+Delta)=q_*(r)-kDelta,                               (9)
```

so the graph of `q_*` is invariant under THM-2605's physical connection

```text
(r,q)->(r+Delta,q-kDelta).                                (10)
```

Choose one positive chamber from (8) in each root and a common-depth
base-13 cylinder inside every chamber.  THM-2605's digit-block specification
then gives, for every finite root word `(r_0,...,r_m)`, one common ancestry of
positive measure carrying

```text
(r_j,q_j)=(r_j,-kr_j)                    at every stop.   (11)
```

For a prescribed target correction `alpha!=0`, choose

```text
r_j=r_0+j k^(-1)alpha.                                   (12)
```

Equations (4) and (12) give

```text
q_j=q_0-j alpha,                  q_7=q_0-7alpha.         (13)
```

Thus the canonical section passes THM-2602's positive twisted-return support
test without choosing a phase or initial shift.

## 4. Which torsor is fixed, and which one is not

Normalize the local shift by

```text
qbar=-k^(-1)q.                                            (14)
```

On the graph (4), `qbar=r`.  Hence the physical target tooth canonically
chooses one fixed-generator equivariant identification between

```text
physical inverse-root digit r
     and
local THM-2599 paired-shift state qbar.                   (15)
```

This is genuine progress beyond an abstract thirteen-choice torsor in the
fixed physical gauge.  It is not invariant under an independent rebasing of
the root and shift torsors, and it supplies no THM-2542 chart gauge.  The
following arrows remain unproved:

```text
local paired-shift q       -/-> THM-2585 next section q',

physical root r            -/-> adjacent THM-2542 chart root,

future local character     -/-> old THM-2334 relation residue.   (16)
```

THM-2608 requires the first arrow together with the adjacent clock; THM-2610
shows positive chronology erases the old `C_13` deck before the third arrow
can be read off.  The THM-2611 proved candidate names a principal ancestry
bibundle as the proposed repair, but is not used as a proved dependency here.
The section (4) supplies none of those missing types merely by sharing the
alphabet `F_13`.

An independent audit checked the legal-role quantifiers, the phase-zero table,
the signs in (9)--(14), and the precise fixed-gauge meaning of "canonical".
No row is excluded, the ledger remains `165`, and LRC(14) remains open.  QED.
