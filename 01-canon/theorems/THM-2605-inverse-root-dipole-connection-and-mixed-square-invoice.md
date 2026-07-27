---
id: THM-2605
title: "Inverse-root dipole connection and uniform fixed-head paired paths"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  THM-2604 all-root role k and its actual paired 13-divisible blocker 13b,
  the physical inverse-root gate has invariant phase v=kr+q under
  (r,q)->(r+Delta,q-kDelta).  The guard null-cell graph is empty.  The
  ordinary null-cell graph has exactly 28 edges, all with b=1 and
  15<=k<=23, and its phase projection has size at most four.  Hence at
  least nine phases make the paired target-danger/blocker-safe cell positive
  on all thirteen physical roots simultaneously.  One phase realizes every
  prescribed finite physical-root word on a common base-13 ancestry, and
  every prescribed initial root has at least nine good initial shifts.
  Choosing Delta=k^(-1)alpha gives q_j=q_0-j alpha and therefore supplies
  the positive physical support required by THM-2602's seven-edge twisted
  return test.  This does not identify the physical root/shift with a
  THM-2542 chart/clock or THM-2585 next-target state, construct a THM-2334
  relation current, preserve a future semantic endpoint, exclude a scalar
  row, or prove LRC(14).
source: common-endpoint-seam-2026-07-28-fixed-head-affine-paired-path
depends_on:
  - THM-2599-rootwise-opposite-shift-paired-slice-law
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
related:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
script: 04-computation/lrc14_fixed_head_affine_paired_path_referee.py
output: 05-knowledge/results/lrc14_fixed_head_affine_paired_path_referee.out
script_sha256: 29be125d8174804cf688e7fd7d6be8ebb70e8ddcf29ea4faf2e89a1ab5df2a02
output_sha256: 3ae8a95514d021a719c0f52c9c7728ef4604061974bd6a094dff4ce7449a7cf1
hash_basis: LF-normalized bytes
---

# THM-2605 -- one affine phase serves every physical inverse root

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2599 proves a positive oppositely shifted target/blocker slice in every
physical root, but its shift may be chosen independently in each root.
THM-2602 shows why unrelated vertex choices cannot supply an ordered
seven-edge transition.  The present theorem classifies the complete failure
graph of the paired physical gate and finds a single affine section which is
positive on all thirteen roots at once.

The mechanism is not a larger marginal census.  The coordinate

```text
v=kr+q                                                     (1)
```

is invariant under the physical inverse-root connection, and the widths of
one complete danger tooth and one complete safe gap force all possible null
cells into a 2,704-cell exact boundary problem.  Its 28 survivors occupy at
most four `v`-rows.  Every other row is a coherent all-root section.

## 1. Physical gate and connection

Put

```text
d_L(x)=1_(||x||<L/14),          u_1=1-d_1,

iota_r(z)=(z+r)/13,             r in F_13,

a=13b,                          b>=1.                     (2)
```

Here `L=1` for an ordinary role and `L=2` for the guard.  The THM-2604
all-root thresholds are

```text
ordinary: L=1, k>=12;

guard:    L=2, k>=10;           13 does not divide k.     (3)
```

For shift `q in F_13`, inverse-root substitution in the THM-2599 paired
gate gives

```text
P_(r,q)(z)
 =d_L((kz+kr+q)/13) u_1(bz-q/13),       0<z<1.            (4)
```

Define

```text
A_v={z:d_L((kz+v)/13)=1},

C_q={z:d_1(bz-q/13)=1}.                                  (5)
```

The cell `(r,q)` is null exactly when

```text
A_(kr+q) subset C_q                    up to endpoints.    (6)
```

Call `(v,q)` **bad** when `A_v subset C_q`.  The affine connection

```text
tau_Delta(r,q)=(r+Delta,q-kDelta)                         (7)
```

preserves the phase:

```text
k(r+Delta)+(q-kDelta)=kr+q=v.                             (8)
```

Thus the correct object is the bad graph in `(v,q)`, not thirteen unrelated
rootwise choices.

## 2. Exact interval geometry

The target-danger components are the intersections with `(0,1)` of

```text
((13m-v)/k-13L/(14k), (13m-v)/k+13L/(14k)),              (9)
```

and the blocker-danger components are the intersections of

```text
((m+q/13)/b-1/(14b), (m+q/13)/b+1/(14b)).               (10)
```

Their complete tooth/gap lengths are

```text
target danger:  13L/(7k),       target safe: 13(7-L)/(7k),

blocker danger: 1/(7b),         blocker safe: 6/(7b).    (11)
```

For ordinary roles every target-safe component, including a component
clipped by an endpoint, has length at most `78/(7k)`.

### 2.1 The guard graph is empty

For `L=2,k>=13`, periodicity gives

```text
mu(A_v)>=(13/k)(2/7)floor(k/13)>1/7=mu(C_q).              (12)
```

Indeed, writing `k=13n+s`, with `n>=1` and `1<=s<=12`, the strict inequality
is equivalent to `13n>s`.  At `k=10,11,12`, exact endpoint clipping gives
the minimum over `v` of the largest target component as

```text
3/35, 6/77, 13/84,                                       (13)
```

all larger than `1/14`.  This rules out `b>=2`; the `507` remaining `b=1`
cells are empty by exact rational containment.  Therefore

```text
B_(2,k,b)=empty                          for all legal k,b. (14)
```

### 2.2 Ordinary badness forces `b=1` and `k<=25`

For the two small legal speeds `k=12,14`, the corresponding component floors
are `13/168` and `13/98`, both larger than `1/14`; hence `b>=2` is impossible.
The `b=1` alignments are included in the exact boundary census below.

For `k>=15`, the interval traversed by `(kz+v)/13` has length
`k/13>1+1/7`, so `A_v` contains one complete ordinary danger tooth.  If
`A_v subset C_q`, its width comparison gives

```text
13/(7k)<=1/(7b),                    so k>=13b.             (15)
```

When `b>=2`, `(0,1)` contains a complete blocker-safe gap.  Complement
containment and (11) give

```text
6/(7b)<=78/(7k),                    so k<=13b.             (16)
```

Equations (15)--(16) force `k=13b`, excluded by (3).  Hence `b=1`.

For `q=0`, the blocker-safe gap is one interior interval of length `6/7`,
so the same comparison forces `k<=13`.  Thus `q!=0`.  Then the blocker-safe
set has two endpoint pieces of total length `6/7`; each lies in an endpoint
target-safe piece.  Consequently

```text
6/7<=2(78/(7k)),                     so k<=26.             (17)
```

The endpoint `26` is thirteen-divisible.  The entire unresolved ordinary
universe is therefore

```text
b=1,             k in {12,14,15,...,25}.                 (18)
```

Together with the `507` guard cells, it contains exactly `2,704` phase/shift
cells.

## 3. The complete bad graph

Two independent exact rational routes give the same table: direct interval
union containment, and midpoint integration over the complete wall atlas.
An arrow denotes `v->q`.

| `k` | complete bad-pair list |
|---:|:---|
| 15 | `5->7, 6->6` |
| 16 | `2->9, 3->8, 7->5, 8->4` |
| 17 | `1->9, 4->7, 5->6, 8->4` |
| 18 | `2->8, 3->7, 5->6, 6->5` |
| 19 | `1->8, 3->7, 4->6, 6->5` |
| 20 | `1->8, 2->7, 4->6, 5->5` |
| 21 | `2->7, 3->6` |
| 22 | `1->7, 3->6` |
| 23 | `1->7, 2->6` |

There are no other bad pairs.  In particular, the graph is empty for every
guard, for ordinary `k=12,14,24,25`, for every `b>=2`, and for every larger
legal speed.  It has exactly `28` edges and satisfies

```text
|projection_v B_(L,k,b)|<=4                              (19)
```

for every legal role.

## 4. One phase realizes every finite root word

Choose

```text
v in F_13 minus projection_v B_(L,k,b).                   (20)
```

There are at least nine choices by (19), and all thirteen choices work for
the guard.  Define the affine section

```text
q_v(r)=v-kr.                                              (21)
```

For every root `r`, `(v,q_v(r))` is not bad, so

```text
mu{z:P_(r,q_v(r))(z)=1}>0                                (22)
```

simultaneously on all thirteen roots.

Choose one positive open chamber for each root, map it through `iota_r`, and
choose a common-depth base-13 cylinder `D_r` inside it.  For
`T(x)=13x mod 1`, every finite root word `(r_0,...,r_m)` and every block
separation `N` at least that common depth `ell` obeys

```text
mu(intersection_(j=0)^m T^(-jN)D_(r_j))
 =13^(-(m+1)ell)>0.                                      (23)
```

Every stop carries the coherent shift `q_v(r_j)`.  Thus one phase embeds the
full one-sided thirteen-shift with the affine cocycle (21), including repeated
roots and arbitrary finite path length.

For a prescribed initial root, `q_0 -> v=kr_0+q_0` is a bijection.  Hence at
least nine initial shifts work for every prescribed finite path.  The lower
bound is sharp in the eight-stop diagnostic: at `(k,r_0,delta)=(16,0,8)` the
good initial shifts are

```text
{0,1,4,5,6,9,10,11,12}.                                  (24)
```

## 5. Exact twisted-return supplier and boundary

Write an affine root path as

```text
r_j=r_0-j delta.                                          (25)
```

Then (21) gives

```text
q_j=q_0+kj delta.                                         (26)
```

For any nonzero target twist `alpha`, choose

```text
delta=-k^(-1)alpha.                                       (27)
```

The eight stops satisfy

```text
q_j=q_0-j alpha,                 q_7=q_0-7alpha.          (28)
```

Equation (23) therefore supplies a positive ordered physical ancestry which
passes THM-2602's necessary seven-edge twisted-return support test on the
actual target-shift alphabet.

This is a transition supplier only at the physical root/shift level.  The
precise connection ledger is

```text
source:       physical inverse-root/shift pairs (r,q),
target:       invariant phase fibres kr+q=v,
map:          tau_Delta(r,q)=(r+Delta,q-kDelta),
preserved:    paired target-danger/blocker-safe positivity,
lost:         chart clock, chamber location, owner/repair endpoint, full X,
needed:       lawful chart/root-to-next-target identification and gluing.
                                                                  (29)
```

In particular, `r`, the coordinate translation `q`, and a THM-2334
left-minus-right relation residue are three different objects.  Neither
their common cardinality nor (28) identifies them.  THM-2608's clock-matched
collapse, THM-2609's external-section diagonal obstruction, and THM-2610's
old/future action-axis boundary remain in force.  No scalar row is removed;
the ledger stays `165` and LRC(14) remains open.

## 6. Exact companion and audit

Run

```text
python3 04-computation/lrc14_fixed_head_affine_paired_path_referee.py
python3 -O 04-computation/lrc14_fixed_head_affine_paired_path_referee.py
```

Both modes byte-match the stored output.  The dependency-free companion uses
`Fraction` throughout and checks:

- both exact interval routes on all `2,704` reduced cells;
- the complete 28-edge table and phase-projection bound;
- a hostile box of `1,774,500` cells with `k<=200,b<=30`;
- all affine-section and twisted-return signs;
- the exact eight-stop histogram
  `{9:86,10:268,11:584,12:372,13:94}`; and
- the common-depth cylinder invoice.

An independent audit replayed normal and optimized modes and hashes,
rederived the tooth/gap reduction including both `b=1` endpoint geometries,
and checked the connection sign and ancestry specification.  It found no
defect.  QED.
