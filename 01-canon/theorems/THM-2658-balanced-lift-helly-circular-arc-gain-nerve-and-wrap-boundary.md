---
id: THM-2658
title: "Balanced lift-Helly circular-arc gain nerve and wrap boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  finitely many proper connected open circular arcs, positive common
  intersection is equivalent to a balanced complete section of the full
  integer lift-gain multigraph.  Balanced sections are canonically in
  bijection with common-intersection components; each component length is
  the minimum selected lifted pair-overlap, and total measure is the sum of
  these minima.  The result is invariant under switching reference lifts.
  After component refinement it gives the exact higher-overlap algorithm for
  finite-union physical charts.  Erasing winding, collapsing disconnected
  components, retaining only endpoint contacts, or using an incomplete graph
  is sharply insufficient.  No physical LRC chart simplex or row exclusion
  is asserted here.
source: wild-holotopy-mining/root-2026-07-28-balanced-lift-helly
depends_on: []
related:
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - HYP-3025-lrc14-closed-arc-cech-nerve-carrier
script: 04-computation/balanced_lift_helly_gain_nerve_thm2658.py
output: 05-knowledge/results/balanced_lift_helly_gain_nerve_thm2658.out
script_sha256: d075eda4493c51355e52b52693d90a5b0a46d911baab375e1bbe95aaa9cff4ca
output_sha256: 17635e7a3e6a72330b5c459b68507561664a2c81c60143e86318b55573af187d
hash_basis: LF-normalized bytes
---

# THM-2658 -- full integer gain balance is the missing circular Helly datum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Pairwise positive overlap fills only the one-skeleton of a chart nerve.  On
the line, pairwise overlap of intervals already forces total overlap.  On the
circle it does not: different pair edges can choose incompatible integer
lifts, and disconnected chart labels can choose incompatible components.
The exact repair is a complete balanced section of the **full integer** gain
multigraph.

## 1. Proper arcs and their full lift gains

Let

```text
T=R/Z,                       pi:R -> T.                    (1)
```

Let `A_1,...,A_N` be proper connected open arcs in `T`.  Choose reference
lifts

```text
I_i=(a_i,b_i) subset R,             0<b_i-a_i<1,
pi(I_i)=A_i.                                                (2)
```

For distinct labels define the finite integer gain set

```text
G_ij={n in Z: |I_i intersect (I_j+n)|>0}.                 (3)
```

Then `G_ji=-G_ij`.  For a label subset `S`, a **complete gain section** is a
choice

```text
n_ij in G_ij,             n_ji=-n_ij                     (4)
```

on every ordered pair in `S`.  It is **balanced** when every oriented
triangle satisfies

```text
n_ij+n_jk+n_ki=0.                                         (5)
```

The gains in (3) retain the integer winding sheet; reducing them modulo a
finite quotient is not the same datum.

## 2. Exact gain-Helly theorem

For every `S` with `|S|>=2`, the following are equivalent:

```text
mu(intersect_(i in S) A_i)>0;

the complete gain multigraph induced by S admits a balanced section.      (6)
```

More strongly, after the reference lifts (2) are fixed, balanced sections on
`S` are in canonical bijection with the connected components of
`intersect_(i in S)A_i`.

For a balanced section `sigma=(n_ij)`, let `C_sigma` be its component.  Its
length is exactly

```text
|C_sigma|=min_(i<j)|I_i intersect (I_j+n_ij)|,            (7)
```

where the minimum uses the **selected lifted pair components**, not the total
circular pair-overlap measures.  Consequently

```text
mu(intersect_(i in S)A_i)
 =sum_(sigma balanced) min_(i<j)|I_i intersect (I_j+n_ij)|.               (8)
```

Thus the positive open Cech nerve of proper connected arcs is exactly the
balance complex of their full integer gain multigraph.  Decorating each
balanced section by (7) retains both component multiplicity and mass.

## 3. Triangle balance produces interval potentials

After relabelling `S`, fix a base label `1 in S` and put

```text
h_1=0,                    h_i=n_(1i).                     (9)
```

The triangle `(1,i,j)` in (5) gives

```text
n_ij=h_j-h_i.                                             (10)
```

Hence the shifted intervals

```text
J_i=I_i+h_i                                               (11)
```

overlap pairwise: translating the selected overlap in (3) by `h_i` gives
`J_i intersect J_j!=empty` in positive length.  Put

```text
L=max_i left(J_i),             U=min_i right(J_i).        (12)
```

Choose `p` attaining `L` and `q` attaining `U`.  Pairwise overlap gives
`L<U`, so one-dimensional Helly yields

```text
K_sigma=intersect_i J_i=(L,U).                            (13)
```

If `p!=q`, then `J_p intersect J_q=(L,U)`.  If `p=q`, then `J_p=(L,U)` is
nested in every other `J_i`, so `J_p intersect J_j=(L,U)` for any `j!=p`.
Every selected pair contains `(L,U)`.  Therefore

```text
|K_sigma|=min_(i<j)|J_i intersect J_j|,                  (14)
```

which is (7).

Because each interval has length less than one, `pi` is injective on
`K_sigma`.  An endpoint of (13) is an endpoint of some proper `J_i`, whose
positive complementary gap prevents the intersection from continuing around
the circle.  Hence `pi(K_sigma)` is one full connected component.

## 4. Every component recovers one balanced section

Conversely, let `C` be a connected component of the common intersection.
There is a unique lift `K` of `C` contained in the base reference interval
`I_1`.  For each `i`, properness gives a unique integer `h_i` with

```text
K subset I_i+h_i.                                        (15)
```

Set `n_ij=h_j-h_i`.  Then (3)--(5) hold, and the construction in Section 3
recovers `K`.  Distinct components cannot yield the same section because
fixed potentials have the unique common interval (13).  This proves the
bijection and, by summing the disjoint component lengths, (8).

## 5. Switching invariance

Replace the reference lift `I_i` by `I_i+s_i`, `s_i in Z`.  The corresponding
gain labels change by

```text
n'_ij=n_ij+s_i-s_j.                                      (16)
```

Triangle sums and every selected overlap length are unchanged.  Thus the
individual integers are a switched presentation, while balance, components,
and the decorated nerve are intrinsic.

## 6. Finite-union physical charts

Suppose a physical chart is a finite disjoint union of proper open arc
components,

```text
U_i=disjoint_union_(alpha in E_i) A_(i,alpha).            (17)
```

Distributivity plus (6) gives the exact higher-overlap test:

```text
mu(intersect_i U_i)>0

iff there is one component alpha_i in E_i for every label i such that
    the complete integer gain multigraph of (A_(i,alpha_i))_i
    has a balanced section.                              (18)
```

For each successful transversal and section, its common mass is again (7).
Because the decompositions in (17) are disjoint, a physical point selects a
unique component at every label, so summing over transversals does not double
count.

Equation (18) is the lawful algorithm for the post-THM-2640 physical carry
atlas: expand every translated chart into individual open components, seek a
balanced transversal clique with full integer gains, then read its mass from
the smallest selected pair component.  A complete label-level one-skeleton is
not enough.

## 7. Sharp failure boundaries

Each hypothesis has a small hostile.

### 7.1 Erasing winding

Take reference lifts

```text
I_1=(-1/10,2/5),     I_2=(1/4,3/4),     I_3=(3/5,11/10). (19)
```

Their directed triangle gains and selected lengths are

```text
G_12={0},        G_23={0},        G_31={1},
lengths=(3/20,3/20,1/5).                                  (20)
```

Every circular pair overlaps positively, but the unique triangle has
holonomy `1`, so the triple intersection is empty.  Any quotient that erases
the integer winding misclassifies (20) as balanced.

### 7.2 Collapsing disconnected components

Let `X,Y,Z` be three disjoint open arcs and set

```text
U_1=X union Z,        U_2=X union Y,        U_3=Y union Z. (21)
```

Every pair overlaps positively—on `X`, `Y`, and `Z`, respectively—and a
coarse label-level gain-zero triangle looks balanced.  But the triple
intersection is empty.  The three edges chose incompatible components.  This
is exactly HYP-3025's reason to use individual physical arc components as
nerve vertices.

### 7.3 Endpoint-only contacts

The closed arcs `[0,1/3]` and `[1/3,2/3]` meet, but in zero mass.  Their open
positive-overlap gain set is empty.  Closed endpoint facets belong to
HYP-3025's boundary-cocircuit sidecar and cannot be inserted into (3).

### 7.4 Incomplete overlap graphs

For

```text
I_1=(0,2/5),       I_2=(3/10,7/10),       I_3=(3/5,1),   (22)
```

the path edges `12,23` both have gain zero.  The connected path is vacuously
balanced, but the triple intersection is empty because edge `13` is absent.
Completeness, equivalently all pair constraints, is what lets interval Helly
fire.

### 7.5 Whole-circle charts

A whole-circle label has no unique proper interval lift and admits ambiguous
sheet choices.  It must be split off separately; it is excluded by (2).

## 8. Relation to the live holotopy frontier

THM-2542 detects a nonzero graph-cycle chart cocycle.  The present theorem is
the complementary sufficiency statement: **zero triangle holonomy on a
complete component-refined lift graph produces an actual positive common
component**.  It does not turn a nonzero THM-2542 cocycle into zero.

THM-2657 proves that the desired carry/root translation quotient is a
nonsplit odometer extension.  That obstruction forbids a literal physical
`C_13` action but does not forbid positive overlaps between individual lifts.
Equations (7)--(18) specify exactly what additional same-base data would turn
those overlaps into a composable simplex.

No assertion is made here that the current thirteen physical LRC charts have
a balanced component transversal.  No owner/root endpoint, positive relation
current, row exclusion, or LRC(14) conclusion follows without that input.

## 9. Exact companion

Run

```bash
python 04-computation/balanced_lift_helly_gain_nerve_thm2658.py
python -O 04-computation/balanced_lift_helly_gain_nerve_thm2658.py
```

Both modes byte-match

```text
05-knowledge/results/balanced_lift_helly_gain_nerve_thm2658.out.
```

The dependency-free referee exhausts all `10,605` multisets of two through
four proper fifth-grid arcs, verifies `3,935` balanced sections and the exact
measure formula with checksum `990`, and performs `1,540` switching checks.
It also verifies the winding, disconnected-component, endpoint-only, and
incomplete-graph hostiles above.  An independent hostile audit rederived the
complete proof, all boundary cases, transcript, and hashes.  The LF-normalized
hashes are declared in the frontmatter.

QED.
