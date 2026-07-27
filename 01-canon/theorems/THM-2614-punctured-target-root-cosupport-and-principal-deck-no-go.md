---
id: THM-2614
title: "Punctured target-root cosupport factorization and principal-deck no-go"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On THM-2600's complete constant-six two-rail bank, every active fine
  (rail,target-section,future-clock) fibre has every nonzero deep-probe root
  and no zero root.  After retaining exactly the globally primitive unit
  target sections, the same-event cosupport in each of the 84 base cells is
  therefore Q_(s,ell) times F13-star, with sizes 156 on 74 cells, 144 on 8,
  and 132 on 2.  Every numerical root-minus-target difference occurs, but
  every diagonal C13 stabilizer is trivial and no total fixed-action or
  semilinear affine graph is contained because its r=0 point is missing.
  The bank contains 970 complete punctured fixed-action graphs: 13 choices
  on each full cell, one on each size-12 cell, and none on the two size-11
  cells.  Thus the retained support is a maximally dense punctured cospan,
  not the principal ancestry bibundle or selected transition demanded by
  THM-2608.  Aggregate unitness is never assigned sheetwise; no row is
  excluded and LRC(14) remains open.
source: wild-holotopy-2026-07-28-punctured-target-root-cosupport
depends_on:
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
related:
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
  - THM-2609-external-target-section-itinerary-saturation-and-root-state-no-go
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
script: 04-computation/lrc14_punctured_target_root_cosupport_thm2614.py
output: 05-knowledge/results/lrc14_punctured_target_root_cosupport_thm2614.out
script_sha256: fb2631e1b482bc1062337976cbbed6952b5f7e0c5fdaf8ce9be1d380edc54c56
output_sha256: f6034c31c6c2957e6a9868c9f89c1f67f639078e85e096b63252c31eefea9e39
hash_basis: LF-normalized bytes
---

# THM-2614 -- the available target/root cospan is a punctured product

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2609 proves that the external target sections have every additive
difference but do not become the constant-six physical future digit.  The
underlying fine bank retains one more same-event coordinate: THM-2600's deep
probe `r`.  Since THM-2608 needs a physical output-root-to-next-target map,
this is the nearest existing candidate for the missing connector.

The exact answer is a sharp no-selection boundary.  The `q` and `r` supports
do coexist physically, but their support is a direct product with the zero
root deleted.  It has every difference and contains many punctured affine
graphs, yet no complete `C_13` graph and no diagonal deck action.  More
incidence is not the missing datum; the missing datum is the thirteenth root
sheet together with a lawful choice of one graph.

## 1. The same-event fine relation

Use THM-2600's canonical typed row and its complete bank of `162` positive
constant-six middle rails.  For a rail

```text
e=(s,ell,t),                 t in {12,0},                 (1)
```

target section `q`, future owner clock `ell5`, and deep probe `r`, let

```text
J_(e,q;ell5,r) >= 0                                      (2)
```

be THM-2600 (8).  Its two target/root coordinates have different types:

```text
q: external absolute target-shift section, s5=-q;
r: same-event translated deep danger Delta_r=d(c3 x-r/13). (3)
```

Every factor in (2) is evaluated at one physical `x`.  Positivity is
therefore genuine fine cosupport, not a coupling of separate marginals.

The complete exact reconstruction gives the stronger fine law

```text
{r:J_(e,q;ell5,r)>0} is either empty or F_13^*.          (4)
```

Across all `162*13*7=14,742` fine fibres, the histogram is

```text
root support 0:       9,638;
root support 12:      5,104.                              (5)
```

Thus factorization occurs before the future-clock marginal and before the
two rails in one base cell are united.  It is not created by a late sum.
The missing root is literal: `r=0` is identically zero in the inherited bank
because the present deepest-safe factor is disjoint from `Delta_0`.

## 2. Unit-section cosupport is exactly `Q times F_13^*`

The unit predicate remains exactly THM-2600's aggregate predicate.  For one
rail `e` and target section `q`, form its globally primitive seven-clock
polynomial `Y_(e,q)(z)` only after the single global content

```text
g0=4,244,240                                               (6)
```

is removed.  Call `(e,q)` unit when `Y_(e,q)` is a unit of
`F_13[z]/(Phi_7)`.  This predicate belongs to the whole `(ell5,r)` aggregate;
it does **not** select or charge one positive entry in (2).

For each base cell `(s,ell)`, define

```text
Q_(s,ell)
 ={q: some rail e over (s,ell) has unit Y_(e,q)},          (7)

R_(s,ell)
 ={(q,r): some such unit rail and some ell5
           have J_(e,q;ell5,r)>0}.                        (8)
```

The fine law (4), checked with the exact unit flags, gives

```text
R_(s,ell)=Q_(s,ell) x F_13^*                              (9)
```

in every one of the `84` cells.  The exact size census is

| `|Q|` | `|R|=12|Q|` | number of cells |
|---:|---:|---:|
| 13 | 156 | 74 |
| 12 | 144 | 8 |
| 11 | 132 | 2 |

The ten deficient cells and their missing target sections are

| `(s,ell)` | missing `q` |
|:---:|:---|
| `(2,0)` | `{8,9}` |
| `(2,2)` | `{11}` |
| `(2,6)` | `{10}` |
| `(6,3)` | `{2}` |
| `(6,5)` | `{2}` |
| `(7,2)` | `{11}` |
| `(7,4)` | `{11}` |
| `(11,0)` | `{4,5}` |
| `(11,1)` | `{3}` |
| `(11,5)` | `{2}` |

The same factorization holds with `ell5` fixed.  Among the `84*7=588`
future-clock relations, the active-`q` counts are

```text
0:204,  9:32,  10:30,  11:226,  12:76,  13:20.          (10)
```

Every one of the `384` nonempty relations has all thirteen differences
`r-q`; the other `204` are empty.

## 3. Full difference support is not a connector

Because every `Q_(s,ell)` has at least eleven elements, (9) immediately
implies

```text
{r-q:(q,r) in R_(s,ell)}=F_13                            (11)
```

for all `84` cells.  Indeed, for a prescribed difference `d`, choose
`q in Q_(s,ell)` different from `-d` and set `r=q+d`, which is nonzero.

This numerical saturation is maximally weak as transition data.  The
literal diagonal `q=r` has

```text
12 points on 74 cells,
11 points on  8 cells,
10 points on  2 cells,                                  (12)
```

but (9) contains all cross-pairs as well.  It is a complete bipartite
cosupport, not the graph of the displayed diagonal or of any other map.

The failure is also equivariant.  Translate both coordinates by `a`:

```text
(q,r) -> (q+a,r+a).                                      (13)
```

For `a!=0`, the second projection `F_13^*` becomes
`F_13\{a}`, so (13) cannot preserve (9).  Hence every one of the `84`
diagonal translation stabilizers is trivial.

## 4. Punctured affine graphs and the exact completion defect

A fixed-action affine identification of two `C_13` torsors has graph

```text
Gamma_c={(c+r,r):r in F_13}.                             (14)
```

Its punctured graph is

```text
Gamma_c^*={(c+r,r):r in F_13^*}.                         (15)
```

Since the first projection of (15) is `F_13\{c}`, equation (9) gives the
exact criterion

```text
Gamma_c^* subset R_(s,ell)
iff F_13\{c} subset Q_(s,ell).                           (16)
```

Consequently the numbers of complete punctured gauges per cell are

```text
13 on the 74 full cells,
 1 on the  8 size-12 cells,
 0 on the  2 size-11 cells.                              (17)
```

There are `74*13+8=970` such graphs in total.  Across all `84*13=1,092`
fixed-action gauges, their numbers of retained punctured points are

```text
12:970,        11:100,        10:22.                     (18)
```

Restoring a **total** affine graph requires its missing `r=0` point and,
on deficient cells, any additional missing positive-root points.  The
minimum completion profiles `(missing points, minimizing gauges)` are

```text
(1,13) on 74 cells,
(1, 1) on  8 cells,
(2, 2) on  2 cells.                                     (19)
```

No total graph (14) is contained in the bank.  More generally every
semilinear affine bijection

```text
r -> kappa r+c,            kappa in F_13^*,              (20)
```

also contains one point over `r=0`; none of the `156` such graphs is
contained in any cell.  The companion checks all of them directly.

Equations (17)--(19) are selectors and completion invoices, not permission
to add the missing physical sheet.  In particular the thirteen choices on
a full cell reproduce the unresolved reference multiplicity rather than
choosing one reference canonically.

## 5. Holotopy interpretation and stopping boundary

THM-2608 asks for a same-carrier next target index `q'` attached to the
outgoing deep root and then for adjacent-clock composition.  The relation
(9) supplies neither:

```text
source:        one physical same-event q/r cosupport;
incidence:     every active q with every r!=0;
missing:       the r=0 sheet;
choice:        no selected q'=phi(r);
composition:   no outgoing-root/next-input gluing.        (21)
```

In path language, a genuinely composable alternating rail/connector loop
would make its rail increments and connector increments telescope to zero.
Thus a nonzero open-rail holonomy requires the connector to carry the
opposite class.  Equation (11) realizes every **numerical** class, but (9)
does not choose one class coherently along a path.  Full difference support
is therefore not transition curvature.

The separately routed principal-bibundle theorem THM-2611 gives the abstract
thirteen-state classification of this failure.  It is not used as a proved
dependency here: (13)--(20) directly prove the physical-bank obstruction.
Likewise THM-2609's chronology and THM-2610's later same-root character graft
do not restore the absent same-event `r=0` point or select one graph in (16).

## 6. Exact evidence and scope

Run

```text
python 04-computation/lrc14_punctured_target_root_cosupport_thm2614.py
python -O 04-computation/lrc14_punctured_target_root_cosupport_thm2614.py
```

The companion rebuilds all `162` rails, `61,248` positive fine entries,
`1,740` globally primitive unit slices, and all `14,742` fine root fibres.
It checks (4) before any marginal, all `588` fixed-future-clock relations,
all `84` products (9), the deficient atlas, all differences and diagonal
stabilizers, all `1,092` punctured fixed-action graphs, and all total
fixed-action and semilinear graphs.  It uses exact integer arithmetic,
exhaustive finite universes, and explicit optimized-mode guards.  Normal and
optimized runs byte-match the stored transcript after LF normalization.

The theorem is confined to support on THM-2600's one canonical typed row.
It does not compare coefficient magnitudes or phases, assign the aggregate
unit Bockstein to a fine sheet, create a principal ancestry fibre, identify
a semantic owner/repair endpoint, compose seven transitions, exclude a
scalar row, or prove LRC(14).  The ledger remains `165`.

QED (candidate; independent hostile audit pending).
