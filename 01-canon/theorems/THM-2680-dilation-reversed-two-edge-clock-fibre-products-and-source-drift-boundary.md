---
id: THM-2680
title: "Dilation-reversed two-edge clock fibre products and source-drift boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the canonical THM-2616 guard-safe raw carrier, D(x)={13x} identifies a
  current rail-owner clock with the following shallow clock and satisfies
  j(Dx)=h(x).  This reverses the printed clock arrows and makes the first
  genuine two-edge object E_0 intersect D^{-1}E_1, with both source labels
  retained.  Exact refined-grid computation finds two positive source-pair
  families whose s_1=6 and s_1=7 zero columns are inherited one-edge
  exceptions.  A sharp all-source hostile has positive one-edge relations
  and formal product but zero D-fibre product.  At source pair (1,1), 146 of
  252 clock triples compose formally, while only 10 survive physically.  The
  formal fixed-offset zero is therefore not a physical
  chronology no-go, but a two-edge support fibre product is not an iterated
  transition: source transport, three-edge coherence, units, endpoints, and
  LRC(14) remain open.
source: root-2026-07-28-dilation-clock-fibre-product
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2670-sharp-graph-clock-incidence-atlas-and-physical-gluing-boundary
related:
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2644-odd-torsor-purity-return-gate-and-nonlinear-fixed-branch-decoder
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
script: 04-computation/lrc14_dilation_reversed_clock_fibre_product_probe.py
output: 05-knowledge/results/lrc14_dilation_reversed_clock_fibre_product_probe.out
script_sha256: 56347893caa028a8d9f1b72c6e886dd38c617019f53e746e966bee65cb91a2ad
output_sha256: 1e7fb8babc74b52a4c18e05ff3604fae3282641df6bb4392a14db5bde266701a
hash_basis: LF-normalized bytes
---

# THM-2680 -- dilation gives a physical two-edge fibre product

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2670 isolates the typing error in the old clock products: different
rail-owner charts are disjoint, so ordinary Boolean multiplication joins
labels but not witnesses.  There is nevertheless a canonical map between
scales.  Pulling the following event back by base-thirteen dilation produces
a genuine common-`x` two-edge object.  It is nonempty, sparse, and sharply
source-sensitive.

## 1. Dilation reverses the clock arrows

Let

```text
D(x)={13x}.                                                (1)
```

Write `O_b` for the THM-2616 rail-owner clock cell labelled `b` and `S_a`
for its shallow clock cell labelled `a`.  Exact half-open interval arithmetic
gives

```text
O_b = D^(-1)(S_b),                  b in F_7.              (2)
```

The seven owner cells are pairwise disjoint.  The companion checks all `42`
ordered disjoint pairs and all seven identities `(2)`, including the wrapped
clock-zero presentation.

The retained numerical states are consecutive base-thirteen digits:

```text
j(x)=floor(13^6 x) mod 13,
h(x)=floor(13 {13^6 x}).                                  (3)
```

Hence, away from the measure-zero endpoints fixed by the half-open convention,

```text
j(Dx)=h(x).                                                (4)
```

Thus a `D`-compatible candidate chronology traverses THM-2670's printed
`(owner,shallow)` arrow backwards.  If its chronological clock path is

```text
a -> b -> c,                                               (5)
```

the current stored atom has `(owner,shallow)=(b,a)` and the following atom
has `(owner,shallow)=(c,b)`.

## 2. The correctly typed fibre product

For a rail occurrence `alpha`, source `s`, binary labels
`epsilon,kappa`, and sharp-graph states `j,h`, let

```text
E^(s)_(b,a;alpha,j,h,epsilon,kappa) subset O_b             (6)
```

be the unintegrated guard-safe atom used in THM-2670.  It is the intersection
of the THM-2616 middle rail, present packet, sharp nonzero deep root,
half-tooth, predecessor digit, carry half, and delayed-word digit.  Every
factor is nonnegative.

The direct intersection of successive stored atoms lies in disjoint owner
cells and is empty.  The `D`-typed two-edge object is instead

```text
F=E^(s0)_(b,a;alpha,j,h,epsilon,kappa)
  intersect
  D^(-1) E^(s1)_(c,b;alpha',h,k,epsilon',kappa').          (7)
```

The intermediate equality in `(7)` is exactly `(4)`.  Both source labels are
retained independently; no source-transition law is assumed.

The exact computation represents `(7)` on the refined integer grid `13T`.
If

```text
y={13^6x},                                                 (8)
```

the two delayed words reduce to the exact interval intersection

```text
Q_A(y) intersect Q_B({13y}).                              (9)
```

The following atom's remaining step profile is pulled back branch by branch
under `D`; the current profile is rescaled to the same grid.  Multiplying the
two positive integer weights preserves support.  The resulting numerator is
therefore an exact nonemptiness certificate, but is **not** asserted to be a
canonical transition mass.  Eight independent one-edge controls verify that
refining both grids multiplies the inherited numerator by exactly `13`.

## 3. Two positive source-pair families

Use the chronological triple notation `(shallow,current owner,next owner)`.
For

```text
(a,b,c)=(3,1,0),                                           (10)
```

the number of positive rail-labelled atom pairs in `(7)` is

```text
20  for 132 of the 144 ordered source pairs,
 0  exactly for (s0,s1)=(s,6),  s=1,...,12.               (11)
```

An exact positive numerator certificate at `(s0,s1)=(1,1)` is

```text
126816337986097204341478787325120                         (12)
```

and current/following labels

```text
(alpha,j,h,epsilon,kappa)=(2,5,2,0,0),
(alpha',j',h',epsilon',kappa')=(0,2,6,1,1).               (13)
```

In particular `h=j'=2`, as required.

For

```text
(a,b,c)=(4,0,1),                                           (14)
```

the corresponding census is

```text
14  for 132 ordered source pairs,
 0  exactly for (s0,s1)=(s,7),  s=1,...,12.               (15)
```

The first `(1,1)` positive numerator certificate is

```text
90407221954729505049108909024000                          (16)
```

and labels

```text
(1,11,3,0,0),                 (3,3,5,1,0).                (17)
```

Again `h=j'=3`.  The exceptional columns in `(11)` and `(15)` are inherited
one-edge zeros: the following banks `(owner,shallow)=(0,1)` at `s1=6` and
`(1,0)` at `s1=7`, respectively, are already empty before the `D`-pullback.
Thus these columns are not a new source-drift obstruction.  Both source
labels remain independently retained, and no source-transport rule follows.

## 4. Sharp hostile and sparse clock atlas

The nearby triple

```text
(a,b,c)=(3,1,2)                                             (18)
```

has zero fibre-product atoms for all `144` ordered source pairs, but for an
inherited reason: its following one-edge cell `(owner,shallow)=(2,1)` is
universally empty.  This is a useful zero-edge control, not a hostile to
`D`-gluing itself.

The sharp hostile is instead

```text
(a,b,c)=(0,1,0).                                           (19)
```

Its `D`-fibre product is empty for all `144` source pairs.  At `(s0,s1)=(1,1)`,
however, both one-edge relations are positive and their formal Boolean
product has support `29`.  Thus the failure in `(19)` occurs only after the
actual witness pullback and cannot be explained by a missing edge or an
incompatible printed state label.

At the generic fixed source pair `(s0,s1)=(1,1)`, exhaust all

```text
7*6*6=252                                                  (20)
```

ordered clock triples with distinct adjacent clocks.  Exactly `146` have
both one-edge relations positive, and all `146` also have positive formal
Boolean product.  Only ten have a positive physical `D`-fibre product.  The
tuple records `(a,b,c,positive rail-labelled atom-pair occurrence count)`:

```text
(3,0,3,6),  (3,0,4,10), (3,0,5,10), (3,0,6,12),
(3,1,0,20), (4,0,1,14), (4,0,2,14), (4,0,3,14),
(4,0,4,12), (4,6,0,18).                                  (21)
```

Only shallow clocks `3,4` occur in `(21)`.  This concentration is a candidate
for an iterated nilpotence mechanism, but no three-edge claim is made here.

## 5. What this changes, and what remains

THM-2670's fixed-step sevenfold Boolean products vanish in the stored arrow
orientation.  Equations `(10)`--`(17)` prove that this cannot be promoted to
a blanket physical chronology no-go: after applying the correctly oriented
`D` handoff, genuine two-edge common-witness fibre products exist.  The
inherited source exceptions `(11)`, `(15)` and the sharp physical/formal gap
`(19)`--`(21)` simultaneously show why one positive overlap is far from a
transition.

The next typed object is

```text
E_0 intersect D^(-1)E_1 intersect D^(-2)E_2,              (22)
```

with a lawful rule for the source labels and retained component data.  A
cycle further requires coherent iteration, not independently chosen source
pairs.  No Bockstein unit, endpoint owner, Perron-sheet selector, global
positive transition, holonomy trivialization, scalar row exclusion, or
LRC(14) conclusion follows.

## 6. Reproduction

Run

```bash
python3 04-computation/lrc14_dilation_reversed_clock_fibre_product_probe.py
python3 -O 04-computation/lrc14_dilation_reversed_clock_fibre_product_probe.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_dilation_reversed_clock_fibre_product_probe.out.
```

The computation uses exact integer interval arithmetic only.  Its finite
universe, positive controls, all-source hostile, refined-grid formula checks,
and current hashes are recorded above.

QED.
