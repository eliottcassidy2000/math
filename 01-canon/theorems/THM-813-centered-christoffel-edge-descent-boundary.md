---
id: THM-813
title: Projected coloured-edge transport first fails from X6 to X7, while reflection-orbit lines are the universal action carrier
status: PROVED (general reflection-orbit action and descent criterion) + FINITE-EXACT (X6 to X7 failure, all 40 natural core lifts, and two-step Mobius audit)
source: codex-2026-07-15-S13
depends_on: [THM-778, THM-796, THM-812]
related: [THM-781, THM-801, HYP-6880]
verification:
  - 04-computation/continued_fraction_edge_descent_boundary_codex_S13.py
  - 05-knowledge/results/continued_fraction_edge_descent_boundary_codex_S13.out
  - 05-knowledge/results/continued_fraction_edge_descent_boundary_codex_S13.json
---

# THM-813 — the first centered-CF edge-descent boundary

Let

```text
E_n=X_n/<kappa_n>
```

be complement lines, let staircase reflection `sigma_n` act on `E_n`, and
let

```text
P_n(e)=(blue/black colour, unordered merged-node endpoint pair). (1)
```

Every complement/reflection-equivariant tiling embedding acts canonically on

```text
Q_n=E_n/<sigma_n>.                                      (2)
```

It need not act on the coarser projected edge cells `P_n(E_n)`.  The centered-
Christoffel embedding `X_5->X_6` in THM-812 descends to `P_5`, but the immediate
next embedding `X_6->X_7` fails on 51 of 187 source cells.  This is the first
exact boundary between an edge identity codec and an action carrier.

## 1. General reflection-orbit action lemma

Suppose `Phi:X_m->X_n` satisfies

```text
Phi kappa_m=kappa_n Phi,       Phi sigma_m=sigma_n Phi.  (3)
```

The first identity induces `bar(Phi):E_m->E_n`.  Since `P_n sigma_n=P_n`,

```text
[e]_sigma -> P_n(bar(Phi)(e))                           (4)
```

is a well-defined function `Q_m->P_n(E_n)`.  It descends further through
`P_m` if and only if (4) is constant on all `Q_m` points contained in each
`P_m` fibre.

This proves the lemma.  Notice what it does not say: the map (4) need not be
injective, and it need not assign a target to a bare merged node.

At `n=5`, every `P_5` fibre is exactly one reflection orbit: there are eight
blue singleton cells and twelve black doubletons.  Hence `P_5=Q_5`, explaining
structurally why THM-812's first action descended.  At `n=6`, this coincidence
ends:

```text
|E_6|=512,             |P_6(E_6)|=187,             |Q_6|=272. (5)
```

## 2. The next centered coordinate copy

For `(q,p,s)=(4,3,0)`, THM-778 gives

```text
F=(1,2,4,5),          d=(1,2,1),          s'=1,     (6)
```

with no midpoint tie.  Replicate the `X_6` high and low legs by this centered
word and its reflection.  Recursively copy the three-bit `X_4` core into the
six-bit target `X_5` core.  In explorer order the target-to-source coordinate
surjection is

```text
rho=(0,1,2,2,3,4,5,6,6,7,8,5,7,8,9).                  (7)
```

It intertwines staircase reflection and therefore defines an injective
`Phi:X_6->X_7` satisfying (3).  Exhaustion of all 1,024 tilings gives zero
injectivity, complement, reflection, or colour failures.

The 187 source `P_6` fibre sizes are

```text
size                 1    2    4    6    8   12   14
number of cells      24  115   31    8    6    2    1. (8)
```

Their numbers of distinct target `P_7` values have histogram

```text
target support        1    2    3    4    6    7
number of cells      136   34    8    6    2    1.      (9)
```

Thus 51 cells fail descent: four blue and 47 black.

## 3. First witness and robustness

The first failure is already blue:

```text
source P cell       B(15,16)
literal lines       132, 370
target P cells      B(143,203), B(27,65)
target masks        4620, 11746.                       (10)
```

The one-tiles of line 132 are `(4,1),(6,3)`.  Those of line 370 are

```text
(5,1),(6,2),(5,2),(4,2),(5,3).
```

Both lines are individually staircase-fixed, so (10) is not the old black
reflection-pair ambiguity.

The core choice is not responsible.  Keep the centered leg map fixed and
enumerate every core-local reflection-intertwining coordinate copy: each
target two-cycle may map in either orientation to the source core two-cycle
or constantly to either fixed bit, and each target fixed bit may map to either
source fixed bit.  Of 64 candidates, exactly 40 are surjective.  All 40 split
the same `B(15,16)` cell; eight lifts fail on 51 `P_6` cells and 32 fail on 52.

## 4. Why `Q` is the correct carrier

Exactly 52 `P_6` cells contain more than one `Q_6` point.  Under (7), 51 of
them fail descent.  The sole accidental success is `K(12,29)`, whose four
lines form two reflection orbits and happen to share target `K(21,84)`.  Every
single-`Q` cell descends.

As the general lemma predicts, all 272 `Q_6` cells act with zero failures.
Their images occupy 268 target `P_7` cells, with four double collisions:

```text
target P       source Q representatives
K(21,35)       13,77
K(21,84)       24,88
K(189,245)     55,325
K(242,260)     198,454.                                 (11)
```

Therefore `Q` is a safe action carrier but is not asserted to be an injective
address.  Equation (4), not accidental injectivity at one size, is the
recursive invariant.

## 5. Composition and the Mobius sector

Composing (7) with THM-812's map gives

```text
rho_(5,7)=(0,1,2,2,2,3,0,4,4,5,4,0,5,4,5).            (12)
```

All 20 `P_5` cells still map to 20 distinct `P_7` cells because `P_5=Q_5`.
Bare nodes remain nonfunctorial: ten source nodes meet 36 target nodes, with
support histogram

```text
{1:3, 3:2, 4:1, 5:2, 7:2}.                              (13)
```

For any coordinate surjection, THM-812's coefficient pushforward composes:

```text
(rho_1)_* (rho_2)_* mu=(rho_1 rho_2)_* mu.               (14)
```

On the original five extension variables, the pulled-back `n=7` node
indicators have truncation cell counts

```text
core 0:  2,5,10,16,18,18; degree-3 collision {73,143,187},
core 1:  2,5,11,17,18,18; degree-3 collision {35,145}.   (15)
```

Degree at most four is again node-identifying.  This two-step stability is
evidence for a quartic node-indicator sector, not a universal closure theorem.
For arbitrary functions, target total degree is not closed: a high-degree
target monomial can collapse to a low-degree source monomial when several
coordinates share one `rho` image.  The exact closed packet for source degree
at most `D` is the `rho`-saturated target family

```text
{B: |rho(B)|<=D}.                                       (16)
```

## 6. Preservation boundary

This result challenges the assumption that colour plus endpoint
isomorphism-class nodes is an intrinsic edge.  It is a static projection of
literal line orbits, and at `n=6` one cell can contain several independently
transported reflection obligations.  The recursive vertices are therefore
`Q`-orbits, core/leg variables, and coefficient packets—not runners and not
bare tournament nodes.

The map preserves complement, staircase reflection, colour, and the target
edge projection of each `Q`-orbit.  It destroys owner clocks, metric walls,
sheet labels, and the LRC loneliness predicate.  Coupling to LRC still
requires THM-808's owner/root stalk and THM-778's phase/tie sidecars.
