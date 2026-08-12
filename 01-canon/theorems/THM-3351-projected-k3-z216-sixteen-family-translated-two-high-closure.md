---
id: THM-3351
title: "Projected-k3 z216 sixteen-family translated two-high closure"
status: >
  PROVED analytic mechanisms + FINITE-EXACT in the declared projected atlas.
  The first sixteen complete intrinsic-cost families after THM-3320, 236
  rows in all, have empty necessary projected screens.  A translation-uniform
  weighted-band lemma repairs the two-high terminal, with a separate exact
  denominator-two measure argument.  The projected ledger falls 373153 to
  372917, wall rows 349 to 113, and complete families 29 to 13; the cap remains
  216.  This does not prove physical entry, the rung, arbitrary k, or LRC(14).
source: root/lrc-math-2026-08-12
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure
related:
  - THM-3308-threshold-chain-modular-multicovers-and-three-layer-status-circuit
script: 04-computation/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.py
output: 05-knowledge/results/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.out
script_sha256: cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62
output_sha256: a88646fbd28d807a0cc9671c509c4424056a539b49d04a2076ba17de57ef5ee4
hash_basis: LF-normalized bytes
---

# THM-3351 -- projected-`k=3`, `z1=216`, sixteen-family translated two-high closure

**PROVED analytic mechanisms + FINITE-EXACT in the declared projected
atlas.**

## 1. Statement

Reconstruct the pinned necessary projected `k=3,z1=216` scalar-body atlas and
the exact queue state after THM-3320.  The atlas has `480` rows: `447` wall
rows and `33` order rows.  THM-3320 leaves `349` live wall rows in `29`
complete intrinsic-ruler families and projected ledger `373,153`.

The first sixteen complete families in the inherited intrinsic-cost order,
from `gcd72/L7056` through `gcd72/L35280`, contain `236` rows.  Every one has
an empty necessary projected screen.  Therefore, within this declared
projected universe,

```text
projected ledger       373153 -> 372917,
z1=216 wall rows          349 -> 113,
complete families          29 -> 13,
projected k=3 cap         216 unchanged.                 (1)
```

The cost boundary is strict: the last selected family has cost `42,124,320`;
the next, `gcd72/L720720` on rows `191,228,332`, has cost `50,450,400`.

This theorem is a necessary-screen closure.  It does not restore physical
entry, endpoint origin, owner, phase, current, arbitrary `k`, the rung, or
LRC(14).

## 2. Exact queue and common-status screen

The selected families are

```text
gcd72/L7056    19 rows       gcd72/L144144   1 row
gcd24/L64680    2 rows       gcd72/L5040    41 rows
gcd24/L240240   1 row        gcd72/L38808    6 rows
gcd72/L55440    9 rows       gcd24/L11760   40 rows
gcd72/L65520   10 rows       gcd72/L17640   33 rows
gcd72/L360360   2 rows       gcd72/L77616    8 rows
gcd36/L97020    7 rows       gcd24/L336336   2 rows
gcd72/L91728   11 rows       gcd72/L35280   44 rows.     (2)
```

The exact ray/common-status screen partitions `227,986` denominator states
as

```text
109,375  crude capacity exclusions,
118,093  exact common-status Farkas exclusions,
    518  residual passports.                             (3)
```

All `118,093` status exclusions are exact-rational checks on all sixteen
Boolean common-table columns: two direct zero-good threshold-row witnesses
and `118,091` inherited full-table witnesses.  The frozen prefix and screen
digests are

```text
7c665728bf21bf641f70d5e683fc615b1a24d3b3cdb0595462cdd17d0daa7721,
d8b0e7e18794bea30c038193db1949585e07e38f3011f311b71570e26f635963. (4)
```

## 3. The high gate and the one-high terminal

By THM-2941 `(25f)`, every projected completion has largest drift

```text
M>13L/132.                                               (5)
```

Put `H=floor(13L/132)+1`.  The exact atlas check gives `216<H` on every
selected wall row, so at least one of the three later drifts is high.  This
uses `(5)`, not an optional printed `HIGH-TAIL` token.  Of the `518` residual
denominator states, `461` pass scalar arithmetic with zero highs; they are
retained as hostile accounting but excluded by `(5)`.

For `6,126` enlarged `(state,high denominator,low-label pair)` cases, let `C`
be the complete cells fixed-safe for the first label and both low labels.
THM-2979's located-torsion certificate finds two cells whose difference has
effective order at most seven modulo the high denominator.  Every unit
numerator separates the two high phases by at least `1/7`; their strict-open
radius-`1/14` danger arcs are disjoint, so one cell survives.  All `6,126`
cases certify.

On thirty of the thirty-five residual bodies, the duplicate-permitting
two-or-more-high scalar bound has positive exclusion gap.  With the at-least-
one-high gate, those bodies therefore have exactly one high and close here.

## 4. A translation-uniform weighted-band lemma

Five `L35280` bodies, atlas rows `146,171,316,350,427`, have nonpositive
duplicate two-high gap.  Exact scalar enumeration leaves twenty-one enlarged
two-high cases and no three-high case.  A common-modulus certificate closes
three.  The remaining eighteen require the following sectionwise lemma.

Let `C` be the multiset of complete cells fixed-safe for the first drift and
the sole low drift.  For a denominator `d|L`, define

```text
B_C(d)=max_(u,J) #{c in C : u*c mod d lies in J},
m_C^tr(d)=|C|-B_C(d),                                   (6)
```

where `u` ranges over `(Z/dZ)^*` and `J` over every translated open cyclic
interval of length `d/7`.  Multiplicity in `C` is retained.

THM-2984 `(13a)--(13h)` supplies the exact translated-band capacity and
finite cyclic-window compression, while `(23b)--(23d)` supplies the
complete-cell/local-coordinate projection law.  For a concrete high drift
`z=(L/d)u+hL` and each fixed local coordinate `y`, its dangerous cells are
exactly the cells for which `u*c mod d` belongs to one such translated band
`J_y`.  Hence for arbitrary high denominators `d_1,...,d_r`, ordinary
inclusion--exclusion gives, pointwise in `y`,

```text
#{c in C safe from all r highs at y}
  >= sum_i m_C^tr(d_i)-(r-1)|C|.                        (7)
```

If the right side is positive, every local coordinate has a complete cell
safe from the first, low, and all high drifts.  Thus the projected safe set
is the whole circle, `P_(E,Z)=T`, contradicting THM-2941 `(25h)`, which
requires `mu(P_(E,Z))<36/91` for a completion.

The integer residues in an open cyclic interval of length `d/7` are
contained in a block of `ceil(d/7)` consecutive cyclic positions; the
intersection can be smaller when a strict endpoint passes through a lattice
point.  Conversely every full block of `ceil(d/7)` positions occurs for a
suitable translate.  Since all cell multiplicities are nonnegative, the
maximum interval weight is therefore exactly the maximum full-block weight
used in `(6)`, including strict endpoints.  Exact enumeration makes `(7)`
positive for fifteen cases; its weakest lower bound is `1,723` cells.

## 5. The three denominator-two equality cases

The other three cases have `d_1=d_2=2`, on rows `171,316,427`, with
fixed-safe cell counts `7,992`, `8,420`, and `9,420`.  Here `(7)` reaches
equality and gives no certificate.  The fact that `(Z/2Z)^*` is a singleton
does not identify the two local safe sets, because their heights and hence
their affine slopes may differ.

For `z=L/2+hL`, on one fixed-safe cell `c` the local phase is

```text
y |-> c/2+(h+1/2)y  (mod 1).                            (8)
```

Across `0<=y<=1`, the image interval has length `h+1/2`.  Its first `h`
complete periods contain danger mass `h/7`; its remaining half-period meets
at most one radius-`1/14` danger interval and contributes at most `1/7`.
Thus one high removes at most

```text
(h+1)/(7(h+1/2))<=2/7                                   (9)
```

of local-coordinate mass.  Two highs remove at most `4/7`, so their common
safe set has mass at least

```text
3/7=39/91>36/91.                                        (10)
```

This again contradicts THM-2941 `(25h)`.  All twenty-one two-high cases and
all `518` passports are now closed.

## 6. Equality boundary and the rejected shortcut

A centered cell count is not enough.  At the distinguished address `c/L`,
the projected local coordinate is zero, and zero lies in every aligned danger
set.  Thus a high-drift-safe grid address can coexist with full projected
containment.  The proof must retain every local coordinate and maximize over
translated bands, as in `(6)--(7)`.  This is recorded as MISTAKE-372.

Likewise, a zero lower bound in `(7)` proves nothing, and a unique residue
unit does not fix the height-dependent slope.  That is why the three
denominator-two rows require the measure argument `(8)--(10)`.

## 7. Exact verification and scope

The companion reconstructs the inherited queue, checks every exact rational
status packet on all sixteen columns, enumerates the residual high masks,
replays all `6,126` located-torsion cases, computes weighted translated cyclic
windows in the fifteen strict two-high cases, and freezes the three
denominator-two rows.  Its terminal and cumulative semantic digests are

```text
10165a2c7a68387bff2be57501b7a2f115126a6855a89a5d9e4841bd071eea95,
3f332da31cc80395396aae575ae8aaffe48e82513a7a3fcf56d72b046ff0b63d. (11)
```

Independent ordinary and `python3 -O` four-worker replays byte-match the
stored output.  The source AST contains no `assert` node and no floating-point
literal.  Dependency files and hashes are pinned inside the transcript.

Conclusion `(1)` is confined to the declared projected `k=3,z1=216` atlas.
Intrinsic cost is only a deterministic queue order; it has no exclusion
force.  The next exact target is the `gcd72/L720720` family.  Physical entry,
arbitrary `k`, the rung, and LRC(14) remain open.
