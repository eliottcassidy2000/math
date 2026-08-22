---
id: THM-3361
title: "Projected k=3 L720720 one-high translated-residue closure"
status: >
  PROVED analytic mechanism + FINITE-EXACT in the declared projected atlas.
  The three gcd72/L720720 rows at k=3,z1=216 have empty necessary projected
  screens.  The inherited high gate and three positive duplicate-two-high
  gaps force exactly one high drift; all 820 one-high cases close because
  their fixed-safe residue support is larger than every translated open
  one-seventh band.  A direct located-torsion route independently closes the
  same 820 cases.  The projected ledger falls 372917 to 372914, wall rows 113
  to 110, and complete families 13 to 12; the cap remains 216.  This proves
  no physical entry, arbitrary-k statement, rung, or LRC(14).
source: codex-2026-08-14-l720720-projected-family
audit: >
  dependency-pinned exact reconstruction, translated-equality hostile,
  independent order-two/three located-torsion implementation, no assert or
  floating-point syntax, and byte-matching ordinary/optimized replays
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-3351-projected-k3-z216-sixteen-family-translated-two-high-closure
related:
  - THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure
script: 04-computation/lrc14_j7_k3_z216_L720720_one_high_translated_residue_thm3361.py
output: 05-knowledge/results/lrc14_j7_k3_z216_L720720_one_high_translated_residue_thm3361.out
script_sha256: 52d70825de6199182c6a8dc3c2674cb0f6f92cfdec8a2b49e22c478307fc3285
output_sha256: 22732fcf2ed2d094834dbbafa9893129856de52922892881f485e2e2b34cdede
hash_basis: LF-normalized bytes
---

# THM-3361 -- projected `k=3`, `z1=216`, `L=720720` closure

**PROVED analytic mechanism + FINITE-EXACT in the declared projected
atlas.**

## 1. Statement

Reconstruct the pinned necessary projected `k=3,z1=216` atlas and the queue
state left by
[THM-3351](THM-3351-projected-k3-z216-sixteen-family-translated-two-high-closure.md).
The next complete intrinsic-cost family is `gcd72/L720720`, of cost
`50,450,400`, on exactly these three wall rows:

```text
row 191   E=(1,5,8,9,11,13)    22 components,
row 228   E=(1,8,9,10,11,13)   26 components,
row 332   E=(2,5,8,9,11,13)    22 components.             (1)
```

Every row in (1) has an empty necessary projected screen.  Therefore, inside
this declared projected universe,

```text
projected ledger       372917 -> 372914,
z1=216 wall rows          113 -> 110,
complete families          13 -> 12,
projected k=3 cap         216 unchanged.                  (2)
```

Intrinsic cost only orders the queue.  It is not an exclusion inequality.
Conclusion (2) does not restore a physical entry, endpoint origin, owner,
phase, current, arbitrary `k`, the rung, or LRC(14).

## 2. Exact inherited screen

The ray/common-status screen partitions the exact denominator states as

```text
16,900 = 8,947 crude capacity exclusions
         +7,234 exact common-status exclusions
         +  719 residual passports.                      (3)
```

The rowwise partitions are

```text
row 191:  8,288 = 4,002 + 3,899 + 387,
row 228:  8,135 = 4,706 + 3,109 + 320,
row 332:    477 =   239 +   226 +  12.                   (4)
```

All `7,234` status exclusions are inherited full-table exact-rational
witnesses.  Every residual denominator tuple contains the first denominator
`10010` exactly once and has lcm `720720`.  The companion pins and checks the
source and semantic hashes of THM-3351 and the earlier atlas, status, and
torsion engines before reconstructing (3)--(4).

## 3. At least one high and at most one high

By THM-2941's strict scalar wall, every projected completion has largest later
drift

```text
M > 13L/132.                                             (5)
```

For `L=720720`, put

```text
H=floor(13L/132)+1=70,981.                              (6)
```

Since the projected first drift is `216<H`, (5) forces at least one of the
three later drifts to be high.  Exact scalar arithmetic still leaves `682`
zero-high hostile records; (5), rather than their scalar score alone, excludes
them.

The duplicate-permitting two-high exclusion gaps are respectively

```text
row 191: 179568616148 / 133555340093175,
row 228: 445315966867 / 211281696826200,
row 332: 8767934779397 / 3096873978798600.               (7)
```

All three fractions in (7) are strictly positive.  Thus even after allowing
the two selected high labels to use the same denominator, no completion has
two or more high drifts.  Combining (5)--(7), every hypothetical completion
has exactly one high drift.

Enumerating the high choice and the two low labels gives `820` one-high cases
over all `719` residual passports and `54` rowwise-distinct low-pair classes
(`26+26+2`):

```text
row 191: 462 cases / 387 passports / 26 low pairs,
row 228: 346 cases / 320 passports / 26 low pairs,
row 332:  12 cases /  12 passports /  2 low pairs.       (8)
```

## 4. The translation-uniform residue lemma

Fix one case from (8).  Let `d|L` be its high denominator and let `C` be the
complete cells that are fixed-safe for the first drift and both low labels.
Discarding multiplicity, define the residue support

```text
R_d(C)={c mod d : c in C} subset Z/dZ.                  (9)
```

Write a concrete high drift in the inherited normal form

```text
z=(L/d)u+hL,               gcd(u,d)=1.                  (10)
```

For every fixed projected local coordinate `y`, THM-2984's complete-cell
projection law identifies the high-danger cells with those satisfying

```text
u c mod d in J_y,                                         (11)
```

where `J_y` is some translated strict-open cyclic interval of length `d/7`.
The height `h` and coordinate `y` change the translation of `J_y`, but not its
length.  Such an interval contains at most `ceil(d/7)` residue classes, and
multiplication by the unit `u` permutes `Z/dZ`.  Consequently

```text
|R_d(C)| > ceil(d/7)                                     (12)
```

implies that, for every `y`, at least one fixed-safe cell is also safe from
the high drift.

Thus (12) makes the projected safe set the whole local-coordinate circle.
That contradicts THM-2941's necessary completion bound
`mu(P_(E,Z))<36/91`.  Notice the quantifier: the surviving cell may depend on
`y`; the conclusion needed here is pointwise coverage of every `y`, not one
globally selected cell.

The exact companion verifies (12) in all `820` cases.  The weakest margin is
one: for a denominator-two case the support has two residues while the danger
capacity is one.  The rowwise fixed-safe cell ranges are `124,168--130,496`,
`130,334--139,040`, and `136,112--138,346` for rows `191`, `228`, and `332`.
Hence all cases in (8), and therefore all passports in (3), close.

The strict inequality in (12) is necessary for this argument.  At `d=14`, the
support `{2,3}` is disjoint from the aligned danger interval, yet the translated
open interval `(3/2,7/2)` contains both residues.  Here
`|R_d(C)|=ceil(d/7)=2`, so aligned safety or equality gives no
translation-uniform survivor.  This hostile control blocks both tempting
shortcuts.

## 5. Independent located-torsion closure

The companion also performs a second terminal search without using the
cardinality conclusion (12).  In every one of the `820` cases it directly
finds two fixed-safe cells whose residue difference has effective order two
or three modulo `d`.  Their high phases differ by at least `1/3`, independently
of the translated local coordinate, so one translated open danger band of
length `1/7` cannot contain both phase values.

The independent direct census is

```text
effective order 2: 758 cases,
effective order 3:  62 cases.                           (13)
```

The inherited located-torsion engine independently reports `752` order-two
and `68` order-three effective witnesses; two of those begin as qualifying
order-six differences and reduce to effective order three.  Both routes close
all `820` cases, although their chosen pairs need not agree.  This is a
cross-check of the terminal, not an additional physical statement.

## 6. Verification, boundary, and conclusion

The exact companion:

1. reconstructs all `16,900` states and all sixteen common-status columns;
2. freezes every residual passport and one-high packet by semantic digest;
3. checks the three exact gaps (7), all `820` residue inequalities (12), and
   all `820` direct located-torsion witnesses (13);
4. includes the translated-equality hostile above;
5. contains no `assert` node or floating-point literal; and
6. pins every inherited source/output dependency.

Ordinary and optimized three-worker replays byte-match the stored output.  Its
semantic digest is

```text
927e52562c64c35833428ea735829ba236b7dcaa5f9e4d2db24f750e1a2db77c. (14)
```

Therefore the three rows (1) have empty necessary screens and (2) follows.
The result is confined to the declared projected `k=3,z1=216` atlas.  It does
not prove physical entry, arbitrary `k`, the rung, or LRC(14).
