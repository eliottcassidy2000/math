---
id: THM-3052
title: "Projected k3 z242 gap241 z240 compositional descent"
status: CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT
source: codex-lrc-z242-z240-composition-2026-08-01
depends_on:
  - THM-3041-projected-k3-z243-first-below-floor-cardinality-descent
  - THM-3033-projected-k3-z246-to-z244-descent-and-z243-high-floor-addendum
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
script: 04-computation/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.py
output: 05-knowledge/results/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.out
script_sha256: 548cb4f0ff095bb73e29cd985b494e80821e0a1b507af302f1cf5729f1282227
output_sha256: c362701c055fd1ed8888f6b0597294c88f72d985913bca880157d47494bb42c7
semantic_sha256: 8e2e6a1ed90033f3d85ec9cd03acf326ffb0e67453bbc206e8f2ece843c7cd68
hash_basis: LF-normalized bytes
---

# THM-3052 -- projected k3 z242 / gap241 / z240 compositional descent

**CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

This file remains outside the proved dependency graph until an independent
hostile audit and explicit status promotion.  No navigation or ledger file is
changed by this candidate.

## Statement

In the lossless projected `k=3` body atlas inherited from THM-3041 and
THM-3033, the exact compositor closes

```text
all 20 occupied z_1=242 rows,
the zero-row atlas gap z_1=241,
all 52 occupied z_1=240 rows.                           (1)
```

The occupied-row decrement in `(1)` is `20+52=72`; the gap costs zero.  Thus,
**if this candidate is promoted**, THM-3041's proved projected necessary-row
ledger and cap would update by the exact arithmetic

```text
374900-72=374828,
projected k=3 cap: z_1<=239.                            (2)
```

Equation `(2)` is an audit target, not a live ledger or navigation mutation
while the theorem has candidate status.  This is not LRC(14), says nothing
about `k<=1` or the final rung, and closes only the named projected necessary
sector.

## 1. Pinned disjoint universe and the `241` gap

The compositor pins LF-normalized THM-3033 source/output/semantic hashes and
the THM-2941 projected-body atlas hash.  It parses every one of the `6,060`
`row=` records by the structured row grammar.  A separate literal `;z1=...;`
scan gives the same neighboring census:

```text
z_1=243: 154 rows = 151 below floor + 3 at/above floor;
z_1=242:  20 rows =  16 below floor + 4 at/above floor;
z_1=241:   0 rows;
z_1=240:  52 rows =  49 below floor + 3 at/above floor;
z_1=239:   4 rows =   2 below floor + 2 at/above floor. (3)
```

The positive `242/240/239` neighbors prevent a parser miss from masquerading
as an empty layer.  The combined occupied bank has exactly

```text
72 rows = 65 wall + 7 order,                         (4)
```

with row-order digest
`c7f760f8571b2510ecd218a0172e2a5ec13b0b58d65c015446b61c19257c99e5`.
The `z_1=241` conclusion is exact absence from the pinned atlas universe, not
an inference from a zero-measure carrier.

## 2. Exact screen and first failed implication

Every occupied row is passed through the unmodified THM-3033 attained-state
generator and its enlarged common-status relaxation.  The exact combined
partition is

```text
states     9449
crude      4783
status     4629
residual     37.                                       (5)
```

The seven order rows contribute `24=7+17+0` states in
`(states,crude,status,residual)` order.  Hence the crude and Farkas screen
closes `9,412/9,449`, not the whole bank.  The first invalid implication would
be to declare the remaining 37 relaxed states empty without a terminal
argument.

The exact level profiles are

| level | rows `(wall,order)` | states | crude | status | residual | residual bodies |
|---:|---:|---:|---:|---:|---:|---:|
| 242 | `20 (16,4)` | 1,430 | 279 | 1,125 | 26 | 4 |
| 240 | `52 (49,3)` | 8,019 | 4,504 | 3,504 | 11 | 2 |

At `z_1=242`, the four order rows have profile `(21,5,16,0)`; at `240`,
the three order rows have `(3,2,1,0)`.  The frozen transcript additionally
prints the exact profile and stage/residual digest of every one of the 72
body rows.

## 3. Terminal closure

For every residual body the inherited duplicate-permitting two-high upper
bound has a strict positive gap.  Thus a genuine residual packet has exactly
one high label.  THM-2984's pointwise complete-cell projection then gives an
actual safe residue at every local coordinate, while an aligned translated
high-danger band has the strict finite cardinality cap checked by the
companion.

The six terminal profiles are:

| `z_1` | body `E` | residual | positive two-high gap | cardinality cases | minimum slack |
|---:|---|---:|---|---:|---:|
| 242 | `(1,8,10,11,12,14)` | 8 | `4936302911/1123781503845` | 8 | 13 |
| 242 | `(2,7,9,11,12,13)` | 12 | `267387356141/62645926164822` | 12 | 367 |
| 242 | `(2,8,10,11,12,14)` | 4 | `53389945853/15058770887160` | 4 | 3360 |
| 242 | `(2,9,10,11,12,14)` | 2 | `92484982943/24066398771415` | 2 | 3960 |
| 240 | `(1,6,8,10,13,14)` | 10 | `9392296609/3809475306000` | 10 | 68 |
| 240 | `(2,9,10,11,12,14)` | 1 | `72297067/19149710580` | 1 | 6 |

The aggregate terminal census is

```text
six residual bodies;
six strict two-high exclusions;
36 zero-high stress controls;
37 one-high cases;
37 cardinality closures;
zero max-gap cases, failures, or unit checks.            (6)
```

The zero-high counts in `(6)` are stress controls: the inherited projected
wall excludes zero-high packets.  The terminal proof is the strict
two-high-plus-one-high direction, not a new scalar conclusion from that
control count.

## 4. Pointwise and evidence boundaries

The common-status screen is a necessary upper relaxation: every literal
pointwise cover maps into it.  At the terminal, every complete-cell coordinate
retains an actual safe residue, so a proper translated open danger union
cannot contain the full projected section.  No step uses the invalid
implication `mu(Safe)=0 => Safe is empty`.

Checkpoint files bind an algorithm schema, all inherited hashes, every frozen
census, and every row-order digest into fingerprint

```text
97c94ba773cdbf3aaadb04eae7e800f319670ee26f2d92595917ec4354de7423. (7)
```

Fresh ordinary and optimized runs use disjoint initially empty checkpoint
directories.  Each rebuilds all `72` screen rows and all six terminal rows;
optimized execution cannot erase truth-bearing checks because the companion
uses explicit `require` calls rather than Python `assert`.

Canonical commands are

```text
python 04-computation/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.py --processes 4 --checkpoint-dir <fresh-normal-dir>
python -O 04-computation/lrc14_j7_k3_z242_gap241_z240_compositional_descent_thm3052.py --processes 4 --checkpoint-dir <fresh-optimized-dir> --output <optimized-output>
```

The exact semantic digest is
`8e2e6a1ed90033f3d85ec9cd03acf326ffb0e67453bbc206e8f2ece843c7cd68`.
Both fresh checkpoint banks contain `72+6=78` envelopes and have the same
semantic digest
`08002b99f7d5e9dde253ed923ae888feb7500906fde75fb6774e58f8533f7598`.
The canonical transcript has `22,683` LF bytes and `96` lines; ordinary,
optimized, and stored bytes agree exactly.
