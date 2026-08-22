---
id: THM-3139
title: "Projected-k3 z225 terminal and z224 screen double-layer descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every projected
  k=3 row at z1=225 and z1=224 is empty.  The projected cap is z1<=223 and
  the scalar ledger is 374090.  No LRC(14) claim is made.
source: root/codex-thm3088-push-2026-08-02
audit: >
  Independent screen and terminal reconstructions verified all 78 z225 rows,
  24 residual bodies, 1,898 residual masks, positive two-high gaps, 3,025
  one-high carrier cases, and the complete four-row z224 screen.  The z224
  audit rebuilt all 10 exact Farkas certificates and found zero residual.
  The canonical companion recomputes both layers and binds only
  solver-witness-free row digests under MISTAKE-331/333.  Fresh normal and
  optimized replays match the stored transcript and declared hashes.
depends_on:
  - THM-3114-projected-k3-z227-screen-and-z226-terminal-double-layer-descent
  - THM-3111-projected-k3-z230-exact-screen-and-compressed-complete-cell-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py
output: 05-knowledge/results/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.out
hash_basis: LF-normalized bytes
script_sha256: 92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36
output_sha256: 4f0ae623d0134406f83f466c0a9f1353525a992acb8ee82c2e92f2f0537f32c5
semantic_sha256: 1eb7e8faefafe41fd6f6cbc108ba09295dbcbb7fe8016d7ec0cbdc258adf358a
---

# THM-3139 -- projected-k3 z225 terminal and z224 screen double-layer descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact adjacent layers

In the pinned THM-2941 projected `k=3` necessary atlas,

```text
z1=225: 78 rows = 78 wall + 0 order,
z1=224:  4 rows =  3 wall + 1 order.                      (1)
```

The complete ordered-row digests are

```text
9a4869fbb3886b058bdb2b60209b1fc88ab4c82ccbc7c56eba42b2d118ee8719,
b42d3f6d42e76a10f89ba7d1fda6f76884081cba9cf32ca4fcedba4757811cb4. (2)
```

The first layer requires the inherited terminal carrier; the second dies at
the exact ray/status screen.

## 2. The z225 screen and terminal

The exact `z1=225` screen gives

```text
30,615 states
 = 15,904 crude + 12,813 exact-status + 1,898 residual.   (3)
```

All `12,813` status exclusions are independently rebuilt exact Farkas
contradictions; none uses the direct normalized branch.  The solver-witness-free
screen digest is

```text
75de10a09defef3d4d388b2eb21511d17f80ddffee274318677bc7cc7bc71cd7. (4)
```

The residual consists of exactly `24` bodies and `1,898` masks.  Every body
has a strictly positive duplicate-permitting two-high gap.  The inherited
wall gate requires at least one high suffix slot, so every actual residual
assignment has exactly one high slot.  The complete terminal census is

```text
zero-high hostile relaxations: 1,754,
one-high cases:              3,025,
body-local low-label sets:     163,
coarse cardinality closes:    2,822,
exact cardinality closes:       203,
max-gap fallbacks / failures:    0 / 0,
minimum support slack:              1.                    (5)
```

The zero-high cases are deliberately enlarged hostile relaxations excluded
by the wall gate.  The one-high bank replaces the actual high band by its ray
supremum, again enlarging the feasible assignment set.  Every one of its
`3,025` cases closes on a complete projected cell, with no max-gap fallback.
The solver-witness-free closure and case-vector digests are

```text
54222f669aae83b2f217b667520ed9b9e27059ebbb8f0b31d177232cbade226b,
151c16219d8b079f910963e19ff3b813f29584fdae2ecc8ec19969c469528cb6. (6)
```

## 3. The z224 exact screen

The complete next layer gives

```text
12 states = 2 crude + 10 exact-status + 0 residual.       (7)
```

This includes all three wall rows and the sole order row.  Every returned
status certificate is checked over exact rationals; the direct/legacy count
is `0/10`.  No high-slot or terminal implication is needed.  The
solver-witness-free screen digest is

```text
f6928f10e464ee694ebce756f7a6a2cb43360a7cd3bd121be4b6f5b98f91b21d. (8)
```

The ray quotient is a necessary upper relaxation: every actual projected
assignment maps to an enumerated state.  Therefore its zero residual proves
the entire four-row layer empty.

## 4. Ledger, evidence, and scope

The two layers contain `78+4=82` rows.  Hence

```text
374172 - 82 = 374090,              z1 <= 223.             (9)
```

The next layer is occupied:

```text
z1=223: 65 rows = 54 wall + 11 order,
row digest 58454a33efe0b74fddeddb498d1518033ccf74064b19ec611d63a223f6065566. (10)
```

All raw Farkas certificates are verified, but under MISTAKE-331/333 each
persisted screen digest binds only `tuple(row[:19])`, never a solver-selected
dual representative or contradiction magnitude.  The terminal closure digest
likewise omits the chosen deterministic duplicate-gap maximizer witness.  The
top-level semantic hash binds these solver-witness-free digests and
basis-invariant counts.

This is only a theorem about the maintained projected `k=3` necessary atlas.
It does not classify physical covers outside that projection, alter `k<=1`,
close the final rung, or prove LRC(14).
