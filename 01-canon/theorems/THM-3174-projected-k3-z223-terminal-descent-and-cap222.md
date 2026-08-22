---
id: THM-3174
title: "Projected-k3 z223 terminal descent and cap222"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/codex-thm3088-push-2026-08-02
audit: >
  An independent hostile audit reconstructed all 65 rows, 1,767 exact Farkas
  exclusions, seven residual bodies, 125 terminal masks, fixed-modulus typing,
  and every carrier closure.  The immutable canonical companion passed normal
  and optimized replay with identical 22-line output, no stderr, and exact LF
  script/output hashes.
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3114-projected-k3-z227-screen-and-z226-terminal-double-layer-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.py
output: 05-knowledge/results/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.out
hash_basis: LF-normalized bytes
script_sha256: 3ac95e58861078828671e606e1c97705ac4def2728019ac21bd78eba0b9f1c18
output_sha256: 4edc6c7efb64ce1062aa863a5f4c17e4735f42bf7b4ac371c6bca660346abe43
semantic_sha256: f6493459d78f6af58b999bd54c8b72e2a8f989c533a1e5a91d5e6bc337352ec7
---

# THM-3174 -- projected-k3 z223 terminal descent and cap222

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact layer and full-order anchor

In the pinned THM-2941 projected `k=3` necessary atlas, the complete next
layer after THM-3139 is

```text
z1=223: 65 rows = 54 wall + 11 order,                    (1)
row SHA256 58454a33efe0b74fddeddb498d1518033ccf74064b19ec611d63a223f6065566.
```

For every row with body `E` and ruler `L`,

```text
L=14 lcm(E),             gcd(223,L)=1,             first_d=L. (2)
```

The inherited evaluator checks each selected denominator tuple `ds`, including
its first label.  Consequently every residual state also satisfies

```text
lcm(ds)=L.                                                   (3)
```

This is an exact full-order anchor, not a variable-modulus or reduced-order
relaxation.

## 2. Complete screen

The exact ray/status screen gives

```text
2,986 states
 = 1,094 crude + 1,767 exact-status + 125 residual.       (4)
```

All `1,767` status exclusions are independently rebuilt legacy exact Farkas
contradictions; no direct normalized certificate is used.  The complete
solver-witness-free screen and residual-bank digests are

```text
88220270f2b1b43aeb0fe237472b31e971ad0a275ef904301f08a89cd748d2eb,
822fcd1cd347acbec8773b95a04bb5e17c598a8225cc9ebc4ab6ce7aeca96b83. (5)
```

All eleven order rows have zero residual.  The seven residual wall bodies are

```text
(1,5,6,9,12,14):       3 masks,
(1,5,9,11,12,14):     26 masks,
(1,6,8,10,13,14):     21 masks,
(1,8,10,11,12,14):    35 masks,
(1,8,10,12,13,14):    29 masks,
(2,6,8,10,12,14):      1 mask,
(2,8,10,11,12,14):    10 masks.                         (6)
```

## 3. Terminal descent

Every body in `(6)` has a strictly positive duplicate-permitting two-high gap;
the smallest is

```text
4592897/1781319540 > 0                                  (7)
```

on `(1,5,6,9,12,14)`.  Because each residual is a wall row, the inherited wall
gate requires at least one high suffix slot.  The positive two-high gap forbids
two or more, so every actual residual assignment has exactly one high slot.
Replacing that slot by its ray supremum enlarges the actual feasible bank.

There is exactly one enumerated one-high case per residual mask.  The complete
terminal census is

```text
residual bodies / masks:      7 / 125,
zero-high hostile relaxations:    115,
one-high cases:                   125,
body-local low-label sets:         13,
coarse / exact closures:       114 / 11,
max-gap fallback / failures:      0 / 0,
minimum support slack:                1.                 (8)
```

Thus all `125` one-high cases close on complete projected cells.  The terminal
semantic, closure, and case-vector digests are

```text
54ca859fd914c26f52b78b93df89c29b2fe85b13983935de372d63fb847ec824,
35b165deda53bd549fec3cf6ad4b70ce8f7c02b26c639d386cde8f80b6dfe071,
11ba29c9a3f933a495709a4449146d9fbce064b9978160724bdf324a52623774. (9)
```

## 4. Consequence and next layer

The layer contains `65` rows, so

```text
374090 - 65 = 374025,               z1 <= 222.           (10)
```

The next layer is much wider:

```text
z1=222: 219 rows = 199 wall + 20 order,
row SHA256 68ea5b7128d3deafd4bb9a14f6d1d09ba867a5f117645a3e1f4f3d03e761de73. (11)
```

It is occupied and is not screened by this theorem.

## 5. Evidence and scope

The canonical companion recomputes all `65` rows directly from the maintained
atlas.  It binds the complete index set, `(1)--(5)`, all seven individual mask
digests, `(8)--(9)`, and the next-layer metadata.  Every raw Farkas certificate
is checked exactly, while the persisted screen digest binds only `row[:19]`.
The terminal closure digest deliberately omits the chosen duplicate-gap
maximizer witness.  These are the solver-witness-free conventions required by
MISTAKE-331/333.

This theorem concerns only the maintained projected `k=3` necessary atlas.
It does not classify physical covers outside that projection, alter `k<=1`,
close the final rung, or prove LRC(14).

Reproduction:

```text
python 04-computation/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.py
python -O 04-computation/lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.py
```
