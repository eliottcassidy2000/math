---
id: THM-3218
title: "Projected-k3 z220 valuation-product terminal descent and cap219"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/projected-k3-z220/2026-08-02
depends_on:
  - THM-3207-projected-k3-z221-coprime-terminal-descent-and-cap220
  - THM-3179-projected-k3-z222-composite-divisor-square-terminal-descent-and-cap221
  - THM-3174-projected-k3-z223-terminal-descent-and-cap222
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
script: 04-computation/lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.py
output: 05-knowledge/results/lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.out
hash_basis: LF-normalized bytes
script_sha256: 97c55e0a9acc9b42c7c7a856889bd2b8fbd49854ac39a1fe974a488f63a3cff8
output_sha256: af7fe69ba68579611f171b9b2252f9a9923bdd2785720a6cb1b208a611ee0565
semantic_sha256: ed67ff55840d07965a6cf8a478c0fc49cab637511ad3d684b2de8831dfcfcd94
audit: >
  An independent hostile audit rebuilt the complete screen from the generic
  pre-candidate z220 checkpoints and recomputed all 31 terminal bodies
  directly through THM-3139 under normal and optimized interpreters.  It
  reproduced every census, Farkas split, residual bank, valuation quotient,
  terminal constant, digest, minimum gap, ledger decrement, and next-layer
  hash.  The canonical heavy-first normal and optimized replays agree
  byte-for-byte with the stored output and the declared LF hashes.
---

# THM-3218 -- projected-k3 z220 valuation-product terminal descent and cap219

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact layer and valuation strata

In the pinned THM-2941 projected `k=3` necessary atlas, the complete layer
after THM-3207 is

```text
z1=220: 289 rows = 249 wall + 40 order,                    (1)
row SHA256 bc75fef8aebdc6b350957ea96842b38a24bad2174200646ebaf2817844515799.
```

For every row with body `E` and ruler `L=14 lcm(E)`, put

```text
g=gcd(220,L),                    first_d=L/g.              (2)
```

The complete row stratification is

```text
g       rows      wall      order
4         72        52         20
20       183       163         20
44         5         5          0
220       29        29          0.                         (3)
```

Thus this layer is neither the coprime reset of THM-3207 nor one fixed
Boolean square.  It is a product-of-valuation-chains layer with four pointed
divisor fibres.

## 2. Complete exact screen

The complete ray/status screen gives

```text
96,780 states
 = 47,054 crude + 49,356 exact-status + 370 residual.     (4)
```

The wall and order pieces are separately

```text
wall:   95,673 = 46,043 + 49,260 + 370,
order:   1,107 =  1,011 +     96 +   0.                   (5)
```

The status bank uses `12` direct normalized and `49,344` legacy exact Farkas
certificates; the split is `12/49,248` on wall rows and `0/96` on order
rows.  All forty order rows therefore close before the terminal step.  The
solver-witness-free screen and residual-bank digests are

```text
f7168e67f3eb21529edcb73cae97815804827cfdb0b2fe7df66dbd37a5e8ba07,
71525cd0a7e90047f6e391e769142b57c37f09506d656a52d0dbde57da1502f7. (6)
```

The residual has exactly `31` bodies and `370` masks.  The companion binds
every body's mask count and mask digest separately.

## 3. The surviving valuation product

For every residual denominator tuple `ds`, define

```text
D=lcm(ds),                       q=D/first_d.              (7)
```

The complete residual quotient census is

```text
(g,q)=(4,4):       2 masks,
      (20,10):    15 masks,
      (20,20):   295 masks,
      (44,44):     5 masks,
      (220,220):  53 masks.                               (8)
```

Hence every residual quotient is the top element of its pointed divisor
lattice except for fifteen `g=20` masks occupying the lower vertex `q=10`.
Together, the occupied values `10` and `20` are the two vertices of the
ambient top `2`-cover; after normalization by `g`, that cover is the
Farey-adjacent pair `1/2,1`.  This is only a poset statement.  No transition,
common carrier, or physical edge between two masks is proved.  There is no
occupied `3`-cover at this layer.

This is the precise survivor of the recurring binary/ternary picture.  The
divisor fibre is a median partial cube under valuation threshold bits, but
its natural restoration maps are joins `q -> lcm(a,q)`.  Every nonidentity
join translation is noninvertible.  Therefore `(8)` is not a `V4` torsor and
does not furnish an action of `C2*C3 = PSL2(Z)`.  To obtain that stronger
object one would need ordered supplier/carry data and actual common-carrier
bijections; the unpointed square or Farey determinant alone forgets exactly
those coordinates.

## 4. Terminal descent

Every residual body has a strictly positive duplicate-permitting two-high
gap.  The smallest is

```text
28369051/33158840715 > 0                              (9)
```

on `E=(1,4,9,10,12,14)`.  Every residual is a wall row, so the wall gate
requires at least one high suffix slot.  The positive two-high gap forbids
two or more; every actual residual assignment therefore has exactly one high
slot.  Replacing that slot by its ray supremum enlarges the actual feasible
bank.

The complete terminal census is

```text
residual bodies / masks:       31 / 370,
zero-high hostile relaxations:      326,
one-high cases:                     644,
body-local low-label sets:          163,
coarse / exact closures:        589 / 55,
max-gap fallback / failures:       0 / 0,
minimum support slack:                 1.                 (10)
```

All `644` one-high cases close on complete projected cells.  The terminal
semantic, tagged closure, and case-vector digests are

```text
4254a3d491e57cbdb3e1a9fda66a002e81bbf3989d39f877e5d2b86357073f11,
afc495948dcc329795cad66cd0aee6a8a616971a99acd67e3b3a477e353d7839,
41500dcc3c8273b9ec8e0d545ca1d6a9f8a52cfb03f15088c143d2d7473e25cc. (11)
```

## 5. Consequence and next layer

The layer contains `289` rows, so the proved consequence is

```text
373716 - 289 = 373427,               z1 <= 219.           (12)
```

The next layer is occupied but unusually small:

```text
z1=219: 16 rows = 15 wall + 1 order,
row SHA256 6a39c83d47b8080fd52e64479b886f1b6a887150f32cfbf30a676fbb0aaa54ba. (13)
```

It is not screened by this candidate.

## 6. Evidence and scope

The checkpoint-independent canonical companion recomputes all `289` screen
rows and all `31` terminal bodies directly from the maintained atlas.  It
pins THM-3207's source, output, and semantic hashes, binds all `31`
individual residual mask banks, checks every raw Farkas certificate exactly,
and verifies `(2)`, `(7)`, and `(8)` on every residual mask.  Under
MISTAKE-331/333, the screen digest binds only `row[:19]`; the tagged terminal
closure digest omits the chosen duplicate-gap maximizer witness.  Both are
solver-witness-free.

This theorem concerns only the maintained projected `k=3` necessary atlas.
The divisor product is not a physical carrier classification or a modular
group action.  The result does not alter `k<=1`, close the final rung, or
prove LRC(14).

Reproduction:

```text
python 04-computation/lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.py
python -O 04-computation/lrc14_j7_k3_z220_valuation_product_terminal_descent_cap219_thm3218.py
```
