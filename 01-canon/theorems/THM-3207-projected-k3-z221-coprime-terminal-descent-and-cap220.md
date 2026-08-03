---
id: THM-3207
title: "Projected-k3 z221 coprime terminal descent and cap220"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/projected-k3-z221/2026-08-02
depends_on:
  - THM-3179-projected-k3-z222-composite-divisor-square-terminal-descent-and-cap221
  - THM-3174-projected-k3-z223-terminal-descent-and-cap222
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3114-projected-k3-z227-screen-and-z226-terminal-double-layer-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
script: 04-computation/lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.py
output: 05-knowledge/results/lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.out
hash_basis: LF-normalized bytes
script_sha256: b3d1f0c451087017c1363d42be8789df78d5eec7db7a05b49dc5ca9e194f2091
output_sha256: 1e4b01f6e5bfb179ad5fe6ba786124ff2c9fdf3833eed4917c0b3ce0abb7b76d
semantic_sha256: aad27ae040933935c6eacfd2801cef1644df3de3a8f031aebf7d1b238ac35e74
audit: >
  An independent hostile audit reconstructed all 90 atlas rows, all 4,420
  exact Farkas exclusions, the 22 residual bodies and 387 individually bound
  masks, the coprime fixed-modulus identities, and every terminal bank.  It
  confirmed the sharp positive gap, all 392 one-high closures, the ledger and
  cap arithmetic, and the occupied z220 metadata.  Canonical normal and
  optimized replays agree byte-for-byte with the stored transcript and the
  declared LF hashes.
---

# THM-3207 -- projected-k3 z221 coprime terminal descent and cap220

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact layer and the coprime reset

In the pinned THM-2941 projected `k=3` necessary atlas, the complete layer
after THM-3179 is

```text
z1=221: 90 rows = 83 wall + 7 order,                      (1)
row SHA256 ab4f5b7fe3b5e1330d6e51dfb150527fa27cac56b755ed5e2ae2522f0f27ceb4.
```

For every row with body `E` and ruler `L`, the exact atlas gives

```text
L=14 lcm(E),                 gcd(221,L)=1,                (2)
first_d=L.
```

Thus every residual denominator tuple `ds` satisfies

```text
D=lcm(ds)=L.                                             (3)
```

This is the full-order coprime anchor of THM-3174, now at a layer whose
screen is much smaller.  It is not the composite `L/6` anchor of THM-3179.

## 2. Complete exact screen

The complete ray/status screen gives

```text
8,705 states
 = 3,898 crude + 4,420 exact-status + 387 residual.      (4)
```

All `4,420` exact-status exclusions are independently rebuilt legacy exact
Farkas contradictions; no direct normalized certificate is used.  The seven
order rows contribute only

```text
32 states = 32 crude + 0 status + 0 residual.            (5)
```

Hence every residual is a wall row.  The solver-witness-free screen and
residual-bank digests are

```text
820ff7d9bc0d593434b672f4bbdbd08d7c7a164b063503e2e0a98d793ce5f82a,
c8b5bb5bcda2034f95968022374e5b455c73e487f81873624a0a4d4c9020f38d. (6)
```

The residual has exactly `22` bodies and `387` masks.  The canonical
companion binds every body's mask count and mask digest separately, as well
as the aggregate digest in `(6)`.

## 3. What happened to the `2`/`3` square

At `z1=222`, THM-3179's quotient `D/(L/6)` retained two independent
idempotent valuation bits and realized an edge of the Boolean join square
`B_2`.  Equations `(2)--(3)` show that at `z1=221` the corresponding quotient
is the singleton

```text
D/first_d = 1.                                          (7)
```

So neither a `V_4` torsor nor an order-sensitive `C_2*C_3` restoration word
survives this level.  The coprime step resets the valuation-carry coordinate
rather than transporting it.  The important hostile boundary is that this
does **not** close the geometry: `387` wall masks remain after the reset.
Coprimality removes the divisor-square sidecar; it does not replace the
screen or the one-high carrier argument.

## 4. Terminal descent

Every residual body has a strictly positive duplicate-permitting two-high
gap.  The smallest is

```text
1931203/776417642 > 0                                   (8)
```

on `E=(4,6,8,10,11,14)`.  Because every residual is a wall row, the wall
gate requires at least one high suffix slot.  The positive two-high gap
forbids two or more, so every actual residual assignment has exactly one
high slot.  Replacing that slot by its ray supremum enlarges the actual
feasible bank.

The complete terminal census is

```text
residual bodies / masks:       22 / 387,
zero-high hostile relaxations:      319,
one-high cases:                     392,
body-local low-label sets:           45,
coarse / exact closures:        318 / 74,
max-gap fallback / failures:       0 / 0,
minimum support slack:                 1.                (9)
```

All `392` one-high cases close on complete projected cells.  The terminal
semantic, closure, and case-vector digests are

```text
93e09fce6debb80536e4c9f40e3c6f6cb79f993f542cd939884757794ac4e924,
6623b2b236e2fb96f003e13e549e66d621b9ad9596c3bcf182ab7dc03b71df03,
6c0bb93eee6fb06bed9ff192b54d3eb7a3d508161d0b0313e3917834ce18fb7b. (10)
```

## 5. Consequence and next layer

The layer contains `90` rows, so the proved consequence is

```text
373806 - 90 = 373716,               z1 <= 220.           (11)
```

The next layer is occupied:

```text
z1=220: 289 rows = 249 wall + 40 order,
row SHA256 bc75fef8aebdc6b350957ea96842b38a24bad2174200646ebaf2817844515799. (12)
```

It is not screened by this candidate.

## 6. Evidence and scope

The checkpoint-independent canonical companion recomputes all `90` screen
rows and all `22` terminal bodies directly from the maintained atlas.  It
pins THM-3174's source, output, and semantic hashes, binds all `22`
individual residual mask banks, checks every raw Farkas certificate exactly,
and verifies `(2)--(3)` on every row and residual mask.  Under
MISTAKE-331/333, the screen digest binds only `row[:19]`; the terminal closure
digest omits the chosen duplicate-gap maximizer witness.  Both are
solver-witness-free.

This theorem concerns only the maintained projected `k=3` necessary
atlas.  Coprimality alone does not close a wall row.  The result does not
classify physical covers outside that projection, alter `k<=1`, close the
final rung, or prove LRC(14).

Reproduction:

```text
python 04-computation/lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.py
python -O 04-computation/lrc14_j7_k3_z221_coprime_terminal_descent_cap220_thm3207.py
```
