---
id: THM-3179
title: "Projected-k3 z222 composite divisor-square terminal descent and cap221"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/codex-thm3088-push-2026-08-02
depends_on:
  - THM-3174-projected-k3-z223-terminal-descent-and-cap222
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3114-projected-k3-z227-screen-and-z226-terminal-double-layer-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
related:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
script: 04-computation/lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.py
output: 05-knowledge/results/lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.out
hash_basis: LF-normalized bytes
script_sha256: f094bb0fa276e1af2e41c2c7db907ad14477c0230ab5a37a7ad9e61c0a350c27
output_sha256: 75801a1dc7ba1086b213fba1630639a54be7c17443742335be07e0db53ccbd0e
semantic_sha256: 8f098a0c1a00f242bc6f3dd2bfa4932a6aff187186d475d889509f3a1b380ad5
audit: >
  An independent hostile audit reconstructed all 219 rows, the exact
  divisor-square strata and coupling census, all 21,285 Farkas exclusions,
  every one of the 37 residual mask banks, the B2 semilattice typing, all
  terminal digests, and the ledger decrement.  Canonical normal and optimized
  replays LF-normalize byte-for-byte to the stored 53-line transcript with
  empty stderr and the declared hashes.
---

# THM-3179 -- projected-k3 z222 composite divisor-square terminal descent and cap221

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact layer and the composite anchor

In the pinned THM-2941 projected `k=3` necessary atlas, the complete next
layer after THM-3174 is

```text
z1=222: 219 rows = 199 wall + 20 order,                    (1)
row SHA256 68ea5b7128d3deafd4bb9a14f6d1d09ba867a5f117645a3e1f4f3d03e761de73.
```

For every row with body `E` and ruler `L`,

```text
L=14 lcm(E),                     37 does not divide L.     (2)
```

Exactly `218` rows have `gcd(222,L)=6`.  The sole exception is atlas row
`8`,

```text
E=(1,2,4,8,10,14),   L=3920,   high=387,
gcd(222,L)=2,         first_d=1960.                       (3)
```

That exceptional wall row and all twenty order rows die in the exact screen.
Thus every residual row belongs to the `gcd(222,L)=6` stratum and has

```text
first_d=L/6.                                              (4)
```

This is a composite fixed-divisor anchor, not the full-order coprime anchor
used at `z1=223`.

## 2. Complete exact screen

The complete ray/status screen gives

```text
49,676 states
 = 27,341 crude + 21,285 exact-status + 1,050 residual.   (5)
```

All `21,285` exact-status exclusions are independently rebuilt legacy exact
Farkas contradictions; no direct normalized certificate is used.  The
solver-witness-free screen and residual-bank digests are

```text
146bae48b8a901cbcc3614dcba88cbeaa87f0361ef0b2a1d9b48d732628a09b5,
2ad57452ad83947b33e3245004dc44be6f7846d6870c2b5a9a3854f18c3b92b5. (6)
```

The residual has exactly `37` bodies and `1,050` masks.  The canonical
companion binds every body's mask count and mask digest separately as well as
the aggregate digest in `(6)`.

## 3. The idempotent `2`/`3` divisor square

For a residual denominator tuple `ds`, put

```text
D=lcm(ds),                         q=D/(L/6).              (7)
```

The inherited evaluator verifies

```text
L/6 divides D divides L,           q in {1,2,3,6}.         (8)
```

The two coordinates of `q` record whether some selected denominator restores
the top `2`-valuation and the top `3`-valuation of `L`.  Because `D` is an
`lcm`, these coordinates combine by idempotent Boolean OR.  The complete
residual census is

```text
q=1: 0,          q=2: 0,          q=3: 38,          q=6: 1,012. (9)
```

Every `q=3` mask occurs on the single body

```text
E=(2,6,8,10,12,14),               D=L/2.                 (10)
```

Among the `q=6` masks, `700` restore the top `2`- and `3`-valuations on
separate selected labels, while `312` contain one label carrying both.
The `38` `q=3` masks have no coupled label.  Thus the ambient quotient in
`(8)` is the Boolean join square `B_2` on the missing valuation bits, while
the realized residual support in `(9)` is only its two-element edge `{3,6}`.

This ambient square is deliberately **not**, under its internal join law, a
`V_4` torsor or an order-sensitive
`C_2*C_3` restoration history.  A join semilattice has no nontrivial internal
units: repeated restoration is absorbed, and the order of the two repairs is
forgotten.  This does not forbid the external `C_2` automorphism swapping the
two atoms of `B_2`; that swap sends the realized edge `{3,6}` to the absent
edge `{2,6}`, so it is not a symmetry of the residual bank.  Equivalently, the
screen forces the top `3`-valuation in every residual and leaves only the top
`2` bit variable.  The bank is a directed one-bit edge, not a two-generator
co-occurrence object.  THM-2197 is the analogous Boolean-deficiency warning;
THM-2597 separates a commuting order-six quotient from the free modular
product.  A restoration word would require a separate carry/path sidecar.

## 4. Terminal descent

Every residual body has a strictly positive duplicate-permitting two-high
gap.  The smallest is

```text
5424109/9842060865 > 0                                  (11)
```

on the body in `(10)`.  Because each residual is a wall row, the inherited
wall gate requires at least one high suffix slot.  The positive two-high gap
forbids two or more, so every actual residual assignment has exactly one high
slot.  Replacing that slot by its ray supremum enlarges the actual feasible
bank.

The complete terminal census is

```text
residual bodies / masks:       37 / 1,050,
zero-high hostile relaxations:       984,
one-high cases:                    1,150,
body-local low-label sets:           161,
coarse / exact closures:        1,075 / 75,
max-gap fallback / failures:        0 / 0,
minimum support slack:                  1.               (12)
```

Thus all `1,150` one-high cases close on complete projected cells.  The
terminal semantic, closure, and case-vector digests are

```text
cfa99f5f95e2e241da6cd2a91a199df21c738645cd7a0689d67dfe9926b26e17,
dc9fd29bc9c1c8525e12913294d312f38f58cdedbac4ffe013bfdddb8d20840c,
cf59598dc5966b961857567d81ea55e3bf0e494ea6af688e0b904ac8142a70c6. (13)
```

## 5. Consequence and next layer

The layer contains `219` rows, so the proved consequence is

```text
374025 - 219 = 373806,               z1 <= 221.           (14)
```

The next layer is occupied:

```text
z1=221: 90 rows = 83 wall + 7 order,
row SHA256 ab4f5b7fe3b5e1330d6e51dfb150527fa27cac56b755ed5e2ae2522f0f27ceb4. (15)
```

It is not screened by this theorem.

## 6. Evidence and scope

The checkpoint-independent canonical companion recomputes all `219` screen
rows and all `37` terminal bodies directly from the maintained atlas.  It
pins THM-3174's source, output, and semantic hashes, binds all `37` individual
residual mask banks, checks every raw Farkas certificate exactly, and verifies
the quotient/coupling census `(9)--(10)`.  Under MISTAKE-331/333, the screen
digest binds only `row[:19]`; the terminal closure digest omits the chosen
duplicate-gap maximizer witness.  Both are solver-witness-free.

This theorem concerns only the maintained projected `k=3` necessary atlas.
It does not classify physical covers outside that projection, recover the
forgotten order of the divisor-square repairs, alter `k<=1`, close the final
rung, or prove LRC(14).

Reproduction:

```text
python 04-computation/lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.py
python -O 04-computation/lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.py
```
