---
id: THM-2837
title: "Divisor/Mordell reduction and polynomial-family exclusions for the ten smallest-open cubic Diophantine equations"
status: >
  PROVED (divisor-pair reformulation; Mordell-family reduction; linear-x
  impossibility; odd-split parity; linear-factor exclusion) + FINITE-EXACT
  (three exhaustive shape scans, boxes 12–14, zero families for all ten
  equations).  The equations' infinitude problems remain OPEN; the residual
  polynomial-family channels are deg x >= 4 even splits and the deg-3
  norm-form channel.
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath "Finiteness Problem for Diophantine Equations" = Grechuk
  catalogue smallest-open cubics)
depends_on: []
related: []
script: 04-computation/grechuk_cubic_dioph_family_exclusion_macmini_S171.py
output: 05-knowledge/results/grechuk_cubic_dioph_family_exclusion_macmini_S171.out
script_sha256: e045c472f68de265486422bb780d3b079ee17fcc361acbe5e43765ebaaade61a
output_sha256: 475dda22f90d5f9f2cdb9044704c4893a0399b2ce0130f2bc04c011bffa0a4a3
hash_basis: LF-normalized bytes
---

# THM-2837 — where solutions of z² + y²z + P(x) = 0 can and cannot come from

## Reformulation (proved)

For `P` one of the ten catalogued cubics (all irreducible over Q):
`z(z + y²) = -P(x)`, so with `A = -z`:

    A·B = P(x),   A + B = y²   (y² - 1 for the -z variant).

Integer solutions = divisors `d | P(x)` with `d + P(x)/d` a perfect square
(then `z = -d`).  Clearing denominators at fixed `d`: `V² = U³ + d³(d² + c)`
with `U = d·x` — a Mordell curve per divisor-parameter; infinitude of
solutions is equivalent to infinitely many `d` carrying an integral point of
the shape `(dx, d²y)`.  This is the clean structural home of the problem and
explains its difficulty (log-Calabi–Yau, cf. sums of three cubes).

## Exclusions

Polynomial families `(x(t), y(t), z(t))` are factorizations
`P(x(t)) = A·B` with `A + B = y²`:

1. **deg x = 1 impossible** (affine substitutions preserve irreducibility).
2. **Odd-degree splits die by parity** (`A+B = y²` has even degree unless
   leads cancel); for deg x = 3 only `(1,8)` (dead: rational root) and the
   non-generic `(3,6)` **norm-form channel** survive scrutiny.
3. **deg x = 2 exhaustively excluded in boxes**: `(2,4)` splits over six
   lead-pairs (box 14), constant divisors `A = d` (box 14), and the
   cancelling `(3,3)` split with `L ∈ {-1,-4,-9}` via the `y⁴ - 4P(x) = W²`
   criterion (box 12): **zero families for all ten equations**.
4. **Pell-shape remark**: the `z = -m²` two-parameter ansatz forces a
   rational root of `P` at the `m⁰` coefficient — impossible.

## Boundary

Boxes are finite; deg x >= 4 shapes and the deg-3 norm-form channel are
unsearched (their parameter counts exceed honest brute force); non-polynomial
(recurrence/Pell-orbit) solution production is untouched.  No claim about
the truth of infinitude — this theorem only maps the terrain.
