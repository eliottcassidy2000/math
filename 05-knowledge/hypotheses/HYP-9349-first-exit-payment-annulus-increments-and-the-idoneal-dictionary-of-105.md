---
id: HYP-9349
title: "Collatz ledgers after THM-4735: the first-exit payment conjecture (D_w < n - h, equivalently P_w < (1 - lambda_w)(n-1), for every first exit), the annulus increment law (E(64N,D) - E(N,D) <= c 153N at every finite arrival horizon with c about 0.6-0.7 on the four certified scales, cross term always negative), and the idoneal dictionary of 105 (the mod-15/mod-35 guards are the genus characters of the one-class-per-genus discriminant -420)"
status: OPEN (pre-audit); finite evidence and first untested cases stated for each item
source: mac-mini-2026-10-20-debtcones
depends_on:
  - 01-canon/theorems/THM-4735-valuation-dependent-debt-cones-optimal-terminal-splits-and-the-debt-two-obstruction-to-universal-cone-coverage.md
  - 05-knowledge/results/native_cancellation_105_20261018g.md
  - 01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md
related:
  - 05-knowledge/results/charge_palettes_boundary_20261020g/charge.md
  - 05-knowledge/results/guarded_children_105_20261018g.md
  - 05-knowledge/results/root_eisenstein_first_arrival_20261018j.md
scripts:
  - 04-computation/experiments/annulus_horizon_20261020/annulus_horizon_increment.py
  - 04-computation/experiments/debt_cones_20261020/genus_420.py
---

# HYP-9349 -- three Collatz ledgers

## (a) Every first exit pays

**Claim.** For every odd `n > 1` whose orbit first drops below `n` after the odd word `w`,
`D_w < n - h`; equivalently `P_w < (1 - lambda_w)(n - 1)`: the descent is never produced by
the debt alone. Evidence: all 500,000 first exits below `1,000,001` (THM-4735 (D)); the
identity of the two forms checked exactly. THM-4735 (C) shows that the uniform certificate
`D < 2` fails for 42 of them, so a proof must use the actual drop. Route: the final step has
`a >= 2`; writing the exit as `lambda_w z + P_w - D_w < z` and the prefix non-exit as
`lambda_u z + P_u - D_u >= z`, the claim is `P_w < (1 - lambda_w) z`; the prefix inequality
gives `(lambda_u - 1) z >= D_u - P_u`, and one needs to compare `P_w = mu P_u + [a = 1]`
with `(1 - mu lambda_u) z`. First untested case: sources above `10^6`; a counterexample would
be a first exit whose prefix is barely non-descending and whose last valuation is large.

## (b) The annulus increment at every horizon

With the inherited checker's own source rows (first-arrival routes to the 627-point frontier,
signed Eisenstein charges `v_s = (-1)^{t_s} A_N(s)`, terminal vectors `W_u`, energy
`E = sum_u Q(W_u)`), define `E(N, D)` using only sources of route depth `<= D`. The exact
scale identity `S(64N) - S(N) = 153N` is reproduced on all four scales (new sources
387,072 / 566,716 / 1,039,689 / 1,548,288 at `N = 4096, 5997, 11002, 16384`). **Observed law:**

```text
N      max_D  [E(64N,D) - E(N,D)] / 153N    at D      full arrival     sum cross(W_u,H_u) at full arrival
4096   0.6782                               149       0.6782           -8,039
5997   0.6502                                76       0.6468           -24,200
11002  0.6153                                76       0.6103           -25,582 (D = 100 shown)
16384  0.5944                                76       0.5933           -45,186 (D = 100 shown)
```

The increment is monotone in `D` up to about `D = 75` and then flat; the cross term
`sum_u cross(W_u, H_u)` is **negative at every horizon** on all four scales (the new annulus
sources interfere destructively with the old terminal vectors), so the increment is below the
new sources' own energy `sum_u Q(H_u)`, itself below `0.7 x 153N`. **Claim.** For every
`N >= 2^12` and every horizon `D`, `E(64N, D) - E(N, D) <= c 153N` with a constant `c < 1`
(the data suggest `c` decreasing with `N`: `0.68, 0.65, 0.62, 0.59`). This is exactly the
increment bound the packet named as sufficient to transport `E(N) <= c S(N)` along the
64-scale ray; it is not proved, and finitely many base scales do not cover all integers.
First test: the next scale `N = 2^15` (routes for `2^23` sources), and a sign proof for the
cross term from the three-shell structure of the annulus (`LLL, LLU, LL, LU, L, U`).

## (c) The idoneal dictionary of 105

`105 = 3 x 5 x 7` is idoneal (one of Euler's 65; so are `15, 21, 165, 385, 1365`, not `35`,
`63`, `153` or `420`). Discriminant `-420` has `h = 8` classes, one per genus; the eight
reduced forms `x^2 + 105 y^2, 2x^2 + 2xy + 53y^2, 3x^2 + 35y^2, 5x^2 + 21y^2, 6x^2 + 6xy + 19y^2,
7x^2 + 15y^2, 10x^2 + 10xy + 13y^2, 11x^2 + 8xy + 11y^2` are separated exactly by the genus
characters `(n/3), (n/5), (n/7)` and the character mod 4, and each represents exactly six
residue classes mod `420` = six `(mod 15, mod 35)` register pairs (genus_420.py). So the
lane's two compatible registers (which agree mod 5 and encode `s mod 105` jointly) together
with the mod-4 phase are precisely the genus characters of the idoneal discriminant `-420`:
"guarded children mod 15 / mod 35" is genus theory for `105`, and idoneality is the statement
that these guards determine the representing form, i.e. that the register pair is a complete
invariant of the class. **Question** (not a claim): does the ROOT first-arrival terminal of a
source depend on its class under the forms of discriminant `-420` beyond its residue mod 105,
and do the 19 idoneal planes of THM-4566 (`m = 105` among them, `chi >= 3`) supply a colouring
structure for the terminals? No transport of energy statements across this dictionary is
asserted; the dictionary itself is exact.
