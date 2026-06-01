---
id: THM-376
name: lrc-first-even-ladder-gap-debt-conservation
status: PROVED (exact finite computation)
date: 2026-06-01
session: codex-2026-06-01-S500
depends_on:
  - THM-357
---

# THM-376: The n=14 and n=18 gate ladders conserve gap-debt product

## Statement

For `n=14`, let the row-parent, gate, and double-gate ladder rows be the speed
sets

```text
L_14(d) = {1} union {d q : 1 <= q < 14, q != 6}
```

for `d=7,14,28`.

For `n=18`, let

```text
L_18(d) = {1} union {d q : 1 <= q < 18, q != 8}
```

for `d=9,18,36`.

Let

```text
gap_ratio = max_gap / (1/n)
debt      = number of unprotected endpoint values.
```

Then:

```text
n=14:
  d=7   gap_ratio=5/924   debt=84    product=5/11
  d=14  gap_ratio=5/1848  debt=168   product=5/11
  d=28  gap_ratio=5/3696  debt=336   product=5/11

n=18:
  d=9   gap_ratio=1/176   debt=176   product=1
  d=18  gap_ratio=1/352   debt=352   product=1
  d=36  gap_ratio=1/704   debt=704   product=1
```

Thus, in these first-even ladder cascades, multiplying the gate depth halves
the visible Archimedean gap and doubles the exposed endpoint debt.

## Verification Record

The exact values are certified by

```text
04-computation/lrc_n14_n18_tournament_pingpong_s481.py
05-knowledge/results/lrc_n14_n18_tournament_pingpong_s481.out
```

using the exact interval-union code from `lonely_runner_residue_probe_s356.py`
and the exact endpoint-protection code from
`lonely_runner_endpoint_protection_s360.py`.

## Proof

For each of the six listed finite speed sets, the scripts compute:

1. the forbidden interval union on `R/Z` exactly over `Fraction`;
2. the largest complementary gap;
3. the finite endpoint set;
4. which endpoints are strictly protected by at least one speed.

The stored rational output gives the six `gap_ratio` and `debt` values listed
above.  Multiplying each pair gives the displayed products.

Because the speed sets are explicit, the interval endpoints are rational, and
the protection predicate is exact, this proves the stated finite certificate.

## Significance

This is a precise version of the "scalar improvement exports endpoint debt"
principle.  The quotient/gate ladders do not approach a counterexample by
erasing obstruction; they translate the same obstruction into deeper endpoint
layers while preserving a gap-debt invariant.

## Related

- HYP-1866
- HYP-1942
- HYP-1950
- THM-357
- `07-reflections/lrc-n14-n18-tournament-pingpong-s481.md`
