---
id: HYP-1952
status: OPEN
source: codex-2026-06-01-S493
related:
  - HYP-1866
  - HYP-1868
  - HYP-1905
  - HYP-1921
  - HYP-1930
  - HYP-1942
  - HYP-1950
---

# HYP-1952: n=18 LRC admits a square-core packet certificate

## Statement

The `n=18` Lonely Runner row should be attacked by a four-part certificate:

```text
local fan lemma
  -> global bridge-charge lemma
  -> dyadic-quotiented square-core packet lemma
  -> pressure-leaf peeling lemma.
```

The point is that the obvious local branch is already exhausted. After the
forced fan

```text
(1,5,7,9,11,13,17)
```

the two bridge choices `6` and `12` cover exactly the same six residual
owner-18 endpoints. The bridge is therefore not a local cover variable. It
is a quotient variable whose sign can only be detected by bringing in global
owner rows and product-depth endpoint packets.

## Evidence

`04-computation/lrc_n18_abstract_certificate_s493.py` finds that the residual
owner-18 endpoints after the forced fan are:

```text
53/324, 55/324, 161/324, 163/324, 269/324, 271/324
```

Both bridges cover all six. The residual symmetric difference is `0`.

The global owner-row shadow does distinguish the bridge sign:

```text
owner endpoints cover_by_6 cover_by_12 delta(6-12) winner
    7        14          6           8          -2     12
    8        16         10           8           2      6
    9        18         12           8           4      6
   10        20         20          18           2      6
   11        22          8          10          -2     12
   12        24         20          14           6      6
   16        32         16          20          -4     12
```

So a dual certificate should weight owner rows `7,8,9,10,11,12,16` so that
both signs pay positive cost unless another proof mechanism fires.

The row-parent/gate/double-gate ladder does not create three separate
obstructions. It is one square-core packet translated in the 2-adic direction:

```text
row-parent 9-ladder skip 8      gap/th=1/176 debt=176 product=1
18-gate ladder skip 8           gap/th=1/352 debt=352 product=1
36-double-gate ladder skip 8    gap/th=1/704 debt=704 product=1
```

After normalizing the scale `9*2^r` by subtracting `r` from the 2-depth and
dividing counts by `2^r`, every row has the same packet:

```text
96*(0,2) + 16*(0,3) + 64*(4,2)
```

This is the `n=18` square-core: a mixed `2 x 3^2` product-building residue
whose apparent gap shrinkage is exactly balanced by endpoint debt.

S493 also reuses the S492 pressure lifts on three `n=18` rows:

```text
n=18 lpd ladder:         max pressure SCC 1, max pressure triangle 0
n=18 gate ladder:        max pressure SCC 1, max pressure triangle 0
n=18 single-gate repair: max pressure SCC 1, max pressure triangle 0
```

This remains true for `k1`, `k2`, and threshold-deficit pressure. Current
rows are therefore pressure-peelable rather than disproof-like pressure cores.

## Predictions

1. A proof of the `n=18` row should not branch deeply on `6` versus `12` at
   the local owner-18 gate; that branch is a local quotient illusion.
2. The finite dual target is a row-weight vector on owner rows
   `7,8,9,10,11,12,16` plus square-core packet rows. The vector should charge
   both bridge signs unless a pressure SCC appears.
3. Nonmultiple repairs can have smaller archimedean gap, but they carry
   endpoint/product debt. Gap alone is not a monotone search objective.
4. If an `n=18` counterexample-like branch exists, it must survive two
   independent peelings at once: endpoint-private peeling and mobile
   pressure-leaf peeling.

## Next Tests

- Solve a small linear program for nonnegative row weights on
  `7,8,9,10,11,12,16` that makes the two bridge signs simultaneously
  expensive.
- Add the normalized packet variables
  `96*(0,2) + 16*(0,3) + 64*(4,2)` to the same linear certificate.
- Turn the pressure-leaf observation into an exact peeling algorithm on
  endpoint rows, not just sampled time slices.
- Search bounded `n=18` perturbations for the first row with
  `k2_largest_scc > 1` or `deficit_largest_scc > 1`; if none appears, use the
  absence as a certificate side condition.

## See Also

- `04-computation/lrc_n18_abstract_certificate_s493.py`
- `05-knowledge/results/lrc_n18_abstract_certificate_s493.out`
- `07-reflections/lrc-n18-abstract-certificate-s493.md`
- HYP-1942
- HYP-1950
- HYP-1930
- HYP-1921
