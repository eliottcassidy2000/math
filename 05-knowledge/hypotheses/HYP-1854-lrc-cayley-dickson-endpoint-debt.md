---
id: HYP-1854
status: EXPLORATORY
source: codex-2026-05-31-S387
related:
  - HYP-1838
  - HYP-1839
  - HYP-1840
  - HYP-1841
  - HYP-1853
---

# HYP-1854: LRC endpoint debt is a Cayley-Dickson-style doubling defect

## Statement

The Lonely Runner denominator tower admits a useful Cayley-Dickson-style
coordinate system:

```text
n = 2^r * odd_core.
```

The `2^r` factor is the doubling row, analogous to the Cayley-Dickson level.
The odd core is the torsion payload carried by that row.  At each gate/quotient
step, LRC behaves like a late Cayley-Dickson algebra: closing one older
freedom exposes a new defect layer.

In this analogy, endpoint debt is the LRC analogue of zero-divisor leakage.
An `n`-gate closes the unit skeleton, but the obstruction reappears in
descendant quotient layers.

## Evidence

`lrc_cayley_dickson_tower_s387.py` compares denominators `14..24` by
2-adic row, odd core, lpd-ladder pressure, endpoint debt layers, and
product-sum seed.

The scan separates the above-14 candidates:

```text
n=16 = 2^4       pure sedenion-row laboratory
n=18 = 2 * 3^2   complex row with square 3-torsion payload
n=20 = 4 * 5     quaternion row with 5-payload
n=21 = 3 * 7     real row with 3-by-7 transfer payload
n=24 = 8 * 3     octonion row stress test
```

For `n=18`, the exact lpd ladder

```text
(1, 9, 18, 27, 36, 45, 54, 63, 81, 90, 99, 108, 117, 126, 135, 144, 153)
```

has:

```text
lpd=9
skip=8
gap/th=0.005682
unprotected endpoints=176
first leak=11/162
endpoint-debt layer histogram=(9:48, 27:16, 99:16, 117:16, 144:64, 153:16)
```

The first exposed layer is exactly the largest-proper-divisor payload `9`.
That is the key structural reason `n=18` is attractive: it is not pure
doubling, but its first defect layer is still crisp.

## Predictions

1. For composite `n`, the first endpoint-debt layer of the lpd ladder is
   generally `lpd(n)`.
2. Pure 2-adic denominators, especially `n=16`, should give the cleanest
   "sedenion zero-divisor" endpoint-debt lemmas.
3. Mixed rows such as `n=18` should be the best place to see how an odd
   torsion payload rides on a doubling row.
4. Stress-test denominators such as `n=24` should not be attacked first; they
   are where the already-proved endpoint-debt mechanism should be validated.

## Sources

- `04-computation/lrc_cayley_dickson_tower_s387.py`
- `05-knowledge/results/lrc_cayley_dickson_tower_s387.out`
- HYP-1838
- HYP-1839
- HYP-1840
- HYP-1841
- HYP-1853
