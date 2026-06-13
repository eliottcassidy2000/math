---
id: HYP-1969
status: OPEN
source: codex-2026-06-01-S505
related:
  - HYP-1967
  - HYP-1968
  - HYP-1942
  - HYP-1952
  - THM-373
  - THM-380
---

# HYP-1969: LRC hard ladders have an H-phase plateau after the first gate

## Statement

For the hard first-even LRC scale ladders, the half-turn phase tournament at
the selected lonely midpoint enters a stable Hamiltonian-path phase after the
first gate:

```text
scale = odd_core       -> high row-parent H phase;
scale = 2*odd_core     -> lower gate H phase;
scale = 4*odd_core,... -> same selected-midpoint H phase.
```

Meanwhile the anchored endpoint ledger keeps translating:

```text
gap halves, endpoint debt doubles, gap*debt stays constant.
```

Thus `H(T)` is a global loneliness entropy for the coarse phase clock, while
endpoint debt is the recursive state variable for the anchored LRC clock.

## Evidence

`04-computation/lrc_h_loneliness_metric_s505.py` computes exact Hamiltonian-path
counts for the half-turn phase tournament through `n=18`.

For the `n=14` hard ladder:

```text
scale 7:   H=22168229, H_ratio=2.0831, gap*debt=5/11
scale 14:  H=17826951, H_ratio=1.6752, gap*debt=5/11
scale 28:  H=17826951, H_ratio=1.6752, gap*debt=5/11
scale 56:  H=17826951, H_ratio=1.6752, gap*debt=5/11
```

For the `n=18` hard ladder:

```text
scale 9:   H=117137481061, H_ratio=2.3981, gap*debt=1
scale 18:  H=102405804217, H_ratio=2.0965, gap*debt=1
scale 36:  H=102405804217, H_ratio=2.0965, gap*debt=1
scale 72:  H=102405804217, H_ratio=2.0965, gap*debt=1
```

The corridor H profiles show small adjacent-cell variation, but the selected
midpoint phase state plateaus exactly in these audited rows.

## Predictions

1. For fixed hard skip in the first-even ladders, all further dyadic row
   doublings after the first gate preserve the selected-midpoint half-turn
   tournament up to isomorphism.
2. The endpoint-depth packet translates in the 2-adic direction while the
   selected H phase remains fixed.
3. H_ratio is a coarse danger filter: rows with tiny scalar gaps but plateaued
   H are endpoint-debt recursions, not new global phase obstructions.
4. A real counterexample candidate should break the plateau by combining
   endpoint-core survival with a new pressure/H phase event.

## Proof Program

The candidate formal lemma is:

```text
For ladder speeds {1} union {2^r*m*q : q != skip}
at the selected hard midpoint t_r,
the signs of all half-turn phase differences stabilize for r >= 1.
```

If true, the LRC recursion decomposes into:

```text
stable phase tournament
+ translated endpoint boundary packet
+ pressure DAG labels.
```

This would explain why H is useful as a loneliness metric without being the
whole LRC proof: it certifies that the global phase state has stopped changing,
so the remaining work is endpoint arithmetic.

## See Also

`07-reflections/lrc-h-loneliness-recursion-s505.md`,
`05-knowledge/results/lrc_h_loneliness_metric_s505.out`,
HYP-1967, HYP-1968, THM-373, THM-380.
