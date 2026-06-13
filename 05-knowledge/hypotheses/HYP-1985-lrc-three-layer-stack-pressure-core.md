---
id: HYP-1985
status: OPEN
source: codex-2026-06-01-S514
related:
  - HYP-1975
  - HYP-1976
  - HYP-1977
  - HYP-1978
  - HYP-1979
  - HYP-1980
  - HYP-1981
  - HYP-1982
  - HYP-1983
  - HYP-1984
  - THM-377
  - THM-380
---

# HYP-1985: The LRC stack obstruction needs an owner-compatible pressure core

## Statement

The S513 add/multiply stack is a useful way to rank LRC hard rows, but a
counterexample-shaped obstruction cannot be certified by coarse row-time
pressure SCCs alone.  The relevant pressure coordinate must be the labelled
owner-compatible endpoint core from the THM-380 certificate trilemma.

Equivalently, a bad LRC row should not merely have:

```text
low additive representation abundance
+ high endpoint debt
+ product-sum/gate pressure
+ marked A000568 projection trouble
+ current pair-cell danger.
```

It must also carry a nonpeeling owner-pressure core after endpoint labels are
kept.  If the coarse pressure tournament peels, then the stack obstruction is
not yet counterexample-shaped, even if the arithmetic and danger coordinates
look bad.

## Evidence

`lrc_three_layer_stack_audit_s514.py` audits 11 selected hard-row families:

```text
n14 initial / row-parent / gate / double-gate
n16 initial / row-parent / gate
n18 initial / row-parent / gate / double-gate
```

It builds Tournament Analysis row tournaments with declared runner, operation,
and denominator observables.  The full stack is not a scalar ledger:

```text
full_stack: H=87, c3=10, SCCs=(8,1,1,1)
top rows: n14 double-gate, n16 gate, n18 double-gate, n18 gate
```

The layer comparison is informative:

```text
runner_dynamic:       H=237, c3=12, SCCs=(8,1,1,1)
operation_weighted:   H=1,   c3=0,  transitive
denominator_static:   H=1,   c3=0,  transitive
```

Edge flips show that the full stack is mostly operation-weighted danger plus a
runner correction:

```text
runner vs denominator:      0.45
operation vs denominator:   0.53
operation vs full stack:    0.07
runner vs full stack:       0.25
```

The row-time state tournament over 38 selected states has:

```text
c3=2, SCCs=(4,1,1,...)
```

and its top states are the expected unsafe gate/double-gate wall states, led
by:

```text
n18 double-gate half-unit
n14 double-gate half-unit
n16 gate unit
n18 gate unit
```

The counterexample-shaped conjunction audit is the decisive refinement.  The
arithmetic coordinates are easy to trigger: n14 and n16 rows score `4/5`,
n18 rows score `3/5` on the median-threshold flags

```text
additive scarcity, endpoint debt, product-sum collisions, A000568 survival,
pressure.
```

But no sampled row triggers the pressure flag:

```text
pressure(nontrivial SCC): absent in every sampled row.
```

So the S513 conjectural obstruction currently fails at the pressure-peeling
coordinate, not at the add/multiply/A000568 coordinates.

## Predictions

1. Replacing the coarse row-time pressure SCC by THM-380 owner-compatible
   endpoint cores should separate true proof pressure from harmless dynamic
   danger.
2. Full-stack row tournaments should continue to rank gate and double-gate
   rows above row-parent rows, because operation-weighted danger tracks the
   current pair-cell deficit.
3. Static denominator labels should remain transitive ledgers on moderate
   selected denominator sets; nontrivial shape should arise from conflicts
   with dynamic runner and endpoint labels.
4. A serious counterexample candidate must satisfy all three layers and retain
   a labelled nonpeeling endpoint-pressure core.  Arithmetic scarcity without
   that core should be treated as hard-row texture, not as obstruction.

## Proof Program

1. Reuse the S514 row/state stack features, but replace
   `pressure_SCC` with THM-380 owner-compatible endpoint-core data.
2. For each selected row, compute the endpoint-owner core after applying
   current pair-cell danger, operation labels, and marked A000568 source data.
3. Build the same row tournament and compare edge flips against the S514
   runner, operation, and denominator layers.
4. Try to prove that the owner-compatible pressure core always peels in the
   first-even hard rows `n=14,16,18`.

## Sources

- `04-computation/lrc_three_layer_stack_audit_s514.py`
- `05-knowledge/results/lrc_three_layer_stack_audit_s514.out`
- `07-reflections/lrc-three-layer-stack-pressure-core-s514.md`
- HYP-1984
- THM-380
