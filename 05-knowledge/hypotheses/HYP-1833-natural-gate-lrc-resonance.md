---
id: HYP-1833
status: EXPLORATORY
source: codex-2026-05-31-S377-natural-gates
related:
  - HYP-1817
  - HYP-1818
  - HYP-1819
  - HYP-1820
  - HYP-1821
  - HYP-1823
  - HYP-1832
  - THM-361
  - THM-362
---

# HYP-1833: Natural product-sum gates mark the fragile LRC quotient coordinates

## Statement

After scalar-ramp quotienting, the hardest one-coordinate Lonely Runner
micro-staircase defects occur at natural-number modes that are product-sum
critical targets.

Equivalently, the same natural modes where the additive and multiplicative
operation cospans collide,

```text
x_1 + ... + x_k = x_1 * ... * x_k,
```

also mark the quotient coordinates where a non-scalar LRC residue vector is
closest to blocking every micro-staircase cell.

## Evidence

S377 natural-gate scan compared all normalized one-coordinate defects for the
full `n=14` and `n=15` cell systems, then tagged a coordinate `i` when `i` is
a product-sum target mode.

For `n=14,k=13`, the product-sum target coordinates inside the quotient are

```text
4, 6, 8, 9, 10, 12.
```

The best collision target misses only `56` cells, while the best non-target
coordinate misses `154`.  The two best coordinates are:

```text
coord 6, residue 7:  missed 56
coord 12, residue 7: missed 84
```

Here `6` is the first visible distinct product-sum resonance
`1+2+3=1*2*3`, and `12` is the first richer target mode with three product-sum
cores through the small arity cutoff.

For `n=15,k=14`, the product-sum target coordinates are

```text
4, 6, 8, 9, 10, 12, 14.
```

The best collision target misses `120` cells, while the best non-target
coordinate misses `260`.  The best one-defects are tied order-3 stencils:

```text
coord 6,  residue 5 or 10: missed 120
coord 14, residue 5 or 10: missed 120
```

Thus the S377 torsion/CRT leak dichotomy is not separate from the natural-mode
graph.  It is concentrated on product-sum collision coordinates.

The same script tested `300` operation-critical speed probes for each frontier
case by inserting product-sum target modes into initial-segment speed sets and
by overloading divisor payloads.  All `600` tested sets had positive complement
gaps; no boundary-only or open-cover disproof candidate appeared.

## Interpretation

The old natural-number graph prompt asked for arrows into `z` when `x+y=z` and
for the sparser multiplication version.  THM-362 says the additive simple
shadow collapses to the transitive order, while multiplication remains the
divisor DAG.  HYP-1833 says the LRC quotient sees the missing labeled
structure: product-sum critical targets are exactly the coordinates where the
cell system is most fragile.

This gives a proof-search order:

1. Remove scalar ramps by THM-363/THM-364.
2. Prove the product-sum target coordinates leak through explicit torsion or
   CRT stencils.
3. Treat non-target coordinates afterwards; they appear to have much larger
   exposed-cell margins.

## Test Plan

1. Prove the coordinate-6 `n=14` eight-stencil lemma as the first visible
   product-sum resonance case.
2. Derive the `n=15` order-3 formulas for coordinates `6` and `14`.
3. Extend the one-defect target/non-target comparison to `n=16` and `n=18`.
4. Feed the product-sum coordinate priority into a branch-and-bound quotient
   certificate for HYP-1823.

## Sources

- `04-computation/lonely_runner_natural_gate_feedback_s377.py`
- `05-knowledge/results/lonely_runner_natural_gate_feedback_s377.out`
- `07-reflections/lonely-runner-natural-gates-s377.md`
- HYP-1820
- HYP-1821
- HYP-1823
- HYP-1832
- THM-361
- THM-362
