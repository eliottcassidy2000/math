---
id: HYP-1832
status: EXPLORATORY
source: codex-2026-05-31-S377
related:
  - HYP-1818
  - HYP-1819
  - HYP-1823
  - HYP-1824
  - HYP-1825
  - HYP-1833
  - HYP-1837
  - THM-363
  - THM-364
---

# HYP-1832: LRC torsion/CRT leak dichotomy

## Statement

After scalar-ramp quotienting, composite Lonely Runner micro-staircase
obstructions leak through small torsion or CRT subgroup layers before they can
become full blockers or open-cover counterexamples.

More concretely:

- For `n=14`, the hard quotient boundary is isolated `2`-torsion.
- For `n=15`, the replacement boundary is order-3 leakage inside `Z/15Z`.
- Gated speed-set disproof attempts inherit the same leaks as positive
  complement gaps.

## Evidence

S377 reused the exact S367 full-cell bitset system.

For `n=14,k=13`, the complete normalized `2`-torsion cube has a unique global
minimum:

```text
vector=(0,0,0,0,0,7,0,0,0,0,0,0,0)
missed=56
```

The missed cells have the same S367/S371 ledger:

```text
7 odd shifts * 8 alpha stencils
gcd shift split: (1,48), (7,8)
all minimum margins: 1
```

Adding any second defect to that coordinate-6 half-turn exposes more cells; the
best such perturbation misses `126` cells.

The existing S372 cycle support-4 scan checked `14,137,695` normalized
support-4 vectors and found best missed `182`, so support `<=4`, full
`2`-torsion, and local search all point back to the same coordinate-6
half-turn.

For `n=15,k=14`, exact small-support scans found:

```text
support 1, all residues: best missed 120
support 2, all residues: best missed 220
support 3, order-3 subgroup {5,10}: best missed 280
support 3, order-5 subgroup {3,6,9,12}: best missed 294
```

The best one-coordinate vectors occur at coordinates `6` and `14` with residue
`5` or `10`.  The displayed extremal has shift histogram on the ten shifts not
divisible by `3`, with gcd split `(1,96),(5,24)`, and all missed cells have
positive margin.

The disproof route also leaked.  Mutated gated interval-cover candidates kept
positive complement gaps:

```text
k=13,n=14 one-gate best max gap: 5/1848
k=14,n=15 one-gate best max gap: 1/495
k=13,n=14 two-gate best max gap: 3/784
k=14,n=15 two-gate best max gap: 37/9405
```

The natural-gate follow-up `lonely_runner_natural_gate_feedback_s377.py`
connects the torsion/CRT leaks to product-sum target coordinates.  For
`n=14`, the best product-sum target coordinate misses `56` cells while the
best non-target coordinate misses `154`; for `n=15`, the corresponding numbers
are `120` and `260`.  The `n=14` extremal coordinate `6` is the first visible
distinct product-sum resonance, and the `n=15` tied coordinates `6` and `14`
are both product-sum targets.

## Why It Matters

HYP-1823 says the fourteen-runner proof should prove that zero is the only
full blocker after quotienting by scalar ramps.  HYP-1832 adds a proof-search
grammar for the next layers: do not search arbitrary quotient vectors first.
Search for the torsion or CRT subgroup through which a defect must leak.

This also gives a cleaner role for disproof attempts.  A genuine LRC
counterexample must defeat not just the scalar-ramp skeleton, but also every
torsion/CRT leak that appears when composite gates are inserted.

## Test Plan

1. Prove the `n=14` coordinate-6 half-turn ledger as an eight-stencil interval
   lemma.
2. Derive a one-defect formula for `n=15` explaining why coordinates `6` and
   `14` with residues `5/10` are extremal.
3. Extend exact subgroup scans to `n=16` and `n=18` to see whether the leak is
   controlled by prime-power torsion or by the squarefree CRT split.
4. Turn endpoint-protection disproof search inside out: prescribe a protection
   cycle first, then solve for speed sets, and check whether the cycle is
   killed by a torsion/CRT leak.

## Sources

- `04-computation/lonely_runner_torsion_crt_feedback_s377.py`
- `04-computation/lonely_runner_cycle_s372.py`
- `05-knowledge/results/lonely_runner_torsion_crt_feedback_s377.out`
- `05-knowledge/results/lonely_runner_cycle_s372.out`
- `07-reflections/lonely-runner-torsion-crt-feedback-s377.md`
- HYP-1818
- HYP-1819
- HYP-1823
- HYP-1824
- HYP-1825
- HYP-1833
- HYP-1837
