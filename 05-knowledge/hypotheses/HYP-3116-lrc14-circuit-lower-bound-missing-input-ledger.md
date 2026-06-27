---
id: HYP-3116
title: LRC14 circuit lower-bound and missing-input ledger
status: RESERVED / active synthesis; not a proof
source: codex-2026-06-27-S266
tangent: T1191
technique: LTI-252
tournament_technique: LTT-150
related:
  - HYP-3115
  - HYP-3114
  - HYP-3113
  - HYP-3112
  - HYP-3111
  - HYP-3108
  - HYP-3107
  - HYP-3098
  - HYP-3083
  - HYP-3074
  - HYP-3054
  - HYP-2997
  - HYP-2991
  - HYP-2963
  - OPEN-Q-108
---

# HYP-3116: LRC14 Circuit Lower-Bound And Missing-Input Ledger

## Claim

Circuit complexity should be used in the LRC14 proof search as a
missing-input discipline, not as a metaphor for hardness.  A proposed shortcut
is proof-facing only after it names the Boolean/proof circuit, its input basis,
which inputs are essential, which certificate minterms close the route, and
which deleted inputs are reconstructible or paid for by sidecars.

The active synthesis connects two older repo threads:

1. **Data circuits for tournament invariants.**  The staircase/tiling Walsh
   work treats Hamiltonian-path counts as a low-degree Boolean/carry circuit
   over fixed-path tiles.
2. **Proof circuits for LRC14.**  HYP-3111/HYP-3115 treat the current frontier
   as a monotone certificate DAG whose outputs are finite-address,
   observer-gluing, root/ear, lattice, and algebraic fold exits.

The live conjectural bridge is that failed scalar proof routes are exactly
low-depth quotient circuits with a nonzero missing-input vector.

## Evidence To Gather

- Mine prior Boolean-circuit, Walsh, cocycle, finite-address, observer-gluing,
  Savitch, and route-center work for reusable gate types.
- Build a small executable ledger of known gates, shortcut attempts, essential
  inputs, and missing-input vectors.
- Run Tournament Analysis on proof-gate families rather than runners, Boolean
  gates, or scalar metrics.
- Use the result to update OPEN-Q-108: any new LRC14 shortcut must either
  provide a certificate minterm or state the first missing input.

## Guardrail

HYP-3115's one-literal finite-bank rule `apex7_error <= 5` is a signal, not a
uniform proof.  This hypothesis exists to prevent that kind of finite fitted
classifier from masquerading as a theorem route.
