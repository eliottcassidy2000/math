---
id: HYP-1828
status: EXPLORATORY
source: codex-2026-05-31-S373
related:
  - THM-357
  - THM-358
  - THM-360
  - HYP-1816
  - HYP-1819
  - HYP-1823
  - HYP-1829
---

# HYP-1828: Fourteen-runner disproof searches should solve protection cycles first

## Statement

For the fourteen-runner reduced problem (`k=13`, threshold `1/14`), a
speed-first counterexample search is structurally backwards.  The necessary
`14`-gate can protect the unit boundary layer, but it repeatedly creates
higher-denominator leaks.  A plausible disproof architecture should first
construct a finite directed endpoint-protection cycle involving the unit layer
and at least one higher quotient layer, then solve for speeds realizing that
cycle.

Equivalently: search over endpoint-protection incidence patterns before
searching over speed sets.

## Evidence

S373 tried five speed-first construction grammars:

- single `14`-gate replacements of the initial segment;
- multi-gate overloads;
- CRT residue lifts, including the S367 half-turn near-blocker residues;
- anchored greedy set-cover pressure on a `12012`-point grid;
- mandatory-gate mutation search.

No exact open-cover candidate appeared among the final exact audits:

```text
open_cover_candidates=0
audited_boundary_only=1
audited_positive_gap=69
```

The best exact positive-gap near-miss was the anchored seven-ladder

```text
(1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)
```

with

```text
forbidden_length = 142/143
max_gap = 5/12936
gap/thresh = 0.005411
unprotected endpoints = 84
first_unprotected = 9/98
```

This is much closer in positive-gap size than simple gate insertion, but it is
not endpoint-protected.  It overcovers the `7`-quotient while exposing many
endpoint residues.

The best single-gate replacement,

```text
(1,2,3,4,5,7,8,9,10,11,12,13,224),
```

has only `12` unprotected endpoints, but its max gap is larger:

```text
max_gap = 11/9408
gap/thresh = 0.016369
first_unprotected = 29/182
```

Thus the two pressures split: small max gaps come from quotient ladders, while
small endpoint exposure comes from near-initial-segment gate replacements.

## Interpretation

The initial segment `(1,...,13)` is already full-measure but boundary-only,
with the six unit witnesses `a/14`.  THM-360 says a primitive open-cover
counterexample must contain a speed divisible by `14` to protect those unit
endpoints.  S373 confirms the catch: adding the gate protects the obvious
layer and reopens leaks at denominators such as `182`, `728`, `9408`,
`12936`, and `35315280`.

This suggests that a genuine disproof, if it exists, is not a random dense
cover.  It must be an engineered endpoint-protection cycle whose protecting
speeds do not merely push the leak to the next quotient layer.

## Next Tests

1. Build the finite endpoint graph on selected quotient layers, starting with
   `a/14`, `29/182`, `15/182`, and `9/98`.
2. Search for small directed protection cycles in that graph before assigning
   speeds.
3. Realize candidate cycles by solving the integer protection inequalities
   from THM-357/THM-360.
4. Use the seven-ladder near-miss as a second model: preserve its tiny gap
   while trying to reduce its `84` exposed endpoints.
5. Run exact one-swap and two-swap neighborhoods of the seven-ladder with a
   cached endpoint classifier rather than raw Fraction recomputation.

## Sources

- `04-computation/lonely_runner_14_disproof_hunt_s373.py`
- `05-knowledge/results/lonely_runner_14_disproof_hunt_s373.out`
- `07-reflections/lonely-runner-fourteen-disproof-hunt-s373.md`
- HYP-1816
- HYP-1819
- HYP-1823
- HYP-1829
