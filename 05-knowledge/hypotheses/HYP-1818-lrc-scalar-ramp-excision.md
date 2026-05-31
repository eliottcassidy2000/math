---
id: HYP-1818
status: EXPLORATORY
source: codex-2026-05-31-S364
related:
  - THM-358
  - THM-360
  - HYP-1813
  - HYP-1816
  - HYP-1817
---

# HYP-1818: LRC scalar-ramp excision before composite micro-staircases

## Statement

For composite `n=k+1`, any micro-staircase tight-lift lemma must first split
off the affine scalar-ramp family

```text
v_i = m i mod n,  i=1,...,n-1.
```

These vectors are not accidental obstructions.  They are the residue shadow of
the Dirichlet equality example `(1,2,...,n-1)`, after shifting and scaling the
time parameter.  Unit multipliers are reindexed initial segments; nonunit
multipliers are quotient/descent cases that should be handled by the
divisibility machinery rather than by a generic cell lemma.

After this scalar-ramp family is excised, the `n=14` and `n=15` full-cell
searches should have no remaining full blockers.  The proof task then becomes
certifying the missed cells of the best non-scalar near-blockers.

## Evidence

`lonely_runner_feedback_loop_s364.py` builds the full discontinuity arrangement
of

```text
floor(n * {i alpha}),  i=1,...,n-1.
```

For `n=14`, the arrangement has `812` patterns and `14*812=11368`
candidate `(s, cell)` pairs.  Every scalar multiplier `m mod 14` blocks all
`11368` candidates.  The gcd histogram among these full blockers is:

```text
gcd(m,14): 1 -> 6, 2 -> 6, 7 -> 1, 14 -> 1.
```

For `n=15`, the arrangement has `960` patterns and `15*960=14400`
candidate pairs.  Every scalar multiplier `m mod 15` blocks all `14400`
candidates.  The gcd histogram is:

```text
gcd(m,15): 1 -> 8, 3 -> 4, 5 -> 2, 15 -> 1.
```

Once scalar ramps are excluded, deterministic local search found no full
blocker:

```text
n=14 best non-scalar: 11312/11368 covered, 56 missed.
n=15 best non-scalar: 14280/14400 covered, 120 missed.
```

The same run tried deterministic and random gated speed-set constructions for
the `14`- and `15`-runner cases.  It found no open-cover candidate; all sampled
hard cases retained either a boundary witness or a positive open gap.

## Interpretation

The scalar ramp identity is the finite cell version of THM-358.  If
`v_i=m i`, then

```text
s v_i + floor(n * {i alpha})
```

is the floor-vector for the initial segment at a shifted time, up to boundary
carry behavior.  The initial segment is tight but only boundary-tight, so every
open cell is blocked.

This explains why the first S363 micro-staircase statement was too optimistic:
it cannot ask every residue vector to have a cell witness.  It must first
remove the Dirichlet equality line, then use quotient descent for nonunit
ramifications, and only then classify genuinely non-scalar vectors.

## Test Plan

1. Prove the scalar-ramp blocking identity for all `n`.
2. Separate unit scalar ramps from nonunit ramps in the lifted prime-grid
   variables and map the nonunit cases to the endpoint/divisibility descent.
3. Enumerate all non-scalar full blockers for `n=14` by exact SAT/backtracking,
   not just local search.
4. For the best `n=14` near-blocker, list the `56` missed cells and search for
   a short mixed-threshold or endpoint-pressure certificate.
5. Repeat the same missed-cell classification for the `n=15` best near-blocker
   to see which features transfer.

## Sources

- `04-computation/lonely_runner_feedback_loop_s364.py`.
- `05-knowledge/results/lonely_runner_feedback_loop_s364.out`.
- HYP-1817 and `04-computation/lonely_runner_k13_microstaircase_s363.py`.
- Sungkawichai and Trakulthongchai, `arXiv:2604.23906`.
- Jensen, `arXiv:2605.27941`.
