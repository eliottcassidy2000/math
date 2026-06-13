---
id: HYP-1818
status: PARTIALLY_PROVED
source: codex-2026-05-31-S364/codex-2026-05-31-S371
related:
  - THM-363
  - THM-364
  - THM-358
  - THM-360
  - HYP-1813
  - HYP-1816
  - HYP-1817
  - HYP-1823
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

THM-363 proves scalar-gauge reindexing, and THM-364 proves the scalar-ramp
blocking identity directly for all `n`, not only for the computed `n=14` and
`n=15` cases.  On every open micro-staircase cell, the shifted scalar-ramp
residues are exactly the initial-segment floor-vector at the shifted time
`alpha+s*m/n`; THM-358 then forces a blocked coordinate.

Once scalar ramps are excluded, deterministic local search found no full
blocker:

```text
n=14 best non-scalar: 11312/11368 covered, 56 missed.
n=15 best non-scalar: 14280/14400 covered, 120 missed.
```

The S371 follow-up reconstructs all `812` representative `n=14` cells and
checks the scalar midpoint identity with zero failures.  It also shows that the
S364 best `n=14` non-scalar vector is not generic:

```text
(8,2,10,4,12,13,0,8,2,10,4,12,6)
```

is the scalar ramp `m=8` with one coordinate changed, `v_6: 6 -> 13`.  Its
`56` missed cells are exactly the cells uniquely blocked by that scalar
coordinate.  An exact scan of all one- and two-coordinate scalar-neighborhood
deformations found no non-scalar full blocker.

The same run tried deterministic and random gated speed-set constructions for
the `14`- and `15`-runner cases.  It found no open-cover candidate; all sampled
hard cases retained either a boundary witness or a positive open gap.

S367 sharpens the excision into an actual quotient.  Adding `m*i` to a residue
vector preserves raw cell coverage by reindexing alpha cells; after normalizing
`v_1=0`, every scalar ramp becomes the zero vector.  In the exact `n=14`
system, the best nonzero quotient search result is

```text
(0,0,0,0,0,7,0,0,0,0,0,0,0),
```

with `56` missed cells.  Exact scans through support `3` and the full
normalized `2`-torsion cube found no full blocker; in that cube, the above
coordinate-6 half-turn is the unique best near-blocker.

## Interpretation

The scalar ramp identity is the finite cell version of THM-358.  If
`v_i=m i`, then

```text
s v_i + floor(n * {i alpha})
```

is the floor-vector for the initial segment at a shifted time, up to boundary
carry behavior.  THM-363 removes the ambiguity: on open cells this is the exact
shifted initial-segment floor-vector.  The initial segment is tight but only
boundary-tight, so every open cell is blocked.

This explains why the first S363 micro-staircase statement was too optimistic:
it cannot ask every residue vector to have a cell witness.  It must first
remove the Dirichlet equality line, then use quotient descent for nonunit
ramifications, and only then classify genuinely non-scalar vectors.

## Test Plan

1. Prove the scalar-ramp blocking identity and scalar-gauge reindexing lemma
   for all `n`.  Done in THM-363 and THM-364.
2. Separate unit scalar ramps from nonunit ramps in the lifted prime-grid
   variables and map the nonunit cases to the endpoint/divisibility descent.
3. Enumerate all non-scalar full blockers for `n=14` by exact SAT/backtracking,
   not just local search.  S371 has completed the radius-1 and radius-2
   scalar-neighborhood scan; a global search remains open.
4. For the best `n=14` near-blocker, list the `56` missed cells and search for
   a short mixed-threshold or endpoint-pressure certificate.  S371 lists the
   cells and identifies them as unique scalar-protection cells.
5. Repeat the same missed-cell classification for the `n=15` best near-blocker
   to see which features transfer.

## Sources

- `04-computation/lonely_runner_feedback_loop_s364.py`.
- `05-knowledge/results/lonely_runner_feedback_loop_s364.out`.
- `04-computation/lonely_runner_k13_scalar_gauge_s367.py`.
- `05-knowledge/results/lonely_runner_k13_scalar_gauge_s367.out`.
- THM-363.
- THM-364.
- HYP-1823.
- `04-computation/lonely_runner_scalar_excision_s371.py`.
- `05-knowledge/results/lonely_runner_scalar_excision_s371.out`.
- HYP-1817 and `04-computation/lonely_runner_k13_microstaircase_s363.py`.
- Sungkawichai and Trakulthongchai, `arXiv:2604.23906`.
- Jensen, `arXiv:2605.27941`.
