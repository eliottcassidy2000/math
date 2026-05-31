---
id: HYP-1817
status: EXPLORATORY
source: codex-2026-05-31-S363
related:
  - THM-363
  - THM-364
  - THM-357
  - THM-358
  - THM-360
  - HYP-1813
  - HYP-1818
  - HYP-1823
---

# HYP-1817: The fourteen-runner case needs a micro-staircase lift

## Statement

For the next open finite Lonely Runner case, `k=13` and `n=k+1=14`, the
right replacement for the prime-modulus tight-class lemma is a
micro-staircase lemma.

Concretely, let

```text
R_alpha(i) = floor(14 * {i alpha}),  i=1,...,13.
```

The proof should not only use the coarse cells `alpha = r/14`.  It should use
the finer cells in the discontinuity arrangement

```text
alpha = a/(14 i),  1 <= i <= 13.
```

For every lifted tight-class residue vector `v in (Z/14Z)^13` that is not
already a scalar-ramp/Dirichlet-equality case and is not handled by the gcd
condition, there should be some `s mod 14` and some cell of `R_alpha` such that

```text
s v_i + R_alpha(i) notin {0,13} mod 14
```

for all coordinates `i`.  If this cell has width `w`, then every prime grid
`1/p Z` with `p > 1/w` contains a usable `r/p` in that cell.

## Evidence

The exact `r/14` analogue of the Sungkawichai-Trakulthongchai prime lemma
fails.  The vector

```text
(7,4,9,6,7,8,5,0,1,12,13,12,7) mod 14
```

blocks all `14*14=196` pairs `(s,r)` with `r/14`.

But it is resolved by the micro-staircase cell

```text
R_alpha = (0,0,0,1,1,1,2,2,2,3,3,3,4)
```

with

```text
alpha in [2/91, 1/42),  width = 1/546.
```

At `s=1`, this sends the vector to

```text
(7,4,9,7,8,9,7,2,3,1,2,1,11),
```

which avoids both forbidden residues `0` and `13`.

The S363 script verifies the same obstruction is resolved on many actual
prime grids, including `p=181,191,199,211,251,307,401,701`.

## Relation To The Current Frontier

The public frontier is `k<=12` in the reduced Wills formulation, i.e. up to
thirteen total runners.  The next case is `k=13`, fourteen total runners.

The recent proof for `k=10` and `k=12` exploits that `k+1` is prime to remove
the tight class `(1,...,k)` analytically.  The case `k=13` has `k+1=14`, so the
prime-field polynomial proof cannot be copied directly.  HYP-1817 says the
correct composite replacement is not a new field argument but a finite
cell-arrangement argument for the floor-vector staircase.

## Computational Note

Running the public `vzsky/13-lonely-runners` verifier locally with its
experimental `LrcVerifier<13>` path on primes `17..101` found:

```text
all final lifted sets empty after Squeeze<2>, Squeeze<3>
```

For every tested prime, once `I(13,p,1)` was computed, the first `c=2` squeeze
already emptied the candidate set.  The cost is almost entirely the initial
cover search.  For example:

```text
p=89:  |I(13,p,1)|=11804159, find_cover=91.812971s, final=0
p=101: |I(13,p,1)|=12697411, find_cover=113.117640s, final=0
```

This suggests that a proof extension needs a stronger understanding of the
initial bad sets `I(13,p,1)`, not a more elaborate lift after they are found.

## S364/S371 Correction

S364 found that the full scalar-ramp family

```text
v_i = m i mod 14
```

blocks every full micro-staircase cell.  This is not a counterexample to the
program; it is the Dirichlet equality spine in residue form.  S371 upgraded
this to THM-364, which proves scalar-ramp cell blocking directly for all `n`.

The same S371 audit shows that the best known `n=14` non-scalar near-blocker
is a one-coordinate defect of scalar ramp `m=8`, and its `56` missed cells are
exactly scalar cells uniquely blocked by the defect coordinate.  HYP-1817
should therefore be read together with HYP-1818: first excise scalar ramps, use
quotient/divisibility descent for nonunit members, and only then demand a
generic micro-staircase witness for non-scalar vectors.

## S367 Quotient Refinement

S367 turns "excise scalar ramps" into a concrete gauge quotient.  Adding
`m*i` to a residue vector reindexes the full alpha-cell system, so we can
normalize by `v_1=0`.  All scalar ramps collapse to zero.  The new finite
micro-staircase target is therefore:

```text
Every nonzero normalized class has a safe (s, alpha-cell).
```

The exact support scans and full normalized `2`-torsion scan found no full
nonzero blocker.  The unique binary extremal is the coordinate-6 half-turn,
which still misses `56` candidates across eight explicit alpha cells and the
seven odd shifts.

## Test Plan

1. Build an exact SAT/backtracking checker for the full `n=14` micro-staircase
   cell arrangement, not merely random search, after excluding scalar ramps.
2. Classify non-scalar residue vectors by the smallest cell width that resolves them.
3. Prove a uniform lower width bound for all non-scalar, non-gcd tight-lift vectors.
4. If the width bound is too small, isolate the exceptional vectors and compare
   them with the `I(13,p,1)` bad sets produced by the public verifier.
5. Add cover-search pruning based on "this partial tuple already has a
   micro-staircase witness" before the DFS reaches depth 13.
6. Prove HYP-1823's quotient lemma, starting with the scalar-gauge identity and
   the exact `2`-torsion interval stencil.

## Sources

- Sungkawichai and Trakulthongchai, `arXiv:2604.23906`.
- Public verifier: `https://github.com/vzsky/13-lonely-runners`.
- `04-computation/lonely_runner_k13_microstaircase_s363.py`.
- `05-knowledge/results/lonely_runner_k13_microstaircase_s363.out`.
- `04-computation/lonely_runner_k13_scalar_gauge_s367.py`.
- `05-knowledge/results/lonely_runner_k13_scalar_gauge_s367.out`.
- THM-363.
- THM-364.
- `04-computation/lonely_runner_scalar_excision_s371.py`.
- `05-knowledge/results/lonely_runner_scalar_excision_s371.out`.
