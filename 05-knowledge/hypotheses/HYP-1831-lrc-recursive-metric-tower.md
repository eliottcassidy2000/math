---
id: HYP-1831
status: EXPLORATORY
source: codex-2026-05-31-S376
related:
  - HYP-1794
  - HYP-1802
  - HYP-1810
  - HYP-1827
  - HYP-1828
  - HYP-1829
---

# HYP-1831: The Lonely Runner frontier has a recursive metric tower

## Statement

The Lonely Runner Conjecture should be studied across the threshold
denominator `n`, not only at one fixed case.  The right analogue of the
repo's tournament structure program is a recursive metric tower:

```text
n
|-> cell-arrangement complexity C(n)
|-> unit-boundary skeleton phi(n)/(n-1)
|-> scalar-puncture moat curvature
|-> endpoint-peeling depth and closure defect
|-> quotient-ladder / gate near-counterexample pressure
```

Prime and composite steps affect different layers.  Prime `n` maximizes the
unit boundary skeleton, while composite `n` introduces quotient layers,
nonunit scalar-puncture deltas, and quotient-ladder near-disproofs.  A proof
strategy should track how these layers transfer from `n` to `n+1`, much as
the tournament program tracks class counts, endpoint-transfer ranks, and
projection defects across vertex insertion.

## Evidence

`lonely_runner_recursive_metrics_s376.py` computes exact metrics for
`n=2..18`.

The micro-staircase arrangement complexity

```text
C(n) = # alpha floor-pattern cells
```

grows with divisor-sensitive jumps:

```text
n:    11   12   13   14   15   16   17   18
C(n): 352  504  598  812  960  1152 1360 1728
dC:    72  152   94  214  148   192  208  368
```

The initial-segment unit skeleton is exactly visible in the endpoint metrics.
For every checked `n`, the initial segment is boundary-only with unprotected
unit endpoints, and the unit density is

```text
phi(n)/(n-1).
```

Primes in the range have density `1.0`; composites expose quotient layers.

The scalar-puncture moat behaves like a local curvature invariant around the
Dirichlet equality spine.  Nonunit best-puncture deltas occur at

```text
n = 4, 6, 8, 9, 10, 12, 14, 15, 16, 18.
```

The familiar frontier values appear in the same metric:

```text
n=14: moat=56,  delta_gcd=7
n=15: moat=120, delta_gcd=5
n=16: moat=128, delta_gcd=8
n=18: moat=198, delta_gcd=9
```

The counterexample side has its own recursive layer.  Quotient ladders with
tiny gap ratio below `0.02` occur at

```text
n = 10, 12, 14, 15, 16, 18.
```

but they carry large endpoint defects.  For example:

```text
n=14: ladder d=7, skip=6, gap/th=0.005411, unprotected=84
n=16: ladder d=8, skip=14, gap/th=0.007576, unprotected=140
n=18: ladder d=9, skip=8, gap/th=0.005682, unprotected=176
```

This extends HYP-1829: the scalar-distance / endpoint-closure interaction is
not merely a fourteen-runner accident; it is a recursive composite-layer
phenomenon.

## Interpretation

The tournament analogy is precise enough to be useful:

- tournament `n` has class counts; LRC has arrangement complexity `C(n)`;
- tournament insertion has endpoint-transfer ranks; LRC has endpoint-peeling
  depth and closure defect;
- tournament special families such as Paley/interval classes expose phase
  behavior; LRC quotient ladders expose tiny-gap but endpoint-defective
  families;
- tournament projection defects warn that scalar counts miss incidence; LRC
  scalar moat counts must be paired with endpoint incidence.

So the LRC object should be treated as a tower of finite incidence systems.
The reduced problem at `n+1` is not just a harder version of the problem at
`n`; it changes the unit group, the divisor lattice, the alpha-cell
arrangement, and the available endpoint-protection gates.

## Predictions

1. The first useful broad `n`-recursive theorem will not be a direct
   monotonicity theorem; it will be a transfer theorem with separate prime and
   composite cases.
2. Composite `n` with a large proper divisor near `n/2` should keep producing
   tiny quotient-ladder gaps and large endpoint defects.
3. Prime steps should look "rounder" on the unit-boundary layer but can still
   have large scalar moats, because there is no quotient ladder to absorb the
   defect.
4. The frontier `n=14` is difficult because it sits at the intersection of a
   small scalar moat density, a strong half-divisor quotient ladder, and a
   nontrivial endpoint-closure defect.
5. A future LRC "complete invariant catalogue" should include at least:
   `C(n)`, denominator spectrum, unit density, scalar moat, puncture
   gcd-profile, gate susceptibility, ladder gap ratio, unprotected endpoint
   count, and peel-depth profile.

## Next Tests

1. Optimize S376's gate channel and extend the atlas to `n=24` or `n=30`.
2. Derive a closed formula for `C(n)` from the denominator arrangement
   `{a/(n i): 1 <= i < n}`.
3. Build an LRC endpoint-transfer matrix from `n` to `n+1`, using unit
   endpoints and first quotient leaks as rows/columns.
4. Compare scalar-puncture gcd profiles to the largest proper divisors of
   `n`.
5. Rank future fourteen-runner disproof attempts by the recursive vector
   `(ladder gap ratio, unprotected endpoints, scalar moat density)`.

## Sources

- `04-computation/lonely_runner_recursive_metrics_s376.py`
- `05-knowledge/results/lonely_runner_recursive_metrics_s376.out`
- `07-reflections/lonely-runner-recursive-metrics-s376.md`
- HYP-1829
