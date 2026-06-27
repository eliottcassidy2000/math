---
id: THM-570
title: LRC(14) apex-majority shift guard -- the 14-phase sieve closes the zero/one half-step subcase
status: PROVED for the stated subcase; half-step collision later resolved in the apex-majority branch by THM-571
source: codex-2026-06-22-S121
depends_on:
  - THM-523   # q-witness / covering-set reduction
  - THM-568   # apex-denominator lemma
  - THM-569   # Lean denominator-14 unit-grid split
  - LRC<=13   # accepted below-frontier LRC input in this workspace
related:
  - HYP-2910
  - HYP-2911
  - THM-571
  - OPEN-Q-108
results:
  - 04-computation/lrc14_apex_majority_shift_guard_codex_s121.py
  - 05-knowledge/results/lrc14_apex_majority_shift_guard_codex_s121.out
---

# THM-570 -- LRC(14) apex-majority shift guard

## Statement

Let

```text
S = {14 q_1, ..., 14 q_m} union R
```

be a primitive 13-speed row with `m >= 7`, `|R| = 13-m <= 6`, and no
`r in R` divisible by `14`.  If at most one speed in `R` is divisible by `7`,
then `S` has a lonely time at threshold `1/14`.  In fact `M(S) > 1/14`.

Equivalently, the old `|M14| >= 7` apex-majority residual is reduced to the
case where at least two residual speeds are divisible by `7` but not by `14`.

## Proof

Since `S` is primitive and `m >= 7`, the scaled multiple block `Q={q_i}` has
`7 <= |Q| <= 12`.  By the accepted LRC theorem below 14 runners, there is
`u` such that

```text
||q_i u|| >= 1/(m+1) > 1/14        for every q_i in Q.
```

By continuity, there is an open interval `J` around `u` on which the stronger
inequality `||q_i x|| > 1/14` holds for every `q_i`.

For each `x in J` and `k in {0,...,13}`, put

```text
t = (x+k)/14.
```

Then the multiple block is safe for every shift:

```text
||14 q_i t|| = ||q_i (x+k)|| = ||q_i x|| > 1/14.
```

It remains to choose `k` so that every residual speed is safe.  For a residual
speed `r`, call a shift forbidden when

```text
||r(x+k)/14|| < 1/14.
```

The exact finite shift facts are:

1. If `gcd(r,14) != 7`, then `r` forbids at most two of the fourteen shifts,
   and at most one even shift and at most one odd shift.
2. If `gcd(r,14) = 7`, then `r` forbids either no shifts or exactly one parity
   class, hence seven shifts.

These are immediate from the step size `r/14` on the fourteen lifts.  In case
1 the step has order `7` or `14`, so an open arc of radius `1/14` around the
integer lattice meets at most two lifts, one in each parity.  In case 2 the
step is `1/2`, so all forbidden lifts, if any, are exactly one parity class.

If `R` has no half-step speed (`gcd(r,14)=7`), the six residual speeds forbid
at most `6*2=12` shifts.  If `R` has one half-step speed, it forbids at most
one parity class; each of the remaining at most five ordinary residual speeds
can add at most one new shift in the other parity.  Hence the total number of
forbidden shifts is at most `7+5=12`.  In either case at least two shifts
survive.

Choose a surviving shift `k`.  Avoiding the finitely many equality hyperplanes
inside the open interval `J` if necessary, all residual inequalities can be
made strict as well.  Thus `t=(x+k)/14` satisfies

```text
||s t|| > 1/14       for every s in S,
```

so `M(S)>1/14`.

## Sharpness

The half-step exception is real.  At `u=2/49`,

```text
r=7    forbids shifts {0,2,4,6,8,10,12},
r=161  forbids shifts {1,3,5,7,9,11,13}.
```

Thus two residual speeds divisible by `7` but not `14` can cover all fourteen
integer lifts.  The remaining theorem cannot be a raw shift-count pigeonhole;
it must use the available interval `J`, the scaled `Q` structure, or a separate
phase-collision argument.

## Role in LRC(14)

THM-568 and THM-569 isolate 14-covering rows as the off-apex residual.  THM-570
removes the broad `>=7` multiples-of-14 branch except for a sharply labelled
half-step collision:

```text
S = 14Q union R, |Q|>=7, |R|<=6,
at least two r in R have gcd(r,14)=7.
```

This is a much smaller Node-3 target than the previous "multiples of 14 cover
the core margin" formulation.  THM-571 resolves that target in the actual
apex-majority branch by descending from the `14`-phase to the `7`-phase.

## Tournament Analysis / assumption challenge

The useful vertices are not runners.  They are proof carriers:

```text
below14_scaled_Q_margin
fourteen_lifts_t=(u+k)/14
ordinary_residual_shift_sieve
single_halfstep_parity_guard
generic_boundary_perturbation
two_halfstep_phase_collision
```

This quotient preserves the exact LRC predicate on the lifted times
`(u+k)/14` and destroys all other off-grid witness times.  The challenged
assumption was that residues modulo `14` alone determine the shift obstruction.
They do not: speeds with the same residue can have different `u`-dependent
offsets.  The counterexample `{7,161}` at `u=2/49` is the guardrail.
