---
id: THM-571
title: LRC(14) apex-majority gamma descent -- every primitive 13-speed row with at least seven multiples of 14 is lonely
status: PROVED modulo the accepted LRC<=13 input
source: codex-2026-06-22-S122
depends_on:
  - THM-570   # 14-phase shift guard for zero/one half-step residual
  - THM-523   # q-witness / covering-set reduction
  - LRC<=13   # recent below-frontier theorem, used only for <=12 scaled speeds
related:
  - THM-568
  - THM-569
  - HYP-2910
  - HYP-2911
  - OPEN-Q-108
results:
  - 04-computation/lrc14_apex_majority_gamma_descent_codex_s122.py
  - 05-knowledge/results/lrc14_apex_majority_gamma_descent_codex_s122.out
---

# THM-571 -- LRC(14) apex-majority gamma descent

## Statement

Let `S` be a primitive 13-speed row and suppose at least seven speeds in `S`
are divisible by `14`.  Then `S` has a lonely time at threshold `1/14`.
In fact `M(S)>1/14`.

This closes the apex-majority branch left open by THM-568/THM-570.

## Proof

Write

```text
S = 14Q union R,        r = |Q| >= 7,  |R| = 13-r <= 6,
```

where no speed in `R` is divisible by `14`.  Since `S` is primitive, `r<=12`.

### Case 1: at most one half-step residual

Assume at most one speed in `R` is divisible by `7`.  This is exactly THM-570.
The below-frontier LRC theorem gives a strict `Q`-safe phase `u`.  All fourteen
lifts

```text
t = (u+k)/14,       k=0,...,13,
```

keep the `14Q` speeds safe.  Ordinary residual speeds forbid at most two lifts
and at most one in each parity; one half-step speed forbids at most one whole
parity class.  Thus at most `12` of the fourteen lifts are forbidden, so one
survives.

### Case 2: at least two half-step residuals

Now let `H` be the full set of speeds in `S` divisible by `7`.  It contains all
`14Q` speeds plus at least two residual half-step speeds, so

```text
|H| >= 9.
```

Because `S` is primitive, not every speed is divisible by `7`, hence
`|H| <= 12`.  Scale by `7`:

```text
H = 7P,       |P| <= 12.
```

By the accepted LRC theorem below 14 runners, there is a phase `v` with

```text
||p v|| >= 1/(|P|+1) > 1/14      for every p in P.
```

Again use continuity to take this strict on a small interval.  For each
`j=0,...,6`, put

```text
t = (v+j)/7.
```

Every speed in `H` is safe at all seven lifts:

```text
||7p t|| = ||p(v+j)|| = ||p v|| > 1/14.
```

The remaining speeds are not divisible by `7`.  There are at most

```text
13 - |H| <= 4
```

of them.  For a speed `b` not divisible by `7`, the seven values
`b(v+j)/7` run through a shifted copy of the seven equally spaced residues.
The forbidden arc around the integers has open length `2/14 = 1/7`, so it
contains at most one of these seven points.  Hence each non-`7` residual speed
forbids at most one lift, and at most four of the seven lifts are forbidden.
Choose a surviving lift.  After a generic perturbation inside the strict
`P`-safe interval, all inequalities are strict, so `M(S)>1/14`.

The two cases cover every primitive row with at least seven multiples of `14`.

## What This Does And Does Not Prove

This theorem eliminates the apex-majority residual.  It does **not** by itself
complete LRC(14), because rows with at most six multiples of `14` still depend
on the separate scale-separated comb/finite-core branch (`S31v` plus its
bounded-core census).  The recent LRC proof through 13 total runners supplies
the below-frontier margin used above, but it does not automatically classify
the 13-moving-speed tight locus.

## Relation To HYP-2911

HYP-2911 found a real obstruction to the raw fourteen-shift sieve: two
half-step speeds can cover both parities.  THM-571 resolves that obstruction
inside the apex-majority branch by changing phase: once two half-step residuals
exist, at least nine speeds are multiples of `7`, so the problem descends to
seven lifts, where the non-`7` residual has size at most four.
