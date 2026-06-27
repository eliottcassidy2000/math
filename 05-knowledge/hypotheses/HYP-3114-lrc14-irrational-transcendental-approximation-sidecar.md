---
id: HYP-3114
title: LRC14 irrational and transcendental approximation sidecar
status: RESERVED / synthesis and executable scout lane; not a proof
source: codex-2026-06-27-S265
tangent: T1190
technique: LTI-251
tournament_technique: LTT-149
related:
  - HYP-3112
  - HYP-3111
  - HYP-3109
  - HYP-3108
  - HYP-3098
  - HYP-3096
  - HYP-3089
  - HYP-3088
  - HYP-3075
  - HYP-3062
  - HYP-2866
  - THM-575
  - THM-565
  - THM-573
  - OPEN-Q-108
---

# HYP-3114: LRC14 Irrational And Transcendental Approximation Sidecar

## Reservation Claim

This lane merges irrational and transcendental approximation into the LRC14
proof frontier without collapsing it to a scalar "good approximation" slogan.

The core transfer rule is the interior-margin lemma:

```text
if t is an LRC14 witness with margin
  delta = min_i(||s_i t|| - 1/14) > 0,
then every rational p/q with
  |t - p/q| < delta / max_i s_i
is also an LRC14 witness.
```

Thus Diophantine approximation is proof-relevant only after the witness
interval, endpoint distance, max-speed scale, and finite-address route are
retained.  Approximation can convert an irrational or transcendental interior
witness into rational grid witnesses, but it cannot replace the proof that the
interior witness interval exists.

## Sidecar Split

```text
continued_fraction_packet
  -> convergents, partial quotients, first denominator entering a witness interval;
algebraic_irrational_packet
  -> Roth/Hurwitz fence after algebraic target and height are named;
transcendental_packet
  -> explicit irrationality-measure or known approximation sequence required;
liouville_spike_packet
  -> infinite-measure/lacunary denominator warning, not a proof shortcut;
finite_interval_packet
  -> exact LRC component, length, margin, endpoint owners, and grid-hit bound;
observer_gluing_packet
  -> HYP-3098/HYP-3112 route legality and destroyed-coordinate audit.
```

HYP-3062 already handles the Roth-Minkowski algebraic fence; HYP-3075 already
handles Hurwitz/Markov/Pell best-approximant walls.  HYP-3114 adds the missing
distinction between:

- algebraic irrational targets, where Roth-type finite-exception fences can
  exist after height data is retained;
- transcendental targets with finite or unknown irrationality measure, where a
  theorem needs an explicit measure sidecar;
- Liouville-type targets, where too-good approximants create sparse denominator
  spikes and cannot be quotient-forgotten;
- the elementary LRC fact that every positive rational interval contains grid
  points at all large denominators.

## Assumption Challenge

Do not assume the useful vertices are runners, rationals, irrationals, or
named constants.  Candidate vertices are witness intervals, endpoint margins,
continued-fraction states, denominator shells, exceptional approximants,
irrationality-measure claims, Liouville spike schedules, ear payloads, and
finite-address proof obligations.

The quotient preserves the predicate "this witness interval can be converted
to legal finite-grid or observer-gluing data."  It destroys raw time, the
chosen irrational representative, partial-quotient spikes, and exceptional
approximants unless they are retained in the sidecar.

## Immediate Work

1. Build a dependency-free scout that computes exact LRC witness intervals for
   named rows, then places algebraic/transcendental/Liouville targets inside a
   positive component and measures the first continued-fraction convergent that
   stays inside the component.
2. Compare that first-convergent denominator with the exact grid-hit bound
   `floor(1/length)+1` and the robustness denominator predicted by the margin
   lemma.
3. Run Tournament Analysis over approximation carriers, with vertices chosen
   from proof obligations rather than from the constants themselves.
