# LRC14 Irrational And Transcendental Approximation Sidecar

Source: codex-2026-06-27-S265  
Status: exact scout evidence and proof-route synthesis; not a proof  
Anchors: HYP-3114, T1190, LTI-251, LTT-149, OPEN-Q-108

## Core Correction

Approximation of irrational or transcendental numbers is not itself an LRC14
proof route.  It becomes useful only after another argument has already
produced an open witness interval or an interior witness margin.

The elementary transfer rule is:

```text
if t is an LRC14 witness with
  delta = min_i(||s_i t|| - 1/14) > 0,
then every rational p/q with
  |t-p/q| < delta / max_i s_i
is also an LRC14 witness.
```

So the proof object is not an irrational constant.  It is the packet

```text
(witness interval, endpoint margin, max speed, robustness radius,
 continued-fraction entry, exceptional approximants, terminal route).
```

## Exact Scout Readout

The dependency-free scout
`04-computation/lrc_irrational_transcendental_approximation_codex_s265.py`
computes exact direct-time witness intervals for named rows, then places
algebraic, transcendental, and Liouville-like targets inside the widest
positive component.

Selected output:

```text
AP13_tight       components=0  measure=0
AP12_tail84      components=8  widest_len=3/1960    grid q>653
loaded_B6        components=64 widest_len=1/5880    grid q>5880
single_tail168   components=8  widest_len=23/11760  grid q>511
```

This confirms two different lessons.

First, AP13 has no positive component, so approximation cannot repair the
boundary-only row.  There is no interior margin for a rational approximation
argument to exploit.

Second, the divisor-loaded THM-575 obstruction has a real positive component,
but in raw time its widest direct component is only `1/5880`.  That agrees
with HYP-3088: raw denominator time is unstable under apex loading.  The
repaired Conjecture 7.1 route must therefore use the normalized THM-565
slow/ruler-coordinate interval, not the original time coordinate.

## Algebraic Versus Transcendental

HYP-3062 already says how to use Roth and Minkowski safely: retain algebraic
degree, height, approximation exponent, exceptional approximants, relation
lattice, covolume, and successive minima.  HYP-3075 adds Hurwitz/Markov/Pell
best-approximant walls.  HYP-3114 supplies the missing split:

- Algebraic irrational targets can have finite-exception fences, but only
  after the height and exceptional approximants are named.
- Transcendental targets do not come with a universal useful bound.  A proof
  needs an explicit irrationality measure or a known approximation sequence.
- Liouville-type targets are warning labels for sparse denominator spikes.
  They are not shortcuts to finite denominator coverage.
- Any positive rational interval already has all sufficiently large
  denominators; this elementary grid fact is the real LRC handoff.

The scout deliberately uses `phi-1`, `sqrt(2)-1`, `e-2`, `pi-3`, and a
truncated Liouville-like target as examples, not as proof carriers.  In
`loaded_B6`, the first continued-fraction hit into the widest component occurs
at denominator `109` or `121`, but the robust radius is only about
`4e-7` to `7e-7`.  That is the scale the proof must retain.

## Tournament Analysis

Tournament vertices were proof carriers, not constants:

```text
finite_interval_margin
continued_fraction_packet
observer_gluing_packet
algebraic_roth_height_fence
transcendental_measure_sidecar
liouville_spike_schedule
raw_named_constant
```

The tournament was transitive:

```text
finite_interval_margin
  -> continued_fraction_packet
  -> observer_gluing_packet
  -> algebraic_roth_height_fence
  -> liouville_spike_schedule
  -> transcendental_measure_sidecar
  -> raw_named_constant
```

The close edges are more important than the ranking.  The interval-margin
packet beats the continued-fraction packet by only one scored axis, and the
continued-fraction packet beats observer gluing by only one axis.  This says
the next useful object is a combined packet, not a choice between arithmetic
approximation and observer charts.

## Next Proof Pull

The next executable step is to replace direct-time widest components by the
normalized THM-565 slow/ruler-coordinate intervals in the HYP-3088/3089 repair.
Each HYP-3098 observer-gluing row and HYP-3112 ear-payload row should then
receive:

```text
witness_interval
endpoint_margin
max_speed
robust_approximation_radius
continued_fraction_first_hit
partial_quotient_spike
irrationality_measure_status
exceptional_approximants
liouville_spike_schedule
terminal_exit
```

That is the way approximation of irrational and transcendental numbers can
advance LRC14: as a controlled-forgetting sidecar converting a proved positive
interval into legal finite-grid evidence.
