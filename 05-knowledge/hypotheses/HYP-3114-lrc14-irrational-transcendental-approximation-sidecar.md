---
id: HYP-3114
title: LRC14 irrational and transcendental approximation sidecar
status: EVIDENCE / exact interval-margin scout and sidecar synthesis; not a proof
source: codex-2026-06-27-S265
tangent: T1190
technique: LTI-251
tournament_technique: LTT-149
related:
  - HYP-3115
  - HYP-3113
  - HYP-3112
  - HYP-3111
  - HYP-3110
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

## Evidence Claim

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

The S265 scout
`04-computation/lrc_irrational_transcendental_approximation_codex_s265.py`
and stored output
`05-knowledge/results/lrc_irrational_transcendental_approximation_codex_s265.out`
make this exact on named LRC14 rows.  The tight AP13 row has no positive
component, so approximation has no place to enter.  The divisor-loaded
THM-575 row `loaded_B6={1,...,11,13,5040}` has `64` positive components,
widest component `(1165/14112, 5837/70560)`, length `1/5880`, and
all-denominator grid bound `q>5880`.  This matches the raw-time denominator
wall from HYP-3088 while also explaining why the repaired route must move to
the normalized THM-565 slow/ruler-coordinate interval rather than the original
time coordinate.

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
root_angle_packet
  -> rational-coefficient PGF root isolation, height, and separation;
bravais_phase_packet
  -> low-denominator resonance walls versus reciprocal-flat phase evidence;
log_gap_packet
  -> Baker-style linear-form gap only after a multiplicative relation exists.
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

## Exact Scout Readout

The computation shows that named constants are not the carrier.  In the
positive rows tested, `phi-1`, `sqrt(2)-1`, `e-2`, `pi-3`, and a truncated
Liouville-like target all enter the widest component through modest
continued-fraction denominators; what matters for proof is the retained packet
containing the component, margin, max speed, first convergent, partial-quotient
spikes, and any irrationality-measure side condition.

Selected exact readout:

```text
AP13_tight       components=0  measure=0
AP12_tail84      components=8  widest_len=3/1960  grid q>653
loaded_B6        components=64 widest_len=1/5880  grid q>5880
single_tail168   components=8  widest_len=23/11760 grid q>511
```

For `loaded_B6`, the first continued-fraction hit into the widest interval
has denominator `121` for `phi_minus_1` and `e_minus_2`, and denominator `109`
for `sqrt2_minus_1`, `pi_minus_3`, and `liouville_10`.  The robust radii are
only about `4e-7` to `7e-7`, so a proof cannot discard the interval-margin
scale even when a low-denominator convergent happens to hit.

## Root, Lattice, And Log-Gap Extension

The incoming HYP-3112/HYP-3113 root-lattice-ear frontier adds three places
where approximation language can become theorem-facing if the right finite
certificate is retained.

First, Lee-Yang root angles are not generic transcendental data.  For a fixed
packet, the miss-count PGF has rational coefficients, so its roots are
algebraic.  A numerical angle close to a `7`th-root direction is useful only
after recording:

```text
root_angle_algebraic_degree
root_angle_height_bound
root_angle_isolating_disk
root_angle_farey_parent_interval
root_angle_separation_certificate
```

Second, Bravais and Minkowski signals should distinguish low-denominator
resonance walls from reciprocal-flat extremal phases.  Current bounded-bank
evidence says high `p0` correlates with residue entropy and against large
Bravais peaks, so the approximation sidecar should measure:

```text
bravais_badly_approximable_score
continued_fraction_word_of_peak_phase
resonance_wall_height
exceptional_low_denominator_peaks
```

Third, Baker-style linear forms in logarithms belong only to lanes that have
already extracted a multiplicative relation, such as a power-resonance gap or a
two-block determinant gap.  Such rows should emit:

```text
linear_forms_log_gap
low_height_resonance_list
multiplicative_relation_lattice
log_gap_exit
```

The rebased HYP-3110 De Moivre/Jacobi/crystallographic lane adds the finite
algebraic/theta version of the same rule.  A De Moivre fold branch, theta tail,
wallpaper quotient, or space-group quotient is approximation-relevant only
after branch id, translation lattice, theta channel, stabilizer word, and
finite-address or observer-gluing exit are named.  In particular, theta tails
belong to finite principal-part / deleted-wall ledgers, not to a generic
"transcendental function" shortcut.

This keeps the pi/e lesson in scope: rationality, irrationality,
algebraicity, and transcendence are not closed-enough proof coordinates by
themselves.  The load-bearing datum is the retained field, height, dependence,
exception list, and finite packet coordinate that survives the quotient.

Post-rebase HYP-3115 supplies the exact anchored-bank
Minkowski/circuit/Ising/De Moivre bridge that can test whether approximation
sidecars align with relation pressure, root strata, and proof-circuit gates.

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

## Tournament Analysis

The scout uses proof-carrier vertices, not constants:

```text
finite_interval_margin
continued_fraction_packet
observer_gluing_packet
algebraic_roth_height_fence
liouville_spike_schedule
transcendental_measure_sidecar
raw_named_constant
```

Its tournament is transitive on the chosen axes:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
priority_path =
  finite_interval_margin
  -> continued_fraction_packet
  -> observer_gluing_packet
  -> algebraic_roth_height_fence
  -> liouville_spike_schedule
  -> transcendental_measure_sidecar
  -> raw_named_constant
```

The close edge flips are also informative: `finite_interval_margin` beats
`continued_fraction_packet` by only one axis, and `continued_fraction_packet`
beats `observer_gluing_packet` by only one axis.  That says the next proof
advance should not choose between interval arithmetic and observer gluing; it
should attach the approximation sidecar to the HYP-3098/HYP-3112 packet route.

## Next Work

1. Replace the direct-time widest-component scout with the normalized THM-565
   slow/ruler-coordinate witness intervals used by the repaired HYP-3088/3089
   route.
2. Attach `witness_interval`, `endpoint_margin`, `max_speed`,
   `robust_approximation_radius`, `continued_fraction_first_hit`, and
   `irrationality_measure_status` to HYP-3098 observer-gluing rows and
   HYP-3112 ear-payload rows.
3. Separate algebraic-height finite-exception fences from transcendental
   irrationality-measure sidecars in any future proof of finite denominator
   coverage.
4. Join the approximation scout to the Lee-Yang ear-payload ledger by adding
   root isolation, Farey parent, separation, and exceptional low-denominator
   resonance fields to rows where roots approach `[-1,0]` or `7`th-root
   directions.
