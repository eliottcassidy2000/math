---
id: HYP-3112
title: LRC14 Lee-Yang ear-payload extremality
status: EVIDENCE / exact scout and proof-route proposal; not a proof
source: codex-2026-06-27-S262b
tangent: T1188
technique: LTI-249
tournament_technique: LTT-147
script: 04-computation/lrc_lee_yang_ear_payload_codex_s262.py
result: 05-knowledge/results/lrc_lee_yang_ear_payload_codex_s262.out
related:
  - HYP-3111
  - HYP-3109
  - HYP-3108
  - HYP-3107
  - HYP-3106
  - HYP-3105
  - HYP-3104
  - HYP-3103
  - HYP-3102
  - HYP-3101
  - HYP-3098
  - HYP-3085
  - HYP-2879
  - THM-577
  - THM-576
  - THM-573
  - OPEN-Q-108
---

# HYP-3112: LRC14 Lee-Yang Ear-Payload Extremality

## Claim

The LRC14 extremality question should not be asked only at the scalar
`p0=G_E(0)`.  It should be asked for the whole miss-count probability
generating function

```text
G_E(z) = sum_t q_t z^t,  q_t = P(number of empty sectors is t),
```

and for the legal extension payloads that move the roots of that function.
This is the exact one-runner payload refinement of HYP-3109's root-curve
ear-map, HYP-3108's broader Lee-Yang/Savitch/Bravais/ear-lattice atlas, and
the HYP-3111 carrier-sidecar reservation, and it uses HYP-3107's Lean
proof-frontier fields as downstream obligations.

If `F=E union {a}` and

```text
A_t(E,a) = P(N_E=t and the new runner a lands in a sector empty for E),
```

then the exact one-runner ear identity is:

```text
q_F[t] = q_E[t] - A_t + A_{t+1}
G_F(z) = G_E(z) + (z^-1 - 1) A_{E,a}(z).
```

Thus `A_t`, not `p0` alone, is the observer-extension / cut payload that
controls root motion.  A quotient that keeps only `q_0`, moments, pair mass,
or the final root count has already forgotten the variable that explains how
the next runner changes the Lee-Yang picture.

## Exact Scout

`04-computation/lrc_lee_yang_ear_payload_codex_s262.py` computes exact
rational miss-count distributions, their roots, a negative-interval
Lee-Yang distance, and one-runner ear payloads.

Stored output:

```text
05-knowledge/results/lrc_lee_yang_ear_payload_codex_s262.out
```

Selected rows:

```text
consec_8/even_AP:
  p0=0.327211, q6=0.020408, p0+q6=0.347619
  real roots=0/6, nearest root modulus=1.4886
  dist(roots,[-1,0])=0.9119, axis gap=25.84 deg

top_cluster:
  p0=0.186703, q6=0.010989, p0+q6=0.197692
  real roots=0/6, nearest root modulus=1.2438
  dist(roots,[-1,0])=0.7179, axis gap=15.25 deg

single_far_21:
  p0=0.217687, q6=0.006803, p0+q6=0.224490
  real roots=0/6, nearest root modulus=1.1467
  dist(roots,[-1,0])=0.2786, axis gap=12.70 deg

break_mid:
  p0=0.045279, q6=0.010989, p0+q6=0.056268
  real roots=2/6, nearest root modulus=0.1212
  dist(roots,[-1,0])=0

random_spread:
  p0=0.011517, q6=0.002646, p0+q6=0.014162
  real roots=4/6, nearest root modulus=0.0472
  dist(roots,[-1,0])=0
```

The AP/consecutive row is not just large at `z=0`.  In this sample it pushes
all roots farthest from the Lee-Yang danger interval `[-1,0]`.

The sampled fugacity curve on `z in [0,3]` is not a single-scalar ranking:

```text
consec_8 wins 52/61 sample points
random_spread wins 6/61
top_cluster wins 2/61
break_mid wins 1/61
winner changes:
  z=0.00 -> consec_8
  z=1.00 -> top_cluster
  z=1.05 -> random_spread
  z=1.35 -> break_mid
  z=1.40 -> top_cluster
  z=1.45 -> consec_8
```

So the proof object is the whole curve and its root motion, not the scalar
`p0` or even the finite moment vector alone.

## Ear-Payload Signal

The exact ear chain separates the AP/consecutive final extension from a
single-far resonance:

```text
nested AP/consec final +7:
  p0=0.327211
  real roots=0/6
  nearest root=1.489
  dist(roots,[-1,0])=0.912
  A_total=0.362585
  A_mean=1.965291
  A_support=1,2,3,4,5,6
  A_even/A_odd=0.110544/0.252041
  reconstruction=ok

single-far final +21:
  p0=0.217687
  real roots=0/6
  nearest root=1.147
  dist(roots,[-1,0])=0.279
  A_total=0.313605
  A_mean=2.993492
  A_support=1,2,3,4,5,6
  A_even/A_odd=0.130612/0.182993
  reconstruction=ok
```

Working readout: good nested AP ears act on a lower miss-count level and push
roots away from `[-1,0]`.  Far resonant ears act on higher miss levels and keep
the root cloud much closer to the danger interval.  Middle-break and spread
rows put roots on the interval itself.

## Two Maps

### Map 1: PGF / Lee-Yang Map

```text
row E
  -> miss distribution q_t
  -> moments S_r and p0
  -> whole curve G_E(z)
  -> root multiset
  -> root metrics:
       real_root_count
       nearest_root_modulus
       lee_yang_negative_interval_distance
       root_axis_gap_deg
       root_modulus_span
       fugacity_winner_profile
  -> gK8 / S2 / convex-order / cap comparison
```

This map preserves the whole analytic object that HYP-3103 exposed.  It
destroys one-step extension information unless `A_t` is retained.

### Map 2: Ear-Payload Decomposition Map

```text
legal growth chain
  -> base E
  -> added runner a
  -> A_t(E,a)
  -> root motion G_E -> G_{E+a}
  -> nested / nonnested / far-resonant ear type
  -> parity and level payload
  -> quotient sidecar or named debt
```

The graph-theory ear facts in the prompt are useful as grammar:

```text
strong directed graphs <-> directed ear decompositions
factor-critical graphs <-> odd ear decompositions
2-vertex-connected series-parallel graphs <-> nested ear decompositions
```

The LRC translation is not a theorem imported from those categories.  It is a
sidecar discipline:

```text
directed ear        -> one-runner extension preserving root motion payload
odd ear             -> parity split of A_t
nested ear          -> AP/consecutive-style legal refinement
nonnested ear debt  -> real-root collision, interval contact, or named route
```

This also reads HYP-2879 in the right direction.  HYP-2879 gives a strong
tournament one-vertex-ear calculus.  HYP-3112 gives the LRC miss-PGF analogue:
the one-runner ear is legal only with its hidden `A_t` payload.

## Proposed New Packet Fields

```text
miss_count_pgf_coefficients
miss_count_pgf_root_multiset
lee_yang_negative_interval_distance
root_axis_gap_deg
root_modulus_span
fugacity_winner_profile
ear_payload_A_vector
ear_payload_mean_level
ear_payload_parity_bias
ear_payload_edge_mass
root_motion_reconstruction_status
nested_ear_status
danger_interval_collision_status
farey_parent_interval
continued_fraction_word
ostrowski_residue_word
root_angle_height_bound
root_angle_separation_certificate
exceptional_low_denominator_resonance
field_of_definition_exit
```

These should be added to the finite-address branch-closure / observer-gluing
ledger whenever the route uses HYP-3103 roots, HYP-3104 maximizer signals,
HYP-3106 perspective functors, HYP-3101 component packets, or HYP-3102
first-obstruction quotients.

## Approximation Merge

The exact ear payload turns irrational/transcendental approximation into a
finite proof-certificate question.  The payload `A_t(E,a)` is computed from
sector incidences on a finite endpoint partition.  Therefore a proposed
irrational witness time, a numerically observed Lee-Yang angle, or a
transcendence analogy must be pulled back to the finite sidecar that decides
which endpoint interval, Farey parent, or root-isolation disk carries the LRC
predicate.

Three guardrails follow.

1. Root angles of `G_E` are algebraic after `E` is fixed, because `G_E` has
   rational coefficients.  Angle proximity to the `7`th-root directions is
   proof-facing only with coefficient height, root isolation, and an explicit
   separation bound.
2. Continued fractions and Ostrowski words are legal schedulers for finite
   root/time localization, not substitutes for endpoint-owner data.  A future
   ledger row should report the Farey parent interval for the last legal ear
   and for any danger-interval contact.
3. Far ears such as `+21` should be read as low-denominator resonance walls:
   they can keep the root set complex while sharply reducing distance from
   `[-1,0]`.  Such rows need `exceptional_low_denominator_resonance` or a
   named debt exit before a quotient is allowed to forget the resonance.

## Bold Conjectures

1. Among binding or low-apex/top-balanced LRC14 packets, the AP/consecutive
   packet maximizes `dist(roots(G_E),[-1,0])` after the expected affine and
   dilation quotients are factored out.
2. The LRC14 cap proof can be reframed as a zero-free strip/interval theorem
   for `G_E`, with the `k=8,9` deviations and known exchange traps appearing
   as root-collision residues rather than as isolated scalar exceptions.
3. Legal proof reductions are root-stability-preserving ear decompositions:
   nested ears are safe, while nonnested or far-resonant ears must emit
   `A_t` sidecars, real-root collisions, first-obstruction syndromes, component
   debt, K33/THM-572 debt, or AP/GW boundary stops.
4. The phi4 density analogy in the prompt should be kept at the correct level:
   `exp(-lambda S^4 - b S^2) dS` is an even distribution where the whole
   transform and zero structure matter, not just one moment.  LRC14 should
   similarly treat `p0` as one evaluation of `G_E`, not as the object.
5. Tournament `H` and LRC `G_E` are parallel partition functions.  Their
   roots should be compared through a transfer ledger, not identified.  The
   comparison belongs downstream of HYP-3105, with the destroyed coordinate
   explicitly named.

## Assumption Challenge

Do not assume tournament vertices are runners or arcs.  For this route the
candidate vertex sets include:

```text
runners
gaps
fixed sectors
sector-boundary events
one-runner ears
payload polynomials A_t
PGF roots
root-motion events
nested-ear obligations
negative-interval collision events
component packets
first-obstruction syndromes
proof obligations
```

Chosen vertices for the next scout should be ear payloads, root-motion events,
and proof obligations.

Preserved LRC predicate: the sector-miss-count PGF under exact one-runner
extension, hence the coverage value `p0` and all factorial moments as
projections.

Destroyed information: labelled interval geometry, endpoint owners, exact
sector incidence, and next-extension root motion unless `A_t`, endpoint
sector data, and owner/cut sidecars are retained.

## Next

Build `lrc14_lee_yang_ear_payload_ledger` over the HYP-2963 packet bank and
the THM-573 `<=6` multiples-of-7 residual.  Start with named rows already in
the observer-gluing, component-bound, cap, O2/O3, K33, and single-far ledgers.
For each row emit `G_E`, root metrics, last legal ear, `A_t`, nested-ear
status, parity/mean payload, destroyed coordinate, and terminal exit.

Then test the concrete conjecture:

```text
if roots approach or meet [-1,0],
then the packet has high-mean ear payload, nonnested ear debt,
component-bound debt, first-obstruction debt, K33/THM-572 debt,
or AP/GW boundary status.
```
