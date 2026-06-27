---
id: HYP-3108
title: Lee-Yang, Savitch, Bravais, and ear-lattice extremality signals for LRC14
status: EVIDENCE / exact signal scout and synthesis maps; not a proof
source: codex-2026-06-27-S262
tangent: T1185
technique: LTI-246
tournament_technique: LTT-144
script: 04-computation/lee_yang_savitch_ear_lattice_extremality_codex_s262.py
result: 05-knowledge/results/lee_yang_savitch_ear_lattice_extremality_codex_s262.out
reflection: 07-reflections/lee-yang-savitch-ear-lattice-extremality-codex-s262.md
related:
  - HYP-3107
  - HYP-3109
  - HYP-3106
  - HYP-3105
  - HYP-3104
  - HYP-3103
  - HYP-3102
  - HYP-3101
  - HYP-3096
  - HYP-3089
  - HYP-3088
  - THM-575
  - THM-573
  - OPEN-Q-108
---

# HYP-3108: Lee-Yang, Savitch, Bravais, And Ear-Lattice Extremality Signals For LRC14

## Reservation Claim

This lane reserves the S262 synthesis after integrating incoming HYP-3107's
Lean proof-frontier ledger.  It merges the HYP-3103 miss-count PGF root signal,
the HYP-3104 maximizer atlas, the HYP-3105 obstruction-transfer atlas, the
HYP-3106 controlled-forgetting perspective groupoid, HYP-3107's formal open
fields, and HYP-3109's completed root-curve/zero-real-ear scout with four
additional lenses:

```text
Lee-Yang extremality: the full PGF zero curve, not only a scalar moment.
Savitch reachability: recursive proof-state reachability, not only a flat
  finite check.
Bravais lattice address: finite signal packets as lattice points with
  covolume/successive-minima sidecars.
Ear decompositions: strong/odd/nested ear growth as a legal gluing grammar for
  proof-obligation carriers.
```

The preserved LRC predicate is: a residual LRC14 packet must either have a
large enough direct lonely interval, pass through the level-7 and dyadic lift
ledger, or emit a named obstruction-transfer sidecar.  The quotient deliberately
destroys raw runner labels, raw times, and most absolute speed magnitudes.  The
required sidecars are PGF root stratum, largest-arc/component estimate, CRT
lift status, endpoint-owner packet, and proof-state reachability certificate.

Challenged assumption: the useful tournament vertices need not be runners,
arcs, or configurations.  In this lane the candidate vertices are signal
families, obstruction-transfer states, PGF root strata, lattice cells, ear
attachments, and proof obligations.

## Work To Complete

Completed in the S262 scout for the first finite map, with HYP-3109 supplying
the exhaustive anchored-bank Lee-Yang root-locus, quartic-profile, and
zero-real one-swap ear-component subscout.  Further work is to lift the
strongest shared signals from named packets to the HYP-2963 residual bank.

## Scout Results

Exact scout:

```text
04-computation/lee_yang_savitch_ear_lattice_extremality_codex_s262.py
05-knowledge/results/lee_yang_savitch_ear_lattice_extremality_codex_s262.out
```

Map 1 measures named LRC packets from S66/HYP-3103, HYP-3104 cap atoms, and
wide binding rows.  It records `q_t`, PGF roots, quartic `phi4` fit, `q0`
component count/largest cell, co-emptiness Perron gap, and the cell-state
transition graph as an ear-decomposition proxy.

Key readouts:

```text
consec_8:    #real roots=0, nearest |z|=1.489, apex7 angle gap=4.04 deg,
             q0 components=10, ear_excess=14, Perron gap=0.6469.
break_8:     #real roots=2, nearest |z|=0.121, apex7 angle gap=20.91 deg,
             q0 components=8, ear_excess=40, Perron gap=0.5508.
random_8:    #real roots=4, nearest |z|=0.047, ear_excess=98.
wide k=12:   #real roots=0, nearest |z|=1.824, apex7 angle gap=2.29 deg.
E_star_k12:  #real roots=0, nearest |z|=1.813, apex7 angle gap=2.47 deg.
```

Thus the S66 Lee-Yang signal survives the wider synthesis: root confinement is
not just a scalar `p0` restatement.  It separates all-complex concentrated and
wide rows from root-collision spread/break rows, while the ear/state sidecar
changes at the same time.

The post-merge supplement adds two reproducible checks.  First, the consecutive
zero arc for `k=8..13` stays in the all-complex stratum; the `k=11` middle pair
lands at `102.9°`, exactly the `2*360/7` angle to one decimal place.  Second, a
deterministic anchored 8-set sample gives
`corr(#real,q0+q6)=-0.372`, `corr(nearest-root-modulus,q0+q6)=+0.952`, and
`corr(phi4-lambda,L_yK8)=+0.690`.  The root confinement signal is robust in
this sample.  The quartic `phi4` coefficient is positive here but remains
noisier and sign-sensitive across named packets, so it stays downstream of the
root-curve sidecar.

The Bravais map gives a full affine rank:

```text
rank(q-vectors across all named packets)=6
covolume proxy=0.0020864042
shortest independent q-differences begin:
  gw12 <-> E_star_k12
  break_8 <-> covering_8
  gw10 <-> gw11
  gw11 <-> gw12
  consec_8 <-> gw10
  random_spread_8 <-> covering_8
```

This is a warning against one-moment scalarization: the named frontier packets
occupy the full six-dimensional `q_t` simplex face.

Map 2 turns HYP-3107's proof frontier into a Savitch-style reachability ledger.
The shortest sidecar routes to the terminal proof interface are:

```text
raw_residual_packet -> ear_state_transition -> observer_gluing_certificate
  -> finite_address_packet -> terminal_lrc14_statement

pgf_zero_curve_signal -> quartic_phi4_energy -> bravais_lattice_address
  -> finite_address_packet -> terminal_lrc14_statement

fine_modp_winding_transfer -> observer_gluing_certificate
  -> finite_address_packet -> terminal_lrc14_statement
```

The recursive midpoint schedule repeatedly lands on
`observer_gluing_certificate`, `bravais_lattice_address`, or
`finite_address_packet`.  That matches HYP-3107's open proof fields: the hard
work is still coverage extremality, reflection-Perron/order-3/order-4
certification, Node-3 peel, finite-ruler glue, fine mod-`p` winding, and a
residual classifier that emits observer-gluing or finite-address packets.

Tournament Analysis uses extremality sidecars as vertices, not runners:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
SCCs are singletons
Hamiltonian_path_count=1
leading path:
  savitch_reachability
  -> ear_decomposition_state_graph
  -> pgf_zero_curve
  -> bravais_lattice_address
  -> coverage_component_profile
  -> reflection_perron_moment
  -> obstruction_transfer_ledger
  -> lean_frontier_open_fields
  -> quartic_phi4_energy
  -> raw_scalar_p0
```

Low-margin edges cluster around Bravais/reflection/coverage/obstruction
sidecars, so those are the next place to test edge flips on the full residual
bank.

The rebased S262 script now also includes the bounded `{0}+7` scan over all
`C(13,7)=1716` rows.  This adds the missing full-bank Bravais and strict-local
descent checks to the named-packet map:

```text
corr(p0,nearest_root)      = +0.899
corr(p0,real_roots)        = -0.483
corr(p0,Bravais_peak)      = -0.430
corr(p0,residue_entropy)   = +0.541
corr(p0,phi4_lambda)       = -0.696

#real=0: count=290,  mean_p0=0.10551, max_p0=0.32721
#real=2: count=1426, mean_p0=0.05732, max_p0=0.20833

j=3,n=13 strict-descent traps: 11, reachable midpoint depth 3
j=4,n=13 strict-descent traps: 14, reachable midpoint depth 3
j=5,n=13 strict-descent traps: 3,  reachable midpoint depth 3
j=5,n=16 strict-descent traps: 14, reachable midpoint depth 3
```

The extra bank-level lesson is sharper than the named-packet Bravais-rank
warning: in this bounded bank, high `p0` prefers reciprocal flatness and high
residue entropy, not large Bragg peaks.  The Lee-Yang root wall is value-facing;
Savitch and ears are proof-facing sidecars that explain when local descent is
compressible and when a named trap needs a nonlocal packet.

## Current Hypotheses

Strengthened:

1. `miss_count_PGF_root_stratum` is a live fine-scale coverage invariant, not
   merely another visualization of `p0`.
2. Bravais/lattice rank should be a named sidecar whenever a proof compresses a
   packet to moment scalars.
3. Ear/state-transition data is a plausible observer-gluing sidecar and should
   be tested on HYP-2963 packets before any raw H=7/H=21 tournament analogy is
   reused.
4. Bounded-bank residue spectra suggest a reciprocal-flat Bravais phase for
   coverage extremizers; Bragg-like residue crystallinity is a warning signal,
   not a maximizer certificate.

Still speculative:

1. The quartic `phi4` coefficient `lambda` may be an extremality curvature
   coordinate, but only `6/15` named rows have positive fitted `lambda` in this
   first scout.
2. The close apex-7 root-angle gaps in wide k=11/12 rows may reflect a real
   level-7 zero arc, but they do not yet prove convergence or a zero-free
   region.
3. Savitch midpoint certificates suggest a finite recursive proof-state search,
   but the midpoint theorem remains the missing observer/coverage certificate.
