---
id: HYP-3530
title: LRC14 multi-denominator rho_D and k<12 attacker-floor scout
status: EVIDENCE / exact finite-denominator and bounded-core attacker scan; not an LRC14 proof
source: codex-2026-06-29 after THM-530 global-witness floor and HYP-3529 R-tail exact-decomposition checkpoint
tangent: T1530
technique: LTI-530
tournament_technique: LTT-430
script: 04-computation/lrc14_multidenom_rhostar_attacker_floor_codex_20260629.py
result: 05-knowledge/results/lrc14_multidenom_rhostar_attacker_floor_codex_20260629.out
reflection: 07-reflections/lrc14-multidenom-rhostar-attacker-floor-codex-20260629.md
related:
  - THM-530
  - HYP-2592
  - HYP-2565
  - HYP-3132
  - HYP-3136
  - HYP-3529
  - OPEN-Q-108
---

# HYP-3530: LRC14 Multi-Denominator rho_D and k<12 Attacker Floor

## Claim

Two proof-frontier objects should now be handled together.

First, replace any single-denominator or via-max target by the finite
multi-denominator global-witness criterion

```text
rho_D(P,E;theta)
  = #{a mod D : a/D in G_P and maxgap({e*a/D:e in E}) > theta}/D,

max_D rho_D(P,E;1/7).
```

This is a finite witness/opportunity bank, not a continuous-measure theorem and
not yet a reconstructed-cover survivor theorem.  It preserves exact rational
candidate witnesses and exposes which denominator families carry the phase
flow.

Second, for the THM-530 union-bound route define an attacker as a normalized
k-shape with

```text
mu_1/7(E) < thr_k := 1 - min_{|P|=13-k} meas(G_P).
```

The exact bounded-core scan through `span <= k+5` for k=8,9,10,11 finds no
attackers, no below-consecutive shapes, and no near-attackers with
`mu < thr_k+1/20`.  Consecutive is the exact bounded-bank minimizer in each
case.  Thus the k<12 union-bound task is now sharply localized: prove every
larger-span row is gentle (`mu_1/7 >= thr_k`) or extend the bank until a known
large-spread lemma applies.

## Exact Readout

Computed by:

```text
04-computation/lrc14_multidenom_rhostar_attacker_floor_codex_20260629.py
05-knowledge/results/lrc14_multidenom_rhostar_attacker_floor_codex_20260629.out
```

Multi-denominator global-witness rows at theta=1/7:

```text
k=8 THM-530 binding:
  rho*=2243/5880, max_D<=196 rho_D=12/13 at D=13

k=9 binding:
  rho*=991/2860, max_D<=196 rho_D=4/5 at D=5

k=10 binding:
  rho*=3803/9555, max_D<=196 rho_D=10/11 at D=11

k=11 binding:
  rho*=278/735, max_D<=196 rho_D=11/12 at D=12

via-max refutation family k=7:
  theta=1/7 rho*=515/1092, max_D<=196 rho_D=10/11 at D=11
  theta=2/7 continuous rho*=0, but D=14 has isolated grid witnesses
```

The last line is a warning, not a contradiction: the finite denominator
criterion can see isolated rational boundary witnesses that have zero Lebesgue
measure.  Therefore `max_D rho_D` should be used as a certificate bank or
route-dispatch layer, not as a stand-alone continuous rho floor.

Exact k<12 attacker scan:

```text
k=8:
  span<=13 primitive_shapes=1716
  min mu_1/7=691/735 = consecutive
  attackers=0, below_consec=0, near_attackers=0
  union floor=1891/5880

k=9:
  span<=14 primitive_shapes=3003
  min mu_1/7=247/294 = consecutive
  attackers=0, below_consec=0, near_attackers=0
  union floor=28117/84084

k=10:
  span<=15 primitive_shapes=5005
  min mu_1/7=38/49 = consecutive
  attackers=0, below_consec=0, near_attackers=0
  union floor=242/637

k=11:
  span<=16 primitive_shapes=8008
  min mu_1/7=1381/2205 = consecutive
  attackers=0, below_consec=0, near_attackers=0
  union floor=10078/28665
```

Large-spread probes were also gentle.  The best probe values were far above
`thr_k` for every k=8..11, with k=11 best probe approximately `0.75745`
against `thr_11=25/91`.

## Proof Pull

The live proof split is:

```text
bounded_core_no_attacker(k<12, span<=k+5)      [exact evidence]
large_span_gentle(k<12)                        [open theorem target]
gp_min_exact_constants                         [known THM-530 data]
---------------------------------------------------------------
unconditional k<12 union-bound floor
```

The multi-denominator route should be formulated as a sidecar-preserving
certificate:

```text
exists D,a with a/D in G_P and maxgap(E*a/D)>1/7
  + denominator family / sidecar alignment
  -> reconstructed survivor candidate
```

The isolated-boundary behavior at theta=2/7 says that the certificate must
retain whether the witness is an open interval point, a boundary point, or a
sidecar-aligned rational.  Raw `max_D` alone forgets that distinction.

## Tournament Analysis

Vertices are proof carriers / denominator families, not runners or arcs.

Pairwise observable:

```text
3*attacker_clearance + 2*witness_support
  + formal_locality - 2*forgetting_penalty
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1]
Hamiltonian path:
  bounded_core_attacker_bank
  -> continuous_union_floor
  -> max_D_rhoD_bank
  -> small_q_survivor_filter
  -> large_spread_gentle_probe
  -> raw_single_D_only
  -> via_max_2over7_shadow
```

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
matroid circuits, and proof obligations.

The selected proof-carrier quotient preserves the LRC predicate "there is a
rational point in the `G_P` plus global-witness carrier."  It destroys
continuous interval volume and reconstructed-cover fast-phase alignment.
Therefore the quotient is useful for routing and certificate discovery, but
the proof still needs a sidecar-alignment lemma or a large-span gentle theorem.
