---
id: HYP-3115
title: Minkowski, circuit-complexity, Ising, and de Moivre sidecars for LRC Lee-Yang extremality
status: EVIDENCE / exact anchored scout plus synthesis hypotheses; not a proof
source: codex-2026-06-27
script: 04-computation/lrc_minkowski_circuit_ising_quintic_bridge_codex_20260627.py
result: 05-knowledge/results/lrc_minkowski_circuit_ising_quintic_bridge_codex_20260627.out
reflections:
  - minkowski-circuit-ising-quintic-bridge-codex-20260627
related:
  - HYP-3116
  - HYP-3111
  - HYP-3114
  - HYP-3113
  - HYP-3112
  - HYP-3110
  - HYP-3062
  - HYP-3109
  - HYP-3108
  - HYP-3107
  - HYP-3106
  - HYP-3105
  - HYP-3104
  - HYP-3103
  - HYP-2963
  - HYP-2614
  - HYP-2613
  - HYP-2612
  - THM-538
  - OPEN-Q-108
external:
  - https://en.wikipedia.org/wiki/Circuit_complexity
  - https://mathworld.wolfram.com/deMoivresQuintic.html
  - https://mathworld.wolfram.com/IsingModel.html
  - https://encyclopediaofmath.org/wiki/Minkowski_theorem
---

# HYP-3115: Minkowski, Circuit, Ising, And de Moivre Sidecars

## Claim

The four new cues should be merged into the LRC14/Lee-Yang proof search as a
single controlled-forgetting packet:

```text
mciq_packet =
  (relation_lattice_short_walls,
   relation_covolume_or_successive_minima,
   root_stratum_ising_domain_wall,
   proof_signal_circuit_basis_size_depth_fanin,
   stationary_de_moivre_quintic_residual,
   uniformity_or_tailored_threshold_debt)
```

The point is not that any one analogy proves LRC.  The point is that each one
names a coordinate a scalar quotient can silently destroy:

- Minkowski's theorem gives a lattice-volume/short-vector gate, but HYP-3062
  already says the relation lattice, covolume, convex body, and low-height wall
  ledger must be retained before volume pressure can become proof currency.
- Circuit complexity asks for the Boolean circuit computing or certifying a
  proof predicate: basis, size, depth, fan-in, and uniformity.  A finite fitted
  threshold is not a proof circuit.
- The Ising model turns the Lee-Yang root stratum into a spin field.  The
  one-swap boundary between `#real=0` and `#real=2` rows is a domain wall.
- de Moivre's quintic suggests looking at the stationary polynomial `G'(z)`.
  Since `G` has degree six, `G'` is a quintic; its distance from translated
  de Moivre normal form is an algebraic-solvability stress signal.

Post-rebase integration note: incoming HYP-3110 owns the De Moivre/Jacobi/
crystallographic frontier, HYP-3111 owns the broad
Minkowski/circuit/Ising/De Moivre carrier atlas, HYP-3112 owns the Lee-Yang
ear-payload extremality ledger, HYP-3113 owns the two-map root-lattice-ear
synthesis, and HYP-3114 owns the irrational/transcendental approximation
sidecar.  This HYP supplies the exact anchored-bank executable subscout that
links those lanes.

## Exact Scout

Artifact:

```text
04-computation/lrc_minkowski_circuit_ising_quintic_bridge_codex_20260627.py
05-knowledge/results/lrc_minkowski_circuit_ising_quintic_bridge_codex_20260627.out
```

The scout recomputes the anchored bank
`0 union A`, `|A|=7`, `A subset {1,...,13}` (`1716` rows), then adds four
signals to the HYP-3109 root-locus metrics.

Named readout:

```text
consec_8:
  p0=0.32721, #real=0, clearance=0.91189, nearest=1.48858
  apex7_error=4.042, pair_relations=22
  minkowski_pressure=1.85934, axis_walls=1
  de_moivre_residual=0.63368

nearest_root_competitor (0,2,4,6,7,8,10,12):
  p0=0.25578, #real=0, clearance=0.80575, nearest=1.55621
  apex7_error=10.176, pair_relations=13
  minkowski_pressure=0.63969, axis_walls=1
  de_moivre_residual=0.66183
```

Minkowski-style finite relation proxy:

```text
corr(p0, pair_relation_count) = +0.427
corr(p0, minkowski_pressure)  = +0.407
top relation leader = consec_8 with pair_relations=22
```

This is consistent with the older support-six lesson: short relations are not
bad by themselves.  The extremal row is relation-rich, but the proof must keep
which relations are AP/consecutive structure, which are low-height walls, and
which are deleted anti-cosets.

Ising root-stratum graph:

```text
one_swap_edges = 36036
spin_counts = {+1: 290, -1: 1426}  (+1 means #real=0)
same_edges = 25952
domain_wall_edges = 10084
wall_fraction = 0.27983
Ising_energy_J1 = -15868
magnetization = -1136
```

This makes the HYP-3109 boundary a literal finite spin interface.  The next
proof target is not just "the zero-real component is connected"; it is to label
which domain-wall edges are legal root-collision exits and which are forbidden
for an extremal packet.

Circuit-complexity sidecar:

```text
target=max_p0_rows, max_p0=0.32721
minimal observed non-p0 conjunction = [apex7_error <= 5]
input_literals = 1
false_positives = 0
false_negatives = 0
```

This is a warning, not a theorem.  On this bounded bank, a single tuned
apex-7 threshold isolates `consec_8`.  A real proof needs a uniform circuit
statement: fixed basis, legal inputs, allowed fan-in/depth, and a reason the
threshold is not just fitted to the finite bank.

de Moivre derivative-quintic sidecar:

```text
corr(p0, -de_moivre_residual) = +0.348
best residual row = (0,1,2,3,4,5,6,9), residual=0.43210, p0=0.18810
consec_8 residual = 0.63368
```

The de Moivre signal is weak and should not be promoted as an extremality
criterion.  Its useful role is different: `G'(z)` controls stationary/root
collision geometry, so a translated de Moivre residual may flag algebraic
branch points or solvable local models near a collision wall.

## Hypotheses To Test

1. **Minkowski plus Lee-Yang, not Minkowski alone.**  Extremality is the
   intersection of high short-relation pressure with zero-free root geometry.
   Relation-rich rows that fail the root sidecar are low-height wall debt, not
   proof exits.

2. **Ising droplet theorem.**  The `#real=0` stratum is a positive-spin droplet
   inside the larger `#real=2` phase.  The root-collision proof should classify
   the `10084` boundary edges by labelled ear move, forbidden wall, or legal
   observer-gluing discharge.

3. **Circuit-uniformity guardrail.**  Any scalar-looking maximizer rule must be
   replaced by a uniform circuit over named sidecars.  A one-literal finite-bank
   classifier is evidence for a signal, not evidence for a proof.

   HYP-3116 sharpens the missing sidecars: the circuit row must retain
   `endpoint_cover_activation_vector`, `phi_gap_sum`, `phi_kernel_status`,
   `P_max_activation`, `endpoint_period_numerator_sidecar`,
   `finite_address_packet`, and `observer_gluing_certificate`, or name the
   first missing input as residual debt.

4. **Stationary quintic stress.**  The de Moivre residual of `G'(z)` is not an
   extremality objective, but low residuals may mark solvable local normal forms
   for root collisions.  Test it on paths crossing from `#real=0` to `#real=2`,
   not only on static rows.

5. **Two-map extension.**  HYP-3108's broad Savitch/Bravais map and HYP-3109's
   anchored root map should be extended with four sidecars:
   `short_relation_wall_class`, `spin_domain_wall_id`,
   `proof_circuit_signature`, and `stationary_quintic_residual`.

## Tournament Analysis

Vertices are proof sidecars, not runners:

```text
root_zero_locus
ising_domain_wall_tension
minkowski_relation_lattice
circuit_depth_sidecar
de_moivre_quintic_residual
bravais_q_lattice_rank
observer_gluing_certificate
single_scalar_p0
```

Pairwise observable: which sidecar better preserves the LRC predicate, retains
finite witnesses, detects phase walls, is certifiable, avoids controlled
forgetting, is computable, carries algebraic structure, and exposes lower-bound
or uniformity debt.

Fingerprint:

```text
score_hist = {0:1,2:3,4:1,5:1,6:1,7:1}
directed_3cycles = 1
Hamiltonian_path_count = 3
SCC with three weakly separated sidecars =
  {bravais_q_lattice_rank, de_moivre_quintic_residual,
   ising_domain_wall_tension}
priority path =
  root_zero_locus > minkowski_relation_lattice > circuit_depth_sidecar
  > observer_gluing_certificate > ising_domain_wall_tension
  > de_moivre_quintic_residual > bravais_q_lattice_rank
  > single_scalar_p0
```

Assumption challenge: considered vertices included runners, gaps,
relation-lattice vectors, Boolean gates, spin domains, stationary quintic
coefficients, PGF roots, Bravais cells, and proof obligations.  The selected
vertices are proof sidecars.  This preserves predicate-retention order and
destroys raw row identity, deliberately.

## Next Work

- Add the four new fields to a HYP-2963 residual sample.
- Replace the pair-relation proxy by an exact relation-lattice/covolume/
  successive-minima ledger for support-six packets.
- Run the Ising wall classifier on actual one-swap paths from `#real=0` to
  `#real=2`, with ears labelled by root motion.
- Re-run the circuit scan with a fixed basis and train/test split over a larger
  row family; reject thresholds that only fit the anchored bank.
- Measure de Moivre residual along root-collision paths and compare it to the
  PGF discriminant/resultant wall.
