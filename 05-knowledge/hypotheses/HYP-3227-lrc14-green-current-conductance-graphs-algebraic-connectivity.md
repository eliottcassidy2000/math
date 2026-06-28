---
id: HYP-3227
title: Green-current conductance graphs and algebraic connectivity extend the LRC14 k=8 normal-fan certificate
status: EVIDENCE / exact bounded-bank scout plus proof-target synthesis; not an LRC14 proof
source: codex-2026-06-28
tangent: T1325
technique: LTI-325
tournament_technique: LTT-225
script: 04-computation/lrc_k8_green_current_conductance_graphs_codex_20260628.py
result: 05-knowledge/results/lrc_k8_green_current_conductance_graphs_codex_20260628.out
reflection: 07-reflections/lrc14-green-current-conductance-graphs-codex-20260628.md
related:
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3214
  - HYP-3213
  - HYP-3212
  - HYP-3211
  - HYP-3210
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3162
  - HYP-3161
  - HYP-3160
  - THM-577
  - OPEN-Q-108
---

# HYP-3227: Green-Current Conductance Graphs And Algebraic Connectivity

## Claim

The Green-current half of HYP-3223 is executable and should be kept as a
sidecar inside the HYP-3224 normal-fan certificate.

The useful object is not a single raw conductance scalar.  It is a pair of
conductance graphs plus a trap-discharge graph:

```text
sector graph:
  vertices = six nonzero sectors
  edge weight = positive empty-sector covariance
  negative covariance = named debt

Green precision graph:
  covariance matrix C_E = grounded Green kernel
  C_E^{-1} = grounded conductance/precision matrix
  edge weight = max(0, -precision_ij)
  positive off-diagonal precision = M-matrix defect

trap-discharge graph:
  vertices = HYP-3202 non-AP exchange traps plus certificate coordinates
  edge weight = normalized coordinate deficit
  algebraic connectivity = how tightly trap debt is connected to legal exits
```

This executes the HYP-3223 request to compute effective-resistance and
Schur-complement style data for the HYP-3202 trap manifold, while preserving
the HYP-3224 warning that the proof object is a normal-fan diagram rather than
an extra scalar.

Post-rebase integration: origin/main claimed HYP-3225 for the local
Green-current / Lorentzian trap-fingerprint scout and HYP-3226 for the
small-pattern motif atlas while this conductance-graph scout was local.  The
roles are complementary:

```text
HYP-3225 = local trap fingerprint classifier on 577 trap-neighborhood rows
HYP-3226 = motif/payload atlas reservation
HYP-3227 = full-bank conductance graph and trap-discharge connectivity scout
```

HYP-3225 says Toeplitz remains the universal first trap discharge and splits
the residual trap mechanisms into Green/Rayleigh/Plucker sidecar classes.
HYP-3227 adds the graph-level fact that the trap/certificate conductance graph
stays connected without Toeplitz and even under Green-only coordinates.

Second rebase integration: HYP-3214 identifies the cyclotomic
Delsarte/Beurling-Selberg magic function with the Fejer kernel `F_7`.  That
gives a plausible global dual above the conductance graph: the next test is
whether Fejer/Toeplitz slack dominates the Green-only trap weights and the
precision-defect island.

## Exact Scout Readout

The scout enumerates the same anchored bounded bank as HYP-3205/HYP-3224:

```text
E = {0} union A
A subset {1,...,14}
|A| = 7
rows_all = 3432
rows_primitive = 3431
```

AP/consecutive and its doubled dilation have no beaters for the main
Green-current capacities:

```text
all_ones_green_energy_MAX:     0 beaters, 2 ties
cov_positive_total_weight_MAX: 0 beaters, 2 ties
cov_positive_lambda2_MAX:      0 beaters, 2 ties
cov_positive_min_degree_MAX:   0 beaters, 2 ties
cov_positive_kirchhoff_MIN:    0 beaters, 2 ties
precision_lambda2_MAX:         0 beaters, 2 ties
precision_kirchhoff_MIN:       0 beaters, 2 ties
precision_killing_abs_MIN:     0 beaters, 2 ties
```

For consecutive:

```text
Sigma_kappa2              = 6237419/8643600 = 0.721622819196
all_ones_green_energy     = 2.611956129391
cov_positive_lambda2      = 0.192033074001
cov_positive_kirchhoff    = 108.654718079151
precision_lambda2         = 2.599903805963
precision_kirchhoff       = 6.605115167554
precision_killing_abs     = 13.895458003620
negative_covariance_debt  = 0
```

The guardrail is equally important:

```text
precision_min_degree_MAX:      3 primitive beaters
precision_mmatrix_defect_MIN:  181 primitive beaters
negative_covariance_debt_MIN:  20 ties
```

So the inverse-Green graph is not a clean "AP minimizes every defect" theorem.
The precision M-matrix defect must be retained as a named sidecar.  This is
consistent with the older raw-conductance failures in the repo: conductance is
a true signal but the wrong terminal order parameter when it is scalarized too
early.

## Trap-Discharge Graph

HYP-3224 found that the `11` non-AP arbitrary-exchange traps all have strict
Toeplitz lambda-min deficits.  HYP-3227 refines that into a conductance graph
on proof obligations.

With all coordinates active, the trap/certificate graph is connected:

```text
nodes = 23
edges = 124
lambda2 = 2.719948208
primary_coordinate_counts =
  D1_cov_layer: 4
  cov_positive_lambda2: 3
  precision_kirchhoff: 2
  Toeplitz_lambda_min: 1
  precision_mmatrix_defect: 1
```

Removing Toeplitz still leaves a connected trap graph:

```text
nodes = 22
edges = 113
lambda2 = 2.537866286
```

Even the Green-only trap graph remains connected:

```text
active coordinates =
  cov_positive_lambda2
  cov_positive_kirchhoff
  negative_covariance_debt
  precision_lambda2
  precision_kirchhoff
  precision_mmatrix_defect

nodes = 17
edges = 58
lambda2 = 1.208613477
```

This does not replace HYP-3224's Toeplitz chart.  It shows the Toeplitz chart
is not the only way the finite trap manifold is visible.  Green-current data
connects every non-AP trap to a legal deficit coordinate, and the Fiedler split
isolates a defect island:

```text
positive Fiedler side =
  (0,2,4,6,7,8,10,12)
  (0,1,2,3,7,8,9,10)
  (0,2,5,7,9,10,12,14)
  (0,1,4,5,7,8,11,12)
  negative_covariance_debt
  precision_mmatrix_defect
```

Those four traps should be treated as the first finite Schur-complement /
M-matrix-defect subcase.

## Proof Target

HYP-3227 upgrades HYP-3223 from "try effective resistance" to a concrete
two-level proof obligation:

```text
Level 1: sector/precision networks
  Prove AP minimizes effective resistance or maximizes capacity in the
  admissible Green-current coordinates that actually have no primitive beaters.

Level 2: trap/certificate conductance graph
  Prove every non-AP exchange trap is connected, by positive deficit weight,
  to a legal certificate coordinate after any forbidden scalar coordinate
  is removed.
```

The theorem-shaped statement is therefore:

```text
Every primitive bounded k=8 row either has an exchange/covariance improvement,
or lies in the finite HYP-3202 trap manifold; every non-AP trap has positive
conductance to the HYP-3205/HYP-3224 certificate boundary; AP and doubled AP
are the only sidecar-free equality rows.
```

The next formal move is not to claim algebraic connectivity proves LRC14 by
itself.  It is to express exchange moves as Schur-complement edits of
`C_E^{-1}` and prove that a non-AP local maximum cannot keep all of:

```text
cov_positive_lambda2
cov_positive_kirchhoff
precision_lambda2
precision_kirchhoff
precision_killing_abs
AP support
Toeplitz margin
D1,D2,D3 covariance layers
```

without becoming the AP/dilation equality row.

## Tournament Analysis

Vertices are proof carriers and certificate graphs, not runners, arcs, or raw
sector labels.

Pairwise observable: AP rank, trap-discharge algebraic connectivity, sidecar
retention, and guardrail risk.

Switch/gauge: orient toward the carrier with higher retained proof payload.

Fingerprint:

```text
score_hist = {9:1, 22:1, 48:2, 96:1, 101:1}
directed_3cycles = 0
SCC_sizes = [1,1,1,1,1,1]
edge_flips_vs_raw_conductance_order = 5
hamiltonian_path_count = 1
priority_path =
  trap_discharge_conductance_graph
  -> normal_fan_dictionary
  -> green_without_toeplitz_sidecar
  -> positive_covariance_sector_graph
  -> green_precision_graph
  -> raw_residue_conductance_guardrail
```

## Assumption Challenge

Alternate vertices considered: runners, gaps, sections, boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
proof obligations, inner sector nodes, trap nodes, and certificate-coordinate
nodes.

Chosen vertices: inner sectors for the response networks, and traps plus
certificate coordinates for the discharge graph.

Preserved predicate: k=8 AP/consecutive coverage/covariance extremality,
HYP-3202 trap identities, HYP-3205/HYP-3224 dictionary coordinates, and
Green-kernel sidecar data.

Destroyed information: literal runner identity, raw residue-conductance
monotonicity, and any config-blind algebraic certificate that HYP-3221 warns
will fail.

Challenged assumption: algebraic connectivity is not a terminal scalar; its
proof value is whether it connects finite trap debt to retained sidecar
coordinates.
