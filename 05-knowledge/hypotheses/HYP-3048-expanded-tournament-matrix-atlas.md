---
id: HYP-3048
title: Expanded tournament-matrix atlas and sidecar observability principle
status: SYNTHESIS / broad matrix-result translation atlas; not a proof
source: codex-2026-06-26-S212
tangent: T1130
script: 04-computation/tournament_matrix_expansion_atlas_codex_s212.py
result: 05-knowledge/results/tournament_matrix_expansion_atlas_codex_s212.out
related:
  - HYP-3047
  - HYP-3046
  - HYP-3043
  - HYP-3039
  - HYP-3040
  - HYP-3042
  - HYP-2121
  - HYP-2120
  - THM-381
  - THM-385
  - OPEN-Q-108
---

# HYP-3048: Expanded tournament-matrix atlas and sidecar observability

## Claim

Tournament matrix methods should be treated as a carrier factory, not as a
single spectral invariant.

A tournament can be represented at least as:

```text
A             0/1 adjacency
S=A-A^T       skew sign matrix
iS            Hermitian spectral operator
L             directed Laplacian
P             stochastic/ranking kernel
boundary      path, arc, Cech, or proof-obligation boundary matrix
incidence     edge/cycle/owner/sidecar incidence matrix
transfer      automaton, Hamilton path, or threshold-state transfer matrix
game          antisymmetric payoff / majority matrix
observability hidden-coordinate sidecar system
```

The safe rule is the same as the controlled-forgetting rule:

```text
a matrix invariant is legal only if its fibers are route-pure,
status-pure, reconstructible, dual-annihilated, descended by a family lemma,
or routed to named residual debt.
```

S212 expands the S210 matrix atlas by `165` additional classic matrix
results/objects across `14` domains.  Together S210 and S212 give `300` named
matrix hooks for tournament and LRC proof work.

## Evidence

The generated S212 atlas covers:

```text
algebra representation: 12
analysis operators: 12
cs machine learning: 9
factorizations canonical forms: 12
finite fields coding: 10
games social choice: 8
graph combinatorics: 14
linear spectral: 18
matrix scaling positivity: 9
number theory arithmetic: 13
optimization algorithms: 12
physics control: 10
probability random: 12
topology geometry: 14
```

The most reusable rows are not the raw spectral ones.  The highest proof-value
clusters are:

1. **Topology/geometry matrices:** boundary, coboundary, cup product,
   intersection, Seifert, Alexander, Hodge, persistence, rigidity, and stress
   matrices.  These preserve cycle/cocycle and owner-essential information.
2. **Graph/combinatorial matrices:** LGV, all-minors matrix-tree, forest,
   Markov chain tree, critical group, BEST, line-digraph, and incidence
   matrices.  These move from vertex quotients to path/edge/cycle carriers.
3. **Optimization/dual matrices:** KKT, Farkas, LCP, Monge, spectral
   sparsification, leverage, multiplicative weights, SOS, Lasserre, and
   copositive matrices.  These point toward proof certificates rather than
   diagnostics.
4. **Arithmetic/coding matrices:** Smith forms, divisor zeta/Mobius matrices,
   GCD/LCM matrices, character circulants, Krawtchouk/MacWilliams/Delsarte,
   Reed-Solomon/Vandermonde, and LDPC/Tanner matrices.  These keep exact
   period, p-adic, and syndrome sidecars.
5. **Control/observability matrices:** controllability and observability
   Gramians, Riccati equations, scattering unitarity, resolvents, and transfer
   matrices.  These formalize when hidden coordinates can be recovered from
   scalar outputs.

## Tournament Analysis

S212 runs a Tournament Analysis whose vertices are matrix-result domains, not
runners.

Pairwise observable:

```text
retained exactness,
incident/cycle payload,
arithmetic hidden-clock payload,
LRC sidecar usefulness,
computability.
```

Switch:

```text
orient toward the domain that keeps more proof-relevant sidecar information
before scalarizing.
```

Fingerprint:

```text
vertices=14
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1}
directed_3_cycles=0
hamiltonian_paths=1
```

Hamiltonian path:

```text
topology geometry
-> graph combinatorics
-> optimization algorithms
-> number theory arithmetic
-> finite fields coding
-> physics control
-> games social choice
-> algebra representation
-> matrix scaling positivity
-> analysis operators
-> linear spectral
-> factorizations canonical forms
-> probability random
-> cs machine learning
```

This is not saying probability or ML matrices are useless.  It says they are
more often compression or scouting carriers, while topology, combinatorics,
optimization, arithmetic, and coding matrices tend to preserve named proof
sidecars.

## Sidecar Observability Principle

The best theorem target emerging from the atlas is:

```text
Given a target residual fiber, build a matrix O whose rows are pairs of packets
that a coarse quotient identifies, and whose columns are candidate sidecar
coordinates.  O_{ij,c}=1 means coordinate c separates pair ij.

A sidecar set is proof-safe for that fiber only if it has full observability:
every route/status-changing pair is separated, reconstructed, annihilated,
descended, or routed to named debt.
```

This turns HYP-3039's hidden-coordinate ledger into a finite matrix problem.
It also extends HYP-3047: directed-edge sectors, cycle chirality, clique
insertion cuts, and endpoint-owner labels should become columns in such an
observability matrix.

## Matrix Pulls For Tournaments

1. **Four-sector edge block matrix.**  For a directed edge `tail -> tip`, build
   the 4-sector matrix on `(tail beats x, tip beats x)` sectors.  Use it for
   the HYP-3047 perspective defect.
2. **Skew-cycle trace sidecar.**  Use `S=A-A^T`, powers of `iS`, and
   nonbacktracking matrices to retain chirality before scalarizing.
3. **Schur-complement deletion rule.**  Deleting a vertex, edge state, proof
   tooth, or endpoint owner must add the correction term; raw deletion is the
   unsafe quotient.
4. **Smith-normal hidden-clock audit.**  When real spectra merge two packets,
   compute integer Smith forms on incidence/boundary sidecars to expose p-adic
   clocks.
5. **Farkas/SOS impossible-packet certificate.**  Failed searches should become
   dual certificates over packet matrices.
6. **Matrix-tree proof-route redundancy.**  Directed arborescence and forest
   counts on proof-carrier tournaments can measure how robust a certificate
   route is after local labels are forgotten.

## Assumption Challenge

The matrix rows and columns need not be runners.  S212 explicitly treats them
as directed edges, cycles, gaps, wall crossings, endpoint owners, denominator
clocks, proof obligations, quotient fibers, cohomology classes, sidecar fields,
or low-rank update directions.  The preserved predicate is whether the matrix
carrier keeps enough information to decide the tournament or LRC proof
obligation.  The destroyed information must be retained, reconstructed,
dual-annihilated, descended, or routed to named residual debt.
