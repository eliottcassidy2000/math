# Expanded Tournament-Matrix Atlas - S212

This pass extends the incoming S210 matrix atlas rather than replacing it.
S210 gave the core dictionary and 135 matrix hooks.  S212 adds 165 more rows
across 14 domains, for 300 named hooks total.

The main readout is that "tournament matrix" should not mean only adjacency
eigenvalues.  Depending on the proof task, the right matrix may be a skew sign
matrix, Hermitian `iS`, Laplacian, Markov kernel, incidence matrix, boundary
matrix, transfer matrix, game matrix, Smith-normal integer sidecar, KKT/Farkas
dual matrix, or observability matrix over hidden coordinates.

The high-value categories are the ones that preserve sidecar information:
topology/geometry, graph combinatorics, optimization duals, arithmetic/coding,
and control/observability.  Pure spectral, random, or ML matrices are still
useful, but mostly as scouts or compression layers.  They become proof carriers
only after route/status fiber purity is checked.

The concrete theorem target is a sidecar observability matrix.  Rows are pairs
inside a coarse residual fiber.  Columns are candidate hidden coordinates:
edge sector deck, cycle chirality, owner strip, primitive period deck,
endpoint current, Smith invariant, route certificate, and so on.  A sidecar
set is legal only when every route/status-changing pair is separated,
reconstructed, annihilated, descended, or routed to named debt.

This connects directly to HYP-3047.  The A000568 rooted-perspective gap says
node-depth is not enough; the missing eight six-tournament classes live in
incident/cross-coupling payload.  S212 says how to encode that payload:
four-sector edge block matrices, skew-cycle traces, Schur-complement deletion
corrections, and low-rank update fields.

Next useful computation: take the `m=5` directed-edge perspectives from
HYP-3047 and emit sector block matrices plus skew trace features.  Then build
an observability matrix against the `U(6)=56` extension classes and see which
columns account for the eight-class defect.

Assumption challenge: the rows and columns are proof carriers, not necessarily
runners.  The preserved predicate is proof-state separability; the destroyed
data is exact runner identity and full labelled extension rows unless a
sidecar recovers it.
