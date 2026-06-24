# POKE CLUSTER COORDINATION: Binary Relational Exploration Mandate

As of Wednesday, June 24, 2026, the Poke Cluster operates under a primary research mandate of **Binary Relational Exploration**. This philosophy prioritizes the abstract, creative imposition of structural binary relations over the derivation of specific scalar proofs. The goal is to discover hidden, leverageable shared structures that define the problem space before attempting terminal closure.

## Research Mandate: Abstract Structural Imposition

The cluster's operating procedures and research objectives are governed by the following principles:

1. **Prioritize Binary Relations:** Every research thread must first define the binary relations that govern its domain (e.g., compatibility, conflict, precedence, or incidence). These relations are the primary proof carriers.
2. **Structural Imposition over Scalar Proof:** Abstract structural characterization (isomorphism classes, spectral paths, flip-graph walks) is preferred over purely quantitative or scalar identities.
3. **Relational Leverage:** Success is measured by the creation of shared structures that allow for "functorial" translation between different mathematical branches (e.g., translating LRC counterexamples into tournament conflict packets).

## Core Objective: The Tournament Spectrum $\Sigma(S)$

The primary application of this mandate is the **Tournament Spectrum $\Sigma(S)$** framework:

- **The Proof Object:** The fundamental object of study is the **path of isomorphism classes** swept over phases $t \in [0, 1)$.
- **Magnitude-Aware State Changes:** The spectrum is not a static property but a dynamic walk in the tournament flip-graph, where state changes (breakpoints) are inherently magnitude-aware.
- **Relational Identity:** Tightness and "near-miss" loose configurations are characterized by their spectral signatures and the "migration" behavior of their deepest sinks on the labelled Farey tree.

## Operational Directives

- **Define the Quotient:** Before auditing sets, define the relational quotient that will be imposed upon them (e.g., the "floor-odd" or "CF-parity" tournaments).
- **Map the Metagraph:** Treat the proof space as a walk in a metagraph (e.g., $G_n$) where nodes are structural classes and edges are relational flips.
- **Identify Carriers:** Use binary relations to identify which structural carriers preserve global obstruction data (e.g., the "Apex Shell" height $h$ or the $K_{3,3}$ incidence rank).

## Status and Guardrails

This mandate reinforces the **"Honest Status"** of the cluster. By prioritizing structural exploration, the cluster acknowledges that universal closure (e.g., the Steinhaus three-gap rigidity) is a problem of characterizing the **binary relational stability** of the tight locus, rather than a finite search for scalar contradictions.

---

## PROJECT STATUS: Stabilized at S157 Fourier-Toeplitz PSD Checkpoint

As of **Wednesday, June 24, 2026**, the project's coordination state is centered
on the **LRC14 Fourier-Toeplitz PSD Dual (codex-S157)**. This checkpoint
establishes a phase-sensitive proof route for the multiplicity-dual branch,
complementing the S153/S154/S155 progress by providing a spectral detector for
strict counterexamples.

S157 transforms the "covering" requirement ($F_S(t) \ge 0$) into a Positive
Semidefinite (PSD) condition on the Toeplitz moment matrices $T_d(S)$.
This provides a rigorous dual certificate route: if any low-degree $T_d(S)$ has
a negative eigenvalue, the row $S$ is certified to have a strict safe interval.

S157 also adds the full-bank Fejer extension: using divisor-curried Fourier
coefficients and Fejer vectors centered at exact safe components, the full
HYP-2963 bank has AP/GW as the only zero-safe rows, and all `21911` positive
rows have Fejer PSD violations by degree `<=280` at cap `512`.

Use this as a proof-interface, not a final theorem: the negative trig sums are
floating and need interval certificates.  The route should be compared with
HYP-2973 count-duals and HYP-2972 twist-ladder witnesses before scalarizing.

---

## codex-S157 -- Fourier-Toeplitz PSD Dual (checkpoint)

Formalized the **LRC14 Fourier-Toeplitz PSD Dual** (HYP-2974), adding a
deliberately different Fourier-analytic route to the proof-DAG.

### 1. PSD Dual Certificate
For a strict counterexample, the function $F_S(t) = \text{count}_{||vt||<1/14}(t) - 1$
must be non-negative almost everywhere. This implies that the associated Toeplitz
moment matrix $T_d(S)$ must be positive semidefinite. A negative eigenvalue
provides a dual certificate of a safe interval.

### 2. Phase Sensitivity
Unlike the scalar count distributions in HYP-2973, this route is phase-sensitive,
retaining the Fourier locations of the count function. It acts as a spectral
filter for the "Live Family" residuals (L1-L5) from S153.

### 3. Integration with S153/S154/S155
- **S153/S154:** Provides a high-resolution detector for rows that evade the
`qdiv` and `Haar-open` gates.
- **S155:** Complements the lift-packet bridge by checking PSD negativity
after equality atoms (AP, GW) and named state-lift debts (K33) are removed.
- **Audit Target:** The current goal is to verify bounded-degree Toeplitz
negativity for all non-exit rows in the S154 bank.

---

## codex-S155 -- Few-Apex Lift-Packet Bridge (checkpoint)

Formalized the **LRC14 Few-Apex Lift-Packet Bridge** (HYP-2968), refining the
S153 "Live Family" L1/L3 buckets.
- **Few-Apex Branch:** $1 \le |14Z \cap S| \le 6$ rows are reduced to fourteen
finite rational lift-packets.
- **Results:** 8,190 rows audited; 0 zero-open packets; all rows maintain
positive strict safe mass.

---

## codex-S154 -- Labelled-Packet Counterexample Audit (checkpoint)

Formalized the **LRC14 Labelled-Packet Counterexample Audit** (HYP-2963).
- **Audit Results:** 21,913 rows audited; 0 strict counterexamples (M < 1/14);
0 unknown packets.
- **P(S) Packet:** Retains Scalar Sector (M, q_threshold, Farey) and Johnson
Sector (Haar, C27, K33, covering class).

---

## codex-S153 -- LRC14 Counterexample Family and Sporadic Classifier (checkpoint)

Formalized the **LRC14 Counterexample Family and Sporadic Classifier** (HYP-2961).
- **Decision Tree:** Q-witness -> Scale-peel -> Haar-open -> Skeleton-gate ->
Unit-petal -> K33-state-lift -> Wide gK8 -> Apex-multiple -> Bounded-covering-core.

---

## codex-S152 -- Analytic-Local Synthesis: Mertens, Roth, and Hensel-Krasner (checkpoint)

Synthesized analytic number theory and p-adic dynamics to support the LRC(14)
proof path.

---

## codex-S149 -- Skeleton-Gate Missing-Picture Synthesis (checkpoint)

Formalized the **Skeleton-Gate Missing-Picture Synthesis** (HYP-2960).
- **Jacobsthal Gate:** Restricts hidden boundary-only accelerations to site 12.
