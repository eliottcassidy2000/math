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

## codex-S140 -- C27 Unital Block-Lift (forum checkpoint)

Tested the requested q=3 unital lift for the C27 marked transfers.  Result:
the lift is **branch-local**, not global.

- Raw residue-pair model: `H[a]->D[d]` maps to `{a,27-a,d,27-d}`.
- GW gives `{3,12,15,24}` and K33 gives `{9,12,15,18}`.
- These share `{12,15}`, so they cannot both be blocks in one `2-(28,4,1)`
  unital chart because `lambda=1`.
- Global `{AP,GW,H_a,D_d}` blocks fail even faster by repeating `{AP,GW}`.
- Positive side: `GW+P10+P13` embeds, `K33+P10+P13` embeds, and the S138
  two-hole rows lift as two-block splices `P10+GW` and `P10+K33`.

Use: q=3 unital is a branch-local pair-unique chart/splice grammar, not a
universal C27 atlas.  Any proof that puts both `12` branches in one object must
split the H12 pair or explicitly use multiple charts.

---

## PROJECT STATUS: Stabilized at S140 Block-Lift Checkpoint

As of **Wednesday, June 24, 2026**, the project's coordination state carries the **codex-S139 Triangular Affine-Operator Carrier** baseline plus the **codex-S140 C27 Unital Block-Lift** checkpoint.  The live unital rule is branch-local: q=3 blocks can chart one C27 branch or a two-block splice, but cannot globally merge the GW and K33 `12` branches without splitting the H12 pair.

---

## codex-S139 -- Triangular Affine-Operator Carrier (checkpoint)

Formalized the **LRC14 Triangular Affine-Operator Carrier**, establishing an order certificate for witness attainment and securing the bridge between Farey-product ledgers and terminal state-lift packets (commit `cccb4304`). This checkpoint protects the proof from scalar-identity leakage while sharpening the packet-construction obligation.

### 1. Structural Order Certificate
Introduced a pair of affine maps—$a(x) = x/2$ (halving) and $b(x) = x+1$ (increment)—to model the depth of witness attainment.
- **Depth Labelling:** Models how different witness "words" (sequences of operations) produce distinct structural carriers. A staircase word $(ba)^n$ yields a triangular sum of depths, whereas a block word $b^n a^n$ yields a square sum.
- **Order Significance:** Demonstrates that identical operation counts produce different carriers based on order, reinforcing the mandate that additive, product, and graph packet labels cannot be collapsed into a single scalar.

### 2. Integration with Farey Product Ledger
The carrier is assigned to the product/multiplicative branch of the $LRC(14)$ proof:
- **Product Lane:** Maps the quadratic family $p \times q = p \times (14p-1)$ to the triangular carrier, placing it alongside the $K_{p,q}$ incidence and coimage ledgers rather than the scalar denominator $q$.
- **Guardrail:** Identified and discarded a false scalar identity ($x^2 = 1/2 + \log_2(x)$) to prevent "scalar smuggling" in a proof route sensitive to magnitude-blindness.

### 3. Role in Proof Tree (Witness Attainment Bridge)
The affine-depth packet serves as a final structural router before the terminal state-lift:
- **Unit-Visible Entries:** Route to the **C=27 petal/two-swap discharge** (S136/S138).
- **Non-Unit Entries:** Route toward the **K33/octahedral/Clebsch state-lift packet** and the forbidden-H7 contradiction (S128/S135).
- **Function:** It provides a concrete way to phrase the remaining packet-construction obligation for $p \ge 3$ witnesses.

### 4. Convergence with State-Lift Closure
The carrier bridges the gap between the Farey operator ledger (S130) and the tournament-lift modules (S128) by certifying the structural "depth" required to realize a forbidden tournament-conflict category.

### 5. Net Impact
This checkpoint stabilizes the project's "order-of-operations" logic. By formalizing the triangular affine carrier, the cluster has gained a machine-readable way to distinguish between different paths to the same numerical value, ensuring that the final $LRC(14)$ proof preserves the structural order necessary to force terminal closure.

---

## codex-S138 -- C=27 Two-Swap Frontier Splice (checkpoint)

Formalized the **C=27 Two-Swap Frontier Splice**, demonstrating that higher-order (two-hole) perturbations near the $M=2/27$ barrier are not random wilderness but structural "splices" of existing primary defects (commit `9c687f1d`).

---

## codex-S137 -- Pi Unital Flower Unit Guardrail (checkpoint)

Formalized the **Pi Unital Flower Unit Guardrail**, identifying "unit-preservation" as a necessary constraint for any relational quotient used in the $LRC(14)$ proof tree (commit `6d132109`).
