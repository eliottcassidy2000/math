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

## PROJECT STATUS: Stabilized at S136 Baseline

As of **Wednesday, June 24, 2026**, the project's coordination state is fully synchronized and indexed. The **codex-S136** documentation—including the **C=27 Shell-Transfer Spectrum** and the **Binary Relational Exploration** mandate—is established as the project's current baseline.

---

## codex-S138 -- C=27 Two-Swap Frontier Splice (checkpoint)

Formalized the **C=27 Two-Swap Frontier Splice**, demonstrating that higher-order (two-hole) perturbations near the $M=2/27$ barrier are not random wilderness but structural "splices" of existing primary defects (commit `9c687f1d`). This checkpoint stabilizes the inductive step for the $LRC(14)$ terminal closure.

### 1. The Frontier Splice
Audited the two-swap perturbation bank (replacements $\le 40$) to test the stability of the $C=27$ marked-transfer protocol.
- **Stability:** The $M \le 3/41$ frontier remains stable, containing only AP, GW ($12 \to 24$), and the near-miss ($12 \to 36$).
- **Splice Rows:** Identified two new $M=2/27$ rows:
    - `drop(10,12)->add(20,24)`
    - `drop(10,12)->add(20,36)`
- **Structural Mapping:** These are precisely splices of a **unit defect** (the $10 \to 20$ petal) with a **non-unit 12-branch defect** (either the GW-tight $12 \to 24$ or the loose $12 \to 36$).

### 2. Relational Identification
The findings reinforce that the shell-transfer quotient correctly classifies complex witnesses:
- **Relational Address:** A two-swap witness is labelled by its component transfers: `[unit hole -> unit double] + [nonunit hole -> nonunit double]`.
- **Simplification:** Rather than treating two-swaps as a new complexity class, the proof can dispatch them by showing their component unit defects already force the second-gap witness.

### 3. Role in Proof Tree (Inductive Closure)
The splice behavior provides a path for closing the witness search:
1. **Filter:** Apply $q \ge 14$ threshold.
2. **Classify:** Identify unit-visible $C=27$ holes.
3. **Dispatch:** Prove that the only two-hole unit petal is $10 \to 20$, and its companions are restricted to the known $12$-branch.
4. **Link:** Connect the $12 \to 36$ component to the $K_{3,3}$/state-lift loose category.

### 4. Integration with Binary Relational Mandate
This checkpoint applies the "Structural Imposition" principle by using the shell-transfer relation to decompose complex perturbations into leverageable, labelled parts. It prevents the "cloud of complexity" that often blocks terminal audits of multiple-replacement sets.

### 5. Net Impact
This checkpoint stabilizes the project's ability to handle multi-replacement perturbations. By formalizing the "frontier splice," the cluster has turned a potentially unbounded search into a structured audit of defect combinations, ensuring that the final $LRC(14)$ closure remains inductive and modular.

---

## codex-S137 -- Pi Unital Flower Unit Guardrail (checkpoint)

Formalized the **Pi Unital Flower Unit Guardrail**, identifying "unit-preservation" as a necessary constraint for any relational quotient used in the $LRC(14)$ proof tree (commit `6d132109`). This checkpoint protects the unital unit-distance geometric branch from leakage into non-unital states.

---

## codex-S136 -- C=27 Shell-Transfer Spectrum (checkpoint)

Formalized the **C=27 Shell-Transfer Spectrum**, specialized the $C=27$ shell carrier by treating it as a relation space rather than a scalar tag (commit `a55d7e5e`). This checkpoint stabilizes the "marked transfer" protocol for identifying candidate tight rows near the AP floor.
