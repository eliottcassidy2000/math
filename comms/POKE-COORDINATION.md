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

## codex-S136 -- C=27 Shell-Transfer Spectrum (checkpoint)

Formalized the **C=27 Shell-Transfer Spectrum**, specialized the $C=27$ shell carrier by treating it as a relation space rather than a scalar tag (commit `a55d7e5e`). This checkpoint stabilizes the "marked transfer" protocol for identifying candidate tight rows near the AP floor.

### 1. The C=27 Shell Relation Space
Redefined the $C=27$ shell (the $p=2$ summand branch) as a space of shell-pair relations $P_a = \{a, 27-a\}$ for $1 \le a \le 13$.
- **Marked Transfers:** Identified that Goddyn-Wong (GW) and near-misses are not random noise but specific transfers across this shell space:
    - **GW:** Transfers a hole at $12$ ($gcd(3)$) to a double at $3$ ($gcd(3)$).
    - **12->36 Near-Miss:** Transfers a hole at $12$ ($gcd(3)$) to a double at $9$ ($gcd(9)$), pushing the defect into a deeper 3-adic layer.
    - **Petals (10->20, 13->26):** Characterized by **unit-visible holes** ($gcd(1)$), mapping them to the second-gap rigidity branch.

### 2. Bounded Frontier Analysis
Audit of the single-replacement atlas through replacement $140$ confirms a remarkably sparse tight/near-tight set:
- **M <= 3/41:** AP, GW, 12->36.
- **M <= 2/27:** AP, GW, 12->36, 10->20, 13->26.
This sparsity reinforces the thesis that low-gap non-AP/GW atoms must exhibit specific shell-transfer signatures.

### 3. Role in Proof Tree (Summand p=2 Branch)
The shell-transfer spectrum provides the structural link between the $p=2$ branch and terminal rigidity:
- **Unit-Visible Holes:** Route to the petal rigidity and two-block shell collapse machinery.
- **Non-Unit Holes (gcd-3, gcd-9):** Route to the $p \ge 3$ $K_{3,3}$ wall or the tournament-conflict state lift via the product-incidence packet.

### 4. Relational Guardrail
Reinforced that the shell quotient is a **carrier**, not a standalone invariant. The transfer labels must be integrated with the exact Farey binding scale ($M=p/q$) to ensure magnitude-awareness and prevent magnitude-blindness leakage.

### 5. Net Impact
This checkpoint stabilizes the project's ability to "label the rooms" in the low-gap witness space. By formalizing shell-transfer signatures, the cluster has moved the $LRC(14)$ proof from a broad search to a targeted audit of finite, marked transitions across the $C=27$ unit/nonunit quotient.
