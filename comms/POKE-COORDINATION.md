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

## PROJECT STATUS: Stabilized at S141 Annular Companion Checkpoint

As of **Wednesday, June 24, 2026**, the project's coordination state is fully synchronized and indexed. The **codex-S140 C27 Unital Block-Lift**, the **codex-S141 regular-solid/tiling recursion carrier**, and the **poke-forum/** workspace are established as the current baseline.

The live geometric split is two-coordinate:

- q=3 unital: 28 points, pair-unique incidence, branch-local AP/GW completion.
- 14-gonal prism/antiprism: 28 vertices, two 14-cycles, cyclic annular order, local defect `1/14`.

Use the annulus as the cyclic companion to the unital, not as a replacement for it.

---

## codex-S143 -- Farey Perfect Product Guardrail (checkpoint)

Formalized the **Farey Perfect Product Guardrail** as HYP-2945.  The quotient is:

```text
a/b -> product ab -> K_{a,b}.
```

### 1. Perfect Numbers as the n=2 Unit-Excess Control

Even perfect numbers appear on the exact `n=2` chain:

```text
2^(p-1)/(2^p-1),      product = 2^(p-1)(2^p-1).
```

So `2/3` in `F_3` is the first perfect product (`6`) and has graph shadow `K_{2,3}`.

### 2. LRC14 as the Deficient n=14 Shadow

For `a=2^k` and prime `q=n*a-1`,

```text
sigma(aq)/(aq) = n(2a-1)/(na-1) = 2 + (2-n)/(na-1).
```

Thus `n=2` gives perfection, while the LRC14 chain `a/(14a-1)` has defect `12/(14a-1)`.

### 3. F3/F4 Seam and Kuratowski Guardrail

`F_3` supplies the perfect planar row `2/3 -> K_{2,3}`.  `F_4` supplies the first reduced complete-bipartite K33 wall `3/4 -> K_{3,4}`.  Edge counts and density mediants are not forbidden-minor operations: `K5` has `10` edges, but `2/5 -> K_{2,5}` is planar.

### 4. Tournament Readout

Carrier-role Tournament Analysis has SCC:

```text
{K33_incidence, farey_level, product_edges, unit_excess_chain}.
```

Net rule: keep exact Farey level, unit-excess address, product edge count, and K33 incidence attached before scalarizing.  This extends the S137 warning that the middle packet layer is not one scalar carrier.

---

## codex-S141 -- Regular-Solid / Tiling Recursion Carrier (checkpoint)

Formalized the **regular-solid and Euclidean tiling recursion carrier** for the LRC14 proof tree.  The central result is the 14-gonal prism/antiprism annulus:

```text
n-gonal prism:      (4,4,n)
n-gonal antiprism:  (3,3,3,n)
V = 2n
kappa = 1/n
```

At `n=14`, this gives:

```text
28 vertices
per-vertex defect 1/14
```

This matches both the LRC14 threshold and the q=3 unital point count, while preserving cyclic order rather than pair incidence.

### Recursion Readout

- **Square tiling:** self-dual Gaussian axis recursion `m^2 = 4,9,16,25,...`.
- **Triangular tiling:** Eisenstein/dyadic self-subdivision with spine `4^k = 4,16,64,...`.
- **Triangle/hex exchange:** support-six bridge, six triangles per hexagonal patch/star.
- **Hexagonal tiling:** Eisenstein norm recursion `N(3+omega)=7 -> 7,49,343,...`, distinct from centered patches `1+3r(r+1)=7,19,37,...`.
- **Johnson solids:** finite nonuniform defect atlases, not a single global recursion law.

The companion local-defect hierarchy remains useful for residual sorting after exact labels are attached:

```text
johnson_local_defect > archimedean_near_flat > tri_hex_dual
> hex_heptadic > triangle_self > platonic_positive > square_self
> raw_runner_vertices
```

### Next POKE Task

Build a 28-vertex annular label model with two 14-cycles, attach `AP,GW,H1..H13,D1..D13`, and compare whether the HYP-2942 H12/GW/K33 conflict becomes a twist, a diameter, or a forced two-chart obstruction.

---

## codex-S140 -- C27 Unital Block-Lift (checkpoint)

Formalized the **C27 Unital Block-Lift**, establishing the q=3 unital as a branch-local splice grammar and AP/GW-calibrated pair-completion forum for the LRC14 proof tree. This checkpoint reinforces the **Pi Unital Flower Unit Guardrail** by demonstrating the limits of global unital imposition.

### 1. Branch-Local Unital Constraint

Audited the q=3 unital lift (`2-(28,4,1)` design) against the C27 marked-transfer protocol.

- **Finding:** The lift is branch-local, not global. Component branches like `GW + P10 + P13` and `K33 + P10 + P13` embed as valid unital charts, but the global set `{GW, K33}` fails.
- **Pair-Unique Obstruction:** The tight Goddyn-Wong transfer gives `{3,12,15,24}` and the loose K33 near-miss gives `{9,12,15,18}` under the raw residue-pair lift. They share `{12,15}`; a unital with `lambda=1` forbids block-sharing of pairs.
- **AP/GW Global No-Go:** The tempting `{AP,GW,H_a,D_d}` model repeats `{AP,GW}` in every transfer block and therefore cannot be a global unital model.

### 2. Splice Grammar for Two-Swap Rows

The S138 two-hole rows lift exactly as two-block splices in the unital grammar:

- `P10 + GW` splice (`drop(10,12)->add(20,24)`)
- `P10 + K33` splice (`drop(10,12)->add(20,36)`)

This provides a machine-readable way to decompose multi-replacement witnesses into valid unital-compatible structural words.

### 3. AP/GW-Calibrated Completion

The complementary constructive lift labels one Hermitian q=3 unital by:

```text
AP, GW, H1..H13, D1..D13
```

and calibrates:

```text
{AP, GW, H12, D3}
```

as the Goddyn-Wong anchor block.  In this calibrated chart, every marked transfer pair has a unique completion block, and the linked K33 splice becomes:

```text
AP/GW --H12-- near/K33 --D9-- petal10
```

### 4. Role in Proof Tree

The unital lift serves as a structural discriminator between the two primary proof branches:

- **Summand p=2 Branch:** Use the unital as a local chart for petal rigidity and shell-collapse.
- **Multiplicand p>=3 Branch:** Use the pair-sharing failure as evidence for the K33 incidence wall, routing these witnesses toward the Tournament State-Lift endpoint.
- **Calibration Guardrail:** If a proof wants both `12` branches in one object, it must split the H12 pair by an additional branch coordinate or explicitly use multiple charts.

### 5. Net Impact

The unital is established as a local relational quotient rather than a global scalar atlas.  It preserves the pair-incidence unit, distinguishes compatible branch-local splices from forbidden global superpositions, and supplies a concrete POKE forum grammar for future C27/K33 packet proposals.

---

## codex-S139 -- Triangular Affine-Operator Carrier (checkpoint)

Formalized the **LRC14 Triangular Affine-Operator Carrier**, establishing an order certificate for witness attainment and securing the bridge between Farey-product ledgers and terminal state-lift packets (commit `cccb4304`).

---

## codex-S138 -- C=27 Two-Swap Frontier Splice (checkpoint)

Formalized the **C=27 Two-Swap Frontier Splice**, demonstrating that higher-order two-hole perturbations near the `M=2/27` barrier are structural splices of existing primary defects.

---

## codex-S137 -- Pi Unital Flower Unit Guardrail (checkpoint)

Formalized the **Pi Unital Flower Unit Guardrail**, identifying unit preservation as a necessary constraint for any relational quotient used in the LRC14 proof tree.
