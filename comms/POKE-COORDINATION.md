# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. The Fiber-Projection & Torsion-Leak Theory

### 2.1 The $n=14$ Half-Turn Leak
At $n=14$, a specific half-turn leak has been identified at **coordinate 6, residue 7**. 
- **Algebraic Alignment:** This point sits exactly at $1 \pmod 2$ and $0 \pmod 7$. 
- **Structural Interpretation:** This represents a 2-torsion subgroup projecting to zero mod 7. This specific alignment suggests a structural vulnerability in the product-sieve cover where the symmetry of the subgroup fails to provide sufficient coverage density.

### 2.2 Comparative Analysis: $n=15$
For $n=15$, we observe order-3 leaks at residues 5 and 10.
- **Structural Interpretation:** These represent a 3-torsion subgroup projecting to zero mod 5. 

### 2.3 General Hypothesis
For a composite $n = p \cdot q$, we hypothesize that the leakage from the product-sieve cover is algebraically constrained to the torsion subgroups of the fiber bundle projections. 
**Corollary:** If the leaks across different divisor fibers are mutually exclusive, a full open cover cannot exist, implying an inherent incompleteness in the covering strategy for composite moduli under certain constraints.

## 3. Concrete Computational Prompts for Cluster Agents

### Task 1: Boundary Mapping ($n=14$)
**Objective:** Map the boundary of the 56 unexcluded cells for the $n=14$ leak. 
- **Verification Requirement:** Determine if these boundaries align with the boundaries of the 7-runner zonotope projection. 
- **Output:** Provide a list of boundary vertices and their distance from the ideal projection boundary.

### Task 2: Rational Interval Probes
**Objective:** Conduct rational interval probes around the 2-torsion and 3-torsion fibers for both $n=14$ and $n=15$.
- **Focus:** Analyze the density of the coverage in the immediate $\epsilon$-neighborhood of these fibers.
- **Output:** Report any gaps found and the maximum width of uncovered intervals.

### Task 3: Protection-Chain Depth Bounding
**Objective:** Attempt to bound the 'protection-chain depth' of the $n$-cover using the known depths of its factor covers.
- **Inquiry:** Can the depth of the composite cover be expressed as a linear combination or a more complex function of the depths of the $p$ and $q$ covers?
- **Output:** Formulate a preliminary bounding equation based on empirical data from previous tasks.
