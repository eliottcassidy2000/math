# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 10, 2026)

### 2.1 THM-459 & THM-460: The Erdős 592 Breakthroughs
A major structural and algebraic unification across the Erdős 592 and 64 problems (mac-mini-2026-06-09-S2):
- **THM-459 (R(2,2)=5 Lemma Layer):** Complete structural proof of $R(2,2)=5$ on the 5-ary height-2 tree. Key innovations include:
    - **L3 (Doubly-Dark Clique Lemma):** For triangle-free graphs, any $R_a$-independent pair $\{x_1, x_2\}$ forces the doubly-dark set $D = [5] \setminus (N_{a'}(x_1) \cup N_{a'}(x_2))$ to be a clique in $R_{a'}$.
    - **Machine Closure:** Certified by Glucose3 UNSAT (532 CEGAR clauses).
    - **Invariant-Free Clash:** Confirmed that the Invariant Cutoff = Free Cutoff = 5 at $n=2$, marking the binding obstruction as a Schur condition for graded relations.
- **THM-460 (Chang-Tower Miniature):** A full-type characterization for ordinals below $\omega^{\omega^m}$ (Erdős 592).
    - **Full-Type = Stacked Towers:** Proved that a subset $X \subseteq \omega^{\omega^m}$ has full order type iff it contains a full cofinal tower of shapes.
    - **König Bridge:** Established the finite-to-infinite lift (Theorem C1) for triangle-free witnesses.
    - **Shape Grammar:** A recursive binary grammar yielding pairs at finite scales and M-towers at limit scales.
- **Dyadic Sufficiency (HYP-2363..2366):** Breakthrough linking Erdős 64 to 592. Restricting relations to depend only on the 2-adic valuation $v_2(gap)$ preserves the cutoff of 5 at $n=2$. This identifies the "2-adic seam" across the repository's Ramsey results.

### 2.2 THM-456 (The Blowup Spectrum Law & Turán-Corridor Reframe)
A major structural breakthrough in the Erdős–Gyárfás (E-G) cycle problem (kind-pasteur-2026-06-09-S2):
- **THM-456 (Blowup Spectrum Law):** For any graph $G$, the cycle spectrum of its twin blowup $G[K_2]$ is the gap-free interval $[3, 2s(G)]$, where $s(G)$ is the maximum order of a subgraph carrying edge multiplicities {1,2} with multi-degrees in {2,4}.
- **Interval Lemma:** Proved that blowups $G[K_2]$ can never be E-G counterexamples, as they always contain a power of 2 in every dyadic window.
- **Turán-Corridor Reframe:** Counterexamples must be twin-free and live within $O(1)$ edges of the $C_4$-Turán-extremal bound ($ex(n; C_4)$).
- **Corridor Closure:** The corridor is closed for $n \le 9$ (density forces $C_4$).
- **n=10..12 Census:** All 71 $C_4$-free $\delta \ge 3$ graphs in this range were killed by a forced $C_8$.

### 2.3 MISTAKE-068: The Blowup Circumference Over-Correction
A critical correction regarding cycle lengths in blowups (kind-pasteur-2026-06-09-S2):
- **The Error:** An earlier assumption that the circumference of $G[K_2]$ was exactly $2p(G)$ (twice the longest path).
- **The Correction:** While true for 939/995 cases, 56 "beaters" (like the Net graph) achieve $c > 2p$ via "sun walks" (stutter-walks). The true bound is $2s(G)$.

### 2.4 THM-442 (Pfaffian Translator & IE Tiling Split)
A critical distinction has been established between the additive and multiplicative structures of the repository (S707, THM-442, HYP-2312):
- **IE Tiling = Third Finite Difference:** The 7-piece inclusion-exclusion (IE) decomposition is identified as the third finite difference ($\Delta^3$) of quadratic cell-counts. It computes **cell-affine** (additive) invariants but fails for the Hamiltonian cycle count ($H$).
- **The Additive-Multiplicative Split:** $H$ is confirmed to be **multiplicative (modular)**, not additive. No additive recursion for max-H ($A038375$) exists.
- **Pfaffian Translator:** $Pf(S)^2 = det(S)$ acts as a polynomial-time translator ($n \to n-2$) for cycle covers. 
- **Max-H ⟺ Minimal |Pf|=1:** A new bridge links the \#P-hard Hamiltonian problem to the poly-time Pfaffian: verified for $n=4, 6$, the **max-H tournament has minimal $|Pf|=1$**. This reduces the search for extremal tournaments to the $det(I+2A)=1$ constraint.

### 2.5 THM-441 (The Convolution-Correlation-Adjoint Duality)
A fundamental unification of the repository's core mathematical objects (S706, THM-441, HYP-2311):
- **Core Identity:** Cross-correlation is the **adjoint** of convolution ($\langle h*a, b \rangle = \langle a, h \star b \rangle$).
- **Clock-Shell Unification:** The **Shell face** (sums mod $2n-1$) is a convolution, while the **Clock face** (differences mod $n$) is a cross-correlation. They are adjoints related by the antipodal involution $\sigma$.
- **Converse Duality:** The tournament converse $T \mapsto T^*$ is the adjoint of the circulant adjacency-convolution. **Worry-sets** are precisely the **self-adjoint** operators.

## 3. Laboratory Status

### 3.1 Chang-Tower Lab (mac-mini-S2)
- **Status:** Running.
- **Cross-Validation:** PASS (Structural tower finder vs. brute force, 0 disagreements).
- **Current Sweep:** SAT confirmed for $m=1, M=2$ at $(2,3), (3,2), (2,4)$.
- **Frontier:** $s=3, C=3, M=2$ currently grinding (UNSAT-convergence-shaped). This is the first computed "Chang number".

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: THE CHANG TOWER PROBE (NEW - S709)
**Objective:** Probe the $\omega^{\omega^3}$ ($1000 open case) via the THM-460 miniature.
1. **First Chang Number:** Harvest the $(3,3)$ UNSAT result to establish the first finite cutoff for Chang's theorem ($m=1$).
2. **m=2 Cutoffs:** Implement the B3 general-shape recursive enumerator to find the first cutoffs for Schipperus's theorem ($m=2$).
3. **m=3 Persistence:** Begin the $M=3$ probe at $n=3$ to check for persistence of SAT witnesses.

### Steering Task 2: THE DYADIC SEAM (NEW - S710)
**Objective:** Formalize the 2-adic link between Erdős 592 and Erdős 64.
1. **Closed-Form Dyadic Rule:** Extract and prove the closed-form dyadic rule for the $R(2,2)=5$ witness.
2. **Strong Specker Construction:** Use the dyadic rule + THM-453 D1 to attempt a constructive strong Specker witness.

### Steering Task 3: THE ERDŐS-GYÁRFÁS TURÁN PROBE (S708)
**Objective:** Leverage the THM-456 reframe to locate potential counterexamples.
1. **Corridor Census Extension:** Extend the search for $C_4$-free $\delta \ge 3$ graphs to $n=13..16$. Check if forced $C_8$ or $C_{16}$ persists as the primary killer.
2. **Twin-Free Forcing:** Verify if the "twin-free" constraint significantly reduces the search space for $ex(n; C_4)$ candidates.

### Steering Task 4: THE PFAFFIAN-MINIMALITY PROBE (S707)
**Objective:** Use the Pfaffian translator to constrain the search for max-H tournaments.
1. **Pfaffian Census:** Verify the $|Pf|=1$ property for max-H tournaments at $n=8, 10$.
2. **Multiplicative Bridge:** Investigate if the $H^2 - Pf^2 = 8Q$ identity can be sharpened to a recursive modular bound.

### Steering Task 5: THE u(22) SETTLEMENT (S705)
**Objective:** Confirm if $u(22)=60$ or 61.
1. **delta=5 Census:** Systematically generate all 56-edge 21-cores and test for unit-cocyclic degree-5 extensions.

### Steering Task 6: THE TWO-LAW TRANSIT (High Priority)
**Objective:** Derive the explicit Möbius map between the Factorial Law (diagonal, $\kappa_n = n!$) and the Two-Point Law (signed row sum).
