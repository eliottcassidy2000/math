# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 8, 2026)

### 2.1 THM-442 (Pfaffian Translator & IE Tiling Split)
A critical distinction has been established between the additive and multiplicative structures of the repository (S707, THM-442, HYP-2312):
- **IE Tiling = Third Finite Difference:** The 7-piece inclusion-exclusion (IE) decomposition is identified as the third finite difference ($\Delta^3$) of quadratic cell-counts. It computes **cell-affine** (additive) invariants but fails for the Hamiltonian cycle count ($H$).
- **The Additive-Multiplicative Split:** $H$ is confirmed to be **multiplicative (modular)**, not additive. No additive recursion for max-H ($A038375$) exists.
- **Pfaffian Translator:** $Pf(S)^2 = det(S)$ acts as a polynomial-time translator ($n \to n-2$) for cycle covers. 
- **Max-H ⟺ Minimal |Pf|=1:** A new bridge links the \#P-hard Hamiltonian problem to the poly-time Pfaffian: verified for $n=4, 6$, the **max-H tournament has minimal $|Pf|=1$**. This reduces the search for extremal tournaments to the $det(I+2A)=1$ constraint.

### 2.2 THM-441 (The Convolution-Correlation-Adjoint Duality)
A fundamental unification of the repository's core mathematical objects (S706, THM-441, HYP-2311):
- **Core Identity:** Cross-correlation is the **adjoint** of convolution ($\langle h*a, b \rangle = \langle a, h \star b \rangle$).
- **Clock-Shell Unification:** The **Shell face** (sums mod $2n-1$) is a convolution, while the **Clock face** (differences mod $n$) is a cross-correlation. They are adjoints related by the antipodal involution $\sigma$.
- **Converse Duality:** The tournament converse $T \mapsto T^*$ is the adjoint of the circulant adjacency-convolution. **Worry-sets** are precisely the **self-adjoint** operators.

### 2.3 THM-438 (The Catalan Law of Paley Cluster Integrals)
Definitive breakthrough in the analytic theory of Paley tournaments:
- **THM-438 (Catalan Law):** Top-order terms $A_{2k}$ are governed by Catalan numbers: $A_{2k} = C_k p^{k+1} + O(p^{k+1/2})$.
- **The Watson Bridge:** Unified the wild diagonal $A088368$ and the density $\rho(x)$ as faces of a single resurgent series with shared coefficients $b = 1, 2, 10, 178/3, \dots$.
- **Bercovici-Pata Partnership:** The wild ($x=0$) and tame ($x=\infty$) endpoints are partners sharing the cumulant signature $\kappa_n = n!$.

### 2.4 THM-440 & HYP-2310 (u(22) Extension Reduction)
- **THM-440:** A 61-edge UDG on 22 vertices must be a unit-cocyclic extension of a 21-vertex core.
- **Moser Ring Caps at 60:** Within $M_L$, $u(22)=60$. Any $u(22)=61$ requires a $\delta=5$ route or a field escape to $\mathbb{Q}(\sqrt{-3}, \sqrt{-d})$.

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)
At $n=14$, identified half-turn leak at **coordinate 6, residue 7** ($1 \pmod 2$ and $0 \pmod 7$).

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: THE PFAFFIAN-MINIMALITY PROBE (NEW - S707)
**Objective:** Use the Pfaffian translator to constrain the search for max-H tournaments.
1. **Pfaffian Census:** Verify the $|Pf|=1$ property for max-H tournaments at $n=8, 10$.
2. **IE Tiling vs. H:** Formally document the failure of additive IE recursions for $A038375$ to prevent further dead-end additive searches.
3. **Multiplicative Bridge:** Investigate if the $H^2 - Pf^2 = 8Q$ identity can be sharpened to a recursive modular bound.

### Steering Task 2: THE ADJOINT SYMMETRY PROBE (S706)
**Objective:** Map every repo object as a convolution operator and identify its adjoint.
1. **Self-Adjoint Census:** Confirm if all "worry-set" configurations are exactly the self-adjoint operators in their respective cyclotomic groups.
2. **Autocorrelation Moments:** Re-examine THM-406 covering-depth moments as autocorrelation integrals.

### Steering Task 3: THE u(22) SETTLEMENT (S705)
**Objective:** Confirm if $u(22)=60$ or 61.
1. **delta=5 Census:** Systematically generate all 56-edge 21-cores and test for unit-cocyclic degree-5 extensions.
2. **Field Escape Check:** Probe for 22-point configurations in $\mathbb{Q}(\sqrt{-3}, \sqrt{-d})$ for new Heegner d.

### Steering Task 4: THE TWO-LAW TRANSIT (High Priority)
**Objective:** Derive the explicit Möbius map between the Factorial Law (diagonal, $\kappa_n = n!$) and the Two-Point Law (signed row sum) to resolve the off-diagonal $P_m$ polynomials.

### Steering Task 5: CYCLOTOMIC DEPTH CENSUS
**Objective:** Run a systematic $q$-census for $n=9$ to $n=14$ to verify if $gcd(q, n) > 1$, testing if $R(14)$ has depth equal to the $3^3$ tower.

### Steering Task 6: TAME HANDOFF PROOFS
**Objective:** Close the remaining gaps by showing $(k)_m | t(k, m)$ and deriving the exact $g_m(-1)$ values.

### Steering Task 7: HAAGERUP-MÖLLER ANALYSIS
**Objective:** Rigorously prove the $1/\sqrt{x}$ edge and the $e \cdot e^{-x}$ tail of the free compound Poisson distribution.
