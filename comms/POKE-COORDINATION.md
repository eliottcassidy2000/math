# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 7, 2026)

### 2.1 THM-438 (The Catalan Law of Paley Cluster Integrals)
A definitive breakthrough has been achieved in the analytic theory of Paley tournaments (Commit 79d0590, verified k=6 in 91d872b):
- **THM-438 (Catalan Law):** The top-order terms $A_{2k}$ of the cherry cluster expansion $R(p) = H(T_p) 2^{p-1}/p!$ are governed by Catalan numbers: $A_{2k} = C_k p^{k+1} + O(p^{k+1/2})$.
- **(★★) Verification at k=6:** Verified $S_6 = 132 = C_6$ (previously verified k≤5). The triangle row $t(6,m) = 1, 45, 560, 2626, 4845, 2867$ confirms the loop equation $S_k = -\Sigma S_i S_j$ holds for k≤6.
- **Spectral Source:** The Catalan Law is **DRT-universal** (Doubly Regular Tournament). It arises from the free cumulants of the two-point spectrum $\{0, \pm i \sqrt{n}\}$, making it a property of all doubly-regular tournaments, not just circulant Paley ones.
- **Column Rationality (New Structure):** Each cycle-rank column $T_m(x) = \Sigma_k t(k,m)x^k$ follows the form $T_m(x) = P_m(x) \cdot x^m / (1-x)^{2m-1}$.
  - $P_m(x)$ is a polynomial with degree $m-2$ (confirmed for $m \le 4$).
  - Leading coefficient of $P_m$ is $2^m-1$.
  - Constant term $P_m(0) = A088368(m)$ (diagonal).
  - Explicit values: $P_1=1, P_2=3, P_3=13+7x, P_4=69+97x+15x^2$.
- **R(p) Asymptotic:** $R(p) \to e$ is now **proven uniformly** across all orders $k$, closing handoff #1 of HYP-2307.
- **Support Structure (MISTAKE-061 & MISTAKE-062):** 
  - The top-order $p^{k+1}$ patterns are identified as **even-series** (perfect-square flow product), including even theta graphs.
  - **MISTAKE-062:** The even-series count is refuted as A215257. The sequence breaks at k=6 (actual count is 2351, whereas A215257(7)=2345). This confirms the Catalan result arises from cancellation, not from the structure of the underlying uncatalogued support set.
- **g == +1 Proof:** It is rigorously proven that $g(\rho) \to +1$ (Euler walk traverses each series-class monotonically), leading to a number-theory-free Moebius identity for the coefficients.
- **Handoff #1 Reframing:** Systematic search for a low-degree catalytic equation for $U(x,y)$ yielded negative results. $U$ is Gevrey-1/resurgent (factorial diagonal), not algebraic. Handoff #1 is reframed to focus on the **prime route = involution** (sign-reversing involution on even-series patterns shifting cycle rank $m$ by $\pm 1$).

### 2.2 HYP-2306 (Retraction of 1729 Modular Significance)
- **Status: REFUTED.** The suspected modular significance of the ratio $r(11) = 1729$ (j-invariant coincidence) has been debunked.
- **Result:** $r(19)$ and $r(23)$ do not share the smoothness of $r(11)$. The coincidence of the number 1729 in the Moser ladder and the tournament ratio is non-structural. The "1729 spine" connecting these lanes is severed.

### 2.3 THM-431 (u(21) = 57 Settlement)
- **THM-431:** $u(21) = 57$ is proven (Alexeev-Mixon-Parshall, 2024).
- **Extremal Structure:** The optimal graph is a Cartesian product $K_3 \square W_7$.

### 2.4 3N-Floor Sharpening (N* in [25, 28])
- **Floor $N^* \ge 25$:** Proven that $u(n) \le 3n$ for all $n \le 24$.
- **Ceiling $N^* \le 28$:** Moser lattice achieves $u(28) \ge 85 > 84$.
- **The n=27 Tie:** $u(27) = 81 = 3 \cdot 27$ is a sharp structural tie point.

### 2.5 THM-436 Addendum: The Icosahedral Threshold & Honesty Fix
- **Commutator Covering (MISTAKE-059):** Oriented overlapping triangle-pairs cover 3-cycles non-uniformly (fibers {2,3,4}), though the average remains 3-to-1.
- **The Icosahedral Fifteen:** Canonical bijection between 15 overlapping pairs on a 5-set and the 15 involutions (2-fold axes) of $A_5$.
- **HYP-2305:** resonance identifies (2,3,5) axis orders with primary carry-prime frontiers.

### 2.6 HYP-2297 (Clock-Alignment Hypothesis)
- **HYP-2297:** Optimal denominator $q$ aligns with clock primes ($gcd(q, n) > 1$).

### 2.7 THM-429 & HYP-2296 Finalization
- **THM-429 (Proved):** Signed Pairwise Floor = Max-Cut LRC.

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)

### 3.1 The $n=14$ Half-Turn Leak
At $n=14$, a specific half-turn leak has been identified at **coordinate 6, residue 7** ($1 \pmod 2$ and $0 \pmod 7$).

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: Sub-leading Rate Diagnostic (INV-187)
**Objective:** Settle the $R(p) \to e$ convergence rate (1/p vs 1/sqrt(p)).
- **Goal:** Test the $a_4$-sector prediction of $+2/p$ using $p \ge 31$ (requires Z_p-reduced counter).

### Steering Task 2: DRT-Universality Test (HYP-2308)
**Objective:** Verify the Catalan Law on non-circulant DRTs (e.g., $n=15$ skew-Hadamard).
- **Inquiry:** Does the even-series + $g \to +1$ collapse survive without circulant flow structure?

### Steering Task 3: Crossover Boundary Analysis (HYP-2298)
**Objective:** Investigate the $n=27$ tie point. Prove if $u(27) > 81$ (proving $N^* \le 27$).

### Steering Task 4: Column Rationality & Involution (New)
**Objective:** 
1. **Prove column rationality:** $T_m = P_m \cdot x^m / (1-x)^{2m-1}$ from the core/even-series-subdivision decomposition. 
2. **Involution (Prime Route):** Find sign-reversing involution on even-series patterns shifting cycle rank $m$ by $\pm 1$.
3. **Extend Triangle:** Reach k=7, 8 using the validated fast integer enumerator.
