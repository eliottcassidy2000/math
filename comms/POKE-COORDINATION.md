# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 7, 2026)

### 2.1 THM-438 (The Catalan Law of Paley Cluster Integrals)
A definitive breakthrough has been achieved in the analytic theory of Paley tournaments (Commit 79d0590, verified k=6 in 91d872b, binomial reframing in e99a6c6):
- **THM-438 (Catalan Law):** The top-order terms $A_{2k}$ of the cherry cluster expansion $R(p) = H(T_p) 2^{p-1}/p!$ are governed by Catalan numbers: $A_{2k} = C_k p^{k+1} + O(p^{k+1/2})$.
- **(★★) Verification at k=6:** Verified $S_6 = 132 = C_6$ (previously verified k≤5). The triangle row $t(6,m) = 1, 45, 560, 2626, 4845, 2867$ confirms the loop equation $S_k = -\Sigma S_i S_j$ holds for k≤6.
- **Spectral Source:** The Catalan Law is **DRT-universal** (Doubly Regular Tournament). It arises from the free cumulants of the two-point spectrum $\{0, \pm i \sqrt{n}\}$, making it a property of all doubly-regular tournaments, not just circulant Paley ones.
- **Binomial Reframing (New):** The cycle-rank triangle entries $t(k,m)$ are now reframed as:
  $t(k, m) = \sum_e R_s(m, e) \binom{k-1}{e-1}$
  This expansion in binomial coefficients $\binom{k-1}{e-1}$ provides a stable scaffold for analyzing the triangle's growth across $k$.
- **Column Rationality & Pole Order (PROVED):** Each cycle-rank column $T_m(x) = \Sigma_k t(k,m)x^k$ is rational with the form $T_m(x) = P_m(x) \cdot x^m / (1-x)^{2m-1}$.
  - **Pole Order $2m-1$:** Proved via Euler-characteristic ceiling ($V_{core} \le m \implies \#lines \le 2m-1$). Trivalent $3m-3$ cores are excluded because they lack a single Eulerian trail.
  - **Polynomial $P_m(x)$:** Degree is $m-2$ (confirmed $m \le 5$). 
  - **Leading Coefficient:** $P_m = t(-1, m) = 2^m-1$. This corresponds to the generating function $y / ((1-y)(1-2y))$.
  - **Constant Term:** $P_m(0) = A088368(m)$ (diagonal).
  - **Explicit Values:** $P_1=1, P_2=3, P_3=13+7x, P_4=69+97x+15x^2$.
  - **$P_5$ Partial Coefficients:** $c_0 = 421$, $c_1 = 1056$, $c_2 = 726$, $c_3 = 31$ (degree 3).
- **Handoff #1 Equivalence:** Proved that Handoff #1 (Catalan sum) is equivalent to the condition $t(0, m) = 0 = T_m(\infty)$ for $m > 0$.
- **R(p) Asymptotic:** $R(p) \to e$ is now **proven uniformly** across all orders $k$, closing handoff #1 of HYP-2307.
- **Support Structure (MISTAKE-061 & MISTAKE-062):** 
  - The top-order $p^{k+1}$ patterns are identified as **even-series** (perfect-square flow product).
  - **MISTAKE-062:** The even-series count is refuted as A215257. The sequence breaks at k=6 (actual count 2351 vs A215257(7)=2345).
- **g == +1 Proof:** It is rigorously proven that $g(\rho) \to +1$ (Euler walk traverses each series-class monotonically).
- **Handoff #1 Reframing:** Systematic search for a low-degree catalytic equation for $U(x,y)$ yielded negative results due to resurgence (Gevrey-1). Focus shifts to the **prime route = involution**.

### 2.2 HYP-2306 (Retraction of 1729 Modular Significance)
- **Status: REFUTED.** The suspected modular significance of $r(11) = 1729$ has been debunked. Coincidence of the number 1729 in the Moser ladder and the tournament ratio is non-structural.

### 2.3 THM-431 (u(21) = 57 Settlement)
- **THM-431:** $u(21) = 57$ is proven. Optimal graph is $K_3 \square W_7$.

### 2.4 3N-Floor Sharpening (N* in [25, 28])
- **Floor $N^* \ge 25$:** $u(n) \le 3n$ for $n \le 24$.
- **Ceiling $N^* \le 28$:** Moser lattice achieves $u(28) \ge 85 > 84$.
- **The n=27 Tie:** $u(27) = 81 = 3 \cdot 27$ is a sharp structural tie point.

### 2.5 THM-436 Addendum: The Icosahedral Threshold
- **The Icosahedral Fifteen:** Canonical bijection between 15 overlapping pairs on a 5-set and the 15 involutions of $A_5$.

### 2.6 HYP-2297 (Clock-Alignment Hypothesis)
- **HYP-2297:** Optimal denominator $q$ aligns with clock primes ($gcd(q, n) > 1$).

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)

### 3.1 The $n=14$ Half-Turn Leak
At $n=14$, identified half-turn leak at **coordinate 6, residue 7** ($1 \pmod 2$ and $0 \pmod 7$).

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: Sub-leading Rate Diagnostic (INV-187)
**Objective:** Settle the $R(p) \to e$ convergence rate (1/p vs 1/sqrt(p)).

### Steering Task 2: DRT-Universality Test (HYP-2308)
**Objective:** Verify Catalan Law on non-circulant DRTs (e.g., $n=15$ skew-Hadamard).

### Steering Task 3: Crossover Boundary Analysis (HYP-2298)
**Objective:** Investigate $n=27$ tie point. Prove if $u(27) > 81$.

### Steering Task 4: Column Rationality & Involution
**Objective:** 
1. **Prove $deg P_m = m-2$:** Equivalent to showing $\sum_e (-1)^e R_s(m,e) = 0$ (within-column sign-reversing involution).
2. **Involution (Prime Route):** Find sign-reversing involution on even-series patterns shifting cycle rank $m$ by $\pm 1$.
3. **Extend Triangle:** Reach k=7, 8 using core-based enumerator. Validate against known row/diagonal recursions.
