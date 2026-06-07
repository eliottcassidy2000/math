# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 7, 2026)

### 2.1 THM-438 (The Catalan Law of Paley Cluster Integrals)
A definitive breakthrough has been achieved in the analytic theory of Paley tournaments (Commit 79d0590, verified k=6 in 91d872b, free-probability unification in 2b24b5e):
- **THM-438 (Catalan Law):** The top-order terms $A_{2k}$ of the cherry cluster expansion $R(p) = H(T_p) 2^{p-1}/p!$ are governed by Catalan numbers: $A_{2k} = C_k p^{k+1} + O(p^{k+1/2})$.
- **The Diagonal is a Free Moment Sequence (PROVED):** The diagonal $t(k,k) = A088368(k)$ is exactly the free moment sequence of the "factorial law" where free cumulants $\kappa_n = n!$.
  - **Proof via Callan's Bijection:** The functional equation $M(z) = F(zM(z))$ for $F = \Sigma n! w^n$ establishes $\sum_{NC(k)} \prod_{B \in \pi} |B|! = A088368(k)$.
  - **Resurgence Unification:** The Gevrey-1 resurgence of $U(x,y)$ (ADD-6) and the free-probability structure (ADD-4) are unified: the resurgence IS the divergence of the $R$-transform $R(z) = \Sigma n! z^{n-1}$.
- **Free-Probabilistic Endpoint Framing:** Both endpoints of the bridge are now framed as free-probabilistic structures of different laws:
  - **Wild Endpoint (x=0):** Free MOMENTS of the factorial law ($t(k,k) = A088368$).
  - **Tame Endpoint (x=∞):** Free CUMULANTS of the two-point law ($\Sigma (-1)^m t(k,m) = (-1)^k C_k$).
  - **Structural Anonymous Path:** The Möbius transit between these laws does not correspond to a standard moment-cumulant map, explaining why off-diagonal sequences are OEIS-negative.
- **Over-count Correction (A000262):** The unconstrained over-count (all partitions) corresponds to the **CLASSICAL moments** of the same factorial law ($\kappa_n = n!$), which is A000262 ("sets of lists"). The gap $A000262 - A088368$ represents the crossing-partition excess.
- **(★★) Verification at k=6:** Verified $S_6 = 132 = C_6$ (previously verified k≤5). The triangle row $t(6,m) = 1, 45, 560, 2626, 4845, 2867$ confirms the loop equation $S_k = -\Sigma S_i S_j$ holds for k≤6.
- **Spectral Source:** The Catalan Law is **DRT-universal**. It arises from the free cumulants of the two-point spectrum $\{0, \pm i \sqrt{n}\}$.
- **Handoffs as Finite Alternating-Binomial Sums:** 
  - **Handoff #1 (Catalan Sum):** $\sum_{k=m}^{2m-1} (-1)^k \binom{2m-1}{k} t(k,m) = 0$.
  - **Handoff #2 (Lead Coefficient):** $\sum_{k=m}^{2m-1} (-1)^{k+1} k \binom{2m}{k+1} t(k,m) = 2^m-1$.
- **Clean Column Form (Falling Factorials):** Each cycle-rank column follows $t(k,m) = (k)_m \cdot h_m(k) = m! \binom{k}{m} h_m(k)$.
- **Column Rationality & Pole Order (PROVED):** $T_m(x) = P_m(x) \cdot x^m / (1-x)^{2m-1}$. Pole order $2m-1$ is an Euler-characteristic ceiling ($V_{core} \le m$).
- **Polynomial $P_m(x)$:** Degree $m-2$ (confirmed $m \le 5$). Leading coefficient $2^m-1$. Constant term $A088368(m)$.
- **Asymptotic Correction (MISTAKE-063):** Retracted the previous refutation. $A088368(m) \sim e \cdot m!$ is correct; the ratio overshoots $e$ (peaking at $m=8$) before descending back.
- **Support Structure (MISTAKE-062):** Even-series count refuted as A215257 (breaks at k=6: 2351 vs 2345).

### 2.2 HYP-2306 (Retraction of 1729 Modular Significance)
- **Status: REFUTED.** Coincidence of 1729 in Moser ladder and tournament ratio is non-structural.

### 2.3 THM-431 (u(21) = 57 Settlement)
- **THM-431:** $u(21) = 57$ is proven. Optimal graph is $K_3 \square W_7$.

### 2.4 3N-Floor Sharpening (N* in [25, 28])
- **Floor $N^* \ge 25$:** $u(n) \le 3n$ for $n \le 24$.
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
1. **Involution (Prime Route):** Find sign-reversing involution on even-series patterns shifting cycle rank $m$ by $\pm 1$.
2. **Diagonal Census:** Prove the bijection *doubled plane trees ↔ Callan noncrossing lists*.
3. **Distribution Probe:** Investigate if the factorial law ($\kappa_n = n!$) is a named distribution.
4. **Tame Handoffs:** Prove Handoff #1 ($t(0,m)=0$) and Handoff #2 ($g_m(-1) = (-1)^m(2^m-1)$).
5. **Extend Triangle:** Reach k=7, 8 using core-based enumerator. Validate against row/diagonal recursions.
