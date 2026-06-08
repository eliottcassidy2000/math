# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 8, 2026)

### 2.1 THM-441 (The Convolution-Correlation-Adjoint Duality)
A fundamental unification of the repository's core mathematical objects has been achieved (S706, THM-441, HYP-2311):
- **Core Identity:** Cross-correlation is the **adjoint** of convolution ($\langle h*a, b \rangle = \langle a, h \star b \rangle$).
- **Clock-Shell Unification:** The **Shell face** (sums $v_i+v_j$ mod $2n-1$, S700) is a convolution $1_S * 1_S$, while the **Clock face** (differences mod $n$, S701) is a cross-correlation $1_S \star 1_S$. They are adjoint operators related by the antipodal involution $\sigma: x \mapsto -x$.
- **Converse Duality:** The tournament converse $T \mapsto T^*$ is the adjoint of the circulant adjacency-convolution. **Worry-sets** (THM-402) are precisely the **self-adjoint** (symmetric) operators.
- **Additive Energy:** Unified as the norm squared of the autocorrelation ($\|1_S \star 1_S\|^2$), which by Wiener–Khinchin is $\sum |\hat{1_S}|^4$.
- **Spectral Rigidity:** Self-adjointness is identified as the root cause of the "worry-set" rigidity and cyclotomic LRC properties.

### 2.2 THM-438 (The Catalan Law of Paley Cluster Integrals)
A definitive breakthrough has been achieved in the analytic theory of Paley tournaments (Commit 79d0590, verified k=6 in 91d872b, free-probability unification in 2b24b5e, Bercovici-Pata partnership in 62f8976, Watson bridge in d65eda8):
- **THM-438 (Catalan Law):** The top-order terms $A_{2k}$ of the cherry cluster expansion $R(p) = H(T_p) 2^{p-1}/p!$ are governed by Catalan numbers: $A_{2k} = C_k p^{k+1} + O(p^{k+1/2})$.
- **The Watson Bridge & Resurgent Asymptotics (PROVED):** The wild diagonal $A088368$ and the density $\rho(x)$ are unified as two faces of a single resurgent series.
  - **Shared Coefficients:** The tail of the density $\rho(x) \sim e^{1-x} \sum b_j x^{-j}$ and the moments $A088368(k)/k! \sim e \sum b_j / (k)_j$ share the **exact same asymptotic coefficients** $b = 1, 2, 10, 178/3, \dots$.
  - **Verification:** Identity verified via 4 independent routes (tail match to $1e-16$, moments up to $k=60$).
  - **Two Humps = One Series:** The spectral density's "two-hump" structure (ADD-15) is resolved into a single resurgent object; the wildness of the moments and the decay of the tail are governed by the same non-perturbative physics.
  - **OEIS Refinement:** The previously observed $a(n) \sim e \cdot n!$ is now rigorously refined to $a(n) \sim e \cdot n! (1 + 2/n + 10/(n(n-1)) + \dots)$.
  - **New Constant:** $A088368(8) = 175769$.
- **The Bercovici-Pata Partnership (PROVED):** The wild (x=0) and tame (x=∞) endpoints are **Bercovici-Pata partners** ($\Lambda(\mu_{classical}) = \mu_{free}$), sharing the cumulant signature $\kappa_n = n!$.
  - **The Borel-Bridge Identity:** $z R(z) = \int_0^\infty e^{-t} C(zt) dt$. This identity unifies the analytic resurgence (ADD-6) with the free-probabilistic structure (ADD-11).
  - **Positive-Definite Measure Family $\mu_q$:** The classical (A000262) and free (A088368) endpoints are slices ($q=1$ and $q=0$) of a crossing-graded measure family $\mu_q$ (Bożejko–Speicher $q$-deformation).
- **The Diagonal is a Free Moment Sequence (PROVED):** $t(k,k) = A088368(k)$ is the free moment sequence of the "factorial law" where free cumulants $\kappa_n = n!$.
- **Free-Probabilistic Endpoint Framing:**
  - **Wild Endpoint (x=0):** Free MOMENTS of the factorial law ($t(k,k) = A088368$).
  - **Tame Endpoint (x=∞):** Free CUMULANTS of the two-point law ($\Sigma (-1)^m t(k,m) = (-1)^k C_k$).
- **Over-count Correction (A000262):** The unconstrained over-count corresponds to the **CLASSICAL moments** of the factorial law, A000262 ("sets of lists").
- **Spectral Source:** The Catalan Law is **DRT-universal**, arising from the free cumulants of the two-point spectrum $\{0, \pm i \sqrt{n}\}$.
- **Handoffs as Finite Alternating-Binomial Sums:** 
  - **Handoff #1 (Catalan Sum):** $\sum_{k=m}^{2m-1} (-1)^k \binom{2m-1}{k} t(k,m) = 0$.
  - **Handoff #2 (Lead Coefficient):** $\sum_{k=m}^{2m-1} (-1)^{k+1} k \binom{2m}{k+1} t(k,m) = 2^m-1$.
- **Clean Column Form (Falling Factorials):** Each cycle-rank column follows $t(k,m) = (k)_m \cdot h_m(k) = m! \binom{k}{m} h_m(k)$.
- **Column Rationality & Pole Order (PROVED):** $T_m(x) = P_m(x) \cdot x^m / (1-x)^{2m-1}$. Pole order $2m-1$ is an Euler-characteristic ceiling ($V_{core} \le m$).
- **Polynomial $P_m(x)$:** Degree $m-2$ (confirmed $m \le 5$). Leading coefficient $2^m-1$. Constant term $A088368(m)$.
- **Free Density & Cauchy Inversion:** Density $\rho$ refines the edge constant to $\approx 0.4-0.6$.

### 2.3 THM-440 (u(22) unit-cocyclic extension reduction)
- **THM-440 (PROVED):** A 61-edge unit-distance graph (UDG) on 22 vertices must delete a degree-4 or degree-5 vertex to a dense 21-vertex core. The deleted vertex's neighbours are a unit-cocyclic delta-set.
- **Moser Ring Caps at 60:** Verification on the 12 densest-21 cores in $M_L$ shows the delta=4 route is empty. Within $M_L$, $u(22)=60$.

### 2.4 HYP-2310 (u(22)=61 and Field Escape)
- **HYP-2310:** $u(22)=61$ is only attainable via the delta=5 route (extending 56-edge 21-cores) or by escaping to $\mathbb{Q}(\sqrt{-3}, \sqrt{-d})$.

### 2.5 HYP-2306 (Retraction of 1729 Modular Significance)
- **Status: REFUTED.** Coincidence of 1729 in Moser ladder and tournament ratio is non-structural.

### 2.6 THM-431 (u(21) = 57 Settlement)
- **THM-431:** $u(21) = 57$ is proven. Optimal graph is $K_3 \square W_7$.

### 2.7 3N-Floor Sharpening (N* in [25, 28])
- **Floor $N^* \ge 25$:** $u(n) \le 3n$ for $n \le 24$.
- **The n=27 Tie:** $u(27) = 81 = 3 \cdot 27$ is a sharp structural tie point.

### 2.8 THM-436 Addendum: The Icosahedral Threshold
- **The Icosahedral Fifteen:** Canonical bijection between 15 overlapping pairs on a 5-set and the 15 involutions of $A_5$.

### 2.9 HYP-2297 (Clock-Alignment Hypothesis)
- **HYP-2297:** Optimal denominator $q$ aligns with clock primes ($gcd(q, n) > 1$).

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)

### 3.1 The $n=14$ Half-Turn Leak
At $n=14$, identified half-turn leak at **coordinate 6, residue 7** ($1 \pmod 2$ and $0 \pmod 7$).

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: THE ADJOINT SYMMETRY PROBE (NEW - S706)
**Objective:** Map every repo object as a convolution operator and identify its adjoint (converse / σ-reflected face).
1. **Self-Adjoint Census:** Confirm if all "worry-set" configurations are exactly the self-adjoint operators in their respective cyclotomic groups.
2. **Autocorrelation Moments:** Re-examine THM-406 covering-depth moments as autocorrelation integrals. Determine if the "Vitali wall" is the autocorrelation's failure to be a finite positive character-combination.

### Steering Task 2: Sub-leading Rate Diagnostic (INV-187)
**Objective:** Settle the $R(p) \to e$ convergence rate (1/p vs 1/sqrt(p)).

### Steering Task 3: DRT-Universality Test (HYP-2308)
**Objective:** Verify Catalan Law on non-circulant DRTs (e.g., $n=15$ skew-Hadamard).

### Steering Task 4: Crossover Boundary Analysis (HYP-2298)
**Objective:** Investigate $n=27$ tie point. Prove if $u(27) > 81$.

### Steering Task 5: THE u(22) SETTLEMENT (NEW - S705)
**Objective:** Confirm if $u(22)=60$ or 61.
1. **delta=5 Census:** Systematically generate all 56-edge 21-cores and test for unit-cocyclic degree-5 extensions.
2. **Field Escape Check:** Probe for 22-point configurations in $\mathbb{Q}(\sqrt{-3}, \sqrt{-d})$ for new Heegner d.

### Steering Task 6: THE TWO-LAW TRANSIT (High Priority)
**Objective:** Derive the explicit Möbius map between the Factorial Law (diagonal, free cumulants $\kappa_n = n!$) and the Two-Point Law (signed row sum) to resolve the off-diagonal $P_m$ polynomials.

### Steering Task 7: CYCLOTOMIC DEPTH CENSUS
**Objective:** Run a systematic $q$-census for $n=9$ to $n=14$ to verify if witness denominators $q$ align with clock primes (HYP-2297: $gcd(q, n) > 1$), specifically testing if $R(14)$ has depth equal to the $3^3$ prime-power tower.

### Steering Task 8: TAME HANDOFF PROOFS
**Objective:** Close the remaining gaps by showing $(k)_m | t(k, m)$ and deriving the exact $g_m(-1)$ values.

### Steering Task 9: HAAGERUP-MÖLLER ANALYSIS
**Objective:** Rigorously prove the $1/\sqrt{x}$ edge and the $e \cdot e^{-x}$ tail of the free compound Poisson distribution.

### Steering Task 10: Column Rationality & Involution (Legacy/Ongoing)
**Objective:** 
1. **Involution (Prime Route):** Find sign-reversing involution on even-series patterns shifting cycle rank $m$ by $\pm 1$.
2. **Diagonal Census:** Prove the bijection *doubled plane trees ↔ Callan noncrossing lists*.
3. **Watson Bridge Analysis:** Use the shared $b_j$ coefficients to derive the off-diagonal $h_m(k)$ polynomials.
4. **Extend Triangle:** Reach k=7, 8 using core-based enumerator. Validate against row/diagonal recursions.
