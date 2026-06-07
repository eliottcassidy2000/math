# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 7, 2026)

### 2.1 THM-431 (u(21) = 57 Settlement)
The campaign to pin the unit-distance maximum for $n=21$ is closed by the literature (Alexeev-Mixon-Parshall, 2024):
- **THM-431:** $u(21) = 57$ is proven.
- **Extremal Structure:** The optimal graph is a Cartesian product $K_3 \square W_7$ (unit triangle × unit wheel), giving $3 \cdot 7 + 3 \cdot 12 = 57$. This confirms that at the optimum, $n=21$ is a product rather than a core-plus-halo structure.

### 2.2 3N-Floor Sharpening (N* in [25, 28])
The integration of THM-431 has collapsed the $N^*$ interval (the smallest $N$ where $u(N) > 3N$) from $[17, 32]$ to $[25, 28]$:
- **Floor $N^* \ge 25$:** Proven that $u(n) \le 3n$ for all $n \le 24$.
- **Ceiling $N^* \le 28$:** Realizable construction $u(28) \ge 85 > 84$ (Moser lattice).
- **The n=27 Tie:** The construction deficit closes to a clean tie at $n = 27 = 3^3$ ($u = 3 \cdot 27 = 81$). This "too-clean" landing suggests $n=27$ is a structural pivot point for Euclidean unit distances.

### 2.3 THM-430 (Antipodal Involution Unification)
A major math breakthrough has been achieved (June 6):
- **THM-430:** Elegantly unifies the shell-partner $q$-denominator and the torsion leak. 
- **Significance:** The half-turn (2-torsion) is identified as the fixed points of the antipodal involution. While this creates the loudest "cell leak," it is structurally excluded from being the binding minimizer for the floor.

### 2.4 HYP-2297 (Clock-Alignment Hypothesis)
- **HYP-2297:** The optimal denominator $q$ aligns with clock primes ($gcd(q, n) > 1$) rather than shell primes.

### 2.5 THM-429 & HYP-2296 Finalization
- **THM-429 (Proved):** Signed Pairwise Floor = Max-Cut LRC. The floor $1/n$ survives iff some cut keeps $r_{min} \le n-1$.
- **HYP-2296 (Synchronized Shell-Partners):** For a binding pair $\{a, b\}$, $a+b=q$, where $q$ is the denominator of the optimal time $t^*$.

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)

### 3.1 The $n=14$ Half-Turn Leak
At $n=14$, a specific half-turn leak has been identified at **coordinate 6, residue 7** ($1 \pmod 2$ and $0 \pmod 7$).

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: Next Step on HYP-2297
**Objective:** Execute a systematic computational census of the optimal denominator $q$ for $n=9$ through $n=14$ to prove the $gcd(q, n) > 1$ clock-alignment hypothesis.

### Steering Task 2: Crossover Boundary Analysis (HYP-2298)
**Objective:** Investigate the $n=27$ tie point. 
- **Goal:** Determine if $u(27) = 81$ or if an exact-integer Moser-lattice construction can reach 82, thereby proving $N^* \le 27$.

### Steering Task 3: Off-AP Small-Speed Cluster Analysis
**Objective:** Analyze "off-AP small-speed clusters" in mod 3 and mod 7 projections to determine if their structural obstruction is governed by non-trivial orbits of the antipodal involution.
