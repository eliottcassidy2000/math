# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 6, 2026)

### 2.1 THM-430 (Antipodal Involution Unification)
A major math breakthrough has been achieved in commit a3e6b37fcfc9917b9d0dab1fd46f5fd1b434ca41:
- **THM-430 (Antipodal Involution Unification):** This theorem elegantly unifies the shell-partner q-denominator and the torsion leak. 
- **Significance:** The half-turn (2-torsion) is identified as the fixed points of the antipodal involution. While this creates the loudest "cell leak," it is structurally excluded from being the binding minimizer for the floor.
- **Odd-Prime Fiber Resolution:** Because -1 = 1 mod 2, the 2-fiber is inert/trivial. Real antipodal structure lives entirely in the odd-prime divisor projections (mod 7, mod 3, etc.).

### 2.2 HYP-2297 (Clock-Alignment Hypothesis)
- **HYP-2297:** The optimal denominator q aligns with clock primes (gcd(q, n) > 1) rather than shell primes.

### 2.3 THM-429 & HYP-2296 Finalization
We have finalized the formalization of the latest breakthroughs regarding the Signed Pairwise Floor:
- **THM-429 (Proved):** Signed Pairwise Floor = Max-Cut LRC. The problem is self-referential: the observer floor $1/n$ survives if and only if some cut keeps the minimum number of distinct relative speeds $r_{min} \le n-1$.
- **HYP-2296 (Synchronized Shell-Partners):** The floor is realized at a 'synchronized shell-partner' exposed by the cut. For a binding pair $\{a, b\}$, we have the relationship $a+b=q$, where $q$ is the denominator of the optimal time $t^*$. This unifies the sign-cut logic with shell theory.

### 2.4 Debunking the '3/19' Artifact
Previous hypotheses regarding the '3/19' floor (HYP-2293) have been debunked. This value was identified as a small-computation depth artifact (Small-B artifact). Computational results for $n=6, 7, 8$ confirm that the true infimum continues to drop as the computation depth $B$ increases, proving that the $1/n$ floor breaks for $n \ge 6$.

### 2.5 Asymptotic Focus Shift
The research focus now shifts to the $n \cdot \inf_S G^*(S)$ asymptotic. We are investigating whether the gap decays to the $1/n^2$ regime or if there exists a true $\Theta(1/n)$ second floor. Current bounds are established as $1/((n-1)(n-2)) \le \inf < 1/n$.

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)

### 3.1 The $n=14$ Half-Turn Leak
At $n=14$, a specific half-turn leak has been identified at **coordinate 6, residue 7**. 
- **Algebraic Alignment:** This point sits exactly at $1 \pmod 2$ and $0 \pmod 7$. 
- **Structural Interpretation:** This represents a 2-torsion subgroup projecting to zero mod 7.

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: Next Step on HYP-2297
**Objective:** Execute a systematic computational census of the optimal denominator q for n=9 through n=14.
- **Goal:** Prove the gcd(q, n) > 1 clock-alignment hypothesis.

### Steering Task 2: Off-AP Small-Speed Cluster Analysis
**Objective:** Analyze the "off-AP small-speed clusters" in the 3-fiber (mod 3) and 7-fiber (mod 7) projections.
- **Inquiry:** Determine if their structural obstruction is governed by the non-trivial orbits of the antipodal involution.
