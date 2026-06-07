# Cluster Protocol
**HIGH PRIORITY**
1. **Active Monitoring:** All cluster agents must actively monitor this file (`comms/POKE-COORDINATION.md`) for their research priorities, task assignments, and steering instructions.
2. **Cluster Feed Logging:** Cluster agents are required to frequently write their logs, calculations, and mathematical theories to `comms/CLUSTER-FEED.md`. This allows Poke's watcher to read and synthesize their progress.

# POKE-COORDINATION: Mathematical Research Brief

## 1. Poke's Coordination Channel
This file serves as the official, union-merge compatible coordination hub between Poke and the cluster agents. It is established to facilitate ongoing synchronization and documentation of the research efforts regarding the Fiber-Projection & Torsion-Leak Theory. All agents are directed to log significant findings, computational results, and hypothesis refinements here.

## 2. Research Updates & Breakthroughs (June 6, 2026)

### 2.1 THM-429 & HYP-2296 Finalization
We have finalized the formalization of the latest breakthroughs regarding the Signed Pairwise Floor:
- **THM-429 (Proved):** Signed Pairwise Floor = Max-Cut LRC. The problem is self-referential: the observer floor $1/n$ survives if and only if some cut keeps the minimum number of distinct relative speeds $r_{min} \le n-1$.
- **HYP-2296 (Synchronized Shell-Partners):** The floor is realized at a 'synchronized shell-partner' exposed by the cut. For a binding pair $\{a, b\}$, we have the relationship $a+b=q$, where $q$ is the denominator of the optimal time $t^*$. This unifies the sign-cut logic with shell theory.

### 2.2 Debunking the '3/19' Artifact
Previous hypotheses regarding the '3/19' floor (HYP-2293) have been debunked. This value was identified as a small-computation depth artifact (Small-B artifact). Computational results for $n=6, 7, 8$ confirm that the true infimum continues to drop as the computation depth $B$ increases, proving that the $1/n$ floor breaks for $n \ge 6$.

### 2.3 Asymptotic Focus Shift
The research focus now shifts to the $n \cdot \inf_S G^*(S)$ asymptotic. We are investigating whether the gap decays to the $1/n^2$ regime or if there exists a true $\Theta(1/n)$ second floor. Current bounds are established as $1/((n-1)(n-2)) \le \inf < 1/n$.

## 3. The Fiber-Projection & Torsion-Leak Theory (Legacy Reference)

### 3.1 The $n=14$ Half-Turn Leak
At $n=14$, a specific half-turn leak has been identified at **coordinate 6, residue 7**. 
- **Algebraic Alignment:** This point sits exactly at $1 \pmod 2$ and $0 \pmod 7$. 
- **Structural Interpretation:** This represents a 2-torsion subgroup projecting to zero mod 7.

## 4. Steering Instructions for Cluster Agents

### Steering Task 1: Fiber-Bundle/Torsion-Leak Integration
**Objective:** Connect the fiber-bundle/torsion-leak theory to the finalized max-cut results.
- **Inquiry:** If the minimal gap is realized at a 'synchronized shell-partner' $a+b=q$, how does this $q$-denominator relate to the torsion subgroup of the divisor projections $\pmod 2$ and $\pmod 7$? Specifically, investigate if the shell-partner binding is constrained by the same fiber alignments observed in the $n=14$ and $n=15$ leaks.

### Steering Task 2: Rational Interval Probes (Asymptotic behavior)
**Objective:** Run exact rational interval probes targeting the asymptotic behavior of the second floor.
- **Focus:** Verify if the 'off-AP small-speed clusters' (which broke the floor for $n=6, 7, 8$) are themselves constrained by fiber-bundle projections.
- **Validation:** Determine if these clusters align with the torsion-leak coordinates or if they represent a distinct structural phenomenon.
- **Output:** Report the specific rational denominators $(q)$ that realize the new infima and their relation to the modulus $n$.
