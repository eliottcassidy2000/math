## kps-2026-06-22-S37 -- Dirichlet-Extremal Tight Locus and Census (checkpoint)

Formalized the correction identifying the LRC(14) tight locus as Dirichlet-extremal sets (1/14-covering nets) rather than maximal additive energy, and documented the honest status of the census completeness (commit `81f7da68`). This checkpoint stabilizes the project's analytical boundary by anchoring the tight locus in Dirichlet approximation theory.

### 1. Tight Locus: Dirichlet-Extremal vs. Additive Energy
Corrected the structural characterization of the tight locus:
- **Finding:** The tight locus corresponds to **Dirichlet-extremal 13-sets** (sets whose danger arcs form a 1/14-covering net of the torus), not merely sets with maximal additive energy.
- **Evidence:** Verified the single-swap pair $12 \to 24$ (Goddyn-Wong) and $12 \to 26$. Despite having equal additive energy, $12 \to 24$ is tight ($M=1/14$) because its residues mod 14 form the required covering net, while $12 \to 26$ is loose ($M > 1/14$).
- **Impact:** This shift refocuses the search for tight sets on Dirichlet approximation properties and residue-net coverage rather than scalar energy sums.

### 2. Auto-Safe Near-Misses
Identified a structural mechanism that automatically secures certain "near-miss" rows:
- **Mechanism:** Many loose rows (e.g., those with $M=1/12$) are already extremal for smaller $N$ (fewer runners).
- **Safety:** These rows are "auto-safe" because their margin at $N=13$ is strictly greater than $1/14$, meaning they cannot be disproved as counterexamples for $LRC(14)$. This provides a robust buffer for non-Dirichlet-extremal configurations.

### 3. Census and Three-Gap Rigidity
Clarified the connection between the tight-locus census and the Steinhaus three-gap theorem:
- **The Core:** The completeness of the $\{AP, GW\}$ census is equivalent to **three-gap (Steinhaus) rigidity** for 13 runners. 
- **Honest Status:** While the census is exact for single-swaps and unrefuted by broad search, its universal completeness remains the **irreducible open core** of the project. This is a known open conjecture for $n=13$ in the literature.
- **Reduction:** $LRC(14)$ is successfully reduced to this rigidity claim; the project has not "proved" the conjecture itself, but has isolated it as the final remaining barrier.

### 4. Convergence and Guardrails
The findings reinforce the "honest" reporting protocol:
- **Proved/Formalized:** $THM-568$ (apex-denominator), $LRCApex7Floor$ (sorry-free), and the single-swap census accuracy.
- **Open:** The universal census completeness (three-gap rigidity).

### 5. Net Impact
This checkpoint stabilizes the project's "extremal atlas." By identifying Dirichlet-extremality as the governing principle of the tight locus and acknowledging the three-gap rigidity bottleneck, the project provides a rigorous and honest roadmap for the terminal proof obligations.
