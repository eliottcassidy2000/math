## mac-mini-2026-06-22-S51 -- Census and Status of Tightness-Star (checkpoint)

Formalized the "forward" proof of the tightness-star theorem, verified the single-swap census of tight rows, and clarified the project's honest status regarding the completeness of the LRC(14) closure (commit `2a852ec4`). This checkpoint defines the "irreducible open core" of the project.

### 1. Tightness-Star Forward: Proved and Formalized
Rigorously proved the forward direction of the (star) crux: $M(S)=1/14 \implies$ the optimum $t^*$ is an apex-7 antipodal point with a binding pair $14 \mid (s_i + s_j)$.
- **Analysis:** Since $1/14 < 1/2$, a local maximum of the minimum-of-sawtooths function requires a "crossing" of an increasing active runner ($\text{frac} = 1/14$) and a decreasing one ($\text{frac} = 13/14$).
- **Arithmetic:** Proved that such a crossing forces the denominator $D$ to satisfy $14 \mid D$, and specifically that $14 \mid (s_i + s_j)$.
- **Lean Formalization:** The arithmetic core is machine-verified and sorry-free in `LRCBindingPair.lean` (`LonelyRunner.BindingPair.binding_pair_dvd`).

### 2. Single-Swap Census: AP and GW
Conducted an exact census of single-speed swaps from the AP $\{1, \dots, 13\}$ for replacements $r \le 300$.
- **Result:** The only tight ( $M=1/14$ ) non-AP set is the Goddyn-Wong-type set $\{1, \dots, 11, 13, 24\}$ (the $12 \to 24$ swap).
- **Consistency:** Every tight row in the census has an optimum at denominator exactly 14, consistent with the binding pairs $\{1, 13\}, \{5, 9\},$ and $\{3, 11\}$.

### 3. Apex Floor Formalized
Formalized the "Apex Floor" lemma in Lean (`LRCApex7Floor.D14_never_certifies`).
- **Theorem:** At any denominator-14 witness point $a/14$, a covering runner (a multiple of 14) necessarily sits on the observer.
- **Implication:** Combined with the forward theorem, this shows that a covering set cannot be tight at a denominator-14 optimum.

### 4. The Irreducible Open Core (Honest Status)
Clarified that while the project has reduced $LRC(14)$ to a single named statement (the census completeness), the proof of that completeness remains an open core.
- **Census Completeness:** The claim that no tight sets exist beyond $\{AP, GW\}$ (e.g., no multi-swap or unbounded exceptions) is equivalent to the Steinhaus/consecutive-maximization conjecture for 13 runners.
- **Status:** This remains an open problem in the literature. The project's bounded search provides strong evidence but does not constitute an exhaustive proof.
- **Disclaimer:** The cluster reports verified reductions and formalizations, not a final proof of $LRC(14)$.

### 5. Net Impact
This checkpoint stabilizes the project's analytical boundary. By proving the forward forcing and the apex-floor contradiction, the $LRC(14)$ problem is successfully focused on the "census completeness" bottleneck. The project now explicitly identifies where the machine-verified and computationally-verified progress ends and where the fundamental open conjecture (consecutive-maximizes) begins.
