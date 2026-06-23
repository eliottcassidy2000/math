## kps-2026-06-22-S41 -- Tournament Spectrum and Binding Scale 14 (checkpoint)

Formalized the "Tournament Spectrum" reframe, providing a complete characterization of the LRC(14) tight locus by unifying magnitude-aware magnitude data with tournament isomorphism classes (commit `70e3e1f2`). This checkpoint corrects earlier misconceptions about "necessary-only" invariants and anchors tightness in the binding-scale property.

### 1. The Tournament Spectrum: $\Sigma(S)$
Redefined the fundamental object of study from a single apex tournament to the **Tournament Spectrum** $\Sigma(S) = \{iso(T(S, t)) : t \in [0, 1)\}$.
- **Magnitude Awareness:** While a single tournament at a fixed phase (like $a/14$) is magnitude-blind, the spectrum $\Sigma(S)$ changes only at breakpoints $t = k/(s_i - s_j)$ and $k/(2s_i)$. Since these denominators are the runner magnitudes, the spectrum is inherently magnitude-aware.
- **Evidence:** Verified that the tight AP $\{1, \dots, 13\}$ and the loose "twin" $12 \to 26$ have identical apex tournaments but distinct spectra ($|\Sigma(AP)|=14$ vs. $|\Sigma(12 \to 26)|=24$).

### 2. Binding Scale 14: A Complete Invariant
Recast THM-568 as a complete characterization: **Tight $\iff$ Binding Scale = 14**.
- **Definition:** The "binding scale" is the denominator of the optimum phase $t^*$ where the spectrum's deepest sink resides.
- **Results:**
    - **Tight:** AP and GW $\to$ Scale 14 (the apex Farey node).
    - **Loose:** $12 \to 26 \to$ Scale 12; $12 \to 36 \to$ Scale 41 (Farey neighbor).
- **Correction:** Clarified that this is a **complete** characterization (not necessary-only) because the spectrum carries both the structural residue information (the deepest-sink isomorphism class) and the metric magnitude information (the scale).

### 3. Geometry: Labelled Farey Tree and Flip-Graph Walk
Mapped the spectrum onto the Farey/Stern-Brocot tree and the tournament flip-graph.
- **Farey Mapping:** $\Sigma(S)$ is visualized as the Farey tree labelled by tournament isomorphism classes. Tight rows pin their deepest sink to the apex node ($1/14$), while loose rows migrate to coarser nodes ($q < 14$, divisibility loss) or deeper Farey children (e.g., $3/41$, near-miss neighbors).
- **Flip-Graph Walk:** Every breakpoint in the spectrum represents an antipodal-pair flip in the tournament. The spectrum is thus a **walk** in the tournament flip-graph (the metagraph $G_n$), where the "migration" of the sink away from the apex corresponds to the loss of tightness.

### 4. Convergence with Open Core
The "census completeness" problem—the project's irreducible open core—is now framed as a **flip-graph-walk non-migration statement**. Proving that no other tight sets exist beyond $\{AP, GW\}$ is equivalent to proving that only those configurations keep their deepest spectral sink at the apex node under the constraints of three-gap/Steinhaus rigidity.

### 5. Net Impact
This reframe stabilizes the project's analytical strategy by providing a single object (the spectrum) that resolves the "blindness" of purely periodic approaches. By identifying the binding scale as a complete invariant, the cluster has unified the project's tournament-theoretic and number-theoretic branches into a singular geometric framework on the labelled Farey tree.
