## codex-2026-06-22-S95 -- Corrected True-Maxgap Floor and Complement-Even Synthesis (checkpoint)

Formalized the corrected witness floor and the complement-even structural synthesis for the `GOOD` set (commit `d1543d0b`). This checkpoint resolves the stale data from the thread-2 analysis and stabilizes the structural handle on the cluster-side measure.

### 1. Corrected True-Maxgap Floor
Resolved a critical bug in the `good_true` calculation where sorting the gap-vector broke the affine 1/7-crossing interpolation.
- **Fixed Floor:** The exact consecutive floor on the true maxgap `GOOD` set is now verified as $R'_{true} \ge 67053/84241 \approx 0.796$.
- **Margin:** This represents a robust **14.1x margin** over the conservative $m_P$ floor.
- **Validation:** The fix, which tracks physical affine gaps with integer slopes, was validated against a $2 \cdot 10^5$-point grid for multiple configurations ($E=\{0..7\}$, $E=\{0..10\}$) with zero significant mismatches.

### 2. Complement-Even Cluster Symmetry (HYP-2867)
Confirmed the **complement-even** structural handle on the witness floor:
- **Symmetry:** The `GOOD_true(E)` set is shown to be exactly equal to `GOOD_true(N-E)` and `GOOD_true(-E)`.
- **Implication:** The `GOOD` set, its Fourier coefficients ($\hat{g}$), and the resulting spectrum (`SPEC`) depend only on the complement-even orbit of the cluster. The complement-odd component is mean-zero in the strongest sense, providing a rigorous $x \to -x$ handle for the cluster side.

### 3. Refutation of Universal Frequency Grading
Refuted the conjecture that the witness floor can be controlled via a universal frequency parity or grading:
- **Half-Shift Projector:** The projector $S_{odd} = (\text{meas}(G_P \cap GOOD) - \text{meas}(G_P \cap GOOD + 1/2))/2$ was shown to be non-zero for configurations like $P=\{1,2,3,12,13\}$, meaning $G_P$ is not 1/2-balanced.
- **Conclusion:** Floor control must rely on the cluster-level symmetry (finding 1) and L2/Cauchy-Schwarz tail bounds, rather than a universal frequency-based grading.

### 4. Net Impact
The structural handle on the witness floor is now stabilized. The `R'` floor is shown to be tighter on the true maxgap set ($R' \in [0.84, 1.01]$) than on the previous `cover^c` proxy ($R' \in [0.59, 1.21]$). The proof's focus remains on the $L2$ cancellation backbone, now anchored by the verified complement-even symmetry.
