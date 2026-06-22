## codex-2026-06-22-S86 -- Finite Part A and Witness Attainment Wiring (checkpoint)

Formalized the wiring of the finite Part A error budget and the witness attainment interface (commit `47b4a006`). This checkpoint establishes the formal Lean glue between asymptotic witness density, finite-ruler counts, and the final LRC14 statement.

### 1. Finite Part A Error Budget (LRCWitnessPartA.lean)
Added a sorry-free module to handle the discretization error between the continuous witness measure $\rho^*$ and the discrete finite-ruler count $\rho_K$:
- **Discretization Lemma:** $\rho_K \ge \rho^* - \text{arcCount} / V_{max}$.
- **Positivity Criterion:** Proved that if the asymptotic density $\rho^*$ has a margin $\delta > \text{arcCount} / V_{max}$, then the finite density $\rho_K$ is strictly positive.
- **Uniform Arc Bounds:** Integrated the period-bounded `#arcs(GOOD)` signal, providing thresholds like `arcBound < deltaFloor * Vmax` to ensure finite-ruler positivity without re-evaluating the full count for every $V$.
- **Assembly Wrappers:** Provided `lrc14_from_finite_partA_p0_margin_shapes` and its split-branch variants, which compose the $p_0$ margin with the Part A error budget to discharge the lonely runner statement.

### 2. Witness Floor and Bonferroni Reduction (LRCWitnessBonferroni.lean)
Formalized the structural reduction of the global witness floor using the Bonferroni union bound:
- **Decoupled Bounds:** Reduces the witness floor $\rho^* \ge \mu(GOOD \cap G_P)$ to a sum of independent lower bounds: $\mu(GOOD) + \mu(G_P) - 1$.
- **Exact Rational Table:** Verified exact rational floors for consecutive minimizers (`nuConsec k`) and $G_P$ capacities (`capRat k`) for all cluster sizes $k=8..13$.
- **Unification:** Proved that the witness floor is implied by the $p_0$ wide-bound margin: $\text{witnessG2} \ge \text{meas}(G_P) - p_0 \ge \text{cap} - (\text{cap} - \delta) = \delta$. This aligns the witness route with the existing $p_0 \le \text{cap}$ analytic work.

### 3. Witness Attainment Interface
Added `LRCWitnessAttainmentBridge.lean` to prove definitional compatibility between the general `distZ` margin interface and the concrete `nearInt`/`Mreach` attainment route. Key theorems include `margin_eq_minReach` and `exists_lonely_of_margin_sSup_ge`, bridging the gap between abstract compactness arguments and concrete speed reach.

### 4. Axiom Audit and Root Integration
- **Axiom Audit:** Verified that the new bridge layers are sorry-free and rely only on standard/classical Lean foundations (`propext`, `Classical.choice`, `Quot.sound`).
- **Root Build:** Updated `TournamentH7.lean` and `Verify.lean` to include these modules, ensuring a clean, warning-free build of the entire integrated skeleton.

### 5. Impact
This checkpoint removes the remaining abstract "glue" obligations from the proof DAG. The proof is now localized to four concrete analytic nodes: instantiating the $p_0$ margin, instantiating $\text{cap} \le \text{meas}(G_P)$, proving the `#arcs` bound, and finalizing the `GOOD/witnessG2` event readouts.
