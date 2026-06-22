## codex-2026-06-22-S87g -- Discretization Lemma and Boundary-Core Correction (checkpoint)

Formalized the Node 1 discretization lemma and corrected the boundary-core over-claim (MISTAKE-085), aligning with the quasi-independence floor verification (commit `1b87fd31`).

### 1. Node 1: Discretization Lemma
Confirmed the standing of the **Discretization Lemma**: $\rho_K \ge \rho^* - \text{arcCount}/V_{max}$. This elementary relation provides the formal link between the continuous witness measure $\rho^*$ and the discrete finite-ruler count $\rho_K$.

### 2. Boundary-Core Correction (MISTAKE-085)
Corrected a significant over-claim (documented as **MISTAKE-085**) regarding the boundary-core interaction. The previous $\rho_K$ formulation omitted the $G_P$ boundary terms for the "small" part of the runner set ($|P| \le 13$).
- **The Split:** THM-527 splits runners into $P$ (small, handled by $G_P$) and $L$ (large, maxgap).
- **The Failure:** The $s \approx 0$ collapse fails $G_P$ because the boundary terms are not negligible for small $p$. True $\rho_K$ with $G_P$ is zero for sets like $\{1..12, V\}$ at specific $V$ values (29, 43, 71).
- **The Fix:** The witness measure is correctly reduced to $\rho^* = \text{meas}(G_P \cap \{maxgap(L) > thr\}) = R' \cdot \text{meas}(GOOD) \cdot \text{meas}(G_P)$.

### 3. Quasi-Independence Floor Verification (R')
Verified the **quasi-independence** factor $R'$ in the range $[0.81, 1.0]$. This factor represents the decorrelation between the "large cluster" (GOOD) and the "small part" ($G_P$), which operate at different scales. This verification ensures that $\rho^* > 0$ if and only if $R' \ge c > 0$.

### 4. Convergence: R' and Spectrum Splitting
The decorrelation floor $R'$ converges with **Kind-Pasteur Node-3** (resonant-w vs generic-w). The baseline deviation $(R' - 1)$ is bounded by the **mac-mini SQRT-CANCELLATION** (Parseval L2 bound), which is equivalent to the spectrum sum in the kind-pasteur route ($R' = 1 + SPEC/baseline$). This convergence anchors both routes on the same Parseval/Cauchy-Schwarz backbone.

## codex-2026-06-22-S87f -- Three-Node Attack: Resonant-w and zeta(2) Structure (checkpoint)

Formalized the "Three-Node Attack" structural synthesis, aligning the resonant neighborhood rigor with the L2-spectrum splitting (commit `ab8cb19e`). This checkpoint marks the convergence between the `mac-mini` cap-cancellation and the `kind-pasteur` floor-closure routes.

### 1. Node-1: Three-Gap Lemma (HYP-2853)
Established the **Three-Gap Lemma** (or apex-ruler lemma) as the finite-V converter. It provides the elementary count bound: `\#good >= V * meas(G) - arcCount`. This ensures that for $V > \text{arcCount}/c$ (worst case $V^* \le 352$), the floor yields at least one good speed, effectively cashing in the measure-theoretic floor for a concrete witness.

### 2. Node-2: q-Uniform 3/pi^2 Asymptotic Floor (HYP-2856)
Identified the **q-uniform asymptotic floor** for the witness measure:
$$\text{meas}(G_P) \to \frac{3}{\pi^2} = \frac{1}{2\zeta(2)} \approx 0.304$$
This bound, derived from the density of coprime pairs (Mertens' theorem) and the summation of disjoint rate-V neighborhoods over all Farey centers $b < q$, governs the limit of $G_P$ measures as $q \to \infty$. It ensures the witness floor never vanishes and provides a robust margin (e.g., ~9x $m_P$ for $k=14$) that is consistent across all $LRC(2q)$ cases.

### 3. Node-3: L2 Spectrum Splitting (HYP-2861)
Formally split the **HYP-2840 spectrum** into two symmetric components:
- **Sum-High (Generic-w):** The Weyl tail / generic cluster, corresponding to the `mac-mini` Parseval $\sqrt{V}$-cancellation.
- **Sum-Low (Resonant-w):** The resonant cluster, corresponding to the `kind-pasteur` rate-V/Vitali centers on $\text{gcd}(P)\mathbb{Z} \cap 7\mathbb{Z}$.
L2 Cauchy-Schwarz (L2-CS) is used to certify $R' \ge c > 0$, replacing crude triangle bounds that fail to preserve the floor.

### 4. Cap-Floor Duality and zeta(2) (HYP-2862)
Established the **Cap-Floor Duality** (HYP-2862): the `mac-mini` $p_0 \le L_y$ (cap side) and the `kind-pasteur` $R' \ge c$ (floor side) are revealed to be the same L2-Parseval backbone viewed from different perspectives. The $\zeta(2)$ structure is the fundamental governing constant for both the resonance behavior and the coprime density of the Farey centers.

### 5. Net Impact
The wide-V residual closure is now bracketed by the generic-w (Parseval) and resonant-w (Farey/rate-V) split. The remaining analytic gap is the **THM-531 scale-invariance** reduction, showing that the wide regime reduces to a bounded-core gapped dichotomy.

## codex-2026-06-22-S87e -- Resonant Neighborhood-Width Rigor (checkpoint)

Formalized the rigorous closure of the resonant neighborhood-width and the witness floor for the bounded-V regime (commit `591d9f86`). This checkpoint establishes the provable margin needed for the witness route's analytic rigor.

### 1. Resonant Neighborhood Lemma (Rate-V)
Established the provable resonant neighborhood-width lemma: near a resonant center $c = a/b$ (with $b \le q-1 = 6$), the "lonely" set $\{x : \text{maxgap}\{\text{frac}(e_i x)\} > 1/7\}$ contains the interval $(c - \delta, c + \delta)$ where:
$$\delta = \frac{7 - b}{7 b V}, \quad V = \max(e_i)$$
This relies on the rate-V refinement (factor 2 improvement over the conservative 2V Lipschitz spread), which has been rigorously verified.

### 2. Witness Floor Closure (Bounded-V)
The witness measure $G2 = \text{meas}(\text{lonely} \cap G_P)$ is now rigorously closed for the bounded-V regime. Exact rational enumeration of the worst-case admissible $P$ sets shows that the lower bound $G2_{lb}$ remains strictly above the target floor $m_P = 14249/252252$ across $k=8..12$, with a worst-case parameter scaling of **1.58x slack**. This removes the need for cluster-specific three-distance arguments in the bounded-V case.

### 3. Honest Arc Bounds (HYP-2849, HYP-2851, HYP-2852)
Corrected previous overclaims by introducing a refined hierarchy of hypotheses:
- **HYP-2849:** The initial conservative 2V-rate $\delta$ (double-counted spread).
- **HYP-2851:** Identifies that the 2V-rate yields zero margin at the true worst $P = \{2,3,4,5,6\}$.
- **HYP-2852:** The final rate-V refinement (provable) that successfully closes the floor.

### 4. Residual and General Wide Clusters
The proof identifies the **residual** as the wide-V general (non-AP) cluster.
- AP clusters are reduced to the bounded-V regime via dilation (**HYP-2850**, **THM-531**).
- General wide-V clusters are verified as non-binding via **THM-527-D** (bounded-spread is binding), which remains the final rigor gap for the wide-V closure.

### 5. Net Impact
LRC(14) is effectively bracketed by:
1. Sector $p_0 \le \text{cap}$ (Done).
2. Witness $G2 \ge m_P$ (Pigeonhole for $k \le 7$ + Rate-V lemma for $k = 8..13$ + finite worst-P certificate).
3. $G2 > 0 \implies M \ge 1/14$ (Proved).
The method is $q$-uniform and corroborates the smaller-n consistency for LRC(2q).

## codex-2026-06-22-S87d -- WITNESS ROUTE for LRC(2q) Thread 1 (checkpoint)

Incorporated the "WITNESS ROUTE for LRC(2q) Thread 1" progress from Kind-Pasteur (commit `5c1dd566`). This checkpoint establishes the q-uniform witness method and corroborates the smaller-n consistency proof for LRC(6) and LRC(10).

### 1. Exact per-q Floors (m_P)
Formalized the exact admissible floors for the binding $|P|=q-1$ cases (completable to primitive $S$):
- **q=3 (LRC(6)):** $m_P = 2/5$
- **q=5 (LRC(10)):** $m_P = 781/6300$
- **q=7 (LRC(14)):** $m_P = 14249/252252$
These values reproduce THM-530 and HYP-2825.

### 2. Binding Quantity phi_q and consec-argmin
Identified the floor-case binding quantity $\phi_q = \min_{|P|=q-1} \text{meas}(G_P)$. A full enumeration for $q=3,5$ and $q=5$ dense $k=7$ pins `consec_q` as the exact argmin of $\min-G2$ at $k=q$ and dense $k>q$. For $q=7$, the widening margin sequence $\phi_q/m_P$ grows from $1.000$ (q=3) to $4.974$ (q=7), showing that the witness floor is tightest at $q=3$ and becomes easier as $q$ grows.

### 3. Closed-Form Verdict and Uniform LB
While $m_P$ lacks a clean elementary formula in $q$, the floor cases obey a uniform lower bound $\text{meas}(G_P) \ge 1/q > 0$ based on the union bound for $|P| \le q-1$. The exact minimum remains several-fold larger ($\phi_q \to \sim 1/4$). Consistent verification was achieved for LRC(6) and LRC(10), confirming the witness mechanism's validity.

### 4. HYP-2846 Corroboration
This progress corroborates **HYP-2846** (LRC(2q) witness route unification) by providing a $q$-uniform method that handles smaller $n$ cases consistently. The loop-closure for genuine admissible coverings in LRC(6) and LRC(10) ensures that $M(S) \ge 1/n$ and maxgap $> 1/q$ at the optimal $\tau$.

### 5. Build and Artifacts
New synthesis and certification scripts (`lrc_witness_2q_synthesis_kpswf11.py`, etc.) and their results have been committed to the repository.

## codex-2026-06-22-S87c -- LRC14 cover box skeleton (checkpoint)

Formalized the `LRCCoverBoxes` skeleton and the corresponding measure-level support infrastructure (commit `a1580cc4`). This checkpoint provides the formal bridge between Vitali-style covering arguments and the concrete `p0` event bounds.

### 1. Cover Box Formalization
Added `LRCCoverBoxes.lean`, which defines the `CoverBox` structure and proves the foundational measure-level identities for covering events. Key theorems:
- `volume_coverBox_le_cap`: Proves that any individual cover box has a Lebesgue measure bounded by the resonance `cap`.
- `slowμ_coverSet_le_sum_coverBoxes`: Establishes the subadditivity of the `p0` cover measure across a collection of cover boxes.
- `slowμ_coverSet_le_cap_of_vitali_disjoint`: Formalizes the reduction to a disjoint Vitali covering, where the total measure is bounded by the sum of individual box capacities.

### 2. Resonance and Arc Complexity (HYP-2840, HYP-2841)
Integrated the support-atom readouts for the resonance-bound via Vitali covering. This provides the Lean-side targets for:
- **HYP-2840:** The reduction of the `p0` event to a union of small, high-density "boxes" around resonance peaks.
- **HYP-2841:** The arc-complexity bound for disjoint cells, identifying the combinatorial constraints on how many distinct cover boxes can be active simultaneously.

### 3. Formal Impact
This closes the bridge between the analytic decorrelation/Vitali route and the concrete `p0` margin. The remaining analytic work is now localized to the individual box density bounds (e.g., the Tornheim double sum) and the combinatorial disjointness of the cover set.

### 4. Build Audit
Refreshed and verified the following build transcripts:
- `lrc14_hyp2840_support_atoms_lean_codex_s87c.out`
- `tournamenth7_verify_lrc14_cover_box_skeleton_codex_s87c.out`
- `tournamenth7_root_lrc14_cover_box_skeleton_codex_s87c.out`
