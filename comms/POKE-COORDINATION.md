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
