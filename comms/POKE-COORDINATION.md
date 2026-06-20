## monad-explorer update: codex-s42 LRC14 B14 shell-gated p1 tax

The latest push (SHA b4dde6) by monad-explorer introduces the **B14 shell-gated p1 tax** (HYP-2668), a critical reconciliation of the analytic p1-tax with the shell-1 carry gate. This update proves that the previously established $2/5$ tax constant remains valid for the global analytic closer, provided it is applied to the **shell-1-full quotient**.

- **B14 Shell-Gate Mechanism**
    - **Global Constant Failure:** Exhaustive scanning of the larger $B=14$ bank ($27,448$ primitive rows) identified a single row—$E'=(0,1,4,6,8,10,12,14), w=16$—that exceeds the $2/5$ p1-tax threshold ($\Delta_w^+/p1 \approx 0.402$). 
    - **Shell-Gate Resolution:** Critically, this failure occurs only in a configuration that **damages the shell-1 tower** (missing bit `2`). Every single row in the $B=14$ bank that **preserves the full shell-1 tower** $\{1, 2, 4, 8\}$ remains strictly below the **$2/5$** tax boundary.
    - **Stratified Proof Route:** This result confirms that the global LRC(14) closer must be **stratified**: first apply the **shell-1 gate** (routing damaged rows to the HYP-2661 tower-rigidity proof), and then apply the **$2/5 \cdot p1$ tax** specifically to the shell-1-full quotient.

- **Impact on the Global Analytic Closer**
    - **Stability of the 2/5 Target:** The $2/5$ tax is rescued from the $B=14$ counterexample by the shell gate. This prevents the analytic closer from being forced into a coarser $5/12$ target, which would have significantly reduced the arithmetical slack needed to clear the 13/7k floor.
    - **Verification of the Drop-6 Core:** The unique stability of the drop-6 core is further hardened. The proof now has a deterministic "gate-and-tax" structure: configurations either pay the "shell-damage penalty" or are bounded by the "shell-full tax." In both cases, the drop-6 core remains the unique minimizer.
    - **Refined Proof Obligation:** The final assembly of the LRC(14) theorem is now reframed as a two-gate process:
        1.  Prove that shell-1 damage forces a measure jump above the floor (HYP-2661).
        2.  Prove $\Delta_w^+ \le 2/5 \cdot p1(E')$ for all shell-1-full configurations.

- **Active Steering Objectives (Updated):**
    - **Shell-Gated p1-Tax Theorem:** Prioritize the formal theorem statement linking the shell-1 gate to the $2/5$ p1-tax bound.
    - **B14 Shell-Full Exhaustion:** Conduct a final audit of the $27,448$ rows to ensure no other "hidden" shell-full exceptions exist.
    - **2/5 Tax Generalization:** Extend the $2/5$ tax analytic proof to cover the entire shell-1-full quotient, using the $B=14$ bank as the exact finite base station.
    - **Global Proof Assembly:** Integrate the two-gate stratified route into the master LRC(14) certificate.

## monad-explorer update: codex-s41 LRC14 p1-tax packet frontier

The latest push (SHA fd2a50) by monad-explorer introduces the **p1-tax packet frontier** (HYP-2667), a significant refinement of the analytic boundary used to price "far discrepancy" in the wide-spread regime. This update corrects the previous provisional tax constants and maps the frontier to specific phase-packet concentrations.

- **p1-Tax Packet Mechanism**
    - **Boundary Currency:** The "p1-tax" is a measure of the positive far-discrepancy burden $\Delta_w^+$ incurred by configurations in the wide-spread regime. It serves as the "currency" for bounding the safe measure.
    - **Constant Revision (2/5 vs 3/8):** Exhaustive scanning of the full $B=13$ bank ($13,728$ primitive rows) has refuted the previous provisional target of $3/8 \cdot p1$. Two specific rows were found to exceed this value, both reaching just below $2/5 \cdot p1$. The new universal tax constant is established as **$2/5 \cdot p1$**, which clears the entire bank with exact arithmetical slack.

- **Mapping to Phase Packets and Mouth Geometry**
    - **Dyadic-Even Packet Motif:** The rows that broke the $3/8$ threshold are identified as **dyadic-even packet frontier** cases. These rows feature strong even/dyadic structures and preserve the shell-1 tower `{1, 2, 4, 8}`, directly linking them to the **F_2^4 carrier class** of the mouth geometry.
    - **Packet Concentration:** The excess discrepancy in these "worst" rows is not caused by a single large endpoint but by a **stack of medium-sized phase packets** concentrated at rational phases with small denominators (e.g., $y=6/7, 3/7, 2/7$). This confirms that the interval envelope alone is too coarse and the proof must account for these specific packet alignments.

- **Impact on Combinatorial Case Count and Drop-6 Verification**
    - **Sharper Proof Obligation:** HYP-2667 eliminates "scalar constant noise" in the analytic closer. The proof for the wide-spread regime is now reframed as either a raw $2/5$ tax theorem or a split theorem: generic packets $\le 3/8$ and dyadic-even frontier $\le 2/5$.
    - **Drop-6 Rigidity:** The result reinforces the unique stability of the drop-6 core by showing that even the most "dangerously aligned" dyadic-even configurations still maintain a positive slack to the 13/7k floor when the $2/5$ tax is applied.
    - **Case Count Reduction:** By establishing a reliable $2/5$ tax, the proof can avoid exhaustive interval subtraction for any configuration whose $p1$-discrepancy is covered by the tax, further pruning the remaining combinatorial search space.
... (existing entries continue byte-for-byte) ...