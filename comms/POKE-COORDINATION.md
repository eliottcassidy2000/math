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

- **Active Steering Objectives (Updated):**
    - **2/5 p1-Tax Theorem:** Formalize the analytic proof of the $2/5 \cdot p1$ far-discrepancy tax.
    - **Phase Packet Alignment Audit:** Conduct a targeted audit of the dyadic-even frontier to characterize the "medium-sized packet stack" mechanism.
    - **Split-Theorem Formalization:** Prepare a split theorem structure for generic vs. dyadic-even phase packets to minimize the tax burden in the global proof.
    - **B=13 Bank Closure:** Finalize the integration of the $B=13$ p1-tax bank results into the master LRC(14) proof.

## kps-S39 update: checkpoint LRC14 three-tail shell-1 frontier

The introduction of the **LRC14 three-tail shell-1 frontier** (HYP-2664, SHA 91a877) establishes a critical arithmetical gate for the three-tail replacement layer. By applying the **shell-1 carry conservation law** (HYP-2661) *before* exhaustive enumeration, the proof reduces the finite residue burden of the three-tail layer by over 50%.

- **Three-Tail Shell-1 Frontier Mechanism**
    - **Carry Gate Priority:** The frontier identifies that the vast majority (87/100) of the most difficult "crude comb" cases—those with high first-tail cutoffs—are configurations that damage the **shell-1 tower {1, 2, 4, 8}**. 
    - **Cutoff Reduction:** By applying the shell-1 gate (HYP-2661), the global maximum first-tail cutoff $R$ for the three-tail layer drops from **308** to **239**. This removes almost all of the "top frontier" cases that previously blocked the combinatorial proof.
    - **Burden Shift:** The total number of tail-triples requiring manual audit below the first cutoff is reduced from **4.19M** to **1.87M**, a significant gain in proof efficiency.

- **Stability of the Drop-6 Core**
    - **Shell-1 as Stability Anchor:** This update reinforces that the drop-6 core's stability is arithmetically anchored by the **full rank of its shell-1 carrier class**. The frontier data proves that any three-tail replacement attempting to undercut the $426/35035$ floor is almost certainly forced to preserve this tower, further narrowing the possible counter-configuration space.

- **Advancing the Combinatorial Tower-Deletion Proof**
    - **Proof Order Optimization:** HYP-2664 redefines the optimal path for the **combinatorial tower-deletion proof**. The proof is now structured as a sequence of quotients:
        1.  **Shell-1 Carry Gate:** Instantly filter packets that delete 1, 2, 4, or 8.
        2.  **Root-Packet Address:** Focus only on shell-1-full packets.
        3.  **Mouth-Owner Ledger:** Separate configurations by mouth-interval retention.
        4.  **Nested Comb:** Apply the finite residue check only as a final closer.
    - **Reframed Theorem Target:** The target for the tower-deletion proof is now reframed: any three-tail AP packet below $426/35035$ *must* be shell-1 full and lie within a small family of old-mouth/root-packet templates. This reduces the "OPEN" combinatorial problem to a bounded, addressed finite residue search.

- **Active Steering Objectives (Updated):**
    - **Shell-1 Deletion Theorem:** Prioritize the independent proof of HYP-2661 to fully unlock the carry gate for the three-tail layer.
    - **Shell-1-Full Nested Comb:** Run the exact nested comb analysis specifically on the shell-1-full packets identified in the frontier atlas.
    - **Mouth-Owner Classification:** Classify the remaining 1.87M triples by their mouth-survivor values to close the mouth-retention rigidity lemma.
    - **Worst-Case Inequality:** Develop a direct inequality for the remaining worst shell-1-full base `holes=(3,5,6,13)`.
... (existing entries continue byte-for-byte) ...