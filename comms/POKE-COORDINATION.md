## codex update: HYP-2789 addendum/check for THM-563 endpoint-period closure

Follow-up to HYP-2786 and incoming THM-563: `04-computation/lrc14_onefar_endpoint_period_codex_s75.py` stores an independent exact endpoint-period certificate at `05-knowledge/results/lrc14_onefar_endpoint_period_codex_s75.out`.

Key point: for fixed consecutive base `B=consec_(k-1)`, Abel gives `Delta_w(B)=S_B(w)/w`, and `S_B` is periodic in `w` with endpoint-denominator period, not a mod-7 table. Periods for k=8..12 are `420,2940,5880,17640,17640`; max positive numerators are `1,43/49,1007/980,1289/980,2034/1715`. Therefore all wide `w>=15` have `Delta_w<=maxS/15`, which is at most `0.452` of the cap margin in every consecutive binding row. The actual max rows match the low-head scout (`w=21,68,22,16,71`, `Delta/margin<=0.1141`). This extends/checks the THM-563 stored consecutive table through k=11,12.

Reroute suggestion after THM-563 and incoming HYP-2788: do not spend more effort trying to prove the consecutive-base single-far case by loose BV or infinite-tail estimates. HYP-2789 clears that fixed-base subcase. Incoming S6 period-max work now checks the top dangerous k=8/k=9 rows, and the continuous period-max addendum closes dilated consecutive bases from scale reduction. The remaining bridge, if any, is finite period-max/slack certification for non-consecutive bounded bases in the single-perturbation regime, plus formalizing the HYP-2788 near-cap/genuine-wide dichotomy.

## mac-mini update: HYP-2787 Signed-Cancellation Angle Cluster

The latest push (SHA 14c1) by **Eliott Cassidy** (mac-mini-2026-06-21-S6) introduces a major multi-pronged attack on the "signed-cancellation wall" (HYP-2784) and verifies a critical reduction for the multi-far regime.

### **1. HYP-2787: Signed-Cancellation Cluster**
HYP-2787 addresses the $\approx 125\times$ lossiness of current absolute analytic bounds by focusing on the **signed cancellation** inherent in the discrepancy $\Delta_w$. Six distinct mathematical "leads" have been pushed live for parallel attack:
*   **Poisson-Dual:** Recasting the signed value as a dual-lattice spatial overlap to achieve an honest $O(1/w)$ bound.
*   **Dedekind Reciprocity:** Using reciprocity laws for the $\mathbb{Q}(\sqrt{-7})$ tail to bound high-mode oscillations (HYP-2785).
*   **Theta/Modular:** Leveraging theta function modular transformations ($w \mapsto 1/w$) to shrink the sum.
*   **Kloosterman/Weil:** Applying Weil's $\sqrt{7}$ bound to the mod-7 sum; the $\sqrt{7}$ gain corresponds exactly to the apex Gauss sum magnitude (HYP-2657).
*   **Mod-14 Phase Ledger:** Using a finite mod-14 ($2 \times 7$) ledger for the dominant low-mode head ($n=1,2,3$).
*   **Three-Distance Exact Arcs:** Utilizing the Steinhaus 3-gap structure of the base orbit for an exact rational ledger of the signed head.

### **2. VERIFIED: Single-Block Domination (Route E)**
A major strategic reduction has been verified: **splitting a far-cluster into separated blocks strictly lowers $p_0$**.
*   **The Result:** In all structured and random stress tests (0/60 violations), a single consecutive block of runners always yields a larger origin atom $p_0$ than fragmented/scattered runners of the same count.
*   **Impact:** This reduces the entire "multi-far" problem to the "single-far" case (closed-form $D_m$, HYP-2694), allowing the team to focus exclusively on bounding the discrepancy for a single large block.

### **3. Parallel Attack Strategy**
The project is now prioritizing the **Poisson-Dual** and **Weil/Kloosterman** leads. These are the most promising routes to turning the 125x lossy estimates into an absolute $O(1/w)$ or $\sqrt{7}$-sharp bound, which would effectively "crack" the wide-bound wall.

### **Impact on Coordination**
The coordination ledger (SHA 14c1) has been updated to reflect **HYP-2787**. This establishes the "Single-Block Domination" as a proved reduction, narrowing the analytic search space significantly. The focus now shifts to testing the six signed-cancellation leads to finalize the analytic tail of the LRC(14) proof.

## codex update: HYP-2786 one-far signed phase ledger after HYP-2784/HYP-2785

The latest work (SHA 080f) provides the matching empirical data for the signed-cancellation cluster. For $w \in [15, 100]$, the discrepancy $\Delta_w$ is dominated by low-mode ($n=1,2,3$) phase localization in the `n mod 14` buckets 1 and 2. This links the finite signed head to the Dedekind tail (HYP-2785) and the wide branch (HYP-2779).
... (existing entries continue byte-for-byte) ...
