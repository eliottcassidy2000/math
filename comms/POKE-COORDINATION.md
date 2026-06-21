## monad-explorer update: HYP-2717 Top-Character Carrier Relation Filter

The latest push (SHA 03cc) by **monad-explorer** (codex-2026-06-21) introduces a critical structural filter for the moderate-span multi-block regime, addressing the "Carrier-Product Bound" gap by accounting for exact carrier resonances.

### **1. Top-Character Carrier Relation Filter (HYP-2717)**
The **Carrier Relation Filter** is an analytic tool designed to bridge the gap between the idealized "decoupled" carrier product and the actual arithmetical coverage of a multi-block row. It expands the coverage probability in carrier Fourier modes and splits the error into two distinct classes based on integer relations between the carrier speeds $M_i$.

*   **Mathematical Formulation:** For a row $E = \{0\} \cup \bigcup_i (M_i + B_i)$, the coverage error $p_0(E) - \text{Product}(E)$ is expanded in carrier characters:
    $$p_0(E) - \text{Product}(E) = \sum_{n \neq 0} \int_0^1 \hat{F}_n(x) \exp(2\pi i (n \cdot M)x) dx$$
    The filter partitions these modes into:
    1.  **Exact Carrier Relations ($\Lambda(M) = \{n : n \cdot M = 0\}$):** These modes survive the line integral regardless of how large the speeds $M$ are. They are the unavoidable "resonances" of the system.
    2.  **Nonresonant Carrier Modes ($n \cdot M \neq 0$):** These modes oscillate and are suppressed by $1/|n \cdot M|$ through integration by parts or bounded-variation estimates.

### **2. Refining the Carrier-Product Bound ($W \approx 26$)**
This filter sharpens the **Carrier-Product Bound** (HYP-2713/HYP-2714) by replacing the assumption of full-torus equidistribution with a rigorous **relation-lattice tail bound**.

*   **Filter Logic:** It proves that while exact relations are unavoidable when multiple blocks are present ($g \ge 2$), they only contribute significantly when their "height" is small. High-height relations have small Fourier coefficients $\hat{F}_n$ and thus do not consume the available proof budget.
*   **W ≈ 26 Regime:** For gapped rows at this span, the filter allows the proof to route low-height resonances and small-denominator nonresonances to a **finite ledger** (HYP-2714), while the analytic lemma (the "comfortable budget" multi-block lemma) only needs to handle the high-height/high-denominator tail.

### **3. Targeting the Moderate-Span Balanced Gap**
The filter provides the missing mechanism to handle the "moderate-span balanced" regime (span 15–300), where gapped "peeling" is too weak.

*   **The "Comfortable" Slack:** The filter leverages the 0.18 slack identified in mac-mini-S8. It shows that the "Top Character" ($\mathcal{M}_6$) from HYP-2716 is the only scalar that must fit under the `cap-product` margin.
*   **Convergence:** By proving that the sum of the exact-relation tail and the nonresonant oscillation stays below the available 0.18 budget, the project can close the multi-block gap without needing an impossibly sharp bound.

### **Impact on Coordination**
This update shifts the "Multi-Block Carrier Margin" from an opaque assumption to a **filtered Fourier problem**. It identifies that the difficulty of the lonely runner problem at $n=14$ is not in the aggregate cover itself, but in the specific **low-height integer relations** between the block speeds. The coordination ledger now reflects that the final proof step is a **signed tail estimate** on the carrier relation lattice.

## Eliott Cassidy update: LRC(14) Sector Route Assembled

The latest push (SHA 585dad) by **Eliott Cassidy** (mac-mini-S8) represents a major structural milestone in the LRC(14) proof, assembling the final **Sector Route** and localizing the remaining gap to a single analytic lemma with a comfortable margin.

### **1. The Four Closed Regimes**
The global proof that the consecutive block `consec_k` maximizes the coverage measure ($\text{meas } S_7$) is now partitioned into four distinct regimes, all of which are considered **closed**:

*   **Bounded Small Span (≤14):** Confirmed via finite exact checks. For $k=12$ and $k=13$, the i.i.d. surjection rate already exceeds the cap, meaning no "wide" regime danger exists.
*   **Single Far Point:** This is the *universal wide supremum*. Since the coverage measure is monotone decreasing in the number of far points, a single far point is the worst-case scenario. It is closed by a "comb bound" ($|\Delta_w| \le 2c_1(E')/7w$), providing a finite verification window ($W^* \le 48$) with margins $\ge 0.12$.
*   **Multiple Far Points / Fully Dissociated:** These cases are strictly easier than the single-far point, as coverage decreases with each added far point, eventually bottoming at the i.i.d. surjection rate which is much lower than the Sector Cap ($cap_k$).
*   **Consecutive Block Tabulation:** The consecutive block itself is verified to satisfy the Sector Cap ($cap_k$) for all relevant $k$, with a remaining slack of 0.023.

### **2. The Comfortable-Budget Multi-Block Lemma**
The remaining gap in the proof has been localized to the **Moderate-Span Balanced** regime (gcd-1 shapes of span between ~15 and a few hundred with no large gaps). 

*   **Localization:** In this regime, the single-far point "peeling" is too weak, and the number of configurations is too large to enumerate. The proof relies on the **Carrier-Product Bound**: well-separated sub-clusters decouple their colorings, and because "splitting strictly costs cover" (codex's Route E), the single block dominates.
*   **The Lemma:** The final proof obligation is to make the finite-separation carrier error rigorous using a multi-dimensional Weyl / Erdős–Turán estimate. 
*   **The "Comfortable Budget":** Crucially, the worst moderate-span shape sampled covers only 0.197 against a cap of 0.381, leaving a **slack of 0.184**. This wide margin means the remaining lemma does not need to be sharp; a relatively crude bound will suffice to close the proof.

### **3. HYP-2713 & HYP-2714**
These hypotheses formalize the assembly of the sector route:
*   **HYP-2713 (Sector Route Assembly):** Bridges the local and global domains by identifying the consecutive block as the global extremizer for coverage measure across all True-Wide regimes.
*   **HYP-2714 (Multi-Block Carrier Margin):** Targets the localized gap in the moderate-span balanced regime. It asserts that the carrier error for multiple sub-clusters is contained within the available 0.18 margin, ensuring the consecutive block remains the global maximum.

### **Status: Finish Line in Sight**
While the Lonely Runner Conjecture at $n=14$ is not yet fully proved, the problem has been reduced from an opaque, irreducible "wall" to a single analytic lemma. The project has moved from a research frontier to a localized "finish line," where the remaining task is to rigorously bound the multi-block error within a generous margin.
... (existing entries continue byte-for-byte) ...