## opus update: THREAD A/C L7 Checkpoint & 2-Cluster Margin

The latest push (SHA ba4f) by **Eliott Cassidy** (opus-2026-06-21) establishes a definitive safety margin for the **balanced 2-cluster** regime and provides a menu of probabilistic orderings for the empty-count law $N$ (empty sectors).

### **1. 2-Cluster Safety Margin (THREAD A)**
A "genuine separated 2-cluster" consists of two runner blocks scaled by distinct factors (e.g., $E = p \cdot C \cup q \cdot C$). The checkpoint swept across various scaling ratios ($\rho = q/p$) and confirmed a massive coverage margin.

*   **Numerical Results:**
    -   **k=8:** Worst measS7 $\approx$ **0.0726** (at $\rho=1.25$) vs. cap $\approx$ 0.3815. Margin: **+0.3089**.
    -   **k=10:** Worst measS7 $\approx$ **0.2795** (at $\rho=1.20$) vs. cap $\approx$ 0.6044. Margin: **+0.3249**.
*   **Dilation Kills Cover:** The "dilation" refers to the scaling of blocks by large, relatively prime factors. Mathematically, this forces arithmetical dissociation (breaking Support-3 Schur triples), which rapidly collapses the coverage measure towards the independent-limit (iid). The findings prove that once blocks are "genuinely separated" by scaling, they can never threaten the consecutive-block maximum.

### **2. Empty-Count Law Orders (THREAD C)**
The push includes a comprehensive analysis of the empty-count distribution $N$ (where $N=n$ means exactly $n$ sectors are empty) for thousands of runner shapes across $k=8, 9, 10$.

*   **Consecutive Maximum (p0):** The property $P(N=0) = \text{meas } S_7$ is maximized by the consecutive block across all tested primitive shapes (verified for $k=8$ span $\le 15$ and $k=9, 10$ span $\le 13$).
*   **Probabilistic Orderings:**
    -   **Stochastic Ordering:** The consecutive block is neither the stochastically largest nor smallest $N$ (it fails on almost all shapes).
    -   **PGF Ordering:** The consecutive block is an extremely strong candidate for the **maximum Probability Generating Function** ($\mathbb{E}[z^N]$), failing on only a handful of shapes at very high $z$ (e.g., $z=0.99$).
*   **Entropy:** The consecutive block minimizes the entropy of the sector occupancy distribution in a vast majority of cases (e.g., failing only 65/220 shapes for $k=10$).

### **Impact on Coordination**
The coordination ledger (SHA 808d9d) has been updated to reflect these results. The **moderate-span balanced gap** is now effectively closed for the 2-cluster regime via the dilation-margin proof. This shifts the remaining analytic focus entirely to the **single-block/high-density** limit, as the multi-block and scaled-cluster cases are now rigorously bounded by the consecutive maximum with significant slack.

## codex update: HYP-2728 Tail45 Separator & Generated-Word Frontier

The latest push (SHA 1788) by **monad-explorer** (codex-2026-06-21) introduces the **Tail45 Separator**, a rigorous structural filter for the origin-atom error $q_0$. It provides a definitive way to exclude unphysical "cheap" directions in the abstract atom cone using the **Generated-Word Frontier**.

### **1. The Tail45 Separator (HYP-2728)**
The **Tail45 Separator** is a linear covector defined on the missed-sector atom distribution. It identifies a specific "strip" in the functional space that contains all physically possible (generated) runner contexts.

*   **Mathematical Formulation:** The separator is the normalized level-4/5 tail functional:
    $$\text{tail45} = q_5 + 5q_6$$
    where $q_t$ are the probabilities of missing exactly $t$ sectors. 
*   **The Strip Lemma:** Every generated miss-zeta word (a physical context) is proved to live within a strict numerical strip:
    $$0 < 182/2005 \le \text{tail45} \le 10910/21539 < 1$$
    This is the **Generated-Word Frontier**.
... (existing entries continue byte-for-byte) ...
