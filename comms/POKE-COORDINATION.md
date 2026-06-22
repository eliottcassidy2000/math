## codex-S114 -- Three-Mode Composition and Committed-Denominator Wall (checkpoint)

Formalized the corrected geometric composition of the three tournament recursion modes and established the "committed-denominator wall" guardrail for the LRC(14) proof (commit `8bc93582`). This checkpoint resolves structural misconceptions about the recursion modes and anchors the proof's split between finite Venn nodes and analytic equidistribution.

### 1. Corrected Three-Mode Composition
Refined the mapping between tournament recursion modes and their geometric/arithmetic counterparts:
- **Möbius Mode ($n \to n-1$):** Functions as the always-on inclusion-exclusion skeleton for all sizes.
- **Eisenstein Mode ($n \to n-3$):** Acts as the even-size fold ($2q \to q$) exposing the apex structure.
- **Legendre Mode ($n \to n-2$):** Functions as a full 3-set Venn diagram for odd sizes.
- **Correction:** Clarified that the Legendre mode is a geometric Venn over corners $A, D, B$ (sizes $N-1, N-2, N-1$) and their overlaps. This distinguishes "corner" from "overlap" labels (e.g., $D$ and $C$ both have size $N-2$ but distinct roles), which is critical for finite realizability in Node 2.

### 2. The Committed-Denominator Wall (LCM Family)
Established a rigorous proof-level guardrail against using fixed finite denominator bases for the LRC(14) closure.
- **Mechanism:** For any bound $X$, the family $S_X = \{1, \dots, 11, 13, \text{lcm}(2, \dots, X)\}$ "kills" every denominator $D \le X$ simultaneously, as the committed speed forces that runner to phase 0 mod 1.
- **Refutation:** This proves that the minimal witness denominator must tend to infinity with $X$, refuting every fixed-denominator atlas shortcut.
- **Opening Behavior:** Probed the first denominators above the wall; the "next-prime" rule is false. The first unblocked witnesses often occur at prime-powers (e.g., $121 = 11^2$ for $X=110, 120$) or specific residue openings.

### 3. Separation of Proof Nodes (2 vs 3)
The three-mode composition provides a clean handoff between the two final proof components:
- **Node 2 (Finite Venn):** Handles the cap/extremality statement within the corrected Legendre Venn geometry, focusing on AP/three-gap rigidity and hull-labels.
- **Node 3 / Part A (Analytic):** Handles the analytic floor through effective equidistribution beyond the committed wall, using exact-period unit packets with prime-power labels.

### 4. Structural Hierarchy
Established the unified composition path:
1. Exact-period packets.
2. Möbius divisor/IE labels ($\phi = \mu \ast id$).
3. Eisenstein even fold ($2q \to q$).
4. Legendre odd Venn (at $q$).
5. Scalar cap/floor comparison (only after labels survive).

### 5. Net Impact
This checkpoint stabilizes the project's proof architecture by replacing abstract scalar recurrences with a labelled geometric and arithmetic hierarchy. By isolating the "committed-denominator wall," it forces the final analytic closure to rely on effective equidistribution rather than heuristic denominator searches.
