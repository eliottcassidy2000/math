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

### **2. Distinguishing the High-Height Tail**
The Tail45 Separator is the "knife" that cuts away unphysical abstract solutions that previous moment-based filters (like $|W_1|+|W_2|$) could not reach.

*   **Excluding Cheap Directions:** The abstract "cheap" atom-cone moves ($r=1..5$) from HYP-2721, which can move the origin atom $q_0$ while preserving low moments, all fall **outside** the Tail45 strip:
    -   Directions $r=2, 4, 5$ are **below** the floor ($\text{tail45} = -1$).
    -   Directions $r=1, 3$ are **above** the ceiling ($\text{tail45} \ge 1$).
*   **Result:** By proving the Tail45 strip lemma, the project can mathematically forbid these dangerous abstract directions without needing to exhaustively check every possible runner configuration.

### **3. Integration with the Delsarte LP (HYP-2726)**
The Tail45 Separator provides the **feasibility filter** that sits between the abstract Delsarte LP and the actual LRC row logic.

*   **Unified Proof Order:**
    1.  **Delsarte LP/Relation Code:** Define the global structure and extremality (consecutive = anti-MDS).
    2.  **Generated-Word Frontier (Tail45 Separator):** Forbid all atom-cone moves that do not satisfy the miss-zeta compatibility strip.
    3.  **Factorial Odd-L1 Envelope:** Use the surviving physical moves to bound the origin-atom error $q_0$ against the remaining margin.
*   **Moderate-Span Balanced Gap:** This closes the gap by showing that the "balanced-wide" configurations have physical atom distributions that are too "expensive" (in terms of tail tax) to ever exceed the sector-cap margin.

### **Impact on Coordination**
The coordination ledger (SHA 366778) has been updated to reflect **HYP-2728**. The "Multi-Block Carrier Margin" route is now structured as a **Filtered Atom-Cone proof**. The Tail45 strip is the official structural guardrail that prevents the analytic lemma from being "hollowed out" by unphysical abstract solutions. The project has also added a self-contained **Lean finite identity module** (`TournamentH7.LRCFactorialAtom`) to formally verify the boundary algebra.

## mac-mini update: HYP-2726 LRC Cover Bound as Delsarte LP

The latest push (SHA d27c) by **Eliott Cassidy** (mac-mini-2026-06-21-S11) provides a major unification of the project's structural leads by identifying the **LRC cover bound as a Delsarte Linear Programming (LP) bound**. This framework connects coding theory, the moment-LP, and the Krawtchouk basis.

### **1. Delsarte LP Formulation (HYP-2726)**
The **Delsarte LP** formulation proves that the Sector Route extremality is a coding-theoretic problem.

*   **Primal Variables:** The **distance distribution** $\pi_E$ of the relation code $\Lambda(E) = \{n : \sum n_i e_i = 0\}$. This represents the distribution of relation supports in the runner set $E$.
*   **Dual Variables:** A function $g(t)$ used to bound the measure ($\text{meas } S_7(E) \le L_y(E)$).
*   **Dual Krawtchouk-Nonnegativity:** The dual feasibility condition is that $g(t)$ expands in the binary **Krawtchouk basis** $K_j(t; 6)$ with all-nonnegative coefficients $c_j \ge 0$.

### **2. Unification of Structural Leads**
The Delsarte LP unifies three previously distinct branches of the proof:
*   **MDS/Arc Coding Lens (HYP-2724):** The relation code is the **Delsarte Scheme**. In this scheme, the consecutive set is the **anti-MDS** code (minimal distance $d=2$, densest in low-weight relations), making it the LP-tight extremal configuration. Sidon sets/arcs are the **MDS** members (LP-slack).
*   **Even-Krawtchouk Structure:** The non-negativity of $c_j$ explains why even-band relations are clean. For $k=8$, the structure is purely even-only ($K_0, K_2, K_4$); for $k > 8$, the robust fact remains Krawtchouk-nonnegativity even as odd coefficients appear.
*   **Moment-LP (HYP-2721):** The THM-534 moment-LP is identified as the Delsarte LP itself. The origin atom $Q_0$ (the primary discrepancy) is the $B_6$ coefficient in this transform.

### **3. Verification and Bounds (k=8–13)**
The Delsarte-positivity of the dual coefficients $c_j$ was verified exactly for $k=8$ through $k=13$.
*   **k=8:** $c = [1/16, 0, 1/40, 0, 3/80, 0, 0]$. This purely even-only structure is unique to $k=8$.
*   **k=9, 10:** $c = [1/12, 1/72, 1/36, 1/48, 0, 0, 0]$.
*   **k=11, 12, 13:** $c = [1/8, 1/24, 1/24, 0, 0, 0, 0]$.
*   **Conclusion:** The consecutive block saturates the Delsarte LP at every binding $k$, confirming it as the LP-tight extremal configuration.
... (existing entries continue byte-for-byte) ...
