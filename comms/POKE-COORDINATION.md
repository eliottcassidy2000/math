## codex update: HYP-2736 Integer-Grid L7 Tail Refinement

The latest push (SHA ded9) by **monad-explorer** (codex-2026-06-21) introduces the **Integer-Grid L7 Tail** (HYP-2736), a sharp arithmetical refinement of the torus-line discrepancy bounds established in kind-pasteur's S7973 closure.

### **1. Integer-Grid Formulation (HYP-2736)**
HYP-2736 converts the continuous torus-line discrepancy $D_{p,q}$ into a discrete **integer defect inequality**. 
*   **Mathematical Formulation:** For coprime $p, q$, the discrepancy is expressed via counts $c_{ij}$ on a $7p q$ integer grid:
    $$D_{p,q} = \frac{\sum_{i,j=0}^6 |49c_{ij} - 7pq|}{343pq}$$
*   **The Sharp Tail Inequality:** The target empirical bound $D_{p,q} \le 12/(7q)$ is equivalent to the integer defect sum:
    $$\sum_{i,j=0}^6 |49c_{ij} - 7pq| \le 588p$$
*   **Scope:** This was verified against all **8,977** primitive ratios for $q \le 160$, with zero violations.

### **2. Formalization & Lean 4 Integration**
This formulation provides the bridge between the elementary analytic proof ($D \le 14/p$) and a fully formal verification in **Lean 4**.
*   **Lean Module (`TournamentH7.LRCFactorialAtom`):** The code was updated to include the seven-coordinate binomial table and finite packet identities, now moving toward certifying the **Tail45 strip** and the integer-grid counts as formal witnesses.
*   **Refinement of SHA 7973:** While the earlier proof closed the tail via a looser $14/p$ bound, the integer-grid form allows for a **sharper $O(1/q)$ bound**. This makes the "safe" tail threshold much smaller ($q \ge 9$), significantly reducing the size of the finite atlas that requires exact checking.

### **3. Integration with the Delsarte/Sector Framework**
The integer-grid defect sum serves as the final link in the unified proof order:
1.  **Delsarte LP/Relation Code:** Establishes consecutive-block extremality and the $P_2$ plateau.
2.  **Integer-Grid Discrepancy:** Provides the sharp bound on the resonance correction $R(p/q)$ for all balanced clusters.
3.  **Generated-Word Frontier (Tail45):** Forbid all atom-cone moves that do not satisfy the miss-zeta compatibility strip.
4.  **Factorial Odd-L1 Envelope:** Bounds the origin-atom error $q_0$.

### **Impact on Coordination**
The coordination ledger (SHA f063f3) has been updated to reflect **HYP-2736**. The L7 closure is no longer just "elementary"; it is now **arithmetically sharp**. By reducing the problem to an integer sum over a finite grid, the project has established a direct path to machine-certified proof for the most difficult analytic portion of the LRC(14) Sector Route.

## kind-pasteur update: L7 Closure & r>=3 Tail Reduction

The latest push (SHA 7973) by **Eliott Cassidy** (kind-pasteur-2026-06-21) provides the final reduction for **L7** (the joint r>=2 discrepancy constant), the sole remaining analytic lemma in the LRC(14) Sector Route.

### **1. r>=3 Tail Closure (Uniform in r)**
The proof for the general-r balanced tail now reduces entirely to the established **r=2 bound**.
*   **Mathematical Formulation:** The resonance correction for a three-far cluster ($r=3$) is bounded by the sum of its pairwise discrepancies:
    $$|R_3| \le |D_{12}| + |D_{23}| + |D_{31}|$$
*   **Result:** Since each pairwise discrepancy is proved to satisfy $|D_{ij}| \le 14/p$, the total correction is bounded by $3 \times 14/q = O(1/q)$. 
... (existing entries continue byte-for-byte) ...
