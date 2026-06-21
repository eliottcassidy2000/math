## codex update: HYP-2750 Signed Tanner Audit Renumbering

The latest push (SHA 3a3f) by **monad-claudebox** (codex-2026-06-21-S74) renumbers the **Signed Tanner/Dessin Audit** (HYP-2750), a critical structural guardrail for the Delsarte dual feasibility proofs.

### **1. Signed Tanner/Dessin Audit (HYP-2750)**
HYP-2750 establishes that the Tanner-like carriers for the THM-534 Delsarte duals are **signed dessins** (bipartite graphs on a surface) with Ferrers/equitable quotients, not sparse LDPC-like graphs or weakly regular unit-distance graphs.
*   **The Audit:** For the binding duals ($k=8 \dots 13$), the carrier graph is defined by moment rows (checks) and missed-depth atoms. The audit confirmed that while the graph has a clean degree quotient, its girth and expansion do not drive the bound.
*   **Result:** The "honestly negative" finding is that the unsigned graph structure forgets the sign/orientation information required for the Delsarte dominance predicate ($g(t) \ge scale \cdot [t=0]$). The useful invariant is the **half-arc sign orbit**, where support automorphisms never mix positive and negative edge classes.

### **2. Relation to Delsarte Weights and Tanner Negatives**
This renumbering consolidates the findings from the **Tanner honest negatives** session:
*   **Expansion vs. Sign:** Unsigned expansion (spectral gap, girth) is not the source of extremality. Instead, the signed weight distribution of the Delsarte quasicode drives the consecutive-set bound.
*   **Doyle-Holt Analogy:** The Doyle-Holt half-arc transitivity survives as a categorical level-up: the support has symmetries, but the orientation (sign) classes are not interchangeable, matching the non-self-converse nature of the extremal tournament.

### **3. Proof Order and Impact**
The audit confirms that the Belyi/dessin language serves as an **address layer** rather than a replacement for the analytic proof stack. It reinforces a strict proof order:
1. **Generated Depth Word**
2. **Signed Delsarte/Quasicode Parity**
3. **Aggregate Consec-Max Ledger**
4. **Cap Scalarization**

### **Impact on Coordination**
The coordination ledger (SHA 406cb1) has been updated to reflect **HYP-2750**. This establishes the Delsarte-dual carrier as a signed combinatorial object, preventing any lossy reduction to unsigned sparse graph theory. The project focus remains on the **full-residue stratum** (HYP-2749) and the machine-certified **Delsarte dual instances** (SHA 1797) as the primary pillars of the LRC(14) proof.

## mac-mini update: HYP-2749 Stratum-Localization

The latest push (SHA e440) by **Eliott Cassidy** (mac-mini-2026-06-21-S14) introduces **Stratum-Localization** (HYP-2749), a powerful reduction technique that narrows the search space for the consecutive-maximum proof to a specific arithmetic stratum.
... (existing entries continue byte-for-byte) ...
