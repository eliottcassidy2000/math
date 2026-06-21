## kind-pasteur update: HYP-2724 MDS/Arc Coding Lens & Support-3 Driver

The latest push (SHA 0e30) by **Eliott Cassidy** (kind-pasteur-2026-06-21) introduces the **MDS/Arc Coding Lens**, a major structural reframing of the LRC(14) proof that confirms the **Support-3 Driver** and adjudicates a critical recurring error in the support-6 floor logic.

### **1. MDS/Arc Coding Lens & Extremal Dichotomy (HYP-2724)**
The **MDS/Arc Coding Lens** treats the relation lattice $\Lambda(E) = \{n \in \mathbb{Z}^k : \sum n_i e_i = 0\}$ as a linear code $[k, k-1, d]$, where $d$ is the minimum support of a relation. This lens establishes a clear **extremal dichotomy** for the proof:

*   **Anti-MDS (The AP/Consecutive Set):** The consecutive set is the "anti-MDS" member of the family, where the minimum distance $d$ collapses to 2. It is the densest configuration in low-weight relations and is confirmed as the global argmax for the coverage measure.
*   **MDS (The Sidon/Arc Set):** The Sidon set is the "MDS" member, representing general position with no low-support relations. It is the easiest configuration to bound, with a correction $\text{corr} \approx 0$.
*   **Significance:** This provides a coding-theoretic mechanism for why the consecutive block is the unique extremizer of the Sector Route.

### **2. Confirmation of the Support-3 Driver**
The session confirmed that **Support-3 relations** (3-APs and Schur triples $a+b=c$) are the primary drivers of the carrier error.
*   **Correlation:** Statistical analysis shows a Pearson correlation of **+0.93** between the count of Support-3 relations ($A_3$) and the total carrier error ($\text{corr}$).
*   **Sign Seam:** It establishes a combinatorial sign seam by support size: Support-2 carries a small/mixed sign, Support-3 carries a large positive mass, and Support-4+ alternates.

### **3. Adjudication of the THM-538 Support-6 Trap**
A major result of this session was the definitive adjudication of the **THM-538/MISTAKE-078 trap**. 
*   **The Trap:** A workflow previously "refuted" the Support-3 framing by claiming a "Support-6 floor" where all relations of support $\le 5$ vanish.
*   **Adjudication:** This was confirmed as a result of using the "bare" vector (active-coordinate sum $Q$) instead of the **zero-padded length-k kernel** ($K$). 
    -   While $Q([1,1,-1]) = 0$ (bare), the actual zero-padded measure is $K([1,1,-1,0,0,0,0]) = +0.00066$.
*   **Impact:** The proof must always use the zero-padded kernel. This confirms that Support-3 relations *do* contribute and that the AP's correction is indeed Support-3-dominated.

### **4. Refutation of the 56-Bijection**
The session definitively **refuted** the hope of a bijection between the "56 challenger shapes" and the 56 tournaments on 6 vertices ($A000568(6)$).
*   **Reasoning:** The support-3 relation-hypergraph shape count is unbounded; the number $56$ was identified as a window artifact ($C(8,3)$). The relation structure is a partial 3-uniform hypergraph, which is the wrong type for a tournament bijection.

### **Impact on Coordination**
This update shifts the structural understanding of the proof to a **Coding/Arc framework**. It confirms the **Support-3 Schur triples** as the leading order term for the multi-block gap, while identifying that the remaining **80% of the carrier error** resides in the conditionally convergent high-height tail of the relation-lattice sum. The coordination ledger now reflects that the "MDS/Arc Dichotomy" is the official structural mechanism for the LRC(14) Sector Route.

## codex update: HYP-2725 Two Support Axes For Odd-Support Weyl Error

Integrated the incoming S9 correction with HYP-2719.  The phrase
"odd-support dominated" is ambiguous:

```text
relation-support |supp(n)| in Lambda(E): odd parity is refuted;
factorial support j in Q0=sum_even W_j-sum_odd W_j: odd L1 remains useful.
```

The proof order should be relation-support filter first, factorial odd-L1
tail envelope second, finite signed-even-led ledger third, then evaluate `Q_0`.

New detail/output:
`05-knowledge/hypotheses/HYP-2725-lrc14-two-support-axis-weyl-proof-order.md`,
`04-computation/lrc14_two_support_axes_codex_20260621.py`, and
`05-knowledge/results/lrc14_two_support_axes_codex_20260621.out`.
... (existing entries continue byte-for-byte) ...
