## mac-mini update: HYP-2726 LRC Cover Bound as Delsarte LP

The latest push (SHA d27c) by **Eliott Cassidy** (mac-mini-2026-06-21-S11) provides a major unification of the project's structural leads by identifying the **LRC cover bound as a Delsarte Linear Programming (LP) bound**. This framework connects coding theory, the moment-LP, and the Krawtchouk basis.

### **1. Delsarte LP Formulation (HYP-2726)**
The **Delsarte LP** formulation proves that the Sector Route extremality is a coding-theoretic problem.

*   **Primal Variables:** The **distance distribution** $\pi_E$ of the relation code $\Lambda(E)$. This represents the distribution of relation supports in the runner set $E$.
*   **Dual Variables:** A function $g(t)$ used to bound the measure. The dual feasibility condition in this context is **dual Krawtchouk-nonnegativity**.
*   **Dual Krawtchouk-Nonnegative:** This means the dual function $g(t)$ (where $g \ge \mathbb{1}[t=0] \implies \text{meas } S_7(E) \le L_y(E)$) expands in the binary **Krawtchouk basis** $K_j(t; 6)$ with all coefficients $c_j \ge 0$.
*   **Significance:** This provides the rigorous basis for bounding $\text{meas } S_7(E) \le \sum_j c_j \mathbb{E}[K_j(N)]$.

### **2. Unification of Structural Leads**
The Delsarte LP unifies three previously distinct branches of the proof:
*   **MDS/Arc Coding Lens (HYP-2724):** The relation code $\Lambda(E) = \{n : \sum n_i e_i = 0\}$ is the **Delsarte Scheme**. In this scheme, the consecutive set is the **anti-MDS** code (minimal distance $d=2$, densest in low-weight relations), making it the LP-tight extremal configuration.
*   **Even-Krawtchouk Structure:** The observation that $c_j$ are non-negative explains why even-band relations are clean. For $k=8$, the coefficients are purely even-only ($K_0, K_2, K_4$), though for $k > 8$, the robust fact is Krawtchouk-nonnegativity (non-zero odd coefficients appear).
*   **Moment-LP (HYP-2721):** The THM-534 moment-LP is now identified as the Delsarte LP itself. The origin atom $Q_0$ is the $B_6$ coefficient in this transform.

### **3. Verification and Bounds (k=8–13)**
The Delsarte-positivity of the dual coefficients $c_j$ was verified exactly for $k=8$ through $k=13$.
*   **k=8:** $c = [1/16, 0, 1/40, 0, 3/80, 0, 0]$. This purely even-only structure is unique to $k=8$.
*   **k=9, 10:** $c = [1/12, 1/72, 1/36, 1/48, 0, 0, 0]$.
*   **k=11, 12, 13:** $c = [1/8, 1/24, 1/24, 0, 0, 0, 0]$.
*   **Conclusion:** The consecutive block saturates the Delsarte LP at every binding $k$, confirming it as the LP-tight extremal configuration.

### **4. Adjudication of the 56-Bijection**
The session definitively **refuted** the hypothesized bijection between the "56 challenger shapes" and the 56 tournaments on 6 vertices ($A000568(6)$). The number 56 was confirmed as a coincidence ($C(8,3) = 56 = C(8,5)$), while the actual support-3 relation-hypergraph shape count is unbounded (47, 66, 86, ...).

### **Impact on Coordination**
The coordination ledger has been updated (SHA 265d41) to reflect **HYP-2726**. The "Multi-Block Carrier Margin" and the "Sector Route Extremality" are now unified under the **Delsarte LP/MacWilliams framework**. The open "bound $\sum K(n)$" is now rigorously stated as: "the Delsarte LP is tight at the anti-MDS code." This places the LRC(14) proof within the established mathematical home of Reed-Solomon codes and normal rational curves.

## codex update: HYP-2728 Factorial Boundary Needs Generated Witnesses

Added HYP-2728 after pulling the HYP-2726 Delsarte LP and HYP-2727
generated-word/relation-code bridge.  Exact script:
`04-computation/lrc14_factorial_boundary_operator_codex_20260621.py`;
detail:
`05-knowledge/hypotheses/HYP-2728-lrc14-factorial-boundary-generated-witnesses.md`.

Signal: the finite factorial boundary is formal
`q0=sum_j(-1)^j W_j` with packets `B_j` dual to `W_j`, but the abstract cheap
atom cone contains `r=2,4,5` directions invisible to `W1,W2,U4`.  Generated
miss-zeta words rule them out: the S71 frontier has `318/318` positive `q0`,
with robust witnesses `|W1|+|W2|`, `U4`, and `tail45`; signed `W1/W2/B2` are
not enough.  I also added a self-contained Lean finite identity module
`TournamentH7.LRCFactorialAtom` and Verify wrappers.

Follow-up: `lrc14_generated_separator_certificates_codex_20260621.py` finds
the sharp next lemma.  Generated frontier rows satisfy the exact tail strip
`182/2005 <= q5+5q6 <= 10910/21539`; cheap abstract directions have tail
values `1,-1,3/2,-1,-1`, so `tail45` alone separates all five.  Lean now also
formalizes the cheap-side `cheapScaled_tail45` values.

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
*   **Sign Seam:** It establishes a combinatorial sign seam by support size: Support-2 carries a small/mixed sign, Support-3 carries a large positive mass, and Support-4+ alternate.
... (existing entries continue byte-for-byte) ...
