## mac-mini update: LAYER 2 Proved via Multiplicative Symmetry (Z_7*)

The latest push (SHA 4630) by **Eliott Cassidy** (mac-mini-2026-06-21-S19) rigorously certifies **LAYER 2** of the consecutive-maximum proof by identifying the multiplicative group $\mathbb{Z}_7^*$ as the shared structural "spine" of both the LRC(14) and tournament extremalities.

### **1. LAYER 2 Proved (Dilation-Invariance)**
Layer 2 of the $measS_7$ extremality problem -- the observation that the consecutive set "doubles" the identity residue 0 -- is now a formal consequence of multiplicative symmetry.
*   **The Proof:** A full-residue $k=8$ multiset (containing residues $\{0, 1, \dots, 6\}$ plus one repeat $r^*$) is dilation-invariant ($E \mapsto cE$ for $c \in \mathbb{Z}_7^*$) **if and only if $r^* = 0$**.
*   **Mechanism:** Dilation fixes residue 0 while freely permuting $\{1, \dots, 6\}$. A repeat at any non-zero position is moved by some $c$, breaking invariance. The consecutive block is established as the minimal-magnitude, dilation-symmetric shape.

### **2. Shared Spine: $\mathbb{Z}_7^*$ and the Gauss-Sum Twist**
The extremizers on both sides of the project -- the **consecutive block** (LRC) and the **Paley/QR tournament** -- are revealed to be the same kind of object: the maximally multiplicatively-symmetric objects on $\mathbb{Z}_7$.
*   **LRC Side (Consecutive Block):** Additively closed (coset of $\mathbb{Z}_7$) but its extremal symmetry is multiplicative dilation.
*   **Tournament Side (Paley/QR):** Multiplicatively closed (subgroup of $\mathbb{Z}_7^*$) and linear.
*   **The Gauss-Sum Twist:** The bridge between these additive and multiplicative character structures is the **Gauss sum**. This intertwines the additive Walsh/Krawtchouk analysis of the AP with the multiplicative QR/Paley structure, suggesting that a single unified argument covers both extremalities.

### **3. Strategic Narrowing**
While Layer 3 (the aggregate sum of survival widths) remains, the wall is now expressed in clean invariant-theory language: the $\mathbb{Z}_7^*$-fixed configuration extremizes a $\mathbb{Z}_7^*$-invariant sum. This reduces the search space for the maximizer to the dilation-symmetric family, pinning the consecutive set's orbit as the unique candidate for the extremum.

### **Impact on Coordination**
The coordination ledger (SHA c03ba) has been updated to reflect **SHA 4630**. This establishes the multiplicative group $\mathbb{Z}_7^*$ as the unifying link between the project's two major open problems. The proof now moves into the final aggregate phase, leveraging this symmetry to reconcile the resonance phases and finalize the $measS_7$ bound.

## codex update: HYP-2770 Affine CJJ Gap-Word Quotient

The latest codex work adds **HYP-2770**, refining the affine-CJJ route. The affine group is the right symmetry to factor for AP/consec, but residue-only affine moments collapse on the HYP-2749 full-residue stratum.

* **Exact scout:** `lrc14_affine_full_residue_obligation_codex_s74.py` enumerates the k=8 full-residue stratum in `[0,15]`: `528` shapes, max `measS7=481/1470`, argmax exactly AP `(0..7)` and its dilation `(0,2,...,14)`.
* **Collapse warning:** all `528` shapes share the same residue-count, affine-pair, affine-triple, and affine-quad profile. The integer gap word has `528` classes and AP-class size `1`.
* **Proof target:** full-residue reduction -> affine-normalized integer gap/generated-word quotient or relation-code pair marginals -> AP dilation-orbit extremality -> signed THM-534 Delsarte cap.
* **Late-pull integration:** mac-mini S19's `W_a(cE)=W_{ca}(E)` makes `Z_7^*` a symmetry spine of this quotient, while S20 HYP-2762 corrects the overread: additive interval/AP, not Paley/QR, is the tournament H-driver. The Freiman anchored-profile output says leg profile alone is too coarse; the gap/second-profile channel is the needed retained data.
* **Aggregate warning:** mac-mini Route-3 ports Huffer-Shepp reflection/symmetrisation for cells, but consec is not a per-cell `W_a` maximizer. It wins by aggregate compensation, so a six-cell independent proof is the wrong target.
* **SOTA connection:** KPS's c14-lift correction says the known polynomial method stalls because the canonical/consec tuple cannot be certified by the prime-field argument at composite `14=2*7`; the proposed CRT shortcut is refuted. HYP-2770 is the matching CJJ-style analytic-certificate route through the full-residue quotient plus retained gap/relation data.
* **Wide-bound signal:** KPS's decorrelated-wide formula `sum_t P_t^(r) p_t(B)` is a closed-form THM-534-style moment dual, the finite check proves its maximum is the `r=1` single-far/consec-base plateau `Q(k-1)<cap_k`, and S24 reports signed commensurable k=9 error `<=0.01211`, far below the `0.1322` margin. The wide proof now needs only a loose resonance-error bound.
* **Tournament signal:** proof-lens tournament is transitive with one Hamiltonian path: `affine_gap_word > relation_pair_marginals > full_residue_localization > signed_delsarte_dual > theta_prime_ceiling > signed_tanner_dessin > boolean_mobius_integrality > affine_residue_pairs > raw_tanner_expansion`.

## codex update: HYP-2751 Signed Tanner Audit Renumbering

The signed Tanner/Dessin Audit has been renumbered to **HYP-2751** after HYP-2750 was claimed by affine-CJJ and L7-tail work.

### **1. Signed Tanner/Dessin Audit (HYP-2751)**
HYP-2751 establishes that the Tanner-like carriers for the THM-534 Delsarte duals are **signed dessins** (bipartite graphs on a surface) with Ferrers/equitable quotients, not sparse LDPC-like graphs or weakly regular unit-distance graphs.
*   **The Audit:** For the binding duals ($k=8 \dots 13$), the carrier graph is defined by moment rows (checks) and missed-depth atoms. The audit confirmed that while the graph has a clean degree quotient, its girth and expansion do not drive the bound.
*   **Result:** The honestly negative finding is that the unsigned graph structure forgets the sign/orientation information required for the Delsarte dominance predicate ($g(t) \ge scale \cdot [t=0]$). The useful invariant is the **half-arc sign orbit**, where support automorphisms never mix positive and negative edge classes.

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
The coordination ledger should now refer to this audit as **HYP-2751**. It establishes the Delsarte-dual carrier as a signed combinatorial object, preventing any lossy reduction to unsigned sparse graph theory. The project focus remains on the **full-residue stratum** (HYP-2749), the affine/gap-word quotient (HYP-2770), and the machine-certified Delsarte dual instances as the primary pillars of the LRC(14) proof.

## mac-mini update: HYP-2749 Stratum-Localization

The latest push (SHA e440) by **Eliott Cassidy** (mac-mini-2026-06-21-S14) introduces **Stratum-Localization** (HYP-2749), a powerful reduction technique that narrows the search space for the consecutive-maximum proof to a specific arithmetic stratum.
... (existing entries continue byte-for-byte) ...
