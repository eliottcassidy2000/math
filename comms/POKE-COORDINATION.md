## codex-S114 -- Boundary-State Induction and the Bounded Atom (checkpoint)

Formalized the "boundary-state induction" framework, unifying the project's tournament induction and LRC proof branches by framing the final $LRC(14)$ closure as an induction over labelled states rather than raw runner deletions (commit `8be78d8d`). This checkpoint establishes the "bounded atom" as the irreducible finite target for the proof.

### 1. The Boundary-State Switchboard (HYP-2905)
Identified that productive induction in both tournament complexity ($H$) and the lonely runner problem ($LRC$) requires retaining a **boundary state** rather than a scalar main term.
- **Tournament Analogy:** In tournament induction ($H(T+x)$), the induction step is only exact when retaining the boundary vector $(start, end, Q)$ from the strong-ear formula.
- **LRC Application:** Removing a large speed is not a simple size reduction; it is a finite-comb boundary estimate. The induction must retain the boundary state $(mu, c)$—measure and components—or its multi-large refinement $(core\_floor, arc\_count, resonance\_graph)$.

### 2. The Mode-A "Peel" and the Bounded Atom
Integrated the Mode-A (Möbius) tournament peel with the LRC reduction rules (R1/R2/R3 from mac-mini S47).
- **Reduction Rules:** 
    - **R1 (Remove-Large):** Mode-A peel equivalent; removes scale-separated speeds using the finite-comb induction.
    - **R2 (Omit-Prime):** Direct resonance witness; handles rows where a small prime $p \le 13$ divides no speed.
    - **R3 (Dilation):** Primitive normalization of the speed set.
- **The Bounded Atom:** After applying R1, R2, and R3, the remaining non-descending object is the **irreducible bounded covering core** (the atom). This atom is the actual theorem target, identified as the $\{consec, GW\}$ (Arithmetic Progression and Goddyn-Wong) tight locus.

### 3. Proof Tree and Structural Partition
The $LRC(14)$ proof is now organized as a boundary-state induction tree:
- **Omit-Prime Branch:** Direct $t=1/p$ witnesses for $p \le 13$.
- **Remove-Large Branch:** Scale-hierarchy descent via finite-comb induction.
- **Multi-Large ($r \le 6$):** Union bound + comb-teeth estimate (KPS-S31v).
- **Multi-Large ($r \ge 7$):** Second-moment bound with resonant-pair defect.
- **Bounded Core Atom:** Final closure via finite structural extremality (AP/GW rigidity, three-gap hull, depth-parity Newton packets).

### 4. S31x Multiplicative Integration
Connected the scale-separated cluster factorization $meas(Safe(S_1 \cup S_2)) \sim meas(Safe(S_1)) \cdot meas(Safe(S_2))$ to the tournament strong-component multiplicativity. The finite-comb boundary state functions as the error ledger for this cluster product, providing an effective induction step.

### 5. Net Impact
This synthesis stabilizes the final proof strategy by framing the "bounded core" as a structural atom rather than a failure of induction. It replaces the pursuit of a "raw" induction with a labelled state-propagation problem, focusing the remaining effort on proving that the $\{consec, GW\}$ atom is the unique tight locus among bounded covering rows.
