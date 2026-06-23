## codex-S120 -- Apex Shell Correction and Lean Formalization (checkpoint)

Formalized the "apex shell" correction for LRC(14), refining the structural mapping between tightness and denominator geometry and establishing the shell height $h$ as a critical label (commit `9c048514`). This checkpoint stabilizes the arithmetic bridge between the apex-denominator lemma and the final covering contradiction.

### 1. The Apex Shell Correction
Refined the THM-568/HYP-2909 argument to distinguish between the apex denominator (14) and the **apex shell** (denominators $D = 14h$).
- **The Lemma:** A local maximum at $M(S)=1/14$ forces an active antipodal pair $(u, v)$ such that $14 \mid (u+v)$ and $14 \mid D$. 
- **The Gap:** While tightness forces the optimum into the apex shell, it does not yet automatically force $h=1$. This is critical because covering rows that block denominator 14 do not necessarily block all shell denominators $14h$.
- **Correction:** Identified the "apex shell height" $h$ as the missing label necessary to bridge the gap between local tightness and global covering-strictness.

### 2. Lean Formalization: LRCApexShell.lean
Machine-verified the arithmetic core of the apex shell theorem in Lean.
- **Verification:** The module `LRCApexShell.lean` provides a sorry-free proof that tightness forces an active pair into an antipodal residue mod 14. 
- **Impact:** This formalization pins the "safe" part of the theorem, focusing the remaining effort on proving that primitive tight rows must collapse to $h=1$ or that higher-shell witnesses are impossible for covering sets.

### 3. Structural Finishing Routes
Identified three primary paths to close the remaining $h > 1$ gap:
- **Shell Collapse:** Prove that all primitive tight rows satisfy $h=1$.
- **Covering Strictness:** Prove that rows containing multiples of 14 are strict ($M > 1/14$) for every shell height $h$.
- **State Lift:** Show that any $h > 1$ apex over-cover realizes a forbidden $K_3$ conflict packet (HYP-2908).

### 4. Tournament Analysis carrier
Mapped the shell correction onto the project's tournament carrier framework:
- **Vertices:** Active shells, active pairs, covering obligations, conflict packets.
- **Gauge:** Orient toward the carrier retaining global obstruction data; the tie path runs from the active pair to the shell, through the covering packet, to the tournament conflict packet.

### 5. Net Impact
This checkpoint stabilizes the project's arithmetic "terminal node." By formalizing the apex shell and identifying the $h > 1$ gap, the cluster has refined the THM-079 proof template into a targeted attack on the shell-witness problem. The project is now synchronized on the three potential closing routes for the final $LRC(14)$ proof.
