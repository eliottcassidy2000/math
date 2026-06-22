## kps-2026-06-22-S42 -- Magnitude-Aware Tournaments and Apex-Twin Separation (checkpoint)

Formalized the "Thread 2" magnitude-aware separation strategy, identifying tournaments that distinguish tight Dirichlet-extremal sets from their "apex-twin" loose counterparts (commit `3118472b`). This checkpoint recovers the metric information necessary to refine the tight-locus census beyond purely periodic winding.

### 1. Apex Winding and Magnitude-Blindness
Proven that the apex/periodic winding tournament $T(a/14)$ is inherently **magnitude-blind**.
- **The Lemma:** The winding tournament is a function of the speed multiset residue mod 14 ONLY. 
- **The Consequence:** The tight AP $\{1, \dots, 13\}$ and loose "twins" like $12 \to 26$ or $12 \to 40$ (which share the same residues mod 14) produce byte-identical apex tournaments. 
- **HYP-2925:** This confirms that the periodic order alone is insufficient to characterize the tight locus, necessitating magnitude-aware discriminators.

### 2. Magnitude-Aware Tournaments (Thread 2)
Successfully identified non-periodic tournaments that separate the tight $\{AP, GW\}$ from the apex-twin loose bank:
- **Floor-Odd Tournament:** Defined by $i \to j$ iff $\lfloor s_i / s_j \rfloor$ is odd.
- **CF-Parity Tournament:** Based on the parity of the continued fraction depth of speeds in the Stern-Brocot tree.
- **Performance:** These tournaments gave distinct isomorphism classes for apex-twins with zero false positives across a 2134-set loose bank.
- **Fingerprint:** Rows are now fingerprinted using $(score-seq, c_3, c_5, H)$, where $H$ is the number of Hamiltonian paths (calculated via bitmask-DP at $n=13$).

### 3. Discriminator Status and Invariants
Clarified that while these magnitude-aware tournaments are powerful **discriminators**, they are not complete invariants.
- **Disjointness vs. Invariance:** The tournaments provide set-disjointness on a finite bank (e.g., AP and GW receive distinct fingerprints under every candidate), rather than a single characteristic constant for all tight sets.
- **Metric Information:** This recovers the metric data that purely cyclic order-based approaches forget, providing the "metric Mac-Mini" signal requested in earlier sessions.

### 4. Honest Status Update
Reinforced the project's honest reporting regarding the "Open Core":
- **Reduction vs. Solution:** The identification of these discriminators does not close the universal $LRC(14)$ residual. 
- **The Barrier:** Characterizing the tight locus remains the irreducible open core of the Steinhaus three-gap rigidity problem. The new tournaments provide a sharper lens for the census but do not replace the fundamental rigidity proof.

### 5. Net Impact
This checkpoint provides the project with a rigorous toolkit for separating tight sets from high-resonance near-misses. By proving apex-blindness and deploying floor-odd/CF-parity tournaments, the cluster has refined the "extremal atlas" and established a formal fingerprinting protocol for candidate tight rows.
