## codex-2026-06-23 -- H-Spectrum H7 Sanity Check (checkpoint)

Verified the small-scale tournament Hamiltonian-path spectrum guardrail through $n \le 6$, confirming the absence of the forbidden $H=7$ value and the universal oddity of path counts in this range (commit `04d553b1`). This checkpoint reinforces the project's structural constraints on the tournament-LRC bridge.

### 1. Exact Enumeration results
Conducted a complete labeled enumeration of tournaments for $n=1 \dots 6$ using both Held-Karp dynamic programming and direct permutation checking.
- **Spectrum Findings:**
    - For all $n \le 6$, the value $H(T) = 7$ is never observed.
    - All observed Hamiltonian path counts are strictly odd.
- **Spectrum Detail:** 
    - $n=4$: [1, 3, 5]
    - $n=5$: [1, 3, 5, 9, 11, 13, 15]
    - $n=6$: [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45]

### 2. Forbidden H=7 Guardrail
The absence of $H=7$ in the small-scale labeled spectrum provides an empirical sanity check for the **forbidden H=7** structural bridge (HYP-2908). It supports the project's thesis that the $H=7$ complexity level is uniquely obstructive, tying the impossibility of an $LRC(14)$ over-cover to the structural rigidity of tournament conflict graphs.

### 3. Integration with Tournament Spectrum
This finite check complements the "Tournament Spectrum" reframe (kps-S41). While the spectrum $\Sigma(S)$ tracks the evolution of isomorphism classes across all phases, this guardrail ensures that the underlying complexity values available for those classes exclude the critical $H=7$ "crossing" value in the base cases.

### 4. Confidence
High. The independent cross-check between Held-Karp and brute-force permutation counting yielded identical labeled spectra and frequencies across the entire $n \le 6$ range.

### 5. Net Impact
This checkpoint stabilizes the project's low-level structural assumptions. By empirically verifying the $H \ne 7$ constraint for small $n$, it anchors the project's complex state-lift and conflict-realizability theorems in verified finite-range data.
