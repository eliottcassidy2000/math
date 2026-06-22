## codex-2026-06-22-S98 -- Finite Residue Bases Are Atlases (checkpoint)

Formalized the "finite residue bases are atlases" framework and resolved the covering atlas problem by refuting global finite closures in favor of scaled residue atlases (commit `f904afc0`). This checkpoint marks the transition from a fixed-denominator shortcut to a moving-moduli spectral control.

### 1. Refutation of Universal Finite Closures (HYP-2876)
Rigorous refutation of the hypothesis that a fixed, finite set of denominators (such as `{83, 89, 21}`) can serve as a universal certificate basis for all covering LRC(14) sets.
- **The Finite-Basis Killer:** For any finite denominator list $B$, a primitive covering row $S_B = \{1, 2, \dots, 11, 13, 84 \cdot \text{lcm}(B)\}$ exists that kills every denominator in $B$ simultaneously. The trailing speed's divisibility forces that runner to the origin for every $D \in B$, zeroing the witness count $N(S,D)$.
- **Synthesis with THM-566:** This refines the previous refutation (which killed all $D \le B_0$) by showing that even sparse, hand-picked bases are vulnerable to divisor-loading. The witness denominator $D$ remains intrinsically unbounded over the class of covering sets.

### 2. Transition to Scaled Residue Atlases
The refutation pivots the proof strategy toward a **moving or scaled residue atlas**.
- **The Main Term:** For a fixed $D$, the unit witness count $N(S,D) = main\_term + resonance\_error$. While individual denominators can dip to zero due to resonance or divisibility, a healthy main term ensures witnesses exist nearby.
- **Adaptive Moduli:** Instead of a fixed basis, the target is a bounded collection of moduli chosen *after* factoring out denominators annihilated by the speed set.
- **Apex Floor strengthened:** The covering condition kills the entire small-denominator floor $\{2, \dots, 14\}$ simultaneously. This explains the sharp difficulty jump at the covering hard core, as all small rational observers are dead by definition.

### 3. Mathematical Framework: Character-Sum Resonance Counts
Established the exact mathematical language for the residue-basis approach:
- **Resonance Packets:** Resonance errors are identified as the additive relations $\sum_s k_s s \equiv 0 \pmod D$.
- **Spectral Alignment:** Coherent divisor/resonance packets are routed to the AP-core good-denominator ladder (scaled finite atlas), while incoherent packets are handled by the HYP-2875 bandlimited $L2$ tail.
- **Verification:** An S98 scout on 602 primitive covering rows confirmed that while $\{83, 89, 21\}$ is an excellent local atlas (certifying 591/602), it is distribution-sensitive and fails on adversarial divisor-loaded rows.

### 4. Structural Guardrails and reflections
- **Minor-Order Challenge:** Reconfirmed that loneliness is not minor-closed under runner deletion; any tournament-style or hole-based analogy must live on residue addresses, not speed subsets.
- **E7/Odd-Hole Analogy:** The $\{7, 21\}$ / $E_7$ diagnostic remains a useful address tool but is not a proof carrier until it preserves the witness predicate.
- **Hammon-Path Selected:** The Hamiltonian path of proof obligations is anchored by `THM566_divisor_loaded_no_go` and terminates at `speed_subset_minor_order`.

### 5. Net Impact
The "finite residue basis" idea is now correctly framed as a sampled atlas for distribution-sensitive regions rather than a universal closure. This aligns the residue route with the spectral/measure-theoretic closure, focusing rigor on the adaptive control of unblocked moduli.
