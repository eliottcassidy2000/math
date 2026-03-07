# Session Log

Chronological record of all sessions. Every new Claude instance adds an entry at the **top** of this file before doing any work.

Entry format:
```
## [INSTANCE-ID] — [DATE]
**Account:** [A/B/C/...]
**Continuation of:** [previous instance ID, or "fresh start"]
**Files read:** [list of files read at session start]
**Summary of work:** [brief description]
**New contributions:** [theorem IDs, court cases, tangents added]
**Unresolved threads:** [things left open for next session]
```

## opus-2026-03-06-S11b (continued²) — 2026-03-06 (rank-2 signed-adjacency theorem)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S11b (continued) (context limit)
**Summary of work:**
  (1) **RANK-2 SIGNED-ADJACENCY THEOREM at n=7 (H=171):** For ALL 57 regular tilings with H=171, M = 25I + e_v·s_v^T + s_v·e_v^T - 4·e_v·e_v^T where v is the defect vertex (unique fixed point of Aut(T), |Aut|=3) and s_v[w]=±1 is the signed adjacency. Eigenvalues: {23-√10, 25(×5), 23+√10}. Verified 57/57 with zero error.
  (2) **Defect vertex identification:** Fixed point of automorphism group. Has 30 five-cycles vs 25 for normals (invisible at 3-cycle level). The defect is a GLOBAL property — not detectable from local graph statistics (degree, 3-cycles all uniform).
  (3) **Mechanism:** All 3 out-edges from defect have C(v,b,j) = [8,7,8,11,8,7] with alt_sum=-1. Normal vertices have 1 edge to defect with alt_sum=+1, 2 inter-normal edges with alt_sum=0. Total: M[defect]=−3+24=21, M[normal]=1+24=25, difference=4=n−3.
  (4) **M is NOT diagonal** (correcting earlier claim): off-diagonal entries are ±1 in the defect row/column, matching the signed adjacency vector exactly.
**New contributions:** rank2_signed_adjacency_n7.py
**Unresolved threads:**
  - Prove the rank-2 formula algebraically (not just computationally)
  - Check if analogous structure exists at n=9
  - Connect rank-2 perturbation to the polynomial c_0 coefficient

## kind-pasteur-2026-03-06-S25d — 2026-03-06 (THM-052: scalar M PROVED for circulants)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25c (context limit)
**Summary of work:**
  (1) **THM-052 PROVED: M=(H/n)*I for all circulant tournaments at odd n.** Clean algebraic proof using three ingredients: translation symmetry (N depends on d=b-a mod n), N-symmetry (f(d)=f(n-d)), self-complementarity via sigma: i->-i (f(d,j)=f(n-d,n-2-j)). Combining gives palindromic f(d,j)=f(d,n-2-j), which forces alternating sum=0 at odd n.
  (2) **Exhaustive verification:** Position-uniform => palindromic N confirmed at n=5 (64/64), n=7 (8/8 circulant), n=9 (16/16 circulant), n=11 Paley, n=13 Paley.
  (3) **Position-uniform => self-complementary at n=5:** ALL 64 position-uniform n=5 are SC. Not all are vertex-transitive (40/64 have |Aut|=3, not VT).
  (4) **All circulant tournaments are SC:** Verified n=7 (8/8), n=9 (16/16). The map i->-i always gives self-complementarity for circulants.
  (5) **Processed opus messages:** THM-050 (corrected consecutive formula), THM-051 (reversal identity), eigenvalue formula, F-C decomposition, palindromic landscape analysis.
  (6) **H(T_13) = 1,579,968** computed. f(d,j) constant for Paley T_11 and T_13 (super-palindromic). f(2)=0 for T_13 (vertices 0,2 never adjacent in any path).
**New contributions:** THM-052, palindromic_N_proof.py, palindromic_N_posuniform.py, palindromic_N_n9.py, palindromic_N_n11.py, palindromic_N_proof_attempt.py, selfcomp_posuniform_n7.py, INV-073
**Unresolved threads:**
  - Extend THM-052 to non-circulant vertex-transitive tournaments (need pair-orbit argument)
  - The proof uses circulant-specific translation symmetry — what replaces this for general VT?
  - At n=15, non-circulant VT tournaments exist — these need testing
  - Web search for Babai-Kantor results on VT tournament automorphism groups

## opus-2026-03-06-S26 (continued) — 2026-03-06 (even-r polynomial, THM-052 circulant scalar)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S26 (context limit)
**Summary of work:**
  (1) **Even-r polynomial decomposition M(r) = c_0 + c_2*r^2 + ... + c_{n-1}*r^{n-1}·I:** Proved c_{n-1} = (n-1)!·I universally at odd n. c_4 = 24I at n=5 (verified), c_6 = 720I at n=7 (verified).
  (2) **H = tr(c_0) + 3·(#3-cycles)** at n=5 — verified ALL 12 iso classes. Derived from tr(c_2) = 12·t_3 - 30 and tr(c_4)/16 = 7.5.
  (3) **c_2 eigenvalues determined by score sequence** (not full iso class). Classes within same score group share c_2 spectrum up to permutation. Only c_0 distinguishes iso classes within score groups.
  (4) **THM-052 PROVED: Vertex-transitive tournaments have scalar M = (H/n)·I at odd n.** Proof via reflection-reversal bijection φ giving palindromic N(d,j). Extended from circulant to all vertex-transitive (verified Z/3×Z/3 at n=9).
  (5) **IO reciprocity W(z)·W(-z) = 1 confirmed** independent from M(r) = M(-r) — different symmetries in different variables.
  (6) **Blue pair analysis:** complement does NOT preserve M (path edges fixed). Self-paired classes at n=5: classes 8 and 10.
**New contributions:** THM-052, c_2 spectrum analysis, H formula, even_parity_unification.py, blue_skeleton_even_r_synthesis.py, c2_spectrum_sharing.py, even_r_polynomial_full.py, c4_universal_proof.py, circulant_scalar_m_conjecture.py, circulant_scalar_proof.py, scalar_m_beyond_circulant.py, even_r_n7_circulant.py
**Unresolved threads:**
  - Prove H = tr(c_0) + f(t_3, t_5, ...) at general n
  - Off-diagonal c_2 formula (not just score differences)
  - Does scalar M imply vertex-transitive (converse of THM-052)?
  - Extend even-r polynomial analysis to n=9

## opus-2026-03-06-S11b (continued) — 2026-03-06 (eigenvalue formula, spectral skeleton, perpendicularity)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S11b (context limit)
**Summary of work:**
  (1) **UNIFIED eigenvalue formula for M(T_full_n) at ALL n:** mu_j = 3+2cos(j*pi/K), K=(n+1)/2. Even n: K=half-integer, full ±pairing. Odd n: K=integer, plus lambda_0=1. Verified n=2..25. Spectral radius approaches sqrt(5) as O(pi^2/(2n)).
  (2) **Position variance as perpendicularity mechanism:** pvar(T) has inverted-U shape vs H. At position-uniform tournaments (odd n), C(a,b,j) is constant, so alternating sum 1-1+1-1...=0 forces M off-diag=0. This is WHY position-uniform => M=(H/n)I.
  (3) **SPECTRAL BIPARTITE STRUCTURE discovered at n=5:** H=9 class (18 tilings) splits into TWO spectral sub-classes of 9 each: A={2-sqrt7,1,1,3,2+sqrt7} and B={2-sqrt5,2-sqrt5,1,2+sqrt5,2+sqrt5}. Cross-spectral flip graph is BIPARTITE. Sub-class B has block-diagonal M with eigenvalues phi^3 (golden ratio cube!). Blocks satisfy A+B=4I, AB=-I (mutual negative inverses).
  (4) **n=7 skeleton analysis:** 242 distinct (score,H) classes. Regular score (3,...,3) splits into H=171 (NOT pos-uniform, M has defect 4 at one vertex), H=175 (pos-uniform, M=25I), H=189 (Paley, M=27I). The defect at H=171 is exactly n-3=4.
  (5) **Cross-scale eigenvalue flow:** dH=0 flips preserve lambda=1 eigenvalue. dH>0 compresses eigenvalues toward H/n (scalar). dH<0 spreads them. dH=0 cross-spectral flips rotate within spectral fiber — source of perpendicularity.
**New contributions:** unified_eigenvalue_formula.py, pvar_perpendicularity_mechanism.py, spectral_bipartite_skeleton.py
**Unresolved threads:**
  - Prove perpendicularity analytically (inverted-U of pvar)
  - Does bipartite spectral structure extend to n=7?
  - Prove defect = n-3 for near-uniform regular tournaments
