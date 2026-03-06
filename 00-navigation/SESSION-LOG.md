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

## kind-pasteur-2026-03-06-S20 — 2026-03-06 (Computational verification: DRT/LTT classification, DC-OCF test)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S19
**Summary of work:**
  (1) **DRT/LTT/Other classification at n=7 — VERIFIED:** The 2640 regular n=7 tournaments split into exactly 3 classes: DRT (Paley, 240), Locally Transitive (720), Other Regular (1680). Directed cycle counts are CLASS INVARIANTS: DRT={3:14,5:42,7:24}, LTT={3:14,5:28,7:17}, Other={3:14,5:36,7:15}. DRT maximizes at EVERY odd length. Confirms Savchenko's "diametrically opposite" characterization.
  (2) **Deletion-contraction does NOT preserve OCF — NEGATIVE RESULT:** At n=4, OCF fails for T\e (39.3%) and T/e (60.7%). OCF is tournament-specific. This blocks naive DC induction but the noncommuting Rédei-Berge framework operates at a richer level.
**New contributions:** T128-T129, INV-051 updated (tested), INV-053 updated (verified), savchenko_cycle_test.py, deletion_contraction_ocf.py
**Unresolved threads:**
- INV-052 (chromatic-Rédei-Berge bridge) still highest priority for next deep read
- Obtain Savchenko's actual polynomial formulas for c_k
- Verify DRT cycle invariance at n=11

## kind-pasteur-2026-03-06-S19 — 2026-03-06 (EXTENSIVE web research: 40+ queries, 16 new investigation leads)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S18h
**Summary of work:**
  Conducted the most extensive web search campaign of the project: 6 parallel search agents + 20+ direct searches covering tournament polynomials, Hopf algebras, extremal theory, graph polynomials, Rédei extensions, OEIS/computational leads.

  **Major new leads discovered:**
  (1) **Mitrovic noncommuting Rédei-Berge function** (arXiv:2504.20968, Apr 2025) — has DELETION-CONTRACTION, enabling inductive proofs. The commutative version lacks this. HIGH PRIORITY for OCF induction. (INV-051)
  (2) **Mitrovic-Stojadinovic chromatic↔Rédei-Berge bridge** (arXiv:2506.08841, Jun 2025) — proves these functions are "almost identical" at poset level. Enables importing chromatic symmetric function results. Proves "converse of Rédei's theorem." (INV-052)
  (3) **Savchenko cycle counting series** (2016-2024) — exact polynomial c_k formulas for regular tournaments. c8(DRT) independent of which DRT. Phase transition at n=39 where DRT advantage reverses. (INV-053)
  (4) **Komarov-Mackey 5-cycle formula** (JGT 2017) — exact c5 from edge score sequence. (INV-054)
  (5) **Linial-Morgenstern cycle minimization** — spectral methods for cycle density extremals. (INV-055)
  (6) **Jerrum-Patel zero-free regions** (JLMS 2026) — characterizes which H-free classes have real-rooted I.P. (INV-056)
  (7) **Herman Terwilliger algebras** (2024) — classifies DRTs, 237 non-isomorphic at n=27. (INV-057)
  (8) **Pantangi critical groups** distinguish Paley from non-Paley DRTs. (INV-058)
  (9) Additional leads: cyclic subtournament Hamiltonicity (INV-059), Eulerian trace formula (INV-060), Hamilton transversals (INV-061), forward arc maximization (INV-062), Paley spectral pseudorandomness (INV-063), Mitrovic Hopf new bases (INV-064), IP root gap (INV-065), low-rank tournament matrices (INV-066).

**New contributions:** INV-051 through INV-066 (16 new leads), T121-T127 (7 new tangents)
**Unresolved threads:**
- INV-051 (noncommuting Rédei-Berge) and INV-052 (chromatic bridge) are the HIGHEST priority for next session
- INV-053 (Savchenko formulas) could immediately verify/extend our cycle maximization theory
- 6 background search agents may have additional results not yet integrated

## opus-2026-03-06-S7 — 2026-03-06 (tiling isomorphism classes: SC+SF kernel discovery, n=7 analysis)
**Account:** Eliott (opus machine)
**Continuation of:** opus-2026-03-06-S6
**Summary of work:**
  Deep investigation of tournament isomorphism classes and their tiling structure, focused on how symmetrical tilings group and how symmetry persists across n values.

  (1) **SC+SF SYMMETRY KERNEL DISCOVERED.** Classes that are both self-converse (T ~ T^op) AND self-flip (flipping all non-path arcs gives isomorphic tournament) form a persistent "kernel." Count: n=4:1, n=5:2, n=6:2, n=7:8. These kernel classes have nearly-regular scores and H values in the upper range.

  (2) **EVEN/ODD PATTERN for H-maximizer.** At even n (4,6), the H-maximizer IS in the SC+SF kernel. At odd n (5,7), the H-maximizer (regular/Paley, high |Aut|) is SC but NOT SF. Proved computationally: Paley(7) fails self-flip for ALL 189 Hamiltonian paths — the flip always changes the score sequence from regular to non-regular.

  (3) **GS ~ 3^floor((n-2)/2) in kernel.** Grid-symmetric tilings per kernel class follow this formula exactly at n=4,5,6 and for most classes at n=7 (6/8 have GS=9=3^2). The exceptions have slightly different score types.

  (4) **n=7 FULLY ENUMERATED.** 456 isomorphism classes (matching OEIS A000568). 88 self-converse, 30 self-flip, 8 both. 3 regular tournament classes: Paley (H=189, |Aut|=21), H=175 (|Aut|=7), H=171 (|Aut|=3). All regular tournaments are SC with anti-auts but none are SF.

  (5) **Self-flip mechanism explained.** Self-flip requires a permutation that is anti-automorphism on non-path arcs but automorphism on path arcs. This MIXED condition is impossible for vertex-transitive tournaments because the base path cannot be distinguished from other arcs by any automorphism. Low-|Aut| tournaments have enough asymmetry to satisfy this condition.

  (6) **SF-only classes come in transpose pairs.** At n=7, all 22 SF-only classes form 11 pairs {T, T^op} with equal H.

**New contributions:** T128, T129, tiling-isomorphism-kernel.md
**Unresolved threads:**
- n=8 kernel count (requires 2^28 = 268M tilings — needs C implementation)
- Prove GS = 3^floor((n-2)/2) formula for kernel classes
- Prove even/odd pattern for H-max in kernel
- Explain why kernel count jumps from 2 to 8 between n=6 and n=7

## kind-pasteur-2026-03-06-S18h — 2026-03-06 (BIBD cycle maximization theorem, directed cycle analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S18g
**Summary of work:**
  (1) **THM-028 PROVED: BIBD Arrangement Maximizes Directed Cycles, Not Disjoint Pairs.** Among all 2640 regular n=7 tournaments, the BIBD (Paley, uniform lambda=2) MINIMIZES alpha_2=7 (disjoint 3-cycle pairs) but MAXIMIZES total directed odd cycles (alpha_1=80). Three rigid classes: (alpha_2=7, H=189, 240 tours), (alpha_2=10, H=171, 1680), (alpha_2=14, H=175, 720). H = 1 + 2*alpha_1 + 4*alpha_2 verified for all 2640.
  (2) **BIBD alpha_2 formula proved:** D = C(b,2) - p*C(r,2) + sum C(lambda_e, 2). The BIBD minimizes the convex sum by Jensen's inequality. Verified at p=3,7,11.
  (3) **Key mechanism identified:** BIBD forces every 5-vertex subtournament to be the regular T_5, which has exactly 2 directed Hamiltonian cycles. Non-BIBD arrangements create less regular subtournaments with fewer directed 5-cycles (28-36 vs 42).
  (4) **CORRECTS previous hypothesis:** Opus's T102 suggested alpha_2 (disjoint pairs) drives H-maximization. WRONG. Alpha_1 (total directed cycles) dominates via the linear term 2*alpha_1. The BIBD's "evenness" creates MORE directed cycles at the cost of FEWER disjoint pairs.
  (5) **Read and analyzed opus's discoveries:** THM-025 (real-rootedness disproved at n=9), transfer matrix symmetry (symbolic proof at n=4-7), Fano-Paley BIBD, determinantal IP approach, Omega_3 complement = matching.
**New contributions:** THM-028, T120, INV-042 updated
**Unresolved threads:**
- Verify BIBD cycle maximization at p=11 (computationally expensive)
- Prove BIBD forces subtournament regularity (algebraic argument needed)
- Prove alpha_1 maximization from BIBD structure at general p
- Transfer matrix symmetry: find conceptual proof (INV-001)

## opus-2026-03-06-S6 — 2026-03-06 (transfer matrix deep analysis: trace formula proved, [[1,0],[0,-1]] disproved)
**Account:** Eliott (opus machine)
**Continuation of:** opus-2026-03-06-S5
**Summary of work:**
  Investigated transfer matrix structure in depth. Three major findings:

  (1) **MISTAKE-011: M = [[1,0],[0,-1]] claim is FALSE.** The 2×2 transfer matrix is NOT always diag(1,-1). Exhaustive check at n=4 shows 2199/2500 failures. M values vary widely (observed -3 to +3 at n=5). The only universal property is SYMMETRY M[a,b] = M[b,a].

  (2) **THM-027 PROVED: Transfer Matrix Trace Formula.** The n×n transfer matrix satisfies tr(M) = H(T) for odd n, 0 for even n. PROOF: Each pair (σ through S ending at a, τ through R starting from a) bijects with a Hamiltonian path P. Contribution is (-1)^{pos(a,P)}. Summing over all a gives (1-(-1)^n)/2 per path. Verified exhaustively n=3,...,7.

  (3) **Off-diagonal sum formula:** sum_{a≠b} M[a,b] = 0 for odd n, 2*H(T) for even n. Verified n=3,...,7.

  Also analyzed the per-subset cancellation structure at n=4,5 looking for an involution proof of symmetry. The Cauchy-Binet decomposition M = E^T·Λ·B shows no simple pairing. The complement pairing D(S)+D(U\S) is constant at n=4 but not at n=5.

**New contributions:** THM-027, MISTAKE-011, T119
**Unresolved threads:**
- Prove transfer matrix symmetry M[a,b] = M[b,a] for general n (INV-001/INV-045)
- Prove off-diagonal sum formula for general n
- The trace formula proof is complete; can it inspire a symmetry proof?

## opus-2026-03-06-S5 — 2026-03-06 (deep web synthesis — Hopf algebra, Feng reversibility, DRT, Irving-Omar)
**Account:** Eliott (opus machine)
**Continuation of:** opus-2026-03-06-S4
**Summary of work:**
  Deep web investigation across all leads with 4 parallel research agents. Major synthesis of 12+ papers and connections.

  KEY FINDINGS:
  (1) **Irving-Omar (2412.10572):** Corollary 20 IS our OCF. ham(D) = Σ_S det(Ā[S])·per(A[S^c]). Walk generating function W=det(I+zXĀ)/det(I-zXA).
  (2) **Hopf algebra route to INV-001 (NEW, INV-045):** Grujić-Stojadinović Hopf comultiplication = our subset convolution. Combined with Feng's dual Burnside reversibility theorem, the tournament constraint T[x,y]+T[y,x]=1 should play role of detailed balance → transfer matrix symmetry.
  (3) **DRT theory (T116, INV-047):** Paley tournaments are doubly regular tournaments ↔ skew Hadamard matrices (Reid-Brown 1972). Nozaki-Suda characterize skew Hadamard via spectra of size n-2 tournaments. Connects our spectral regularity finding to maximizer theory.
  (4) **Asymptotic convergence (T117, INV-048):** Adler-Alon-Ross proved max H(T) ≥ (e-o(1))·n!/2^{n-1}. Our Paley ratios converge toward e. Paley = quasi-random explains this.
  (5) **El Sahili-Ghazo Hanna (2023):** T and T^op have same oriented Ham path TYPE distribution. Our transfer matrix identity M_{T^op}=(-1)^{n-2}M_T is stronger.
  (6) **Our BIBD discovery appears NOVEL:** No prior work found connecting Paley tournament 3-cycles to BIBDs or Fano plane decomposition.
  (7) **Pantangi (2019):** Critical groups distinguish Paley from other DRTs — potential algebraic invariant for H-maximization.
  (8) **Satake (2025):** New cyclotomic nearly-doubly-regular tournaments. Savchenko's conjecture.

**New contributions:** T114-T118, INV-045 through INV-050, web-synthesis-opus-S5.md
**Unresolved threads:**
- Formalize Feng reversibility → transfer matrix symmetry proof (INV-045, highest priority)
- Check if non-Paley DRTs exist at small p and compare H values (INV-047)
- Compute H(T_p)/(p!/2^{p-1}) at p=31 (INV-048)
- Read Ai (2025) on new digraph polynomials (INV-049)

## kind-pasteur-2026-03-06-S18g — 2026-03-06 (hereditary CORRECTION, deletion ratio formula, R-minimization refuted)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S18f
**Summary of work:**
  (1) **HEREDITARY CORRECTION (MISTAKE-010)**: Previous claim "all maximizers at odd n are hereditary" was WRONG. Exhaustive check: n=5 has 64 maximizers, only 24 (regular, score (2,2,2,2,2)) are hereditary. The 40 non-regular (score (1,2,2,2,3)) have del_hs=[3,5,5,5,3], NOT all max H(4)=5. Correct statement: only REGULAR maximizers at odd n are hereditary.
  (2) **THM-026 PROVED: Deletion-Sum Ratio Formula.** R(T) = sum_v H(T-v)/H(T) = n - E_weighted[|U(S)|] where E is over independent sets of Omega(T) weighted by 2^{|S|}, and |U(S)| = total vertices covered. Equivalent: U_sum = 2*sum_C |V(C)|*mu(C). Verified exhaustively n=3-6. (Renumbered from THM-025 to avoid collision with opus's THM-025 disproof.)
  (3) **R-MINIMIZATION CONJECTURE REFUTED at n=7.** Conjectured H-maximizer minimizes R. Holds at n=3-6 but FAILS at n=7: tournaments with H=123 have R=1.585 < R(max)=5/3. The formula R = n - E[|U|] is proved but is NOT a variational characterization of maximizers.
  (4) **Score rigidity at n=4**: ALL 24 tournaments with score (1,1,2,2) have H=5=max. This means hereditary n=5->n=4 is trivially forced — any regular deletion lands on a maximizer by score rigidity.
  (5) **Non-trivial hereditary at n=7->n=6**: Score (2,2,2,3,3,3) has 2640 tournaments with H in {41,43,45}; only 480 achieve max. Paley T_7 deletion specifically lands on Type B (IP=[1,20,1]).
  (6) **U_sum decomposition**: For Type A (n=6): 2*54 + 4*24 = 204. For Type B: 2*84 + 4*6 = 192. Type A wins via more disjoint pairs despite fewer total cycles.
**New contributions:** THM-026, MISTAKE-010, T113, OPEN-Q-017 (refuted), OPEN-Q-018
**Unresolved threads:**
- Why does Paley T_7 deletion land on a maximizer? (score rigidity fails at n=6, so structural argument needed)
- Prove regular hereditary maximizer for general odd n
- Find weaker variational principle that DOES characterize maximizers

## opus-2026-03-06-S18 — 2026-03-06 (DISPROVE real-rootedness of I(Omega(T), x) at n=9)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S17 (multi-session overnight)
**Summary of work:**
  (1) **DISPROOF OF REAL-ROOTEDNESS AT n=9 (THM-025)**: Found explicit counterexample tournament with score sequence [1,1,3,4,4,4,6,6,7]. Full Omega has 94 directed odd cycles, I(Omega, x) = 1 + 94x + 10x^2 + x^3. Newton's 2nd inequality FAILS: a_2^2=100 < a_1*a_3*(3/2)=141. Two complex roots. H(T)=I(Omega,2)=237 (OCF unaffected). This disproves OPEN-Q-015 conjecture for general n.
  (2) **Omega_3 also fails**: For same tournament, I(Omega_3, x) = 1 + 12x + 6x^2 + x^3 with discriminant = -1323 < 0. Newton's inequalities both HOLD for Omega_3 (a subtlety: disc < 0 is strictly stronger than Newton failure for degree-3).
  (3) **n <= 8 remains valid**: THM-020 (claw-free => real roots via Chudnovsky-Seymour) confirmed. 200 random n=8 tournaments with full Omega: 0 failures. The claw-free bound is SHARP.
  (4) **Failure is extremely rare**: 0 failures in 100k random n=9 tournaments (Omega_3) and 0 in 1323 high-variance n=9 tournaments (full Omega). The counterexample requires very specific structure (11/12 three-cycles share vertex 3). Arc-flip analysis: 3/36 neighbors also fail, all involving vertex 3.
  (5) **Cut-flip analysis, determinantal IP**: Explored cut-flip delta preservation (trivial since OCF proved), Omega_3 complement = matching at n<=6 (new, T109), Turan proves Newton's 1st inequality for n<=11 (T110).
**New contributions:** THM-025 (disproof), T109-T112, updated OPEN-Q-015 (resolved/disproved), updated THM-020
**Scripts:** full_omega_n9_correct.py, real_roots_n9_verify.py, newton_n8_sample.py, real_roots_failure_rate.py, and others
**Unresolved threads:**
- What fraction of n=9 tournaments fail real-rootedness? (appears < 1 in 100k)
- Does failure occur at n=10 more frequently?
- New question: what structural property of Omega(T) ensures real-rootedness for "most" tournaments?

## kind-pasteur-2026-03-06-S18f — 2026-03-06 (hereditary maximizer chain + spectral analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S18e
**Summary of work:**
  (1) **HEREDITARY MAXIMIZER CHAIN (exhaustive n=3-7)**: At odd n=3,5,7: EVERY vertex deletion from the global H-maximizer gives the (n-1)-maximizer. At even n=4,6: NO deletion gives the (n-1)-maximizer. The deletion spectrum is CONSTANT for all maximizers (vertex-transitive). Two types at n=6: Type A (del=11, IP=[1,14,4], more disjoint pairs) and Type B (del=13, IP=[1,20,1], more total cycles). Both perfectly vertex-homogeneous (4 c3 + 5 or 10 c5 through every vertex).
  (2) **SCORE OBSTRUCTION ANALYSIS**: The even-n hereditary failure is NOT due to score incompatibility — deletion score at n=6 IS a maximizer score at n=5. The subtournament is simply suboptimal within that score class. At odd n, regular score + Paley structure forces deletion to be optimal. Score change: regular n → delete v → SC score (d-1,d-1,...,d,d) = maximizer score at n-1.
  (3) **AA^T SPECTRAL CORRELATION**: corr(H, lambda_1(AA^T)) = -0.97 at n=5 and -0.96 at n=6. H-maximizers have smallest leading eigenvalue of AA^T (most spectrally regular). Spectral gap also strongly anti-correlated with H.
  (4) **Cayley transform test**: tr(S^k) for skew-symmetric S is always 0 (trace of skew-symmetric is 0). The Irving-Omar arctanh approach uses adjacency A not skew S. tr(A^3)/3 = c3 confirmed. Pfaffian not directly related to H. The determinantal representation exists iff roots are real negative (circular), so the key is finding a NATURAL matrix from tournament structure.
  (5) **n→n-2 pair deletions**: From T_7, ALL 21 pair-deletions give score (1,2,2,2,3) with H=13 (not max 15). The n→n-2 chain does not preserve maximality.
  (6) **SC MAXIMIZER AT n=8 CONFIRMED**: SC tournaments with score (3,3,3,3,4,4,4,4) achieve H=661 = global max (OEIS A038375). Generated via fpf involution (2^16 per sigma). 19 SC score classes, all tested.
  (7) **HEREDITARY BREAKS AT n=8** (as predicted): No deletion from n=8 maximizer gives max H(7)=189. Deletion H-values: 131, 133, 151 (not constant, unlike n=6). Even-odd dichotomy confirmed n=3-8.
**New contributions:** T104 (hereditary maximizer), T105 (AA^T spectral), INV-044, OPEN-Q-016 extended to n=8
**Unresolved threads:**
- Prove hereditary maximizer algebraically (score compatibility proved, optimality within class needs proof)
- Find natural PSD matrix M_T such that det(I+xM_T) = I(Omega(T), x)
- Prove SC maximizer within each score class (not just global max) at n=8

## opus-2026-03-06-S18 — 2026-03-06 (cross-refs complete, THM-024 proved, Paley deletion verified)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S17
**Files read:** All theorem files, TANGENTS, INVESTIGATION-BACKLOG, OPEN-QUESTIONS, SESSION-LOG
**Summary of work:**
  (1) **CROSS-REFERENCING COMPLETE**: Added "Verification Scripts" sections to all 15 remaining theorem files (THM-003, 005-010, 012-014, 021, 023, PROP-001, LEM-001, LEM-002). All 30 theorem files now have script references.
  (2) **THM-024 PROVED: Every SC tournament has an involution anti-aut.** Corrects T095/T093 which claimed ALL anti-auts are involutions (false at n=6). Proof uses Moon's theorem (|Aut(T)| odd) + Cauchy's theorem (⟨Aut(T), σ⟩ has even order → order-2 element in anti-aut coset). Clean 5-line group theory proof.
  (3) **PALEY DELETION VERIFIED at p=11**: H(T_11 − v) = 15745 = a(10) = max H at n=10 (OEIS A038375). Extends the Paley deletion pattern: T_3→a(2), T_7→a(6), T_11→a(10). All confirmed.
  (4) **QUASI-REGULARITY EXPLAINED (T103)**: Omega_3(T) is quasi-regular because adjacency depends on vertex-set intersection (Johnson graph structure). CV of degree = O(1/√m) → 0, forcing λ_max/avg_deg → 1. Verified n=5-20.
  (5) **WINDOWS PATH FIXES**: Fixed 7 more scripts with Windows paths (kind-pasteur S18e additions).
  (6) **CLAIM A DECOMPOSITION for T_7**: sum_mu = 6×3 + 30×1 + 24×1 = 72. All 3-cycle complements in Paley T_7 have a 3-cycle (mu=3). Corrected c_5 count (42 directed, not 21).
**New contributions:** THM-024, T102, T103, INV-042, INV-043
**Unresolved threads:**
- SC maximizer at n=8 (computationally expensive, OPEN-Q-016)
- Paley deletion at p=19 (need H(T_19−v) = a(18)?)
- Full involution existence test at n=7 (too slow for brute force; proved theoretically via THM-024)
## opus-2026-03-06-S4 — 2026-03-06 (Fano-Paley, transfer matrix polynomial proof, Paley deletion p=19)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S3 (context continuation, ran overnight)
**Files read:** Full warm-up sequence, kind-pasteur S18e messages, sc-maximizer-mechanism.md, interlacing-clique-deletion.md
**Summary of work:**
  (1) **PROVED: Anti-aut sigma induces Omega(T) automorphism.** Clean proof: sigma maps directed cycle C to reversed cycle on sigma(V(C)), preserving vertex-sharing. At even n, sigma* is fixed-point-free on cycles.
  (2) **PROVED: Path reversal identity** M_{T^op}[i,j] = (-1)^{n-2} M_T[j,i] via subset complement.
  (3) **BREAKTHROUGH: Transfer matrix symmetry is a POLYNOMIAL IDENTITY** (T103). M[a,b]-M[b,a]=0 after tournament substitution T[j,i]=1-T[i,j]. Proved symbolically at n=4,5,6,7. With independent variables, nonzero (12-48 terms). Tournament constraint is essential. This is the key to INV-001.
  (4) **DISCOVERED: Fano-Paley connection** (T102). The 14 cyclic triples of T_7 = TWO copies of the Fano plane PG(2,2). Each of 7 disjoint pairs takes one line from each Fano. alpha_2=7 perfectly explained. Generalizes: cyclic triples of T_p form 2-(p,3,(p+1)/4) BIBD.
  (5) **PALEY DELETION p=19**: H(T_19-v)=117,266,659,317 for all v. Self-comp scores (8^9,9^9). Conjecture a(18)=117266659317.
  (6) **T_11 design**: 55 cyclic triples form a 2-(11,3,3) BIBD. 495 disjoint pairs, 550 disjoint triples. Aut-transitive on triples.
  (7) **Omega clique analysis**: At n=5, Omega is always a clique (theta=1, alpha_2=0). At n=6, H=45 achievers have two routes: alpha=[1,14,4] or [1,20,1].
**New contributions:** T102 (Fano-Paley), T103 (transfer polynomial), T104 (Paley p=19), symbolic_symmetry_proof.py
**Unresolved threads:**
- Prove transfer symmetry conceptually (not just by symbolic computation)
- SC maximizer proof: Fano structure might give algebraic handle
- n=8 SC maximizer test
- Transfer matrix symmetry: M_{T^op} = (-1)^{n-2} M_T^T is interesting but doesn't directly give symmetry

## opus-2026-03-06-S17 — 2026-03-06 (comprehensive verification + spectral analysis of Omega)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S18
**Files read:** Full warm-up sequence, all theorem files, all computation scripts
**Summary of work:**
  (1) **COMPREHENSIVE THEOREM AUDIT**: Verified all 30 theorem files. Found and resolved: duplicate THM-013 (renamed to THM-012b), duplicate THM-016 (removed hamiltonian-alternating-sum, kept claim-b), wrong MISTAKE-006 cross-ref in THM-012b, outdated claw-free question in THM-019, LEM-001/LEM-002 certainty mismatch (2->5), THM-001 Route D outdated conditionality.
  (2) **31 PATH FIXES**: Fixed 6 Windows paths and 25 non-portable relative paths across computation scripts.
  (3) **PROVED: No blueself at ANY odd n** -- algebraic proof via endpoint score multiset analysis. THM-022 Theorem 5 upgraded from computational (n<=7) to universal.
  (4) **DISPROVED: Blueself = SC maximizer** at n=6. Blueself H=41 is not SC max (SC max=45).
  (5) **DISCOVERED: Through-v cycles form clique** in Omega(T) (proved). At n=5, 100% of remaining cycles adjacent to some through-v cycle.
  (6) **SPECTRAL ANALYSIS OF OMEGA_3(T)**: lambda_max/avg_degree ~ 1.005-1.011 for all n=5-15. Omega is quasi-regular. 0 real-root failures in 700 samples (n=5-10).
  (7) **VERIFICATION SCRIPTS**: Created verify_all_theorems.py, blueself_odd_n_proof.py, interlacing_verify.py, interlacing_structure.py, blueself_sc_maximizer_connection.py, omega_spectral_fast.py.
  (8) **CROSS-REFERENCING**: Added verification script paths to 5 theorem files.
**New contributions:** THM-012b rename, THM-016 dedup, T096, T097, T098, INV-038-041, interlacing-clique-deletion.md draft
**Unresolved threads:**
- verify_all_theorems.py stuck on n=7 OCF (find_odd_cycles slow at n=7)
- Quasi-regularity of Omega: needs theoretical explanation
- 20 of 30 theorem files still lack verification script references

## kind-pasteur-2026-03-06-S18e — 2026-03-06 (SC maximizer n=7, anti-aut involution theorem)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S18d
**Summary of work:**
  (1) **SC MAXIMIZER CONFIRMED AT n=7 (exhaustive)**: All 15 self-complementary score classes at n=7 have SC tournament achieving max H. Global max H=189 at regular (3,3,3,3,3,3,3) = Paley T_7. Key refinement: only self-complementary score sequences (s_i+s_{n-1-i}=n-1) can contain SC tournaments; the conjecture is vacuous for non-SC score sequences. 647s computation over 2^21 = 2M tournaments.
  (2) **ANTI-AUTOMORPHISM INVOLUTION THEOREM**: All anti-automorphisms of SC tournaments at n=5,6 are involutions (sigma^2=id). At even n=6: sigma is fixed-point-free with 3 two-cycles, creating 4 vertex-disjoint 3-cycle pairs. At odd n=5: 1 fixed point + 2 two-cycles. The involution structure is the algebraic mechanism behind the SC advantage.
  (3) **DISJOINT PAIR MECHANISM**: SC tournament at n=6 has 8 three-cycles all selecting one vertex from each sigma-orbit, forming 4 complementary pairs covering all 6 vertices. NSC tournament with same score has only 1 disjoint pair. This produces alpha_2=4 vs alpha_2=1 in the independence polynomial, explaining H=45 vs H=43.
**New contributions:** T091 (updated), T093 (updated), T095, sc_maximizer_n7_fast.py, anti_aut_analysis.py
**Unresolved threads:**
- Prove SC maximizer algebraically (involution mechanism is clear, need formal proof)
- Test at n=8 (even n, can have fpf involution, expect strong result)
- Is the anti-automorphism ALWAYS an involution? Or only at small n?

## opus-2026-03-06-S3 — 2026-03-06 (Deep census, blueself parity theorem, POS uniformity)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S2 (context continuation)
**Files read:** Full warm-up sequence, inbox messages from kind-pasteur S17/S18, new computation scripts
**Summary of work:**
  (1) **THM-023 PROVED: Blueself requires even n.** Flip changes endpoint scores by score'(0) = n - score(0), so same-score requires score(0) = n/2 (integer only at even n). Verified computationally at n=3,...,8.
  (2) **Deep census analysis (exhaustive n=3,...,6):**
    - POS orientation perfectly UNIFORM across all tilings AND within GS tilings
    - SC always maximizes H within each score sequence class (confirms kind-pasteur S18 finding)
    - Blueself achieves max or near-max H (rank 1 at n=4, ranks 1,5 at n=6)
    - Blackself at even n exclusively in NSC (paired) classes; at odd n in SC classes
    - SF tilings come in flip-pairs; count per class varies (2 at n=6, 4 at n=5)
    - Self-flip fraction decreasing: 25%, 12.5%, 1.56% at n=4,5,6
  (3) **Score(0) distribution in GS tilings:** Binomial-like, centered at (n-1)/2. At n=8: 1280/4096 have score(0)=4=n/2 (blueself-eligible).
  (4) **Background tasks launched:** n=7 full census (deep_census_analysis.py), n=8 approximate census (census_n8.py), n=7 POS census (pos_tiling_census.py). All still running.
**New contributions:** THM-023, INV-038, deep_census_analysis.py, census_n8.py, pos_tiling_census.py
**Unresolved threads:**
- n=7 deep census (background, very slow canonicalization ~hours)
- n=8 approximate census (background, ~2 hours ETA)
- Prove why blueself maximizes H (structural reason needed)
- Count blueself at n=8 (1280 eligible, need canonicalization)

## kind-pasteur-2026-03-06-S18 — 2026-03-06 (SC maximizer theorem, interlacing exhaustive n=5,6)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S17
**Files read:** THM-022, opus inbox MSG-015, TANGENTS.md, OPEN-QUESTIONS.md
**Summary of work:**
  (1) **SC MAXIMIZER THEOREM (exhaustive n=4,5,6)**: Within each score sequence class, the isomorphism class with maximum H is ALWAYS self-converse. At n=6 score (3,3,3,2,2,2): SC max=45 > NSC max=43. At (4,3,3,2,2,1): SC max=37 > NSC max=31.
  (2) **INTERLACING EXHAUSTIVE**: I(Omega(T-v),x) interlaces I(Omega(T),x) — 0/5120 at n=5, 0/196608 at n=6. EXHAUSTIVE.
  (3) **Cycle analysis**: Within score class (3,3,3,2,2,2), all classes have c3=8 but SC has more disjoint pairs (alpha_2=4 vs 1).
  (4) **Score regularity**: Low score variance correlates with high H. SC -> regular scores -> high H.
  (5) **Paired class FALSE**: "Paired minimize H" is false (transitive H=1 is SC). But "SC max within score class" is TRUE.
**New contributions:** T091, T092, paired_class_H_analysis.py, sc_maximizer_proof.py, interlacing_fast.py
**Unresolved threads:**
- Interlacing at n>=7 (FAILS — 36/1113 failures, corrected in T084)
- Prove SC maximizer algebraically
- Blueself = highest H among SC at even n?

## opus-2026-03-06-S16 — 2026-03-06 (Tiling model ↔ independence polynomial: blueself/blackself theorems)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S15 (context continuation)
**Files read:** Full warm-up sequence, tournament-tiling-explorer.html, definitions.md, MISTAKES.md
**Summary of work:**
  (1) **Resolved merge conflict** from previous session's push (agents/.session-state.json).
  (2) **Deep investigation** of relationship between tiling model structures (flip, grid-symmetry, transpose) and the independence polynomial I(Omega(T), x). Wrote 4 computational scripts totaling ~600 lines.
  (3) **PROVED Theorem 1 (THM-022):** T^op preserves I(Omega(T), x) — cycle reversal argument.
  (4) **PROVED Theorem 2:** Flip and transpose commute (pointwise ops commute with permutations).
  (5) **PROVED Theorem 3:** Grid-symmetry is a flip-invariant (corollary of Thm 2).
  (6) **PROVED Theorem 4:** For grid-symmetric tilings, k_0 + k_{n-1} = n-2 (endpoint constraint).
  (7) **PROVED Theorem 5 (odd-n obstruction):** No blueself tilings exist at odd n. Score-sequence product constraint: endpoint products differ by n-2-2k_0, which is nonzero at odd n. Exhaustive verification: 0/16 (n=5), 0/512 (n=7).
  (8) **VERIFIED Theorem 6:** Blueself/blackself mutual exclusivity at class level for n<=7. All 30 self-flip classes at n=7 are pure blackself. At n=6: 2 blueself, 6 blackself, 0 mixed.
  (9) **DISCOVERED:** flip(T) is NOT isomorphic to T^op in general (only 3.9% at n=6). Flip does not preserve H or I.P.
  (10) **Filed THM-022** with 7 theorems/observations. Added tangents T087-T090.
**New contributions:**
- THM-022 (tiling-indpoly structure theorems)
- T087 (flip-invariance of grid-symmetry), T088 (odd-n obstruction), T089 (mutual exclusivity), T090 (flip ≠ T^op)
- 04-computation/tiling_indpoly_investigation.py (main investigation)
- 04-computation/tiling_proofs.py (focused proofs)
- 04-computation/mutual_exclusivity_n7.py (n=7 verification)
- 04-computation/blueself_parity.py (parity analysis)
**Unresolved threads:**
- Mutual exclusivity at n>=8 (even n, non-trivial)
- Algebraic proof of odd-n obstruction for ALL odd n (current proof is exhaustive for n<=7)
- Connection between blueself classes and H-maximization (blueself classes have high H, regular scores)
## opus-2026-03-06-S2 — 2026-03-06 (Pin-grid sigma structure: POS-free identity, two-sigma analysis, induction framework)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S1 (context continuation)
**Files read:** SESSION-LOG.md, INVESTIGATION-BACKLOG.md, symmetry_check.py, definitions.md
**Summary of work:**
  (1) **PIN-GRID SIGMA VERIFIED**: Pin-grid sigma (r,c)->(c,r) acts WITHIN strips (same end vertex). Tournament sigma (i,j)->(n-1-j,n-1-i) acts ACROSS strips. They agree only on diagonal r=c.
  (2) **POS-FREE IDENTITY PROVED**: free(strip k) = cumul_POS(k) = floor(k/2). Algebraically: delta_free(k) = POS(k) = [k even]. This confirms user's observation that per-strip sigma-free bits grow at the same rate as POS.
  (3) **n->n+2 INDUCTION STRUCTURE**: Adding two strips adds exactly n sigma-free bits and exactly 1 POS. The single POS arc is always the midpoint arc: from vertex floor(k/2)-1 to vertex k-1 (0-indexed). Identity: floor(n/2) + floor((n+1)/2) = n.
  (4) **TWO-SIGMA ANALYSIS**: Tournament sigma ALWAYS preserves H (100%, verified n=3,...,7). Pin-grid sigma does NOT preserve H (only 5% at n=7). The two sigmas do NOT commute; their composition has order 3, generating a group isomorphic to S_3.
  (5) **MOD-4 ANALYSIS**: Pin-grid sigma preserves H mod 2 trivially (Redei). Does NOT preserve H mod 4 (~50-58%). Flipping sigma pairs does NOT preserve H mod 4. The POS/parity connection operates through tournament structure, not through simple modular arithmetic.
  (6) **TRANSFER MATRIX**: The M[a,b] transfer matrix symmetry (INV-001) is NOT the raw path count — it involves inclusion-exclusion: M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S). Sum over extensions does not preserve this symmetry.
**New contributions:**
- 04-computation/sigma_structure.py (verified all structural claims n=3,...,19)
**Unresolved threads:**
- How does the POS structure connect to OCF proof? Pin sigma doesn't preserve H, so simple induction via sigma pairs fails
- The group generated by pin sigma and tournament sigma (order 6) — does it have algebraic significance?
- Transfer matrix symmetry (INV-001) remains the key target

## kind-pasteur-2026-03-06-S17 — 2026-03-06 (Integrate background agents: H(T_19), web research, interlacing discovery)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S16 (context continuation)
**Files read:** SESSION-LOG.md, TANGENTS.md, OPEN-QUESTIONS.md, THM-020, background agent outputs
**Summary of work:**
  (1) Resolved git merge conflicts from S16 push (TANGENTS.md, SESSION-LOG.md, inbox files, rename conflict) — renumbered kind-pasteur T074-T078 to T079-T083 to avoid collision with opus T074-T075
  (2) Integrated 4 background agent results:
    - **H(T_19) confirmed**: 1,172,695,746,915, |Aut|=171, H/|Aut|=6,857,869,865 (exact integer). All 19 per-endpoint counts equal (vertex-transitive).
    - **Web research**: arXiv:2505.22766 polymer-gas representation new; Sokal conjecture (Peters-Regts) only helps at Delta<=3; Prakash-Sharma root gap doesn't give real-rootedness; no published work on Omega(T) independence polynomial.
    - **Independence poly analysis**: Full Omega(T_11) has 21169 cycles, I=[1,21169,10879,1155], degree 3, all real roots. INTERLACING under vertex deletion: 0/1324 failures — major new lead!
    - **J(n,3) large n**: 0/50 failures at n=12,15 for random induced subgraphs (failures only appear at n=9 with small subsets).
  (3) Added tangents T084 (interlacing), T085 (full Omega T_11), T086 (polymer-gas)
  (4) Updated OPEN-Q-015 with interlacing lead, web research results, full Omega data
**New contributions:** T084, T085, T086; updated OPEN-Q-015
**Unresolved threads:**
- INTERLACING is the most promising lead for proving real roots by induction. Need to formalize and test more.
- Can the Omega(T) -> Omega(T-v) relationship be made explicit enough to prove interlacing?
- H(T_31) computation (next Paley prime)

## opus-2026-03-06-S15 — 2026-03-06 (Elementary real-rootedness proof via Turán's theorem + discriminant)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S14 (continuation after merge conflict resolution)
**Files read:** Full warm-up sequence, all inbox messages, INVESTIGATION-BACKLOG.md, human inbox (tournament-tiling-explorer.html)
**Summary of work:**
  (1) **RESOLVED MERGE CONFLICTS** from previous session's failed push. Three files (TANGENTS.md, SESSION-LOG.md, .session-state.json) had conflict markers from concurrent kind-pasteur-S15 and opus-S13 commits. All resolved: kept both sides, renamed my session S12→S14 to avoid collision. Successfully pushed.
  (2) **INTEGRATED HUMAN INBOX**: Tournament Tiling Explorer (interactive HTML visualization tool). Archived to 03-artifacts/visualizations/. Extracted perspective conjecture (T075): #vertex-orbits at n predicts #classes at n+1, holds for n=3→4 and n=4→5, fails at n=5→6.
  (3) **NEW THEOREM THM-021: ELEMENTARY PROOF OF REAL-ROOTEDNESS** for n≤8 via discriminant analysis + Turán's theorem. Key insight: for n≤8, alpha(Omega)≤2, so I(Omega,x) = 1+a₁x+a₂x² has real roots iff a₁²≥4a₂. The "disjoint 3-cycle graph" is triangle-free (three disjoint triples need 9>n vertices), so by Turán: a₂≤c₃²/4, giving discriminant≥0. At n=8, 3-5 pairs handled by AM-GM ((c₅-1)²≥0). This is an INDEPENDENT proof of THM-020, using only classical combinatorics instead of Chudnovsky-Seymour.
  (4) **EXHAUSTIVE VERIFICATION**: Discriminant bound verified for ALL tournaments at n=5 (1024), n=6 (32768), n=7 (2,097,152). Sampled at n=8 (2000) and n=9 (500). Zero failures.
  (5) **TIGHTNESS ANALYSIS**: The bound is TIGHT at n=6,7 — tournaments with exactly 2 complementary 3-cycles achieve disc=0, giving I=(1+x)² (double root at -1). These form a single isomorphism class at n=6 with score sequence (4,4,4,1,1,1).
**New contributions:**
- THM-021 (discriminant real-rootedness proof)
- T074 (discriminant approach tangent), T075 (perspective conjecture from HTML tool)
- 04-computation/discriminant_real_roots.py (exhaustive verification script)
- 03-artifacts/visualizations/tournament-tiling-explorer.html (human's visualization tool)
- OPEN-Q-015 updated with alternative proof reference
**Unresolved threads:**
- For n≥9, degree≥3: need Newton/ULC bounds to extend real-rootedness proof
- The perspective conjecture (T075): why does it fail at n=5→6? Is there a corrected version?
- Turán bound is tight: can the double-root tournaments (disc=0) be characterized algebraically?

## opus-2026-03-06-S1 — 2026-03-06 (Deep tiling geometry investigation)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S13 (backlog investigation)
**Files read:** INVESTIGATION-BACKLOG.md, SESSION-LOG.md, OPEN-QUESTIONS.md, definitions.md, backlog-investigation-results.md, knuth-connections-research.md, tiling-class-structure.md, parity_tournaments_fixed.tex
**Summary of work:**
  (1) **SIGMA SYMMETRY**: Proved sigma (converse + relabel) acts cleanly on tiling classes. sigma PRESERVES bits (permutes, no complement) and weight. Sigma-fixed tilings = 2^floor((n-1)^2/4). Self-converse classes: 2,2,8,12 at n=3,4,5,6.
  (2) **COMPLEMENT DOES NOT RESPECT CLASSES**: Unlike sigma, bit-complement (flip all non-path arcs) does NOT map classes to classes. This asymmetry is fundamental.
  (3) **TRIANGLE 3-CYCLE FORMULA**: P(3-cycle) = 1/2 for consecutive triples, 1/4 for all others. E[c3] = (C(n,3) + n-2)/4 over random tilings.
  (4) **STRONG H~c3 CORRELATION**: r = 0.956 at n=5,6. H = 1+2c3 exact at n<=4, breaks at n>=5 due to higher odd cycles.
  (5) **CLASS TRANSITION GRAPH**: Always connected. ΔH always even (Rédei). E[ΔH] = 0 for every arc position. Self-loop fraction: 42%, 10%, 5% at n=4,5,6.
  (6) **INDISTINGUISHABLE CLASSES**: At n=6, classes sharing all tournament invariants are either sigma pairs (converse twins), self-converse coincidences, or cross-score partners. Weight distribution can distinguish beyond tournament invariants.
  (7) **BIT-POSITION VARIANCE**: Longest arc (gap=n-1) is most predictive of class. Middle arcs vary most within classes.
**New contributions:**
- INV-036 (tiling grid geometry), tiling-symmetry-analysis.md, 5 computation scripts
**Unresolved threads:**
- Find grid-local rules predicting class membership beyond c3
- Connect sigma reduction to arc-flip proof strategy (INV-004)
- Higher odd cycle contribution as correction to H ≈ 1+2c3

## kind-pasteur-2026-03-05-S16 — 2026-03-05 (Real roots to degree 6, perpendicular bell curve, J(n,3) NOT hereditary)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S15 (context continuation)
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, INVESTIGATION-BACKLOG.md, THM-020-real-roots.md, knuth-connections-research.md, tournament_lib.py, Irving-Omar paper (arXiv HTML)
**Summary of work:**
  (1) **RESOLVED GIT REBASE**: Merged S15 commit with opus-S12, pushed successfully.
  (2) **PERPENDICULAR BELL CURVE**: H(T) forms a perfect symmetric bell curve vs tiling Hamming weight at n=5 (exhaustive). Avg H by HW: 1.0, 3.2, 5.3, 6.9, 8.0, 8.4, 8.0, 6.9, 5.3, 3.2, 1.0. Ratio perp/diag = 2.79. NSC avg H is CONSTANT (~4.5) across all HW, while SC avg varies 1.0-9.7. NSC H values severely constrained: {3,5} at n=5 vs SC {1,3,9,11,13,15}.
  (3) **REAL ROOTS EXTENDED TO DEGREE 6**: Verified 0 failures for I(Omega_3(T), x) at n=5 (exhaustive, deg 1), n=6 (500, deg 1-2), n=7 (200, deg 1-2), n=8-10 (50 each, deg 2-3), n=12 (30, deg 3-4), n=15 (10, deg 5), n=20 (3 confirmed, deg 6). All roots real and negative.
  (4) **ULC CONFIRMED**: Ultra-log-concavity (alpha_k/C(m,k) log-concave) holds with 0 failures at ALL tested n.
  (5) **J(n,3) HEREDITARY REAL-ROOTEDNESS IS FALSE**: Arbitrary induced subgraphs of Johnson graph J(n,3) do NOT always have real-rooted I.P. Counterexample: 7 triples at n=9, alpha=[1,7,5,1], discriminant=-44. About 1.1% of random 7-element subsets fail. BUT tournament-realizable triple sets ALWAYS give real roots. This proves real-rootedness is tournament-specific, not a general J(n,3) property.
  (6) **HOOK EXPANSION ANALYSIS**: Irving-Omar Prop 26 hook coefficients [s_{(i,1^{n-i})}]U_T are palindromic for all tournaments. Transitive hooks = binomial coefficients C(n-1,k). At n=4: exactly 3 hook patterns (H=1→[1,3,3,1], H=3→[3,3,3,3], H=5→[5,3,3,5]).
  (7) **IRVING-OMAR PROOF CHAIN**: Fully understood Corollary 20 proof. S(D) = permutations whose nontrivial cycles are D-cycles. 2^psi = 2^(#nontrivial cycles). Each directed cycle gives one permutation cycle. So sum = I(Omega(D), 2).
  (8) **SCHWESER-STIEBITZ-TOFT**: Only mod-2 results. No new mod-4 information beyond what we have.
  (9) Processed inbox: tournament_n5_v5.html (interactive tiling explorer visualization tool, no new theorems).
**New contributions:**
- T074 (perpendicular bell curve), T075 (NSC constrained), T076 (hook palindrome), T077 (degree-6 real roots at n=20), T078 (J(n,3) NOT hereditary)
- Updated THM-020 verification table with polynomial degrees
- 04-computation/perpendicular_propagation.py, real_roots_deep.py, hook_expansion_test.py
**Unresolved threads:**
- H(T_19) computation (agent launched, may be running)
- Why does tournament structure prevent non-real roots? The constraint is: cyclic triples of a tournament form a specific subset of J(n,3). What algebraic property of this subset forces real-rootedness?
- The perpendicular bell curve: does it hold for all n? What is the closed-form relationship between HW and avg H?
- Web research agent for recent real-roots results (launched, may be running)

## kind-pasteur-2026-03-05-S15 — 2026-03-05 (Perpendicular maximizer, Claude's Cycles, ultra-log-concavity)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S14b (same session)
**Files read:** n8-anomaly-deep-dive.md, tiling-class-structure.md, INVESTIGATION-BACKLOG.md, parity_tournaments_fixed.tex (Open Problems section), Irving-Omar paper (web)
**Summary of work:**
  (1) **PERPENDICULAR MAXIMIZER PRINCIPLE**: Self-converse tournaments (T ~ T^op) have significantly more Hamiltonian paths: SC/NSC mean ratio 1.97 at n=5 (exhaustive), 1.59 at n=7 (random). The global H-maximizer is ALWAYS self-converse. In tiling space, SC tournaments sit perpendicular to the transitive<->full diagonal. This is the geometric reason Paley tournaments (always SC) maximize H.
  (2) **ULTRA-LOG-CONCAVITY**: I(Omega(T), x) coefficients satisfy ULC (alpha_k/C(m,k) is log-concave) with 0 failures across n=5 (exhaustive), n=6-9 (random). This is STRONGER than Newton inequalities. Despite Omega NOT being a matroid (62% non-matroid at n=6), ULC holds.
  (3) **CLAUDE'S CYCLES CONNECTION**: Knuth's paper (Feb 2026) on m^3-vertex Hamiltonian cycle decomposition. Works for ODD m (Gray code construction), fails for EVEN m. Parallels our Paley maximizer working for odd primes and n=8 being the first even anomaly. Three coordinate axes (i,j,k) mirror three directions in tournament space.
  (4) **IRVING-OMAR BARRIER**: Schur expansion of Redei-Berge function is NOT positive for tournaments (counterexample in paper). Rules out Schur-positivity route to real-rootedness. But Corollary 20 confirms our OCF formula.
  (5) Revisited n=8 anomaly: BlackSelf(8) = T_657, which extends P(7), has uniform D_v=54, and (H-Fix)/2=312 (even). The mu explosion at n=8 is the tournament analog of Claude's Cycles failing at even m.
**New contributions:**
- T070 (perpendicular maximizer), T071 (ultra-log-concavity), T072 (Claude's Cycles), T073 (Irving-Omar barrier)
**Unresolved threads:**
- Can ULC be proved from tournament structure? (Would imply log-concavity, a step toward real roots)
- What is the precise relationship between odd/even dichotomy and perpendicular propagation?
- Does the "three-axis" analogy (Hamming weight / self-converse / automorphism) have algebraic content?

## opus-2026-03-06-S14 — 2026-03-06 (Rigorous verification of all claims, OPEN-Q fixes, Irving-Omar matrix algebra exploration)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S10 (continuation)
**Files read:** Full warm-up sequence, all inbox messages, INVESTIGATION-BACKLOG.md, all new session logs
**Summary of work:**
  (1) RIGOROUS VERIFICATION: Re-derived all critical computational claims from scratch.
    - OCF exhaustive n=3,4,5: PASS (0 failures with correct all-cycle enumeration)
    - Paley H values P(3)=3, P(7)=189, P(11)=95095: PASS
    - H(P(19))=1172695746915, H(P(23))=15760206976379349: PASS
    - Max H at n=7 = 189 achieved by 240 tournaments (exhaustive): PASS
    - T_11 OCF decomposition 1+2*21169+4*10879+8*1155=95095: PASS
    - No claws in Omega at n=6 (exhaustive): PASS
    - Claw-free argument for n<=8: VERIFIED (sound vertex-counting proof)
  (2) BUG FOUND in many auxiliary scripts: cycle-finding with `break` after first cycle per vertex set misses multiple directed cycles on same vertex set (e.g., two directed 5-cycles on 5 vertices). The CORE verification script (verify_core_results.py) is CORRECT (uses fix-minimum-vertex without break). This bug does NOT affect any verified theorems.
  (3) AUDIT AGENTS found 3 issues in OPEN-Q-013: H(T_19) still marked "unknown" (computed in S5/S10), attribution incomplete, next target obsolete. ALL FIXED.
  (4) AUDIT AGENTS found theorem files healthy with 2 minor issues (THM-018 scope clarification, THM-013 filename). No critical errors.
  (5) IRVING-OMAR PAPER VERIFIED: Corollary 20 says exactly what repo claims. Proof uses matrix exponential + Cayley transform of tournament adjacency matrix. CREATIVE LEAD: odd-cycle extraction via arctanh splitting connects to orthogonal/unitary spectral theory, potentially explaining real roots.
  (6) NUMBER-THEORETIC OBSERVATION: alpha_3(T_11) = 1155 = |Aut|*21. alpha values of Paley tournaments have divisibility by |Aut| for alpha_3 (but not alpha_1, alpha_2).
  (7) c_3(P(7)) = 14 (not 7). The value 7 was for undirected cyclic triples with a different normalization. 14 cyclic + 21 transitive = 35 = C(7,3). Consistent with p(p^2-1)/24 = 14.
**New contributions:**
- OPEN-Q-013 updated with H(T_19), H(T_23), a(8)=661 discovery, Szele ratio data
- T065 (Irving-Omar Cayley transform connection to real roots) in TANGENTS.md
- 04-computation/verify_all_claims.py (comprehensive verification script)
**Unresolved threads:**
- Irving-Omar Cayley transform: can I(Omega(T), x) be expressed as a determinant?
- Alpha_k divisibility by |Aut| for Paley tournaments — systematic study needed
- Real roots: algebraic path via symmetric functions remains the key open problem

## opus-2026-03-05-S13 — 2026-03-06 (Rigorous backlog investigation: 2-adic tower, full Omega real roots, coefficient analysis)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S12 (continuation)
**Files read:** INVESTIGATION-BACKLOG.md, OPEN-QUESTIONS.md, symmetry_check.py, paper-deep-connections.md
**Summary of work:**
  (1) **INV-014 (2-adic tower): PARTIALLY RESOLVED.** v_2(H(T)) = 0 for ALL tournaments (exhaustive n≤6, sampled n=7). This IS Redei's theorem. H mod 4 ≡ 1+2·alpha_1 (mod 4) via OCF. At n=3,4 this equals 1+2·c_3 exactly, but breaks at n≥5 when 5-cycles contribute. H mod 2^k approaches uniform on odd residues as n grows.
  (2) **OPEN-Q-015 (real roots): FULL OMEGA VERIFICATION.** Tested I(Omega(T), x) with ALL odd cycles (3+5+7+9-cycles), not just 3-cycles. 185 samples at n=5..9: 100% OCF pass, 100% real roots, 100% log-concave, 100% unimodal. Paley p=7: full Omega has 80 cycles, I.P.=[1,80,7], roots both real negative.
  (3) **Coefficient structure discovered.** Alpha_k/alpha_{k-1} ratios drop super-exponentially (~0.3 for k=1→2, ~0.03 for k=2→3). Root ratio ≈ alpha_1²/alpha_2 from Vieta formulas. Independence polynomial degree = floor(n/3) for full Omega.
  (4) **Source cone identity: H(source_cone(T')) = H(T').** Confirmed by INV-035 background agent. Trivial from path structure but important for R-cone proof strategy (INV-004).
  (5) Launched 3 background agents for INV-014, INV-013, INV-035. Two produced partial results; INV-014 agent got stuck on long computation.
**New contributions:**
- 03-artifacts/drafts/backlog-investigation-results.md (comprehensive results document)
- OPEN-Q-008 partially resolved (v_2=0 always)
- OPEN-Q-015 strengthened (full Omega verification)
- INV-014 status updated in backlog
**Unresolved threads:**
- Background agents may still be running (INV-013 realizable Omega, INV-035 tournament families)
- Full Omega at n≥10 needs more efficient cycle enumeration
- Algebraic explanation for real-rootedness remains the key open problem
- Irving-Omar symmetric function approach to real-rootedness not yet explored

## opus-2026-03-05-S12 — 2026-03-05 (Deep web research: Knuth connections to OCF — 7 channels identified)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S11 (continuation)
**Files read:** MSG-014 (kind-pasteur-S14b), SESSION-LOG.md, TANGENTS.md
**Summary of work:**
  (1) **KNUTH FASCICLE 8a** (Dec 2025): Section 7.2.2.4 of TAOCP on Hamiltonian paths and cycles. "Path flipping" = our arc-flip approach. Directly relevant to OCF proof strategy. Priority: read for tournament-specific exercises.
  (2) **CLAUDE'S CYCLES** (Knuth, Feb 2026): Hamiltonian cycle decomposition of Cayley digraphs solved by Claude Opus 4.6 via modular Gray code. 760 valid decompositions. Parallel: algebraic structure (Cayley/Paley) enables exact solutions.
  (3) **RÉDEI REVISITED** (arXiv:2510.10659, Oct 2025): Schweser-Stiebitz-Toft exhibit stronger Rédei theorems by Dirac and Berge. Dirac's theorem constrains parity of Hamiltonian paths in mixed graphs — could give mod-4 information beyond OCF.
  (4) **ZDD TECHNOLOGY**: Knuth's Zero-Suppressed Decision Diagrams could efficiently compute I(Omega(T), x) for moderate n, exploiting sparse complement of dense Omega.
  (5) **ALGORITHM X / DLX**: Exact cover formulation of odd-cycle packing. DLX could enumerate all vertex-disjoint odd-cycle collections, providing alternative OCF verification.
  (6) **PERMANENT / CYCLE COVERS**: Ryser's inclusion-exclusion for permanents is structurally parallel to independence polynomial evaluation. Bridge to Björklund's GF(2) method.
  (7) **PÓLYA ENUMERATION**: Aut(T_p) acting on odd-cycle collections. Cycle index of affine group could give structured formula for I(Omega(T_p), 2).
  (8) Read kind-pasteur-S14b findings: S_{2,1,1}-free at n<=9, quasi-line fails at n=8, line graph at n=6. Full hierarchy: line(n<=5) < quasi-line(n<=7) < claw-free(n<=8) < S_{2,1,1}-free(n<=9) < S_{1,1,1}-free(n<=11).
**New contributions:**
- 03-artifacts/drafts/knuth-connections-research.md (7 connections with priority ranking)
**Unresolved threads:**
- Read Knuth fasc8a.pdf for tournament exercises
- Read arXiv:2510.10659 (Rédei revisited) for stronger parity theorems
- Implement ZDD-based I(Omega(T), x) for n=8-12
- Formalize Cayley/Paley algebraic structure parallel

## opus-2026-03-05-S11 — 2026-03-05 (Continued research: Szele-Adler-Alon e-convergence, real roots at n=15, Paley asymptotic)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S9 (continuation)
**Files read:** All new messages (MSG-019/S10), TANGENTS.md, n8-anomaly-deep-dive.md, SESSION-LOG.md
**Summary of work:**
  (1) **SZELE BOUND DISCOVERY:** Adler-Alon-Ross (2001) proved max H(T) >= (e-o(1))·n!/2^{n-1}. Paley ratios H(P(p))/(p!/2^{p-1}) = 2.00, 2.40, 2.44, 2.53, 2.56 converge exactly toward e=2.718. Paley tournaments achieve the asymptotically optimal constant.
  (2) **REAL ROOTS EXTENDED to n=15:** 0/95 tournaments at n=9..15 have non-real roots for I(Omega_3, x). Crucially, this holds at n=12-15 where S_{1,1,1}-freeness FAILS. Alpha(Omega_3) grows as ~floor(n/3).
  (3) **FRIEDGUT-KAHN UPPER BOUND:** Improved Alon's n^{3/2} to n^{3/2-0.25} ≈ n^{1.25}. If Paley achieves ~e·n!/2^{n-1}, the polynomial factor may be entirely unnecessary.
  (4) **LINE GRAPH HYPOTHESIS ALREADY DISPROVED** in S9. K_5-e in Omega_3 at n=6.
  (5) Updated zero-free-subdivided-claw-research.md with Connections 7 and 8.
**New contributions:**
- Extended real-roots verification data (n=11..15)
- Szele-Adler-Alon analysis connecting Paley ratios to e
- Updated 03-artifacts/drafts/zero-free-subdivided-claw-research.md
**Unresolved threads:**
- Does H(P(p))/(p!/2^{p-1}) → e exactly? Needs H(P(31))
- Can Alon upper bound be improved to O(1)? (major open problem)
- Real roots: algebraic explanation needed (symmetric function path)
- Irving-Omar paper Corollary 20: exact statement still unverified

## kind-pasteur-2026-03-05-S14b — 2026-03-05 (Structural hierarchy: subdivided-claw-freeness, line graph refutation)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S14 (context continuation)
**Files read:** THM-020, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md
**Summary of work:**
  (1) T054 REFUTED: Omega(T) is NOT a line graph at n=6. K5-e (Beineke forbidden subgraph) in 53% of n=6 tournaments. Heilmann-Lieb cannot explain real roots.
  (2) NEW FINDING: Omega(T) is S_{2,1,1}-free (subdivided-claw-free) at n<=9 (0/100 failures), but FAILS at n>=10 (92%). Combined with opus-S9: S_{1,1,1}-free through n=11, fails at n=12. Each subdivision level buys ~3 vertices.
  (3) Quasi-line graph property FAILS at n=8 (49%), n=9 (10%), n=10 (0%).
  (4) Full structural hierarchy: line graph (n<=5) < quasi-line (n<=7) < claw-free (n<=8) < S_{2,1,1}-free (n<=9) < S_{1,1,1}-free (n<=11) < ??? Real roots hold at ALL tested n.
  (5) LITERATURE: Alon (1990) + Adler-Alon-Ross (2001): max H(T) = Theta(n!/2^{n-1}). Paley maximizer conjecture consistent.
  (6) LITERATURE: Jerrum-Patel (2026) proves real roots for bounded-degree H-free (H = subdivided claw). Key insight from opus-S9: every fixed subgraph eventually appears in Omega(T), so explanation must be algebraic.
  (7) Verified all real roots + log-concavity + unimodality at n=9 (10 random).
**New contributions:**
- T061-T064 (line graph refuted, S_{2,1,1}-free, quasi-line fails, Alon bounds)
- Updated THM-020 with expanded graph property table (including line graph, quasi-line, S_{2,1,1})
- Updated OPEN-Q-015 with comprehensive structural analysis
**Unresolved threads:**
- What structural property explains real roots at n>=10? (opus-S9 says: must be algebraic)
- Real roots: Irving-Omar/Grinberg-Stanley symmetric function framework is most promising path

## opus-2026-03-05-S10 — 2026-03-05 (Paley maximizer verification, n=8 H-maximizer discovery, Omega structure analysis)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S8 (context continuation after summary)
**Files read:** SESSION-LOG.md, TANGENTS.md, OPEN-QUESTIONS.md, THM-020, MSG-013, n8-anomaly-deep-dive.md, verify_n8_claims.py, INVESTIGATION-BACKLOG.md
**Summary of work:**
  (1) DISCOVERED: a(8)=661 from OEIS A038375. The H-maximizer at n=8 is a SC tournament with |Aut|=1, NOT T_657 (H=657, the Paley extension). T_661 does NOT contain P(7) as a vertex-deletion and has non-uniform D_v=(67,63,68,72,72,68,63,67).
  (2) CONFIRMED: P(7) is the GLOBAL H-maximizer at n=7, verified by exhaustive enumeration of all 2,097,152 tournaments on 7 vertices. 240 tournaments achieve H=189.
  (3) CONFIRMED: Paley maximizer conjecture at p=11: H(P(11))=95095=a(11).
  (4) NEW COMPUTATIONS: H(P(19))=1,172,695,746,915 and H(P(23))=15,760,206,976,379,349. These are candidate values for a(19) and a(23) in OEIS.
  (5) RATIO ANALYSIS: H(P(p))/(p!/2^{p-1}) = 2.00, 2.40, 2.44, 2.53, 2.56 for p=3,7,11,19,23 — converging toward e=2.718, consistent with Szele-Alon-Friedland asymptotic bound.
  (6) OMEGA STRUCTURE: Full Omega at n=8 has 76-78 vertices (density 0.98). The 3-cycle subgraph accounts for only 30% of H(T); 70% comes from 5-cycle and 7-cycle contributions. This is the structural signature of the n=8 transition.
  (7) IMPORTANT NEGATIVE: The Paley maximizer conjecture applies ONLY at Paley primes (p=3 mod 4). At non-Paley n=8, the maximizer has trivial automorphism group and no Paley substructure.
**New contributions:**
- T059 (n=8 maximizer and 5-cycle dominance), T060 (Paley ratio convergence to e) in TANGENTS.md
- T053 updated with p=19,23 data and n=8 negative result
- INV-024 extended through p=23, INV-032 updated with perfectness failure
- 04-computation/check_max_h8.py, paley_maximizer_test.py, paley_h_sequence.py, n8_omega_structure.py
**Unresolved threads:**
- Submit H(P(p)) values to OEIS (a(19), a(23) are new)
- Test line graph hypothesis for Omega(T) (Beineke forbidden subgraphs)
- Test subdivided-claw-freeness at n=9
- Compute H(P(31)) if feasible

## opus-2026-03-05-S9 — 2026-03-05 (Deep web research: zero-free regions, subdivided claws, real-roots mystery)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S8 (context continuation after compaction)
**Files read:** MEMORY.md, all inbox messages (MSG-010 through MSG-013, MSG-015), THM-020-real-roots.md, n8-anomaly-deep-dive.md, web-research-connections.md, tribonacci-web-research.md, grinberg-stanley-connection.md, MISTAKES.md, THM-002, THM-019, INVESTIGATION-BACKLOG.md, hard-core-statistical-physics-connections.md
**Summary of work:**
  (1) SYSTEMATIC AUDIT: Found and fixed errors in tribonacci-web-research.md (4 sections with false claw-free/perfect claims), web-research-connections.md (Connection 4), INVESTIGATION-BACKLOG (INV-032).
  (2) BUG FOUND: hard_core_fast.py chordality check (lines 157-168) only checks 1 of 3 C4 orderings. Reports 0% non-chordal when correct is ~14% at n=6.
  (3) **NEW COMPUTATIONAL: S_{1,1,1}-freeness** — Omega_3(T) is S_{1,1,1}-free (subdivided claw with each edge subdivided once) through n=11, FAILS at n=12 (100%). Each subdivision level buys ~3 more vertices. Pattern: claw-free n<=8, S_{1,1,1}-free n<=11.
  (4) **LINE GRAPH HYPOTHESIS DISPROVED:** K_5-e (Beineke forbidden subgraph) found in Omega_3(T) at n=6 (45%). Not a line graph.
  (5) **JERRUM-PATEL 2026 ANALYSIS:** Subdivided claws are the EXACT boundary for zero-free regions of I(G,x) in bounded-degree H-free graphs. Best possible: non-subdivided-claw H allows zeros on positive real axis. BUT: requires bounded degree AND fixed H, both of which fail for Omega(T) at large n.
  (6) **BEZAKOVA et al. 2024:** Matching mixing-time dichotomy — O(n log n) for S_{a,b,c}-free bounded-degree, exponential otherwise. Same limitations.
  (7) **KEY INSIGHT: Real-rootedness cannot be explained by any forbidden subgraph property.** Every fixed subgraph eventually appears in Omega(T). Explanation must be algebraic, likely through Irving-Omar/Grinberg-Stanley symmetric function framework.
  (8) **AUTHOR CORRECTION:** arXiv:2412.10572 is by Irving & Omar, not Grinberg & Stanley. They build on G-S's framework.
  (9) **SZELE BOUND:** H(T) <= c·n^{3/2}·n!/2^{n-1} (Alon 1990, solving Szele 1943 conjecture). Lower bound n!/2^{n-1} achieved probabilistically. Paley maximizer conjecture (kind-pasteur-S14b) says T_p achieves the max; ratio H(T_p)/(p!/2^{p-1}) ≈ 2, 2.4, 2.44 for p=3,7,11.
**New contributions:**
- 03-artifacts/drafts/zero-free-subdivided-claw-research.md (comprehensive research document)
- Fixed 4 errors in existing documentation files
- S_{1,1,1} computational data for n=9..12
**Unresolved threads:**
- Prove real roots for n>=9 — algebraic approach needed (symmetric functions?)
- Rate of growth of max subdivided claw in Omega(T) as function of n
- Is there a weaker graph property ("quasi-claw-free"?) that holds universally?
- Irving-Omar framework: does it say anything about I(Omega(T), x) beyond x=2?
- Compute H(T_19) for Paley maximizer conjecture
- Szele bound tightness: does H(T_p) ~ c·p!/2^{p-1} for some constant c?

## kind-pasteur-2026-03-05-S14 — 2026-03-05 (Deep research: real roots theorem, Paley maximizer, literature survey)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S13 (context continuation)
**Files read:** MEMORY.md, THM-019, OPEN-QUESTIONS.md, TANGENTS.md, SESSION-LOG.md
**Summary of work:**
  (1) NEW THEOREM (THM-020): All roots of I(Omega(T), x) are real and negative for n<=8. PROVED via claw-freeness + Chudnovsky-Seymour (2007). Real-rootedness conjectured for all n.
  (2) COMPUTATIONAL: Verified real roots at n=5 (exhaustive), n=6 (500), n=9 (50), n=10 (30) — 0 failures anywhere, even beyond claw-free range.
  (3) TESTED & DISPROVED: Omega(T) is NOT always a comparability graph (fails n=7: 1%, n=8: 58%).
  (4) TESTED: det(I+2A) =/= H(T) in general. No simple spectral formula.
  (5) FIXED: Cycle enumeration bug — multiple directed cycles per vertex set.
  (6) **MAJOR: Paley tournaments MAXIMIZE H(T)!** OEIS A038375 confirms: a(3)=3=H(T_3), a(7)=189=H(T_7), a(11)=95095=H(T_11). New conjecture: T_p maximizes Hamiltonian paths for Paley primes.
  (7) LITERATURE (5 background agents, all completed):
    - 7 papers cite Grinberg-Stanley, including chromatic-Redei-Berge connection (arXiv:2506.08841)
    - Bezakova et al. 2024 dichotomy for H-free graphs and hard-core model
    - Omega(T) NOT interval graph at n=6 (13.9% fail)
    - Line graph hypothesis: if Omega(T) = L(H), Heilmann-Lieb gives real roots
    - H(T_p) and H/|Aut| sequences NOT in OEIS
    - El Sahili-Abi Aad (2019) extends Forcade parity; Grunbaum conjecture proved
    - Our conflict graph / independence polynomial formulation is NOT in any G-S paper — genuinely original
    - Cycle divisibility threshold corrected: k >= (p+3)/2, not (p+1)/2
  (8) FINDING: lambda=2 is OUTSIDE hard-core uniqueness for Omega(T). Real roots are special.
  (9) Graph property hierarchy: interval (n<=5) < chordal (n<=5) < comparability (n<=6) < perfect (n<=7) < claw-free (n<=8).
**New contributions:**
- THM-020-real-roots.md (new theorem)
- OPEN-Q-015 (real roots conjecture)
- T048-T058 in TANGENTS.md (Paley maximizer, line graph, interval graph, LGV connection, subdivided claw, chromatic-RB connection)
- Corrected cycle divisibility threshold in OPEN-Q-013
**Unresolved threads:**
- Test line graph hypothesis (Beineke's 9 forbidden subgraphs)
- Test subdivided-claw-freeness at n=9
- Compute H(T_19) — would test both Paley maximizer conjecture and extend the sequence
- Submit H(T_p) to OEIS
- Bijective proof of OCF still open

## kind-pasteur-2026-03-05-S13 — 2026-03-05 (Web research: hard-core model, independence polynomial at lambda=2, statistical physics connections)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S12 (new session)
**Files read:** All warm-up files, inbox messages (MSG-007, MSG-013, MSG-014 x2), tribonacci-web-research.md
**Summary of work:**
  (1) Comprehensive web research on connections between H(T) = I(Omega(T), 2) and statistical physics.
  (2) FINDING: lambda=2 is in the NON-PERTURBATIVE regime of the hard-core model. The uniqueness threshold lambda_c(Delta) < 2 for all Delta >= 5. Cluster expansion diverges; approximate counting is #P-hard. OCF is a rare exact evaluation in this regime.
  (3) FINDING: I(G, 2) counts 2-labeled independent sets: pairs (S, f) where S is independent and f: S -> {1,2}. For OCF, each Ham path corresponds to a colored collection of vertex-disjoint odd cycles.
  (4) FINDING: Scott-Sokal (2005) connects LLL to independence polynomial zero-free regions. The Shearer zero-free disc is far too small to reach lambda=2 for any nontrivial max degree. Positivity of I(Omega(T), 2) cannot be deduced from LLL; it requires the combinatorial identity H(T) >= 1.
  (5) FINDING: Chudnovsky-Seymour (2007) gives all-real-roots for claw-free graphs, applicable to Omega(T) only for n <= 8. Jerrum-Patel (2026) extends zero-free regions to H-free classes.
  (6) Confirmed Grinberg-Stanley proof (arXiv:2307.05569 + arXiv:2412.10572) uses matrix algebra, NOT statistical mechanics. No bijective interpretation exists yet.
  (7) Created comprehensive research report: hard-core-statistical-physics-connections.md with 11 references.
**New contributions:**
- 03-artifacts/drafts/hard-core-statistical-physics-connections.md (comprehensive research report)
**Unresolved threads:**
- Bijective proof of OCF: each Ham path <-> 2-colored independent cycle collection
- Whether Omega(T) has special structure beyond claw-freeness that explains positivity at lambda=2
- Potential use of transfer matrix methods from stat mech for tournament problems

## opus-2026-03-05-S8 — 2026-03-05 (Deep n=8 anomaly investigation: BlackSelf, mu explosion, Paley connection)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S7
**Files read:** THM-019, tiling-class-structure.md, FINAL_FINDINGS.md, THM-016, INVESTIGATION-BACKLOG.md, omega_claw_fast.py, q009_adj_decomposition.py, tex paper (oq:n8 definition)
**Summary of work:**
  (1) Investigated the n=8 anomaly at the user's request, connecting blueself/blackself/pos formulas.
  (2) Exhaustive enumeration of all 65536 self-converse n=8 tournaments with reversal anti-automorphism.
  (3) Found 2560 with |Aut|>1, surveyed all: 10 distinct (H, Fix(beta), |Aut|) types.
  (4) DISCOVERED: T_657 (H=657, |Aut|=3) contains Paley tournament P(7) as vertex-deletion — the Hadamard connection for INV-012.
  (5) DISCOVERED: T_657 has perfectly uniform D_v = 54 for all vertices (unique among specimens).
  (6) DISCOVERED: mu values at n=8 reach {3, 5, 7, 9, 11} for 3-cycles and {1, 3} for 5-cycles.
  (7) DISCOVERED: Full Omega at n=8 has 76 vertices (20+48+8 cycles of length 3,5,7).
  (8) IDENTIFIED: T_A (|Aut|=9, H=621) is the ONLY SC+Aut>1 specimen with no C5 in Omega_3.
  (9) CORRECTED: Signed position identity requires T and T' (flipped tournament), not same T on both sides.
  (10) Web research: Chudnovsky-Seymour, El Sahili parity results, doubly regular tournaments ↔ skew Hadamard.
**New contributions:**
- `03-artifacts/drafts/n8-anomaly-deep-dive.md` — comprehensive analysis
- `04-computation/blackself8_*.py` — four investigation scripts
- INV-012 updated with detailed findings
**Unresolved threads:**
- Confirm T_657 is isomorphic to a known construction
- Resolve BlackSelf(8) definition ambiguity
- Full Omega independence polynomial at n=8 (76 vertices, needs specialized algorithm)
- Investigate whether uniform D_v characterizes Paley extensions

## opus-2026-03-05-S7 — 2026-03-05 (Omega perfectness DISPROVED at n=8; claw-freeness trivial at n<=8)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S6 (context continuation x3)
**Files read:** THM-019, grinberg-stanley-connection.md, tiling-class-structure.md, OPEN-QUESTIONS.md, SESSION-LOG.md, all inbox messages, omega_perfectness_implications.py
**Summary of work:**
  (1) Reviewed all new material from git pull: OCF proved by Grinberg-Stanley, THM-019, Tribonacci, tiling class structure.
  (2) Committed endpoint decomposition / linearity work (q009_endpoint_decomp.py, q009_G_identity.py, q009_G_in_test.py, q009_linearity_proof.py) as historical contributions.
  (3) PROVED: Omega(T) is trivially claw-free for n<=8 by vertex counting (3 pairwise v.d. odd cycles + 1 touching all three requires >= 9 vertices).
  (4) DISPROVED: Omega(T) is NOT always claw-free at n=9 (90% of random tournaments have claws).
  (5) DISPROVED: Omega(T) is NOT always perfect at n=8. 53.8% of random n=8 tournaments have a C5 (5-hole) in the 3-cycle conflict graph. Explicit counterexample constructed with 13 forced arcs.
  (6) VERIFIED: Omega(T) appears perfect at n<=7 (0/1000 failures).
  (7) CONFIRMED: OCF holds regardless of Omega perfectness (H=I(Omega,2)=417 at the n=8 counterexample).
**New contributions:**
- THM-019 CORRECTED: perfectness fails at n=8, holds for n<=7
- OPEN-Q-014 RESOLVED (disproved)
- 04-computation/omega_claw_fast.py (claw-freeness analysis)
- 04-computation/omega_c5_test.py (C5 construction at n=8)
- 04-computation/omega_claw_free_test.py (comprehensive test, slow)
- Vertex counting theorem: claw impossible for n<=8
**Unresolved threads:**
- Do all-real-roots / log-concavity of I(Omega,x) survive at n>=8 despite imperfection?
- What structural property of Omega(T) DOES explain OCF? (not perfectness, not claw-free)
- Prove Omega(T) perfectness for n<=7 (clean vertex-count argument for C5?)

## kind-pasteur-2026-03-05-S12 — 2026-03-05 (CRITICAL: OCF proved by Grinberg-Stanley; comprehensive audit)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S11 (context continuation)
**Files read:** All warm-up files, all theorem files (THM-001 through THM-018, CONJ-001, CONJ-002, PROP-001, LEM-001, LEM-002), MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, tiling-class-structure.md, signed-adjacency-identity.md, paper-deep-connections.md
**Summary of work:**
  (1) COMPREHENSIVE AUDIT of all mathematical results in the project. Read every theorem file, checked proof chains, cross-referenced claims.
  (2) FOUND AND FIXED: THM-016 significance section incorrectly claimed even-odd split implies OCF (contradicts MISTAKE-008). Fixed.
  (3) FOUND AND FIXED: signed-adjacency-identity.md status line said "equivalent to OCF". Fixed.
  (4) FOUND: THM-013 numbering collision (THM-013-insertion-decomposition.md says "THM-012" inside). Fixed header.
  (5) FOUND: Two duplicate THM-016 files (claim-b and hamiltonian-alternating-sum). Noted for future merge.
  (6) FOUND: LEM-001 and LEM-002 outdated after CONJ-002 refutation. Updated with resolution notes.
  (7) VERIFIED numerically (via background agents): THM-016 inductive proof (all steps correct, no gaps), Tribonacci proof (all 5 components verified), OCF at n=3,4,5 (independent recomputation), H(T_11)=95095, Claim B at n=6.
  (8) **CRITICAL DISCOVERY**: Web search found that OCF (H(T) = I(Omega(T), 2)) is ALREADY PROVED in the literature by Grinberg & Stanley (arXiv:2307.05569, 2023; arXiv:2412.10572, 2024). Their Corollary 20 states ham(D̄) = Σ 2^{ψ(σ)} over permutations with all odd D-cycles. For tournaments, D̄ = D^op and ham(D^op) = ham(D), so this gives H(T) = I(Omega(T), 2) exactly.
  (9) Updated CONJ-001 (now PROVED), THM-002 (now PROVED for all n), PROP-001 (now PROVED), OPEN-Q-002 and OPEN-Q-009 (RESOLVED). Created grinberg-stanley-connection.md documenting the full equivalence.
  (10) Updated LEM-001, LEM-002 with computed values (h_QR=h_NQR=201).
  (11) Web searches also confirmed: Tribonacci connection is NEW (not in OEIS for tournaments), H(T_19) value is NEW, Omega(T) perfectness is NEW, the H/|Aut| sequence 1,9,1729 is NOT in OEIS.
**New contributions:**
- CONJ-001 → PROVED (via Grinberg-Stanley OCF + Claim B)
- THM-002 → PROVED for all n (three independent proofs documented)
- PROP-001 → PROVED
- OPEN-Q-002, OPEN-Q-009 → RESOLVED
- 03-artifacts/drafts/grinberg-stanley-connection.md (full analysis of equivalence)
- Error fixes: THM-016, signed-adjacency-identity.md, THM-013 header, LEM-001, LEM-002
- 04-computation/verify_core_results.py (independent OCF verification)
- 04-computation/verify_tribonacci_proof.py (tribonacci proof verification)
**Unresolved threads:**
- Merge duplicate THM-016 files
- Canonize Omega(T) perfectness (potentially new theorem)
- The project's independent contributions (THM-016/017 even-odd split, THM-018 coefficient identity, tiling class structure, Tribonacci theorem) remain valuable and NEW
- Consider submitting Tribonacci result to OEIS
- The Grinberg-Stanley proof technique (symmetric functions / matrix algebra) is entirely different from our combinatorial methods — worth studying for further applications

## opus-2026-03-05-S6 — 2026-03-05 (Tribonacci web research: interval graph structure, Chudnovsky-Seymour connection)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S5 (context continuation)
**Files read:** tiling-class-structure.md, full_class_analysis.py, web-research-connections.md, hard_core_fast.py, paper-deep-connections.md
**Summary of work:** (1) Deep web research with Tribonacci structure in mind. (2) DISCOVERED that all odd cycles of T_full_n are consecutive intervals — Omega(T_full_n) is an INTERVAL GRAPH. (3) PROVED I(Omega(T_full_n), 2) satisfies Tribonacci recurrence independently of H(T_full_n) via weighted interval packing DP with telescoping. (4) Connected to Chudnovsky-Seymour: Omega(T) claw-free implies I(Omega(T), x) has ALL REAL (negative) ROOTS — log-concavity, positivity at x=2. (5) Connected to Chvátal-Sbihi: claw-free perfect graphs decompose via clique cutsets into elementary (line graphs of bipartite graphs) and peculiar atoms. (6) Identified Jerrum-Patel 2026 paper extending zero-free regions. (7) Found multiple combinatorial interpretations of Tribonacci (tilings, compositions, binary sequences).
**New contributions:**
- 03-artifacts/drafts/tribonacci-web-research.md (comprehensive synthesis)
- INV-035 added to INVESTIGATION-BACKLOG (Tribonacci interval graph structure)
- INV-032 updated with Chudnovsky-Seymour connection
- Verified OCF for T_full family n=3,...,8 via interval graph computation
**Unresolved threads:**
- Prove Omega(T) is always claw-free (INV-032 — would unlock Chudnovsky-Seymour)
- Direct bijection between run decompositions and weighted interval packings
- Extend Tribonacci analysis to other tournament families
- Check transfer matrix eigenvalues for T_full

## opus-2026-03-05-S5b — 2026-03-05 (THM-018: alpha_w^H = alpha_w^I PROVED symbolically at n=4,...,7)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4d (context continuation x2)
**Files read:** All warm-up files, previous computation files, THM-013, THM-002
**Summary of work:** (1) Reformulated alpha_w^H as insertion decomposition: for each perm pi of B_w, insert w at each position with interface correction factors. (2) Discovered base = -2*H(W) and correction = 2*(H(W)-H(B_w)) + cycle_derivs. (3) Proved key structural insight: T(I,v) = T(J,v) at s=0 (i,j interchangeable). (4) PROVED alpha_w^H = alpha_w^I as polynomial identity at n=4,5,6,7 using SymPy symbolic computation (diff = 0 exactly). (5) Identified that n>=8 needs additional Delta(alpha_2) terms from THM-013 (VD cycle pairs). (6) Created THM-018 documenting the proof.
**New contributions:**
- THM-018: coefficient identity alpha_w^H = alpha_w^I proved at n<=7 (symbolic)
- 04-computation/q009_alpha_identity.py (insertion decomposition and factor analysis)
- 04-computation/q009_insertion_decomp.py (defect = cycle derivatives verification)
- 04-computation/q009_symbolic_proof.py (symbolic proof n=4,5,6)
- 04-computation/q009_symbolic_n7.py (symbolic proof n=7)
- 04-computation/q009_symbolic_n8.py (n=8 attempt — needs full delta_I formula)
**Unresolved threads:**
- General proof of alpha_w^H = alpha_w^I for all n (key remaining step for OCF)
- Need full delta_I formula (including Delta(alpha_2)) for n>=8 symbolic verification
- Inductive proof approach: use OCF at sub-tournament level to simplify Delta(alpha_k)

## kind-pasteur-2026-03-05-S11 — 2026-03-05 (Tribonacci theorem; tiling class structure deepened)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S10 (context continuation)
**Files read:** All warm-up files, full_class_analysis.py, tiling-class-structure.md, tiling_structure_analysis.py, OPEN-QUESTIONS.md, INVESTIGATION-BACKLOG.md, paper-deep-connections.md, broadcast MSG-011
**Summary of work:**
  (1) DISPROVED the 2^(n-2)+1 conjecture for the full tiling class at n=7,8 (H=31,57 vs expected 33,65).
  (2) DISCOVERED the correct formula: H(T_full_n) = Tribonacci(n) = H(n-1)+H(n-2)+H(n-3), OEIS A000213. Verified through n=12.
  (3) PROVED the Tribonacci recurrence via a combinatorial argument: Ham paths biject to "run decompositions" (ordered sequences of consecutive intervals with gap condition). Conditioning on first interval size yields f(n) = S(n-2) + S(n-3) where S is the partial sum; telescoping gives f(n) = f(n-1)+f(n-2)+f(n-3). Complete proof written up.
  (4) Analyzed self-paired class "perpendicular" geometry: confirmed all self-flip pairs have Hamming midpoint exactly at m/2 (hypercube center). Self-paired classes straddle the transitive<->full diagonal. All self-paired tournaments are self-converse with odd H values.
  (5) Pulled and reviewed opus-S5 findings: H(T_19) computed, alpha_1 conjecture disproved, Omega(T) always perfect, hard-core non-perturbative.
**New contributions:**
- Tribonacci theorem for full tiling class (complete proof in tiling-class-structure.md)
- Updated tiling-class-structure.md with F3 correction, F11 (perpendicular geometry), F12 (Tribonacci)
- Self-paired class H values all odd: 5, 11, 13, 15, 25, 41, 43, 45
**Unresolved threads:**
- Omega perfectness implications for OCF (agent running but incomplete)
- Formula for number of self-paired classes: 1, 2, 8, ...
- Bridge even-odd split to OCF
- Transfer matrix symmetry proof (INV-001, highest priority per opus-S5)

## opus-2026-03-05-S5 — 2026-03-05 (Priority C+E investigation)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4b (continuation after context compaction)
**Files read:** INVESTIGATION-BACKLOG.md, paper-deep-connections.md, hard_core_lattice_gas.py, hard_core_fast.py, ballot_sequence_test.py, compute_H_T19.py, mod4_score_test.py, conflict_graph_catalog.py, SESSION-LOG.md, TANGENTS.md, OPEN-QUESTIONS.md
**Summary of work:** Continued Priority C and E investigation from previous context. Ran all computational scripts created in prior context and gathered results:
  - H(T_19) = 1,172,695,746,915 (INV-024 DONE). H/|Aut| = 6,857,869,865. All 19 endpoints have equal count (Paley symmetry). Computed via C DP in 0.5s.
  - alpha_1 ≡ c_3 (mod 2) conjecture DISPROVED (INV-026). Counterexamples at n=3,4,5. All have alpha_1=0 but c3 odd.
  - Hard-core lattice gas (INV-028): lambda=2 above ALL convergence thresholds. OCF is non-perturbative.
  - Omega(T) structure for n=4,5: max degree grows, independence number = 1 at n=4.
  - mod4_score_test.py has a bug: reports OCF mismatches at n=5, but OCF is correct (bug in cycle-finding code). Alpha_1 result is valid since it's independent of OCF.
  Updated INVESTIGATION-BACKLOG.md with all findings: INV-015, INV-016, INV-019 (paper connections deepened), INV-024 (H(T_19) computed), INV-026 (conjecture disproved), INV-028 (non-perturbative confirmed), INV-029 (bijective proof found).
**New contributions:**
- INVESTIGATION-BACKLOG.md updated with 7 investigation results
- paper-deep-connections.md (created in prior context, verified this session)
- compute_H_T19.py output: H(T_19) = 1,172,695,746,915
- ballot_sequence_test.py: bijective proof of C(L-2, 2k-1) confirmed
- hard_core_lattice_gas.py + hard_core_fast.py: non-perturbative regime confirmed
- mod4_score_test.py: alpha_1 ≡ c_3 (mod 2) disproved (NOTE: has OCF mismatch bug at n=5)
**Unresolved threads:**
- INV-001 (transfer matrix symmetry proof) — highest priority, Feng connection is the best lead
- INV-004 (R-cone + cut-flip proof strategy) — enabled by Rajkumar Prop 1
- INV-027 (conflict graph catalog) — script exists but too slow, needs optimization
- INV-020/021 (Striker/Chapman papers) — correct arXiv IDs not found
- mod4_score_test.py OCF mismatch bug needs fixing
- conflict_graph_catalog.py needs optimization for n>=5

## kind-pasteur-2026-03-05-S10 — 2026-03-05 (Claim B PROVED; gap in proof chain identified)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S9 (context continuation)
**Files read:** All warm-up files, opus-S4c session entry, q009_claim_b.py, q009_direct_proof.py, q009_claim_b_proof.py, signed-adjacency-identity.md, MISTAKES.md
**Summary of work:** (1) Ran q009_claim_b_proof.py — all steps of the inductive proof verified at m=1,...,8 (max error 1.9e-12 at m=8, well within floating-point tolerance). (2) Created THM-016: formal write-up of the inductive proof of Claim (B) for all m. The proof uses h_start expansion, exchange of summation, inductive hypothesis on |W|-1 vertices, and path extension. Tournament condition T(a,b)+T(b,a)=1 is essential (Step 5). (3) IDENTIFIED GAP in opus-S4c proof chain: the step "B(Li,Rj)=B(Lj,Ri) -> OCF" is INVALID per MISTAKE-008. Even-odd split is a consequence of OCF, not equivalent. Corrected signed-adjacency-identity.md. (4) Updated OPEN-Q-009 to reflect that even-odd split is now proved for all n, but OCF remains open.
**New contributions:**
- THM-016: Claim (B) alternating sum identity — PROVED for all m
- q009_claim_b_proof.py verified at m=1,...,8
- q009_creative_hypotheses.py (10 hypotheses tested, dual identity confirmed)
- Corrected signed-adjacency-identity.md (removed false equivalence claim)
- Updated OPEN-Q-009 (even-odd split proved, gap documented)
**Unresolved threads:**
- Bridge the gap: even-odd split -> OCF needs additional identity
- The odd-S sum of Delta(S,R) needs to be connected to the cycle formula
- Investigation backlog items remain (Forcade 1973, Chapman ASM, transfer matrix proof)

## opus-2026-03-05-S4d — 2026-03-05 (THM-016 PROVED + THM-017 PROVED — even-odd split for all n)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4c (context continuation)
**Files read:** All previous session files, OPEN-QUESTIONS, theorem files
**Summary of work:** PROVED THM-016 (Claim B path identity) by induction on m using first-step decomposition. Key steps: (1) separate boundary term, (2) first-step decomposition of h_start, (3) exchange summation order, (4) recognize inner sum as THM-016 at size m-1, (5) use T(v,u)=1-T(u,v) to get H(W')-h_end(W,v), (6) cancellation (-1)^{m-1}+(-1)^m=0. Then PROVED THM-017 (B(Li,Rj)=B(Lj,Ri)) for all n by induction on |W|, using THM-016 for the d-independent part and the induction hypothesis for the d-dependent part. Both proofs verified numerically at all sizes. CORRECTED: the even-odd split is necessary but NOT sufficient for OCF (kind-pasteur-S8 was right). Investigated the gap: F(x)=(1+x)Q(x) where F(1)=delta_H and F(-1)=0, but Q doesn't have universal structure.
**New contributions:**
- 04-computation/q009_claim_b_proof.py (THM-016 proof verification)
- 04-computation/q009_full_proof.py (THM-017 full proof verification)
- 04-computation/q009_bridge_gap.py (investigating gap to OCF)
- 01-canon/theorems/THM-016-hamiltonian-alternating-sum.md (PROVED)
- 01-canon/theorems/THM-017-even-odd-split-proof.md (PROVED)
- Updated OPEN-QUESTIONS, signed-adjacency-identity.md with honest assessment
**Unresolved threads:**
- OCF for all n still OPEN — even-odd split insufficient
- Need delta_H = delta_I for all n (currently proved n<=8)
- The Q polynomial in F(x)=(1+x)Q(x) may yield a path forward

## opus-2026-03-05-S4c — 2026-03-05 (Claim B discovery — key reduction for OCF proof)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4
**Files read:** All warm-up files, previous session's exploration scripts
**Summary of work:** Major progress on OPEN-Q-009. (1) Corrected the s-evenness reduction (C_w+D_w=0 was wrong; correct statement is s-linear coefficients vanish). (2) Verified B(Li,Rj)=B(Lj,Ri) as polynomial identity at n=3,...,7 with real-valued arcs in [-2,3]. (3) Discovered the PAIRING STRUCTURE: for each u!=v, left(u)=right(u) exactly, meaning the d-dependent part of alpha_v vanishes by the INDUCTIVE HYPOTHESIS. (4) Isolated CLAIM (B): a standalone identity about tournament Hamiltonian path weights: sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S,v) = (-1)^{m+1} h_end(W,v). Verified at m=1,...,8. Proved algebraically at m=1,2,3. (5) Established full proof chain: Claim B -> B(Li,Rj)=B(Lj,Ri) -> OCF -> Claim A.
**New contributions:**
- 04-computation/q009_cw_dw_proof.py (showed C_w+D_w!=0, corrected reduction)
- 04-computation/q009_algebra_n4.py (exact n=4 algebra, s-coordinate formula)
- 04-computation/q009_s_quadratic.py (quadratic s-structure analysis)
- 04-computation/q009_direct_proof.py (alpha_v decomposition, pairing discovery)
- 04-computation/q009_claim_b.py (Claim B verification m=1,...,8)
- Updated signed-adjacency-identity.md with proof structure and Claim B
**Unresolved threads:**
- Prove Claim (B) for all m — the ONLY remaining step for OCF proof
- Claim (B) is an alternating-sum identity about internal tournament Hamiltonian paths
- Proved at m=1,2,3 by direct algebra; need general proof (possibly by induction on m)

## kind-pasteur-2026-03-05-S9 — 2026-03-05 (tex deep analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S8 (context continuation)
**Files read:** Full parity_tournaments_fixed.tex (2189 lines), all opus-S4b new files (paper-connections.md, INVESTIGATION-BACKLOG.md, 4 new computation scripts), inbox messages MSG-004 and MSG-005
**Summary of work:** Pulled opus-S4b contributions (transfer matrix symmetry discovery, investigation backlog, paper connections). Then performed line-by-line analysis of entire tex file. Found 5 issues: (1) DR mod-4 proof (Thm 7.4) has broken arithmetic. (2) SE-SYT formula (Thm 7.3) produces non-integer. (3) Transitive uniqueness proof incomplete. (4) Verification record outdated. (5) Rajkumar et al. missing from bibliography.
**New contributions:**
- 03-artifacts/drafts/tex-deep-analysis.md (comprehensive analysis report)
- INV-028b, INV-029b, INV-030b added to backlog
**Unresolved threads:**
- Forcade 1973 GF approach (INV-023)
- Chapman 2001 ASM bijection (INV-021)
- Transfer matrix symmetry proof (INV-001)

---

## opus-2026-03-05-S4b — 2026-03-05 (repo scour + backlog creation)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4b (context continuation — second part)
**Files read:** Full tex paper (open problems, bibliography, literature), TANGENTS, OPEN-QUESTIONS, MISTAKES, MASTER_FINDINGS, proof-landscape, signed-adjacency-identity, all code file headers, paper-connections.md
**Summary of work:** Completed the user's directive to establish a systematic practice of scouring the repo for leads. Read every significant file in the repo. Created INVESTIGATION-BACKLOG.md with 31 cataloged leads across 5 priority tiers (A: blocks OCF proof, B: structural understanding, C: references to investigate, D: computational targets, E: tangent investigations). Updated CLAUDE.md with mandatory Step 5b (repo scouring every session). Notable findings: Rajkumar paper (2110.05188) is NOT in the tex bibliography despite being connected; 6 bibliography entries have NOT been investigated for OCF connections (Striker, Chapman, Eplett, El Sahili x2, Forcade GF approach); 3 open problems from the paper have zero investigation (mixed graphs, Striker-Chapman equivariance, realizable conflict graphs).
**New contributions:**
- `00-navigation/INVESTIGATION-BACKLOG.md` — comprehensive lead catalog (31 items)
- `CLAUDE.md` Step 5b — mandatory repo-scouring practice for all agents
- Memory file created at `.claude/projects/.../memory/MEMORY.md`
**Unresolved threads:**
- INV-001 (transfer matrix symmetry proof) — highest priority
- INV-008 (Striker-Chapman equivariance) — completely uninvestigated
- INV-015 (add Rajkumar to bibliography)
- INV-024 (compute H(T_19)) — feasible, not done
- n=8 C verifier may still be running from earlier in this session

## kind-pasteur-2026-03-05-S8 — 2026-03-05 (error audit)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S7 (context continuation)
**Files read:** All canon files, all opus-S3/S4 new files, all proof scripts, TANGENTS, OPEN-QUESTIONS
**Summary of work:** Pulled opus-S3 and opus-S4 contributions, then audited the entire codebase for errors. Found and fixed 4 issues: (1) BUG in sympy_proof_n8.py — used simplified n<=7 formula which THM-013 says FAILS at n=8; rewrote to use full A-clique formula. (2) LOGICAL ERROR in even-odd split lemma — claimed "equivalent to OCF" but it's only a consequence; the odd-S sum of Delta(S,R) differs from the cycle formula [g-l]*H(R) at the per-subset level. (3) Stale claim in proof-landscape saying n<=7 instead of n<=8. (4) Duplicate tangent number T040. Verified all SymPy proofs (n=4,5,6) and opus n=7 proof pass correctly. NOTE: opus-S4b ran in parallel and found additional structure (Signed Position Identity, tournament-specificity, bracket analysis) — see their session entry.
**New contributions:**
- MISTAKE-008: even-odd split equivalence claim corrected
- MISTAKE-009: sympy_proof_n8.py formula bug documented and fixed
- sympy_proof_n8.py rewritten with correct full A-clique formula
- Tangent numbering fixed (T040 duplicate -> T044-T047)
- Even-odd-split-lemma.md corrected, OPEN-Q-009 corrected, TANGENTS corrected
**Unresolved threads:**
- Run sympy_proof_n8.py overnight (corrected version, full A-clique formula)
- Prove OCF for all n — central open problem
- Reconcile even-odd split equivalence claim with opus-S4b's Signed Position Identity analysis

---

## opus-2026-03-05-S4b — 2026-03-05 (Transfer matrix symmetry + paper connections)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S3 (parallel session to S4, then extended)
**Files read:** All warm-up files, inbox snippet from S4 diff, all n=8 computation files, arXiv:2110.05188, arXiv:2510.25202
**Summary of work:** Three phases: (1) Verified Even-Odd Split Lemma, proved tournament-specific/polynomial, built C verifier for n=8. (2) Analyzed two external papers: "Tournament Representations" (flip classes, locally transitive, R-cones) and "Dual Burnside Process" (Q=AB factorization, spectral correspondence). (3) MAJOR DISCOVERY: the transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is ALWAYS SYMMETRIC (verified 7500+ tests, n=4..8). This is STRONGER than Even-Odd Split — each cross-term individually agrees. Connected to Burnside detailed balance. Also found locally transitive tournaments DO have 5-cycles and 7-cycles (conjecture "only 3-cycles" fails). All OCF tests pass for R-cones, automorphism-symmetric tournaments, locally transitive tournaments.
**New contributions:**
- 03-artifacts/drafts/paper-connections.md (10 connections to external papers)
- 04-computation/locally_transitive_test.py (cycle structure of rank-2 tournaments)
- 04-computation/burnside_connection_test.py (character/Fourier decomposition)
- 04-computation/transfer_matrix_test.py (transfer matrix structure)
- 04-computation/symmetry_check.py (VERIFIES transfer matrix symmetry, 7500+ tests)
- T045 (Transfer matrix symmetry discovery)
- T046 (Tournament representations / flip class connection)
**Unresolved threads:**
- PROVE transfer matrix symmetry for all n (would prove OCF)
- Explore detailed-balance / reversibility interpretation
- R-cone simplification strategy (prove OCF for R-cones, extend via flip)

---

## opus-2026-03-05-S4 — 2026-03-05 (n=7 and n=8 PROVED, even-odd split discovery)
**Account:** Eliott (opus machine)
**Continuation of:** opus-2026-03-05-S3 (context continuation)
**Files read:** TANGENTS.md (resolved merge conflict), THM-015, PROP-001, OPEN-QUESTIONS.md, SESSION-LOG.md, symbolic_proof.py, symbolic_proof_fast.py, kind-pasteur inbox messages
**Summary of work:** Resolved git merge conflict from S3 rebase (renumbered my T032 to T039). Then wrote optimized n=7 exhaustive verifier using numpy-vectorized bitmask approach — PROVED OCF at n=7 (2^20 = 1,048,576 configs in 4 seconds). Extended to n=8 with chunked processing — PROVED OCF at n=8 (2^27 = 134,217,728 configs, 57 minutes). Discovered the Even-Odd Split Lemma: the adj decomposition delta = sum_S Delta(S,R) splits equally between even-|S| and odd-|S| terms, so delta = 2*(odd-S sum). This connects directly to the cycle formula (only odd cycles). The alternating sum sum(-1)^|S| Delta(S,R) = 0 is equivalent to OCF but provides a clean algebraic reformulation. Verified n=5,...,8.
**New contributions:**
- THM-015 updated: n=7 and n=8 PROVED
- 04-computation/q009_prove_n7.py (numpy-vectorized exhaustive verifier)
- 04-computation/q009_prove_n8.py (chunked n=8 verifier)
- 04-computation/q009_even_odd_split.py (even-odd split discovery)
- 04-computation/q009_alternating_sum.py (alternating sum analysis)
- 03-artifacts/drafts/even-odd-split-lemma.md (new proof angle documentation)
- Tangent T040 (even-odd split)
- OPEN-Q-009 updated with proof frontier and even-odd split
**Unresolved threads:**
- n=8 exhaustive verification COMPLETE (134M/134M, 57 min)
- Prove the alternating sum identity for ALL n (equivalent to OCF)
- The even-odd split is a clean algebraic reformulation but doesn't simplify the proof

---

## kind-pasteur-2026-03-05-S7 — 2026-03-05 (n=7 proof + structural analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S6 (context continuation)
**Files read:** SESSION-LOG.md, OPEN-QUESTIONS.md, TANGENTS.md, MISTAKES.md, definitions.md, THM-013, THM-014, THM-015, all S6 code files
**Summary of work:** Extended OCF proof to n=7 via two independent methods: (1) SymPy polynomial identity verification (20 arc variables, 1858 monomials, 77s) and (2) exhaustive {0,1} enumeration (1,048,576 cases, 775s). Both confirm delta_H = delta_I as polynomial identity at n=7. Also explored structural decomposition for general proof: proved delta_H = -sum s_x*(R1+R2) algebraically, showed R1+R2 = 2*H(B_x) at n=4, found excess terms at n=5 satisfy sum s_x*exc(x) = -2*(D5-C5) on {0,1} but NOT as polynomial identity (non-multilinear). Investigated f(S)+f(S^c) pairing: universal at n=4, fails at n>=5 due to 5-cycle corrections. Created complete hand proof for n=4 via f(S) decomposition.
**New contributions:**
- THM-015 updated to PROVED at n<=7 (was n<=6)
- 03-artifacts/code/sympy_proof_n5.py (SymPy polynomial identity proof at n=5)
- 03-artifacts/code/sympy_proof_n6.py (SymPy polynomial identity proof at n=6)
- 03-artifacts/code/sympy_proof_n7.py (SymPy polynomial identity proof at n=7)
- 03-artifacts/code/algebraic_proof_n4.py (complete hand proof at n=4)
- 03-artifacts/code/pairing_proof.py (f(S)+f(S^c) analysis, shows n>=5 failure)
- 03-artifacts/code/proof_structure_analysis.py (pred/succ decomposition, s-signature analysis)
**Unresolved threads:**
- Prove polynomial identity for ALL n — central remaining open problem
- n=8 SymPy verification feasible (user offered overnight compute)
- Structural proof approaches: transfer matrix, induction on n, permanent expansion
- Non-multilinear excess identity blocks naive decomposition approach

---

## opus-2026-03-05-S3 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S2 (context overflow, resumed)
**Files read:** All warm-up + THM-012, THM-013, TANGENTS, OPEN-QUESTIONS
**Summary of work:** Extended tournament_lib.py with arc-flip functions. Independently verified THM-013 adjacency identity. Extended OCF verification to n=10 (first time n>=8 verified).
1. Added adj_count(), flip_arc(), verify_thm013(), verify_ocf(), independence_poly_at_fast() to tournament_lib.py
2. THM-013 adjacency identity: VERIFIED n=5 exhaustive (10240/10240), n=6 sampled (3000/3000)
3. OCF H(T)=I(Omega(T),2): PROVED n<=7 exhaustive (1,048,576 arc assignments at n=7, 27min), n=8 (500 random), n=9 (100), n=10 (30) -- all 0 failures
4. Algebraic analysis: decomposed adj(i,j)-adj'(j,i) as subset convolution. Per-vertex LHS/RHS decomposition does NOT match -- proof must work globally.
5. Key insight: at n<=8, max independent set size in Omega(T) is 2; at n=9, it's 3. Fast I computation exploits this.
6. Explored proof strategies: per-vertex nadj decomposition (dead end), transfer matrix/permanent (no clean formula found), bijection search (suggestive structure but no proof), inclusion-exclusion over blockers (correct framework but doesn't simplify).
7. Key structural finding: s=0 vertices never block swaps (cleanly separates U_T from U_T'), but nadj_sum/H(B_x) ratio is not constant -- the identity is fundamentally global.
**New contributions:**
- tournament_lib.py: 6 new functions (adj_count, flip_arc, compute_s_x, count_directed_5_cycles_through_arc, verify_thm013, verify_ocf, independence_poly_at_fast)
- verify_ocf_sweep.py: clean verification script
- symbolic_proof_n7.py: exhaustive n=7 OCF proof (1,048,576/1,048,576 PASS)
- OCF PROVED through n=7, verified through n=10
- Updated THM-013, THM-015, OPEN-QUESTIONS.md with n=7 exhaustive results
- 4 exploratory scripts: unmatched_decomposition_by_vertex.py, inductive_structure.py, transfer_matrix_approach.py, bijection_search.py, proof_strategy_analysis.py
**Unresolved threads:**
- PROVE the identity for all n (not just n<=7): the central open problem
- Subset convolution framework is the right language but per-vertex doesn't match; need global argument
- The "2-colored independent cycle set" interpretation of I(Omega,2) may lead to a bijective proof
- n=8 proof feasible with optimization (2^27 ~134M cases) but would take hours

---

## kind-pasteur-2026-03-05-S6 — 2026-03-05 (polynomial identity proof)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S5 (context continuation)
**Files read:** SESSION-LOG.md, OPEN-QUESTIONS.md, TANGENTS.md, THM-013, THM-014, all S5 code files
**Summary of work:** Resolved 6-file git merge conflict from S5 rebase. Then discovered and proved a NEW PROOF METHOD for OCF: the swap involution's unmatched path count U_T'-U_T can be expressed as a polynomial identity in arc variables, and this polynomial equals delta_I from THM-013. Proved by hand at n=4 (clean formula: U_T'-U_T = 4-2(p+q+r+t) = 2*sum(s_x)). Verified exhaustively at n=5 (512 cases) and n=6 (16384 cases). This constitutes a PROOF of OCF for n<=6. Found and fixed a sign convention error in THM-013's D/C labels.
**New contributions:**
- THM-015 in 01-canon/theorems/ (swap polynomial identity, proves OCF at n<=6)
- 03-artifacts/code/symbolic_proof.py (polynomial identity verifier for n=4,5,6)
- 03-artifacts/code/symbolic_proof_fast.py (optimized n=6 version)
- 03-artifacts/code/unmatched_decomposition.py (structural analysis of blocking)
- 03-artifacts/code/unmatched_vs_formula.py (aggregate identity tests)
- Tangents T037-T038 (polynomial identity proof, sign convention warning)
- Resolved S5 merge conflicts (TANGENTS renumbered T028-T036, all files clean)
**Unresolved threads:**
- Prove the polynomial identity for ALL n (not just n<=6) — this would prove OCF/Claim A
- The n=4 hand proof suggests a transfer matrix / generating function approach
- n=7 verification feasible (2^19 = 524K cases) but needs optimization
- OPEN-Q-009 partially resolved: proof method identified, works at n<=6

---

## kind-pasteur-2026-03-05-S5 — 2026-03-05 (arc-flip session)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S4
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, CONJ-001, THM-003 (Claim B proof), LaTeX paper claim_strategies section
**Summary of work:** Attacked OPEN-Q-009 (arc-reversal invariance) with new angle. Discovered that E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips (exhaustive n<=5, random n=6). Derived algebraic formula for delta_I via A-clique argument adapted to arc flips. Showed the simple formula delta=2*(gained-lost) works at n<=5 but fails at n>=6. Tested and ruled out contiguous-block decomposition. Formulated PROP-001 (arc-flip identity) as a new proof strategy for OCF/Claim A equivalent to proving a combinatorial identity about Ham path changes under arc flips weighted by complement H-counts.
**New contributions:**
- PROP-001 in 01-canon/theorems/ (arc-flip identity, new proof strategy)
- 03-artifacts/code/arc_reversal_study.py (delta_H decomposition under arc flips)
- 03-artifacts/code/ocf_arc_flip_study.py (E(T) invariance, joint delta_H/delta_I distribution)
- 03-artifacts/code/delta_I_correct_formula.py (algebraic delta_I formula verification)
- 03-artifacts/code/delta_I_fast_test.py, paths_via_cycles_test.py, contiguous_block_test.py
- OPEN-Q-009 reformulated with cleaner E(T) version
- Tangents T028-T031 (arc-flip approach, dead ends)
**Unresolved threads:**
- PROP-001: prove the arc-flip Ham path identity combinatorially
- OPEN-Q-009: the E(T) formulation is cleaner but still needs proof
- Three possible proof approaches: transfer matrix, generating function, involution

---

## kind-pasteur-2026-03-05-S3 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S2
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, DISC-001, CONJ-001, theorem files
**Summary of work:** Rebased kind-pasteur onto main (resolved 6-file conflict). Closed DISC-001 (mu bug does not contaminate verification). Investigated OPEN-Q-010/011 at n=7: wrote test_n7_ABD.py computing A (TypeII sum), B (mu-weighted), D (cycle mu sum) for 1050 pairs. Finding: A=/=D at n=7 in general (only 5.9% exact match). Near-cancellation is statistical (mean A-D~0.1) not algebraic.
**New contributions:**
- 03-artifacts/code/test_n7_ABD.py (correct A-B-D computation at n=7)
- 03-artifacts/code/test_perpath_n7.py (initial wrong test -- documents why naive approach fails)
- DISC-001 moved to resolved/; OPEN-Q-010/011 updated with n=7 findings
- TANGENT T027 (n=7 A-B-D near-cancellation is statistical not algebraic)
**Unresolved threads:**
- OPEN-Q-009: arc-reversal invariance -- the key step, still untouched
- OPEN-Q-012: tower hypothesis -- not yet investigated
- OPEN-Q-013: H(T_p) formula for Paley primes

---

## kind-pasteur-2026-03-05-S2 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S1
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, DISC-001, CONJ-001, CONJ-002, LEM-001, LEM-002, THM-002
**Summary of work:** Resolved OPEN-Q-009 (kind-pasteur numbering) by direct computation. h_QR=h_NQR=201, c_9(T_11)=11055. Direct enumeration gives H(T_11)=95095=55*1729, refuting CONJ-002 for p=11. Added finish_session.py enforcement tooling.
**New contributions:**
- H(T_11)=95095 computed directly (03-artifacts/code/compute_H_T11.py)
- CONJ-002 REFUTED for p=11; OPEN-Q-013 opened (correct formula for H(T_p))
- MISTAKE-006: ratio-coincidence c_k/C(11,k) has no basis
- TANGENTS T022-T026 (Paley structure, 1729, symmetry)
- agents/finish_session.py and agents/check_session_closed.py (end-of-session enforcement)
**Unresolved threads:**
- OPEN-Q-013: What is the correct formula for H(T_p)?
- OPEN-Q-009 (opus): arc-reversal invariance — the key proof step
- DISC-001: needs formal close

---

## kind-pasteur-2026-03-05-S1 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** fresh start — first session on this machine
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, inbox/PROCESSING-REPORT.md, agents/REGISTRY.md, all theorem files, DISC-001
**Summary of work:** Registered as new agent. Processed PALEY_T11_c9_ANALYSIS.md inbox doc: extracted LEM-001, LEM-002, CONJ-002. Argued Position A in DISC-001.
**New contributions:** LEM-001, LEM-002, CONJ-002, MISTAKE-006 (ratio-coincidence), DISC-001 Letter 2
**Unresolved threads:** DISC-001 response needed; h_QR/h_NQR computation

---

## kind-pasteur-2026-03-05-S4 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S3
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, opus broadcast messages
**Summary of work:** Structural fix session. kind-pasteur was pushing to claude/kind-pasteur branch (invisible to opus and CI). Merged claude/kind-pasteur into main. Fixed Windows encoding issues in processor.py (CP1252 fails on Unicode). Updated finish_session.py to push HEAD:main always.
**New contributions:**
- All kind-pasteur S1-S3 research now on main (Paley computation, n=7 A-B-D, DISC-001 resolved)
- agents/processor.py: Windows UTF-8 fix (stdout.reconfigure + write_text encoding)
- agents/finish_session.py: push to origin HEAD:main (not current branch)
**Unresolved threads:**
- OPEN-Q-009 (arc-reversal invariance): top priority, untouched
- OPEN-Q-012 (tower hypothesis): untouched

---

## opus-2026-03-05-S2 — 2026-03-05
**Account:** Claude Opus 4.6 (e's MacBook)
**Continuation of:** opus-2026-03-05-S1 (parallel session)
**Files read:** All navigation, canon, court files; LaTeX paper; FINAL_FINDINGS.md; file.txt; tournament_lib.py
**Summary of work:** Proved THM-008b (general mu triviality), THM-012 (partial mu invariance under arc flips). Disproved MISTAKE-004 (OCF IS a valid closed form — verified H(T)=I(Omega(T),2) for all n<=6). Resolved OPEN-Q-011 (near-cancellation is statistical, not structural). Verified adjacency formula H(T)-H(T')=adj(i,j)-adj'(j,i). Established arc-reversal decomposition framework for Q-009.
**New contributions:**
- THM-008b: general mu triviality bound L >= n-2
- THM-012: mu invariant under arc flips when at least one endpoint in V(C)\{v}
- DISC-002: MISTAKE-004 is wrong — OCF is a valid closed form (verified 33,864 tournaments)
- OPEN-Q-011 resolved (near-cancellation is statistical artifact)
- Adjacency formula verified; arc-reversal decomposition framework built
- 6 computation scripts in 04-computation/
**Unresolved threads:**
- DISC-002 needs formal resolution (retract MISTAKE-004)
- Arc-reversal proof strategy (Q-009) partially developed — key obstacle is tracking cycle creation/destruction and mu changes simultaneously
- Claim A proof for general n remains open (OPEN-Q-002)

---

## opus-2026-03-05-S1 --- 2026-03-05
**Account:** opus (Claude Opus 4.6)
**Continuation of:** SYSTEM-2026-03-05-S1 (initial setup)
**Files read:** All navigation files, all canon files, DISC-001, agents/processor.py, inbox/processor.py, README.md
**Summary of work:** Built the complete tournament computation library (tournament_lib.py) implementing all core objects from definitions.md. Independently verified Claim A at n<=6 (196,608 pairs, 0 failures) with a MISTAKE-001-compliant implementation. Also verified Claim B at n<=5, Redei at n<=6, and Claim A at n=7 by random sampling (100 tournaments, 0 failures). Ingested two contributed files (FINAL_FINDINGS.md, file.txt) containing substantial new results. Extracted 4 new theorems (THM-008 through THM-011), 2 new mistakes (MISTAKE-004, MISTAKE-005), resolved OPEN-Q-001 and OPEN-Q-003, added 4 new open questions (Q-009 through Q-012), and documented 3 dead ends (T016-T018).
**New contributions:**
- 03-artifacts/code/tournament_lib.py (CODE-000): complete computation library
- 03-artifacts/code/verify.py (CODE-000v): CLI verification runner
- THM-008: mu=1 trivially for n<=5 (resolves OPEN-Q-001)
- THM-009: per-path failure characterization at n=6 (resolves OPEN-Q-003)
- THM-010: n=4 block-counting theorem (exactly 3 Ham paths)
- THM-011: general block-counting formula H_C^+(T)
- MISTAKE-004: OCF is recursive, not closed-form over all cycles
- MISTAKE-005: cycle bijection under arc reversal fails
- OPEN-Q-009: arc-reversal invariance (key unproved step)
- OPEN-Q-010 through Q-012: new questions from ingested files
- Dead ends T016-T018 documented
- Registered agent: opus
**Unresolved threads:**
- Claim A proof for general n (CONJ-001) -- arc-reversal invariance (OPEN-Q-009) is the key step
- DISC-001 can be moved to resolved
- The near-cancellation phenomenon (OPEN-Q-011) is a promising proof strategy lead
- Per-path formula with 3+5 cycles at n=7 (OPEN-Q-010) needs computational testing
- Computation library could be extended: inshat, per-path identity, arc-reversal testing

---

## SYSTEM-2026-03-05-S1 — 2026-03-05
**Account:** System (initial setup)
**Continuation of:** fresh start — first Cowork session
**Files read:** MASTER_FINDINGS.md (uploaded), parity_tournaments_fixed.tex (uploaded)
**Summary of work:** Built the full research directory system. Processed two uploaded source files as the first contributions. Extracted theorems F1–F5 and Claim A to canon. Logged MISTAKE-001. Opened DISC-001 as a potential court case seed (μ computation bug vs. paper's 0-failure verification claim).
**New contributions:**
- THM-001 through THM-007, CONJ-001 in 01-canon/theorems/
- MISTAKES.md: MISTAKE-001 (μ computation bug in scripts 6-9)
- 00-navigation/TANGENTS.md, OPEN-QUESTIONS.md populated from source files
- 03-artifacts/drafts/parity_tournaments_fixed.tex archived
- 02-court/active/DISC-001-mu-bug-vs-verification.md opened
**Unresolved threads:**
- Claim A proof (the central open problem)
- Why does per-path identity hold at n=5 despite 5-cycles? (OPEN-Q-001)
- μ computation bug in scripts 6-9 needs investigation (DISC-001)
- Formal proof of Fix(σ) = 2^{m²} for self-evacuating SYT (OPEN-Q-007)
