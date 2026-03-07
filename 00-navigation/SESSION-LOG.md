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

## opus-2026-03-07-S32 — 2026-03-07 (DEFORMED EULERIAN NUMBERS — complete a_k(T) formula)
**Account:** opus
**Continuation of:** opus-2026-03-07-S31 (context limit)
**Summary of work:**
  (1) **THM-062: Deformed Eulerian Numbers (PROVED).** Complete closed-form:
      a_k(T) = A(n,k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T)
      where c_k^{(f,d)} = sum_j A(f+1,j) * C(d-f, k-j) * (-1)^{d-f-k+j}
  (2) **Key insight:** Transitive tournament gives standard Eulerian numbers A(n,k).
      All invariant corrections sum to zero (preserving n! total).
  (3) **Bivariate formula:** Phi_T(x,y) = A_n(x,y) + sum_I 2^parts * A_{f+1}(x,y) * (x-y)^{n-1-f} * I(T).
      Reduces to W(r) on line x-y=1.
  (4) **Derivative formulas:**
      - F'_f(1/2) = 2^{f+1} - 2 (Mersenne-like)
      - a_{n-2}(T) = OCF polynomial. At n=7: 120 + 48*t3 - 12*t7
      - F''_f(1/2) = f(f-1) + 2(f-1)(2^{f+1}-f-2) + 2*A(f+1,2)
      - a_{n-3}(T) at n=7: 1191 + 30*t3 - 18*t5 + 30*t7 - 36*bc (constant = A(7,2))
  (5) **Negative results:** Tournament Eulerian poly E_T(x) NOT real-rooted (74-100% complex).
      W(r) roots also generally off real/imaginary axes.
  (6) **Verification:** 260+ coefficient checks (n=5: 75, n=7: 140, n=9: 45). All exact.
**New contributions:** THM-062-forward-edge-distribution.md, 7 computation scripts
**Unresolved threads:** General F^{(k)}_f(1/2) closed form; a_k at n=9 with full invariants; combinatorial interpretation of c_k^{(f,d)}; positivity constraints from a_k >= 0

## kind-pasteur-2026-03-07-S26 — 2026-03-07 (THM-061 Anti-Evaluation + W(r) Flip Analysis)
**Account:** kind-pasteur
**Continuation of:** kind-pasteur-2026-03-06-S25h (context limit)
**Summary of work:**
  (1) **THM-061 (PROVED): W(-1/2) = (-1)^{n-1} * H(T) for all tournaments.** Elementary proof: W(-1/2) counts backward Ham paths with sign (-1)^{n-1}; reversing backward paths bijects with forward paths. Verified exhaustive n<=6, sampled n=7.
  (2) **F_f(-1/2) = (-1)^f PROVED** via Eulerian formula: at r=-1/2, only k=f term survives since (r+1/2)=0, giving A(f+1,f)*(-1)^f = (-1)^f.
  (3) **W(0) vanishing at n=5:** W(0) = 1 - t3 + 2*t5 = 0 for ALL n=5 tournaments with odd t3 (100%). This follows from the constraint: odd t3 forces t5 = (t3-1)/2, hence H = 3*t3 when t3 is odd.
  (4) **Palindromic Eulerian distribution:** a_k(T) = a_{n-1-k}(T) for ALL tournaments (forward edge distribution is symmetric). Proof: reversal bijection. This is EQUIVALENT to r-parity of W(r) (THM-059 property iii).
  (5) **W(r) flip analysis:** W(T) - W(flip(T)) = C_t3*dt3 + C_t5*dt5 (verified all 8 GS pairs at n=5). W(-r) ≠ W_flip(r) in general.
  (6) **Skeleton spectral + Eulerian:** Silver ratio (1+sqrt2) eigenspaces are t3-independent. H projects strongly onto the 1+sqrt2 eigenspace. Degree = class size for skeleton adjacency.
  (7) **Total sum formula:** sum_T W_T(0) = 2^{m-n+2} at odd n, 0 at even n. Only backbone-only permutations survive averaging.
**New contributions:** THM-061-anti-evaluation.md, 8 computation scripts in 04-computation/
**Unresolved threads:** W(0) vanishing pattern at n>5; algebraic proof of odd-t3 => t5=(t3-1)/2 at n=5; skeleton spectral structure at n=7; free position universality (opus's open question)

## opus-2026-03-06-S11b (continued^8) — 2026-03-07 (OCF PROVED at n=7 via Fourier decomposition)
**Account:** opus
**Continuation of:** opus-2026-03-06-S11b (continued^7, context limit)
**Summary of work:**
  (1) **DEGREE-4 IDENTITY PROVED**: The last remaining identity at n=7. Discovered that degree-4 Fourier space is exactly 2-dimensional: Type P (5-vertex spanning paths) and Type Q (6-vertex disjoint P₂ pairs), with DISJOINT supports.
  (2) **EXACT RATIONAL CONSTANTS**: [deg-4 of t₇] = (1/2)·[deg-4 of t₅] + [deg-4 of α₂]. And w₂/4 = 3·[deg-4 of t₅] + 6·[deg-4 of α₂]. Both verified analytically over all 5985 four-edge monomials with zero error.
  (3) **COUNTING LEMMAS**: For type P: c₅=2, c₇=4, paths=12. For type Q: c₇=8, paths=24, α₂=1. All match clean combinatorial arguments (vertex insertion).
  (4) **OCF AT n=7 COMPLETE**: All 4 degree-homogeneous identities proved: degree 0 (trivial), degree 2 (proportionality), degree 4 (counting lemmas), degree 6 (path-cycle bijection).
  (5) **n=9 PRELIMINARY**: 9-cycles introduce new graph types (7v, 8v). Degree-2 c_{a2}=5 confirmed. t₅_d4 and α₂_d4 remain nearly uncorrelated (corr ≈ -0.05), suggesting 2D structure persists.
**New contributions:** `degree4_identity_n7.py`, `degree4_proof_n7.py`, INV-050 updated
**Unresolved threads:** Extend to n=9 (degrees 4 and 6). General-n proof of middle-degree identities. Connection to master polynomial / Eulerian numbers.

## opus-2026-03-07-S31 — 2026-03-07 (Master Polynomial PROVED: Eulerian Numbers + EGF)
**Account:** opus
**Continuation of:** opus-2026-03-06-S30 (context limit)
**Summary of work:**
  (1) **RECURRENCE FIX**: THM-059 had j^2 instead of (j+1)^2 in the b-triangle recurrence. Fixed and logged as MISTAKE-016.
  (2) **F_7 COMPUTED**: F_7(r) = -124r + 3024r^3 - 20160r^5 + 40320r^7. Matches n=8 background polynomial.
  (3) **n=10 VERIFICATION**: All 25 W-coefficients verified with zero error. New invariant bc55 (disjoint 5-cycle pairs) enters. F_9 confirmed.
  (4) **PERMUTATION FORMULA PROVED**: F_f(r) = sum_{sigma in S_{f+1}} prod (r + sign(sigma(i+1)-sigma(i))/2). Verified f=0..6.
  (5) **EULERIAN NUMBER DECOMPOSITION**: F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k. This identifies the master polynomial with the Eulerian polynomial evaluated at the "ascent/descent weighting."
  (6) **EXPONENTIAL GENERATING FUNCTION PROVED**: sum F_j(r) x^j/(j+1)! = (e^x-1)/(x((r+1/2)-(r-1/2)e^x)). Derived from classical Eulerian polynomial EGF.
  (7) **TANGENT NUMBERS PROVED**: At r=0, the EGF becomes (2/x)tanh(x/2), giving F_{2k}(0) = (-1)^k T_{k+1}/4^k. Previously verified; now PROVED algebraically.
  (8) **CENTRAL FACTORIAL CONNECTION EXPLAINED**: The b_{k,j} numbers are exactly the Eulerian polynomial A_{2k+1}(t) expanded in the u=pq basis. The (j+1)^2 factor in the recurrence comes from this change of basis.
**New contributions:** THM-059 upgraded to PROVED, MISTAKE-016, master_poly_n10_verify.py
**Unresolved threads:** Algebraic proof of FREE POSITION UNIVERSALITY (why F_f is independent of which positions are free); deformed independence polynomial properties

## kind-pasteur-2026-03-06-S25h — 2026-03-06/07 (Bipartite Skeleton + t3 Parity)
**Account:** kind-pasteur
**Continuation of:** kind-pasteur-2026-03-06-S25g (context limit)
**Summary of work:**
  (1) **THM-060 (PROVED): Blue line skeleton is bipartite at odd n, with t3 parity determining the bipartition.** For GS tilings, t3(T)+t3(flip(T)) = n-2 + even = ODD at odd n. Verified n=3,5,7,9.
  (2) **Proof mechanism:** Consecutive triples contribute exactly 1 each (total n-2, odd). Non-consecutive triples contribute even total (Types A and B). GS constraint is essential for Type B evenness.
  (3) **Even/odd dichotomy:** At even n, t3 sum is EVEN, skeleton is NOT bipartite (has 3-cycles at n=6), self-flips occur. At odd n: 100% cross-class, bipartite, no self-flips.
  (4) **NSC sidedness:** NSC pair members always have same t3 (since T and T^op have same t3), so both sit on same side of bipartition.
  (5) **Diff distribution = scaled Binomial(n-2, k)**: The n-2 consecutive triples act as independent ±1 coins.
  (6) **Spectral analysis (n=5):** Skeleton eigenvalues {±(1+√2), ±1, ±1, ±(√2-1)} — silver ratio appears! K^2 diagonal = GS class sizes. Antiferromagnetic interpretation.
  (7) **GS cube geometry:** Flip = antipodal map on {0,1}^k. Cube adjacency (local, Hamming-1) ≠ skeleton (global, antipodal flip). Two fundamentally different graph structures.
**New contributions:** THM-060-bipartite-skeleton.md, bipartite-skeleton-synthesis-S25h.md, 12 computation scripts in 04-computation/
**Unresolved threads:** Algebraic proof of Type B evenness; spectral structure at n=7; connection between skeleton eigenvalues and tournament invariants; Möbius strip interpretation at even n

## opus-2026-03-06-S30 — 2026-03-06 (Universal Master Polynomial + Central Factorial Numbers)
**Account:** opus
**Continuation of:** opus-2026-03-06-S29 (context limit)
**Summary of work:**
  (1) **UNIVERSAL MASTER POLYNOMIAL THEOREM (THM-059)**: The per-invariant r-polynomial C_I(r,n) = 2^{parts(I)} * F_f(r) where f = free position count. The F_j are determined by the CENTRAL FACTORIAL NUMBER TRIANGLE (OEIS A036969) via the recurrence b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}. VERIFIED 22/22 cases across n=4..9.
  (2) **SHIFT PRINCIPLE**: C_{t_{2j+1}}(r) at n = C_{t_{2j-1}}(r) at n-2. Follows as corollary of (1).
  (3) **COMPLETE n=8 COEFFICIENT TABLE**: w_7=40320, w_5=-20160+1440t3, w_3=3024-480t3+48t5+96bc, w_1=-124+34t3-8t5+4t7-16bc+8bc35. All EXACT (0 error).
  (4) **PENALTY VERIFIED at n=9**: penalty = -30 + (21/2)*t3 + 3*t7 + 6*bc35 + 12*a3. Zero error.
  (5) **PREDICTIONS for n=11**: F_8(r) = 31 - 2640r^2 + 40320r^4 - 211680r^6 + 362880r^8. C_{t3}(r) at n=11 = 2*F_8(r). Testable computationally.
  (6) **Even-n W-polynomial**: Has only ODD powers of r. Top coefficient = n! (same as odd n). Same hierarchical structure.
**New contributions:** THM-059 (Universal Master Polynomial), w1_n8_complete.py, w_even_n_hierarchy.py, w_generating_function_test.py, universal_master_polynomial.py
**Unresolved threads:** Algebraic proof of central factorial recurrence; C_0(r) structure; F_8/F_9 predictions awaiting n=11 verification

## opus-2026-03-06-S11b (continued^5) — 2026-03-06 (Cross-scale perpendicularity + n=9 formulas)
**Account:** opus (overnight)
**Continuation of:** opus-2026-03-06-S11b (continued^4, context limit)
**Summary of work:**
  (1) **General tr(c_{n-5}) formula verified at n=9**: The decomposition into (4,) and (2,2) position patterns extends to all odd n, with coefficients following clean combinatorial formulas.
  (2) **OCF cycle counting convention clarified**: alpha_k in I(Omega(T),2) counts independent sets of DIRECTED cycles weighted by Hamiltonian cycle counts on each vertex set. At n=9, a 5-vertex set can have MULTIPLE directed 5-cycles (critical correction).
  (3) **EXACT c_2 formula at n=9**: c_2 = 462*t_3 - 60*t_5 + 12*t_7 - 120*bc33 + 24*bc35_w + 48*alpha_3 - 2640. Verified with zero error over 30 random tournaments.
  (4) **Coefficient hierarchy at n=9**: c_8 universal, c_6 depends on t_3, c_4 depends on (t_3,t_5,bc33), c_2 depends on (t_3,t_5,t_7,bc33,bc35_w,alpha_3). Each level introduces new OCF invariants.
  (5) **Perpendicularity explained**: At n=7, residual correlations between hierarchical invariants are EXACTLY ZERO (algebraic, not statistical). This is because t_5 enters via (4,) patterns and bc enters via (2,2) patterns — different position topologies.
  (6) **Complete OCF decomposition at n=7**: H = 1 + 2*(t_3+t_5+t_7) + 4*bc, where t_7 (7-cycle count) was previously missing. c_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc + 253/4.
**New contributions:** Updated THM-055 with n=9 formulas, trc2_exact_n9.py
**Unresolved threads:** c_0 formula at n=9 (needs more data), algebraic proof of singleton cancellation

## kind-pasteur-2026-03-06-S25g — 2026-03-06 (W-Hierarchy Spectral Decomposition + THM-058 PROVED)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25f (context limit)
**Files read:** SESSION-LOG.md, simplex-cube-perspectives-S25g.md, opus S28 scripts
**Summary of work:**
  (1) **w_0 = H - 3*t_3 VERIFIED** at n=5 exhaustively. FAILS at n=7 (needs more invariants).
  (2) **EXACT w_0 at n=7**: w_0 = 2*t_3 - t_5 + 2*t_7 - 2*alpha_2 - 17/4. Zero error over 15 samples.
  (3) **COMPLETE W-coefficient hierarchy at n=7** (ALL formulas exact, 0 error):
      w_6 = 5040, w_4 = 240*t_3 - 2100, w_2 = -60*t_3 + 12*t_5 + 24*alpha_2 + 231,
      w_0 = 2*t_3 - t_5 + 2*t_7 - 2*alpha_2 - 17/4.
      Each level adds one new cycle complexity. This is a spectral decomposition.
  (4) **THM-058 PROVED**: w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2] for ALL tournaments.
      Algebraic proof: non-adjacent sigma = 0, sigma_adj = (n-3)!*[2t_3-C(n,3)/2].
      Verified at n=5 (exhaustive), n=7 (20 samples), n=9 (15 samples, via DP).
  (5) **w_{n-5} decomposition**: Only (1,g,1) patterns contribute (massive vanishing).
      sigma_1111 = -12*t_3 + 4*t_5 + 42 (sees t_5 via consecutive positions).
      sigma_{1g1} = -8*t_3 + 8*alpha_2 + 35 for g>=2 (sees alpha_2 via separated pairs).
      This explains the hierarchy mechanism: t_5 enters via consecutive, alpha_2 via separated.
  (6) **P(n) = OEIS A093934**: Rooted tournament count. P(n) = 2*(n-1)! for n<=5, FAILS at n=6.
  (7) **Penalty shift principle**: H-w_0 penalizes t_3 at n=5 but t_5+alpha_2 at n=7.
  (8) **Creative synthesis**: Connection map linking tournaments to 6 mathematical categories.
**New contributions:** THM-058 (PROVED), W-coefficient-hierarchy-S25g.md, creative-synthesis-S25g.md, simplex-cube-perspectives-S25g.md, INV-082/083/084, T173-T176
**Unresolved threads:**
  - Prove sigma_1111 and sigma_{1g1} formulas algebraically (extend THM-058 approach)
  - Complete w_{n-5} proof for general n
  - Is there a topological interpretation of w_0 (Euler characteristic of some space)?
  - Check sigma_{1g1} independence of g at n=9

## kind-pasteur-2026-03-06-S25f — 2026-03-06 (Grand Synthesis + W(r) Stratification + Pfaffian Duality)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25e (context limit)
**Files read:** grand-synthesis-S25f.md (prior context), SESSION-LOG.md, INVESTIGATION-BACKLOG.md, opus S27 scripts
**Summary of work:**
  (1) **Grand Synthesis document** (`03-artifacts/drafts/grand-synthesis-S25f.md`): Comprehensive map of all mathematical structures — 5 equivalent algebraic perspectives (independence polynomial, transfer matrix, symmetric functions, Hopf algebra, tiling geometry), 6 novel creative connections (Mobius inversion, perpendicular plane/Mobius strip, Steiner systems, Cayley transform, groups/SC boundary, 1729 mystery).
  (2) **W(r) coefficient stratification** (`04-computation/W_coefficient_stratification.py`): Verified w_0 = -t_3 + 2*t_5 + 1 at n=5 (EXACT for all 11 iso classes). All odd-indexed W coefficients are exactly 0. W(r) stratifies tournament invariants by odd-cycle complexity.
  (3) **W(r)/OCF connection** (`04-computation/W_ocf_connection_v4.py`): At n=5, H = 1 + 2*(total directed odd cycles) since a_2 = 0 (no room for disjoint cycle pairs). Fixed cycle enumeration bug (canonical form filter was wrong).
  (4) **Recursive Hopf structure**: Discovered that overlap=3 contribution to w_{n-5} at n=7 uses OCF at n=5 as subroutine via the Hopf algebra coproduct. This creates a recursive hierarchy: W(r) coefficients at n use OCF at smaller n as building blocks.
  (5) **Pfaffian-path duality** (`04-computation/pfaffian_path_duality.py`): At n=4, det(S) is exactly determined by t_3 (det=9 iff t_3 odd). At n=6, needs finer invariants. Path-cycle duality: odd n has H (paths survive), even n has Pf (cycles survive).
  (6) **Paley tournament eigenvalue analysis**: Fixed Paley construction (requires p ≡ 3 mod 4). Paley T_7 has W(r)/7! = [1/320, 0, 1/80, 0, 1/4, 0, 1]. All non-trivial eigenvalues degenerate → scalar M.
  (7) **Deep connections document** (`03-artifacts/drafts/deep-connections-S25f.md`): Extended analysis of W(r)=tr(M(r)), Cayley transform, Pfaffian structure, BIBD embedding, transfer matrix as random walk, Hopf algebra recursion, 1729 number theory.
  (8) **Integrated opus S27 results**: THM-055 (coefficient hierarchy theorem) connects perfectly — e_{2k}(s_P) is polynomial in f_P, moments of f_P determine W coefficients.
**New contributions:** INV-079 (W(r) stratification), INV-080 (Pfaffian-path duality), INV-081 (Paley W(r) structure), grand-synthesis-S25f.md, deep-connections-S25f.md
**Unresolved threads:**
  - Explicit formula for w_0 at n=7 (depends on finer invariants than (t_3, t_5))
  - Prove Hopf algebra recursion algebraically
  - Does the 4th moment of f_P have a cycle-theoretic interpretation?
  - Compute W(r) for Paley T_11 and T_3
  - Formal path-cycle duality between H and Pf(S)?

## opus-2026-03-06-S11b (continued^4) — 2026-03-06 (Coefficient Hierarchy Theorem)
**Account:** opus
**Continuation of:** opus-2026-03-06-S11b (continued^3)
**Summary of work:**
  (1) **THM-055: Coefficient Hierarchy Theorem** — Proved that tr(c_{n-1-2k}) = sum_P e_{2k}(s_P) where e_{2k} is a polynomial of degree 2k in f_P (forward arc count). Via Newton's identity reduction, all power sums p_j of s_i = ±1/2 collapse to functions of p_1 = f - (n-1)/2.
  (2) **Explicit formulas computed** via sympy for e_0, e_2, e_4, e_6 at n=5,7,9.
  (3) **Moment hierarchy discovered:** sum_P f^j for j≤3 depends only on t_3; sum_P f^4 depends on MORE than t_3. This is why tr(c_{n-5}) cannot be expressed via simple cycle counts.
  (4) **Position decomposition at n=7:** f_(4) splits into 4 types by position adjacency pattern: "4 consec" (3 subsets, 5 vertices = subtournament H sums), "3+1" (6 subsets, 6 vertices), "2+2" (3 subsets, 6 vertices), "2+1+1" (3 subsets, 7 vertices = irreducible whole-tournament invariant).
  (5) **At n=5:** f_(4) = H directly (only 1 subset of 4 positions), giving tr(c_0) = H - 3*t_3.
  (6) Read algebraic proof of c_{n-3} from opus-S27. Read grand synthesis from S25f.
**New contributions:** THM-055, coefficient_hierarchy_proof.py
**Unresolved threads:**
  - Prove sum_P f^3 depends only on t_3 algebraically (currently only verified computationally)
  - What is the explicit invariant that sum_P f^4 captures at n≥7?
  - Can the hierarchy be used to bound or relate c_0 to H?

## kind-pasteur-2026-03-06-S25e — 2026-03-06 (THM-052 DISPROVED for non-SC VT; McKay database)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-06-S25d (context limit, then continuation)
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, INVESTIGATION-BACKLOG.md, frobenius21_palindromic_N.py
**Summary of work:**
  (1) **THM-052 DISPROVED for non-SC VT tournaments:** Computed N(0,1,j) for the F_21 non-normal Cayley tournament at n=21 via prefix+suffix bitmask DP (1075s total). N is NOT palindromic (all 20 values differ from mirror), alternating sum = M[0,1] = 45,478,409 != 0. M is NOT scalar.
  (2) **McKay database validation:** Downloaded McKay's VT tournament data for n=21 (users.cecs.anu.edu.au/~bdm/data/). 88 circulant + 22 non-circulant = 110 total VT tournaments. ALL 88 circulant are self-converse, ALL 22 non-circulant are NOT self-converse. n=21 is the SMALLEST order with non-circulant VT tournaments.
  (3) **digraph6 decoder bug fixed:** Initial decoder skipped diagonal entries in digraph6 format, corrupting adjacency matrices. Fixed by advancing bit index on diagonal. Verified: all circulant tournaments now correctly decode as regular (all out-degree 10).
  (4) **3-cycle counts:** All 22 non-circulant VT tournaments have identical 3-cycle count (385), same as our F_21 construction. All are regular.
  (5) **Literature search:** El Sahili-Ghazo Hanna (2023): T and T^op have same oriented Hamiltonian path type distribution. Ai et al. (2025): converse-invariant digraph polynomial. Neither directly addresses position distributions.
  (6) **MISTAKE-013, MISTAKE-014 logged:** Self-converse assumption false for non-abelian VT; THM-052 scope must be restricted.
**New contributions:** mcKay_vt21_selfconverse.py, MISTAKE-013, MISTAKE-014, INV-077, INV-078
**Unresolved threads:**
  - Which of McKay's 22 non-circulant VT tournaments corresponds to our F_21 construction?
  - What is the exact M matrix structure for the non-SC tournament? (only M[0,1] computed)
  - Verify Aut+Anti characterization at n=7 exhaustively
  - Does any structure theorem hold for non-SC VT M matrices?

## opus-2026-03-06-S26 (continued²) — 2026-03-06 (position-uniform, converse disproof, even-r hierarchy)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S26 (continued) (context limit)
**Summary of work:**
  (1) **THM-052 CONVERSE DISPROVED at n=5:** Scalar M does NOT imply vertex-transitive. Counterexample: "cone over C_3" with scores (1,2,2,2,3), |Aut|=3, but M=3I. Has anti-automorphisms swapping poles.
  (2) **POSITION-UNIFORM <=> SCALAR M at n=3,5 (exhaustive):** N[v,j] = H/n for all v,j is exactly equivalent to M=(H/n)*I. 64/1024 labeled tournaments on 5 vertices satisfy this.
  (3) **Aut+Anti transitivity analysis:** Aut(T) union Anti(T) transitive => scalar M (= VT for tournaments), but NOT conversely. Non-VT scalar-M tournaments have 2 orbits under Aut+Anti.
  (4) **Even-r polynomial for non-VT PU:** VT tournaments have EACH c_k individually scalar. Non-VT PU tournaments have non-scalar c_0 and c_2 that cancel at r=1/2. Key: c_0 off-diagonal = 0.5 cancelled by c_2*0.25 = -0.5.
  (5) **Cone construction does NOT generalize to n=7:** H=135, H%7=2. No non-circulant PU found at n=7 (1000 random + 84 flips tested).
  (6) **FORMULA: tr(c_{n-3}) = 2*(n-2)!*t_3 + const(n).** Verified at n=5 (12*t_3-30) and n=7 (240*t_3-2100). Constants: -0.5, -30, -2100.
  (7) **tr(c_2) at n=7 depends on full iso class,** not just cycle counts. Varies from -57 to +63 within same score sequence.
  (8) **n=9 even-r: top 3 coefficients universal for regular.** tr(c_8)=362880=9!, tr(c_6)=90720, tr(c_4)=6480 all universal. Only tr(c_2), tr(c_0) vary.
  (9) **Even-r vs OCF connection:** H = tr(c_0) + tr(c_2)/4 + ... refines the OCF H = 1 + 2*alpha_1 + 4*alpha_2. At n=5: tr(c_0) = H - 3*t_3.
**New contributions:** THM-052 converse disproof, position-uniform characterization, cone construction analysis, c_{n-3} formula, even-r hierarchy
**Unresolved threads:**
  - Prove tr(c_{n-3}) = 2*(n-2)!*t_3 + const algebraically
  - What is const(n)? Ratios are 60, 70 — pattern unclear
  - At n=7 (prime), is PU = VT? (evidence: no non-VT PU found)
  - What determines tr(c_2) at n>=7? (finer than score sequence)
  - Extend OCF connection: express tr(c_k) in terms of independence polynomial coefficients

## opus-2026-03-06-S11b (continued³) — 2026-03-06 (diagonal signed position + perpendicularity)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-06-S11b (continued²) (context limit)
**Summary of work:**
  (1) **DIAGONAL SIGNED POSITION THEOREM:** M[v,v] = sum_P (-1)^{pos(v,P)} where sum is over all Ham paths. Verified all 12 iso classes at n=5. M[v,v] can be NEGATIVE (not a path count). "Defect vertex" = vertex with position bias far from H/n.
  (2) **c_0 RANK STRUCTURE across H-classes at n=7:**
    - H=189 (Paley): c_0 = 2.25·I (SCALAR, rank 0)
    - H=175: c_0 = 0.25·I (SCALAR, rank 0)
    - H=171: c_0 has rank-2 perturbation, 1 defect vertex with diag=-3.75 vs 0.25
  (3) **PERPENDICULARITY CONFIRMED at n=7 (790 iso classes):**
    Mean cosine of M-directions = -0.0485 (near 0 = perpendicular).
    Low-H classes have positive cosine (aligned), high-H negative (anti-aligned).
    Crossover (true perpendicularity) at H≈95-105 (near median).
  (4) **ALL 43 H-maximizers at n=7 have M = 27·I** (scalar). Supports conjecture: max-H => VT => scalar M.
  (5) **Defect count decreases with H** within score class (1,2,2,2,3) at n=5:
    H=11: 2 defects, H=13: 1 defect, H=15: 0 defects (scalar).
  (6) **M(r) symmetry for ALL r** verified: each c_{2k} individually symmetric (polynomial identity). Already proved by THM-030.
  (7) **MISTAKE-012 corrected:** Blue pair ≠ tournament complement. Complement formula: M(T^c) = diag(M(T)) - offdiag(M(T)).
**New contributions:** diagonal_signed_position_theorem.py, perpendicularity_cosine_n7.py, c0_concentration_theorem.py, rank2_signed_adjacency_n7.py
**Unresolved threads:**
  - Prove diagonal signed position theorem from IE formula
  - Why does perpendicularity occur near median H? Is there a spectral explanation?
  - Extend rank-2 structure to n=9

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
