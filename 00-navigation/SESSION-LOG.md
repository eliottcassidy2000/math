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

## opus-2026-03-07-S35 (continued^5) — 2026-03-07 (Skeleton structure, Walsh-SC connection)
**Account:** opus
**Continuation of:** opus-2026-03-07-S35c4
**Summary of work:**
  (1) **SC/NSC skeleton analysis at n=5 and n=7.** Complete classification of tournament isomorphism classes into self-complementary (SC) and non-self-complementary (NSC) complement pairs.
  (2) **Signed S_n symmetry of Walsh spectrum (T181).** Proved: hat{H}[S] = (-1)^{|S cap F_sigma|} * hat{H}[sigma^{-1}(S)] for any vertex permutation sigma. Verified at n=5 (0 mismatches over 1024 monomials).
  (3) **Even-degree Walsh invisibility of SC/NSC (T182).** Proved that even-degree Walsh coordinates cannot distinguish SC from NSC — NSC pairs are "chiral dimers" invisible to H. At n=5: 3 H=3 classes share identical (C2,C4)=(-6,0) but 1 is SC, 2 form NSC pair.
  (4) **Palindromic score necessary but not sufficient (n>=7).** At n=5 palindromic score <=> SC. At n=7 this breaks: ~124 palindromic-score NSC classes exist.
  (5) **Orbit coordinate decomposition (T183).** H = 15/2 + 3/4*C2 + 1/8*C4 at n=5. t3 shares same degree-2 Walsh signs with ratio 3. C2 = 4*t3 - 10.
  (6) **SC constructive enumeration at n=7.** Using Moon's theorem (SC witness has type (2,2,2,1)), enumerated all 4096 labeled SC tournaments with fixed sigma. Found 65 distinct (H,t3,score) signatures, H ranging 1-189.
  (7) **Position function symmetry.** For SC tournaments: f_v(T) = f_{sigma(v)}(T^op) where sigma is SC witness. NSC pairs have different position multisets.
**New contributions:** T181, T182, T183
**Unresolved threads:** Complete SC enumeration at n=7 (need all 88 classes); palindromic-score NSC characterization; connection between skeleton and transfer matrix symmetry

## opus-2026-03-07-S39 (continued) — 2026-03-07 (H=21 structural impossibility analysis)
**Account:** opus
**Continuation of:** opus-2026-03-07-S39
**Summary of work:**
  (1) **THM-078: u_T size-weighted independence polynomial** (renamed from THM-077 to avoid collision with S35c4's Walsh OCF proof).
  (2) **THM-079: H=21 component reduction PROVED.** Part A: disconnected Omega ruled out (requires I(component)=7, blocked by THM-029). Part B: P_4 realization ruled out (two 3-cycles sharing vertex on 5 vertices always force 3rd cycle). Corollary: H=21 requires connected Omega with >=6 vertices.
  (3) **I(P_4, 2) = 21 discovery.** Path graph on 4 vertices gives independence polynomial value exactly 21 at x=2. Key target structure for H=21.
  (4) **Graph classification for I(G,2)=21.** Exhaustive enumeration: v=4 only P_4, v=5 none, v=6 only K_6-minus-2-edges (13 edges, 2 non-edges).
  (5) **n=6 exhaustive (alpha_1, alpha_2) analysis.** All 6 decompositions of alpha_1+2*alpha_2=10 have zero tournaments. Missing sums {3,10,17,19} correspond to gaps {7,21,35,39}.
  (6) **Score sequence analysis for t_3=4 at n=9.** 26 score sequences found; 0 of 1M random samples achieved exactly 4 three-cycles in P_4 arrangement.
  (7) **t_3=4 forces t_5>0 at n=7.** 0 of 17,837 tournaments with t_3=4 had t_5=0.
**New contributions:** THM-079, graph I=21 classification, P_4 impossibility proof
**Unresolved threads:** K_6-minus-2-edges route still open; need n>=9 analysis; H=21 full proof remains conjectural

## opus-2026-03-07-S35 (continued^4) — 2026-03-07 (THM-076 general-r + THM-077 new OCF proof)
**Account:** opus
**Continuation of:** opus-2026-03-07-S35 (continued^3)
**Summary of work:**
  (1) **THM-076 general-r proof COMPLETE.** Used EGF on set partitions: F(t,x)=(1+bt)/(1-at), giving G_r(x) = 2r!/((1-x)^r(1+x)). The half-weighted partial sum S(j)=1/2 for ALL j (alternating sums telescope). Hockey-stick identity gives Sigma=(m+r)!/m!. Crucially, Sigma depends ONLY on m and r, NOT on individual path sizes a_i. Verified r=1..5, n up to 20.
  (2) **THM-077: NEW ELEMENTARY PROOF OF OCF.** Direct Walsh analysis proves H(T) = I(Omega(T), 2) without using GS machinery. H side: even-length path components force contiguity in any HP (each internal vertex has degree 2). Sign product (-1)^{2a}=1, so block count = 2^r*(n-2k)!. I side: THM-076 gives same amplitude. Sign: epsilon_S = (-1)^{asc(S)} on both sides independently. hat{H}[S] = hat{I}[S] for ALL S. QED.
  (3) **Sign formula discovered.** hat{H}[S] = (-1)^{asc(S)} * 2^r*(n-2k)!/2^{n-1}, where asc(S) counts ascents in path components. Consistent with THM-068 (PCD descent-sign).
  (4) **Hypercube/4D connection clarified.** At n=7, H lives in effective 4-parameter Walsh space: degree 0, degree 2, degree 4 (two amplitude types: P_4 and P_2+P_2), degree 6. This is the "4D structure" noted by previous sessions.
  (5) **Skeleton analysis.** GS flip = toggle non-backbone bits on hypercube. Complement invariance = even-power vanishing = H(T)=H(T^op). NSC classes are those with non-palindromic score sequences.
**New contributions:** THM-076 (complete general-r proof), THM-077 (new OCF proof), T179, T180
**Unresolved threads:** Write up THM-077 for paper; formal proof that contiguity extends to multi-component (currently only verified); extend OPCD to general degree

## opus-2026-03-07-S38 — 2026-03-07 (H=21 gap proved n<=7, tangent number proof, Mitrovic DC explored)
**Account:** opus
**Continuation of:** opus-2026-03-07-S37
**Summary of work:**
  (1) **THM-075: H=21 permanent gap PROVED through n=7.** Exhaustive computation of all 2,097,152 tournaments at n=7 confirms H=21 never occurs. Complete H-spectrum at n=7 has 77 distinct values (all odd, range 1-189). OCF constraint analysis: none of the valid (alpha_1, alpha_2) decompositions for H=21 are achievable at n=6 or n=7. Fixed cycle counting bug (must count ALL directed cycles per vertex set, not just one).
  (2) **Tangent number connection PROVED (INV-093).** P_n(0,0) = 2^{(n-1)/2} * T_n follows from evaluating Eulerian EGF at t=i. Agent found Hetyei (2017) paper connecting tournaments to median Genocchi numbers.
  (3) **Type-count = A000009 (INV-092).** OCF cycle types at size n biject with partitions into distinct parts via removing 1's. Null-dim sequence not in OEIS (novel).
  (4) **Mitrovic noncommutative deletion-contraction explored (INV-094).** W_X = W_{X\e} - W_{X/e}^up. Edge contraction T/e is NOT a tournament. OCF fails for T\e and T/e. Not useful for proving OCF.
  (5) **Bags-of-sticks for OCF: DEAD END (INV-095).** Under OCF specialization, every bag of sticks contributes 1. Decomposition gives no new info.
  (6) **Even-cycle cancellation explained.** Factor 2^k in OCF = orientation multiplicity. Each odd cycle has 2 directed orientations (T and T^op), both contributing sign +1.
  (7) **GS-OCF bridge to P(u,x) analyzed.** U_T and P(u,x) are different objects encoding overlapping info. Neither is a specialization of the other. They share the independence polynomial as common shadow.
  (8) **Degree-4 Fourier at n=9 re-analyzed.** Agent found 2D structure (same as n=7), contradicting S37's rank>>200 finding. Needs reconciliation.
**New contributions:** THM-075, INV-091-095, complete n=7 H-spectrum, OPEN-Q-019 update
**Unresolved threads:** Prove H=21 at general n; reconcile degree-4 n=9 dimensionality; is 63 a permanent gap?; THM-063 vs THM-074 contradiction on G_T(t,2)=E_T(t)

## opus-2026-03-07-S35 (continued^3) — 2026-03-07 (THM-076: Walsh-OCF Factorization)
**Account:** opus
**Continuation of:** opus-2026-03-07-S35 (continued^2)
**Summary of work:**
  (1) **THM-076: Walsh-OCF Factorization Identity (PROVED).** The Walsh amplitude formula for I(Omega(T), 2) factorizes: each Walsh monomial's coefficient decomposes into "covering configurations" — ways to assign the monomial's path components to disjoint odd cycles. The remaining vertices contribute E[I(Omega,2)] = f(n-used). Each term telescopes via the constant-term identity C(m,j)*j!*(m-j)! = m!.
  (2) **Single-component proof complete.** For any P_{2a} path, the factored sum gives exactly 2*(n-2a)!/2^{n-1}. Verified for a=1,...,5 and n up to 19.
  (3) **Multi-component proof verified.** For r=2 (P2+P2, P2+P4) and r=3 (P2+P2+P2), the factorization works by enumerating covering types (one big cycle, separate cycles, mixed). All verified at multiple n values including brute-force cycle counting at n=11.
  (4) **Cycle polynomial structure discovered.** t_k (directed k-cycle count) is a polynomial of degree k-1 in edge variables (odd k: top degree cancels between forward/reverse). Walsh degrees 0,2,...,k-1 only (complement invariance kills odd degrees).
  (5) **Walsh-orthogonal OCF invariants.** At n=5: q2 = t3 (centered), q4 = t5 - (1/2)*t3 (centered). The mixing coefficient t5_hat/t3_hat = 1/2 exactly at degree 2. OCF invariants have Walsh content at multiple degrees (t5 contributes to both degree 2 and 4).
  (6) **Key insight: higher-degree OCF follows from degree-0 OCF.** The factorization uses E[I(Omega,2)] on remaining vertices, which is degree-0 OCF. This creates a self-consistent tower: degree 0 implies degree 2, which (with appropriate bookkeeping) implies degree 4, etc.
**New contributions:** THM-076, extensive verification scripts
**Unresolved threads:** Full general-r algebraic proof; formal proof that covering types exhaust all contributions; connection to Grinberg-Stanley proof structure

## kind-pasteur-2026-03-07-S30 — 2026-03-07 (Factor-2 explained, u_T polynomial, tangent numbers, agent investigations)
**Account:** kind-pasteur
**Continuation of:** kind-pasteur-2026-03-07-S29
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, INVESTIGATION-BACKLOG.md, THM-065, THM-074, THM-002, THM-073, gs_specialization_check.py, typed_GT.py, bags_of_sticks.py, opus MSG-041
**Summary of work:**
  (1) **Factor-2 in OCF EXPLAINED.** ps_1(U_T)(1) = H(T) because: even-cycle permutations cancel (-1 from T-direction + 1 from T^op), odd-cycle perms contribute +1 from BOTH directions. The 2^psi in OCF = T + T^op contributions. Verified n=3,5.
  (2) **u_T(m) = ODD polynomial in m.** Only odd powers appear (m, m^3, m^5,...). u_T(m) = m*Q_T(m^2). Coefficients = f-level weighted sums S_f from THM-065. Verified exhaustive n=5.
  (3) **P_n(-2,0) = T_n exactly.** Base polynomial evaluates to tangent number at u=-2 (in addition to 2^m*T_n at u=0). P_n(u,0) interpolates: T_n at u=-2, 2^m*T_n at u=0, n! at u=2. Updated THM-074.
  (4) **G_T is NOT a specialization of U_T** — confirmed by computation (transitive gives g1^5 = A_5(t), not a polynomial root) and background agent analysis.
  (5) **H=21 still absent at n=8** — 200k random samples, gap between H=19 and H=23 persists.
  (6) **Schweser-Stiebitz-Toft (2510.10659)** — agent found expository/historical only, no connections to our framework. One new paper to catalog: Ai et al. 2407.17051.
  (7) **Bags-of-sticks/deletion-contraction** — agent confirmed dead end for alternative OCF proof. Edge-DC produces non-tournaments where OCF fails.
  (8) **THM-065 updated** with A078408 OEIS identification and extended null_dim table to n=21.
**New contributions:** factor2_explained.py, UT_to_GT_bridge.py, THM-065 update, THM-074 tangent number update
**Unresolved threads:** H=21 permanent gap needs exhaustive n=8 check; u_T(m) root structure; Irving-Omar matrix algebra (agent hit rate limit); Ai et al. 2407.17051 to catalog

## kind-pasteur-2026-03-07-S29 — 2026-03-07 (THM-074 master decomposition, deep core review)
**Account:** kind-pasteur
**Continuation of:** kind-pasteur-2026-03-05-S3
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, INVESTIGATION-BACKLOG.md, all opus inbox messages (MSG-039, MSG-040), typed_GT.py, gs_specialization_check.py, bags_of_sticks.py, P_hierarchy_general.py, ocf_factor2_investigation.py, degree4_n9_fast.py, THM-073-gs-typed-independence.md
**Summary of work:**
  (1) **THM-074: Master Decomposition of P(u,x) — PROVED.** Created and verified the complete factorization: g_I(u) = P_{n-S}(u,0) * (u-2)^{S/2}. Verified exhaustive at n=5, random at n=7. Key: corrections depend only on S = sum(l_i-1), not individual cycle lengths.
  (2) **Tangent number connection discovered:** P_n(0,0) = 2^{(n-1)/2} * T_n where T_n is n-th tangent number. Verified n=1,3,5,7,9,11,13. Connects tournament base polynomial to alternating permutations.
  (3) **Null space dimension formula verified:** null_dim = #{cycle types} - #{distinct S values}. Gives sequence 0,0,1,3,6,11,19,29,44,65. Types grow as partitions into odd parts ≥3, S-values grow linearly as floor(n/2)+1.
  (4) **G_T(t,x) ≠ E_T(t) CORRECTED.** typed_GT.py line 183 claim "G_T^typed(t;2,2,...) = E_T(t)" is FALSE. G_T is the "inflated independence polynomial": G_T(0,x)=I(Omega,x), G_T(1,x)=n!, G_T(t,0)=A_n(t). NOT the HP descent polynomial.
  (5) **THM-073 number collision resolved:** Renamed kind-pasteur's THM-073 to THM-074 to avoid collision with opus's THM-073 (GS typed independence).
  (6) **Investigated opus messages thoroughly:** Read MSG-039 (THM-070 clean Claim A proof, GS specialization, Mitrovic-Stojadinovic analysis) and MSG-040 (THM-068 PCD proved all degrees, THM-072 OPCD discovered). Ran and verified gs_specialization_check.py and bags_of_sticks.py outputs.
**New contributions:** THM-074 (master decomposition), tangent number connection, null space formula verification, G_T≠E_T correction
**Unresolved threads:** Null space dimension OEIS lookup; H=21 permanent gap; GS-OCF bridge to P(u,x); degree-4 Fourier at n=9

## opus-2026-03-07-S37 — 2026-03-07 (THM-073, degree-4 rank explosion, opp_edges bugfix)
**Account:** opus
**Continuation of:** opus-2026-03-07-S36 (context limit continuation)
**Summary of work:**
  (1) **THM-073: GS generator expansion = typed independence polynomial (PROVED).** The Grinberg-Stanley expression U_T = poly(p₁, 2p₃, 2p₅, ...) with non-neg integer coefficients directly encodes the typed independence polynomial. Coefficients α_S count collections of s₃ disjoint 3-cycles, s₅ disjoint 5-cycles, etc. ps₁(U_T)(1) = I(Ω(T), 2) = H(T).
  (2) **BUGFIX: opp_edges in bags_of_sticks.py and gs_specialization_check.py.** Both scripts had `opp_edges = {(j,i) : (i,j) ∉ T}` which equals T itself (not T^op!). Fixed to `{(j,i) : (i,j) ∈ T}`. This caused prior "OCF specialization" checks to compute a T-only sum rather than the actual specialization of U_T. The correct U_T has ALL-ODD partition support.
  (3) **Degree-4 Fourier at n=9: dim >> 200.** P4 monomials NOT proportional (117/126 anomalous). Fourier proof strategy is INFEASIBLE for n≥9 middle degrees.
  (4) **Mixed-direction cycle analysis at n=9.** Mixed T/T^op disjoint cycle pairs exist abundantly at n=9 (100+ per tournament) but not at n=5,7. The full U_T sum (T+T^op, signed) still gives H(T) via cancellation.
  (5) **Three equivalent H(T) formulas clarified:** (a) Held-Karp, (b) I(Ω(T),2), (c) ps₁(U_T)(1). The OCF specialization p₁→1, p_odd→2, p_even→0 applied to full U_T does NOT give H(T); it's ps₁ at m=1 that works.
  (6) **THM-065 updated with null space formula, Irving-Omar attribution fixed, new papers cataloged.**
  (7) **Typed independence polynomial more refined than standard:** 203 distinct types vs 166 standard signatures at n=7 (50k samples). 34 I_std values map to multiple I_typed versions.
**New contributions:** THM-073, INV-088-090, THM-065 update, attribution fixes, bugfixes
**Unresolved threads:** bags-of-sticks deletion-contraction OCF path; U_T at n=9 (verify all-odd support); typed independence polynomial applications

## opus-2026-03-07-S35 (continued^2) — 2026-03-07 (PCD General Proof + OPCD)
**Account:** opus
**Continuation of:** opus-2026-03-07-S35 (continued)
**Summary of work:**
  (1) **THM-068 PCD: PROVED AT ALL DEGREES, ALL ODD n.** Complete algebraic proof via block-placement. Key insight: all macro-items have ODD width, so start-position parity = macro-position parity. Constant descent sign cancels in ratio. Signed position sum factors through (-1)^{internal_offset}. Alternating sum = 1 for odd n-2k.
  (2) **THM-072: Off-diagonal PCD (OPCD).** Mc_hat[a,b,S] characterized at even Walsh degrees. Interior vertices (deg_S >= 2) have zero rows/cols. Cross-block endpoint, endpoint-free, free-free entries follow 1:2:4 ratio pattern. Formula: 1/((n-2k)(n-2k-1)) scaled by block orientation factors.
  (3) **Mc symmetry at even Walsh degrees PROVED.** Block-placement bijection (reverse macro-perm + flip orientations). Antisymmetric part lives ONLY at odd Walsh degrees.
  (4) **Dimension reduction:** At n=5, H(T) on 1024-point hypercube has only 3 independent amplitude parameters (constant, t3 coeff, t5 coeff). Walsh spectrum: 91 nonzero coefficients but only 3 distinct amplitudes.
  (5) Updated THM-068 from CONJECTURED to PROVED status.
**New contributions:** THM-068 (upgraded to full proof), THM-072 (new, was THM-070 before renumbering)
**Unresolved threads:** Full off-diagonal formula at general degree; rank = n-2k proof; relationship to Grinberg-Stanley OCF proof

## opus-2026-03-07-S36 — 2026-03-07 (THM-070 clean Claim A proof, Mitrovic-Stojadinovic, GS specialization)
**Account:** opus
**Continuation of:** opus-2026-03-07-S34 (context limit, then S36 continuation)
**Summary of work:**
  (1) **THM-070: Clean Claim A from OCF (PROVED).** 4-step proof:
      Through-v cycles form clique in Omega(T) → each indep set has at most one through-v cycle → graph equality gives f(C)=2mu(C) → Claim A. NO inclusion-exclusion needed.
  (2) **Through-v clique theorem.** Any two cycles through v share vertex v, hence adjacent. This trivially makes higher-order IE terms zero.
  (3) **THM-069 renumbering.** Fixed collision: graph equality = THM-069, Walsh-Fourier = THM-071.
  (4) **GS OCF specialization verified.** p_1→1, p_{odd≥3}→2, p_even→0 gives H(T) from U_T coefficients. Confirmed for C_5 (H=15) and Paley T_7 (H=189).
  (5) **Mitrovic-Stojadinovic paper analyzed (arXiv:2506.08841).** Key connections:
      - X_{inc(P)} = omega(U_P) (chromatic = omega of Redei-Berge)
      - Bags-of-sticks decomposition for U_X
      - Their phi(sigma) = our S = sum(l_i-1)
      - Noncommutative deletion-contraction for W_X
      - Connects Stanley-Stembridge conjecture to h-positivity of U_P
  (6) **Typed independence polynomial at n=7.** S=4 collision between (5,) and (3,3) types confirmed computationally. Only size-2 independent sets at n=7 are disjoint 3-cycle pairs.
  (7) **INV-033 updated** with Mitrovic-Stojadinovic findings.
**New contributions:** THM-070, typed indep scripts, GS specialization scripts, bags-of-sticks analysis
**Unresolved threads:**
  - Can bags-of-sticks decomposition give direct OCF proof?
  - Degree-4 Fourier structure at n=9 (agent launched)
  - Does noncommutative deletion-contraction on W_T prove Claim A directly?
  - General null space dimension formula (agent working)

## kind-pasteur-2026-03-07-S28 (continued) — 2026-03-07 (P(u,x) hierarchy, Mersenne vanishing, H=21 gap, THM-065 n=9)
**Account:** kind-pasteur
**Continuation of:** kind-pasteur-2026-03-07-S28 (context continuation)
**Summary of work:**
  (1) **THM-065 verified at n=9:** Null space dimension = 3, all 3 predicted null vectors confirmed (250 random tournaments, max |C*v| < 1e-11).
  (2) **THM-067: c_1 Formula and Mersenne Vanishing (PROVED).** c_1^(f,d) = 2^{f+1} - d - 2. Vanishes iff n = 2^{f+1}-1 (Mersenne numbers: 1, 7, 31, 127, ...). Explains why p_{m-1}(x) is linear at n=7.
  (3) **P(u,x) coefficient hierarchy.** Complete formulas at n=5 (exhaustive) and n=7. p_m(x) = I(Omega,x), p_{m-1}(x) linear at n=7 (bc coefficient = 0 due to c_1^(2,6)=0).
  (4) **H=21 permanent gap evidence.** Absent at n=7 (100k samples), n=8 (200k), n=9 (50k). At n=7, no valid (alpha_1, alpha_2) decomposition gives H=21: alpha_1=10 forces alpha_2=2, alpha_1=8 forces alpha_2=0, etc. All paths blocked.
  (5) **Investigated opus's THM-068/069/070.** Clean 4-step Claim A proof, position character decomposition, graph equality.
  (6) **GS typed IP synthesis.** U_T = p_1^n * I_typed(Omega; 2p_k/p_1^k), G_T(t,x) is genuinely new beyond U_T.
**New contributions:** THM-067, P_hierarchy_general.py, null_space_n9.py, h21_alpha_constraints.py
**Unresolved threads:** H=21 proof (rigorous impossibility), alpha constraint structure for general n, typed G_T ↔ P(u,x) connection

## opus-2026-03-07-S35 (continued) — 2026-03-07 (Position Character Decomposition)
**Account:** opus
**Continuation of:** opus-2026-03-07-S35 (context continuations)
**Summary of work:**
  (1) **THM-068: Position Character Decomposition (PCD).** M[v,v]_hat[S] = (-1)^{[v in N(S)]} * H_hat[S] / (n-2k) for Walsh degree 2k. |N(S)| = k, C(n,k) patterns per degree.
  (2) **Degree-2 PROVED ALGEBRAICALLY for all odd n.** Key: 2*(n-3)! valid permutations per block position, descent count always 1, alternating position sum gives -1/(n-2) ratio for shared vertex.
  (3) **Even-length component rule PROVED.** H_hat[S] ≠ 0 iff every path component of S has even length. Odd-length components have opposite descent parity → exact cancellation. Verified at n=7.
  (4) **Walsh degree ↔ OCF invariant correspondence.** P_{2j} ↔ (2j+1)-cycle, P_{2j1}+P_{2j2} ↔ disjoint odd-cycle pairs. Partition count = p(k) matches invariant count.
  (5) **Off-diagonal rank pattern.** M1[a,b]_hat has rank n-2k at degree 2k. M1 (as defined) is NOT symmetric; OCF transfer matrix is different.
  (6) **Fixed w_6 constant** in THM-069: 5040 = 7! (not 720 = 6!), w_0 constant = -17/4.
**New contributions:** THM-068, pcd_verification.py
**Unresolved threads:** Algebraic proof of PCD at degree 4+; OCF proof via Walsh framework; exact off-diagonal transfer matrix definition

## opus-2026-03-07-S34 — 2026-03-07 (GS-OCF Bridge, Graph Equality THM-067, f(C)=2mu(C))
**Account:** opus
**Continuation of:** opus-2026-03-07-S33 (context limit, then S34)
**Summary of work:**
  (1) **GS-OCF Bridge (PROVED).** Verified exhaustively at n=5: Grinberg-Stanley's U_T in power-sum basis equals the TYPED independence polynomial of Omega(T). U_T = sum_{S indep} 2^|S| * prod p_{len(c)} * p_1^{n-sum_len}. Each directed k-cycle contributes exactly 2 to the p_k coefficient.
  (2) **G_T(t,x) is genuinely new beyond U_T.** The power-sum expansion LOSES descent information. Our G_T requires the quasisymmetric structure, not just power sums.
  (3) **Null space explained via S = sum(l_i - 1).** The correction function g(t) = A_{f+1}(t)*(t-1)^{d-f} depends ONLY on S, not individual cycle lengths. This explains why the forward-edge distribution can't distinguish cycle types with the same S. Null dims: n=3:0, 5:0, 7:1, 9:3, 11:6, 13:11. Confirms and extends THM-065.
  (4) **THM-067: Graph Equality (PROVED).** For any cycle C through v: Omega(T)\N[C] = Omega(T-v)|_{avoid C\{v}}. Trivial proof: v in C implies any cycle disjoint from C avoids v. Immediately gives f(C) = 2*mu(C) for ALL tournaments.
  (5) **f(C) = 2*mu(C) verified at n=5 (exhaustive), n=7 (2664 pairs), n=9 (30511 pairs).**
  (6) **Claim A structural mechanism clarified:** alpha_0 = 1 (empty set) cancels in deletion, leaving only cycles through v. Each cycle contributes independently via f(C) = 2*mu(C). The higher-order inclusion-exclusion terms vanish (implied by Claim B).
**New contributions:** THM-067, 12+ computation scripts, GS-OCF bridge analysis
**Unresolved threads:**
  - Connection between THM-067 graph equality and THM-068 position character decomposition
  - Does the S = sum(l_i-1) grouping connect to the Walsh-Fourier degree structure?
  - Can typed G_T with separate y_k weights reveal more structure?
  - General n formula for null space dimension (related to odd partitions)

## opus-2026-03-07-S35 — 2026-03-07 (Walsh-Fourier Diagonalization PROVED at n=7)
**Account:** opus
**Continuation of:** opus-2026-03-06-S11b (many continuations)
**Summary of work:**
  (1) **THM-069 (formerly THM-067): Walsh-Fourier Diagonalization (PROVED at n=5,7).** The OCF H(T) on the Boolean hypercube {0,1}^m decomposes into Walsh degree components D_{2j}, each with Walsh degree EXACTLY 2j. The Fourier coefficients w_k = 2^k * D_{n-1-k} have pure Walsh degree.
  (2) **Analytical proofs of all Walsh coefficients at n=7:** t3_hat=-1/4, t5_hat=-3/4, t7_hat=-3/8, bc_hat=-1/4 at degree 2. t5_hat=1/16, t7_hat=1/32, bc_hat=0 (path) and t5=0, t7=1/16, bc=1/16 (double-fan) at degree 4.
  (3) **Generalized telescoping:** t7_hat = (t5+2bc)_hat/2 at degree 4, extending t5_hat = t3_hat/2 at degree 2. This connects to THM-065 f-level grouping.
  (4) **Key insight:** The Walsh basis diagonalizes the Fourier decomposition because f-level sums are the natural variables, and each cycle's contribution at lower Walsh degrees telescopes via the factor-of-2 identity.
**New contributions:** THM-067 (originally THM-066, renamed due to collision)
**Unresolved threads:** General n proof; degree-6 structure of w_0 at n=7; connection to skeleton eigenvalues

## opus-2026-03-07-S33 — 2026-03-07 (Trivariate GF, empty set insight, reduced polynomial P(u,x))
**Account:** opus
**Continuation of:** opus-2026-03-07-S32 (context limit, then S33 continuation)
**Summary of work:**
  (1) **THM-063: Trivariate GF G_T(t,x) (PROVED).** The generating function
      G_T(t,x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t)(t-1)^{d-f}
      has clean special evaluations:
      - G_T(t,0) = A_n(t) [Eulerian poly, T-independent]
      - G_T(0,x) = I(Omega(T), x) [INDEPENDENCE POLYNOMIAL!]
      - G_T(1,x) = n! for all x [T-independent]
      - G_T(t,2) = E_T(t) = sum a_k t^k
      Proof: (-1)^{d-f} = 1 always (d-f always even for cycle collections).
  (2) **Palindromic symmetry:** G_T(t,x) = t^{n-1} G_T(1/t, x) for ALL x.
      Generalizes the known a_k = a_{n-1-k} palindromy.
  (3) **Reduced polynomial P(u,x):** Writing u = t+1/t, G_T = t^m P(u,x).
      The LEADING coefficient p_m(x) = I(Omega(T), x) (independence poly!).
      P(2,x) = n! for all x. Explicit OCF formulas at n=5 (exhaustive) and n=7.
  (4) **Null space discovery:** The forward-edge distribution a_k(T) has a rank
      deficiency: changing (t5, bc) by (-2, +1) leaves all a_k unchanged.
      The forward-edge distribution cannot distinguish "2 five-cycles" from
      "1 disjoint 3-cycle pair." Verified computationally.
  (5) **E_T(-1) = deformed zigzag number.** Verified at n=3,5,7,9.
      Coefficient pattern: 2^{parts+2S} * E_{n-2S} where E_j = tangent numbers.
  (6) **E_T(i) alternates purely real/imaginary** by n mod 4. Verified.
  (7) **Bilinear structure:** G_T is NOT a product A_n(t)*I(Omega,x).
      Delta = G_T - A_n*I vanishes on both axes, equals -n!*(I-1) at t=1.
      N_T(t,2) = E_T(t)/A_n(t) smoothly interpolates from H to 1.
**New contributions:** THM-063-trivariate-gf.md, 10+ computation scripts
**Unresolved threads:**
  - General n formula for p_j(x) coefficients in the reduced polynomial
  - Does the null space structure persist/grow at n=9? (more null vectors?)
  - Combinatorial interpretation of p_j for j < m
  - Connection between P(u,x) structure and Claim A
  - Literature search for "deformed Eulerian" or bivariate Eulerian-independence GF

## kind-pasteur-2026-03-07-S27 — 2026-03-07 (W(i/2)=0 characterizes H-maximizers; deep core review)
**Account:** kind-pasteur
**Continuation of:** kind-pasteur-2026-03-07-S26 (context limit)
**Summary of work:**
  (1) **THM-064: W(i/2) Evaluation (NEW).** Discovered that W(i/2)=0 characterizes
      H-maximizers at n=3,5,7. Equivalent to (x^2+1)|A(x) where A(x) is the
      forward-edge generating polynomial. Fails at n=11 (Paley has W(i/2)=-10010).
      At even n, W(i/2) is pure imaginary, never zero.
  (2) **Free position universality:** Tested two hypotheses (random tournament avg,
      backbone model avg) — both WRONG. F_f != avg W in either model. The
      universality is deeper: it involves permutation sums, not just edge sums.
      Binary sign sums give (2r)^f (trivially position-independent), but F_f
      involves Eulerian weighting from permutation structure.
  (3) **Agent results integrated:**
      - Positivity: a_k=0 only at n=3 (3-cycle). For n>=4, always a_k>0.
      - n=7 constraints: n=5 identity (odd t3 => W(0)=0) doesn't generalize.
        W(0) always has frac part 3/4 at n=7.
      - W(i/2) is real at odd n, pure imaginary at even n.
      - W(i/2) = H - 6*t3 + 15 exactly at n=5.
  (4) **H=21 gap:** Confirmed absent at n=6 (exhaustive), n=7 (100K samples),
      n=8 (20K samples). H=7 and H=21 are the only gaps <= 25 through n=8.
  (5) **Tangent number connection verified:** sum(-1)^k A(n,k) = (-1)^{(n-1)/2} T_{(n+1)/2}
      connecting W(0) of transitive tournament to OEIS A000182.
  (6) **Synced with opus-S33:** Integrated Fourier OCF decomposition, bivariate
      Phi_T, reduced polynomial P(u,x), k-colored independence polynomial.
**New contributions:** THM-064, 10+ computation scripts
**Unresolved threads:**
  - W(i/2)=0 at n=9: does it still characterize maximizers?
  - Free position universality proof (permutation-level, not edge-level)
  - H=21 permanent gap proof
  - Background agents still running (n=9 W(i/2), Fourier patterns, deformed Eulerian properties, H=21 search)

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
