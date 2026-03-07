# Mistakes Log

**Purpose:** Every error that has been made and corrected — with enough context that no Claude instance ever repeats it. Read this before doing any computational or proof work.

Format per entry:
- What was assumed / done
- Why it was wrong
- The correct framing
- Impact on existing results
- Source (who found it, when)

---

## MISTAKE-001: μ Computation Bug in Scripts 6-9

**Date discovered:** ~2026-02-26 (pre-Cowork sessions); logged 2026-03-05
**Found by:** Claude instance (Account unknown, pre-Cowork era), reported in MASTER_FINDINGS.md
**Affects:** Scripts 6-9 (sum_mu() function); does NOT affect scripts 1-5

### What was assumed
`ind_poly_at_2_restricted()` correctly computes μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2) by iterating over cycles using `T[perm[i]][perm[(i+1)%len]]`.

### Why it was wrong
`perm` is a permutation of T−v vertices that does NOT include v. The cycle-checking code uses the **full T matrix** instead of T−v. Specifically, for a cycle of length L:
- `T[perm[i]][perm[(i+1)%len]]` checks arcs in the full tournament T
- But cycles of T−v (cycles not using v) should be checked in T−v
- The wrap-around `perm[(i+1)%len]` may also have indexing issues for general L

### The correct framing
The restriction to T−v must be done BEFORE checking cycles. Any cycle-finding in `ind_poly_at_2_restricted()` must operate on T−v's adjacency matrix exclusively.

### Impact assessment
- **Scripts 1-5 (inshat analysis):** UNAFFECTED (do not use μ computation)
- **Paper's Claim A verification:** STATED UNAFFECTED by MASTER_FINDINGS (verification runs used a separate code path)
- **Paper's Claim B verification:** Needs confirmation — Claim B involves I(Ω,2) directly
- **The verification table in the paper (0 failures for Claim A at n≤6):** Stated to be valid

⚠️ **UNRESOLVED:** It is not fully clear which code path the paper used for Claim A verification at n=6. If sum_mu() was used, the 0-failure count may not be reliable. This is the subject of DISC-001 in 02-court/active/.

### Lesson
When computing μ(C) for any cycle C in T:
1. First restrict to T−v (remove v and all incident arcs from the adjacency matrix)
2. Then find all odd cycles in T−v that are vertex-disjoint from C\{v}
3. Build Ω(T−v)|_{avoid C\{v}} from these cycles
4. Evaluate I(·, 2) on this restricted conflict graph

Never use the full T matrix when computing anything about T−v.

### Resolution
**Independent verification completed** (opus-2026-03-05-S1): tournament_lib.py implements mu(C) correctly per the above steps. Exhaustive verification of Claim A at n<=6 (196,608 pairs, 0 failures) confirms the paper's results are valid. See DISC-001 for the formal resolution.

---

## MISTAKE-004: RETRACTED — OCF IS a Valid Closed Form Over All Odd Cycles

**Date originally entered:** During file.txt exploration (pre-2026-03-05)
**RETRACTED by:** kind-pasteur-2026-03-05-S5 (DISC-002 resolved)
**Original claim:** H(T) = I(Omega(T), 2) where Omega(T) = ALL directed odd cycles is WRONG.

### Why the original claim was itself wrong

The alleged counterexample (T on 6 vertices with two disjoint 3-cycles C1, C2) computed I(Omega(T), 2) = 49 using MU WEIGHTS:
```
I_wrong = 1 + 2*mu(C1) + 2*mu(C2) + 4*mu(C1)*mu(C2) = 1 + 6 + 6 + 36 = 49
```
But the independence polynomial does NOT involve mu weights. The correct computation:
- alpha_0 = 1, alpha_1 = 2 (either C1 or C2), alpha_2 = 1 ({C1, C2} — vertex-disjoint, so non-adjacent)
- I(Omega(T), 2) = 1 + 2*2 + 1*4 = 9 = H(T) CORRECT

### The correct framing (replaces the false framing)

**H(T) = I(Omega(T), 2) IS a valid closed-form identity**, where:
- Omega(T) = conflict graph on ALL directed odd cycles of T
- Two cycles are adjacent in Omega(T) iff they share a vertex
- I(G, x) = sum_{k>=0} alpha_k * x^k is the plain independence polynomial (no mu weights)

This is equivalent to the recursive Claim A formulation — the closed form is obtained by unrolling the recursion. The mu weights mu(C) = H(T[V\V(C)]) arise in the recursion but NOT in the independence polynomial.

### Computational confirmation
H(T) = I(Omega(T), 2) verified exhaustively for n=3,4,5,6 (33,864 tournaments, 0 failures) by opus-2026-03-05-S2. Further confirmed by T_11 with H=95095 matching exact OCF calculation.

### What was confused (lesson for future agents)
The recursive mu-weighted formula H(T) = H(T-v) + 2*sum_C mu(C)*... uses mu weights at each step. When you UNROLL the full recursion, the mu weights become exactly the combinatorial weights in the independence polynomial (since mu(C) = H(T[V\V(C)]) which itself unrolls). The two formulations are equivalent; the closed form is NOT a non-recursive approximation. Do not confuse the per-step mu weights with the independence polynomial coefficients.

---

## MISTAKE-005: Cycle Bijection Under Arc Reversal Fails

**Date discovered:** During file.txt exploration
**Found by:** Claude instance (Account unknown)
**Affects:** Proof strategy for arc-reversal invariance D(T,v) = D(T',v)

### What was assumed
When T' is obtained from T by flipping arc i->j to j->i (with i,j != v), odd cycles through v containing i->j in T biject with odd cycles through v containing j->i in T', preserving V(C).

### Why it was wrong
A cycle C through v containing i->j can be written as v ->^{P1} i -> j ->^{P2} v. The "conjugate" C' would need j->i in T', requiring the path segments to be REVERSED. But reversing a directed path P2 does not generally give a valid directed path in the same tournament (arcs may point the wrong way).

### The correct framing
There is NO individual bijection C <-> C' preserving vertex sets. The arc-reversal invariance, if true, must hold as a SUM equality: sum_{C: i->j in C} H(T[V\V(C)]) = sum_{C': j->i in C'} H(T'[V\V(C')]). This is a weaker statement that doesn't require a cycle-by-cycle bijection.

### Impact
The arc-reversal invariance D(T,v) = D(T',v) remains the key unproved step for a general proof of Claim A. The cycle bijection approach is a dead end; a sum-level argument is needed.

---

## MISTAKE-002: Exact Path Formula H(T) = B_v + S_v + R_v

**Date discovered:** During script verification at n=4
**Affects:** Any approach trying to decompose H(T) exactly as B_v + S_v + R_v

### What was assumed
H(T) = B_v + S_v + R_v (exact equality)

### Why it was wrong
Verified FALSE: 96 failures out of 256 pairs at n=4.

### The correct framing
Only the PARITY version holds: H(T) ≡ B_v + S_v + R_v (mod 2). The exact identity is wrong.

### Impact
Similarly, S_v + R_v = 2Σμ(C) is FALSE (144/256 failures at n=4). These exact decompositions are not the right approach.

---

## MISTAKE-003: Per-Path Identity as a Path to Claim A for n≥6

**Date discovered:** During n=6 verification
**Affects:** Any proof strategy that relies on the per-path identity to prove Claim A for general n

### What was assumed
The per-path identity (inshat−1)/2 = Σ_{3-cycle embeddings} μ(C) could serve as an intermediate step toward proving Claim A for all n.

### Why it was wrong
The per-path identity fails for 2,758/9,126 ≈ 30% of (T,v,P') triples at n=6. It is provably a 3-cycle-only formula (by THM-005). Longer cycles contribute to Claim A's RHS but are invisible to the per-path identity.

### The correct framing
The per-path identity is valid and useful for n≤5 (where it implies Claim A). For n≥6, a different per-path formula is needed — one that accounts for contributions from all odd cycles, not just 3-cycles. See OPEN-Q-004.

### Impact
Any proof of Claim A for n≥6 must go beyond the per-path identity. The five proof strategies in the paper are the current best alternatives.

---

## MISTAKE-006: c₉/C(11,9) = c₇/C(11,7) "Ratio Coincidence" Has No Structural Basis

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S1 (from inbox document PALEY_T11_c9_ANALYSIS.md)
**Affects:** Any earlier estimate of c₉(T₁₁) ≈ 220 derived from this ratio

### What was assumed

The ratio c_k / C(11,k) is constant across odd k, noting c₇/C(11,7) = 1320/330 = 4 and (if so) c₉/C(11,9) = 4 → c₉ = 4 · 55 = 220.

### Why it was wrong

The ratio is NOT constant:
- k=3: c₃/C(11,3) = 55/165 = 1/3
- k=5: c₅/C(11,5) = 594/462 = 9/7 ≈ 1.286
- k=7: c₇/C(11,7) = 1320/330 = 4 ← coincidence only
- k=9: unknown

The pattern 1/3, 9/7, 4 is not monotone, not constant, has no arithmetic progression structure. The k=7 ratio is accidental. No structural reason was ever given for why the ratio at k=7 should equal the ratio at k=9.

### The correct framing

c₉(T₁₁) is determined by c₉ = (55/2)(h_QR + h_NQR) where h_QR = h({0,1}) and h_NQR = h({0,2}) are ham-cycle counts in explicit 9-vertex sub-tournaments (LEM-001). It must be computed directly.

### Impact

Any theorem or conjecture that relied on c₉ = 220 (or any ratio-derived value) should be flagged as unverified. **UPDATE (kind-pasteur-S2):** c₉ = 11055 (computed directly), and H(T_11) = 95095. CONJ-002 is fully refuted for p=11. The ratio estimate of 220 was off by a factor of 50.

---

## MISTAKE-007: Trace-Method Cycle Count Errors for c_6, c_7 in T_11

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S5 (from inbox documents more.txt, other.txt, stuff.txt)
**Affects:** Any hand computation of c_k(T_11) via tr(A^k) minus non-simple walk corrections

### What was assumed

The hand computation (inbox: more.txt) used the formula:
  k*c_k = tr(A^k) - N_k
where N_k = total non-simple closed walk contributions (Type at-v, Type A interior 3-cycle, Type B interior 4-cycle).

This gave c_6=1375 and c_7=1320.

### Why it was wrong

The non-simple walk corrections were computed incorrectly -- specifically the Type A and Type B contributions had arithmetic errors. The correct values, verified by direct DFS enumeration (other.txt), are:
  c_6(T_11) = 1595 (not 1375)
  c_7(T_11) = 3960 (not 1320)

### The correct framing

Direct DFS enumeration is more reliable than the trace correction method for large k. The correct complete cycle count table for T_11:
  c_3=55, c_4=165, c_5=594, c_6=1595, c_7=3960, c_8=7425, c_9=11055, c_10=10681, c_11=5505

These values are confirmed by the OCF identity: H(T_11) = 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155.

### Impact

The "corrected conjecture" H(T_11)=4455 (inbox: stuff.txt) was derived from the wrong c_7=1320 and was also false. The ratio 1320/330=4 was a coincidence. The actual sequence H(T_p)/|Aut(T_p)| = 1, 3, 9, 1729 for p=3,7,11 has no 3^k pattern.

### Lesson

For computing cycle counts in specific tournaments, use direct enumeration (DFS/backtracking) rather than eigenvalue-trace corrections. The trace method is valid in principle but requires extremely careful non-simple walk accounting. The LEM-001 formula c_9 = (55/2)(h_QR+h_NQR) is the correct approach (kind-pasteur-S2).

---

## MISTAKE-008: Even-Odd Split Claimed "Equivalent to OCF"

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S8 (reviewing opus-S4 output)
**Affects:** even-odd-split-lemma.md, OPEN-Q-009, TANGENTS T040

### What was assumed

The even-odd split (sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R)) was claimed EQUIVALENT to OCF in multiple places.

### Why it was wrong

The odd-S sum of Delta(S,R) = L_j(S)*R_i(R) - L_i(S)*R_j(R) is NOT the same as the cycle formula [g(S)-l(S)]*H(R). Specifically:
- L_j(S) = sum_a h_end(S,a)*T[a][j] does not include the T[i][first] factor
- g(S) = sum_{perms} T[i][first]*...*T[last][j] does include it

So even=odd gives delta_H = 2*(odd sum of Delta), but this odd sum is not the cycle formula. The even-odd split is a CONSEQUENCE of OCF, not equivalent to it. Proving even=odd would not prove OCF.

### The correct framing

The even-odd split is a necessary condition for OCF. It is a valid structural observation (verified n=5,...,8) but strictly weaker than OCF. To prove OCF via this route, one would additionally need to show that the odd-S sum of Delta(S,R) equals the cycle formula — which is essentially the full OCF identity.

### Impact

The even-odd split cannot serve as a standalone reformulation of OCF. It may still be useful as a sanity check or structural constraint, but claims of "equivalence" in the documentation have been corrected.

---

## MISTAKE-009: sympy_proof_n8.py Used Simplified n<=7 Formula

**Date discovered:** 2026-03-05
**Found by:** kind-pasteur-2026-03-05-S8
**Affects:** 03-artifacts/code/sympy_proof_n8.py (original version)

### What was assumed

The formula -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7) applies at n=8.

### Why it was wrong

THM-013 explicitly states this formula FAILS at n=8 due to VD 3-5 cycle pairs. At n=8, the complement of a 5-cycle has 3 vertices, and H(3-vertex) can be 1 or 3 (not always 1 as at n<=7). The simplified formula doesn't weight 5-cycle contributions by their complement's H value.

### The correct framing

At n=8, must use the full A-clique formula: delta_I = 2*sum_C [gained-lost]*H(comp(C)). The script has been rewritten to use this formula.

### Impact

The original script would have produced WRONG results. Fixed immediately upon discovery.

---

## MISTAKE-010: Hereditary Maximizer Chain Claimed for ALL Maximizers at Odd n

**Date discovered:** 2026-03-06
**Found by:** kind-pasteur-2026-03-06-S18g
**Affects:** INV-044, T104 (hereditary maximizer chain)

### What was assumed

"At odd n, every vertex deletion from ANY maximizer gives the (n-1)-maximizer."

### Why it was wrong

Exhaustive check at n=5 shows: 64 maximizers (H=15), of which only 24 (regular score (2,2,2,2,2)) are hereditary. The remaining 40 (score (1,2,2,2,3)) have del_hs including H=3 values, NOT max H(4)=5.

Full data:
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary (all non-regular)
- n=5: 24/64 hereditary (only regular ones)
- n=6: 0/480 hereditary (all non-regular)
- n=7: 240/240 hereditary (all regular)

### The correct framing

**Only REGULAR maximizers at odd n are hereditary.** At n=3,5,7, the regular maximizers (those with score (k,k,...,k)) have ALL vertex deletions giving max H(n-1). Non-regular maximizers at n=5 do NOT have this property.

The pattern: regular maximizers exist at odd n and are vertex-transitive, so all deletions are isomorphic. One deletion being optimal implies all are.

### Impact

The investigation INV-044 and tangent T104 need correction. The hereditary chain is: T_7 -> T_6_max -> (chain breaks at n=5 for non-regular) or T_7 -> T_6_max -> T_5_regular -> T_4_max -> T_3. The chain from Paley T_7 goes through regular maximizers at each odd step.

---

## MISTAKE-011a: Transfer Matrix M = [[1,0],[0,-1]] Claim Was Wrong

**Date discovered:** 2026-03-06
**Found by:** opus-2026-03-06-S6
**Affects:** `04-computation/transfer_matrix_test.py`, INV-001 documentation

### What was assumed

The 2x2 transfer matrix M indexed by endpoints {i,j} of a fixed arc always equals [[+1, 0], [0, -1]]. This was stated as a verified fact in `transfer_matrix_test.py`.

### Why it was wrong

Exhaustive check at n=4 (64 tournaments x 5 arc pairs = 320 tests) shows 2199/2500 failures at n=4 with the original test parameters. Example: the transitive tournament on 4 vertices gives M[0,1] = 1, not 0. The values of M[a,b] vary widely depending on the tournament (observed values include -3, -2, -1, 0, 1, 2, 3 at n=5).

The original test used random tournaments with seed=42 and checked only specific arc pairs. The claim was never properly validated.

### The correct framing

The transfer matrix M[a,b] = sum_{S ⊆ V\{a,b}} (-1)^|S| E_a(S) B_b(V\{a,b}\S) satisfies:
1. **Symmetry**: M[a,b] = M[b,a] for all a,b (PROVED symbolically n<=7, verified numerically n<=8)
2. **Trace formula**: sum_a M[a,a] = H(T) for odd n, 0 for even n (THM-027)
3. **Off-diagonal sum**: sum_{a!=b} M[a,b] = 0 for odd n, 2*H(T) for even n

The individual entries M[a,b] are NOT always 0 or ±1. They are integers whose distribution depends on the tournament structure.

### Impact

INV-001 (proving transfer matrix symmetry) is still the correct goal — the symmetry M[a,b] = M[b,a] IS always true. But the diagonal values are NOT fixed constants, which changes the picture for any proof attempt that assumed diag(1,-1).

---

## MISTAKE-011b: Paley "tournament" at p = 1 mod 4 is NOT a tournament

**Date discovered:** 2026-03-06 (kind-pasteur-S25d)
**Found by:** kind-pasteur
**Affects:** nonham_vanish_n13_paley.py, palindromic_N_proof.py (n=13 and n=17 entries), paley_super_uniform.py, any script using QR mod p for p = 1 mod 4

### What was assumed

Computations labeled "Paley T_13" and "Paley T_17" used the quadratic residue set QR mod p as the circulant generator set. This was assumed to give a tournament.

### Why it was wrong

Paley tournaments exist ONLY at p = 3 mod 4 (where -1 is a quadratic NON-residue). At p = 1 mod 4, -1 is a QR, so QR is closed under negation: d in QR implies -d in QR. This means both T[i,j]=1 and T[j,i]=1 for QR-separated pairs, giving bidirectional edges. The resulting digraph is NOT a tournament.

Valid Paley tournament primes: 3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, ...
INVALID: 5, 13, 17, 29, 37, 41, 53, 61, 73, ...

### The correct framing

- "H(T_13) = 1,579,968" is the path count of a non-tournament circulant digraph, not a tournament.
- "H(T_17) = 5,587,473,776" is similarly invalid.
- The "Paley T_5" with S={1,4} is also NOT a tournament (p=5 = 1 mod 4).
- THM-052 proof is unaffected (algebraic, works for any circulant tournament).
- n=13 re-verified with S={1,2,3,4,5,6} (valid tournament): H=3,711,175, palindromic N confirmed.

### Impact

- THM-052 verification examples at n=13 and n=17 need correction: use valid circulant tournament generator sets, not QR mod p for p=1 mod 4.
- The earlier NONHAM vanishing verification at n=13 (nonham_vanish_n13_paley.py) was on a non-tournament. Needs re-running with valid tournament.
- No impact on the algebraic proof of THM-052.
- The memory entry "H(T_p) = 3, 189, 95095" is still correct (those are p=3 mod 4 primes).

---

## MISTAKE-012: Blue Pair / Tournament Complement Confusion

**Date discovered:** 2026-03-06 (opus-S11b continued²)
**Found by:** opus
**Affects:** blue_skeleton_even_r_synthesis.py (opus-S26)

### What was assumed

The script blue_skeleton_even_r_synthesis.py claimed that "BLUE PAIR HAS SAME TRANSFER MATRIX at odd n", conflating the tiling blue pair (flipping only non-path arcs) with the tournament complement (flipping ALL arcs).

### Why it was wrong

1. The tiling blue pair only flips arcs (a,b) with |a-b| >= 2, keeping path arcs (i,i-1) fixed. So s → -s only for non-path pairs, NOT all pairs. M(T') ≠ M(T) in general.

2. The tournament complement T^c (A^c[i][j] = 1-A[i][j]) does flip ALL arcs, giving s → -s everywhere. But M(T^c) ≠ M(T) in general either!

3. The correct formula: M(T^c) = diag(M(T)) - offdiag(M(T)). The diagonal is preserved and off-diagonal is NEGATED. This is because M[a,b] for a≠b involves n-2 edge weights (odd s-degree at odd n), while M[a,a] involves n-1 edge weights (even s-degree at odd n).

### The correct framing

**THM (verified n=5 exhaustive, n=7 spot):** For any tournament T at odd n:
  M(T^c)[a,a] = M(T)[a,a]   and   M(T^c)[a,b] = -M(T)[a,b] for a≠b.

COROLLARY: M(T^c) = M(T) iff M(T) is diagonal (scalar M at odd n).

### Impact

- The claim in blue_skeleton_even_r_synthesis.py is WRONG. M(T) ≠ M(T') for tiling blue pairs, and M(T) ≠ M(T^c) for tournament complement (unless M is diagonal).
- The analysis of which iso classes are "self-paired" in that script is unaffected (it's about iso class mapping, not about M equality).
- The complement formula is a NEW theorem that should be recorded.

---

## MISTAKE-013: "All VT tournaments are self-converse" assumption is FALSE

**Date discovered:** 2026-03-06 (kind-pasteur-S25e)
**Found by:** kind-pasteur (confirmed by independent background agent + McKay database)
**Affects:** opus-2026-03-06-S26 THM-052 extension claim (commit 8ed2516)

### What was assumed

opus-S26 claimed THM-052 extends to ALL vertex-transitive tournaments via the reflection-reversal bijection phi = sigma * rev, stating "for any vertex-transitive tournament T, there exists an anti-automorphism tau." This implicitly assumes all VT tournaments are self-converse.

### Why it was wrong

At n=21, ALL 22 non-circulant VT tournaments (from McKay's database) are NOT self-converse. These are Cayley tournaments on F_21 = Z/7 x| Z/3 (Frobenius group, smallest non-abelian group of odd order) with non-normal connection sets. Exhaustive backtracking search confirms no anti-automorphism exists for any of them. All 88 circulant tournaments at n=21 ARE self-converse.

The inversion map g -> g^{-1} gives an anti-automorphism iff S is a union of conjugacy classes ("normal" S). For non-abelian groups, non-normal S creates VT tournaments without any anti-automorphism.

### The correct framing

- THM-052 is PROVED for circulant tournaments (all abelian Cayley tournaments at odd n)
- THM-052 covers all VT tournaments at n <= 19 (all circulant by McKay data)
- THM-052 **FAILS** for non-abelian-group VT tournaments at n=21 (PROVED by computation)
- M[0,1] = 45,478,409 for the F_21 non-normal tournament (H = 123,522,430,238,361)
- N(0,1,j) is NOT palindromic: all 20 values N[j] != N[19-j]

### Impact

- opus's general VT theorem claim is FALSE and must be retracted
- THM-052 must specify "circulant" or "self-converse VT" as hypothesis
- Self-converse is the RIGHT boundary: SC VT => palindromic N => scalar M; non-SC VT => non-palindromic N => non-scalar M
- The per-orbit structure (N depends on vertex-pair orbit) still holds, but palindromicity requires SC

---

## MISTAKE-014: THM-052 extension to all VT tournaments is DISPROVED

**Date discovered:** 2026-03-06 (kind-pasteur-S25e)
**Found by:** kind-pasteur
**Affects:** THM-052 scope, opus-S26 proof

### What was assumed

THM-052 was claimed to hold for ALL vertex-transitive tournaments at odd n: M = (H/n)*I.

### Why it was wrong

Computational proof at n=21: the Cayley tournament on F_21 with non-normal connection set has M[0,1] = 45,478,409 != 0. The N(0,1,j) sequence is NOT palindromic:
  N[0] = 581,223,220,317 vs N[19] = 581,314,958,778 (differ by 91,738,461)
  Alternating sum = 45,478,409

### The correct framing

THM-052 holds for self-converse VT tournaments at odd n. Self-converse is necessary (not just sufficient) for palindromic N, and palindromic N at odd n is equivalent to scalar M for VT tournaments.

Hierarchy of scalar M:
1. Circulant (always SC) => scalar M [PROVED]
2. Abelian Cayley (always SC via -x) => scalar M [PROVED]
3. Non-abelian Cayley with normal S (SC via inversion) => scalar M [PROVED by same argument]
4. Non-abelian Cayley with non-normal S (NOT SC) => M NOT scalar [DISPROVED at n=21]

### Impact

- THM-052 scope must be restricted
- Opens question: which specific non-abelian VT tournaments have scalar M?
- Answer: exactly the self-converse ones (those with normal connection sets)

## MISTAKE-015: THM-055 coefficient table has wrong c_6 and c_0 at n=7

**Date discovered:** 2026-03-06 (opus-S29, independently confirmed by kind-pasteur-S25g)
**Found by:** opus + kind-pasteur
**Affects:** THM-055 coefficient table, c_0 formula

### What was assumed

THM-055's coefficient table at n=7 claimed:
- tr(c_6) = 720 = (n-1)! (universal)
- tr(c_0) = H - 6*bc - 3*t_5 + 249/4

### Why it was wrong

tr(c_{n-1}) = sum_P e_0(s_P) = sum_P 1 = n!, NOT (n-1)!. Direct polynomial fitting of W(r) = sum_P prod(r+s_i) confirms the leading coefficient is 5040 = 7! at n=7.

The verification script `trc2_exact_formula.py` actually produces max error c_0 = 67.5, NOT 0.0. The error was visible in the output (`c0= 1.8/ 69.2`) but was not caught. The root cause: c_0 was derived from `H - c2/4 - c4/16 - 720/64` using the WRONG c_6 value, hiding the error.

### The correct framing

The correct W(r) = tr(M(r)) coefficients at n=7 are:
- w_6 = 5040 = n! (universal)
- w_4 = 240*t_3 - 2100 (unchanged)
- w_2 = -60*t_3 + 12*t_5 + 24*bc + 231 (unchanged)
- w_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc - 17/4

### Impact

- THM-055 coefficient table corrected: c_6 = n! and c_0 constant = -17/4 (not 249/4)
- The c_0 formula constant changed from +253/4 to -17/4 (difference = 270/4 = (5040-720)/64)
- The n=9 claim c_8 = 362880 = 9! = n! was already correct
- All middle coefficients (c_2, c_4) are unaffected

---

## MISTAKE-016: THM-059 recurrence had j^2 instead of (j+1)^2

**Date discovered:** 2026-03-07
**Found by:** opus-2026-03-07-S31
**Affects:** THM-059 central factorial recurrence statement (table and formulas were correct, only the recurrence formula was wrong)

### What was stated
b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}

### Why it was wrong
Plugging in: b_{2,1} should be 5, but j^2 * b_{k-1,j} = 1^2 * 1 = 1, giving b_{2,1} = 1+1 = 2 (not 5).

### The correct formula
b_{k,j} = b_{k-1,j-1} + (j+1)^2 * b_{k-1,j}

This was confirmed by checking all 15 entries of the b-triangle for k=0..4. The correct recurrence is equivalent to the standard central factorial number recurrence with shifted column indices.

### Resolution
THM-059 corrected. The (j+1)^2 factor now has a combinatorial explanation via the Eulerian polynomial decomposition: F_f(r) = sum_k A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k, where the central factorial structure emerges from expanding in u = (r+1/2)(r-1/2) = r^2-1/4.

### Impact
- The numerical table and all computed F_j values were always correct
- Only the stated recurrence formula was wrong
- The OEIS A036969 identification may need clarification (different column conventions)

---

## MISTAKE-010: Missing 2^s Factor in M[a,b] Walsh Formula (THM-080)

**Date discovered:** 2026-03-07 (S35c7)
**Found by:** opus-2026-03-07-S35c7
**Affects:** THM-080 Walsh formula for M[a,b]

### What was assumed
The Walsh coefficient hat{M[a,b]}[S] = (-1)^{asc(S)} * (n-2-d)!/2^{n-2}, with NO dependence on the number of components. This was described as a "fundamental simplification" compared to H.

### Why it was wrong
The formula was verified exhaustively only at n=5, where ALL valid monomials have s=0 (zero unrooted components). At n=5, the maximum degree is 3, and with only 3 interior vertices, there's no room for unrooted even-length components to coexist with rooted ones. So the 2^s factor was always 1, making it invisible.

At n=7, degree-3 monomials like P1(a-rooted) + P2(unrooted) have s=1, and the formula without 2^s gives wrong reconstruction (16/20 failures).

### The correct formula
hat{M[a,b]}[S] = (-1)^{asc(S)} * **2^s** * (n-2-d)!/2^{n-2}

where s = number of unrooted (even-length) components. Each unrooted component has 2 valid orientations in the HP (both giving the same chi_S sign), contributing a factor of 2. Rooted components have only 1 valid orientation (pinned at a or b).

### Impact
- THM-080 formula corrected with 2^s factor
- Walsh proof of M[a,b]=M[b,a] symmetry still holds (2^s is symmetric in a,b)
- H-M comparison now shows PARALLEL structure: H has 2^r (all components unrooted), M has 2^s (only unrooted components contribute orientations)
- The "no r-dependence" claim was wrong; M DOES depend on component structure via s
- n=7 reconstruction: 20/20 match with corrected formula

### Lesson
Always verify formulas at the NEXT size up before claiming generality. n=5 was too small to expose the s-dependence.
