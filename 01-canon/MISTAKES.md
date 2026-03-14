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

**RESOLVED** (DISC-001, closed by kind-pasteur-2026-03-05-S3): Independent verification using tournament_lib.py confirms the paper's results are valid. See 02-court/resolved/DISC-001-mu-bug-vs-verification.md.

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

---

## MISTAKE-017: "Non-Paley DRT at n=11" from invalid tournament connection set

**Date discovered:** 2026-03-07
**Found by:** kind-pasteur-2026-03-07-S39b
**Affects:** INV-068, MEMORY.md DRT analysis section, TANGENTS.md DRT entry

### What was assumed
A "non-Paley DRT at n=11" was constructed using connection set {1,2,3,5,8} in Z_11 (circulant digraph). Claims: c3=44, c5=407, H=69311, |Aut|=11. "Paley strictly dominates in ALL cycle counts."

### Why it was wrong
The connection set {1,2,3,5,8} does NOT give a tournament. For a circulant tournament on Z_p, the connection set S must satisfy S ∩ (-S) = ∅ (so each pair {i,j} has exactly one directed arc). But {1,2,3,5,8} contains BOTH 3 and 8=11-3, and BOTH 1 and 10=11-1... wait, 10 is NOT in S. Let me re-check: -S = {11-s : s ∈ S} = {10, 9, 8, 6, 3}. S ∩ (-S) = {3, 8} ≠ ∅.

So for any pair (i,j) where (j-i)%11 ∈ {3, 8}: BOTH T[i][j]=1 AND T[j][i]=1. The resulting digraph has bidirectional edges and is NOT a tournament. All computations (c3, c5, H, is_doubly_regular) were performed on a non-tournament digraph and are MEANINGLESS.

### The correct framing
An exhaustive search of all 32 valid tournament connection sets in Z_11 (choosing one from each pair (d, 11-d)) found exactly 2 that are (11,5,2)-difference sets: {1,3,4,5,9} (QR) and {2,6,7,8,10} (NQR). These give ISOMORPHIC tournaments (both Paley T_11). There is NO non-Paley circulant DRT at n=11.

Whether a non-circulant DRT exists at n=11 remains an open question. At prime order p, all groups are Z_p, so all Cayley tournaments are circulant. A non-circulant DRT would need a different construction.

### Impact
- ALL claims about "non-Paley DRT at n=11" are INVALID
- INV-068 "Paley dominance" finding needs complete re-evaluation
- The claimed c3=44 was wrong — Moon's formula gives c3=55 for ALL regular n=11 tournaments
- The claimed "Paley strictly dominates in all cycle counts" is unverifiable since no valid comparison tournament exists
- MEMORY.md entry on DRT analysis at n=11 needs correction

### Lesson
When constructing a circulant tournament from a connection set S ⊂ Z_p^*, ALWAYS verify S ∩ (-S mod p) = ∅. A (v,k,λ)-difference set is NOT automatically a valid tournament connection set.

---

## MISTAKE-016: Wrong formula for ker(d_2^rel) in relative homology

**Date discovered:** 2026-03-08 (kind-pasteur-S41)
**Found by:** kind-pasteur-S41, via manual computation contradicting script output
**Affects:** beta2_relative_homology.py, beta2_relative_correct.py; HYP-213 verification

### What was assumed
The script `beta2_relative_homology.py` computed ker(∂_2^rel) as:
  `(ker ∂_2 + V_2) / V_2`
where V_2 = Ω_2(T\v) (non-v subspace of Ω_2).

### Why it was wrong
The correct formula for ker(∂_2^rel) in the quotient complex Ω_*/V_* is:
  `∂_2^{-1}(V_1) / V_2`
where V_1 = Ω_1(T\v) (non-v arcs). This is the preimage of V_1 modulo V_2.

The wrong formula misses elements x ∈ Ω_2 whose boundary ∂_2(x) is NONZERO but lies entirely in V_1 (non-v arcs). Such elements are relative 2-cycles but NOT absolute 2-cycles.

Concretely: P_v ∘ ∂_2(x) = 0 (projection of boundary onto v-arcs vanishes), but ∂_2(x) ≠ 0.

### The correct framing
ker(∂_2^rel) = dim(Ω_2) - rk(M) - dim(V_2), where M = ∂_2|_{Ω_2} restricted to rows of v-arcs.

This correctly counts the preimage of V_1 in Ω_2.

### Impact
- **HYP-213 is REFUTED**: H_2(T, T\v) > 0 for many (T,v) pairs at n ≥ 4.
  - n=4: 16/256 pairs (6.25%)
  - n=5: 840/5120 pairs (16.4%)
  - n=6: 35,328/196,608 pairs (18%)
- The proposed inductive proof of β_2 = 0 via H_2(T,T\v) = 0 does NOT work.
- However, β_2 = 0 itself is NOT affected — it remains computationally verified.
- The connecting map δ: H_2(T,T\v) → H_1(T\v) is always injective (verified n=4,5), consistent with β_2 = 0 via the long exact sequence.

### Lesson
When computing relative homology H_*(X, A) via quotient complexes:
1. ker(∂_p^{rel}) is NOT (ker ∂_p + C_*(A)) / C_*(A).
2. ker(∂_p^{rel}) = ∂_p^{-1}(C_{p-1}(A)) / C_p(A).
3. These differ whenever there are elements whose boundary is nonzero but lands in the sub-complex.
4. Always verify relative homology against the long exact sequence.

---

## MISTAKE-018: beta_3 <= 1 Assumed for All Tournaments

**Date discovered:** 2026-03-09 (kind-pasteur-S48)
**Found by:** kind-pasteur-S48 via extended sampling at n=8 (5000 random tournaments)
**Affects:** THM-123 (was THM-110) proof architecture, HYP-371b, HYP-375, HYP-342, HYP-380, HYP-393 scope

### What was assumed
Multiple hypotheses and proof strategies assumed beta_3 <= 1 for ALL tournaments:
- HYP-371b: "beta_3=2 impossible"
- HYP-375: "beta_3 <= 1 at n=9"
- THM-123 proof architecture: Claims I, II, III designed to prove beta_3 <= 1
- The opus exhaustive proof at n=7 was incorrectly assumed to generalize

### Why it was wrong
beta_3 = 2 DOES occur at n=8. Four examples found in 5000 random tournaments (rate ~0.08%):
- Profile: (1, 0, 0, 2, 0, 0, 0, 0) — two independent H_3 generators
- Scores: (2,3,3,3,4,4,4,5) and (3,3,3,3,4,4,4,4) — near-regular
- Confirmed by BOTH max_p=5 and max_p=7 in full_chain_complex_modp (mod-p exact)
- All b3=2 tournaments have good vertices (b3(T\v)=0 for some v)

Previous sampling (200 at n=9, 100 at n=8) was insufficient to detect 0.08% rate.

### The correct framing
- beta_3 <= 1 is proved ONLY at n <= 7 (exhaustive, HYP-393)
- beta_3 = 2 at n=8 (confirmed, 4/5000)
- beta_3 may grow further at n >= 9

### Impact
- THM-123 proof architecture is valid ONLY at n <= 7
- Claims I (i_*-injectivity) also FAIL at n=8 (13 violations in 5000 trials, even with b4=0)
- Claim III (consecutive seesaw) FAILS at n=8 (beta_3+beta_4 coexistence)
- The beta_3 <= 1 bound is a SMALL-n PHENOMENON, not a universal property

### Lesson
1. Small sample sizes (100-200) cannot detect 0.1% phenomena. Use 5000+ for rare events.
2. Properties proved exhaustively at n<=7 do NOT automatically extend to n>=8.
3. n=8 is a critical threshold where many path homology structural properties break down:
   consecutive seesaw, i_*-injectivity, beta_3<=1, bad vertex acyclicity.

## MISTAKE-019: Int64 Overflow in Chained Numpy Matrix Multiplication

**Date discovered:** 2026-03-10 (kind-pasteur-S50)
**Found by:** kind-pasteur-S50 via comparison of two K_tv computation methods
**Affects:** opus-S59's tv_cycle_structure.py (Ghost Cycle "failures" are spurious), potentially any script using `A @ B @ C % PRIME` pattern

### What was assumed
`D3 @ (tv_omega @ ob3).T % PRIME` safely computes the matrix product mod PRIME.

### Why it was wrong
With `RANK_PRIME = 2^31 - 1`:
- `tv_omega @ ob3` produces entries up to ~4.6 × 10^18 (fits in int64, max 9.2 × 10^18)
- BUT these are NOT reduced mod PRIME — entries can be >> PRIME
- `D3 @ X.T` then involves products up to `(2^31) * (4.6e18) = 9.9e27`, massively exceeding int64 max
- Result: silent int64 overflow → wrong matrix entries → wrong rank → wrong K_tv

This caused opus-S59's tv_cycle_structure.py to report Ghost Cycle failures in 14/504 pairs at n=7 and 11/304 at n=8. ALL of these "failures" are arithmetic artifacts.

### The correct framing
ALWAYS reduce mod PRIME between chained matrix multiplications:
```python
# WRONG (can overflow):
result = A @ B @ C % PRIME

# RIGHT:
temp = A @ B % PRIME
result = temp @ C % PRIME

# BEST: use the new safe utility:
from tournament_utils import matmul_mod
result = matmul_mod(matmul_mod(A, B), C)
```

The `matmul_mod()` function in tournament_utils.py automatically chunks the inner dimension to prevent overflow, even for single multiplications with large entries.

### Impact
- Ghost Cycle (K_tv = B_tv) HOLDS universally at n ≤ 8 (0 real failures in 1000+ tests)
- HYP-408 (codim-1 universality) remains computationally verified at n ≤ 8
- No real mathematical result is affected; only the false "counterexamples" are invalidated

### Lesson
1. NEVER chain numpy `@` without intermediate `% PRIME` when PRIME ≈ 2^31
2. Use `matmul_mod()` from tournament_utils.py for all modular matrix arithmetic
3. When two equivalent computations disagree, suspect numerical issues before mathematical failure

## MISTAKE-020: Truncated Chain Complex Gives False Betti Numbers at Top Degree

**Date discovered:** 2026-03-10
**Found by:** kind-pasteur-S50

### What was assumed
Using `full_chain_complex_modp(A, n, max_p=6)` for n=8 tournaments, opus-S59 reported β_6 nonzero for 89.8% of tournaments (HYP-420), with values ranging 1-25.

### Why it was wrong
With `max_p=6`, the computation gives β_6 = ker(d_6) - ranks.get(7, 0). Since degree 7 is not computed, `ranks.get(7, 0)` returns 0. The reported "β_6" is actually just dim(ker d_6), NOT the true Betti number.

With `max_p=7` (full complex): d_7 is injective on Omega_7, and rk(d_7) = ker(d_6) EXACTLY for all 50 tested tournaments. True β_6 = 0 always.

### The correct framing
The Betti number at the highest computed degree is always an UPPER BOUND (missing the image from the next degree). For n-vertex tournaments, always use `max_p=n-1` to get correct Betti numbers, especially at degrees n-2 and n-1.

Correct results: β_{n-1} = β_{n-2} = 0 for ALL tournaments at n=3-8 (HYP-423, HYP-424). The top boundary map d_{n-1} is always injective.

### Impact
- HYP-420 is FALSE. β_{n-2} is NOT generically nonzero at n=8.
- The "β_6 among β_3=1" distribution in opus's beta4_at_n7.out is entirely artifactual.
- All lower-degree Betti numbers (β_0 through β_5) from that computation are correct.

### Lesson
1. ALWAYS use max_p=n-1 when computing Betti numbers to avoid truncation artifacts
2. Betti at the max computed degree is an UPPER BOUND (Betti at max_deg-1 and below are exact)
3. When β at max_deg seems surprisingly large/nonzero, check if im(d_{max_deg+1}) is missing

---

## MISTAKE-019: THM-136 Sign Convention Error

**Date discovered:** 2026-03-12
**Found by:** kind-pasteur-S57
**Affects:** THM-136 formula statement (not the verbal description or proof mechanism)

### What was assumed
The trace alternation sign formula was stated as:
`sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-3)/2}`

### Why it was wrong
Direct computation at p=7,11,19,23 shows:
- k=5: Delta > 0 (positive), but (-1)^{(5-3)/2} = (-1)^1 = -1 (WRONG)
- k=7: Delta < 0 (negative), but (-1)^{(7-3)/2} = (-1)^2 = +1 (WRONG)

The formula gives the OPPOSITE sign for every k.

### The correct framing
`sign(tr(A_P^k) - tr(A_I^k)) = (-1)^{(k-1)/2}`

Equivalently: positive for k = 1 mod 4, negative for k = 3 mod 4.
Verified with 1218+ individual (k,p) tests, zero failures.

Note: the VERBAL description in THM-136 was always correct ("k=1 mod 4: Paley wins").
Only the symbolic formula was off by one power.

### Impact
- Formula in THM-136 theorem file CORRECTED
- No downstream impact: all proofs used the verbal description, not the formula
- The algebraic proof (kind-pasteur-S57) uses the correct convention throughout

## MISTAKE-021: S70 "GLMY Betti Numbers" Use Wrong Chain Complex

**Date discovered:** 2026-03-13
**Found by:** opus-2026-03-13-S71
**Affects:** ALL scripts from S70 session: betti_omega_connection.py, betti_divisibility.py, per_eigenspace_betti.py, per_eig_betti_n9.py, and all results/theorems derived from them (THM-154, eigenspace Betti uniformity)

### What was assumed
The S70 scripts computed "GLMY path homology Betti numbers" using:
- Allowed paths = "regular paths" (v_i→v_{i+1} AND v_{i-1}→v_{i+1})
- Boundary = interior-only deletion (indices 1 to m-1)
Results were called "GLMY Betti numbers" and compared with GLMY literature.

### Why it was wrong
The actual GLMY path homology uses:
- Allowed paths = directed paths (v_i→v_{i+1} only, NO skip-one requirement)
- Boundary = full vertex deletion (indices 0 to m), but restricted to Ω_m subspace
- Ω_m = {u ∈ A_m : ∂u has all components in A_{m-1}}

**These give DIFFERENT chain complexes with different Betti numbers:**
- Paley P_7 GLMY: β = [1,0,0,0,6,0,0], dim(A_2)=63
- Paley P_7 S70:  β = [7,0,0,21,21,21,21], dim(A_2)=21

The "regular path + interior boundary" complex IS a valid chain complex
(d²=0 verified), but it is NOT standard GLMY path homology.

### The correct framing
There are TWO distinct valid chain complexes for tournaments:

1. **GLMY Path Homology** (standard): directed paths, full boundary on Ω_m.
   Implemented correctly in path_homology_v2.py.
   β_0 = 1 for all tournaments. β_2 = 0 for all tested tournaments (n≤8).

2. **Tournament Regular Homology (TRH)** (novel?): regular paths, interior boundary.
   Used in S70 scripts. β_0 = n for all tournaments on n vertices.
   Has eigenspace Betti uniformity and divisibility by n for circulants.

Both are valid mathematical objects. But they should not be conflated.

### Impact
- THM-154 (Betti divisibility) applies to TRH, not GLMY
- Eigenspace Betti uniformity applies to TRH, not GLMY
- β_2=0 for all tournaments holds for BOTH (GLMY verified n≤8, TRH verified n≤8)
- The S70 "per-eigenspace Betti" results are self-consistent but not GLMY
- The S38-S41 β_2=0 results (from path_homology_v2.py) are correct GLMY
- circulant_homology.py implements yet another convention (full boundary on regular paths) which is NEITHER GLMY nor TRH

### Lesson
ALWAYS verify which chain complex you're computing. The three ingredients
(allowed paths, boundary convention, Ω subspace) must be consistent.
When reading "path homology" results, check which convention is used.

---

## MISTAKE-019: TWO bugs in independent set backtracking algorithm

**Date discovered:** 2026-03-13, kind-pasteur-S60
**Found by:** kind-pasteur
**Affects:** alpha3_p7_only.py, alpha3_moment_analysis.py, overlap_weight_analysis.py, H_energy_decomposition.py, cycle_walsh_decomposition.py, moment_cancellation_mechanism.py, overlap_gauss_bridge.py, alpha_directed_p11.py, alpha_full_p11.py, alpha2_direct_verify.py, backtrack_debug.py (ALL files with independent set enumeration)

### Bug 1: Missing vertex 0 (`backtrack(0,0,0)` should be `backtrack(-1,0,0)`)

The backtracking function `backtrack(v, mask, size)` iterates `for w in range(v+1, n)`. When called with `v=0`, the loop starts at `w=1`, SKIPPING vertex 0 entirely. This undercounts all alpha_j.

**Fix:** Call `backtrack(-1, 0, 0)` so the loop starts at `w=0`.

### Bug 2: Skipping consecutive indices (`backtrack(w+1, ...)` should be `backtrack(w, ...)`)

The recursive call `backtrack(w + 1, mask | nbr[w], size + 1)` passes `v = w+1`. At the next level, the loop starts at `w' = v+1 = w+2`, SKIPPING index `w+1`. This means any independent set containing cycles with consecutive indices is missed.

**Fix:** Change to `backtrack(w, mask | nbr[w], size + 1)`. Then the next level's loop starts at `w+1`, correctly considering all higher indices.

### Concrete example

At p=7, Interval tournament S=[1,2,3]:
- 59 directed cycles, 14 disjoint (3,3)-pairs (correct)
- Bug 2: Pair (5,6) = ({0,3,6}, {1,2,5}) has consecutive indices and was SKIPPED
- Backtracking gave alpha_2=13 instead of 14, H=171 instead of 175
- Held-Karp gives H=175 (correct)

### Impact
- All previous alpha_j values from backtracking are SUSPECT
- THM-027 alpha_2 values at p=7 need recheck (Paley alpha_2=7 was coincidentally correct because no consecutive disjoint pairs)
- Any H derived from backtracking alpha may be wrong

### Lesson
When implementing independent set enumeration via backtracking, the recursive call after selecting vertex w should pass v=w (NOT v=w+1). The `range(v+1, n)` in the next level already excludes w.

---

## MISTAKE-022: Sparse Gaussian Elimination Fill-In Bug

**Date discovered:** 2026-03-13, opus-S71c (9th context window)
**Found by:** opus, when k=0 eigenspace Betti numbers came out negative
**Affects:** p19_omega5_sparse.py, p23_omega5_sparse.py, p31_omega5_sparse.py, p43_omega5_sparse.py (ALL scripts using the sparse Gaussian pattern with single-pass row iteration)

### What was assumed
The sparse Gaussian elimination iterated over `sorted(col.keys())` once, subtracting each matching pivot. This should correctly eliminate all pivot contributions.

### Why it was wrong
When subtracting a pivot at row `r`, the pivot vector has entries at rows `r' > r` (fill-in). These new entries at rows NOT in the original column are never checked against existing pivots at those rows, because the `sorted(col.keys())` list was computed BEFORE the subtraction and doesn't include fill-in entries.

Concrete example: column has entries at rows {3, 7}. Pivot at row 3 has entries at {3, 5, 7}. After subtracting pivot 3, the column has entries at {5, 7}. Row 5 was NOT in the original sorted list, so even if there's a pivot at row 5, it's never subtracted. This causes the rank to be OVERCOUNTED (some columns that should reduce to zero don't).

### The correct framing
After any pivot subtraction, restart the row scan from the beginning (or at least from the newly-created entry). A simple fix: wrap the elimination loop in `while changed: ... break after subtraction`.

### Impact
- **P_19 Omega_5 was 12602 (WRONG), correct is 23832**
- **P_23 Omega_5 was 50715 (WRONG), correct is 78430**
- **P_31 Omega_5 was 252065 (WRONG), correct is 456330**
- **P_43 Omega_5 was 1429652 (WRONG), correct is 2865660**
- P_7 and P_11 were unaffected (small enough that fill-in didn't change rank)
- HYP-790 ("Omega_5 not polynomial in m") was based on wrong data — **RETRACTED**
- **CORRECTED**: Ω_5 = m(m-1)(m³-6m²+10m-2) — a **clean integer polynomial** in m!
- All formulas Ω_d for d ≤ 5 are now proven/verified

### Lesson
In sparse Gaussian elimination, fill-in from pivot subtraction can create new entries at rows that were not in the original column. These MUST be processed against their pivots. Always use a while loop that restarts after each subtraction, or maintain a priority queue of unprocessed rows.

---

## MISTAKE-023: α₁ Counts DIRECTED Odd Cycles, Not Vertex-Sets

**Date discovered:** 2026-03-14
**Found by:** opus-2026-03-14-S71d
**Affects:** two_and_three_universality.py, i3_mod3_proof.py, vandermonde_sigma_connection.py, jacobsthal_23_deep.py (first version), and any script computing I(CG, x) by counting cycle vertex-sets

### What was assumed

The independence polynomial I(Ω, x) was computed by enumerating odd cycle **vertex-sets** (frozenset of vertices), counting each set once regardless of how many distinct directed cycles it supports.

### Why it was wrong

The conflict graph Ω(T) has vertices = **directed odd cycles** (definition in definitions.md line 37). For 3-cycles in tournaments, each vertex triple supports at most 1 directed 3-cycle, so vertex-set counting is correct. But for 5-cycles and above, a single vertex-set can support **multiple** distinct directed cycles:

- Example: bits=40 at n=5, the 5-vertex set {0,1,2,3,4} supports **3** distinct directed 5-cycles
- Vertex-set method gives α₁=5 → I(2)=11, but H=15
- Directed-cycle method gives α₁=7 → I(2)=15 = H ✓

### The correct framing

When computing I(Ω, x):
1. For each vertex-set of size k, enumerate ALL distinct directed k-cycles (normalize by fixing start vertex and direction)
2. Each distinct directed cycle is a SEPARATE vertex of Ω(T)
3. Two vertices of Ω are adjacent iff the underlying vertex-sets intersect

**Exhaustive verification at n=5:**
- Vertex-set method: 184/1024 mismatches with H
- Directed-cycle method: 0/1024 mismatches with H

### Impact

- All α₁ values from scripts using vertex-set counting at n≥5 are WRONG (undercounted)
- The Vandermonde extraction results (HYP-867, HYP-868) were based on the wrong α values
- The 3/2 ratio result may still hold (it was measured within lambda fibers, not from α directly)
- The structural insights about 7→8 transition are UNAFFECTED (vertex-set counting is correct for α₂ when cycles have different sizes)
- Scripts need to be updated to use directed cycle enumeration

### Lesson

The definition says "vertices are **directed** odd cycles." For 3-cycles in tournaments, vertex-set = directed cycle (1-to-1). For k≥5 cycles, a k-vertex tournament subtournament can have multiple Hamiltonian cycles. Always enumerate directed cycles explicitly.
