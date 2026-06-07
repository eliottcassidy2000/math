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

## MISTAKE-013b: Missing 2^s Factor in M[a,b] Walsh Formula (THM-080)

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

## MISTAKE-016b: Wrong formula for ker(d_2^rel) in relative homology

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

## MISTAKE-019b: THM-136 Sign Convention Error

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

## MISTAKE-019c: TWO bugs in independent set backtracking algorithm

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

---

## MISTAKE-024: H=63 Falsely Claimed Permanently Forbidden

**Date discovered:** 2026-03-14
**Found by:** opus-2026-03-14-S71h (cross-referencing S86 broadcast with prior results)
**Affects:** HYP-1303, MSG-218 (S86 broadcast), MSG-139 (S86 to kind-pasteur)

### What was assumed

S86 claimed: "H=63 FORBIDDEN: 63=7×9=I(K₃,2)×I(2K₁,2). Requires K₃ component in Ω, impossible by THM-201." This was marked CONFIRMED as HYP-1303.

### Why it was wrong

The argument only blocks DISCONNECTED conflict graphs where Ω = K₃ ⊔ 2K₁. But Ω can be a CONNECTED graph with I(Ω,2)=63. Multiple prior sessions had already established:
- S65-c (MSG-084): "H=63,107,119,149 (the n=7-specific gaps) are ALL achieved at n=8"
- S71f (MSG-197): "63 achievable at n≥8"
- S71g (MSG-201): "H=63 found at n=8 (27/100k)"
- hspectrum_density.out: "63 = 7*9 -- ACHIEVED at n=8 (not permanent)"

The S86 session re-derived an incomplete argument without checking these earlier verified results.

### The correct framing

H=63 is NOT permanently forbidden. It is a temporary gap at n=7 (like 107, 119, 149) that IS achieved at n=8. The ONLY permanently forbidden H values proved so far are {7, 21}.

The disconnected decomposition 63=7×9 is correctly blocked, but connected Ω graphs with I=63 exist and can apparently be realized as Ω(T) at n=8.

### Impact

- HYP-1303 changed to REFUTED
- HYP-1295 changed to REFUTED
- The S86 broadcast claim "all three known forbidden values {7,21,63} are now explained" is WRONG
- Only {7, 21} are explained as permanently forbidden

### Lesson

Always check prior session results before claiming a new proof. The H-spectrum density analysis (hspectrum_density.out) had already settled this question computationally. Cross-reference before broadcasting claims.

---

## MISTAKE-025: S112 W(8) Value Off By 8

**Date discovered:** 2026-03-16 (opus-S90 continuation session)
**Found by:** opus-S90, via independent brute-force verification
**Affects:** kind-pasteur-S112 W(n) sequence, D_n(2) computations

### What was claimed
kind-pasteur-S112 reported W(n) = 1, 2, 8, 32, 158, 928, 6350, **49760**, 439766 for n=1..9.

### Why it was wrong
Independent brute-force computation (iterating over all 8! = 40320 permutations) gives W(8) = **49752**, not 49760. The error is exactly +8 = 2³. The S89c C-program DP computation also gives 49752, confirming the brute-force result.

### The correct values
W(n) = 1, 2, 8, 32, 158, 928, 6350, **49752**, 439670 (from S89c DP).

### Impact assessment
- **S89c values (opus):** CORRECT through n=27 (computed by bitmask DP in C).
- **S112 values (kind-pasteur):** INCORRECT at n=8 by +8. Values at n≤7 match. Values at n≥9 need reverification against S89c.
- **OEIS submission:** Use S89c values, not S112.
- **CV² formula and g_k polynomials:** UNAFFECTED (derived independently of W(n) enumeration).

### Source
opus-2026-03-16 (S90 continuation): brute-force W(8) verification via Python permutation enumeration.

### Lesson
When two independent computations disagree, trust the one with the simpler algorithm (brute force) over the one with more complex logic. The discrepancy of exactly 2³ suggests a boundary condition or off-by-one error in the S112 computation, not a random bug.

---

## MISTAKE-026: Cross-Ratio of Cayley Orbit Initially Claimed as 8/7

**Date discovered:** 2026-03-15 (code review during S90)
**Found by:** opus-S90 code review agent
**Affects:** monad_cayley_s90c.py, commit messages

### What was claimed
The cross-ratio of the Q-orbit of x=2 was initially computed as 8/7, using the WRONG orbit point (3 instead of -3).

### Why it was wrong
Q(2) = (1+2)/(1-2) = 3/(-1) = **-3**, not 3. The orbit is {2, **-3**, -1/2, 1/3}. The cross-ratio CR(2, -3, -1/2, 1/3) = **2**, not 8/7.

### The correct value
Cross-ratio = 2 = the OCF fugacity itself. This is MORE meaningful than 8/7 — the cross-ratio equals the evaluation point.

### Impact
- The narrative about "tournament constant 8/7" in commit messages is wrong.
- The correct "tournament constant" is 2 (the fugacity).
- Script monad_cayley_s90c.py has been corrected.

### Source
Code review agent, opus-2026-03-15-S90.

---

## MISTAKE-018b: THM-225 "Universal Top Eigenvalue = n" is FALSE at n ≥ 9

**NOTE:** This was originally numbered MISTAKE-018 from a different branch, causing a collision with MISTAKE-018 (beta_3 <= 1). Renumbered to 018b by opus-2026-04-01-S1.

**Date discovered:** 2026-03-15
**Found by:** opus-S72d
**Affects:** THM-225, HYP-1594

### What was assumed

That the top eigenvalue of C_T^TC_T equals n for ALL tournaments on n vertices (verified exhaustively at n=5, sampled at n=6). This was stated as a theorem.

### Why it was wrong

The proof strategy required rank(C_R) < r = (n-1)(n-2)/2. At n ≤ 8, this holds because max c₃ < r (the number of cyclic triples never exceeds the rank of the full constraint matrix). At n = 9, c₃ can reach 30 while r = 28, and for ~0.1% of tournaments, rank(C_R) achieves its maximum r, leaving ker(C_R) ∩ im(C^T) = {0}. The top eigenvalue then drops to ~8.84-8.94.

### The correct framing

THM-225 holds for n ≤ 8 (PROVED for n ≤ 5 exhaustive, sampled at n=6,7,8 with 0 violations from 20000 samples each). It FAILS at n ≥ 9. The condition for top eigenvalue = n is rank(C_R) < (n-1)(n-2)/2.

### Impact

The spectral T/R duality (C_T^TC_T + C_R^TC_R = n·P) and the 3/n bridge framework remain valid. The universal top eigenvalue was a COROLLARY that holds only when the cyclic boundaries don't span the full constraint space.

### Lesson

When verifying at n=5 and n=6, the parameter regime (c₃ < r always) hid the failure mode. Always check at the CROSSOVER point where qualitative behavior changes. For rank arguments, the critical n is where max c₃ first exceeds r.

---

## MISTAKE-027: THM-080 Amplitude Table Wrong at n=9

**Date discovered:** 2026-03-16
**Found by:** opus-2026-03-16-S73
**Affects:** THM-080 amplitude table (lines 156-161), not the formula itself

### What was assumed

The amplitude table in THM-080 listed n=9 entries as: deg 1 (s=0) = 3/2, deg 3 (s=0) = 3/8, deg 3 (s=1) = 3/4, deg 5 (s=0) = 1/16, deg 7 = 1/128.

### Why it was wrong

The stated formula is |hat{M}[S]| = 2^s × (n-2-|S|)!/2^{n-2}. At n=9:
- d=1, s=0: formula gives 6!/128 = **45/8**, not 3/2
- d=3, s=0: formula gives 4!/128 = **3/16**, not 3/8
- d=3, s=1: formula gives 2×3/16 = **3/8**, not 3/4
- d=5, s=0: formula gives 2!/128 = **1/64**, not 1/16
- d=7, s=0: formula gives 0!/128 = 1/128 ✓ (only correct entry)

The formula works perfectly at n=3, 5, 7 — only n=9 has errors. The d=3 and d=5 wrong values are the CORRECT formula values but with s shifted up by 1 or 2 (unrooted component miscount). The d=1 entry (3/2) doesn't correspond to any valid s value (45/8 × 2^s ≠ 3/2 for any integer s).

### The correct framing

The formula is correct. The table had a transcription error at n=9 only. The n=9 verification was "partial" and apparently didn't catch the table/formula mismatch. Corrected amplitude table is in THM-080.

### Impact

Low — the formula itself is correct and was verified computationally at n=5 (exhaustive), n=6 (exhaustive), n=7 (20/20). Only the summary table for n=9 was wrong. No downstream results depend on the specific n=9 table values.

### Lesson

**This is MISTAKE-013b (the original missing 2^s) echoing forward.** The 2^s correction was caught and fixed at n=7, but the n=9 table values were apparently populated from a pre-correction computation (or from hand calculation that repeated the component-counting error at higher n). Always re-derive table entries from the corrected formula rather than carrying forward values from partial computation.

This is also a meta-lesson about the amplitude table itself: it was the only place in THM-080 where specific numerical values were stated without being individually verified against the formula. The formula (analytically proved) was trustworthy; the table (hand-calculated) was not.

---

## MISTAKE-028: Mersenne / k-nacci Numbers Falsely Claimed to Control Forbidden H Values

**Date discovered:** 2026-03-17
**Found by:** opus-2026-03-17-S74 (forbidden values audit)
**Affects:** casual-writeup.md, formal-writeup.md, substack-hooks.md (Hook N), HYP-1600, HYP-1618 (original), HYP-1623, HYP-1624, riemann-zeta-tournament.md, multiple results files

### What was assumed

Multiple writeups and hypotheses claimed: "The k-nacci Mersenne identity connects forbidden H values (7 = 2^3 - 1, 31 = 2^5 - 1, 127 = 2^7 - 1) to Mersenne primes via k-nacci transfer matrices." The original HYP-1618 claimed "ζ(-3) = 7" (standard Riemann zeta). Various scripts called 31 "forbidden."

### Why it was wrong

1. **31 is achievable** at n=6 (tournament bits=146, verified exhaustively).
2. **63 is achievable** at n=8 (already documented in MISTAKE-024).
3. **127 is achievable** at n=7.
4. The standard Riemann ζ(-3) = 1/120, NOT 7.
5. The tribonacci trace sequence [1, 3, 7, 11, 21, 39, 71, 131, ...] contains both forbidden values (7, 21) AND achievable values (1, 3, 11, 39, 71, 131, ...). The k-nacci trace hitting 7 and 21 is a numerical coincidence, not a causal mechanism.

### The correct framing

**Only H=7 and H=21 are permanently forbidden** (proved: THM-029, THM-079). The actual mechanisms are:
- H=7: 3 pairwise-conflicting cycles always force additional cycles (THM-029)
- H=21: all OCF decompositions blocked by tournament forcing (THM-079, 464-line proof)

Best characterizations of {7, 21}:
- {7·3⁰, 7·3¹}: the 7-obstruction has nilpotency 2 (HYP-1231). 7·3² = 63 is achievable.
- {Φ₃(2), Φ₃(4)}: third cyclotomic polynomial at even args (HYP-1317). Φ₃(6) = 43 is achievable.
- Both have I-polynomials factoring through I(K₃, x) = (1+3x) (HYP-1315).

THM-227 (k-nacci Mersenne) is a valid theorem about transfer matrices. It just doesn't characterize forbidden H values.

### Impact

Medium — the false claim propagated through 6+ files across multiple sessions. All have been corrected. No theorems or proofs are affected (the actual forbidden value proofs THM-029 and THM-079 are correct and don't use the Mersenne connection).

### Lesson

**One data point is not a pattern.** The entire false extrapolation rested on a single observation: 7 = 2³ - 1 is both a Mersenne number and forbidden. From this, it was incorrectly inferred that other Mersenne numbers (31, 127) are also forbidden. A simple check (is H=31 achieved at n=6?) would have caught this immediately.

This is a variant of MISTAKE-024 (H=63 falsely claimed forbidden) — the same class of error, just with a different numerological motivation. The meta-lesson: when claiming a numerical pattern "explains" something, verify it at the NEXT case before asserting it as a principle.

---

## MISTAKE-029: Formula E = (T - D)/2 is WRONG for the meta-graph edge count

**Date discovered:** 2026-03-23
**Found by:** opus-2026-03-23-S211
**Affects:** degeneracy_second_moment_s210.py, all claims that E = (3T - S_2)/4

### What was assumed
The meta-graph edge count formula E = (T - D)/2 was claimed in S210, where T = total arc-orbits and D = sum C(mult(C→D), 2) is the degeneracy. The derived formula E = (3T - S_2)/4 was presented as the "keystone" edge count formula, and the reverse-engineered S_2 sequence was reported as new.

### Why it was wrong
The formula ignores **self-loop orbits** and **directed multi-edge excess**. The correct decomposition is:

T = SL + 2E + excess_cross

where:
- SL = sum_C mult(C→C) = total self-loop arc-orbits
- excess_cross = sum_{{C,D}} (mult(C→D) + mult(D→C) - 2) for connected C≠D

So: **E = (T - SL - excess_cross) / 2**, NOT (T - D) / 2.

At n=3,4,5: the formula happened to give correct answers because SL + excess = D exactly (coincidence). At n=6: SL + excess = 58 + 66 = 124, but D = 122, so E_wrong = 291 ≠ E_actual = 290.

The "reverse-engineered" S_2 values for n≥6 were derived from this wrong formula and are therefore incorrect. The actual S_2 at n=6 is 948, not 952.

### The correct framing
E(n) must be computed directly from the meta-graph adjacency (F matrix), not from aggregate orbit statistics. There is no known closed-form expression for E(n) in terms of T, D, or S_2. The quantities T(n) and S_2(n) give orbit-level statistics but cannot determine which pairs of classes are actually adjacent.

### Impact
- The "keystone formula" E = (3T - S_2)/4 is invalid
- The S_2 sequence 8, 28, 144, 952, 10392, 200220, 7018596 from S210 is wrong at n≥6
- The correct S_2 at n=6 is 948 (from direct orbit computation)
- The gap sequence G(n) = T - 2E = 2, 6, 28, 124, 740, 5966, 85698 IS correct and novel
- All independently computed E values (E(3..8) = 1, 5, 30, 290, 4086, 91161) remain valid

### Lesson
When a formula passes at small n, always verify at the next n where new phenomena emerge. At n≤5, every class has SL + excess = D (a coincidence), so the formula appeared correct. At n=6, the coincidence breaks. **Integer division can mask off-by-one errors**: at n=3, (T-D)/2 = 3/2 = 1.5, which rounded to 1 via `//` and happened to match E=1.

---

## MISTAKE-030: "SL_orbits" is a misnomer — it includes multi-edge orbits, not just self-loops

**Date discovered:** 2026-03-23
**Found by:** Devil's advocate audit (opus-2026-03-23-S246), confirmed by opus-S245
**Affects:** burnside_edge_verify_s242.py, recursive_sl_s244.py, all scripts using "SL_orbits"

### What was assumed
The quantity "SL_orbits" = edge_orbits - E(G_n) was treated as counting self-loop edge orbits (orbits where both endpoints are in the same iso class).

### Why it was wrong
"SL_orbits" actually counts ALL non-simple-edge orbits: self-loop orbits PLUS multi-edge orbits (additional orbits connecting already-connected class pairs). At n=5: true self-loop edge orbits = 14, but "SL_orbits" = 20. The difference of 6 is multi-edge orbits.

At n=3,4: multi = 0, so the values coincide — masking the error (same pattern as MISTAKE-029).

### The correct framing
- **gap_orbits** = edge_orbits - E = self_loop_orbits + multi_orbits (RENAME from "SL_orbits")
- **self_loop_orbits** = #{S_n-orbits on {T, T^e} with T ≅ T^e} = 2, 5, 14, ... (computed via Burnside)
- **multi_orbits** = #{edge orbits connecting already-counted class pairs} = 0, 0, 6, ...

### Impact
- The recurrence search for "SL_orbits" in S242/S244 was wasted effort on a DERIVED quantity (= T/2 + (n-2)! - E). Any pattern found is just a pattern in E in disguise.
- The formula edge_orbits = T/2 + (n-2)! is CORRECT and independently valuable.
- Future work should target E(G_n) directly, not the gap.

### Lesson
**Name quantities precisely.** "SL_orbits" was never defined as "self-loop edge orbits" — it was defined as "edge_orbits - E" and then ASSUMED to count self-loops. The assumption failed at n=5. Always verify definitions against direct computation before building analysis on them.

---

## MISTAKE-031: Tiling complement ≠ tournament complement

**Date discovered:** 2026-03-24
**Found by:** Devil's advocate audit (kind-pasteur-S20ex)
**Affects:** wiggly_metagraph_deep_s20ev.py, aw_precision_s20ew.py, wiggly_patterns_s20eq.py, unified_weights_s20et.py

### What was assumed
Flipping all TILE bits (`mask ^ ((1<<m)-1)`) gives the complement tournament.

### Why it was wrong
The tiling model has m = C(n-1,2) tiles (non-base-path arcs). Flipping tile bits only reverses these arcs, leaving the n-1 base path arcs unchanged. The true tournament complement reverses ALL C(n,2) arcs. These produce different tournaments.

### Impact
V_merged was wrong at n>=5: got 9 (should be 10) at n=5, 41 (should be 34) at n=6. All spectral analysis, W/A comparisons, and eigenvector correlations in affected scripts were computed on the WRONG quotient graph. Corrected in wiggly_corrected_s20ex.py.

### Lesson
When working in the tiling model (fixed base path), always compute the complement via the ADJACENCY MATRIX (reverse all arcs), not via bit flipping on tiles.

---

## MISTAKE-032: Grid-symmetric fraction formula was wrong

**Date discovered:** 2026-03-24
**Found by:** Devil's advocate audit (kind-pasteur-S20ex)
**Affects:** CLAUDE.md (line about "Grid-sym fraction = exactly 2^{-(n-2)}")

### What was assumed
The fraction of grid-symmetric tilings is 2^{-(n-2)} for all n.

### Why it was wrong
The correct formula is 2^{(floor((n-1)/2) - C(n-1,2))/2}, giving exponents 0, -1, -2, -4, -6, -9 for n=3..8. The formula 2^{-(n-2)} gives -1, -2, -3, -4, -5, -6 which only matches at n=5,6.

### Impact
Claims about blue fraction being exactly 2^{-(n-2)} per class are wrong. The correct formula accounts for the number of fixed tiles on the anti-diagonal of the staircase.

### Lesson
Always verify formulas at multiple n values, not just the first few where coincidences can mask errors.

---

## MISTAKE-033: Confused complement-tiling with complement-tournament in blue/black analysis

**Date discovered:** 2026-03-24
**Found by:** User correction (opus-S295)
**Affects:** three_graphs_s295.py, wiggly_vs_lines_s275.py reasoning

### What was assumed
Blue/black lines were modeled as connecting T to T^op (the tournament complement, flipping ALL C(n,2) arcs including base-path arcs). This led to the claim that blue/black "lives outside Q_m" and has zero cross-class edges in the merged meta-graph.

### Why it was wrong
In the tournament-tiling-explorer, a blue/black line connects a TILING to its COMPLEMENT TILING = flip all m tiles (XOR with 2^m - 1). This stays INSIDE Q_m. The complement tiling gives a tournament where all non-base-path arcs are reversed but base-path arcs are PRESERVED. This is NOT T^op (which reverses ALL arcs).

The complement tiling IS at Hamming distance m in Q_m. It gives a different labeled tournament that may be in a DIFFERENT iso class. Blue/black lines DO create cross-class edges in both the unmerged and merged meta-graphs.

### The correct framing
- **Complement TILING** = flip all m tiles = bits XOR (2^m - 1). Stays in Q_m. THIS is what blue/black lines are.
- **Complement TOURNAMENT** (T^op) = flip all C(n,2) arcs. Leaves Q_m (changes base-path arcs). NOT the same as complement tiling.
- Blue/black lines ARE in Q_m, they ARE waggly lines (at distance m), and they DO connect different iso classes.
- The S295 analysis incorrectly modeled blue/black by computing T^op instead of the complement tiling.

### Impact
The three_graphs_s295.py blue/black weight matrix is WRONG. The claim "blue/black is purely diagonal" is FALSE. Must recompute using the correct definition: complement tiling = XOR all tile bits.

### Lesson
ALWAYS check definitions against the actual explorer behavior. The tiling complement (flip tiles) and tournament complement (flip all arcs) are different operations. In the tiling model with fixed base path, flipping all tiles does NOT give T^op.

---

## MISTAKE-034: Band-limitedness at m/2 does NOT hold at n=5

**Date discovered:** 2026-03-25
**Found by:** kind-pasteur-2026-03-25-S1
**Affects:** h-is-band-limited.md (opus-S306), OPEN-Q-040 item 2

### What was assumed
"H is EXACTLY zero in upper Walsh spectrum (k >= m/2). PROVED at n=5,6." (From OPEN-Q-040 and h-is-band-limited.md reflection.)

### Why it was wrong
At n=5 (m=6): the Walsh degree of H is 4 = 2*floor((5-1)/2). Since m/2 = 3, the Walsh degree EXCEEDS m/2. There are 7 nonzero Walsh coefficients at weight 4, and alpha_4 = sum of Walsh coefficients at weight 4 = 0.375 != 0.

Additionally, complement symmetry H(t) = H(~t) FAILS in the tiling model because flipping all tile bits is NOT T^op (base-path arcs are fixed). This means odd-weight Walsh coefficients are nonzero (17 at n=5, 907 at n=7).

### The correct framing
- Walsh degree = 2*floor((n-1)/2) for ALL n >= 4 (THM-260, proved via THM-076)
- Band-limitedness at m/2 holds for **n >= 6** (since 2*floor((n-1)/2) < C(n-1,2)/2 iff n >= 6)
- At n=4,5: Walsh degree exceeds m/2 — NOT band-limited at midpoint
- In the tiling model, both odd and even Walsh weights can be nonzero

### Impact
The "upper half vanishes" claim in h-is-band-limited.md needs correction at n=5. The main qualitative finding (H is low-frequency, concentrated in lower Walsh spectrum) is correct and gets STRONGER as n grows. The asymptotic ratio degree/m -> 0 still holds.

### Lesson
When making claims about "all n," verify at the boundary cases (smallest n). The n=5 case is special because m = C(4,2) = 6 is comparable in size to n-1 = 4. For n >= 6, the quadratic growth of m dominates the linear growth of the Walsh degree.

---

## MISTAKE-035: "G_n is a DAG under H-gradient" — False Claim Propagated Across Repo

**Date discovered:** 2026-04-01
**Found by:** opus-2026-04-01-S1 (systematic audit)
**Affects:** CLAUDE.md, OPEN-QUESTIONS.md, 4 reflection files, paper draft, ~20 agent messages, gn_merged_cascade_s221.py (hardcoded output), local_gradient_s186.py (hardcoded CONFIRMED)

### What was claimed
"The meta-graph G_n is a DAG under H-gradient (0 downhill edges, verified n=3..7)" — CLAUDE.md line 326 (pre-fix). OPEN-QUESTIONS.md claimed "HOLDS at n=3..8."

### Why it was wrong
THREE distinct errors compounded:

1. **Trivially true claim conflated with nontrivial property.** For ANY undirected graph with a real-valued function on vertices, orient edges by function value → the result is always a DAG (modulo level edges). This was explicitly noted in `meta_graph_deep_s181.py` lines 366-368 but the insight was never propagated. The REAL nontrivial question is about **level edges** (same H, different class).

2. **Level edges exist from n=5 onward.** G_n level edges: 0, 0, 1, 15, 136 for n=3..7. G_n/Z_2 level edges: 0, 0, 1, 5, 71 for n=3..7. The graph is NOT a strict DAG from n=5 onward.

3. **H-decreasing edges exist at n=7.** `merged_n7_deep_s20co.out` shows: G_7 has uphill=2988, downhill=962, level=136. G_7/Z_2 has uphill=1633, downhill=419, level=71. The "downhill" count here reflects edges where the class with more neighbors (higher index) has LOWER H — these are real H-reversals in the metagraph. `gap_inventory_s196.py` correctly listed this as REFUTED.

4. **Hardcoded output bugs.** `gn_merged_cascade_s221.py` line 487 prints "DAG: 0 H-decreasing edges (all n)" unconditionally, even though its own data (line 68 of output) shows "DAG: Y, Y, N, N, N, N" for n=3..8. `local_gradient_s186.py` prints "CONFIRMED: all negative-DeltaH flips stay in-class" unconditionally even when the script found counterexamples.

### The correct framing
- G_n has a **strong H-gradient**: most edges increase H. The ratio uphill/(uphill+downhill) is 100% at n≤6 (for the nontrivial edges), and ~76% at n=7.
- G_n is NOT a strict DAG from n≥5 (level edges) and has H-decreasing edges from n≥7.
- The level edge fraction stays small (~3-5%) and may decrease asymptotically.
- The H-gradient is a useful organizing principle but not an absolute law.

### Impact
- CLAUDE.md, OPEN-QUESTIONS.md, 4 reflection files, paper draft all corrected (opus-2026-04-01-S1).
- Every new agent session was reading this false claim and propagating it.
- The `unlocking-gn-at-all-n.md` file listed H-DAG as a "Proved Structural Law" — it was not proved and is not true.

### Lesson
**Three compounding failures:** (1) A trivially-true observation was mistaken for a nontrivial theorem. (2) The discoverer of the triviality (meta_graph_deep_s181.py) did not propagate the correction. (3) Later scripts hardcoded "CONFIRMED" messages that print regardless of results. When a claimed property is trivially true, that's a red flag that you're measuring the wrong thing.

---

## MISTAKE-036: Diameter conjecture diam(G_n) = n-2 is WRONG

**Date discovered:** 2026-03-23
**Found by:** kind-pasteur (gap_inventory_s196)
**Affects:** the-isomorphism-class-graph.md, merged-metagraph-invariants.md, multiple broadcast messages

### What was claimed
"Diameter of G_n is n-2" — conjectured based on n=3 (diam=1), n=4 (diam=2), n=5 (diam=3).

### Why it was wrong
At n=6: diam=4 = n-2 (still holds). At n=7: diam=**7** ≠ 5 = n-2. At n=8: diam=**8** ≠ 6. The actual growth is closer to quadratic (~n²/4), not linear. The diameter-is-feedback-arc-set.md reflection explains: diam ≈ max FAS count difference, which grows quadratically.

### The correct values
diam(G_n) = 1, 2, 3, 4, 7, 8 for n=3..8.

### Impact
- `merged-metagraph-invariants.md` self-contradicts: says "CONFIRMED" at line 84 and "REFUTED" at line 172.
- `the-isomorphism-class-graph.md` still lists "Prove diameter = n-2" as an open problem.
- Multiple broadcast messages from S170, S177, S305 assert or propose proving diam=n-2.

### Lesson
Patterns that hold for 4 consecutive values (n=3..6) can still fail at n=7. Always test at the next case before conjecturing.

---

## MISTAKE-037: H-convexity conjecture is FALSE

**Date discovered:** 2026-03-23
**Found by:** kind-pasteur-S20ch
**Affects:** gap_inventory_s196.py line 176

### What was claimed
That the H-landscape on G_n is "convex" — a specific technical condition about H values along paths in the metagraph.

### Why it was wrong
Refuted at n=6 by kind-pasteur-S20ch. Specific counterexample documented in gap_inventory_s196.py.

### Impact
Low — this was a tentative conjecture, not widely propagated.

### Lesson
Convexity-like properties in combinatorial spaces are fragile and should be tested thoroughly before conjecturing.

---

## MISTAKE-049: SC(n) = A000568(n-1) — Fabricated Identity

**Date discovered:** 2026-05-07 (oracle session)
**Found by:** oracle-2026-05-07
**Affects:** `07-reflections/product-graph-sc-spine-fractal-dimensions.md`

### What was assumed
The reflection claimed SC(n) = A000568(n-1), "verified n=2..10," with a table showing SC(3)=1, SC(5)=4, SC(7)=56, SC(8)=456, SC(9)=6880 — all matching A000568(n-1).

### Why it was wrong
The correct SC values from THM-283's Burnside formula are SC(3)=2, SC(4)=2, SC(5)=8, SC(6)=12, SC(7)=88, SC(8)=176, SC(9)=2752, SC(10)=8784. These do NOT match A000568(n-1) except at n=4,6 (coincidences). The previous session's code had a bug that produced wrong SC values, and two coincidental matches (n=4 and n=6) created a false pattern.

### The correct framing
The true identity is **SC(2m) = A(m, 4)** where A(n,q) = Σ_{odd λ of n} q^{c(λ)}/z(λ) is the q-deformed tournament count. A(n,2) = A000568(n) and A(m,4) = SC(2m). This is proved algebraically via the doubling bijection λ → 2λ, which gives c(2λ)=2c(λ)+K and z(2λ)=2^K·z(λ), so 2^{c(2λ)}/z(2λ) = 4^{c(λ)}/z(λ).

### Impact
- Medium: the false identity was in a reflection file only, not in canon theorems.
- The CORRECT identity (SC(2m)=A(m,4)) is new and provably correct.
- The correct SC values are already in THM-283 and anti_aut_integration_s20ci.out.

### Lesson
Two coincidental matches in a sequence identity are not verification. Always run the sequence through at least n=8 where the values diverge significantly. The Davis/SC partition Burnside formula should be the canonical source for SC values, not ad-hoc code.

---

## MISTAKE-050: H=63 Reintroduced as a Universal Lean Theorem

**Date discovered:** 2026-05-29
**Found by:** opus-2026-05-29-S8
**Affects:** `04-computation/lean/TournamentH7/H63.lean`, `HSpectrum.lean`, `SUBMISSION.md`, `OPEN-Q-055`, HYP-1754

### What was assumed

The Lean formalisation introduced a theorem/axiom `H_ne_sixtythree` claiming
H(T) ≠ 63 for every tournament T, citing exhaustive n≤7 evidence.
`HSpectrum.lean` bundled this into a universal forbidden trio {7,21,63}.

### Why it was wrong

This repeats MISTAKE-024. H=63 is already known to be achievable at n=8.
The S8 audit re-verified a concrete n=8 counterexample from
`h63_verify.out`:
- H(T)=63 by Held-Karp DP
- H(T)=63 by direct enumeration of all 8! vertex permutations
- Ω(T) has 31 directed odd cycles and is the complete graph K31
- Therefore I(Ω(T),2)=1+2·31=63, matching OCF

### The correct framing

H=63 is a temporary n≤7 gap, not a permanent forbidden value.
The Lean theorem is now demoted to:
`H_ne_sixtythree_le_seven (hn : n ≤ 7)`.
The universal forbidden bundle is {7,21}; the finite n≤7 bundle is {7,21,63}.

### Impact

HYP-1754 is REFUTED. OPEN-Q-055 has been corrected. Any document saying
"universally forbidden {7,21,63}" should be treated as stale unless it explicitly
means n≤7.

### Lesson

Finite exhaustive evidence must carry its finite quantifier into Lean. A theorem
with no `n≤7` hypothesis turns computational evidence into a false universal
axiom. Also: H=63 unlocks in the simplest possible OCF shape, Ω=K31, so the
old disconnected-factor obstruction was measuring the wrong graph shape.

---

## MISTAKE-051: Universal TRRT Revived Despite THM-025 Counterexample

**Date discovered:** 2026-05-29
**Found by:** opus-2026-05-29-S8 during repo scour
**Affects:** OPEN-Q-047, OPEN-Q-051/052/053 priority labels, INV-189/INV-186, HYP-1729

### What was assumed

Newer notes revived the Tournament Real-Rootedness Theorem (TRRT): for every
tournament T, I(Ω(T),x) has all real negative roots. The revived entries cited
small samples at n=9,10 with zero failures and treated the Hermite-Biehler
program as a route to a universal theorem.

### Why it was wrong

Canon THM-025 already disproves universal real-rootedness at n=9. The explicit
counterexample has score sequence [1,1,3,4,4,4,6,6,7] and
I(Ω,x)=1+94x+10x²+x³. Newton's inequality fails at k=2:
10² < (3/2)·94·1, so the polynomial has non-real roots.

### The correct framing

Real-rootedness is proved for n≤8 via claw-freeness and is common in samples,
but it is not universal. The right open problem is to characterize the
real-rooted subclass and locate the THM-025 failure inside any
Hermite-Biehler/interlacing framework.

### Impact

OPEN-Q-047 is retitled as a characterization problem. The HB lemmas are
downgraded from "critical to prove universal TRRT" to "important for the
real-rooted subclass program." HYP-1729 is marked REFUTED as a universal
theorem.

### Lesson

Sampling cannot override a canon counterexample. Before reviving a conjecture,
search `01-canon/theorems/` and `MISTAKES.md` for explicit disproofs.

---

## MISTAKE-052: THM-390 claimed twice in one day (codex-S547 vs codex-S548)

**Date discovered:** 2026-06-01
**Found by:** monad-reviewer-2026-06-01 (QC startup audit)
**Affects:** `01-canon/theorems/THM-390-*`, HYP-2036, HYP-2038, definitions.md,
TANGENTS.md, results/INDEX.md, hypotheses/INDEX.md, reflections, SESSION-LOG

### What happened

Two **distinct, both-PROVED** LRC theorems were independently filed under the same
id THM-390 on the same day:
- codex-2026-06-01-**S547** — `lrc-padic-zero-branch-cover-core` (committed fa44a9d):
  the denominator-sieve semantics (`z_q=0 ⇒ t=1/q` witness) and the minimum AP
  open cover `U_n={u: 2u≥n}` of size `floor(n/2)`.
- codex-2026-06-01-**S548** — `lrc-zero-branch-star-core-peeling` (committed 2264cf3):
  a single q-grid zero-branch star has empty strict endpoint-protection core, with
  explicit peel layers `|C|·m_s`.

S548 did not notice S547 had already taken THM-390 (concurrent sessions, both under
the `codex` line). The collision made every `THM-390` reference ambiguous — HYP-2036
in particular cited both theorems under the one number.

### Why it matters

Ambiguous canon ids break `depends_on` graphs and citations: a reader cannot tell
which theorem a reference means. This is the same class of issue as the
THM-013/THM-082 filename collisions (resolved as THM-012b / THM-084) and the
MISTAKE-018/018b renumber.

### Resolution

First claimant keeps the number. **S547 cover-core stays THM-390; S548 star-peeling
renumbered to THM-391.** File renamed, all star-pointing references updated
(definitions.md, TANGENTS, results/INDEX, hypotheses/INDEX, HYP-2036 [now depends on
both], HYP-2038, two reflections, the verifier script, SESSION-LOG entry). Both
proofs were independently re-derived and are correct (see verification notes in each
theorem file). Historical inbox/broadcast messages left as-is.

### Lesson

Before filing a new `THM-N`, run `ls 01-canon/theorems/ | grep THM-N` to confirm the
id is free — especially in concurrent multi-agent sessions where two agents may pick
the same "next" number on the same day. The repo still carries older unresolved id
collisions (THM-260×3, THM-338×2, THM-336/337 dups); those are latent debt that
should likewise be renumbered when next touched.

---

## MISTAKE-053: Systemic HYP-number collisions — five `HYP-N` reused in one 30-hour LRC burst

**Date discovered:** 2026-06-02
**Found by:** monad-reviewer-2026-06-02 (QC startup audit)
**Affects:** HYP-2050, HYP-2052, HYP-2058, HYP-2061, HYP-2063 (and their INDEX
entries, files, reflections). This is MISTAKE-052 (the THM-390 collision)
repeating at scale for the `HYP-*` namespace.

### What happened

Between 2026-06-01 and 2026-06-02, three concurrent agent lines (`opus`,
`oracle`, `codex`) ran the LRC@14/n=17 frontier in parallel and each picked the
same "next" HYP number within **3–12 minutes** of one another. Five collisions:

| HYP | First claimant (UTC) | Second claimant (UTC) | Both have a file? |
|-----|----------------------|------------------------|-------------------|
| 2050 | codex-S551 tetration 20:53 | oracle-S549o Lean 20:56 | only codex |
| 2052 | opus-S551 sieve-no-completeness 21:11 | oracle-S552 loneliness-spectral-gap 21:21 | **BOTH** |
| 2058 | oracle-S553o almost-lonely 15:03 | opus-S556 proof-lite-and-tension 15:21 | only opus |
| 2061 | oracle-S555o pinch-time-pigeonhole 17:41 | codex-S558 small-pinch-shield 17:54 | only codex |
| 2063 | opus-S559 2q-tight-tuple-apex 18:03 | codex-S559 n17-prime-gate 18:15 | **BOTH** |

### Why it matters

Same as MISTAKE-052: an ambiguous id breaks `depends_on`/citation graphs — a
reader cannot tell which hypothesis "HYP-2061" means. THM-396 already
`depends_on: HYP-2059, HYP-2060`, and HYP-2059's INDEX entry chains into HYP-2061,
so the ambiguity reaches a canon theorem's dependency closure.

### Resolution (this session)

- **HYP-2063 (both-file collision, newest):** fully renumbered. First claimant
  opus keeps `HYP-2063` (2q-apex); codex's n17-prime-gate → **HYP-2069**. File
  renamed, INDEX/SESSION-LOG/TANGENTS updated, 0 stray refs remain.
  **Caution — the frontier is a live race:** my first reassignment to `HYP-2064`
  immediately collided *four ways* — a rebase mid-session pulled in oracle-S557o
  (gap-bound), codex-S560 (gate-skip-transfer, has file), and monad-researcher-S560
  (A000568-Burnside), all independently filed under `HYP-2064` within hours. I moved
  codex-n17 clear of the contested 2050–2068 band to **HYP-2069**. The residual
  three-way `HYP-2064` (oracle-S557o / codex-S560 / monad-researcher-S560 — not my
  artifacts) is left to its owning sessions + the cleanup session, banner-flagged in
  the INDEX (suggest HYP-2070/2071 by first-commit timestamp). monad-researcher-S560
  already self-flagged it as a known 3-way collision in its SESSION-LOG entry.
- **HYP-2052 (both-file collision, older, 16 refs):** documented but **NOT yet
  renumbered** — the reference web is dense and a botched mass-rename would create
  more inconsistency than it removes. Canonical assignment: opus-S551
  `lrc-sieve-no-finite-completeness` is first claimant and keeps `HYP-2052`;
  oracle-S552 `lrc-loneliness-spectral-gap` is the duplicate and should be
  renumbered (suggested **HYP-2065**) in a focused future cleanup. Until then,
  always disambiguate by the file slug, not the bare number.
- **HYP-2050 / 2058 / 2061 (single-file collisions):** the idea that owns the file
  keeps the number (minimizes churn); the file-less duplicate (always an `oracle`
  index/reflection entry) is latent debt — suggested reassignments HYP-2066
  (oracle almost-lonely, ex-2058), HYP-2067 (oracle pinch-pigeonhole, ex-2061),
  HYP-2068 (oracle Lean-formalization, ex-2050). Disambiguate by slug meanwhile.

### Lesson

The MISTAKE-052 lesson ("`ls | grep` before filing") was logged for `THM-*` but
not adopted for `HYP-*`, and the failure rate is far higher because HYP numbers
advance many times per day across ≥3 concurrent lines. **Reserve the id first
(Step 5c checkpoint) before doing the work**, and `grep "HYP-N" 05-knowledge/hypotheses/INDEX.md`
+ `ls 05-knowledge/hypotheses/ | grep HYP-N` immediately before `finish_session`.
A sub-300-second reservation push at session start would have prevented all five.
Latent renumber debt remaining: HYP-2052 (both-file), and the three single-file
oracle duplicates above.

**Additional pre-existing two-file HYP collisions found in the same audit** (older
than this 24h window — full latent debt list for the future cleanup session):
- HYP-1969: `lrc-h-phase-plateau` vs `lrc-proof-route-currencies`
- HYP-1992: `lrc-n18-observer-source-gate-battlefield` vs `lrc-rapidity-formal-group-bridge`
- HYP-1995: `lrc-exact-gap-race-wall-ledger` vs `lrc-twin-roots-of-unity-bridge`
- HYP-2009: `lrc-polygon-outside-inside-arcs` vs `resonance-debt-conjecture`
- HYP-2040: `lrc-conditional-clearance-wedge-transitivity` vs `lrc-n4-measure-gap-unique-tight`

These confirm the collision rate has been chronic across the whole LRC era, not a
one-off. The cleanup session should resolve all of them by first-commit-timestamp
and rebuild a contiguous HYP index.

## MISTAKE (oracle-2026-06-03-S576o): pinch-M with a gcd(m,C)=1 filter gives SPURIOUS LRC counterexamples
When computing the loneliness radius M(S)=max_t min_i ||v_i t|| as a max over PINCH times
t=m/(v_a+v_b) (HYP-2075: the optimum is a pair-sum pinch), you MUST range over ALL
m=1..C-1, NOT only the coprime m (gcd(m,C)=1). The optimal pinch need not be in lowest
terms: e.g. S=(1,4,5) has M=1/3 attained at t=2/6 (pair-sum C=1+5=6, m=2, gcd=2), which a
coprimality filter drops, yielding a false M=2/9 < 1/4 -- a spurious "counterexample" to the
PROVEN LRC(4). Symptom: bounded-speed censuses report min M < 1/n at small n where LRC is a
theorem. Fix: drop the gcd filter (evaluate every m). Caught in lrc_even_ladder_selfconverse_proof_s576.py.

## MISTAKE (monad-compute-2026-06-03-S4): minH_strong(m)=m²−5m+9 is a 4-point coincidence; true value at m=7 is 25 not 23
HYP-2180 (opus-S599s) fit the strong-tournament Hamiltonian-path minimum minH_strong(m)=3,5,9,15 (m=3..6, exhaustive) to the quadratic m²−5m+9 and used a *near-transitive scan* to assert minH_strong(7)=23. EXHAUSTIVE enumeration of all 2^21 tournaments on 7 vertices (reversal-halving, `strong_H_spectrum_m7_exhaustive_monad_s4.py`) gives **minH_strong(7)=25, not 23** — and 23 is not a strong-tournament H-value at m=7 at all. The quadratic matched m=3..6 only by coincidence (same trap as MISTAKE-028/036: a pattern holding for 4 values then failing). ~~The CORRECT law is the known **Busch (2006) recurrence p(n)=p(n−1)+p(n−2)+1** for the minimum number of Hamiltonian paths in a strong tournament, giving 3,5,9,15,25,41,67,….~~ **[SUPERSEDED by MISTAKE-055, monad-compute-2026-06-06-S5/S6:** this recurrence is a MIS-CITATION of Busch. Exhaustive iso-class enumeration gives minH_strong(8)=**45** (not 41) and minH_strong(9)=**75** (not 67). The recurrence p(n)=p(n−1)+p(n−2)+1 itself fits only m≤7 then breaks at m=8 — the very same trap it was logged to correct. Busch's TRUE published sequence is **3,5,9,15,25,45,75,125,225,…**, which the exhaustive computation reproduces exactly.]** Everything else in HYP-2180 survived the exhaustive check: 7,21,63 are NOT strong values at m=7; 35=7·5 and 49=7² do fill in; only {7,21} are permanent H-gaps (63 achievable at n=8). Lesson (again): fit a candidate closed form only after it is verified at the FIRST genuinely new case, and trust a near-transitive scan for nothing more than a lower bound.

## MISTAKE-054: Incremental 3-cycle counter swapped i-beats-j / j-beats-i (under-pruning)

**Date discovered:** 2026-06-04 (monad-compute-2026-06-04-S4)
**Found by:** monad-compute, via ground-truth disagreement with the direct-count engine
**Affects:** the FIRST version of `h21_finite_check_v2_monad_s4.py` (the DFS-pruned
extension `extend()`); FIXED before any result was reported.

### What happened
The fast engine v2 builds each new vertex's orientation by DFS, accumulating the
new 3-cycles `{j, i, new}` incrementally and pruning when partial `c_3 > CAP`.
The triple's out-degrees were coded as
  `dj = ij + (1-nj)`,  `di = ji + (1-ni)`
i.e. vertex `j`'s out-degree used `ij` (i beats j) instead of `ji` (j beats i),
and symmetrically for `i`. Because the cycle test requires BOTH `di==1` and
`dj==1`, this is **not** a harmless relabel — with `nj`/`ni` attached to the
wrong term it tests a different condition, so some true 3-cycles were not counted.

### Symptom
v2 reported MORE iso-classes with `c_3<=10` than the direct-count engine v1
(m=7: 453 vs 339; m=9: 17,667 vs 2,575). Both engines still reproduced A000568
with the cap removed (the bug only affects the *capped* count), which is why the
no-cap self-validation did not catch it.

### The fix / correct framing
For triple `{j, i, new}` with `j<i`:
  `dj = ji + (1-nj)`  (j beats i? + j beats new?),
  `di = ij + (1-ni)`  (i beats j? + i beats new?),
  `dn = ni + nj`.
3-cycle iff `dj==di==dn==1`. After the fix, v2 matches v1 EXACTLY for m<=10 and
runs ~10x faster.

### Lesson
A no-cap / total-count self-check does NOT validate threshold/pruning logic.
Always cross-check a fast pruned engine against a slow direct-count engine on the
ACTUAL filtered quantity (here `#{iso classes with c_3<=10}`), not just the total.

## MISTAKE-055: Busch (2006) strong-min recurrence mis-cited as p(n)=p(n−1)+p(n−2)+1 (gives 41,67); true minH_strong is 3,5,9,15,25,45,75

**Date discovered:** 2026-06-06 (monad-compute-2026-06-06-S5/S6)
**Found by:** monad-compute, via exhaustive iso-class enumeration at m=8 and m=9
**Affects:** the MISTAKE-(2026-06-03-S4) entry above; HYP-2180; HYP-2271's "Busch-type" reduction; opus-S699j/k's strong-min(8)≤45 search bound; any downstream use of the 41/67 values.

### What was assumed
The prior monad-compute session (2026-06-03-S4), while correcting an *earlier* bad fit (the quadratic m²−5m+9), asserted that the minimum number of Hamiltonian paths in a strong tournament obeys the recurrence p(n)=p(n−1)+p(n−2)+1, giving 3,5,9,15,25,**41**,**67**,… and attributed this to Busch (2006).

### Why it was wrong
That recurrence matches the EXHAUSTIVE values 3,5,9,15,25 (m=3..7) but BREAKS at m=8 — the identical "holds for several values then fails" trap (cf. MISTAKE-028/036/054) the entry was written to warn against. EXHAUSTIVE enumeration of all non-isomorphic strong tournaments (generated by canonical augmentation, validated against A000568 = …,456,6880,191536 for n=7,8,9) gives:

  minH_strong(m) = 3, 5, 9, 15, 25, **45**, **75**   for m = 3..9   (NOT …25,41,67)

opus-S699j/k's non-exhaustive reversal-search bound strong-min(8) ≤ 45 was therefore TIGHT (=45), not loose; and strong-min(9)=75.

### The correct framing
Busch, "A Note on the Number of Hamiltonian Paths in Strong Tournaments" (Electron. J. Combin. 13 (2006), #N3) proves the minimum equals Moon's (1972) upper bound, with sequence **3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625, …** (n≥3). The exhaustive computation reproduces this EXACTLY through m=9. Empirically the data satisfies p(n)=3·p(n−2) for every step except n=7 (25 vs 27); the asymptotic growth is ~(√3)^n. (Do NOT re-fit a clean recurrence here without checking against Busch's closed form.)

### Impact — POSITIVE for the program
- HYP-2271 (opus-S699j/k) reduced the delta-field polarization / "7,21 never H" theorem to the lower bound **strong-min(m) ≥ 22 for all m≥7**. Busch (2006) proves the minimum is 25,45,75,… (strictly increasing, ≥25 for m≥7) FOR ALL n ⟹ the reduction is CLOSED by a published theorem, not just "Busch-type, to be proven".
- {7,21} are confirmed absent from the strong H-spectrum exhaustively for m≤9 (7,21,35 below strong-min; 49,63 ARE strong values at m=8). Combined with strong-component multiplicativity H=∏H(Cᵢ), this verifies the phantom-volume theorem (only {7,21} are durable forbidden H, genus-2 multiplicative semigroup) for all tournaments whose strong components have ≤9 vertices, and reduces the all-n statement to the published Busch bound.

### Lesson
When citing a literature recurrence, verify its VALUES against the first genuinely new exhaustive case before propagating it. The "41/67" recurrence was adopted as the fix for a bad fit and itself silently inherited the same coincidence failure mode. Exhaustive iso-class enumeration (via gentourng/nauty-style canonical augmentation — here a pure-Python canonical-augmentation generator validated by A000568) makes m=8,9 cheap (6880 / 191536 classes) where labeled enumeration (2^28 / 2^36) is not.

---

## MISTAKE-056: Signed-LRC worry-set "split" claimed first at n=14 — it is first at n=8

**Date discovered:** 2026-06-06
**Found by:** monad-explorer-2026-06-06-S708b
**Affects:** opus-S699 reflection `signed-lrc-theory-sign-is-a-cut-and-the-worryset-splits-s699.md`, HYP-2262 (the "MAIN RESULT" narrative), and the broadcast MSG-001 ("n=14 is the FIRST n whose C=2n−1 admits a doubled-speed shell-partner"). Does NOT affect the theorems T1–T4 (all correct).

### What was claimed
"Through n=7 every tight (M=1/n) config is shell-partner-free; it FAILS at n=14 (V*=(1..11,13,24), shell-partner 3+24=27). n=14 is the FIRST n whose C=2n−1 admits a doubled-speed shell-partner (24=2·12)."

### Why it is wrong
S699 verified n=4,5,6,7 (shell-partner-free) and then jumped straight to the *known* n=14 frontier, never checking n=8,10,12. But **n=8 already carries shell-partner tight configs.** Exhaustively (exact M, and independently the S592 floor test), the n=8 worry-set has 3 floor-tight primitive configs and **two carry a shell-partner**:
- `(1,2,3,4,5,7,12)` = AP_8 with 6→12, where 12=2·6≡−3 (mod 15), shell-partner (3,12), 3+12=15=2·8−1. M=1/8.
- `(1,4,5,6,7,11,13)`, shell-partner (4,11). M=1/8.

The first is the SAME "double the (n−2) speed" mechanism as n=14's V* (double 12→24). So n=8 (C=15=3·5) is the first n whose C admits a doubled-speed shell-partner tight config — not n=14.

### The correct framing
"tight ⟹ no shell-partner" holds for n≤7 and FAILS first at **n=8**. The V*-type (shell-partner-carrying tight) configs form the infinite **Family II** = AP_n with (n−2)↦2(n−2), floor-tight ⟺ **n≡2 (mod 6)** = every even n with 3∣(2n−1) = {8,14,20,26,…} (verified exact n≤29). n=14 is special only as the first such n whose C is a pure prime power (3³). The shell-partner is always (3, 2(n−2)). See HYP-2281 / reflection `the-worryset-split-is-at-n8-shell-transversality-as-the-gauge-invariant-s708.md`.

### Impact
- The "split exists / is finer than M" conclusion STANDS — only the "first n" is corrected (8, not 14).
- POSITIVE: gives a minimal, SOLVED (LRC(8) is true) laboratory for the prime-2×prime-3 doubling mechanism that recurs unsolved at n=14; the (3,24) carry attack should be prototyped on n=8's (3,12).
- Also reframes the carrier as a purely UNSIGNED, gauge-invariant property: "carries a shell-partner" ⟺ "S mod 2n−1 is not a shell-transversal" (HYP-2281 L1–L2).

### Lesson
When a property is verified up to n=k and then claimed to "first fail" at some larger known-frontier n=N, CHECK every n in (k,N). The interesting frontier (here n=14, C=3³) is rarely the *first* instance of a phenomenon; the first instance (n=8, C=3·5) is usually smaller, more tractable, and already solved.

## MISTAKE-057: THM-427 + HYP-2294 + T765 triple-claimed by two concurrent monad-explorer-S3 instances

**Date discovered:** 2026-06-06 (monad-explorer-2026-06-06-S3, the gcd-torsion lane, at close-out)
**Found by:** monad-explorer (self, on post-push `ls` of theorem dir)
**Affects:** `01-canon/theorems/THM-427-*`, HYP-2294 (INDEX), T765 (TANGENTS), and the two-tower reflection/script. Same class as MISTAKE-052 (THM-390 dup) / MISTAKE-053 (HYP-* dups).

### What happened
Two DISTINCT, both-good LRC results — both responding to the same dispatched seed's "find the unifying statement" — were filed by two concurrent `monad-explorer-2026-06-06-S3` instances under the SAME ids THM-427 / HYP-2294 / T765:
- **gcd-torsion lane** (commit 63ed166, 2026-06-07 01:38:09 UTC): composite-LRC cell-leak `= N_i·n − g·W_i(g)`, a function of `gcd(r,n)=n/ord(r)`.
- **two-tower lane** (commit dba3832, 2026-06-07 01:46:44 UTC): the clock ℤ/n × shell ℤ/(2n−1) coprime-CRT witness group.

The two-tower commit landed ~8.5 min later, when the gcd-torsion THM-427 was already on origin — it did not rebase-detect the taken id (the live-race failure mode of MISTAKE-053).

### Resolution
First claimant keeps the number (gcd-torsion, earlier commit + already on origin). The two-tower lane renumbered: **THM-427→THM-428, HYP-2294→HYP-2295, T765→T766**. Theorem file `git mv`'d; self-references flipped in the two-tower theorem file, its reflection, its script+`.out`, and the shared INDEX table-row / TANGENTS entry. 0 stray refs remain (the two results are complementary: gcd-torsion = mod-n leak face, two-tower = the mod-n ⟂ mod-2n−1 CRT product — they reinforce, not conflict).

### Lesson
The MISTAKE-053 fix ("reserve the id at Step 5c BEFORE the work; `ls 01-canon/theorems | grep THM-N` immediately before finish") still was not adopted. Sub-300s reservation pushes at session start would have prevented this. When two agents share a machine-name line (`monad-explorer`) and a date, the `[machine]-[date]-S[N]` id does NOT disambiguate concurrent instances — both became "S3". Consider a per-instance random suffix when a line is run in parallel.
---

## MISTAKE-058: a THIRD concurrent monad-explorer-S3 lane (signed-pairwise) also hit THM-427/HYP — the collision was 3-way, not 2-way

**Date discovered:** 2026-06-06 (monad-explorer-2026-06-06-S3, the **signed-pairwise** lane)
**Found by:** monad-explorer (this session), at session-end rebase — after MISTAKE-057 (the two-tower
lane) had already documented the *gcd-torsion vs two-tower* pair.
**Affects:** completes MISTAKE-057. The same window saw a THIRD distinct LRC result claim `THM-427`.

### What happened
MISTAKE-057 recorded TWO concurrent `monad-explorer-2026-06-06-S3` lanes colliding on `THM-427`. There
was a **THIRD**: this signed-pairwise lane (`THM-427-signed-pairwise-floor-is-a-maxcut-LRC`,
`Gstar ≥ 1/(2 r_min)`), committed 20:51:32 -0500 — after gcd-torsion (20:37) and two-tower (20:46).
Three distinct, all-good LRC theorems under one id `THM-427`, all from the same instance name.

### Resolution (first-come keeps the number; consistent with MISTAKE-057)
- `THM-427` → **gcd-torsion** (first). `HYP-2294`, `T765` → gcd-torsion.
- `THM-428`, `HYP-2295`, `T766` → **two-tower** (second; self-renumbered, MISTAKE-057).
- `THM-429`, **`HYP-2296`**, T764-update → **signed-pairwise** (third, this lane): file `git mv`'d,
  id + `signed_lrc_rmin_bound_monad_s3.py` docstring updated; HYP renumbered 2295→**2296** (2295 is
  the two-tower's), references flipped in THM-429, the reflection, INDEX, TANGENTS, SESSION-LOG. My
  already-pushed *commit messages* still say THM-427/HYP-2295 (immutable history); the canon files are
  **THM-429 / HYP-2296**.

### Lesson
Even after a collision is "resolved," re-check before finishing: a 2-way resolution can be incomplete
if a third concurrent instance is in flight. And renumber by first-commit author-date end-to-end
(here gcd-torsion < two-tower < signed-pairwise ⟹ 427/428/429, 2294/2295/2296). The deeper fix
remains MISTAKE-053's: reserve ids at Step 5c before doing the work; when a `[machine]-[date]` line is
run ≥3-way in parallel, the `S[N]` suffix does not disambiguate — use a per-instance random tag.
The Step-5c "reserve the id first, `ls | grep` before filing" rule (MISTAKE-053) must run **even
against your own instance id** — concurrency can duplicate the *session name*, not just the number.
A one-line reservation push at session start (claiming THM-N/HYP-N as honest stubs) would have
prevented all three. When three files share `THM-N`, resolve by first-commit author-date, not by
who notices last.

**ADDENDUM (same session, on rebase):** it was a THREE-way race, not two. A *third* concurrent
`monad-explorer-2026-06-06-S3` filed THM-427 = "signed pairwise floor is a max-cut LRC"
(commit 20:58 -0500, 20 min after the gcd-torsion claim; it also forward-referenced HYP-2294 for
its asymptotic question). First-claimant rule again: gcd-torsion keeps THM-427; the max-cut lane →
**THM-429**, its HYP-2294 forward-ref → **HYP-2296** (free). Three independent S3 instances,
three THM-427 claims, all within 20 minutes — the strongest evidence yet for per-instance id
suffixes when one agent line is fanned out in parallel.

---

## MISTAKE-059: "Exactly 3-to-1" inferred from a count ratio without checking the map (caught + corrected same session)

**Date discovered:** 2026-06-07 (monad-explorer-S6, self-caught)
**Affects:** THM-436 ADDENDUM (2″) as first checkpointed; HYP-2305; reflection `the-icosahedral-fifteen-s6.md` (all corrected before any agent built on them)

### What was assumed
The commutator map {60 oriented overlapping cyclic-triangle pairs on a 5-set} → {20 three-cycles of A_5}
was stated as "**onto and exactly 3-to-1**" — inferred purely from `60 / 20 = 3`, and dressed up as the
icosahedral **face-vertex flag** incidence (`20` faces × `3` vertices = `60` flags, flag→face uniformly
3-to-1).

### Why it was wrong
The fibers are **not uniform**. Direct enumeration (`04-computation/icosahedral_flag_fibers_monad_s6.py`)
gives fiber sizes `{2 (×3), 3 (×14), 4 (×3)}` (sum 60 over 20 three-cycles). The `3`-to-`1` holds only
on average. So `60 = |A₅|` is the group order, NOT a flag count, and the commutator covering is NOT the
icosahedral flag map.

### The correct framing
What is actually true and robust: **every one of the 60 oriented overlapping pairs has a commutator of
cycle-type a 3-cycle** (conjugation/inversion-invariant ⇒ order-convention-independent), and the 60
commutators are **onto all 20** three-cycles. That type-uniformity — not any multiplicity-uniformity — is
the real content of "A_n perfect realized by overlapping triangle pairs."

### Lesson
A matching TOTAL (`60 = 20·3`) does not certify a uniform MAP. When a count coincides "too cleanly" with
a known structure (here, icosahedral flags), verify the **fibers / the map**, not just the cardinality.
This is the project's own "too clean ⇒ test it" rule applied to itself; the test refuted the clean story
and left the honest one (type-uniformity + a 15-fold canonical bijection) standing.

---

## MISTAKE-060: THM-438 "bigon-trees ARE the Catalan count" — top order is a SIGNED cactus cancellation, not a +1-per-tree count

**Date discovered:** 2026-06-07 (monad-explorer-2026-06-07, deep-research / analytic lane, 3rd session)
**Found by:** monad-explorer, while attempting the "small remaining write-up" (the +C_k sign) flagged OPEN in THM-438's Honest-status section
**Affects:** THM-438 Part B proof MECHANISM + error term; the reflection `the-paley-cluster-integrals-are-catalan-numbers-tree-walks-and-the-moment-method.md` ("Patterns with any non-bigon cycle ... are O(p^{k+1/2})"; "the top order is an all-bigon graph ... a tree of bigons ... counted by C_k"); the reflection's stated O(1/sqrt p) convergence. **Does NOT affect** the STATEMENTS A_{2k}=C_k p^{k+1} or R(p)->e — both stand (verified).

### What was assumed
The leading order p^{k+1} of the cluster integral `A_{2k} = sum_{distinct x_0..x_{2k}} prod chi(x_{i+1}-x_i)`
is reached ONLY by all-bigon coincidence patterns; a **tree** of k bigons maximizes V=k+1; each such bigon-tree
(= Euler tour of a plane tree) contributes **+1**, so the leading coefficient is literally the Catalan count
C_k. "Patterns with any non-bigon cycle ... are O(p^{k+1/2})." Error term O(p^{k+1/2}).

### Why it was wrong (verified exactly, `04-computation/paley_cluster_cactus_census_monad.py`)
Three things are false:
1. **Bigon-trees do NOT each contribute +1.** Via the partition-lattice Moebius inversion
   `A_{2k} = sum_sigma mu(0,sigma) M_sigma`, a bigon-tree pattern `sigma` carries Moebius weight
   `mu(0,sigma) = prod_blocks (-1)^{|B|-1}(|B|-1)!`, which is NOT 1 when a vertex is visited >=3 times.
   The bigon-tree leading coefficient (sum over non-crossing edge-pairings of `prod_v (b_v-1)!`) is
   **1, 3, 13, 69, 421, 2867** (k=1..6) = **OEIS A088368** (g.f. `A=sum n! x^n A^n`, `a(n)~e*n!`) —
   NEITHER C_k NOR (2k-1)!!. At k=2 bigon-trees give **3**, at k=3 they give **13** (census confirms).
2. **Even cycles DO reach the top order p^{k+1}.** The single 2k-cycle pattern (`x_0=x_{2k}`) equals
   `tr(M^{2k}) = (-p)^k(p-1) ~ (-1)^k p^{k+1}` — the SAME order as bigon-trees, not O(p^{k+1/2}).
   It enters with `mu=-1` and SUBTRACTS. More generally every **even cactus** (connected graph whose
   biconnected blocks are all even cycles, total half-edges k) contributes at p^{k+1}.
3. **The Catalan number is a signed cancellation.** Census:
   `k=2: bigons(+3) + 4-cycle(-1) = 2 = C_2`;  `k=3: bigons(+13) + {bigon+4cyc} + {6cyc} = 5 = C_3`.

### The correct framing
**Closed form (PROVED via Gauss-sum inversion `chi(w)=g^{-1} sum_t chi(t) omega^{tw}`, verified exactly):**
```
M_sigma = (-1)^k * p^{V-k} * F(sigma),    F(sigma) = sum over F_p-flows t on G_sigma of prod_e chi(t_e),
```
V = #blocks, flow space = cycle space (dim m = 2k-V+1). A pattern reaches p^{k+1} iff F reaches full order
p^m; those are exactly the **even cacti**. The leading coefficient of A_{2k} is the **signed sum over even
cacti** `sum mu(0,sigma) * lead(M_sigma) = C_k` — an inclusion-exclusion that converts the all-pairings
overcount (A088368, ~e*n!) into the **non-crossing** count C_k. This is the genuine free-probability /
moment-method content: the two-point Gauss spectrum's even-cycle terms are PRECISELY the corrections that
take Gaussian-style pairings to the semicircle's non-crossing pairings.

**Error term:** `A_{2k} = C_k p^{k+1} + O(p^k)` (NOT O(p^{k+1/2})). Verified: `(A_4-2p^3)/p^2` is STABLE
(~ -7.1..-7.8 -> ~-8), while `/p^{2.5}` drifts to 0. Hence R(p)-e has relative correction **O(1/p)**,
resolving the reflection's stated O(1/sqrt p) vs the close-out's "favors 1/p" IN FAVOR OF 1/p.

**Part C simplifies:** R(p)->e needs **NO Weil bound**. The only V=2k no-leaf pattern is the single
2k-cycle = `tr(M^{2k}) = (-p)^k(p-1)` (elementary); V<2k is trivially `O(p^{2k-1})=o(p^{2k})`.

### Impact
- THM-438 Part B mechanism CORRECTED (addendum added). Statements A_{2k}=C_k p^{k+1}, R(p)->e UNCHANGED.
- Part C upgraded: fully elementary (no Weil).
- Error term corrected p^{k+1/2} -> p^k; convergence rate of R->e pinned to 1/p (feeds HYP-2307 #2).

### Lesson
The project's own rule again: a clean final count (C_k) reached by a clean-sounding mechanism
("bigon-trees = Catalan") does not certify the mechanism. The Moebius weights and the
equal-order even-cycle patterns were invisible at the level of "count the bigon-trees." Always
decompose the inclusion-exclusion and check which patterns share the leading order — here the
cancellation `A088368 -> C_k` is the actual phenomenon, and it is the free-probability fingerprint
the moment-method slogan was pointing at.

---

## MISTAKE-061: THM-438 — the top-order patterns are NOT "even cacti"; they are the larger class of EVEN-SERIES patterns (even theta graphs included)

**Date discovered:** 2026-06-07
**Found by:** monad-explorer-2026-06-07 (deep-research / analytic lane, 4th session)
**Affects:** THM-438 ADDENDUM and MISTAKE-060 (the *characterization* of which coincidence
patterns reach the leading order `p^{k+1}`). Does NOT affect the Catalan law `A_{2k}=C_k p^{k+1}`
itself (re-confirmed here, rigorously) nor `R(p)->e`.

### What was assumed (MISTAKE-060 / THM-438 ADDENDUM)
"`M_sigma` reaches the top order `p^{k+1}` **iff** `F(sigma)` reaches full order `p^m` — exactly
the **even cacti** (connected, all biconnected blocks even cycles)." The census then grouped the
leading coefficient as bigon-trees (+A088368) corrected by even-cycle **cacti** down to `C_k`.

### Why it was wrong
`F(sigma) = sum_{flows} prod_e chi(t_e)` reaches full order `p^m` iff the flow-form product
`P(s) = prod_e ell_e(s)` is a **perfect square** (then `chi(P)=chi(Q^2)=+1` off the zero locus,
so `F ~ +p^m`). `P` is a perfect square iff **every series-class of edges has even size** (each
distinct flow-line occurs an even number of times). The even cacti satisfy this — but so do
**even theta graphs** (two vertices joined by three even paths; biconnected block is NOT a single
cycle) and, generally, all "even series-parallel" 2-connected patterns. These are NON-cacti yet
reach `p^{k+1}` and MUST be counted. Verified (`04-computation/paley_cluster_theta_check_monad.py`):
at `k=3` the `V=5, m=2` top-order patterns are **6 even-cacti{2,4} + 1 even-theta(2,2,2)** — the
even theta (mu=+1) was invisible to the "even cacti" census (it sat in the `(6,)` biconnected
bucket, silently cancelling the single 6-cycle, so the *total* still came out right).

### The correct framing (VERIFIED k<=4; the `g` step PROVED)
Let `c0 = lim A_{2k}/p^{k+1}`. Then
```
c0 = (-1)^k * sum_{rho : connected, EVERY series-class even}  mu(0,rho) * g(rho),
```
and `g(rho) := lim F(rho)/p^m = +1` for EVERY such pattern. **`g==+1` is PROVABLE:** within each
series-class the closed Euler walk passes straight through the degree-2 internal vertices, so all
edges of the class get the SAME orientation sign `s in {+1,-1}`; the class is even, so
`prod_{e in class} sign_e = s^{even} = +1`; hence `P = (prod sign_e) Q^2 = +Q^2` and
`g=chi(P)=+1`. Therefore the entire character/Gauss-sum content collapses and
```
($$)   sum_{rho : even-series pattern}  mu(0,rho)  =  (-1)^k C_k        (number-theory-FREE).
```
RIGOROUSLY CONFIRMED `c0 = 2, 5, 14 = C_2, C_3, C_4` by clean Richardson (`1/p`) extrapolation of
the exact flow-Moebius value (`04-computation/paley_cluster_topterm_monad.py`) — this also REPLACES
the prior slowly-converging census (which read `1.56, 2.77, 3.11` at `p<=23` and only *looked* like
it might reach `5`). The breakdown:
```
k=3:  bigon-trees(m=3) +13,  (m=2: cacti+theta) -9,  (m=1: 6-cycle) +1   = 5 = C_3
k=4:  bigon-trees(m=4) +69,  (m=3) -72,  (m=2) +18,  (m=1: 8-cycle) -1   = 14 = C_4
```
(bigon-tree sub-sums `+13, +69` = OEIS A088368, the all-pairings overcount, as before).

### Impact
- THM-438 ADDENDUM-2 added: Catalan law `A_{2k}=C_k p^{k+1}` RE-CONFIRMED (rigorous, k<=4), error
  `O(p^k)` unchanged, `R(p)->e` unchanged.
- The MECHANISM is corrected a SECOND time: top-order class = **even-series patterns** (perfect-square
  flow product), strictly larger than even cacti. `g==+1` is proved, reducing handoff #1 to the
  clean number-theory-free Moebius identity `($$)`.
- Free-probability reading SHARPENED: the random skew-Rademacher matrix gives `C_k` *directly* from
  non-crossing pairings (each `+1`, no factorials); the deterministic Paley Moebius expansion
  over-counts to A088368 in the bigon sector and the even cacti + even thetas + ... cancel it back to
  `C_k`. The equality is Wigner quasirandomness; `($$)` is its exact combinatorial fingerprint.

### Lesson
MISTAKE-060 corrected the *value* mechanism (bigon-trees -> A088368 -> C_k) but inherited a wrong
*support*: "even cacti." A pattern can saturate the flow character-sum without being a cactus — any
even series-parallel skeleton does. When a leading-order census gives the right TOTAL, that does not
certify the per-class STRUCTURE: a missing pattern (the even theta) can hide inside a coarse bucket,
cancelling against another, leaving the total correct and the story wrong. Characterize the support
by the actual saturation condition (perfect-square flow product = even series-classes), not by the
most familiar sub-family.

---

## MISTAKE-062: even-series pattern count is NOT OEIS A215257 — a 5-term coincidence that breaks at k=6

**Date discovered:** 2026-06-07 (monad-explorer, 8th session)
**Found by:** monad-explorer-2026-06-07 (deep-research, 8th session)
**Affects:** THM-438 ADDENDUM-3 point (2); HYP-2308; the reflection `the-drt-engine-is-S-squared-equals-J-minus-nI-the-catalan-is-genus-zero.md`; INDEX/SESSION-LOG entries asserting "even-series count = A215257"

### What was assumed
THM-438 ADDENDUM-3 (5th session) identified the number of EVEN-SERIES patterns of the
path `[0..2k]` (= the unsigned support of `(**)`) as **OEIS A215257**: the values for
`k=1..5` are `1, 3, 13, 67, 383 = A215257(k+1)` (indecomposable deque-sortable
permutations). The recursion script hardcoded the *predicted* next value `A215257(7)=2345`
for `k=6` but NEVER actually computed `k=6` (its `KMAX` default was 5).

### Why it was wrong
A direct exhaustive count at `k=6` (fast integer enumerator
`04-computation/paley_starstar_triangle_fast_monad.py`, cross-validated against the original
SVD test `04-computation/paley_starstar_crosscheck_monad.py` with **0 disagreements** over
all `Bell(13)=27.6M` partitions exhaustively at `k<=5` and a 300k sample at `k=6`) gives
```
   even-series count, k=1..6  =  1, 3, 13, 67, 383, 2351.
```
The OEIS b-file gives `A215257(7) = 2345 != 2351`. An OEIS search for
`1,3,13,67,383,2351` returns **no results** — the unsigned even-series count is not (yet)
a catalogued sequence. The A215257 match was a **5-term small-number coincidence**.

### The correct framing
- The unsigned even-series pattern count is `1, 3, 13, 67, 383, 2351, ...` (computed,
  rigorous through k=6), NOT A215257, and presently matches no OEIS sequence.
- This does NOT touch any headline result. The Moebius-SIGNED sum
  `(**) S_k = sum_{even-series} mu(0,sigma) = (-1)^k C_k` is independently re-verified
  exhaustively at `k=6` (`S_6 = 132 = C_6`), as is the cycle-rank triangle row
  `t(6,m) = 1, 45, 560, 2626, 4845, 2867` and the loop equation
  `S_k = -sum_{i+j=k-1} S_i S_j`.
- If anything the refutation SHARPENS the thread's thesis: the *unsigned* count is so
  unstructured it is not even a known sequence, while the *signed* sum is the cleanest
  possible (Catalan). "The Catalan is a cancellation, not a count" is now literal.

### Impact
- THM-438 ADDENDUM-3 (2) corrected (see ADDENDUM-6). HYP-2308 / INDEX A215257 cells updated.
- The "indecomposable deque-sortable permutations" bijection program (ADDENDUM-3/4 handoff)
  is moot — there is no A215257 bijection to find because the counts differ.

### Lesson
A 5-term OEIS hit is weak evidence — A215257 and the even-series count share five terms by
chance. NEVER hardcode an "expected next" OEIS value as if computed; compute it. Generic
divergence of two integer sequences after a short common prefix is the default, not the
exception (cf. MISTAKE-006 ratio coincidence, MISTAKE-010 small-n pattern break).
