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
