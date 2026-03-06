# Open Questions

**Status codes:** 🔴 CRITICAL (blocks main proof) | 🟡 IMPORTANT (needed for paper) | 🟢 INTERESTING (worth exploring)

---

## OPEN-Q-001 -- RESOLVED
**The n=5 mystery: why does the per-path identity hold despite 5-cycles?**

**RESOLVED by THM-008:** The per-path identity holds trivially for n<=5 because mu(C) = 1 for ALL odd cycles C through v. For 3-cycles, the complement V\{v,a,b} has at most 2 vertices, which cannot form an odd cycle. For 5-cycles, C\{v} exhausts all of T-v, leaving 0 available vertices. The identity reduces to #TypeII = #TypeII. There is no "delicate balance" -- the identity is vacuous at n<=5.

**Additional detail (opus-S2):** More generally, mu(C) = 1 whenever cycle length L >= n-2 (THM-008 mu triviality bound). At n=6, mu(3-cycle) is in {1, 3}: mu=1 (76.7%) when 3 available vertices form transitive subtournament, mu=3 (23.3%) when cyclic.

---

## OPEN-Q-002 🔴
**Prove Claim A: H(T) − H(T−v) = 2Σ_{C∋v} μ(C)**

The central open problem. Equivalent to OCF: H(T) = I(Ω(T), 2). Claim B (the algebraic companion with I(Ω,2) instead of H(T)) is proved. The gap is the combinatorial identity connecting H(T) to I(Ω(T),2).

**Verification record (opus-2026-03-05-S2+S3):**
- n≤6: EXHAUSTIVE (33,864 tournaments, 0 failures) via verify_ocf_sweep.py
- n=7: EXHAUSTIVE (1,048,576 arc assignments, 0 failures) via symbolic_proof_n7.py (27min, opus-S3)
- n=8: 500 random, 0 failures
- n=9: 100 random, 0 failures (first verification at this n; max indep set size=3)
- n=10: 30 random, 0 failures (first verification at this n)

**Current approaches (from paper §claim_strategies):** Five routes documented. None yet complete.

---

## OPEN-Q-003 -- RESOLVED
**Characterize when the per-path identity holds at n=6**

**RESOLVED by THM-009:** The per-path identity fails for path P' iff some Type-II position (a,b) in P' has mu(v,a,b) > 1, which happens iff the 3 vertices V\{v,a,b} form a directed 3-cycle in T-v. This is a perfect binary separation: mu>1 at any TypeII position => always fails; all mu=1 => always holds.

---

## OPEN-Q-004 🟡
**Find a correct per-path formula for all n**

The 3-cycle-only formula (per-path identity) fails for n≥6. The natural generalization summing over all odd cycles overcounts. The maximal-embedding-only formula also fails. Is there a formula of the form Σ_{cycles C, (non-v consecutive in P')} f(C, P') = (inshat−1)/2 that works for all n?

---

## OPEN-Q-005 🟡
**Combinatorial proof of the C(L-2, 2k-1) distribution (THM-007)**

The number of internal signature patterns giving exactly k Type-II positions within an L-cycle embedding window is C(L-2, 2k-1). This is proved (verified computationally). A bijective proof connecting this to ballot sequences or Dyck paths is missing. See TANGENT T001.

---

## OPEN-Q-006 🟢
**Asymptotic formula for Σ_C μ(C)**

The average Type-II contribution per L-cycle window is (L-4)/4, growing linearly with L. Does this yield an asymptotic formula for Σ_C μ(C) as a function of the cycle-length distribution of T? What happens for random tournaments as n→∞?

---

## OPEN-Q-007 🟡
**Full proof of Fix(σ) = 2^{m²} for self-evacuating SYT**

Verified for n=5 and n=7 (m=2 and m=3 respectively, giving 4 and 512 self-evacuating SYT). Full proof is conditional on a precise classical reference not yet pinned down. The identification with TSSCPPs may provide the reference.

---

## OPEN-Q-008 🟢
**2-adic tower: what is the 2-adic valuation of H(T)?**

The OCF formula gives H(T) = I(Ω(T), 2). Evaluating at x=4,8,... gives higher-order mod-2^k invariants. Is there a combinatorial characterization of v_2(H(T)) (the 2-adic valuation)?

---

## OPEN-Q-009 🔴 → MAJOR PROGRESS (Claim B path identity PROVED)
**Prove arc-flip identity: E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips**

The key unproved step for a general proof of Claim A. Equivalently: prove OCF (H(T) = I(Omega(T),2)) by showing delta_H = delta_I for each arc flip from transitive.

**BREAKTHROUGH (opus-2026-03-05-S4/S4c):** The Hamiltonian path alternating sum identity
(THM-016, "Claim B path identity") is now PROVED by induction on m:
  sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) = (-1)^{m+1} h_end(W, v)
This was the LAST MISSING PIECE for the inductive proof of B(Li,Rj)=B(Lj,Ri).
THM-016 → THM-017: B(Li,Rj)=B(Lj,Ri) for ALL n — the even-odd split is PROVED.
**CAVEAT:** The even-odd split is necessary but NOT sufficient for OCF (kind-pasteur-S8).
OCF for all n still requires proving delta_H = delta_I (unsigned sum = cycle formula).
**Remaining:** Find a route from the even-odd split to OCF, or prove delta_H=delta_I directly.

**GENERAL FORMULA (opus-S2, THM-013; independently derived by kind-pasteur-S5 via A-clique):**

  delta_I = sum_{k>=1} 2^k * Delta(alpha_k)

where Delta(alpha_k) depends recursively on alpha_{k-1} of cycle complements. Equivalently (kind-pasteur-S5): delta_I = 2 * [sum_{gained C'} I(R_{C'}, 2) - sum_{lost C} I(R_C, 2)] where R_C = Omega(T[V\V(C)]). By strong induction, I(R_C, 2) = H(complement). This is the A-clique argument adapted to arc flips (same technique as Claim B proof). Verified n=4,...,9 (opus) and n=5 exhaustive + n=6 random (kind-pasteur).

**Independent verification (opus-2026-03-05-S3):** Adjacency identity adj(i,j)-adj'(j,i) = RHS verified:
- n=5: 10,240/10,240 (exhaustive, all tournaments x all arcs)
- n=6: 3,000/3,000 (random sampling)

**Simplified forms by n (opus-S2):**
- n<=5: delta_H = 2*sum_L(DL-CL) [simple formula]
- n=6,7: delta_H = -2*sum s_x*H(B_x) + 2*sum_{L>=5}(DL-CL)
- n=8: needs correction for VD 3-5 pairs (simplified formula fails 8/30)
- n=9+: needs alpha_3 terms (three mutually VD 3-cycles)

**EXHAUSTIVE PROOF FRONTIER (opus-S4, extending kind-pasteur-S6/S7):**
- n<=6: proved by kind-pasteur-S6 (symbolic polynomial identity, THM-015)
- n=7: PROVED by opus-S3/S4 (2^20 = 1,048,576 configs) and kind-pasteur-S7 (SymPy + enumeration)
- n=8: PROVED by opus-S4 (2^27 = 134,217,728 configs, 57 minutes, ALL passing)
  Also independently verified by opus-S4b C implementation (0 fails through 3M+ configs)
- n=9 would need 2^35 ~ 34B configs — infeasible without new approach

**Even-Odd Split Lemma (opus-S4/S4b, new discovery):**
The adj decomposition delta = sum_S Delta(S, V\{i,j}\S) satisfies:
  sum_{|S| even} Delta(S,R) = sum_{|S| odd} Delta(S,R)
i.e., the alternating sum vanishes: sum (-1)^|S| Delta(S,R) = 0.
Equivalently (Signed Position Identity, opus-S4b): sum_{P: i->j} (-1)^{pos(i)} = sum_{P': j->i} (-1)^{pos(j)}.
This is a POLYNOMIAL identity requiring T[a][b]+T[b][a]=1 (FAILS for general digraphs: 282/500 fail at n=4).
Verified n=4,...,8 exhaustive, continuous weighted tournaments pass.
**CAVEAT (kind-pasteur-S8):** This is a CONSEQUENCE of OCF, not equivalent to it. The odd-S sum
of Delta(S,R) differs from the cycle formula [g-l]*H(R) per-subset; see even-odd-split-lemma.md.
See 03-artifacts/drafts/even-odd-split-lemma.md.

**PROVED (kind-pasteur-S10):** The signed position identity / alternating sum identity / even-odd
split is now proved for ALL n, via Claim (B) (THM-016). See THM-016 for the inductive proof.

**Remaining task:** Bridge the gap between even-odd split and OCF. The even-odd split gives
delta_H = 2 * sum_{|S| odd} Delta(S,R), but this does NOT equal the cycle formula delta_I
(see MISTAKE-008). An additional identity is needed to connect these two expressions.

**Key structural facts (both agents independently confirmed):**
- All affected cycles contain {i,j} (complement unchanged by flip)
- At most one affected cycle in any independent set (A-clique)
- Complement values computable by inductive OCF
- The swap involution (opus THM-014) gives adj(i,j)-adj'(j,i) = #U_T - #U_T'
- Even-odd split: delta decomposes equally between even-S and odd-S terms (T040/T044)
- Bracket B(u,w) = T[u][i]*T[j][w] - T[u][j]*T[i][w] has block structure by s-types (T043)
- When all s_v=0: brackets AND boundary terms vanish (potential induction base)

See PROP-001 (kind-pasteur), THM-013, THM-014, THM-015 (both agents).

**Source:** opus-2026-03-05-S2/S3/S4/S4b; kind-pasteur-2026-03-05-S5/S6/S7

---

## OPEN-Q-010 🟡
**Per-path formula including 3-cycles AND 5-cycles at n=7**

At n=7, mu(5-cycle) = 1 always (V\{v + 4 cycle vertices} has 2 vertices, no odd cycles). So 5-cycle contributions are "trivially weighted" just like 3-cycles at n<=5. A per-path formula summing over both 3-cycle and 5-cycle embeddings (each with their mu weights) might work at n=7. Test computationally.

**Status (kind-pasteur-2026-03-05-S3):** NEGATIVE RESULT. The per-path formula does NOT simplify at n=7. The algebraic identity (inshat-1)/2 = #{TypeII} = #{3-cycle embeddings} is trivially satisfied, but this just restates THM-004+005 -- it does not encode 5-cycle information. Computing the actual A, B, D quantities (see test_n7_ABD.py) shows A=/=D in general. A=/=D means: total TypeII count (A) does not equal total odd-cycle mu sum (D). The 5-cycles contribute non-trivially even when mu=1. See T027 and OPEN-Q-011.

**Source:** FINAL_FINDINGS.md, Q3; kind-pasteur-2026-03-05-S3

---

## OPEN-Q-011 -- RESOLVED (statistical artifact, not structural)
**Near-cancellation of two error effects at n=6**

**Resolved by:** opus-2026-03-05-S2, confirmed by kind-pasteur at n=7

**Answer:** The near-cancellation is a statistical observation, NOT an exact identity. Computational verification (3000 pairs at n=6, opus-S2) shows:
- A = D exactly for only 836/3000 (28%) of pairs
- A - D ranges from -12 to +9 (mean ≈ 0)
- Mean(A-B) ≈ -Mean(B-D) is approximate, not exact

The decomposition A-B = -(B-D) does NOT hold pair-by-pair. The two effects cancel on average but not structurally. This is NOT a viable proof strategy for Claim A.

**Status (kind-pasteur-2026-03-05-S3):** PARTIAL ANSWER. At n=7, tested 1050 (T,v) pairs: mean A-D = 0.097 (near zero), but NOT zero in general (range -39 to 26). Mean A-B = -73.78, mean B-D = +73.88 (near-cancellation on average). The near-cancellation is STATISTICAL, not algebraic. The per-pair |A-D|<=1 holds only 13.1% of time. The decomposition Claim A = (A=B) + (B=D) does NOT yield two tractable sub-identities. The near-cancellation at n=6 was likely a low-n coincidence.

**Source:** FINAL_FINDINGS.md; kind-pasteur-2026-03-05-S3

---

## OPEN-Q-012 🟢
**Tower hypothesis: L-cycle corrections from (L+2)-cycles**

At n=2k, the first cycle whose mu can exceed 1 has length 2k-1. The "excess" mu from shorter cycles may be exactly compensated by contributions from cycles 2 vertices longer. Is there a recursive structure where L-cycle corrections are expressed in terms of (L+2)-cycle contributions, creating a tower that sums to Claim A?

**Source:** FINAL_FINDINGS.md, Q5

---

## OPEN-Q-013 🟡
**Correct formula for H(T_p) for Paley primes p ≡ 3 (mod 4)**

Both conjectures are FALSE for p=11:
- Original conjecture H(T_p) = p * 3^((p-1)(p-3)/8) gives 649,539 for p=11, not divisible by 55.
- Revised conjecture H(T_p) = |Aut(T_p)| * 3^((p-3)/2) gives 4455 for p=11 (off by factor 21.4).

**Known values (all confirmed by direct computation):**
- p=3: H=3, |Aut|=3, H/|Aut|=1
- p=7: H=189, |Aut|=21, H/|Aut|=9=3^2
- p=11: H=95095, |Aut|=55, H/|Aut|=1729=7*13*19 (Hardy-Ramanujan taxicab number)
- p=19: H=unknown, |Aut|=171, H/|Aut|=unknown

**Sequence:** 1, 9, 1729 — no obvious pattern. 3^k pattern breaks catastrophically at p=11.

**Complete cycle count table for T_11** (confirmed kind-pasteur-S5 from inbox/other.txt, all consistent with H=95095):
| k | c_k(T_11) | C(11,k) | c_k/C(11,k) | integer? |
|---|-----------|---------|-------------|----------|
| 3 | 55 | 165 | 1/3 | no |
| 4 | 165 | 330 | 1/2 | no |
| 5 | 594 | 462 | 9/7 | no |
| 6 | 1595 | 462 | 145/42 | no |
| 7 | 3960 | 330 | 12 | YES |
| 8 | 7425 | 165 | 45 | YES |
| 9 | 11055 | 55 | 201 | YES |
| 10 | 10681 | 11 | 971 | YES |
| 11 | 5505 | 1 | 5505 | YES |

**OCF verification:** 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155 EXACT

**Integrality observation:** C(11,k) | c_k(T_11) for ALL k >= 6 = (p+1)/2. This appears to be a genuine structural theorem; see T028.

**Next computational target:** H(T_19) (requires computing Ham paths of 19-vertex tournament). |Aut(T_19)|=171=9*19.

**Source:** kind-pasteur-2026-03-05-S2 computation; cycle table from kind-pasteur-2026-03-05-S5 (inbox/other.txt)

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery. Independently confirmed by opus-S2 with explicit mu triviality bound L >= n-2.
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
- **Paley computation (kind-pasteur)**: h_QR=h_NQR=201, c_9(T_11)=11055, H(T_11)=95095. CONJ-002 refuted for p=11 (kind-pasteur-2026-03-05-S2). Code: 03-artifacts/code/compute_H_T11.py.
- **OPEN-Q-011**: Near-cancellation is statistical, not structural. Not a viable proof strategy. (opus-S2, confirmed kind-pasteur n=7)
