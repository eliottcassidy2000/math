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

The central open problem. Verified exhaustively for n≤6 (196,608 pairs) and by random sampling for n=7 (3,500 pairs, 0 failures). Claim B (the algebraic companion with I(Ω,2) instead of H(T)) is proved. The gap is the combinatorial identity connecting H(T) to I(Ω(T),2).

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

## OPEN-Q-009 🔴
**Prove arc-reversal invariance: D(T,v) = D(T',v) for flips not involving v**

The key unproved step for a general proof of Claim A. Equivalently: prove OCF (H(T) = I(Ω(T),2)) by showing ΔH = ΔI for each arc flip from transitive.

**MAJOR PROGRESS (opus-S2, THM-013):** At n=6, derived explicit formula:
  ΔI = -2·Σ_x s_x·H(B_x) + 2·(D5-C5)
where s_x = 1-T[x][i]-T[j][x], B_x = V\{i,j,x}. Verified 2216/2216 flips.

At n=5: even simpler: Ω(T) is always complete, ΔI = 2·(#destroyed - #created). Verified 732/732.

**Remaining task:** Prove the adjacency identity:
  adj(i,j) - adj'(j,i) = -2·Σ_x s_x·H(B_x) + 2·(D5-C5)

Key structural facts already proved:
- All destroyed/created cycles have {i,j} ⊆ V(C) (complement unchanged)
- THM-012: mu invariant when ≥1 flip endpoint in V(C)\{v}
- Complement adjacency formula for persist-change cycles (1896/1896 at n=6)

**Source:** file.txt exploration, opus-2026-03-05-S2

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

CONJ-002 refuted for p=11: H(T_11) = 95095 = 55*1729, not 4455. Known values:
- p=3: H=3, H/|Aut|=1
- p=7: H=189, H/|Aut|=9=3^2
- p=11: H=95095, H/|Aut|=1729=7*13*19

The sequence 1, 9, 1729 has no obvious pattern. Is there a number-theoretic formula fitting these values?

**Source:** kind-pasteur-2026-03-05-S2 computation

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery. Independently confirmed by opus-S2 with explicit mu triviality bound L >= n-2.
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
- **Paley computation (kind-pasteur)**: h_QR=h_NQR=201, c_9(T_11)=11055, H(T_11)=95095. CONJ-002 refuted for p=11 (kind-pasteur-2026-03-05-S2). Code: 03-artifacts/code/compute_H_T11.py.
- **OPEN-Q-011**: Near-cancellation is statistical, not structural. Not a viable proof strategy. (opus-S2, confirmed kind-pasteur n=7)
