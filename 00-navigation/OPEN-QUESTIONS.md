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

The key unproved step for a general proof of Claim A. D(T,v) = H(T) - H(T-v) - 2*sum_{C through v} H(T[V\V(C)]). If D is invariant under arc flips not involving v, then D=0 follows from a base case (transitive tournament). The naive cycle bijection FAILS (MISTAKE-005). A sum-level argument is needed: show sum_{C: i->j in C} H(T[V\V(C)]) = sum_{C': j->i in C'} H(T'[V\V(C')]).

**Source:** file.txt exploration

---

## OPEN-Q-010 🟡
**Per-path formula including 3-cycles AND 5-cycles at n=7**

At n=7, mu(5-cycle) = 1 always (V\{v + 4 cycle vertices} has 2 vertices, no odd cycles). So 5-cycle contributions are "trivially weighted" just like 3-cycles at n<=5. A per-path formula summing over both 3-cycle and 5-cycle embeddings (each with their mu weights) might work at n=7. Test computationally.

**Source:** FINAL_FINDINGS.md, Q3

---

## OPEN-Q-011 🟡
**Near-cancellation of two error effects at n=6**

Define A = sum_P' #TypeII(P'), B = sum_P' sum_{TypeII at j} mu(v,P'[j],P'[j+1]), D = sum_{all odd C through v} mu(C). Then A - B approx -5.88 (mu>1 inflation) and B - D approx +5.88 (5-cycle contributions). These are NEARLY EQUAL AND OPPOSITE. Is there an exact algebraic identity relating (A-B) + (B-D) = 0 to the structure of 5-cycles? Can the proof of Claim A be decomposed into proving A-B = -(B-D)?

**Source:** FINAL_FINDINGS.md

---

## OPEN-Q-012 🟢
**Tower hypothesis: L-cycle corrections from (L+2)-cycles**

At n=2k, the first cycle whose mu can exceed 1 has length 2k-1. The "excess" mu from shorter cycles may be exactly compensated by contributions from cycles 2 vertices longer. Is there a recursive structure where L-cycle corrections are expressed in terms of (L+2)-cycle contributions, creating a tower that sums to Claim A?

**Source:** FINAL_FINDINGS.md, Q5

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery. Independently confirmed by opus-S2 with explicit mu triviality bound L >= n-2.
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
