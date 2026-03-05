# Open Questions

**Status codes:** 🔴 CRITICAL (blocks main proof) | 🟡 IMPORTANT (needed for paper) | 🟢 INTERESTING (worth exploring)

---

## OPEN-Q-001 🔴
**The n=5 mystery: why does the per-path identity hold despite 5-cycles?**

At n=5, 5-cycles through v DO exist (they use all 5 vertices: v→a→b→c→d→v requires 4 other vertices, which T-v has exactly). Yet the per-path identity (inshat−1)/2 = Σ_{3-cycle embeddings} μ(C) holds with 0 failures at n=5.

At n=6 it fails 2,758/9,126 times. The per-path identity only sums over 3-cycle embeddings (by the bijection in THM-005), but Claim A's RHS includes all odd-cycle contributions. At n=6 the 5-cycles add uncaptured weight. Why doesn't this break things at n=5?

**Hypothesis:** At n=5, every 5-cycle embedding is always "accompanied" by exactly the right set of 3-cycle embeddings to account for its μ contribution. This delicate balance breaks at n=6.

**What's needed:** A theorem explaining why the 5-cycle μ contribution is always "absorbed" by 3-cycles at n=5.

---

## OPEN-Q-002 🔴
**Prove Claim A: H(T) − H(T−v) = 2Σ_{C∋v} μ(C)**

The central open problem. Verified exhaustively for n≤6 (196,608 pairs) and by random sampling for n=7 (3,500 pairs, 0 failures). Claim B (the algebraic companion with I(Ω,2) instead of H(T)) is proved. The gap is the combinatorial identity connecting H(T) to I(Ω(T),2).

**Current approaches (from paper §claim_strategies):** Five routes documented. None yet complete.

---

## OPEN-Q-003 🟡
**Characterize when the per-path identity holds at n=6**

~70% of (T,v,P') triples at n=6 satisfy the per-path identity. ~30% fail. What structural property of T distinguishes the two groups? Is it the number of 5-cycles through v? The structure of Ω(T-v)?

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

## Resolved Questions (moved here when answered)

*None yet.*
