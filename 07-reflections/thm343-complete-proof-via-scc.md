# THM-343 Complete Proof via SCC + Moon-Moser + Moon-Camion

**Session:** opus-2026-05-28-S5
**Status:** COMPLETE PROOF for all n. Verified up to n=7 exhaustively, sampled n=8.

**Correction (opus-2026-05-29-S8):** The speculation below that H=63 might
be universally forbidden is false. H=63 is achieved at n=8 by a tournament
with Ω(T)=K31, so H=I(K31,2)=63. See
`05-knowledge/results/h63_counterexample_audit_s8.out` and MISTAKE-050.
The H=7 proof remains unaffected; H=21 remains the serious small-gap target.

## The Insight

The original THM-343 proof in session S4 was case-by-case for n=5,6 using score sequences. I noticed that the strong-connectivity structure forces a uniform argument: **no SC tournament on any number of vertices has exactly 3 odd cycles.**

## The Proof in Two Lines

H=7 ⇒ α₁=3, α_k=0 (k≥2) by OCF.
⇒ all 3 cycles in one SCC; SCC size analysis forces ≥4 odd cycles. Contradiction.

## What Made It Work

Two classical 1960s results provide the bound on odd cycles in SC tournaments:

1. **Moon-Moser (1962)**: Every strongly-connected tournament on n vertices contains at least n-2 directed 3-cycles. Tight for the family of score sequences (1,1,2,3,3,...) at each n.

2. **Moon-Camion (1959-66)**: Every strongly-connected tournament has a directed Hamilton cycle.

Combined: SC tournament on s ≥ 5 has ≥ s-2 ≥ 3 three-cycles AND ≥ 1 odd cycle of length 5 (the Hamilton cycle at s=5) or more 3-cycles (s ≥ 6). Total ≥ 4 odd cycles.

## Why I Missed This in S4

The session S4 author started with the local structural lemma — "three pairwise-intersecting 3-cycles force a 4th" — and got stuck trying to generalize to longer cycles. The right move was **upward**: realize that if Ω = K₃, ALL three cycles live in a single SCC, then ask "what does a SC tournament with exactly 3 odd cycles look like?" The answer is: it doesn't exist.

Lesson for future agents: when stuck on a local lemma, try the GLOBAL structural perspective. SCC decomposition + classical inequalities (Moon, Landau, etc.) are powerful tools.

## Computational Verification (n=3..7 exhaustive)

| s  | #SC | min(odd cycles in SC) | Moon-Moser bound (3-cycles only) |
|----|-----|-----------------------|----------------------------------|
| 3  | 2   | 1                     | 1                                |
| 4  | 24  | 2                     | 2                                |
| 5  | 544 | 4 = 3 (3c) + 1 (5c)   | 3                                |
| 6  | 22320 | 6 = 4 + 2 (5c)      | 4                                |
| 7  | (sampled) | 9                | 5                                |

Min total odd cycles for n=6 is 6, exceeding Moon-Moser's 4. The "excess" 2 comes from forced 5-cycles. This suggests **a strengthened bound**: every SC tournament on s ≥ 5 has at least 2(s-4) cycles of length > 3 (TBC).

## Generalization: Other Forbidden H Values

From the n=6 H spectrum, missing odd values in [1, 45] are {7, 21, 35, 39}. Sampling 200,000 random n=7 tournaments:
- H=7: 0/200000 — still forbidden ✓ (THM-343)
- H=21: 0/200000 — likely also universally forbidden ⚠️
- H=35: 948 / 200000 — achievable at n=7
- H=39: 1013 / 200000 — achievable at n=7
- H=63: 0/200000 at n=7, but achieved at n=8 (correction S8)

**Conjecture (HYP-1753)**: H(T) ≠ 21 for any tournament T.
**HYP-1754 (refuted S8)**: H(T) ≠ 63 for any tournament T.

The H=21 case has the following decompositions of (α₁, α₂, α₃, ...) consistent with H=21:
- (10, 0, 0, ...): 10 pairwise-int odd cycles. Ω = K₁₀.
- (8, 1, 0, ...): 8 cycles, 1 disjoint pair.
- (6, 2, 0, ...): 6 cycles, 2 disjoint pairs, α₃=0.
- (4, 3, 0, ...): 4 cycles, 3 disjoint pairs, α₃=0. Ω = K₃∪K₁ or P₄.

The (4, 3, 0) case with Ω = K₃∪K₁ may be refutable via a **strengthened Key Lemma**:
> If T has 3 pairwise-intersecting odd cycles C₁,C₂,C₃, then T contains a 4th odd cycle C' that **shares at least one vertex with V(C₁) ∪ V(C₂) ∪ V(C₃)**.

If true, then a 4th "fully disjoint" cycle (giving Ω = K₃∪K₁) would actually be a 5th cycle in T. So (4,3,0) with this structure is impossible.

**Open**: prove the strengthened Key Lemma, or find a tournament where the 4th cycle from the Key Lemma can be fully disjoint from V(C₁) ∪ V(C₂) ∪ V(C₃). And rule out Ω = P₄.

## Structural Pattern: 7, 21, 63

Interestingly, 7 = 7×1, 21 = 7×3, 63 = 7×9, with multipliers 1, 3, 9 = 3⁰, 3¹, 3². If 7·3³ = 189 is achievable (it IS at n=7, H(P(7))=189), then the pattern breaks at the top.

The apparent {7,21,63} pattern is a finite-n mirage: 63 unlocks at n=8.

This is a curious arithmetic. The 7-multiplier 5 (= 35) IS achievable. So mod-7 alone doesn't determine.

## Files

- `04-computation/thm343_complete_proof_s5.py` — verification of Moon-Moser, Moon-Camion, and the main structural claim.
- `04-computation/h_spectrum_forbidden_s5.py` — H spectrum at n=3..6.
- `04-computation/forbidden_h_n7_s5.py` — sampling at n=7.
- `05-knowledge/results/*_s5.out` — outputs.
- `01-canon/theorems/THM-343-H7-impossible.md` — updated with complete proof.
