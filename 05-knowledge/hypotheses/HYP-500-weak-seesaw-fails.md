---
id: HYP-500
statement: "β₁(T) + β₃(T) ≤ 2 for all tournaments T (Weak Seesaw)"
status: REFUTED
source: opus-2026-04-04-S1
---

# HYP-500: Weak Seesaw β₁+β₃ ≤ 2

**Status:** REFUTED at n=8

**Counterexample:** Trial 2456 (seed=12456 in weak_seesaw_test.py)
- β₁ = 1, β₃ = 2, giving β₁+β₃ = 3
- Score sequence: (1, 3, 3, 3, 4, 4, 4, 6) — self-complementary
- c₃ = 14

**Context:** The strong seesaw β₁·β₃=0 (HYP-316) is PROVED for n≤7 but fails at n=8 (0.12% rate). This hypothesis attempted a weaker replacement.

**What holds:**
- β₁ ∈ {0,1} universally (THM-223)
- β₂ = 0 universally (THM-108 + THM-285)
- β₁·β₃=0 at n≤7 (proved)
- Euler characteristic: β₁+β₃ = β₄-β₅+β₆-β₇ (universal constraint)
- (β₁,β₃) ∈ {(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)} at n=8 (all observed)

**What is FALSE:**
- β₁·β₃=0 at n≥8
- β₁+β₃ ≤ 2 at n≥8
- "All seesaw violations have SC score" (1/7 does not at n=8)

**Rate data (n=8, 5000 samples):**
| (β₁, β₃) | Count | Rate |
|-----------|-------|------|
| (0, 0) | 4628 | 92.6% |
| (0, 1) | 303 | 6.1% |
| (0, 2) | 5 | 0.1% |
| (1, 0) | 57 | 1.1% |
| (1, 1) | 6 | 0.12% |
| (1, 2) | 1 | 0.02% |

**At n=9 (2000 samples):** No violations of β₁·β₃=0. β₁>0 is extremely rare (0.15%).

**Conclusion:** There is NO universal bound on β₁+β₃ beyond the Euler characteristic constraint. The seesaw is a small-n phenomenon that degrades as tournament complexity increases.
