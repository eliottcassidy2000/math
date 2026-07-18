---
id: THM-1026
title: THE FIVE-SLOT LEDGER OF LRC(14), AND WHERE THE PAIR FLOOR LANDS — pairs provably CANNOT close 13 runners (S₁ = 13·2λ = 13/7 ≈ 1.857, so Hunter must recover 6/7 ≈ 0.857 but 12 tree edges at ~1/49 recover only ~0.245); the level-5 truncation B₅ = 1 − S₁ + S₂ − S₃ + S₄ − S₅ is what LRC(14) actually needs, it clears at equidistribution (+0.1221), and it requires LOWER bounds on S₂, S₄ and UPPER bounds on S₃, S₅; THM-1025's all-regime pair floor fills the S₂ slot to 99.9% of the equidistributed value (median 1.5896 vs 1.5918, verified ≤ exact throughout) — one slot of five, and the other three are open
status: verified exactly — 60 comparable 13-packets, S₂ floor vs exact vs equidistributed; the Hunter-budget arithmetic is decisive; NO completed proof of LRC(14) is claimed or implied
source: opus-2026-07-17-S358 (owner: work the remaining open mathematics vigorously; formalize only if a proof is complete)
depends_on: [THM-1025 (the all-regime pair floor this places), THM-1021 (the negative that motivated it), THM-964/897 (the Bonferroni ladder), the fleet's B5 / level-5 machinery]
scripts: 04-computation/bonf5_pairfloor_opus_S358.py -> 05-knowledge/results/bonf5_pairfloor_opus_S358.out
---

# THM-1026 — the five-slot ledger, and an honest placement

**Pairs cannot close 13 runners.** At λ = 1/14 each comb has measure
2λ = 1/7, so S₁ = 13/7 ≈ 1.857 and the union bound overshoots by 6/7.
Hunter's second-order correction along a spanning tree recovers at most
12 edge-overlaps; at the independence value 1/49 each that is ≈ 0.245,
against a deficit of 0.857. The gap is structural, not a matter of
sharper constants: **the 7-wall closed because 7·(1/7) = 1 exactly, and
no such coincidence is available at 13.**

**What LRC(14) actually needs** is the odd level-5 truncation

> uncovered ≥ B₅ = 1 − S₁ + S₂ − S₃ + S₄ − S₅,

which at equidistributed values S_k = C(13,k)(1/7)^k gives **+0.1221 > 0**
— the level-5 wall clears, exactly as the ladder (THM-897) predicted.
Because the signs alternate, certifying B₅ > 0 requires **lower** bounds
on S₂ and S₄ and **upper** bounds on S₃ and S₅.

**Where the pair floor lands.** THM-1025's floor is precisely an S₂
lower bound, and it is valid in every regime. Measured on 60 comparable
13-packets (each floor verified ≤ the exact value):

| quantity | range | median |
|---|---|---|
| S₂ floor (THM-1025) | 1.5570 – 1.6414 | 1.5896 |
| S₂ exact | 1.5857 – 1.6757 | 1.5923 |
| S₂ equidistributed | — | 1.5918 |

The floor recovers **99.9%** of the equidistributed S₂. That slot is
effectively closed.

**The honest ledger.**

| slot | need | status |
|---|---|---|
| S₁ | fixed | 13/7, exact |
| S₂ | LOWER | **supplied** (THM-1025; kernel-pure for coprime pairs) |
| S₃ | UPPER | open (the fleet's T2 / THM-926 addresses the triple deviation) |
| S₄ | LOWER | open — the quadruple analogue of THM-1025 |
| S₅ | UPPER | open — the fleet's level-5 / B5 machinery |

**No completed proof is claimed.** The pair floor closes one slot of
five. Saying otherwise would be false, and the arithmetic above is the
reason: pairs are provably insufficient at 13 runners, and three of the
five slots have no bound yet.
