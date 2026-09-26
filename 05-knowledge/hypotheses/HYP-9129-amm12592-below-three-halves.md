---
id: HYP-9129
title: "AMM 12592: the owner's question -- is C* < 3/2? (the optimal handoff limit lies in [1.3775, 1.570])"
status: >
  OPEN. Proved: C* >= 11/8 (THM-4467; reach of the Polya method 1.3775) and
  C* <= C_* = 1.598 (constructions). Numerically, the best super-blocks
  reach about 1.570. Deciding C* < 3/2 is a weighted-potential optimisation
  over handoff zero-measures mu (mass <= (B-2)/2, S(0) = +-1): minimise the
  asymptotic mountain-pass fold threshold. Candidate sharpness question: is
  C* = 1 + gamma_1 = 1.3775...? This is probably false; it would need
  integer spine series with capacity exactly 1 on W_(gamma_1).
source: collatz-procgen-20260923 AMM lane (candidates H2, H3), owner's question 2026-09-23
depends_on:
  - 01-canon/theorems/THM-4467-uniform-polya-capacity-gap-amm12592.md
  - 05-knowledge/hypotheses/HYP-9128-amm12592-superblocks-beat-golden.md
---

# HYP-9129 -- below 3/2?

**The question.** A deterministic fair extractor with constant `1 + eps`
is impossible for `eps <= 3/8` (THM-4467). Is a constant below `3/2`
possible?

**State.**
* Lower bounds: `C* >= 1.375` proved, `1.3775` numerically.
* Upper bounds: `C* <= 1.598` proved, about `1.570` numerically, via
  HYP-9128.
* The two sides meet in a potential-theoretic optimisation: the best
  handoff measure against the Pólya capacity of the quotient domain.

## 2026-09-23 update
* The proved window is now `1.377 <= C* <= 1.59`: THM-4467 extended, and THM-4468.
* Numerically, realizable handoff states reach `1.567`.
* The unrestricted real zero-measure relaxation reaches about `1.50`. But its zeros sit on an arc of `|w| = 1`, and an integer discriminant/Kronecker/parity argument (only sketched) excludes that for integer states.
* The necessary-condition LP stays feasible down to `1.435` (`B = 4`), so it cannot prove `C* >= 3/2`.

## 2026-09-26 update (opus, THM-4488)
* The proved window is `1.377 <= C* <= 197/125 = 1.576`: THM-4468's bottom-regime majorant replaced by the exact binomial ratio, everything else re-run. The certificate of this family stops at about `1.575` (level rates), consistent with the measure-level threshold `1.578`.

Still OPEN.
