---
id: LEM-044
title: THE CONSECUTIVE CLOSED FORM, THE c = 8 THEOREM, THE EDGE CAP, AND THE TREE-HUNTER CEILING (the wall program past c = 7). (A) THE CONSECUTIVE CLOSED FORM: μ(D_k ∩ D_{k+1}) = 1/49 + r(6−r)/(49k(k+1)), r = k mod 7 — the LEM-042 integer sum collapses to 49Σ − 14k(k+1) = 14r(6−r) exactly; zero excess iff 7 | k, max at r = 3; specializes to 1/(7(k+1)) for k ≤ 6. Refereed exact k = 1..400. (B) THE c = 8 CONSECUTIVE THEOREM: every consecutive block [v..v+7] crosses the wall — the seven path credits sum to 1/7 + Σ_{r_i≠0} r_i(6−r_i)/(49(v+i)(v+i+1)), and a 7-window of k's contains EXACTLY ONE r = 0, so the excess is STRICTLY positive uniformly in v: good ≥ excess > 0 (referee: min excess 1.7e-5 at v = 200; three blocks vs true good). THE c = 9 BOUNDARY: the pure consecutive route needs excess > 6/49 and dies for ALL v (referee: no crossing v ≤ 100). (C) THE EDGE CAP: μ(D_a ∩ D_b) ≤ 1/14 for every distinct pair, equality iff reduced (1,2) — proof: each LEM-042 term ≤ 2a and there are ≤ (a+b)/7 + 1 terms, giving μ ≤ (a+b+7)/(49b) ≤ 1/14 outside the finite set a, b ≤ 4, checked directly; scan to 300 confirms uniqueness. (D) THE TREE-HUNTER AND ITS CEILING: Hunter's inequality holds for ANY spanning tree T (leaf-plucking induction; klein's kernel-pure Lean lemma is the path case, the tree case named): μ(∪A_i) ≤ Σμ − Σ_{T} μ(pairs) — refereed exact on 60 random family/tree instances. Crossing at block size c needs tree credits > (c−7)/7 = (2c−14)/14 while (C) caps them at (c−1)/14: the tree-hunter can cross ONLY for c ≤ 12 — AT c = 13 IT IS IMPOSSIBLE FOR EVERY TREE AND FAMILY (12/14 = 12/14, not strict). Feasibility is real: doubling families {3·2^i} cross at c = 9, 10, 12 with credits exactly (c−1)/14 (brute-verified against true good)
status: PROVED ((A) algebraic collapse, verified; (B) one-zero-per-window; (C) finite-check + termwise bound; (D) leaf-plucking + (C)) + REFEREED EXACT throughout
source: boxeph-2026-07-17-S71 (owner directive: more similar LRC progress; extends LEM-042/043, klein's hunter, opus's block reduction)
depends_on: [LEM-042 (formula), LEM-043 (cone floor/dichotomy), klein path_hunter_add_le]
script: 04-computation/lrc14_consecutive_treehunter_boxeph_S71.py -> 05-knowledge/results/lrc14_consecutive_treehunter_boxeph_S71.out
---

# LEM-044 — past the wall: c = 8 proved, the ceiling at 13

The wall program now has exact boundaries: c = 7 falls to the dichotomy
(LEM-043), c = 8 falls to the consecutive closed form (every consecutive
8-block, uniformly), c = 9..12 fall per-family to tree-hunter credits
(decide; doubling families demonstrate), and c = 13 is provably beyond ANY
tree-hunter — the edge cap 1/14 meets the requirement (2c−14)/14 exactly at
c = 13. Larger blocks need non-pair structure (the deep-well/citation
routes), which is where the program's other tools already live.

## Evidence log
- [x] closed form k ≤ 400; c = 8 theorem + boundary; edge cap ≤ 300 + finite proof
- [x] tree-hunter 60 instances; doubling demos c = 9, 10, 12; ceiling at 13
- [ ] named: Lean — the tree case of hunter (leaf-plucking), the closed form
      (integer decide per k), the c = 8 window argument (7-residue pigeonhole)
