---
id: LEM-045
title: THE MULTI-PARENT (ANCESTOR-SET) HUNTER, WITH THE HONEST c = 9 NEGATIVE — and the mechanical remainder CASHED IN LEAN. (A) THE INEQUALITY: for arbitrary ancestor sets P(i) ⊆ {0..i−1}: μ(∪_{i<n} A_i) + Σ_{i≥1} μ(A_i ∩ ⋃_{j∈P(i)} A_j) ≤ Σ μ(A_i) — the same leaf-peeling disjointification (⋃_{P(n)} A_j ⊆ S); interpolates union bound (P = ∅) → tree/arborescence Hunter (|P| = 1) → the EXACT disjointification identity (P = all). Kernel-pure in Lean (multi_parent_hunter_add_le), Python-refereed exact on 25 random family/ancestor-set instances. (B) THE HONEST NEGATIVE: two-parent credits do NOT revive the consecutive c = 9 route — measured exactly: credits 0.2351 vs 2/7 = 0.2857 needed (v = 50); three-parent at c = 10: 0.3609 vs 0.4286 — consecutive triples are strongly correlated (the triple correction eats the union gain); the LEM-044 boundary stands. The tool remains available for MIXED blocks (a far parent adds nearly-independent credit) — per-family decide. (C) THE MECHANICAL REMAINDER CASHED (LRCTreeHunter.lean, all kernel-pure [propext, Classical.choice, Quot.sound]): consecutive_closed_form — the ALL-k induction DONE in-kernel (sum_shifted evaluation + the 686q² + 196qr + 98q + 98r polynomial collapse via linear_combination); good_pos_of_path_credits — the general assembly skeleton (n events of measure ≤ 1/7 with path credits C and n/7 < 1 + C leave a good set of positive measure; at n = 8 the hypothesis is exactly 1/7 < C = the c = 8 theorem's shape); lt_one_of_add_le (ENNReal cancellation). Residual named: the arc-measure rendering of LEM-042(A) (μ(D_a ∩ D_b) = the muNum value) to feed concrete credits into the skeleton
status: PROVED + LEAN kernel-pure (build green 8475 jobs) + Python referee (25 instances; the c = 9/10 negatives exact)
source: boxeph-2026-07-17-S73 (owner directive: the mechanical remainder + arborescences/hunter generalizations)
depends_on: [LEM-044 (the boundary being tested), klein path hunter (the p = pred instance), LRCTreeHunter.lean S72]
script: 04-computation/lrc14_multiparent_hunter_boxeph_S73.py -> 05-knowledge/results/lrc14_multiparent_hunter_boxeph_S73.out
---

# LEM-045 — the multi-parent hunter; the remainder cashed

One induction now covers the whole credit spectrum: union bound, paths,
stars, arborescences, ancestor-set (DAG-credit) unions, up to the exact
identity. The consecutive c = 9 hope dies honestly (correlated triples);
the all-k closed form and the assembly skeleton live in the kernel. What
separates the c = 8 theorem from a fully formal statement is exactly one
named item: rendering μ(D_a ∩ D_b) = muNum/(14ab) as a Lean measure fact.

## Evidence log
- [x] multi-parent inequality: Lean kernel-pure + 25 exact Python instances
- [x] c = 9 (2-parent) and c = 10 (2/3-parent) consecutive negatives, exact
- [x] all-k closed form in-kernel; assembly skeleton in-kernel
- [x] arc-measure rendering (S74) + cast bridge (S75) + THE PER-RUNNER BOUND (S76, danger_measure_le kernel-pure): μ(dangerR v ∩ window) ≤ 1/7 every v ≥ 1 (tooth cover Icc(−⌈v/2⌉,⌈v/2⌉); interior teeth 1/(7v); extreme teeth ≤ ofReal(1/2 − (M−ρ)/v) — 1/(14v) even, nonpositive odd). ALL c = 8 INGREDIENTS IN-KERNEL; the end-to-end composition (restrict-measure plumbing) is the one named final item
