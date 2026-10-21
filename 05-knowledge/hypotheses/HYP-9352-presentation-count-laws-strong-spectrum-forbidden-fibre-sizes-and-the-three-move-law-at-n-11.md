---
id: HYP-9352
title: "After THM-4741: tc >= 5 for every strong tournament other than K1, C3, T5 (so tc = 3 classes are exactly the one-T5 stacks for all n); fibre sizes 7 and 21 never occur; the three-move law at n = 11 (random sample); growth of the strong spectrum minimum; the pendant nodes of the fundamental black graph (converse pairs of non-palindromic {K1,C3}-stacks plus the self-converse classes with tc - bm = 2) and their count law"
status: "PARTIAL: (1) basic strong floor/tc=3 classification PROVED; (4) pendant classification PROVED, quasi-polynomial expectation REFUTED; sharper spectrum, forbidden fibres, detours and descent OPEN"
session: opus-2026-10-20-S48
source: 05-knowledge/results/presentation_counts_single_view_classes_rotation_inequality_20261020.md
related:
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 05-knowledge/hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md
  - 05-knowledge/hypotheses/HYP-9346-self-lines-connectivity-and-growth-laws-of-the-tiling-line-metagraph.md
  - 05-knowledge/results/black_descent_20261020c/cycles.md
---

# HYP-9352: presentation-count laws after THM-4741

**Current scoped update, 2026-10-21.** The independently audited
[three-presentation classification](../results/metagraph_classification_20261021a/README.md)
proves the basic floor in (1): the only strong `tc=3` type is the cyclic
five-vertex tournament, and all ordinary `tc=3` types are ordered
`{K1,C3,T5}` stacks with exactly one `T5`. Their series is
`z^5/(1-z-z^3)^2`. The stronger monotonicity and `2n-1` clauses remain OPEN.

The pendant-center classification in (4) was already proved by the audited
[deficit-two center packet](../results/deficit_two_centers_20261020/README.md):
the only strong self-converse center of deficit two is `S4`. Thus its
self-converse pendant series is `z^4/(1-z^2-z^6)`, and the complete pendant
series is `(2z^4-z^5-z^7)/((1-z-z^3)(1-z^2-z^6))`. The expectation that every
fixed-degree count is eventually quasi-polynomial is REFUTED: mirroring
arbitrary words of `k` blocks chosen from `(K1,C3)` and `(C3,K1)` gives at
least `2^k` distinct degree-zero nodes of order `8k`. The repaired direction
is rational generating functions with possible exponential growth. This
repaired fixed-degree statement is now **PROVED** for every fixed degree in
the fundamental converse-merged paired-black graph by the audited
[fixed-degree rationality packet](../results/fixed_degree_metagraph_rationality_20261021/README.md).
The [final Session22 overview](../results/session22_tournament_tail_synthesis_20261021/SESSION22_SYNTHESIS.md)
records the scope through incoming `359ddd20` plus Session22 results: the degree-three formula still contains unknown
strong deficit-six center coefficients at seven possible orders
`{15,18,20,21,24,28,35}`, with degree at most35.
This does not complete adjacency classification. The remaining sharper
spectrum, forbidden-fibre, detour, descent, wiring, bridge and connectivity
questions remain **OPEN**.

The numbered original proposals below are retained for lineage; the scoped
update overrides their older status phrases.

1. **Strong spectrum floor.** Every strong tournament other than `K1`, `C3`, `T5` has presentation count `tc >= 5`; equivalently `H(T) > 3 |Aut(T)|` for strong `T` of order `>= 6`. With the stack multiplicativity (THM-2450 section 2) this makes the `tc = 3` classes exactly the stacks over `{K1, C3, T5}` with one `T5` for every `n` (exhaustive at `n <= 9`: `1, 2, 3, 6, 11`). Sharper: the minimum of the strong spectrum `min S(n)` (`1, 5, 3, 9, 9, 17, 25` at `n = 3..9`) is non-decreasing for `n >= 5` and at least `2n - 1` for `n >= 8`. What THM-4741 (C) gives for `tc = 3`: along every Hamiltonian cycle `T` is a circulant, or `|Z_C| = n/3`, or `n = 2 mod 4` with `|Z_C| = n/2`. Hostile: a strong tournament with a large automorphism group and few Hamiltonian paths (the wreath products `C3[C3[...]]` have `|Aut| = 3^((n-1)/2)` but far more paths).

2. **Forbidden fibre sizes.** No node of the line metagraph has `tc = 7` or `tc = 21` at any order (THM-1370 proves it for rigid nodes, where `tc = H`; a non-rigid counterexample needs `H = 7 |Aut|` or `21 |Aut|` with `|Aut| >= 3`, e.g. `H = 35, 49, 63`). Exhaustive at `n <= 9`. Companion (corrected by the audit): at `n = 9` every odd value below 2883 other than 7 and 21 occurs as a presentation count, with gaps only between 2883 and the maximum 3299; conjecture: for every `n` the presentation-count spectrum is an initial segment of the odd numbers minus `{7, 21}` followed by a sparse tail.

3. **Three-move law at `n = 11`** (HYP-9347 (1)): every closed blue line between non-pure-blue classes has a literal three-black-move detour at `n = 11`. Evidence: 50,190 random closed blue lines (drawn with replacement, about 80 duplicates), all layers, all with three-move detours and no exception (THM-4741's note, section 6); the full census (`2^24` lines, about 160,000 core-seconds) was running at close-out.

4. **Pendant nodes of the fundamental black graph** (corrects HYP-9347 (4), which called them "exactly the non-self-converse single-presentation classes"; MISTAKE-648). The merged black graph `B` of THM-4731 has degree `tc` at a merged non-self-converse node and `(tc - bm)/2` at a self-converse node, so its pendant nodes are (i) the converse pairs of non-palindromic `{K1, C3}`-stacks, `(a(n) - p(n))/2` of them (`a` = A000930, `p` the palindromic compositions into parts 1 and 3), and (ii) the self-converse classes with `tc - bm = 2`; at `n <= 9` the latter are exactly the two `S4`-centred palindromic stacks `(1, S4, 1)` (`n = 6`) and `(1, 1, S4, 1, 1)` (`n = 8`), both with `(H, |Aut|, tc, bm) = (5, 1, 5, 3)`, so the pendant counts are `1, 3, 3, 6, 7` at `n = 5..9`. Conjecture: the self-converse pendant nodes are exactly the palindromic stacks with a central `S4` and `{K1, C3}` wings (`tc - bm = 2` because `S4`'s five tilings contain three blue ones), and the degree-`d` nodes of `B` for fixed `d` are eventually a quasi-polynomial in `n`; the degree-0 nodes are exactly the pure-blue classes (`tc = bm`).

5. **Descent potentials.** `black_descent_20261020c/cycles.md` proves that the cyclic-triangle count descends along a black line from every strong black-active node; THM-4741 (D) shows that `H` has no such property. Conjecture: every non-pure-blue class with `tc >= 3` has a black neighbour with strictly smaller `tc`, except the strong classes with `tc = min S(n)`; untested (the census column for `tc`-descent is the next cheap addition).

Non-consequences: (1)'s strong floor does not bound `tc` for arbitrary non-strong classes (products); (3)'s listed evidence is a sample, not a certificate. The former beyond-order-nine uncertainty in (4)(ii) is superseded by the audited center classification linked above. Full metagraph adjacency, black connectivity, and the remaining spectrum clauses are not settled here.
