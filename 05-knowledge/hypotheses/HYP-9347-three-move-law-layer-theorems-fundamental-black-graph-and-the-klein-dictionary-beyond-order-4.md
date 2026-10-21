---
id: HYP-9347
title: "After THM-4732: the three-move law (every closed blue line between non-pure-blue classes has black distance exactly 3 for n >= 8), the (n-5,4) layer and a uniform layer theorem, bounded black distance (<= 5 for all n, the two-arc family F_n as the only distance-5 lines), the fundamental black graph's pendant nodes, and the chirality character on black twin pairs"
status: "NUMERICAL (n <= 10 exhaustive for (1), (3); n <= 9 for (4)); (2) OPEN (structural); (5) OPEN"
session: opus-2026-10-20-S47
source: 05-knowledge/results/black_detours_closed_blue_lines_klein_dictionary_20261020.md
related:
  - 01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
  - 05-knowledge/hypotheses/HYP-9346-self-lines-connectivity-and-growth-laws-of-the-tiling-line-metagraph.md
  - 05-knowledge/results/endpoint_three_repairs_20261020h/README.md
---

# HYP-9347: the three-move law and what follows it

Conventions as in THM-4732 (literal black detour, black distance, layer `k`, pure-blue = PB).

1. **Three-move law.** For every `n >= 8`, every closed blue line between two non-PB classes has black distance exactly 3 (a chain of three black K moves to its literal labelled partner). Evidence: exhaustive at `n = 8, 9, 10` (2048, 32768, 524288 lines), every layer; the endpoint lanes' structural theorems prove it for the layers `k = 1, 2` (`n >= 6`; at `n = 5` layer 1 holds `T5`'s pure-blue-touching line), `k = 3` (`n >= 11`, and `n = 10` by two concurrent certificates) and most of `k = 4` for `n >= 15` (endpoint_four_repairs_20261021). Hostiles: `F_5` and `F_7` (distance 5) show the law fails below `n = 8`. Next test: `n = 11` (`2^24` closed blue presentations, 55-bit tilings; feasible layer by layer in C++).

2. **Layer theorems.** (a) The `(n-5, 4)` layer admits a structural three-move proof for every `n >= 8` (the concurrent endpoint_four_repairs_20261021 proves `n >= 15` except the closed middle with `I = {y,c,h}`, `c -> x` or `x -> c` with `h in {x, n-4}`; `F_7` sits in that sub-case at `n = 7`, and in THM-4734's bank-miss pattern `{2, n-3, n-2}`; THM-4732 (4) guarantees the layer is PB-free for `n >= 8`; mac-mini's HYP-9348 names the same structural obligation; the open sub-case at `n = 11` is FINITE-EXACT: all 196,608 presentations have three-move detours, THM-4732 (3)); (b) a uniform proof over all layers exists in which the only orders needing finite control are those where the open member is path-poor (THM-4732's `F_n` phenomenon: black-move degree `6, 20` at `n = 5, 7`), i.e. the thresholds `n >= 11` of the lanes are regularity-coincidence artefacts and the true threshold of the three-move law is `n = 8`.

3. **Bounded distance.** For every `n` and every closed blue line between non-PB classes the black distance is at most 5, and equals 5 only for `F_5`, `F_7`. (Exhaustive to `n = 10`; `F_9`, `F_11` have distance 3.) A weaker form: the black distance is bounded by an absolute constant.

4. **Fundamental black graph.** The degree of a merged node of THM-4731's `B` is `tc` (non-self-converse) or `(tc - bm)/2` (self-converse), so its pendant nodes are exactly the non-self-converse single-presentation classes (`tc = 1`, i.e. `H = |Aut|`; THM-2444's rigidity lemma makes the self-converse ones pure-blue, of degree 0): PROVED, trivial. FINITE-EXACT at `n = 5..9`: the `tc = 1` classes are exactly the stacks over `{K1, C3}` (`sum_j C(n-2j, j) = 4, 6, 9, 13, 19` classes, all with `H = |Aut| = 3^j`), the palindromic ones pure-blue and the others pendant in `B` (14 at `n = 9`: six with `H = 3`, eight with `H = 9`); the unique degree-2 node at `n = 9` is the self-converse class `(H, |Aut|, tc, bm) = (15, 3, 5, 1)`. Conjecture: for every `n` the single-presentation classes are exactly the stacks over `{K1, C3}` (so `T5` never has `tc = 1` in any stack, and the pendant nodes of `B` are the non-palindromic such stacks), and `B` minus the PB nodes is connected for all `n` (equivalent to HYP-9346 (2)).

5. **Chirality beyond order 4.** THM-4732 (5)'s chirality character (compare `t` with `sigma t`; a black line takes the chirality of its closed member) splits every merged black multiplicity `2B` into `B + B` (a line and its sigma-twin have opposite chirality: PROVED, trivial from the sigma-twin theorem); conjecture: the two chirality halves of the black graph on classes are isomorphic as multigraphs for every `n` (true at `n <= 8` by the twin bijection; the content is whether a chirality-preserving relabelling exists at the class level beyond the twin map). Also: for `N = p^k q r` with odd `k` an equivariant embedding of `D(N)` into the order-5 tilings can be chosen to read `sign f` as closedness and `sign g` as chirality (untested beyond the eight-element case of `k = 2`).

Non-consequences: none of (1)-(5) implies black-core connectivity for all `n` without (a) of THM-4732 Lemma 2; (3) does not bound the length of class-level black walks through other classes (those can be shorter). Hostile for (2b): any closed blue line at `n = 11` without a three-move detour would refute (1) and the `n = 8` threshold at once.
