---
id: THM-5266
title: "The extension grammar of tournament classes, curried. (A) For a tournament T with count h, end-counts e_x, root responses B_x (THM-4865) and pair-cover numbers M(y, x) (ordered pairs of nonempty paths covering T, the first ending at y, the second starting at x), a new vertex v with out-neighbourhood X has H(T + v_X) = h + sum_{x in X} (B_x - e_x) - sum_{y != x in X} M(y, x): the count of a one-vertex extension is a quadratic form on the cube, and S(n+1) is the union of the images of the forms of the order-n classes. (B) The uniform contacts - sink (B doubles, e vanishes, B_s = h), source (B_src = P_T(2)), sole exit ((h, e, B) -> (h - e + B, h - e, h + B)), ordinal sum (THM-5260 (B)) - are closed on (h, e, B, P) and compose to the near-join eps_t o sigma^k (h + 2^k B_t), the tower and the midpoint of THM-5264; a general contact needs M, two of them the 3-path tensor. (C) The additive sub-grammar (near-joins, towers, dead-end appends) never outputs 7 or 21: the response sets are B(1) = {2^j}, B(3) = {3 2^j}, B(5) = {5, 6, 7, 9, 32} 2^j (tournaments of count 1, 3, 5 are padded chains, padded 3-cycles, padded strong order-4), 12 is not in B(9), and B >= h; so 7 is not 2^m + 1 and 21 is not 3(2^m + 1) - the holes are the Mersenne shadows of 9 and 27. The full grammar encodes them through multiplicativity and Busch's floors (THM-1370). (D) Extension graph to order 9 (nauty canonical forms): contacts 2, 7, 29, 215, 2591, 50357, 1656619 at orders 3-9; mean distinct parents 1, 1.75, 2.42, 3.84, 5.68, 7.32, 8.65 (of n); children of an order-(n-1) class from 22 to 248 at n = 9; the transitive tournament is the most fertile only through order 4; in the canonical-parent tree the branching at 8 -> 9 has maximum 227 (the transitive trunk) and 2,977 of 6,880 order-8 classes are childless; the toothpick diagrams draw one toothpick per class crossing its least canonical deletion. (E) First creators (THM-4922) are 1, 1, 6, 33, 328, 5457, 163676 of the classes of orders 3-9 (85 % at order 9, 1,520 values); the out-degree of the deleted vertex over the genealogy edges is close to binomial with mode at the middle (35176, 144462, 296548, 372259, ... at order 9), so the genealogy is nearly degree-blind; first creators with a digit-0 vertex are 21 %, 7.3 %, 1.7 %, 0.37 % at orders 6-9; of the orphan pair creating 63 one member has a digit-0 vertex (a dead-end append over an order-7 base) and its converse none."
status: >
  PROVED (klein-2026-11-16-tooth): (A) by splitting a Hamiltonian path of T + v at v (replayed on all 62,422 extensions
  of orders <= 7, 0 mismatches); (B) from THM-5260 (B), THM-4865, THM-4862 (a), MISTAKE-776; (C) by the classification
  of the tournaments of count 1, 3, 5 (counts multiply over strong components; Busch's floors 5 and 9 at orders 4, 5), the
  ordinal-sum laws, and the census of the two strong order-5 classes of count 9 (replayed over all rooted tournaments of
  orders <= 8 with h < 21). FINITE-EXACT: (D), (E) (extgraph.py with labelg; creators9.c; spectrum9.c and spectrum9z.c
  agree class by class). Not independently audited. CITED: THM-4922, THM-4889, THM-1370, THM-5260, THM-4865, THM-5264,
  THM-5265, A000568, A002785, A139250.
session: klein-2026-11-16-tooth
date: 2026-11-16
source: 05-knowledge/results/toothpick_grammar_creators_20261116/GRAMMAR.md
scripts:
  - 05-knowledge/results/toothpick_grammar_creators_20261116/scripts/quadform.py (+ receipts/quadform.out)
  - 05-knowledge/results/toothpick_grammar_creators_20261116/scripts/holes_grammar.py (+ receipts/holes_grammar.out)
  - 05-knowledge/results/toothpick_grammar_creators_20261116/scripts/extgraph.py (+ receipts/extgraph.out)
  - 05-knowledge/results/toothpick_grammar_creators_20261116/scripts/creators9.c (+ receipts/creators9.out)
  - 05-knowledge/results/toothpick_grammar_creators_20261116/scripts/tree_creators.py (+ receipts/tree_creators.out)
  - 05-knowledge/results/toothpick_grammar_creators_20261116/scripts/toothpick.py (+ figures/toothpick_7.svg, toothpick_8.svg)
related:
  - 01-canon/theorems/THM-4922-first-creators-of-hamiltonian-path-counts-converse-fibers-endpoint-genealogy.md
  - 01-canon/theorems/THM-5264-the-out-degree-1-law-for-digit-0-vertices-is-proved-at-every-order-two-path-lemma-unique-digit-0-vertex-in-a-strong-tournament-and-the-converse-digit-0-implies-a-near-join-count-below-22228.md
  - 01-canon/theorems/THM-5265-the-open-half-of-the-converse-is-a-density-statement-extended-near-join-family-same-core-obstruction-15-37-reach-398129-and-63-divides-mersenne-iff-p-1-mod-6.md
  - 01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci.md
  - 05-knowledge/results/toothpick_grammar_creators_20261116/TOOTHPICK.md
  - 05-knowledge/results/toothpick_grammar_creators_20261116/CREATORS.md
---
# THM-5266 — the extension grammar, curried

**(A)** A Hamiltonian path of `T + v_X` has `v` first (`Σ_{x∈X} s_x` paths), last (`Σ_{y∉X} e_y`), or inside as `P_1 v P_2`
with `last(P_1) = y ∉ X`, `first(P_2) = x ∈ X` (`Σ_{y∉X, x∈X} M(y, x)`). With `Σ_{y∉X} e_y = h − Σ_{y∈X} e_y`,
`Σ_{y∉X, x∈X} M = Σ_{x∈X} Σ_y M(y, x) − Σ_{y, x∈X} M(y, x)` and `B_x = s_x + Σ_y M(y, x)` the form follows. ∎

**(C)** Counts multiply over strong components and strong tournaments of order ≥ 4 (≥ 5) have count ≥ 5 (≥ 9), so the
tournaments of count 1, 3, 5 are `M_n`, `M_a ⇒ C_3 ⇒ M_b`, `M_a ⇒ S_4 ⇒ M_b`. `B_t(M_n) = 2^{n−1−i}`, `B_t(C_3) = 3`,
`B(S_4) = {5, 6, 7, 9}`, `P_{C_3}(2) = 12`, `P_{S_4}(2) = 32`, and the ordinal-sum laws give the three response sets. For
`k ≥ 1`, `7 − h ∈ {6, 4, 2}` and `21 − h ∈ {20, 18, 16}` for `h = 1, 3, 5` are not in `2^k 𝔅(h)`; `h ≥ 9` gives `h + 2B ≥ 27`.
For `k = 0` (dead-end appends) the only candidate is `21 = 9 + 12`, and the responses of the count-9 tournaments (padded
`C_3 ⇒ C_3`: `9·2^j, 36·2^j`; the two strong order-5 classes: odd parts `1, 5, 9, 11, 15, 21, 23, 25`) exclude 12. ∎

**(B), (D), (E)** are in GRAMMAR.md §2, TOOTHPICK.md and CREATORS.md with their receipts. The grammar of first creators is
the full extension grammar, not its additive part: the holes are encoded intrinsically only in the additive sub-grammar
(C), and in the full grammar by THM-1370's floors. Chan–Pak OPEN.
