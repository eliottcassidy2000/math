---
id: THM-5271
title: "The metagraph, the microcosm and the selector. (A) On the flip metagraph of order n the Hamiltonian count is a potential whose edge differences are the antisymmetric pair-cover numbers M(b, a) - M(a, b) (THM-5267 (A)), even everywhere, with H = 1 only at the transitive corner; if every even c >= 24 is a rung (a one-arc seam centre: classes of counts c - 1 and c + 1 joined by a flip) then every odd number >= 23 is a Hamiltonian-path count (induction from 23); the rungs are the one-contact seam centres, all even c in [264, 10000] are rungs from bases of order <= 9, and the only forced jumps of the ladder are the holes 7 and 21 (HYP-9864). (B) The first 32 odd numbers split 3 | 1 | 6 | 1 | 21 around the holes with all block sizes triangular, least orders 1, 3, 4 | - | 5, 5, 5, 5, 6, 6 | - | 6 to 8 (63 last, at 8), and sum 32^2; the tripling chain 7, 21, 63, 189, 567, 1701 is a count from 63 on; every odd number in [65, 2881] is a count of order <= 9; the near-join misses of bases <= 9 above 2^19 are flat across residue classes mod 64, 8 and 7 (classes 7 and 21 mod 64 at relative rates 1.01 and 0.99) and 5 % lower on multiples of 3. (C) The square-sum graph on 1..n has Hamiltonian paths 1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 392, 669, 4041 for n = 15..34 and Hamiltonian cycles 0 for n <= 31, exactly 1 at n = 32 and 33, 11 at 34: the first 32 nodes carry the first cycle. (D) Image selector: the 1,200 counts first born at order 9 are values of the quadratic forms of 31 classes of order 8 (greedy), the best class (count 393) giving 226, the six maximal-count classes ranking about 4,843 of 6,880, and the first creators of order 8 covering all; near-join values from all rooted classes of orders <= 10 realise every odd number <= 19,190,825 except 1, 7, 11, 13, 21, 47, 67 (first miss 19,190,827), so the composite selector (census below 2,882, near-joins above) is exact for every odd N <= 19,190,825 except 7 and 21."
status: >
  PROVED (klein-2026-11-16-micro): (A) from THM-5267 (A) and the trivial induction; the ladder facts from THM-5267 (C).
  FINITE-EXACT: (B) (microcosm.py, residues.py on the order-9 census and bitmap), (C) (hamcyc.c), (D) (cover.py on
  profiles8.npz; nj10.c over all 9,733,056 classes of order 10 plus the smaller orders). ANALOGY, typed: the square-sum
  and count microcosms. Not independently audited. CITED: THM-5263, THM-5264, THM-5265, THM-5266, THM-5267, THM-5268,
  THM-4094, THM-4505, HYP-9453, HYP-9863, HYP-9864, THM-282, THM-1965, THM-4730 (the metagraph), A071983.
session: klein-2026-11-16-micro
date: 2026-11-16
source: 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/METAGRAPH.md
scripts:
  - 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/scripts/nj10.c (+ receipts/nj10.out)
  - 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/scripts/cover.py (+ receipts/cover.out)
  - 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/scripts/hamcyc.c (+ receipts/hamcyc.out)
  - 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/scripts/microcosm.py (+ receipts/microcosm.out)
  - 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/scripts/residues.py (+ receipts/residues9.out, residues10.out)
related:
  - 05-knowledge/hypotheses/HYP-9864-the-seam-ladder-every-even-centre-from-264-on-is-a-one-arc-seam-centre-which-implies-chan-pak.md
  - 01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md
  - 01-canon/theorems/THM-5265-the-open-half-of-the-converse-is-a-density-statement-extended-near-join-family-same-core-obstruction-15-37-reach-398129-and-63-divides-mersenne-iff-p-1-mod-6.md
  - 01-canon/theorems/THM-282-blue-lives-on-spine.md
  - 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/SELECTOR.md
---
# THM-5271 — the metagraph as a potential, the microcosm, the selector

**(A)** `H` is odd (Rédei) and changes by `M(b, a) − M(a, b)` along a one-arc flip, an even number; `H = 1` iff the tournament
is transitive. A rung at `c` carries a count from `c − 1` to `c + 1`; 23 is a count; induction. The forced jumps are at
the holes: no rung at 6, 8, 20, 22. ∎ The rest of (A), and (B)–(D), are the censuses of MICROCOSM.md, SQUARESUM32.md,
METAGRAPH.md and SELECTOR.md with their receipts.

The holes do not recur in the macrocosm as residue classes; they recur only as the two forced jumps of the seam ladder
and as the two members `7·1, 7·3` of the lane `7p`. Chan–Pak OPEN.
