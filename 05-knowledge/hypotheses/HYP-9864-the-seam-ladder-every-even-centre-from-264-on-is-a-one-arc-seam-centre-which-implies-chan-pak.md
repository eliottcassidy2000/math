---
id: HYP-9864
title: "The seam ladder. Call an even number c a rung if some tournament T and arc a -> b have H(T) = c - 1 and H(T with the arc reversed) = c + 1 (a one-arc seam at centre c, THM-5267; equivalently a one-contact seam one order down). Conjecture: every even c >= 264 is a rung. Consequence (PROVED, trivial induction from the count 23): the conjecture implies Chan-Pak completeness, every odd number >= 23 being a count; the two solo separators 7 and 21 are the only forced jumps of the ladder (5 -> 9 and 19 -> 23), and below 264 the non-rungs other than 6, 8, 20, 22 are 36, 40, 54, 56, 60, 64, 84, 94, 108, 148, 188, 196, 262, each jumped by other realisers."
status: >
  OPEN (klein-2026-11-16-micro). FINITE-EXACT evidence: THM-5267 (C) — from bases of order <= 9 every even c in
  [264, 10000] is a one-contact seam centre, hence a rung; the misses below 264 are exactly the seventeen listed. The
  implication to Chan-Pak is proved in METAGRAPH.md section 2 of the packet. The conjecture is a density statement
  about the antisymmetric part of the pair matrix M(b, a) - M(a, b) (THM-5267 (A)); no structural mechanism is known,
  and the near-join analogue (THM-5265) is likewise density-only. A proof would need, for every large even c, a tournament
  of count c - 1 carrying an arc whose reversal adds exactly 2.
session: klein-2026-11-16-micro
date: 2026-11-16
source: 05-knowledge/results/microcosm_squaresum32_metagraph_selector_20261116/METAGRAPH.md
related:
  - 01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md
  - 01-canon/theorems/THM-5265-the-open-half-of-the-converse-is-a-density-statement-extended-near-join-family-same-core-obstruction-15-37-reach-398129-and-63-divides-mersenne-iff-p-1-mod-6.md
  - 05-knowledge/hypotheses/HYP-9453-the-near-join-sweep-toward-chan-pak-completeness-density-of-b-sets-two-strong-realisers-for-every-attained-odd-count-at-least-15.md
  - 05-knowledge/hypotheses/HYP-9863-heaviest-deletions-of-first-creators-are-first-creators-and-the-finitely-many-even-centres-without-a-one-contact-seam.md
---
# HYP-9864 — the seam ladder

A rung at `c` is an edge of the flip metagraph across which the count potential steps by exactly 2, from `c − 1` to
`c + 1`. If `c − 1` is a count and `c` is a rung, `c + 1` is a count; 23 is a count; so rungs at every even `c ≥ 24` give
every odd number from 23 on. The holes force the two jumps `5 → 9` and `19 → 23`; the thirteen other small non-rungs are
jumped because their neighbours are counts by other realisers. From 264 to 10,000 every even number is a rung
(THM-5267 (C)). The conjecture that this continues is equivalent in spirit, not in statement, to HYP-9453's density of
response sets: both ask that a one-parameter arithmetic family (here `M(b, a) − M(a, b) = ±2`) never leaves a gap.
