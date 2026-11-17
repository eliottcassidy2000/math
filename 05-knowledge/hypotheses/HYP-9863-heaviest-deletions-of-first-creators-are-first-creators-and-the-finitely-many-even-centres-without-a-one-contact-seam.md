---
id: HYP-9863
title: "(a) For every first creator F of order n >= 4 (a class of least order with its count, THM-4922) other than the two first creators of 63, every vertex deletion of F of maximum count is a first creator of order n - 1; so the first creators form a tree under the heaviest-deletion map with 63 the only break. (b) Only finitely many even centres c have no one-contact seam (a pair T + v_X, T + v_{X u x} with counts c - 1, c + 1): below 10,000 and from bases of order <= 9 they are exactly 6, 8, 20, 22, 36, 40, 54, 56, 60, 64, 84, 94, 108, 148, 188, 196, 262, and none of them is 2 mod 8; conjecturally these seventeen are all. (c) The +2 window (3H + 2, 3H + 4) is a one-contact seam for every odd count H other than 1, 11, 17, 19, 27, 35."
status: >
  OPEN (klein-2026-11-16-seam). FINITE-EXACT evidence: (a) all 169,501 first creators of orders 4-9, every tied heaviest
  deletion checked (heavy_ties_ic.out); (b) centres <= 10,000 from bases <= 9 (seams9.out), the 20 misses from bases <= 8
  show the same residues; (c) the 498 odd counts <= 1000 realised as contact values. No mechanism is known for the residue
  pattern in (b); 6, 8, 20, 22 are forbidden by the holes 7 and 21 (THM-1370). A proof of (a) would need an inequality
  tying a new count at order n to a new count at order n - 1 through the heaviest deletion; THM-4889's
  H(T) <= sum_v H(T - v) gives only the weak form in THM-4922 (2).
session: klein-2026-11-16-seam
date: 2026-11-16
source: 05-knowledge/results/seams_plusminus_ext10_20261116/EXT10.md
related:
  - 01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md
  - 01-canon/theorems/THM-4922-first-creators-of-hamiltonian-path-counts-converse-fibers-endpoint-genealogy.md
  - 05-knowledge/results/seams_plusminus_ext10_20261116/SEAMS.md
---
# HYP-9863 — the heaviest-deletion tree of first creators, and the seam centres

**(a)** THM-4922 proved that some deletion of a first creator is a first creator when `h > n(h)·max H(n(h) − 2)` and
found it true at orders 4–10 except for 63. The census here says more: the deletion of *largest* count is always one,
ties included, at orders 4–9, again except for 63 (whose heaviest deletions have count 37, born at order 6). Deletion
monotonicity `H(T − v) ≤ H(T)` (THM-5267 (B)) makes the heaviest-parent map a tree on all classes with no childless node;
(a) says the first creators are a subtree of it.

**(b), (c)** Seams move a count by one arc (THM-5267 (A)); the quadratic form of THM-5266 (A) realises every even centre in
`[264, 10000]` from bases of order ≤ 9. Whether the thirteen small misses other than the hole-neighbours are reached by
larger bases, and why no miss is `≡ 2 (mod 8)`, is open; the window failures `11, 17, 19, 27, 35` are the familiar hard
small counts (11 sterile-only, 17 and 27 fertile-only, 35 an atom).
