---
id: THM-5268
title: "Neutral choices of the core. (A) Padding is transparent: for T = M_a => S => M_b with S strong of order >= 3, the balanced arcs of T (M(a, b) = M(b, a), whose reversal preserves the count, THM-5267 (A)) are the adjacent-rank arcs of the paddings, whose flips give isomorphic tournaments, and the balanced arcs of S, trivial or nontrivial exactly as in S; so the count-preserving moves of a fibre F_h are those of its strong realisers, and for h < 125 (Busch's floor f(10) = 125) the fibre's neutral structure is determined at orders <= 9. (B) The fibres with no nontrivial count-preserving flip at any order, for h < 125, are exactly 1, 3, 5, 11, 13, 17, 19, 31, 35, 39, 47, 53, 55, 57, 63, 67, 85, 95, 113, 115, 121; they include the permanent near-join exceptions 11, 13, 47, 67, the +2-window failures 11, 17, 19, 35, the atom 63 and the tower-blocked 53, 67, 113; 9 and 27 are not rigid. (C) Order 9 has 82,384 balanced arcs, 54,928 trivial and 27,456 nontrivial; the neutral-richest strong classes have counts 37 (order 6, seven balanced arcs), 51 and 33 (order 7, six), 657 (order 8, seven), 495 (order 9, ten); the fibre with the most nontrivial moves is 37 at orders <= 8 and 547 at order 9; Spearman's correlation between nontrivial neutral edges and the size of the response set B(h) over the 61 odd h <= 125 is 0.613. (D) The holes are the k = 1 near-joins of the rigid fibres 1 and 3 with the next fibre's response: 7 = 1 + 2*3 and 21 = 3 + 2*9 with 3 in B(3) but not in B(1) and 9 in B(9) but not in B(3); the pencils 1 + 2^k*3 and 3 + 2^k*9 miss only their first terms; in the lane 7p = p + 6p the holes are the two members whose shift 6p is twice a wrong-fibre response; the even prime 2 is the doubling both need and the centre of the first seam (1, 3)."
status: >
  PROVED (klein-2026-11-16-core): (A) by the path-insertion arguments of NEUTRAL.md section 1 and the canonicity of the
  condensation. FINITE-EXACT and, by (A) with Busch's floor, PERMANENT: (B) (neutral.py on the pair matrices of all
  classes of orders <= 8; neutral9.c + labelg on all 191,536 classes of order 9). FINITE-EXACT: (C). (D) is arithmetic
  on THM-5266 (C) and the census of counts to order 9 (holes_arith.out). Not independently audited. CITED: THM-5267,
  THM-5266, THM-5264, THM-5263, THM-5260, THM-4094 (the lanes p and 7p), THM-1370, Busch's floors, and the concurrent
  packets ready_core_growth_and_mod18_20261117 and one_port_neutral_continuations_20261117 (their "neutral" is a
  property of an appended continuation; here it is a property of the core).
session: klein-2026-11-16-core
date: 2026-11-16
source: 05-knowledge/results/holes_neutral_core_20261116/NEUTRAL.md
scripts:
  - 05-knowledge/results/holes_neutral_core_20261116/scripts/neutral.py (+ receipts/neutral.out)
  - 05-knowledge/results/holes_neutral_core_20261116/scripts/neutral9.c (+ receipts/neutral9.out)
  - 05-knowledge/results/holes_neutral_core_20261116/scripts/neutral9.py (+ receipts/neutral9_py.out)
  - 05-knowledge/results/holes_neutral_core_20261116/scripts/holes_arith.py (+ receipts/holes_arith.out)
related:
  - 01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md
  - 01-canon/theorems/THM-5266-the-extension-grammar-curried-one-vertex-extension-count-is-a-quadratic-form-the-additive-sub-grammar-never-outputs-7-or-21-extension-graph-and-toothpick-tree-to-order-9-first-creators-are-85-percent-with-degree-blind-genealogy.md
  - 01-canon/theorems/THM-5263-digit-0-vertices-backward-junctions-sole-exits-and-the-near-join-the-primes-11-13-47-67-are-permanent-near-join-exceptions.md
  - 01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md
  - 05-knowledge/results/holes_neutral_core_20261116/HOLES.md
  - 05-knowledge/results/holes_neutral_core_20261116/CORE.md
---
# THM-5268 — neutral choices of the core

**(A)** Inside a padding, a non-adjacent flip creates a strong block of count ≥ 3 and multiplies the count; an
adjacent-rank flip returns transitive padding. For `x` in `M_a` and `y ∈ S`: if `x` is not last in `M_a`, all old paths
survive and (`M_a` before `x`, a path of `S` ending at `y`, `x`, the rest of `M_a`, `M_b`) is new; if `x` is last, the
`s_y(S)` old paths through `x → y` are lost but the pairs (path of `S` ending at `y`, the rest of `S`) give `E_y(S) ≥ h(S)`
new ones, and `s_y(S) < h(S)` since `y` is not a source of a strong `S`. For `x ∈ M_a`, `z ∈ M_b`: (`M_a − x`, `S`, `M_b` to
`z`, `x`, the rest of `M_b`) is new and nothing is lost. Flips inside `S` give `M_a ⇒ S′ ⇒ M_b`, isomorphic to `T` iff
`S′ ≅ S`. ∎ Strong realisers of `h < 125` have order ≤ 9, so (B) is permanent.

**(B), (C)** are the censuses of NEUTRAL.md. **(D)** is the arithmetic of HOLES.md on THM-5266 (C): `𝔅(1) = {2^j}` and
`𝔅(3) = {3·2^j}` contain neither 3 nor 9, so the two `k = 1` near-joins with the next fibre's response are the holes.

The rigid fibres cannot be re-arranged without a free path run, and a run only doubles; targets in the lanes `p`, `7p`
should be written over flexible fibres (CORE.md). Chan–Pak OPEN.
