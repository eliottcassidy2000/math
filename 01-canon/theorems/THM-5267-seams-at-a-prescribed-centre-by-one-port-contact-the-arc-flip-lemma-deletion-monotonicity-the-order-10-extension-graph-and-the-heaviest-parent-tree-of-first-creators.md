---
id: THM-5267
title: "Seams and the plus-or-minus-one digits. (A) Arc-flip lemma: reversing the arc a -> b of a tournament T changes the count by M(b, a) - M(a, b), where M(y, x) is the number of ordered pairs of nonempty paths covering T with the first ending at y and the second starting at x (THM-5266); sum over the arcs a -> b of M(a, b) = (n - 1) h; a contact seam (T + v_X, T + v_{X u x}) is the arc-flip seam of the arc between v and x one order up, so the contact digit (B_x - e_x) - sum_{y in X} (M(y, x) + M(x, y)) and the flip digit M(b, a) - M(a, b) are the symmetric and antisymmetric parts of one matrix. (B) Deletion monotonicity H(T - v) <= H(T) (insertion of v at the front, the end, or before the first vertex it beats); hence every class is the heaviest parent of its sink extension. (C) One-contact seams from bases of orders <= 9 place every even centre c in [4, 10000] as counts c - 1, c + 1 except 6, 8, 20, 22, 36, 40, 54, 56, 60, 64, 84, 94, 108, 148, 188, 196, 262 (all below 263; none congruent to 2 mod 8); the +2 window (3H + 2, 3H + 4) around 3(H + 1) is a contact seam for 492 of the 498 odd counts H <= 1000 realised as contact values, failing only at H = 1, 11, 17, 19, 27, 35; the contact digits take every value in [0, 3638]; only 2 of 980 first-witness seams change the number of strong components. (D) Order-10 extension graph: 9,733,056 classes, 95,441,218 contacts, mean 9.806 distinct parents of 10, 8,150,261 classes with ten distinct deletion classes and 3 with one; children per order-9 class from 23 to 503 of 512; the transitive M_9 has 475, not the maximum, but is the largest branch of the least-string tree, in which 78,498 of 191,536 order-9 classes (41.0 %) are childless. (E) Heaviest-parent tree: no childless class at any order; every heaviest deletion (all ties included) of every first creator of orders 4-9 is a first creator, with the single exception of the converse pair creating 63 (heaviest deletions of count 37); branching 3, 8, 13, 44, 115, 247. (F) Moon's interval chains have H(IC(p_1, ..., p_k)) = product of (2^{p_i - 1} + 1), so Busch's strong floors 3, 5, 9, 15, 25, 45, 75, 125 are products of the Fermat-type numbers 3, 5, 9."
status: >
  PROVED (klein-2026-11-16-seam): (A) by counting the Hamiltonian paths through an arc as pair covers (replayed on 10,566
  flips of orders 3-7, 0 mismatches); (B) by the insertion argument (0 violations in 198,944 classes); (F) KNOWN (Moon
  1972, as cited in THM-4744), checked on all 33 compositions of N <= 8. FINITE-EXACT: (C) (seams.py by the quadratic
  form of THM-5266 (A) on orders <= 8; seams9.c on orders <= 9, centres <= 10,000), (D) (extgraph10.c, nauty labelg,
  86 s), (E) (heavy_tree.py, heavy_ties_ic.py on the extension graph of orders <= 9). (C)'s residue pattern and (E)'s
  genealogy are conjectured for all orders in HYP-9863. Not independently audited. CITED: THM-5266, THM-5264, THM-5260,
  THM-4922, THM-4744, THM-4889, A000568, center_126_new_interface_20261116, ready_restarts_rounding_20261116.
session: klein-2026-11-16-seam
date: 2026-11-16
source: 05-knowledge/results/seams_plusminus_ext10_20261116/SEAMS.md
scripts:
  - 05-knowledge/results/seams_plusminus_ext10_20261116/scripts/seams.py (+ receipts/seams.out)
  - 05-knowledge/results/seams_plusminus_ext10_20261116/scripts/seams9.c (+ receipts/seams9.out)
  - 05-knowledge/results/seams_plusminus_ext10_20261116/scripts/extgraph10.c (+ receipts/extgraph10.out)
  - 05-knowledge/results/seams_plusminus_ext10_20261116/scripts/heavy_tree.py (+ receipts/heavy_tree.out)
  - 05-knowledge/results/seams_plusminus_ext10_20261116/scripts/heavy_ties_ic.py (+ receipts/heavy_ties_ic.out)
related:
  - 01-canon/theorems/THM-5266-the-extension-grammar-curried-one-vertex-extension-count-is-a-quadratic-form-the-additive-sub-grammar-never-outputs-7-or-21-extension-graph-and-toothpick-tree-to-order-9-first-creators-are-85-percent-with-degree-blind-genealogy.md
  - 01-canon/theorems/THM-4922-first-creators-of-hamiltonian-path-counts-converse-fibers-endpoint-genealogy.md
  - 01-canon/theorems/THM-5260-forcing-the-midpoint-but-not-the-plus-2-child-history-blind-gadgets-are-ordinal-sums-the-ordinal-sum-laws-for-root-response-end-count-and-charge-and-the-tower-criterion-for-3h-plus-2.md
  - 01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md
  - 05-knowledge/hypotheses/HYP-9863-heaviest-deletions-of-first-creators-are-first-creators-and-the-finitely-many-even-centres-without-a-one-contact-seam.md
  - 05-knowledge/results/seams_plusminus_ext10_20261116/EXT10.md
  - 05-knowledge/results/seams_plusminus_ext10_20261116/PLUSMINUS.md
---
# THM-5267 — seams, the arc-flip lemma, and the heaviest-parent tree

**(A)** A Hamiltonian path of `T` through `a → b` is a pair (path ending at `a`, path starting at `b`) covering `V`, so there
are `M(a, b)` of them; after the flip the paths through `b → a` are the `M(b, a)` pairs (path ending at `b`, path starting at
`a`), which never use the arc between `a` and `b`; all other paths survive. Summing `M(a, b)` over the arcs counts each
Hamiltonian path once per arc. ∎ The contact statement is the flip of the arc `v–x` in `T + v_{X∪x}`.

**(B)** For a Hamiltonian path `p_1 … p_{n−1}` of `T − v`: if `v → p_1` insert `v` first; else if `p_{n−1} → v` insert it
last; else `p_1 → v` and `v → p_{n−1}`, and before the first `p_i` with `v → p_i` we have `p_{i−1} → v → p_i`. Distinct
paths give distinct extensions. ∎

**(C), (D), (E)** are censuses (SEAMS.md §2, EXT10.md §1–2). **(F)** The back arc of block `i` contributes the factor
`2^{p_i − 1} + 1` (Moon); checked for every composition of `N ≤ 8`.

The suffix word moves a count along a pencil `h + 2^k B`, the quotient multiplies it, and the port contacts move it by
the quadratic form; only the last places an arbitrary centre (every even number in `[264, 10000]` from bases ≤ 9). The
centres 6, 8, 20, 22 are forbidden by the holes; the other thirteen misses are small and unexplained. Chan–Pak OPEN.
