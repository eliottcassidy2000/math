---
id: HYP-9664
title: "Signatures, carriers of both signs and extremal sign patterns. (A) Carriers of both signs: for every odd d >= 3 there are two-port kernels with at most two visits of the forms (1, 2^u, 0, d 2^w) and (-1, -2^u, 0, d 2^w), so that the plus grid 1 + 2^s + d 2^(r+s) and the minus grid d 2^(r+s) - 2^s - 1 are path counts for all large r, s. (B) The share of hyperbolic kernels among the two-port kernels with at most two visits tends to a positive limit. (C) For every n >= 4 the Cayley-Dickson sign tournament minimises the number of even zero-sum 4-sets among the tournaments on PG(n-1, 2) with all lines cyclic, the minimum being 42 [n-1 choose 3]_2, and at most [n choose 3]_2 - 7 [n-1 choose 3]_2 planes are odd at once. (D) Module-side sources for a fixed grid exist at finitely many levels only: for fixed odd d there are finitely many s with a tournament X, H(X) = 2^s + 1 and odd part of Pi_2(X) dividing d."
status: >
  OPEN (all four parts). PRE-AUDIT DRAFT.
  (A) FINITE-EXACT evidence: plus carriers d = 3, 21 (order 8; KNOWN, OFFSET_TRANSPORT.md), d = 35 (KNOWN, same
  source), d = 43 (THM-5130 (F), a port of a unit root of order 9), d = 145 (unaudited, S22); minus carrier d = 31
  (THM-5133 (C), order 8); the census of all two-port kernels of order 9 adds none on either side beyond the known
  d = 35. No carrier for d = 9 among kernels of order at most 9 nor among unit roots of order at most 10 (S22). The hypothesis does not imply
  Chan-Pak: a carrier covers its grid only from some level and some r on. For d = 1 the ranges must avoid
  N+(1,1,1) = 7, N+(1,2,2) = 21 and N-(1,1,3) = 7.
  (B) FINITE-EXACT: the share is 0, 0, 2.5, 9.8, 18.7, 27.2, 34.8 per cent at orders 3..9 (THM-5133 (C)). No model.
  (C) FINITE-EXACT for n = 4 (THM-5133 (F), exhaustive linear algebra): the minimum is 42 and the bound on odd planes
  is 8 = 15 - 7. For n = 5, 6 only the values of CD_n are known (630 and 6510 even 4-sets; 50 odd planes of 155 at
  n = 5). Without the condition on lines the minimum at n = 4 is 21, not 42.
  (D) FINITE-EXACT evidence for d = 9: among all classes of order at most 9 only the levels s = 1, 2, 3, 7; among
  one-apex words of chain length at most 22 only s = 1, 2 (pre-registered predictions F1, F2 of the opposite sign
  FAILED). For closings with a transitive module the statement is an S-unit equation in 2 and 3 (finitely many
  solutions for each pair, CITED, not effective). Orders 10 and more are unexplored for s >= 7.
session: kind-pasteur-2026-11-11-S23
source: 05-knowledge/results/kps_s23_flip_spin_kernels_20261111/README.md
related:
  - 01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md
  - 05-knowledge/hypotheses/HYP-9662-two-sided-seeds-and-the-reach-of-sign-flip-towers.md
  - 05-knowledge/hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md
---

Notation as in THM-5133. A two-port kernel `(K; x, y)` with at most two visits has a form
`al h h' + be pi h' + ga h pi' + de pi pi'` in the coordinates `(h, pi) = (H, Pi_2)` of its two modules.

## (A) Carriers of both signs

Every odd `N` has a plus address `N = 1 + 2^s + d 2^(r+s)` (unless `N = 2^s + 1`) and a minus address
`N = d' 2^(r'+s') - 2^(s') - 1`; exactly one of the two has level 1. HYP-9662 (A) asks for a seed and an anti-seed
for every odd `c`, one line at a time. A carrier is a uniform family of them: one kernel, two free transitive
modules, a whole grid. The hypothesis is that the supply of carriers is complete on both sides. What would refute
it: a proof that some `d` has no carrier (for instance through the finiteness of the admissible visit matrices).
What would make it useful: carriers whose ranges start at level 2, since every odd `N` has one address of level at
least 2.

## (B) Signature statistics

The sign of `det G` separates kernels whose path count is a sum of two products from those where it is a
difference. The data do not yet distinguish a limit below 1/2 from 1/2.

## (C) Extremal sign patterns

On `PG(3, 2)` the sedenions are one of two extremal families (THM-5133 (F2): 120 sets of type "a plane and the
planes through a point off it", 15 of type "the planes off a point"). The hypothesis says that the doubling stays
extremal in every dimension. A first test is `n = 5`: decide by linear algebra over `F_2` whether more than 50 of
the 155 planes of `PG(4, 2)` can be odd at once.

## (D) Fermat modules

By THM-4756, `Pi_2 = (-1)^n H (mod 3)`, so for `d = 9` a module at an even level must have `Pi_2` a power of two.
The classes of order at most 9 with `Pi_2` a power of two have `H = 1, 5, 19, 25, 217, 223, 259, 599, ...`; none of
these is `2^s + 1` beyond `s = 2`.
