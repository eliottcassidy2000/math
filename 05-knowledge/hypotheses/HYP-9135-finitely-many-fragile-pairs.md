---
id: HYP-9135
title: "Collatz has only finitely many fragile pairs: flipping one pair {2i-1, 2i} creates a new cycle only for the 24 known i <= 2308"
status: >
  OPEN HYPOTHESIS; FINITE-EXACT to i <= 10^7 (lane), and independently to
  i <= 20000 (orchestrator). In the pairing family of THM-4470, flipping
  pair i sends 2i-1 down to i-1 and 2i up to 3i. A new cycle appears
  exactly for i = 1, 4, 5, 10, 11, 13, 20, 22, 40, 61, 84, 122, 126, 167,
  189, 217, 244, 325, 334, 433, 445, 577, 1154 and 2308. Equivalently,
  either the T-orbit of 3i contains 2i, or the orbit of i-1 contains 2i-1
  (with the gluing described in THM-4470 §6). This is a cycle-gate
  statement: i (2^(K+1) - 3^(a+1)) = B_w.
source: collatz-procgen-20260922 session, brackets/pairings lane candidate (P3), 2026-09-24
depends_on:
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
---

# HYP-9135 -- finitely many fragile pairs

**Statement.** The set of `i` such that the one-flip mutant of Collatz at
pair `i` has a cycle other than `{1,2}` (or `{0}` when `i = 1`) is finite.
The strong form is that it is exactly the 24 values listed above.

**Refutation form.** A fragile pair `i > 2308`, i.e. a path `3i ~> 2i` or
`i-1 ~> 2i-1` of the Collatz map, of a cycle-gate type.

**Evidence and caution.**
* The new cycles of `T` sit at upper approximations of `log_2 3`:
  `(a,K) = (3,5), (5,8), (10,16), (17,27), (29,46), (34,54), (46,73)`.
  A late fragile pair would need a later good upper approximation, such as
  the convergent `485/306`.
* The `3n-1` analogue has a late member, `i = 12029` at `84/53`, after a gap
  from 410. So late members *can* occur, and the strong form is less secure
  than the finite form.
* Heuristically, fragility is a cycle problem for the one-flip mutants. The
  same Diophantine rarity that makes nontrivial Collatz cycles implausible
  should make fragile pairs finite.

**Why it matters.** Fragile pairs are Collatz's "microcosm" in the pairing
family: the only places where the truth of tree-ness is sensitive to a
single local change. They are the pairing-ladder counterpart of the no-cycle
half (NC). HYP-9135 is a Pi^0_1-type statement, given Collatz below the
relevant heights. It can be attacked by the Simons–de Weger / Hercher
machinery on the mutant cycle equation.
