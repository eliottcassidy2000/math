---
source: collatz-procgen-20260922 (mac-mini)
status: REFLECTION (provenance, not truth). The results live in 05-knowledge/results/collatz_procgen_20260922_*.md; statuses there govern.
tags: [collatz, procedural-generation, exceptional-set, choice, relaxation-ladder, barriers, e-scc, applegate-lagarias, sheets, trunk, althofer]
---

# The missing Collatz insight lives where choice would be

**Prompt (user):** "procedurally generate new possible approaches for open
problems like collatz ... test your hypotheses and explore past work
extensively ... surprisingly unexpected fringe ideas we have brushed" (with
a pasted snippet on unit-gap clocks, signed sheets and a level-11 `1/15`
density).

## What the procedure was

Past sessions had audited pasted blueprints one claim at a time. This
session instead built a *generator* and an *instrument*:

* **The generator** (the foundry) types every approach by the controls it
  is blind to (`3n-1` sheet, `5n+1` drift, planted defects, rational
  2-adic cycles, undecidability) and by whether it needs a thin exceptional
  set. It turns "is this idea promising?" into "which control must this
  idea fail on, and does it?"
* **The instrument** (the exceptional-set profiler) counts residue classes
  with no descending certificate. It applies to Collatz, to both sheets, to
  `5n+1`, and to any relaxation with choice.

Neither is deep. What they buy is comparability: a hundred-odd
generated variants can be put on one scale.

## What the scale showed

Collatz's non-descending set has dimension `h(log_3 2)=0.95`. Three kinds
of freedom move it:

* branch choice (graph `E`) makes it three orders of magnitude thinner;
* sign choice (`3n+-1`) makes it empty;
* wild multipliers (Applegate--Lagarias) leave a single class.

Partial choice pins the effect to the entries of rising runs (`6 mod 8`;
the greedy fingerprint leads with `54 mod 64`). The equivalent
choiceful reformulation of Collatz, the undirected game, does *not*
collapse.

So the "missing insight" has a job description:

* **Where.** It must do, for positive integers and without choice, what
  choice does in the relaxations: control rising runs.
* **Which half.** It concerns the divergence half. The cycle half already
  has live mechanisms (Baker-type), and the literature lane confirms that
  no published all-orbits mechanism overcomes both the sheet and the
  drift.
* **Where the sign lives.** In the sign of `2^K-3^L`. Positive plus-sheet
  cycles contract multiplicatively, so class-level certificates never see
  them. Positive minus-sheet cycles expand, and their minima are hostile.

## How the fringe ideas connected

* **Signed sheets (the snippet):** all class-level data are sheet-blind,
  because odd `b` are 2-adically conjugate. The negative cycles are exactly
  the hostile points of the relaxed plus-sheet problem.
* **The trunk `(4^i-1)/3`** (earlier sessions' Cipolla towers): it is the
  exit set of the 3-adic hostile point `1/2`.
* **Applegate--Lagarias's lone class `-1 mod 2^j`:** the one-point end of
  the ladder. Their near-free escape `1+2^(-j)` is exactly what `E` lacks
  near `1/2` (cost at least `32/27`). That is why E-SCC keeps a
  Collatz-type core.
* **Althöfer's `3n+-1` game:** the two-player version of the sign
  choice. The one-player version is trivial; the two-player version is
  open, with a prize.
* **The `1/15` density:** its exact Collatz echo, density `2/15` among odd
  `n`, is a density statistic, blocked by DEFECT and SHEET.

## What surprised me

1. **E is already in the literature** (Le--Smith's Loosened Collatz Graph).
   The reachability questions, Q1 and Q2, were not. A cheap search before
   the method saved a false novelty claim.
2. **Q2's hard core collapsed to a mod-27 statement.** Every class but two
   descends within two moves. It is short enough to kernel-check in core
   Lean in a second.
3. **Chained hostile landings are typical** (31% of shortest descents from
   the `1/2` class). A thin exceptional set does not mean an easy
   all-integers statement.
4. **My own claim that "dimension 0" holds for `E` did not survive.**
   Local-minimum growth of about `0.1` per level fits a small positive
   dimension just as well. I withdrew it within the session.

## Candidate card (not promoted)

**"Measure the exceptional set along a relaxation ladder."**
* Trigger: an all-orbits statement whose certificates are residue-class
  based.
* Action: compute exceptional-set growth for the system, a partial-choice
  family and a full relaxation. Read the missing freedom off the drop, and
  treat the hostile points of the relaxation as its escape obligations.
* Counterindication: equivalent reformulations with choice can keep the
  full dimension (the undirected game). A thin exceptional set can still
  carry an unbounded chain of escapes (Q2 near `1/2`).
* Evidence: this session; Applegate--Lagarias's single-class obstruction;
  Caraiani's `5x+1` semigroup (per the atlas).
