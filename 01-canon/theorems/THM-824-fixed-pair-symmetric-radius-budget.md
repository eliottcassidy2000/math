---
id: THM-824
title: Fixed-pair folded target and symmetric Minkowski radius budget
status: RESERVED — exact all-size fixed-(13,5) factorization proved; canonical proof and replay in progress
source: codex-2026-07-15-S10 deep sum-arc continuation
depends_on: [THM-774, THM-789, THM-821]
related: [THM-803, THM-807, THM-817, HYP-6820]
---

# THM-824 — Fixed-pair symmetric radius budget

Namespace reservation after a live-main check that THM-824 is free.

For `Q(t)=||9t||+||4t||`, the exact superlevel target
`{Q>=11/13}` is the union of the two radius-`2/169` balls centred at
`5/13` and `8/13`.  If `E` is nonempty compact and `R=-R` is compact with
`0 in R`, the current proof shows that `E+R` lies in this target exactly when
the maximum displacement of `E` from those two centres plus the circular
radius of `R` is at most `2/169`.  Applied to the fixed folded branch, this
replaces all component-by-cell containment obligations by one exact additive
radius budget.

This stub reserves only the theorem and artifact namespace.  The no-switch
lemma, necessity/sufficiency proof, scaled `13`-phase corollary, exact replay
on the THM-821 bank, owner-transport guardrail, and Tournament Analysis are
being written before promotion.  The result is fixed to `(x,y)=(13,5)` and
does not by itself prove that every core violates the budget or that the
global sporadic branch is empty.
