# HYP-2189: Cauldron Game As Schur Coloring With Removal

**Status:** OPEN

## Claim

The cauldron game is an online additive-coloring problem whose first-boil
version is governed by Schur-type finite obstructions, while the "boil and
remove the cauldron" version is a distinct sacrifice/removal dynamic over the
same additive predicate.

For the base rule with `k=3` cauldrons and a two-term boil condition
`A+B=C`, the first-boil variant should have safe prefix length `13` and force
a boil at `14`, matching the classical Schur number `S(3)=13`. The removal
variant should be treated separately because a player may sacrifice one color
class to preserve cleaner active cauldrons later.

## Current Evidence

Reserved during codex-2026-06-03-S618. The computation is still pending:
`04-computation/cauldron_game_s618.py` will compare first-boil, all-boiled
removal, fixed term-count, and finite-sums variants, then store output in
`05-knowledge/results/cauldron_game_s618.out`.

## Assumptions To Challenge

- Values are placed once, so the default boil relation uses distinct summand
  values. The repeated-summand variant is stronger and should be reported
  separately.
- "Unlimited terms" is interpreted as some finite subset of already placed
  values summing to another value, not as an infinite series.
- A removed cauldron is unavailable for later placements.

## Tournament Analysis Plan

Tournament vertices should be proof routes and quotient states, not just
cauldrons or natural numbers. Candidate vertices include first-boil Schur
coloring, removal/sacrifice strategy, fixed arity Rado route, finite-sums
route, greedy residues, active-state DP, subset-sum obligations, and proof
obligation quotients.

The preserved predicate is survival of the relevant boil condition over a
prefix of the natural numbers. The quotient destroys label history of
interchangeable cauldrons and, in finite-sums variants, may destroy which
specific subset witnessed the first boil unless witnesses are retained.

