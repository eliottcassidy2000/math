# HYP-2189: Cauldron Game As Schur Coloring With Removal

**Status:** CONFIRMED exact small scout; OPEN higher-`k` and higher-arity variants

## Claim

The cauldron game is an online additive-coloring problem whose first-boil
version is governed by Schur-type finite obstructions, while the "boil and
remove the cauldron" version is a distinct sacrifice/removal dynamic over the
same additive predicate.

For the base rule with `k=3` cauldrons and a two-term boil condition
`A+B=C`, the literal "three values" reading uses distinct summands. Thus the
first-boil variant is the weak Schur problem: safe prefix `23`, forced boil at
`24`. The classical repeated-summand comparison has safe prefix `13`, forced
boil at `14`.

The removal variant should be treated separately because a player may
sacrifice one color class to preserve cleaner active cauldrons later. For
`k=3` and the literal distinct two-term rule, S618's exact active-state search
gives optimal last boil at `27`; under repeated summands, the last boil is
`20`. Under the distinct finite-sums rule, where any previous subset of size at
least two may sum to the new value, the exact `k=3` last boil is `25`.

## Current Evidence

`04-computation/cauldron_game_s618.py` stores exact small searches in
`05-knowledge/results/cauldron_game_s618.out`.

Computed first-boil thresholds for `k=3`:

- distinct two-term cauldron: safe `23`, forced `24`;
- repeated two-term classical comparison: safe `13`, forced `14`;
- distinct two-or-three-term cauldron: safe `22`, forced `23`;
- distinct finite-sums cauldron: safe `21`, forced `22`.

Computed all-boiled/removal thresholds:

- distinct two-term cauldron: `k=1,2,3` last boils at `3,10,27`;
- repeated two-term comparison: `k=1,2,3` last boils at `2,7,20`;
- distinct finite-sums cauldron: `k=1,2,3` last boils at `3,10,25`.

The exact-three-summand, three-color removal threshold is not brute-forced
here; it is a higher generalized Schur/weak-Schur problem.

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

S618's route tournament is transitive with one Hamiltonian path. It ranks
`removal active-state DP > finite-sums subset obligations > weak Schur
first-boil > fixed-arity generalized Schur > classical Schur comparison >
greedy residue coloring > raw cauldron labels`.
