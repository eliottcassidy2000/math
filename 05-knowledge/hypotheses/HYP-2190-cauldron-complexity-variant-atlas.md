# HYP-2190: Cauldron Last-Boil Is Schur-Frontier Hard, While Enriched Variants Import The Harder Repo Structures

**Status:** SUPPORTED by S619 finite computation and variant analysis; OPEN as a formal complexity classification.

## Claim

The base one-player last-boil cauldron problem for fixed `k` is difficult in a
Schur/Rado frontier sense, but it is not yet as structurally difficult as
A000568, LRC, or the unit-distance frontier.

The reason is negative and useful: the base state is only a canonical tuple of
live sum-free color classes. It has no `S_n` orbit/Burnside quotient like
A000568, no continuous `p_0`/all-orders cover coimage like LRC, and no
geometric embedding or unfaithfulness side-channel like unit distance. It is
still hard: as `k` grows, the safe-coloring frontier explodes and last-boil
computation inherits Schur-number style obstruction search.

The richer cauldron variants become more connected exactly when they add one of
the missing structures:

- adversarial play adds alternating quantifiers and residue strategy;
- hidden labels/orbit views add A000568-style quotient side data;
- modular/CRT rules add Collatz/LRC residue channels;
- LRC depth cauldrons add all-orders coverage and the Vitali wall;
- unit-distance cauldrons add geometric embedding side-channels.

## Evidence

`04-computation/cauldron_complexity_variants_s619.py` stores the finite scout
and `05-knowledge/results/cauldron_complexity_variants_s619.out` stores output.

One-player exact last-boil state counts:

- weak two-term: `k=1,2,3` last boils `3,10,27`, states `3,80,184948`;
- repeated/classical: `2,7,20`, states `2,20,6020`;
- finite-sums: `3,10,25`, states `3,77,62554`.

Bounded pressure test:

- weak two-term `k=4` removal hits the `250001` state cap without solving;
- weak two-term first-boil frontier for `k=4` crosses `250000` canonical states
  already at `n=13`, while `k=3` dies exactly at `n=24`.

Attack-only adversarial cauldrons:

- if A places odd `n` into B's cauldrons and B places even `n` into A's, then
  two-term weak sums give B a residue advantage: `odd+odd` cannot hit an odd
  target, while B's even stream can boil via `2+4=6` and later witnesses;
- allowing three/finite terms lets A use `1+3+5=9`, flipping the asymmetric
  `A=2,B=1` case.

## Variants Worth Pursuing

The S619 variant tournament ranks attack-only adversarial parity cauldrons
first, followed by Collatz residue cauldrons, LRC depth cauldrons,
gift-or-poison adversarial cauldrons, modular CRT cauldrons, and unit-distance
geometric cauldrons. The ranking is not a theorem; it is a proof-route atlas.

The most promising next concrete game is attack-only/gift-or-poison adversarial
cauldrons with parity-residue normalization, because it already shows the
Collatz/two-block lesson: density is blind to correlated residue streams.

## Tournament Analysis

Vertices are cauldron variants and proof routes, not cauldrons. The pairwise
observable is a weighted score combining repo connection, proof payload,
difficulty, novelty, and computability. The switch/gauge favors variants that
retain side-channel payload while still permitting exact finite computation.
The S619 route tournament is transitive, has score histogram `{0:1,...,11:1}`,
zero directed `3`-cycles, and one Hamiltonian path.

The quotient preserved by the base last-boil DP is active survival of the
selected additive boil predicate. It destroys cauldron labels and, unless
witnesses are retained, the reason a cauldron is near boiling. The challenged
assumption is that cauldron difficulty comes from the same place as tournament
or geometric difficulty; in the base game it does not.
