# HYP-2192: Block-Turn Adversarial Cauldrons Are A Shared-Removal Minimax

**Status:** SUPPORTED by exact S621 minimax for `k<=3`; OPEN for larger `k`,
other start-player conventions, and complexity classification.

## Claim

The user's adversarial cauldron variant, where one player moves twice and then
the other moves twice, is best modeled as alternating control over the shared
S618 all-boiled/removal DP.

The payoff is literal final-boil time.  The Delayer chooses placements to
maximize the last boil; the Spoiler chooses placements to minimize it.  A
placement into an active cauldron that creates the selected additive relation
removes that cauldron, and play ends when every cauldron has boiled.

This is not the same game as S619's parity cauldrons or S620's two-block
schedule-channel scout, where players attack separate cauldron pools and the
main invariant is each player's residue stream.  S621 has one shared active
resource pool and alternating minimax quantifiers.

## Evidence

`04-computation/cauldron_block_turn_minimax_s621.py` stores the exact scout and
`05-knowledge/results/cauldron_block_turn_minimax_s621.out` stores output.

Schedules:

- `single-start`: `D,S,S,D,D,S,S,...`, matching the first-turn-only-one-move
  normalization;
- `double-start`: `D,D,S,S,D,D,...`, the unnormalized comparison.

Exact `k=3` final-boil values:

| Variant | One-player removal | Single-start block minimax | Double-start block minimax |
|---------|--------------------|----------------------------|----------------------------|
| distinct two-term | `27` | `13` | `14` |
| repeated/classical two-term | `20` | `11` | `10` |
| finite-sums distinct | `25` | `13` | `14` |

Thus the block-turn rule adds genuine adversarial content.  In the literal
weak two-term rule, the first full Spoiler block cuts the `k=3` final boil from
`27` to `13`.  Giving the Delayer the first full two-move block changes the
same value to `14`.

Representative optimal single-start weak two-term `k=3` boil line:

- `S@7`: cauldron `(1,4,6,7)` boils via `1+6=7`;
- `S@11`: cauldron `(2,5,8,9,11)` boils via `2+9=11`;
- `D@13`: cauldron `(3,10,12,13)` boils via `3+10=13`.

The scout sidesteps resonance-energy scoring entirely; the only scalar value is
the exact minimax final-boil time.

## Tournament Analysis

Vertices are schedules/proof routes, not runners, arcs, or raw cauldrons.  The
pairwise observable is which route preserves the requested adversarial block
predicate.  The switch/gauge is

`2*prompt_match + 2*strategic_payload + exactness + repo_connection + computability`.

The resulting route tournament is transitive:

`single-start block minimax > double-start block minimax > attack-only parity cauldrons > gift-or-poison cauldrons > one-player removal DP > raw first-boil Schur > raw cauldron labels`.

Fingerprints: score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, zero directed
`3`-cycles, singleton SCCs, one Hamiltonian path, and `11` edge flips versus a
computability-only gauge.

## Assumption Challenge

Alternate vertex sets considered: active cauldron states, controllers,
two-move blocks, schedules, boil witnesses, parity residue streams, proof
routes, and raw cauldron labels.

The quotient preserves exact minimax final-boil time under the selected block
schedule and additive boil relation.  It destroys interchangeable cauldron
labels and all off-equilibrium non-optimal continuations; the stored trace keeps
one deterministic optimal line.

The challenged assumption is that "adversarial cauldrons" denotes one game.
It does not: S619/S620 are separate-pool residue/schedule-channel games, while
S621 is a shared-pool alternating-quantifier removal game.

## Next Checks

- Solve the opposite-start schedules where the Spoiler receives the initial
  one-move or two-move block.
- Add alpha-beta pruning and policy extraction for `k=4` pressure tests.
- Test whether block schedules with period `(a,b)` have monotone value
  envelopes between one-player removal and immediate-spoiler control.
