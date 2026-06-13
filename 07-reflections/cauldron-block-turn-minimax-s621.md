# Cauldron Block-Turn Minimax (S621)

The useful correction is that the user's "one player goes twice, then the other
goes twice" variant has at least two faithful formalizations.  S620 treats
separate attack streams as schedule-channel sumsets.  S621 treats the same
turn-order idea as shared-pool minimax over the all-boiled removal game.

S619/S620 are separate-pool games: A attacks B's cauldrons, B attacks A's
cauldrons, and schedule words create correlated residue streams.

S621 is shared-resource minimax.  There is one pool of active cauldrons.  The
controller of the next natural number chooses which active cauldron receives it.
If the cauldron boils, it is removed.  The Delayer maximizes the time of the
final boil; the Spoiler minimizes it.

The user's first-turn normalization is `D,S,S,D,D,S,S,...`.  That matters.  In
the literal weak two-term rule with `k=3`, one-player removal lasts to `27`, but
single-start block minimax lasts only to `13`.  The unnormalized
`D,D,S,S,...` comparison lasts to `14`.  The rule is not just a stylistic turn
order; it changes the quantifier geometry of the removal DP.

The exact `k=3` table:

| Variant | One-player | Single-start | Double-start |
|---------|------------|--------------|--------------|
| distinct two-term | `27` | `13` | `14` |
| repeated/classical | `20` | `11` | `10` |
| finite-sums distinct | `25` | `13` | `14` |

For the main single-start weak case, one optimal boil line is:

- `S@7`: `(1,4,6,7)` via `1+6=7`;
- `S@11`: `(2,5,8,9,11)` via `2+9=11`;
- `D@13`: `(3,10,12,13)` via `3+10=13`.

The result is a good small laboratory for alternating quantifiers.  The
one-player DP asks for an existential delaying policy.  The block-turn game
asks for a policy robust against locally concentrated minimizing control.  That
is exactly the difference between "can I keep all resources alive?" and "can I
keep them alive through the next adversarial block?"

No resonance-energy surrogate is needed here.  The scalar is just final-boil
time.  The retained payload is the minimax value, the schedule, and one optimal
witness trace.

Tournament Analysis deliberately uses schedules/proof routes as vertices:
single-start block minimax, double-start block minimax, attack-only parity,
gift-or-poison, one-player removal, first-boil Schur, and raw labels.  The route
tournament is transitive with one Hamiltonian path and `11` edge flips versus a
computability-only gauge.  That flip count is the warning: the easiest exact
game is not necessarily the game that preserves the user's predicate.

Next concrete work: add opposite-start schedules, alpha-beta policy extraction,
and a `k=4` pressure test with pruning.  If the exact value is capped, record
frontier profiles rather than pretending the cap is evidence for a value.
