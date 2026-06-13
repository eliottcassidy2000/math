# Cauldron Two-Block Adversarial Schedules (S620)

The user's suggested rule is not just a cosmetic turn-order variant.  It
changes the residue channel each player uses to attack.

I modeled the attack-only version: on turn `n`, the scheduled player places
`n` into one of the opponent's active cauldrons; a cauldron that contains an
allowed additive witness boils and is removed.  The two-block schedule is

`A, B, B, A, A, B, B, ...`

so A owns residues `{0,1}` mod `4` and B owns `{2,3}` mod `4`.

The old parity schedule is lopsided under weak two-term sums.  A attacks with
odd numbers, and odd plus odd is even, so A's stream has no pair-sum
self-hit.  B attacks with evens, and even plus even is even, so B's stream
does self-hit.  That is why S619's parity game let B win the weak small cases.
The advantage is not density; each player still moves half the time.  It is
closure.

The one-handicap two-block word changes the closure:

- A stream `{0,1}` has pair sums `{0,1,2}` and hits both of its target lanes;
- B stream `{2,3}` has pair sums `{0,1,2}` and hits lane `2`, but misses lane
  `3`;
- strict `A,A,B,B` has the complementary asymmetry: A misses lane `1`, B hits
  both lanes.

The exact minimax results match the residue story.  In weak `2v1`, parity is
B at `14`, while one-handicap two-block is A at `5` and strict two-block is A
at `6`.  In weak `2v2`, parity is B at `14`, one-handicap two-block is A at
`13`, and strict two-block is B at `15`.  So the first-turn handicap is not
an incidental detail; it selects which period-four channel has the fuller
early closure.

This suggests a clean abstraction:

1. A schedule word induces a stream set `S_P` for each player.
2. A boil relation induces an arity-limited sumset closure.
3. The strategic bias begins where `closure(S_P)` intersects future targets in
   `S_P`.
4. Extra cauldrons are reserve capacity that delay the opponent's first useful
   self-hit.

That is a tiny analogue of the Collatz/two-block issue.  Density says the
players split turns evenly.  Residue-channel closure says one player may be
able to make legal witnesses while the other is feeding values into a blind
lane.  This also feels adjacent to the LRC lesson: the retained object is not
the scalar amount of overlap, but which orders and which residue channels of
overlap survive.

The concepts worth pursuing:

- **Additive combinatorics:** classify schedules by `S+S`, `S+S+S`, and
  finite-sum reachability in `Z/mZ`.
- **Automatic sequences:** period-four is the first square-wave case; more
  serious schedules are morphic, Sturmian, Beatty, or low-discrepancy words.
- **Maker-Breaker Schur games:** the board is the hypergraph of triples
  `(a,b,a+b)`; cauldrons are bins, and a boil is a completed edge.
- **Combinatorial game theory:** gift-or-poison play turns this into a finite
  partizan game over active hypergraph resources.
- **CRT/Collatz cauldrons:** make boiling modular, then ask whether residue
  channels can simulate the two-block `2`-adic/`3`-adic pressure from the
  Collatz thread.

The next computation should not merely increase cauldron counts.  It should
enumerate schedule languages and build the residue-channel graph for each
language, then use exact minimax only on representative channel types.

Artifacts: `04-computation/cauldron_two_block_adversarial_s620.py` and
`05-knowledge/results/cauldron_two_block_adversarial_s620.out`.
