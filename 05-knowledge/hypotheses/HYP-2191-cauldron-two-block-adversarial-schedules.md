# HYP-2191: Two-Block Adversarial Cauldrons Are Schedule-Channel Sumset Games

**Status:** SUPPORTED by S620 finite minimax scout; OPEN for general schedule
classification.

## Claim

The adversarial cauldron variant with schedule

`A, B, B, A, A, B, B, A, A, ...`

should be analyzed as an additive channel problem.  A periodic turn word gives
each player a subset `S` of residue classes.  For the weak two-term boil rule,
the first obstruction is whether `S+S` returns to `S`: if the active player's
attack stream has no pair-sum self-hit, then density is misleading because the
stream cannot create a future target in its own channel.

Ordinary parity alternation is therefore biased under weak two-term sums.  A
attacks with odd residues `{1,3}` mod `4`, and `{1,3}+{1,3}={0,2}`, so A's
odd stream cannot boil by two-term self-hits.  B attacks with even residues
`{0,2}`, and `{0,2}+{0,2}={0,2}`, so B has a self-returning channel.

The one-handicap two-block schedule changes the channels to `{0,1}` for A and
`{2,3}` for B.  Both streams now have pair-sum access to at least one of their
own target residues, and exact small minimax flips many weak two-term outcomes
from B-favored parity to A-favored block play.

## Evidence

`04-computation/cauldron_two_block_adversarial_s620.py` stores the finite
scout and `05-knowledge/results/cauldron_two_block_adversarial_s620.out`
stores output.

Residue-channel diagnostics mod `4`:

- parity alternation: A stream `(1,3)` has pair sums `(0,2)` and no self-hit;
  B stream `(0,2)` has pair sums `(0,2)` and full self-hit;
- one-handicap two-block: A stream `(0,1)` has pair sums `(0,1,2)` and full
  self-hit; B stream `(2,3)` has pair sums `(0,1,2)`, self-hit at `2`, and a
  missing `3` lane;
- strict two-block: A stream `(1,2)` self-hits only `2`; B stream `(0,3)`
  self-hits both `0` and `3`.

Exact attack-only minimax highlights:

- weak `2v1`: parity gives B at `14` using `164` cached states; one-handicap
  two-block gives A at `5` using `9` states; strict two-block gives A at `6`
  using `10` states;
- weak `2v2`: parity gives B at `14`, one-handicap two-block gives A at
  `13`, while strict two-block gives B at `15`;
- with two-or-three-term or finite-sums rules, one-handicap two-block keeps A
  winning `2v2` at `13`; strict two-block flips `2v2` from B in weak two-term
  to A in the richer rules.

## Mathematical Concepts

The natural concepts around this variant are:

- additive combinatorics of schedule channels: `S+S`, higher sumset closures,
  and missing residue lanes;
- automatic and Walsh/Rademacher sequences: `A,BB,AA,...` is a period-four
  square wave, and generalized schedules become finite automata or symbolic
  dynamics;
- Maker-Breaker/positional games on the Schur hypergraph whose edges are
  triples `(a,b,a+b)`;
- Collatz two-block residue blindness: scalar move density is not enough when
  residue streams are correlated with the update rule;
- LRC/Helly analogies: retaining the right order of overlap matters more than
  low-order averages, but the current cauldron scout is still a finite
  hypergraph game rather than a full `p_0` Vitali-wall problem.

## Tournament Analysis

The S620 route tournament uses proof routes and schedules as vertices rather
than raw cauldrons.  The pairwise observable scores residue balance,
small-game depth, repo connection, computability, and novelty.  The score
quotient is transitive with one Hamiltonian path:

`residue-channel sumset calculus -> automatic/Walsh schedule analysis ->
Maker-Breaker Schur hypergraph game -> Collatz two-block residue analogy ->
adaptive adversarial cauldrons -> LRC depth/Helly cauldron analogy -> raw
cauldron-label minimax`.

The quotient preserves forced-boil outcome, schedule channel, and explicit
additive witness type.  It destroys most cauldron-label history and exact
near-boil provenance.  The challenged assumption is that "fair" turn density
implies fair cauldron power; residue closure shows it does not.

## Next Problems

1. Classify periodic schedules by the directed graph on residue classes
   `r -> r'` when `r'` is reachable by an allowed sum from the same stream.
2. Add gift-or-poison choices where a player may place `n` into either own or
   opponent cauldrons, and compare normal-play complexity with attack-only
   play.
3. Replace period-four schedules by Sturmian or Beatty schedules to test
   whether low-discrepancy density still hides additive-channel bias.
4. Add CRT cauldrons where boiling is modular subset-sum, making the Collatz
   residue connection literal.
