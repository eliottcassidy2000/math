# Session39 final addendum audit

Status: PROVED scoped arguments; FINITE-EXACT producer replays.
No earlier frozen audit or producer artifact was changed.

## Native antiisometry obstruction

The proof in CONVERSE_IS_NOT_TIME_REVERSAL.md passes. Using the previously
certified orthogonal decomposition <1/256> plus a group of exponent at
most4, every element of order256 has form ag+h with a odd. Its self-pairing
numerator modulo256 is a^2 plus a multiple of64, hence1 modulo8.
An antiisometry would send the distinguished generator to an order256
element with numerator-1=7 modulo8, a contradiction. An antiisometry of
the full pairing restricts to the characteristic2-primary subgroup, so
unknown odd factors cannot repair this obstruction.

The operation typing is also correct. Converse fixes every converse-merged
vertex and consequently its vertex-divisor critical-group element. On
literal complement edges, reflection can exchange black twins. Because
the closing free cell is reflection-fixed but complemented, reflection
cannot interchange the endpoints of a complement pair; fixed pairs are
exactly blue. None of this changes the induced identity on vertex divisors.

Thus this converse cannot implement sign reversal of the declared K=A
linking construction, even after changing discriminant labels. The note
properly confines 'time reversal' to the formal operation K->-K and makes
no unsupported physical realization or central-charge claim. Its argument
requires the stated exponent gap and does not cover arbitrary orders or
other constructions. The positive order-two self-negative control guards
against overgeneralizing the obstruction.

Both normal and optimized check_converse_obstruction.py runs PASS136816
checks with identical outputs, saved here. They replay the frozen receipt
without rewriting it. The inherited full native group decomposition is
not recomputed by this small addendum; the first audit independently
replayed its terminal linking generator from the original edge ledger.

## Static rooted-cutoff lower bound

ROOTED_CUTOFF_BOUNDARY.md passes. For an occupied proper subset of the
cycle's cells, the support is a forest. Gauge equations can be solved with
rotation zero in the component containing the chosen root, so the actual
root stays fixed. For all g occupied cells at at mostg vertices, the subset
is a transversal. The map y->-y+2x_root fixes that root, interchanges the
two nonzero phases, and preserves every induced arc: no internal-cell arc
occurs. Occupancy-wise bijections therefore give pointwise equality of
all singly-rooted profiles throughg vertices. No global isomorphism is
asserted or follows from these occupancy-dependent local bijections.

Global C3 rotations make each profile cellwise constant. The cell partition
is equitable with phase-independent directed counts. Hence every ordinary
refinement record remains identical between the two graphs by induction.
Rooted two-vertex counts already include the distinct cell outdegrees,
so the full initial profile for cutoff at least2 actually gives the
intrinsic triple partition and is stable immediately.

Nonisomorphism, strong connectivity, and self-converseness are established
in the earlier approved family. Thus a separator on ALL strong self-converse
tournaments using only this static profile plus ordinary refinement must
have k(3g)>=g+1. This is not a prime-cyclic improvement, a sufficiency result,
a complexity lower bound, or an obstruction to actual individualization.
An actual global mark is a different observer and can repair this family.

The standalone rooted producer was read and replayed. It compares complete
canonical rooted profiles, not just their hashes, on every actual vertex
for g3,4,5, and explicitly checks1620 root-fixed full transversals.
All3612 checks PASS and the printed receipt matches its frozen receipt.
The finite test corroborates the all-g bijective proof rather than replacing
it. The exact replay output is saved as rooted_replay.json.

## Reproduction and pinning

PRODUCER_PINS.json gives the approved exact bytes. Replay the producers:

- python3 check_converse_obstruction.py
- python3 -O check_converse_obstruction.py
- python3 check_rooted_cutoff.py

Neither default replay mutates its receipt. The source scopes and dependency
hashes are retained in the producer manifests. No Git mutation was performed
by this auditor.
