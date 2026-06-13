# HYP-2194: Unit-Distance Impairment Spectroscopy Exposes Load-Bearing Carrier Channels

**Status:** SUPPORTED by S622 small-size Moser-carrier lab; OPEN as a proof
method for `n=22`.

## Claim

Unit-distance construction and proof methods should be stress-tested by
controlled impairments before being scaled to the `n=22` one-bit frontier.
The productive small-size question is not only "how dense can the beam get?"
but "which retained channels break when the method is deliberately damaged?"

The S622 impairment lab tests three such ablations on the rank-4 Moser carrier:

1. **Direction dropout:** remove one antipodal unit-direction pair and measure
   the edge loss.
2. **Observed direction-support Helly:** read exact small witnesses and count
   how many unit-direction pair types they actually use.
3. **Gain-cap throttling:** forbid frontier extensions above a local gain cap
   and measure the first small size where exactness fails.

The guiding analogy is from the recent LRC/cauldron/tournament work.  LRC
Helly scale says higher overlap orders must sometimes be retained; cauldron
two-block schedules say density is blind without the right channel closure;
tournament coimage work says scalar shadows are weak unless side channels are
reattached.  Incoming S621's cauldron block-turn minimax adds a second
impairment template: give an adversary control of the removal pool and measure
how much the one-player optimum collapses.  Incoming HYP-2193's H=21 Moon
reduction adds the obstruction-first template: keep the right local load until
the global obstruction becomes finite.  For unit distances, the corresponding
side channels are unit direction pairs, high-gain extension events, dense
deletion cores, automorphism quotients, and totally-unfaithful obstruction
labels.

## Evidence

`04-computation/unit_distance_impairment_lab_s622.py` stores the lab and
`05-knowledge/results/unit_distance_impairment_lab_s622.out` stores output.

The full Moser-carrier beam with width `260` recovers the exact small values
through `n=14`, including:

```text
n=8:  14 edges
n=10: 20 edges
n=12: 27 edges
n=14: 33 edges
```

The best `n=14` witness uses direction-pair counts

```text
(0, 1, 6, 4, 0, 6, 7, 3, 6)
```

so two direction pairs are unused in that witness.  Dropping pair `0`, pair
`1`, or pair `4` still reaches `33` edges at `n=14`; dropping any of pairs
`2,3,5,6,7,8` falls to `32`.  This is the first useful "shadow-price" ledger:
usage and necessity are not identical.  Pair `1` is used once in the chosen
witness but is replaceable, while several heavily used coordinate directions
cost one edge when removed.

Observed direction support in exact small witnesses grows gradually:

```text
n=6:  support pairs 3, mask {2,3,5}
n=8:  support pairs 4, mask {2,3,5,6}
n=10: support pairs 6, mask {2,3,5,6,7,8}
n=12: support pairs 6, mask {2,3,5,6,7,8}
n=14: support pairs 7, mask {1,2,3,5,6,7,8}
```

This is not a minimality theorem.  It is a constructive upper bound saying
which direction orders a known exact witness actually needs.  The next exact
subset search should start from those masks rather than from all `2^9`
direction subsets.

Gain caps expose the local overlap-order budget:

```text
cap 2: first misses exact at n=7, reaches only 25 edges at n=14
cap 3: first misses exact at n=9, reaches 31 edges at n=14
cap 4: recovers exact through n=14
```

Thus small exactness already needs gain-4 extension events by `n=9`.  For
`n=22`, this suggests that gain-4 and gain-5 extension events should be
isolated as proof obligations, not merely left as raw beam children.

## Novel Techniques

- **Impaired-carrier spectroscopy:** deliberately drop one direction, one gain
  level, or one side channel and measure the exact small-size failure mode.
- **Direction Helly certificates:** prove a candidate dense core needs only a
  small named subset of unit-direction pairs, then reattach omitted directions
  as repair lanes.
- **Shadow-price ledgers:** rank directions by dropout loss, not just by edge
  usage in one witness cluster.
- **Gain-order sieves:** treat gain-4 and gain-5 extension events as higher
  overlap orders, analogous to LRC Helly depth; isolate them before raw
  enumeration.
- **Obstruction-first beams:** score children by edge gain plus proximity to
  known totally-unfaithful filters, not by edge gain alone.

## Tournament Analysis

The S622 route tournament uses impairment techniques as vertices, not points.
The transitive score quotient ranks:

`observed direction-support Helly -> direction-dropout spectroscopy ->
gain-cap ladder -> obstruction-first extension solver -> automorphism-canonical
beam -> raw wider Moser beam -> triangular-only construction -> graph-only
dense enumeration`.

The quotient preserves small exactness, diagnostic value, proof transfer,
computability, novelty, and risk.  It destroys exact embedding provenance
unless the technique explicitly reattaches side-channel labels.  The challenged
assumption is that improving a unit-distance method means only widening the
beam or counting faster; S622 says controlled failure modes are more valuable
than undifferentiated search width.

## Next Problems

1. Run an exact subset search around the observed masks `{2,3,5}`,
   `{2,3,5,6}`, `{2,3,5,6,7,8}`, and `{1,2,3,5,6,7,8}`.
2. Build a gain-4/gain-5 extension solver over `21`-vertex `56/57` cores.
3. Add Moser-shell automorphism canonicalization so dropout losses are measured
   modulo carrier symmetries, not one arbitrary coordinate presentation.
4. Turn totally-unfaithful subgraphs into a pre-beam penalty/obstruction score.
5. Compare the Moser small-size direction masks with triangular/Eisenstein and
   CM/class-field carriers to see which unit directions survive projection.
