# Unit-Distance Impairment Spectroscopy (S622)

The useful move this session was to stop asking only for a stronger unit
distance construction and instead ask how our current construction method
fails when we deliberately damage it.

The inspiration came from several nearby threads:

- LRC Helly scale: a finite low-order summary can miss the real obstruction,
  so retain exactly the overlap orders that matter.
- Cauldron two-block schedules: density is blind if it forgets residue channel
  closure.
- Incoming cauldron block-turn minimax: adding an adversary to the shared
  removal pool is itself a controlled impairment of the one-player cauldron
  optimum; the useful datum is how far the optimum collapses.
- Tournament/coimage work: raw scalar shadows are weak; side-channel retention
  is the proof object.
- Incoming HYP-2193 H=21 reduction: a hard global obstruction becomes finite
  after retaining the right low-order load, here Moon's 3-cycle pressure.
- S617 frontier-gain: once a carrier is chosen, count child states locally by
  gain rather than recomputing everything.

For unit distances, the Moser carrier gives nine antipodal unit-direction
pairs.  S622 asks what happens when we remove directions, cap local frontier
gain, or inspect how many directions small exact witnesses actually use.

The full width-260 Moser beam still recovers the exact values through `n=14`.
The best `n=14` witness uses seven of the nine direction pairs.  Its direction
usage is

```text
(0, 1, 6, 4, 0, 6, 7, 3, 6)
```

so pairs `0` and `4` are unused in that particular witness.  But the dropout
audit is subtler: removing pair `0`, `1`, or `4` still reaches `33` edges at
`n=14`, while removing any of `2,3,5,6,7,8` drops to `32`.  Pair `1` is used
once but replaceable.  Usage is not necessity.

The direction-support ladder is also encouraging:

```text
n=6  uses 3 direction pairs: {2,3,5}
n=8  uses 4 direction pairs: {2,3,5,6}
n=10 uses 6 direction pairs: {2,3,5,6,7,8}
n=12 uses 6 direction pairs: {2,3,5,6,7,8}
n=14 uses 7 direction pairs: {1,2,3,5,6,7,8}
```

This is a small-size analogue of a Helly number.  It does not prove minimality,
but it gives a constructive upper bound on how many direction orders a known
exact witness needs.  The next exact search should start from those masks and
ask which directions are genuinely forced.

The gain cap audit gives the other half.  If children with gain above `2` are
forbidden, the beam first misses exactness at `n=7` and reaches only `25` at
`n=14`.  Cap `3` first misses at `n=9` and reaches `31`.  Cap `4` recovers
exactness through `n=14`.  So high-gain extension events are not optional
noise.  They are the construction-side analogue of higher overlap order in
LRC.

The novel techniques I would carry forward:

1. **Impaired-carrier spectroscopy.**  Before widening a beam, remove one
   direction, one gain level, or one side channel and record exactly what
   breaks.
2. **Direction Helly certificates.**  Certify dense candidates by a small named
   subset of unit directions, then reattach omitted directions as repair lanes.
3. **Shadow-price ledgers.**  Score directions by dropout loss and repair cost,
   not just by frequency in a witness.
4. **Gain-order sieves.**  Treat gain-4/gain-5 extension events as proof
   obligations.  Enumerate them directly over dense `21`-cores.
5. **Obstruction-first beams.**  Fold totally-unfaithful filters into the beam
   score before spending time on edge-rich but doomed children.

The guiding taste is now: make the method fail in a controlled way at small
sizes, then turn the failure into a certificate vocabulary.  For `n=22`, that
means a good future run is not just a larger Moser beam.  It is a beam whose
states carry direction support, gain-order obligations, dense-core deletion
context, automorphism class, and unfaithful-subgraph risk.

Artifacts: `04-computation/unit_distance_impairment_lab_s622.py` and
`05-knowledge/results/unit_distance_impairment_lab_s622.out`.
