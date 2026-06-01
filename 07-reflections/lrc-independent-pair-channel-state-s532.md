# Independent-pair channel state for the n=14 multi-channel generalization (S532)

The user's refinement is the right next coordinate.  "Support width" says how
many harmonic/chord channels are active.  But it does not say how much
independent binary state those channels can carry once the surrounding scaffold
is fixed.

The four-vertex toy model is exact.  Keep the scaffold

```text
0 -> 1 -> 2 -> 3, plus 0 -> 3
```

fixed, and let only the two independent crossing pairs `(0,2)` and `(1,3)`
vary.  Those two bits hit all four unlabeled tournament isomorphism classes:

```text
state  scores        H  c3  SCCs
00     (1,1,2,2)    5   2  (4)
01     (0,2,2,2)    3   1  (3,1)
10     (1,1,1,3)    3   1  (3,1)
11     (0,1,2,3)    1   0  (1,1,1,1)
```

So in this scaffold, two independent pair-arcs are a complete coordinate for
the K4 ordinary isomorphism class.  Across all K4 choices of two disjoint free
arcs and all fixed orientations of the other four arcs, the audit finds:

```text
1 class reached: 12 scaffolds
2 classes reached: 12 scaffolds
4 classes reached: 24 scaffolds
```

That is the useful lesson: the metric is not merely "number of free arcs."  It
is the amount and state of an independent-pair matching inside a fixed scaffold.

For n=14 the connection is immediate.  The clasp-deleted regular 14-gon has
seven chord channels.  Channels `d=1..6` each have 12 edges, maximum matching
size 6, and 7 maximum matchings.  The diameter channel `d=7` has 6 edges,
maximum matching size 6, and exactly 1 maximum matching; all of its edges are
mutually independent.  It is literally six independent pairs plus one singleton
runner.

This sharpens the multi-channel parity lead:

```text
support width:       which channels can carry inside debt
independent rank:    how many disjoint pair states a channel can carry
matching shape:      whether the independent states are unique or scaffolded
state vector:        the orientation/tie/export state on those pairs
compatibility:       which state vectors can coexist across channels
```

The diameter channel is therefore not just "the last channel" in the
`12,12,12,12,12,12,6` inventory.  It is the unique pure independent-pair channel.
If n=14 is hard because hidden diameter ties compactify the wall, then the proof
should track the six diameter pair states explicitly.

Possible proof target:

```text
Fix the non-diameter scaffold around a candidate wall-only n=14 state.
If the six diameter independent-pair states are not the AP-balanced state,
then either the outside clasp opens or endpoint/resonance debt exports to a
deeper labelled layer.
```

This also explains why the n=4 parity law is a genuine local model rather than
just a small-n curiosity.  K4 has two diameter independent pairs; n=14 has six.
The multi-channel generalization should be a matching-state law, not only a
residue-support law.

Near-term computational tests:

1. Add a `diameter_pair_state` vector to n=14 hard-row audits.
2. Compare diameter-pair state with source measure, endpoint debt, and
   resonance-debt ratio.
3. For finite-check survivors, quotient first by non-diameter scaffold, then
   classify the six-bit diameter state.
4. Search for K4-like local windows inside the n=14 permutohedral handoff graph:
   fixed scaffold plus two independent pair toggles determining the marked
   local class.
5. Test whether AP/scalar-ramp wall classes are exactly the scaffold states in
   which every maximum independent-pair matching is balanced.

Artifacts:

- `04-computation/lrc_independent_pair_channels_s532.py`
- `05-knowledge/results/lrc_independent_pair_channels_s532.out`
- HYP-2015
