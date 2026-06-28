# n=4 Tournament Subbasis Reflection

Source: codex-2026-06-28-S276  
Hypothesis: HYP-3143  
Script: `04-computation/tournament_n4_dual_basis_erdos870_codex_s276.py`

The useful surprise in the prompt is that the two n=4 models differ by
representation discipline.

The Hamiltonian-path tiling model keeps the geometric feeling of the three
remaining chords.  It gives exactly the requested class-valued table:

```text
* | E | a | b | c
E | T | + | - | S
a | + | T | S | S
b | - | S | T | S
c | S | S | S | T
```

But the full 3-bit cube has class fiber sizes `{T:1, +:1, -:1, S:5}`.  This
means `S` is not a true order-2 obstruction: it can be obtained as the
one-flip state `c`, as two-flip states, or as the three-flip state.  That is a
controlled-forgetting failure in miniature.

The alternate model fixes four arcs with partial score `0,1,1,2`.  The exact
search found `12` witnesses, all with free edges a perfect matching.  The
clean witness fixes `01,03,12,23` and leaves `x=02`, `y=13`, producing
`E->T`, `x->+`, `y->-`, `xy->S`.  This is the local packet I would hand to
the next LRC proof agent: a deterministic boundary cycle plus two diagonal
obstruction bits.

The Erdős-870 connection gives the right language.  The solved problem is not
just "everything is representable"; it is "representable at the declared order
and not below it."  The LRC analogue is:

```text
fixed filler + free obstruction basis + lower-order exclusion
```

If a quotient breaks this exact-order condition, it can still be useful as a
visualization or search coordinate, but it must carry a collision sidecar.
This is the same warning now appearing in HYP-3138 odd leakage, HYP-3139
center/boundary leakage, HYP-3141 tail/tip edge payloads, and HYP-3142's
terminal k=8 U4 sidecar.

Assumption challenged: tournament vertices need not be runners, and even
"remaining arcs after a Hamiltonian path" need not be the correct proof
coordinates.  In this packet the live vertices are representation words,
fixed filler arcs, matching free bits, class fibers, and proof obligations.

Bucket for future pulls:

- Search n=5 and n=6 for partial arc assignments whose free bits represent
  A000568 classes with exact-order separation.
- Add `first_packet_order` and `lower_order_leakage` fields to edge-witness
  ledgers.
- Treat q=3 unital four-point blocks as `C4 + matching` packets after
  AP/Goddyn-Wong labels are attached.
- Recheck the k=8 bounded-core exit: the U4/resolvent sidecar should be a
  terminal packet, not a state reachable by a lower-order scalar quotient.
