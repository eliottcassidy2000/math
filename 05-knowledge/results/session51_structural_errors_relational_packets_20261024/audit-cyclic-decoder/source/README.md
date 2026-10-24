# Cyclic prime-cell automorphism decoder

An all-size, domain-checked decoder for tournaments whose actual stable
equitable cells are singletons or consecutive cyclic prime tournaments.
The complete permutation action is recovered from gain edges and rigid-edge
pins. This extends the inherited triangle-only and two-prime-cell mechanisms;
it does not classify arbitrary tournaments or solve native routing.

- THEOREM.md gives the proof, scope, provenance, and lost-coordinate ledger.
- decoder.py performs domain recognition and returns literal generators,
  all vertex orbit sizes, group order, and exact refinement error.
- check.py independently tests internal actions by full permutation search,
  all 2160 labelled degree-one/two five-cell crossings, all 216 triple
  matchings, fifty three-five-cell controls, converse, relabelling, and a
  regular-seven out-of-domain hostile. It compares the full generated
  permutation group with literal rotation actions, not only the group order.
- sources/ contains the inherited proof inputs named in SOURCE_PINS.json.

Run python3 -B verify_packet.py --full for normal and optimized tests in a
fresh manifest-only copy. No external package, repository, or network needed.

The main conceptual correction: a tree of interactions can freeze symmetry
through a rigid crossing, before cycle holonomy is even relevant. At p=3
such crossings cannot occur; at p=5 they already do.
