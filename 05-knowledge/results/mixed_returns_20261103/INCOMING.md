# Reconciliation with concurrent sessions 155 and 156

Status: **CITED proofs / VERIFIED literal return examples.** Read at incoming
commit `3f58f67cee` after the main packet was completed; this is an integration
note, not a claim to have discovered the incoming results.

Session155's [mixed-return theorem](../session155_spectrum_20261103/returns/MIXED.md)
uses a closed 6,400-element monoid modulo30 to exclude formal integer resets
for seven macros and a specified seed residue class. Its separate
[depth theorem](../session155_spectrum_20261103/returns/NOTE.md) characterizes
the fixed D lattice at every depth. Our positive-cone test instead treats
786 blocks and the seven supplied targets, excludes physical predecessors
even when integer gates pass, and needs no finite word-length cutoff. These
are complementary tests with different universes, not interchangeable claims.

Session156's [physical predecessor classification](../session156_spectrum_20261103/returns/NOTE.md)
independently uses the same physical cone. It proves an all-word trailing-b
criterion for its seven-macro target family and excludes53 last atoms in a
broader source test. Its other actual sources give

    v(q1) P_b = 3 v(q5 a) - 6t,
    v(b a) P_b = 3 v(q5 a) - 7t.

Here v(q5 a)=(13,5,9), v(q1)=(3,3,3), v(b a)=(5,2,3), and t=(2,-1,0).
Appending D, whose gain is -4, gives corrections24t and28t. The resulting
24-vertex target has H=12784595 and the two26-vertex outputs have
H=38353833 and38353841. Appending another b reverses the correction signs;
its actual targets and outputs retain ports26 and28. `incoming_audit.py`
replays their literal-graph checks normally and under optimization through
our reproduction driver, with all six computational inputs hash-pinned.

These positive identities do not contradict our finite-bank no-cycle theorem.
They are transported seed relations along unbounded suffixes, not a return
to a fixed finite bank. They also show why our seven-target obstruction must
not be stated for all actual sources. A graph word, port, correction gain
and vertex clock remain necessary even when two sources share a scalar H.

The incoming [two-apex state language](../session156_spectrum_20261103/states/NOTE.md)
offers a concrete alternative to the saturated return catalog: a finite
weighted transfer preserving a chosen apex order, intersected with an exact
strongness language. Its witnesses and bounded-modulus coverage do not imply
exact scalar coverage or a source selector. The next experiment can search
there for the missing H295/H883 branches while retaining the marked-port
observer. Before using it as a replacement rule, certify the literal graph,
named port, complete state, legal conversion program and both successor
branches. The incoming [clone-family bridge](../session156_spectrum_20261103/interfaces/BRIDGE.md)
also has its own mod49 exclusion; that is a different condition on a
different family, not our integral return lattice.

No incoming claim is needed for the cone, determinant or carry proofs in this
packet. The new integration supplies positive controls and a better next
search language, without weakening the existing stopping reason.

## Later incoming projective classification

Commit `ebdd5e4325` arrived during publication. Its
[eight-macro projective decoder](../spectrum31_return_20261103/PROJECTIVE.md)
proves disjoint ratio images on [1,2], reconstructs every source word on
the supplied seed ray, and classifies its exact ternary edges. It independently
rules out iteration of the95/97 companions inside that infinite family.
This is **CITED**, not an independently rerun certificate here. Our786-block
interval union need not have disjoint images and is not claimed to decode
words; our expansion/determinant argument instead excludes finite exact
seed cycles without restricting all seed bulks to one ray. Neither statement
replaces the other's hypotheses.

That packet's [joint F7 action monoid](../spectrum31_return_20261103/TRANSPORT.md)
retains the F49 field state and transverse F7 charge. It is a different
representation from our48-element GL2(F3) action on normalized bulk phase.
The cardinalities are not interchangeable carrier sizes. Its
[two-live-port branch failures](../spectrum31_return_20261103/LIVE_FIBRES.md)
also reinforce the next obligation: a new source selector must certify both
successor branches, not just scalar coverage or one successful return.
