# Independent audit of the m=12 fixed-color counterexample

**Result: PASS.** The explicit integer-difference witness has exactly five
graceful tree completions, and its at-most-two-color reconfiguration graph is
the disjoint union of a triangle and an edge.

## Independent computation

`independent_m12_audit.py` transcribes the nine fixed edges directly and uses a
disjoint-set cycle/connectivity test. It does not use the discovery search,
its component mapping, its stored support lists, or the supplied verifier's
graph traversal. It compares the transcription with `witness.json` only after
the mathematical checks pass.

There are 7 choices for difference 6, 4 for difference 9, and 2 for difference
11. All 56 triples were tested; exactly 5 are trees. Their edges agree exactly
with the five rows of `ARITHMETIC_FIXED_COLOR_COUNTEREXAMPLE.md` in the listed
order. In zero-based state indices the width-two graph has exactly the edges
(0,1), (0,2), (1,2), (3,4). Every pair within a component differs in one color;
every pair across components differs in all three free colors. Thus even the
width-one graph of this particular fixed-color fiber has the same components.

The audit also verifies the four fixed-forest components, the complete actual
edge lists for every supported contact type, and all contracted loops. Its
full machine-readable output is `independent_m12_audit_receipt.json`.
Both the existing `m12_exact_receipt.json` and a fresh run of the supplied
`verify_m12_witness.py` agree on the complete state list and component list.

## Proof review

The fixed edges are acyclic. Contracting their four components converts any
completion into a three-edge tree on those four components. A contracted loop
would already make a cycle with the fixed forest, so loops are inadmissible.
The union of all nonloop supports is precisely the path B--A--D--C, hence all
three of its edges must occur. Assigning the three free colors to distinct
path edges admits exactly the two assignments stated in the proof. The first
has multiplicity 1 x 3 x 1 = 3; the second has multiplicity 2 x 1 x 1 = 2.
All three assignments by color differ between these two cases, proving the
width-two disconnection without reliance on computation.

Therefore the claimed counterexample to connectivity of arbitrary nonempty
fixed-color fibers, and the impossibility of a width-two path between the
two groups that retains all nine shared edges, are valid.

## Scope and verifier note

This audit does not independently establish the separately referenced
upper-tail theorem or the negative search results for m=9,10,11. The sentence
asserting existence of an unrestricted global width-two path relies on that
external theorem; the local counterexample and shared-edge obstruction do not.
No minimality or global one-color disconnection conclusion follows.

The supplied verifier's output is correct, but its assertions only require
more than one state and exactly two components. The independent audit adds
explicit assertions for all five states, the exact simple-graph edge set,
component sizes 3 and 2, and all cross-component distances. Self-adjacencies
in the supplied verifier are harmless for its component computation; the
independent receipt records the loop-free reconfiguration graph.

Reproduce with:

    python independent_m12_audit.py

No supplied files were changed and no Git actions were performed.
