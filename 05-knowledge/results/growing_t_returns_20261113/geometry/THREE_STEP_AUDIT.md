# Independent audit of the physical881-to6551 program

**FINITE-EXACT, accepted:** the actual finite program, all name maps,
patches, births, root states, and final neutral sink deletion. This is a
three-return witness, not an all-k return theorem or a closed family.

The independent checker imports no counting code from the repair lane.
It uses Python unbounded-integer endpoint subset dynamic programming,
reverses the adjacency to count starts, and directly sums ending-root
two-path covers over all subsets. It separately verifies tournament
validity and strong connectivity for every literal stage.

| stage | order | root | state(H,s,q,F) | strong? |
|---|---:|---:|---|---|
| helper from A2 at old3 |11|10|(1589,698,881,1772)|yes|
| first repaired parent |11|0|(351,1,3,1341)|yes|
| first ready target |12|11|(1691,350,351,1692)|yes|
| second repaired parent |12|0|(853,1,5,2459)|yes|
| second ready target |13|12|(3311,852,853,3312)|yes|
| third repaired parent |13|0|(495,1,33,6057)|no|
| third ready target |14|13|(6551,494,495,6552)|no|
| delete old sink8 |13|12|(6551,494,495,6552)|yes|

The first helper is exactly the canonical birth at old vertex3. The
three parent patches contain respectively5,6,12 distinct old arcs and
reproduce the supplied parent rows exactly. Their subsequent births
retain every old arc and add a vertex beaten only by the selected port.
For each birth H_new=H_parent+F_parent-s_parent, and readiness is exactly
s_parent=1. Thus the full Hamilton count and readiness are both checked.

The canonical-to-actual name map

    (0,2,3,10,9,4,5,6,1,7,8,11)

consists of the unique Hamilton path beginning at0 in the first parent,
followed by its newborn. It is a permutation, is checked against every
arc of the named alternate1691 input, and transports both the second
patch and its selected port. Appending the next newborn name12 similarly
transports the third patch and port. The checker compares both repaired
graphs in canonical coordinates against their original raw records;
it does not merely compare their scalar counts.

At the last stage old vertex8 has outgoing row0. It is a global sink,
the final active root13 is retained, and deleting8 gives exactly the
listed final rows after reindexing. Direct counts verify all four root
coordinates are unchanged, and the remaining13-vertex graph is strong.
The intermediate third parent and14-vertex target are weak; no strongness
claim is made for those two stages.

This chain evades the proved all-birth obstruction by changing old arcs
and lowering the parent count before each ready newborn. It uses the
351->853->495 parent chain. The separate183->205 module-family branch
is not a dependency of this witness.

`THREE_STEP_AUDIT.json` contains the recounted stages and SHA256 pins for
the chain program, both original receipts, canonical input and all three
raw patch records. Normal and optimized audit receipts match. Reproduce
`python -B .../geometry/three_step_audit.py`, then with `-O`.
