# Independent referee: forbidden-wedge topology consumer

**Verdict: PASS, with no repair.** The complete eleven-word exclusion in
`continuing10_20260907_lrc_wedge_topology.md` follows from the two audited
native-wedge theorems and the inherited complete minimum-tree/profile
certificate. The conclusion is restricted to actual connected seven-label
complements of a selected six-label body of gcd 7200 in a primitive,
distinct thirteen-speed row. The four surviving words are necessary-condition
survivors only; no unsafe realization or full 7200 closure is claimed.

## Exact role cover and direction of optimization

Under hypothetical failure, no actual margin-24 vertex can have both an
18-neighbor and a neighbor of margin 4 or 16: those two actual edges would
be one of the newly closed native wedges. Therefore each positional
24-vertex satisfies the exact disjunction

    all 18-edges absent OR all 4/16-edges absent.

If both groups are absent, both roles are allowed. Equal-margin labels
remain distinct. Taking the product of these local role choices covers
every actual wedge-free graph; no consistency of independently chosen
ratio labels is assumed. The referee checks the equivalence for every
possible local incident-edge subset at every relevant vertex of every word.

For a role assignment, the producer's permitted graph contains all locally
possible atlas edges not forbidden by that role. Its edge weight is the
inherited minimum credit over every compatible actual ratio. Any actual
graph in that branch is a subgraph, and each actual edge has at least this
weight. If the actual graph is connected, choose any of its spanning trees.
Its actual overlap credit is at least the minimum weighted spanning-tree
cost of the generous permitted graph. Thus a branch cost greater than the
full-profile excess gives a safe lift. A disconnected permitted graph cannot
contain an actual connected graph.

This verifies the required MINIMUM direction; using a maximum tree of the
generous graph would be invalid. The earlier maximum-forest argument applied
instead to credits on known actual labels and ratios. The two optimizations
have different quantifiers, and the present proof keeps them separate.

## Independent complete verification

The standalone referee imports no producer and invokes no spanning-tree
optimization routine. It enumerates all 7^5=16,807 labelled seven-vertex
trees through the elementary Prüfer bijection, checks that every decoded
tree is connected with six distinct edges, and checks uniqueness of the
complete tree universe. It then filters and evaluates every tree on each
of the 29 permitted branch graphs and on each of the 15 unrestricted
inherited graphs. These 44 exhaustive minima independently recover all
displayed values and saved attainers.

The complete branch minima, in the producer's word order, are

    0:216,114   1:114,266   2:114,168   3:114,302   4:114,190
    5:162,180   6:114,92    7:114,148   8:114,92    9:114,226
    10:114,112  11:162,66   12:102      13:162,208  14:162,66.

Consequently the deleted indices are exactly
1,2,3,4,5,6,7,8,9,10,13. The remaining indices and best branch cost/excess
are 0:114/116, 11:66/79, 12:102/103, and 14:66/87. The word without a
24-vertex correctly has one branch, giving 29 branches in total. Every
word, role, removed edge, permitted edge weight, threshold, and attainer is
checked against the full saved certificate. Strict comparison cost>E is
retained; equality is not silently promoted to closure.

This is a deletion of eleven whole projected word classes in the connected
actual-graph scope. It is stronger than checking eleven selected ratio
assignments. It has no implication for a disconnected actual complement
that does not itself contain a closed wedge.

## Reproduction and pins

    python 04-computation/continuing10_20260907_lrc_wedge_topology_audit.py
    python -O 04-computation/continuing10_20260907_lrc_wedge_topology_audit.py

Both runs pass **51,576 always-active exact gates**. Captured stdout is raw
LF and byte-identical, without text normalization. The referee supports the
outside directory and filed repository layout and pins its inherited inputs.

    producer source 6fc998e7596bfa2c4182bfe41423f3d4cbdd9d9a4f7ce77724afebdfdb653e6b
    producer JSON   4008a6ae441672637becabb0f0ba17f8e0fe9e5975657fe9dddefe2bb5766e79
    audit source    ce2e0976fcb4a08138377b19e3216cf7acc9708e4211c06090c6822756a9ada2
    audit output    9e67c262a701453f7e66bbe08250bcc2567f07b540fa3546432140ef0b2696c3

No producer files were modified. Parent owns promotion and filing. A separate
stronger third-wedge exploration is not a dependency of this frozen audit.
