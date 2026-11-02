# Escaping every six-vertex target17 trap with the source retained

**FINITE-EXACT complete six-vertex result, plus a PROVED ordinal-context
extension.** This resolves the first test proposed in the incoming Boolean
realization handoff; it is not an all-order descent theorem.

The incoming `boolean_realization_71_17_20261102/HANDOFF.md` records768
isolated sublevel sources and asks whether two-edge moves suffice. They do
not:288 require three. All counts below are raw directed Hamiltonian words;
vertices retain their original labels.

## 1. Complete result

On all32768 labelled six-vertex tournaments, let e(T)=(H(T)-17)^2.
There are1440 target sources and768 strict local minima under single-edge
flips. Every strict minimum has H15. They form exactly four isomorphism
classes, represented by the following lexicographic edge masks (bit1 points
from the smaller vertex to the larger):

| mask | labelled orbit size | Aut order | scores | minimum flips to H17 |
|---:|---:|---:|---|---:|
|72|240|3|1,2,2,2,3,5|2|
|140|144|5|2,2,2,2,2,5|3|
|178|240|3|0,2,3,3,3,4|2|
|314|144|5|0,3,3,3,3,3|3|

All four are nonstrong: a strong five-vertex H15 block with one source or
sink. Both target17 classes have scores1,1,2,3,4,4. A flip moves one unit of
score between two vertices. Half the minimum l1 distance of the score
multisets therefore bounds the number of flips below by2 or3, respectively.
The exact target packets attain these bounds while preserving an explicitly
recorded Hamilton path of the source. The score calculation explains why
one of the two symmetry types needs the additional edit.

Every strict minimum also has a **three-flip** path with count sequence

    15 -> 13 -> 23 -> 17,                                (1)

preserving one original Hamiltonian path throughout. Its largest squared
error is36. This is optimal: no strict minimum lies in a target-containing
component of the e<=16 sublevel graph, whereas all do at e<=36. There is no
intermediate squared-error level between16 and36, since all H are odd.
All three edits in(1) are distinct, outside the displayed preserved path.
Thus it is a source-preserving packet with a certified temporary budget;
it is not monotone descent in the original error.

For mask72 one such packet flips bit positions0,14,1, giving masks

    72 -> 73 -> 16457 -> 16459,

and preserves the path(5,0,4,1,3,2). The other three complete certificates
are in `escape_audit.json`. Applying relabellings transports them to all768
sources without changing the source's labels in the actual execution.

## 2. Native path complements give the same distance partition

For a current Hamilton path P, define K_P to keep its five arcs and reverse
all ten other arcs. P remains Hamiltonian, so the same operation is an
involution. Build the complete **labelled** guarded graph with these edges;
it has737280 path-labelled directed incidences. Parallel incidences are kept.

No strict minimum reaches H17 by one K_P. Exactly480 require two such
complements and288 require three. A legal path is supplied at each step.
For example:

    mask72, H15 --P=(5,2,1,0,4,3)--> mask26526, H31
                --P=(4,3,1,0,5,2)--> mask48, H17.

This uses actual guarded re-marking between the two complements. It does
not assert that the first named path is preserved by the second operation;
the three ordinary-edge packet in section1 is the version retaining one
common path throughout. The two operation graphs have different transitions
and error budgets even though their2/3 distance partitions agree here.
Distances are computed before any isomorphism or converse quotient; no
shortest-path lift is inferred from a scalar H quotient.

## 3. An all-size port-preserving use of the finite repair

Let L and R be arbitrary tournaments and C any one of the768 trapped
six-vertex presentations. In the ordinal join L->C->R, put
m=H(L)H(R), with the empty tournament count1. Execute one of the certified
three-edge packets **inside C**, leaving all outside arcs fixed. Then

    H: 15m -> 13m -> 23m -> 17m,

with exactly the same vertex set and clock. Choose a Hamilton path in L and
R and concatenate them with the packet's common path in C; that full path
is preserved at every stage. The maximum squared error toward17m is36m^2.
This follows from SCC/ordinal multiplication and the actual inside packet.
Optimality is only claimed among moves confined to C, not against outside
edits. Nor must the padded initial graph remain a local minimum.

This illustrates the useful repair to a failed energy argument: retain the
actual source, name a bounded temporary budget, and certify the whole legal
word. Existence of small packets in these four classes does not prove an
analogous bound for arbitrary sources or arbitrary targets.

## 4. Independent calculations

`local_escape.cpp` counts all sources by endpoint subset DP, then uses
multi-source cube BFS and minimax Dijkstra. `escape_audit.py` independently
constructs the same H table by summing the720 Hamilton-path Boolean
subcubes, each containing1024 completions. It checks all32,768 sources,
floods the two sublevels without Dijkstra, enumerates the four complete
relabeling orbits, checks the preserved path for each packet, and constructs
the entire path-labelled complement graph. Tests use explicit exceptions
and agree under optimized Python. This is independent algorithmic checking
by the same agent, not an independent-person audit or a formal proof.
