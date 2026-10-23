# A sharp five-vertex obstruction: full source plus endpoints does not determine the complement target

Status: PROVED explicit hostile, with exhaustive controls at smaller orders. This is not a hostile to the endpoint-marked *target-class* decoder: it distinguishes that decoder's input from a source-only observation.

Let T have vertices0,1,2,3,4 and arcs

1→0, 2→0, 0→3, 4→0, 2→1, 3→1, 4→1, 3→2, 4→2, 4→3.

It has five directed Hamiltonian paths. In particular

P=(4,0,3,2,1),     Q=(4,2,0,3,1)

are paths with exactly the same actual ordered endpoints(4,1).

For a directed Hamiltonian path R, define C_R(T) by retaining its four arcs and reversing every other arc. Then

H(C_P(T))=15,     H(C_Q(T))=11.

Since H is invariant under isomorphism and converse, these two targets belong to different converse-merged classes. Their sorted outdegree lists nevertheless agree, both(1,2,2,2,3).

Thus even the ENTIRE source tournament together with its two actual endpoints does not determine its target class under non-path arc flipping. Any source-only statistic, including any collection of endpoint-marked source incidences or source electrical invariants, also fails on this pair. This says nothing about additional measurements made on the target or measurements retaining the actual path's internal order.

The obstruction is order-minimal. On at most four vertices, fixing ordered endpoints leaves at most two internal vertices. If there are two, the two possible path orders would require both directions of their mutual arc. Hence there is at most one directed Hamiltonian path for each ordered endpoint pair. The five-vertex example is the first possible failure, not merely the first one a search happened to find.

## Exactly what is missing

The path-position sidecar differs. For example, the successor of4 is0 in P and2 in Q. Both candidates use0→3, so an arbitrary shared-edge observation is not necessarily discriminating. Preserving the full labelled Hamiltonian path makes C_R(T) deterministic by definition; the research question is which smaller shared-incidence bank suffices for the intended quotient.

The Session40 marked-quotient/chirality theorem reconstructs a tournament supplied as a marked object. Applied to the actual complement target, it can authenticate that target's ordinary class. It does not generate the target from the source endpoints alone. Full class-wiring reconstruction additionally requires the complete actual source/path inventory and the proved transport of observations along each such path. Literal partner authentication remains stronger still.

## Exact reproduction convention

List unordered pairs lexicographically01,02,03,04,12,13,14,23,24,34, with bit1 meaning the smaller vertex points to the larger. The source mask is4; C_P(T) is871 and C_Q(T) is733. The checker directly enumerates directed Hamiltonian paths and independently applies arc-by-arc complementing. It also checks every tournament and every endpoint pair through order4. No isomorphism oracle is needed to separate the targets, because their path counts already differ.

Run python3 check_endpoint_wiring_hostile.py and python3 -O check_endpoint_wiring_hostile.py.
