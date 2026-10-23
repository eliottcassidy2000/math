# Independent audit of the parent odd-cell network theorem

23 October 2026. `seven_cell_budget` read `NETWORKS.md` and
`check_networks.py` in full, then implemented a separate literal-hostile
checker in `seven/audit_network.py`. No proof repair was needed.

## Combined marking and chart rules

For cells of sizes1,3,5,7, the coprime block argument correctly forces
informative support to separate by size. Root costs1,1,2 apply per
connected nonuniform component, including isolated non-singleton cells.
Each propagation uses a special small-cell equitable-partition lemma;
the parent correctly does not infer refinement-discreteness merely from
a trivial automorphism group. Original cell names, actual owner labels,
and complete cross blocks remain input data.

The source-target chart bound3,5,21 is per component and the records retain
shared choices on all incident edges. Matching the quotient and every
non-tree edge is essential. The theorem is not a scalar signature decoder.

## Network Burnside formula

The fixed-point sum is correct for the stated full internal group C_p
and connected support. A group vector with a zero and a nonzero coordinate
has an informative edge crossing its zero/nonzero cut, so fixes no
configuration. Every all-nonzero vector gives p positional orbits per
edge, each meeting every row/column exactly once; choosing d_e orbits
therefore gives binom(p,d_e) fixed regular blocks.

The stabilizer restriction to one root is injective, so its order is1 orp.
Double-counting nonidentity stabilizers, followed by division by the
correct orbit size p^(v-1), gives the displayed symmetric count. Subtracting
it from the total yields the stated rigid-count numerator. Thus the
shared-rotation rule is not the product of independent edge orbit counts.

## Nine-vertex module construction and its scalable lower bound

The parent common-outneighbour module detector is correct: inside one
inner C3 pair the count is3, across two different triples it is1.
Consequently every internal automorphism preserves the intrinsic triple
partition, permutes its three parts by an outer C3, and rotates each inner
triple. This proves the group order81 without an assumed decomposition
algorithm.

For m named R9 cells and connected support, each equality block forces
agreement of the outer rotations but does not constrain any inner
rotation. Hence the exact group order is3^(3m+1). Every one of the3m
inner triples supports a nontrivial automorphism fixing all other
vertices. Any set of fewer than3m owners misses an entire triple and is
not a base. Choosing one owner in each triple fixes all inner rotations
and the common outer rotation, proving the exact3m base bound.

This is a lower bound for all mark choices, not just the six test marks
used by the two-cell numerical fixture. It also explains why connected
coarse support cannot supply a universal two-owner budget.

## Prime-eleven hostile: independent construction and checks

The four-vertex tournaments have score pattern2,2,1,1. Each high-side
score-one vertex gets its one H-to-L matching arc, and each low-side
score-two vertex receives exactly one such arc. Together with H->M->L
this gives all eleven vertices outdegree5.

The separate checker constructs all adjacency entries directly from that
description and independently enumerates automorphisms. The actual
regular-eleven tournament has group C3, not merely a known C3 subgroup.
Its two nonidentity maps fix exactly the eight vertices outside M.
The degree-three J3 direct-sum circulant block has exactly the same C3
as its one-sided internal kernel. It agrees entrywise with the parent's
frozen fixture.

Therefore prime size alone demonstrably fails. The C_p network theorem's
full-group hypothesis cannot be replaced by 'regular tournament of prime
order'. No contradiction with the seven-cell lemma results: its elementary
cycle-shape exclusion stops before an internal C3 with eight fixed
vertices becomes possible.

## Reproduction

    python -B 05-knowledge/results/odd_cell_completion_20261023/seven/audit_network.py
    python -O -B 05-knowledge/results/odd_cell_completion_20261023/seven/audit_network.py

Both checks pass. Output: R9 internal81, pair2187, target fibre27 and
minimum pair base6; R11 internal3 and one-sided kernel3. The proof audit
does not re-run the parent's full120^3 R5 matching-triangle enumeration,
nor assert an unproved universal observer or enumeration improvement.
