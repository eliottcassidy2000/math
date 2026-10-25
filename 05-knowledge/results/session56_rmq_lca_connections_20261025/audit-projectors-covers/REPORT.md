# Independent audit: signed projectors and hierarchy covers

2026-10-25. Status: proof audit PASS; exact checks must accompany this report as matching results.json and results.optimized.json. No arithmetic existence, novelty, external complexity classification, or LRC consequence is asserted.

## Proof findings

The signed-projector statement is sound with the stated finite nonempty meet-semilattice and Q/R coefficients. Refinement order gives upper zeta Zf(r)=sum(x>=r)f(x). Its inverse column for p is e_p=sum(x<=p)mu(x,p)delta_x, not the opposite orientation. Multiplication becomes pointwise, so all signed idempotents are precisely binary zeta-coordinate vectors. The identity is sum(e_p); identifying it with delta_top requires a top. A bottom with two incomparable atoms supplies the independent hostile identity -delta_bottom+delta_a+delta_b.

The nonnegative-real classification is correctly stronger than the integer-count result. Mass forces zero or probability mass one; strictly increasing downset rank, equality in expectation, symmetry, and independent identical draws force X=Y almost surely and hence a Dirac law. This is a statement about independent meet convolution. Coupled repetition Y=X preserves every probability law and is not a counterexample. Uniform C2 is idempotent under group convolution, not meet convolution. Signed projectors and symmetry-orbit sums cannot be interpreted as physical negative choices. The physical provenance restriction survives normalization, and a Dirac normalized state does not preserve unnormalized array counts when repeated.

The tree-cover proof correctly identifies minimum same-label unrestricted rooted-tree covers with the chromatic number of the crossing graph. Inclusion parents exist uniquely for a nonempty laminar cluster; adding full set and singleton clusters ensures actual descendant labels. Binary resolution preserves clusters. Empty requested family has minimum zero; a nonempty trivial family requires one tree. Empty literal clusters have no ordinary tree representation. Empty leaf universe is excluded, with the zero-vertex graph realization assigned a dummy label.

Private vertex labels prevent containment in the arbitrary-graph realization, while distinct edge labels realize adjacency exactly. This is an explicit polynomial-size reduction, not a cited external hardness theorem. C5 proves clique number need not equal the cover minimum. Four affine directions on F3² give precisely K3,3,3,3 and four trees. A synchronous common-refinement CRT filtration is laminar and has one tree; its nodes need not include individual prime blocks. Numerical labels or prime-partition sidecars can restore independent queries externally, so no claim of total numerical information loss is justified.

## Recovered Session55 source

The pinned Session55 CRT theorem gives A={{0,4,8},{2,3,7},{1,5,6}} and B={{2,4,6},{0,5,7},{1,3,8}}. Every A/B intersection has one cell. On the actual nine named cells the six requested blocks therefore realize K3,3, with nine edges and exactly two trees. The producer's six-label mod-2/mod-3 instance is a separate valid K2,3 example, not this nine-cell witness. Both are checked independently here.

## Verification integrity and repair

The projector producer checker uses explicit exceptions and remains substantive under python -O. Its source snapshots are pinned and were additionally matched against their original source paths at audit time; original_source_verification.json records that check.

The hierarchy producer checker uses removable assert statements. Consequently its optimized run omits core obligations, including the direct-cover comparison and hierarchy validations. Matching optimized output alone cannot certify those checks. Its mathematical proofs remain valid and its normal replay is useful; this audit does not modify its frozen files. The independent verify.py here replaces that missing optimized assurance with explicit RuntimeError gates throughout. It imports neither producer checker.

The independent algebra path constructs partitions as sets of blocks, computes meet via nonempty block intersections, and obtains projectors by backward triangular solve rather than the producer's equality masks and Möbius interval recurrence. It tests all projector products and basis diagonalization pairs for partitions on 1 through 5 labels; all 32 signed Pi3 idempotents; exact probability vectors of denominators 1 through 8; coupled repetition; C2; and a no-top semilattice.

The independent hierarchy path compares canonical exhaustive graph coloring to laminar-subfamily cover dynamic programming and reconstructs actual descendant clusters through inclusion-parent recursion. It checks all 128 families on three labels, every labelled graph on up to five vertices, affine lines, six-label CRT, nine-cell Session55 A/B, diagonal refinement, C5, empty family, singleton universe, and rejection of an empty cluster. Bounds are explicit finite controls supplementing all-size proofs.

Reproduction:

    python3 verify.py > results.json
    python3 -O verify.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 verify.py --self-test-failure
    python3 -O verify.py --self-test-failure

Both last commands must exit nonzero with RuntimeError. SHA256SUMS pins the entire delivered audit except itself. SOURCE_PINS.json pins all frozen input files. Producer checksum transcripts and normal/optimized replays are retained. This is a full overlay artifact; producer files and remotes are untouched. No blocked mathematical action remains.
