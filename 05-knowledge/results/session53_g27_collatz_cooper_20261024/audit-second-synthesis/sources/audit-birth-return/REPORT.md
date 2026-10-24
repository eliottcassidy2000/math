# Independent audit: graph birth and Collatz return packets

Verdict: PASS for the stated elementary theorems and finite computations, with scope limitations below. No repository edits. Producer packets were reported frozen by the parent; every manifested byte was independently hashed before testing. This audit does not certify novelty or external literature priority.

## Mathematical review

The graph polynomial follows from the disjoint degree-zero, degree-one and degree-two cases at the new vertex. A singleton is one available vertex, not two choices. A degree-two attachment requires distinct old components; otherwise it closes a cycle. The three terms have component counts k+1, k, k-1 respectively. Taking the lowest nonzero coefficient proves the stated trichotomy; in particular a J_(m+1) term can preserve m when E_m=0. The P4 middle-vertex attachment is a valid witness. Endpoint-pair states retain precisely the component relation and exposed vertices needed for subsequent births. Fresh global vertex labels and valid input states are preconditions; the state-only API cannot detect reuse of a label hidden inside an old path.

For G27, new neighbors are exactly 9 and 22. Deleting 27 from each unoriented Hamiltonian path gives either one of three endpoint extensions or one of 32 compatible two-covers. Full compatible-cover set equality is checked against an independent selected-edge forest search. The latter search has degree caps and rollback union-find cycle exclusion; the former search enumerates vertex orders, fixing a leaf as start when available and canonicalizing reversal. The 1,116 total G26 two-covers are not all usable. This is a finite decomposition, not a Hamiltonicity threshold or all-N recurrence.

The Collatz parity-cylinder induction is valid because the two source lifts have endpoints of opposite parity, their difference being the odd numerator 3^s. Prefix affine formulas retain the source rather than comparing only with a reset child. All 58 pre-descent slopes exceed one, and the 59th endpoint difference is -4+(3^37-2^59)t<0. Hence first shortcut descent is exactly 59 throughout n=27+2^59 t, t>=0. There are 37 odd shortcut steps, giving ordinary clock 96; inserted odd peaks cannot be earlier descents. Independently generated prefix_certificate.json retains all 59 slopes, denominators, carries, baseline values and ordinary clocks. Direct ordinary trajectories additionally check t=0,1,2,7,100,2^200+1.

The exact full U word needs modulus 2^71, whereas the 70-letter T word only needs 2^70. Its endpoint 1+3^41 t is not a common total-stopping certificate. CRT gives exact finite-word residue counts over complete periods, not convergence or short-window/prime independence. Fixed-word mod-7 affine transport has difference-character sign (-1)^s. The explicit same-residue pair 27,34 disproves an autonomous mod-7 quotient. No bridge to square-sum graph dynamics or s7 supercongruences follows.

## Exact verification routes

Graph check.py: paths and edge_covers provide separate enumeration routes; birth and main compare the complete deletion-cover sets. state_birth/state_transfer and responses implement the general transfer. test_rejections.py compares all 1,099 labeled graphs of orders 1 through 5, totaling 33,646 spanning linear forests, and exercises the P4 near-minimum channel, component pairing, singletons and higher-degree fibres.

Collatz check.py: orbit and oddorbit directly compute the clocks; prefix constructs affine words. run verifies the k=59/70 cylinders, every pre-descent numerator inequality, all 2,046 words of lengths 1..10, 14,322 CRT lifts, 85,932 directed Paley-pair checks, 129 lifts per return word and 32 exact-U lifts. Checks use explicit exceptions and remain active under optimization.

Both packages were copied using only the producer manifest entries, plus their manifests, into separate directories. No unlisted two_covers.json, cache, repository source file or external dependency was needed. Normal and -O check outputs match each other and the producer's recorded bytes; the graph test outputs also match in both modes.

Nine deliberately false mutations fail in both normal and optimized modes: same-component joins, doubled singleton ports, deleted join transitions, omitted near-minimum middle case, wrong odd affine multiplier, omitted U final bit, universally positive Paley orientation, wrong ordinary clock, and false common total stopping. audit_results.json records rejection messages.

## Packaging and limitations

Graph MANIFEST.json SHA256: 22e395ef8e62cec05be58d3ed681d227d02607157a6d241399d3afda41ccd12f (9 payload files).
Collatz manifest.json SHA256: 2640877e5482ff3a93821b237451d146f8aa1b96b74d9be4eed9c99fbb795ab6 (12 payload files, including 8 source snapshots).

Graph SOURCE.json contains inherited source paths and hashes, not source snapshots. Thus graph runtime/proof closure is verified, but the packet does not itself carry the complete provenance documents. No claim that those references were packaged is made. Collatz source snapshots are included and hash-verified; review of the inherited family section confirms its distinction between reaching 47 and source-relative descent for the exceptional source 27.

This audit's MANIFEST.json is a closed payload list, excluding itself. Its inputs contain only the producer-manifest file sets plus their manifests. Run python3 audit.py from any directory after relocating that closed set. It regenerates disposable relocated/mutants work directories and output certificates. The independent proofs, not the finite test universes, support the all-size and all-height conclusions.
