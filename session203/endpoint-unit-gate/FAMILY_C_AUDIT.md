# Independent audit of the family-C rooted correction

**VERIFIED all-order tensor identity.** This review is independent of the author's endpoint-state skeleton implementation in `session203-edge-extractor/rooted_skeleton_candidate.py`.

## Literal data

Parent core P=(2,4,41,3,15,91,31), module mask27; intermediate child core C=(0,5,41,131,15,219,159,23), module mask155. Both use copied root q=5. Insert V=U minus its root as a module toward every other core vertex, but retain its actual input incidences to q. Outgoing bitmask conventions agree with the earlier physical core repair.

The physical operation reverses core arc0->1 and adds vertex7. Old vertex0 is then a universal sink, including toward the inserted vertices, and may be deleted without changing H. This verifies the same-order replacement interpretation.

## Different tensor basis and exact check

For r nonempty runs partitioning V, let A_(r,i,j) count ordered run assignments with i in {0,1} indicating that one designated run ends at a vertex beating the copied root, and j in {0,1} indicating that another designated run starts at a vertex beaten by the root. When both conditions occur, they apply to two distinct runs. The other runs are unconditioned. Their order is fixed by their corresponding core gaps. Relabeling bins gives identical counts for different choices of designated positions.

Enumerate every permutation of the core. For each linear gap independently, allow empty exactly when the neighboring core arc is forward (boundaries are allowed), and allow occupied exactly when the left core vertex beats the module and the module beats the right core vertex. Adjacency to q imposes the corresponding endpoint condition instead of a uniform mask test. The resulting integer coefficient vector counts literal path words by their unique core order and occupied gaps. It is valid for every rooted U, without interpolation.

The independent code `check_known.py` produces exactly:

    coefficient(child) - 2 coefficient(parent)
      = A_(0,0,0) + A_(1,0,1) - A_(1,1,0) - A_(2,1,1).

For nonempty V, H(U)=A_(1,0,1)+A_(1,1,0)+A_(2,1,1) and start_root(U)=A_(1,0,1). For V empty, A_(0,0,0)=1, H(U)=start_root(U)=1. Thus in every case

    H(child(U))-2H(parent(U)) = 2start_root(U)-H(U).

The complete integer vectors are in KNOWN_C_TENSOR.json. Every other coefficient cancels exactly. This agrees with the author's structurally different tensor representation.

## Scope

This is an actual endpoint-sensitive replacement law, not a unit-correction or eventual-spectrum theorem. A physical start-one input gives correction2-H(U). In particular the directed triangle gives correction-1, rather than1. The simpler edge/triangle contexts in PROOF.md realize the same correction and explain its compatibility with the mod-four obstruction.

Reproduce with `python3 check_known.py`; all values are exact integers, and the proof is the path-skeleton bijection above.
