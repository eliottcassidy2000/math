# Correction to Session 133 checkpoint 5811bb4

## Packing polynomials require normalization (2026-11-01)

The published checkpoint **5811bb4** contained a hostile-control assertion in
`spectrum/ATLAS_PROOF.md` and `spectrum/build_atlas.py` that H15 had distinct packing profiles
among classes through order six. The check compared variable-length tuples
`(1,7)` and `(1,7,0)`, which are **the same polynomial**. Its apparent positive
result was therefore a representation-padding error, not a mathematical
collision. That claim and check have been removed.

The repaired control strips trailing zero coefficients before comparison and
uses H9: a strong five-tournament has polynomial `1+4x`, whereas the ordered
sum of two directed triangles has polynomial `1+2x+x^2`. Both evaluate to nine
at x=2. Normal and optimized receipts have been regenerated.

This correction is scoped to the finite test, not a denial of all-order H15
collisions: `C3 ordinal-sum S4` has polynomial `(1+x)(1+2x)` and H15,
whereas a strong five-vertex H15 tournament has `1+7x`. The former has seven
vertices and was outside the published order-at-most-six atlas.

No Hamiltonian count, strong inventory, self-converse flag, forbidden basis,
obstruction-sequence term, or structural theorem changes. The H15 symmetry
control remains correct: its strong order-five classes are self-converse,
whereas its two strong order-six classes are a converse pair.
