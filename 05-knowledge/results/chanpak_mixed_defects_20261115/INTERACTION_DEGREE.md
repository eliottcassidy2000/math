# The degree of an exact repair is controlled by its contacted vertices

**Interface theorem.** Fix a tournament and independent binary controls
which reverse disjoint sets of arcs. Suppose every controlled arc has
at least one endpoint in a set R of r vertices. In the Boolean control
coordinates the literal Hamiltonian path count is an integer multilinear
polynomial of degree at most 2r. Its nonconstant coefficients are even.

**Proof.** For each permutation of the vertices, its contribution is the
indicator that all successive arcs exist in the requested orientation.
Arcs outside the control support either reject the permutation or impose
no condition. Each controlled arc used imposes a frozen bit value on its
control; if a permutation requires both settings of the same control its
indicator is zero. Otherwise the indicator is a product of distinct
factors z_i or (1−z_i). At most 2r successive arcs touch the r vertices,
since each such vertex has at most two incident arcs in the permutation.
Sum the products. The coefficient of a nonempty subset S is the exact
alternating sum of counts over all settings in S, with other settings
zero. The parity theorem makes this alternating sum even. No actual
Hamilton paths are identified across different orientations here.
The degree bound is sharp: an ordinal sum of r disjoint three-vertex
one-apex blocks has H equal to the product of their counts. On each
block, toggle the two apex marks independently from 00; the highest
coefficient is −2. The product has degree exactly 2r.

An exact promised repair on a downward-closed family of legal sets
(for example matchings) holds precisely when the higher coefficients
vanish on every legal face of sizes 2 through 2r. Larger faces have
degree too high; for a legal set, all present monomials other than
the constant and singles must vanish. Conversely subtract the affine
part and invert the alternating-sum formula on the legal faces.
Thus one apex only requires pair interactions, whereas multiple
apices may hide a higher interaction from every two-control observation
at the baseline. This is a concrete specialization of the
[older Boolean Hasse interaction language](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends.md).
It suggests how to organize a multi-apex search without treating
all interaction tests as independent of the geometry.

**Caveats.** Independent refers to which graphs exist for each setting,
not to the resulting counts. Degree bounds do not provide coefficients
of the desired magnitude. Low-degree interaction and exact additivity
at the baseline do not by themselves give a step at every prescribed N.

## Literal witness: three-way defect hidden from every base square

Code a 7-vertex tournament by lexicographic vertex pairs i<j,
with bit 1 for i→j. Start with code 863667 and use three disjoint
xor masks 24, 2, 768, touching only vertices 0 or 1.
In order by subsets ∅, A, B, C, AB, AC, BC, ABC, their literal path
counts are

    79, 133, 77, 15, 131, 69, 13, 65.

All three pairwise dividends at the base vanish. The one-way
steps are +54, −2, −64, which predict 67 if added. The actual
triple value is 65, so the third dividend is exactly −2. This
makes an elementary counterexample to combining multiple
apices by checking pairwise mixed changes at a single base.
The companion counts all eight tournaments by direct
enumeration of their vertex permutations.
