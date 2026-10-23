# An actual kernel that is normal only after restricting the chart domain

23 October 2026. **FINITE-EXACT**, with a direct proof of the mechanism.
This is an independent hostile control for [KERNELS.md](../KERNELS.md).

Take the lexicographic regular nine-tournament `R9=C3[C3,C3,C3]`, labeling
the outer triples `0,1,2`, `3,4,5`, and `6,7,8`. Its full internal group
has order 81. Between two copies put

    B = diag(J3,C6),

where `C6(i,j)=1` exactly when `j-i mod6` belongs to `{0,1,3}`.
Every row and column has degree three. The first three rows coincide, and
the other six rows are pairwise distinct; the analogous column statement
also holds. These are actual biregular tournament attachment blocks.

Let `H <= G x G` be the actual block stabilizer. Full enumeration gives

    |G|=81, |H|=9,
    P=projection_row(H)=K_row=C3,
    Q=projection_column(H)=K_column=C3,
    H=K_row x K_column.

Both kernels rotate only the first triple and fix all other owners. In
particular the accepted quotient has one element, with exactly nine actual
lifts. The extra six-by-six block eliminates every other simultaneous
internal rotation; this last finite assertion is checked against every
pair of the 81 full internal automorphisms.

The outer rotation `sigma(i)=i+3 mod9` is an internal tournament
automorphism but is not in `P`. It conjugates the first-triple rotation
to a second-triple rotation, which does NOT preserve the cross block
one-sidedly. Thus `K_row` is not normal in `G`, while it is normal in `P`.
The quotient `G/K_row` is not a group; the restriction to the projected
domain in KERNELS.md is essential, not cosmetic.

## Source and target kernels act on opposite sides

Construct a literal target block by transporting the row cell with
`sigma` and the column cell with `sigma^2`. Let `r=(sigma,sigma^2)`.
The source-to-target comparison relation is exactly `R=rH`, with nine
pairs. After translating by those actual projected representatives, it
returns to `H` and hence to the domains `P,Q` above.

The SOURCE kernels act on the RIGHT of comparison maps: `f -> f*k`.
Applying the source first-triple rotation on the LEFT of the reference
row map produces a pair outside `R`. Left action becomes valid only for
the conjugated TARGET kernel `sigma*K_row*sigma^(-1)` (and analogously
on the column side). The checker constructs the target adjacency block
and verifies all of these statements from its entries.

There is also a literal observation loss. Across the nine maps in `R`,
the image of source row owner zero is respectively 3,4,5, each three
times. The quotient's scalar multiplicity nine is correct for unweighted
counting, but cannot determine that actual owner image. Its exact formal
enumerator is `3*z_3+3*z_4+3*z_5`. This supplies a concrete hostile to
discarding a kernel fibre before evaluating a port-sensitive observation.

## Reproduction and exact universe

    python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/lattice/audit_non_normal_kernel.py
    python -O -B 05-knowledge/results/attachment_lattices_square_sum_20261023/lattice/audit_non_normal_kernel.py

The script checks all `9!` permutations to recover the full group, all
`81^2` chart pairs for both source and target block comparisons, exact
conjugation and multiplication, and the owner-image histogram. It has no
external dependencies or removable assertions. Literal source and target
matrices and the failing left-action witness are retained in
`non_normal_kernel_receipt.json`.
