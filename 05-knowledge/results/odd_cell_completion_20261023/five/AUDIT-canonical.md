# Independent audit of source-only canonical charts

23 October 2026. **PASS**, within the supplied named equitable
size1/3/5/7-cell carrier, retaining all actual blocks. Read
[CANONICAL.md](../seven/CANONICAL.md) and inspected canonical.py.

## Proof review

The deterministic tree is chosen from the support graph on named cells,
not from raw vertex labels. Internal minimal templates are invariant
under vertex relabelling. For a fixed parent chart, equality of two child
block words would make their quotient a forbidden one-sided kernel
automorphism. Thus the greedy minimizer is unique, and the construction
commutes with each cell-preserving isomorphism.

It is correct that this need not give the absolute least matrix over all
chart tuples: equivariance plus a complete final matrix is enough. The
final code retains all non-tree blocks, so it cannot lose a global cycle
constraint. Equality of two full codes explicitly constructs the complete
cell-preserving isomorphism. Uniform inter-component blocks impose no
additional chart coupling but their actual quotient bits are still output.

The minimum-fibre count equals the full cell-preserving automorphism
order. The group sends any minimizing chart to another, root restriction
is faithful, and two minimizing complete chart systems differ by an
actual automorphism. A lexicographic choice by input owner labels is
correctly presented only as a gauge; it is not an intrinsic pointwise
labelling when a residual automorphism moves vertices.

For a full actual Hamiltonian path, every vertex appears once. Therefore
two distinct canonicalizer charts yield distinct component-restricted
coordinate words. Their minima are unique. Independent minimization is
valid because positions belonging to each named component are fixed;
the first difference of full path words belongs to one component and is
its own first restricted-word difference. This does not discard shared
ports or construct a hypothetical path from a quotient.

The group-gcd saving is sound: faithful restriction to every cell makes
the global group order divide each internal automorphism order. A
component containing group3 and group7 cell types is rigid. The owner's
base cost and ordinary refinement cost remain carefully distinguished.

## Separate exact controls

The audit does not repeat the already complete order31 run. It targets a
different universe: **all4320 regular nonuniform five-by-five blocks**.
An independent C5 x C5 rotation-double-orbit calculation supplies ground
truth. The canonicalizer produces exactly the same192 equivalence classes,
without merging different local orbits or splitting any one orbit. Its
minimum fibre equals the independently calculated stabilizer on every
block. All192 representatives receive a fresh global vertex relabelling
and a literal Hamiltonian path/complement transport check.

The parent's new fixed-local-orbit loss formula receives a separate
hostile test. There are four rigid degree-one matching block orbits. Each
of their16 ordered pairs is put on a three-cell tree. Varying the phase
at the shared middle cell gives exactly five different full canonical
codes for each local signature, hence80 in all. Every whole network is
rigid, and the third pair of cells is joined uniformly. This confirms
that cycle data alone do not explain all shared-incidence loss: a tree
already loses relative phase when two rigid attachments are independently
canonicalized. It also independently verifies the claimed16 x5 case of
the parent's p^(beta+r-1+a) formula.

These controls use actual supplied cells and binary blocks, not an
anonymous margin relaxation. They are not an all-tournament census or a
proof of faster unrestricted A000568 enumeration. The theorem's all-size
claim is bounded-cell, with input ownership and equitability supplied.

## Reproduction

    python -B 05-knowledge/results/odd_cell_completion_20261023/five/audit_canonical.py
    python -O -B 05-knowledge/results/odd_cell_completion_20261023/five/audit_canonical.py

The script imports the canonicalizer being audited, but its independent
orbit truth is computed by direct literal rotation actions. Frozen output
is in audit_canonical.json. No git mutations occur.
