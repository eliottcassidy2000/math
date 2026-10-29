# Recovered ordered-run witness: all four counts replayed

**FINITE-EXACT independent replay; PROVED explanation of the one-bit fiber.**
Source: commit `eea15f5d968c8e493a809862bbfec3989d09c436`,
`05-knowledge/results/moving_selector_atoms_20261028/THEORY.md`, C2, and
`frames.py`, `parity_fiber`. That source arrived during Session92 publication.
It closes the explicit source limitation in the earlier frozen
`SEAM_CARRIER.md`; the earlier caution was appropriate to its then-available
source pin, not a refutation of these counts.

## Exact objects and observation language

Outgoing-neighbour bitmask rows, vertices numbered from zero:

    A = (10,12,9,16,7),
    D = (2,4,9,3),
    T4 = (14,12,8,0).

A is the strong five-vertex tournament C3[C3,1,1]. Its automorphisms are exactly
identity, (0 1 2), and (0 2 1), each fixing3 and4. This is checked against all
120 permutations. D is strong; T4 is transitive.

The two **ordered all-singleton spanning path covers** of A are

    U=(3,0,1,2,4), V=(3,0,2,1,4).

They are not asserted to be Hamilton paths inside A. Each singleton is a
legal run. The forest complement reverses every arc, hence is A^op.

At cursor index j, observe the whole child rooted at the ordered roles
(first outer port3, last outer port4, current vertex W_j), on A or A^op.
The index is known. There is no retained previous interior bookmark and no
persistent writing operation. Cursor motion is along the ordered run list;
complement exchanges the two sheets. Boundary behavior is identical.

Explicit automorphisms mapping U_j to V_j, while fixing the outer roles, are
identity, identity, (0 1 2), (0 2 1), identity, at indices0,...,4 respectively.
The same maps work on A^op. These pointwise indexed matches imply equality
of every finite adaptive read-only scan/complement transcript: equal observed
histories induce equal next choices, and the same index and sheet again have
equal observations. They do **not** give a single consistent transporter.
The unique map carrying all of U to V is (1 2), which is not an automorphism.
This is not a counterexample to the different single-Hamilton-presentation
K_scan completeness question in source C1.

## Exact assembly and native partner operation

Use quotient C3 with arcs A->B->C->A. Give B vertices5,...,8 and C vertices9,...,12.
For either W=U,V, the global path is

    (W0,5,9,W1,6,10,W2,7,11,W3,8,12,W4).

Its quotient walk is (A,B,C)^4 A. Thus every successive global arc is supplied
by the quotient, irrespective of internal run order. It is a legal Hamilton
path in the same original substituted tournament for either W.

The observed target is the **native partner**: reverse every arc except the
consecutive arcs of this global path. It is not the unchanged ordinary
substitution. This resolves the apparent conflict with THM-3121's ordinary
path-cover substitution law.

| B and C child | even middle order (0,1,2) | odd middle order (0,2,1) | gap |
|---|---:|---:|---:|
| T4, T4 | 2,483,993 | 2,491,877 | 7,884 |
| D, D | 2,899,675 | 2,906,091 | 6,416 |

The checker covers **all six** middle orders in each row, not just the two
quoted examples. All three even permutations give the first value; all three
odd permutations give the second. Both original substituted tournaments and
all twelve partners are strong. In the second row all three children are
strong as well.

## Why one bit is exactly enough on this fiber

The six orders of the three interior vertices are acted on freely by the
three-element rotation group. They therefore form exactly two orbits, the
even and odd orders relative to(0,1,2). Both outer ports remain fixed.
The direction between the first two interior ports detects the two cases:
0->1, while2->0. Cyclic rotations preserve this direction and account for
every order in each sheet. Consequently one bit is both necessary and
sufficient to classify this six-element ordering fiber, and the differing
Hamilton counts prove necessity even for the specified partner-count observer.
This is not a universal one-bit repair for arbitrary run families.

A bookmark at the first interior port also repairs the fiber because its
stabilizer inside the residual C3 is trivial. It is retained information,
not information that can be discarded between operations.

The source proves five is the smallest strong child for this particular
all-singleton, two-distinct-fixed-outer-ports obstruction: C3 is regular and
has trivial vertex stabilizers; a nontrivial order-four tournament symmetry
would cycle three vertices while fixing one, forcing that fixed vertex to be
a source or sink. This statement does not minimize the total assembled order.

## Independent reproduction and audit

    python3 check_recovered_runs.py
    python3 -O check_recovered_runs.py

The standalone script imports no producer or repository code. It reconstructs
both source tournaments, all paths and partners, enumerates Aut(A), verifies
the cursor maps and the absent common transporter, and counts H using both
forward endpoint DP and backward memoized suffix DP. It checks strongness of
all sources and targets. Normal and optimized receipts are identical; they
contain every source row, target row, global path, orbit and count.

An independent Session92 audit uses edge sets and a6/7 meet-in-the-middle
Hamilton counter rather than these recurrences. It reproduced all twelve
cases under normal and optimized Python, as well as anchored cursor profiles
on both sheets and all120 candidate child automorphisms. See the parent audit
lane's `check_recovered_runs_independent.py` and receipts.

## Consequence for the earlier seam theorem

These two paths share the first vertex3. The proved shared-start padding in
`SEAM_CARRIER.md` therefore now applies to each recovered pair: add a new
vertex beating only the common first vertex, prefix it to both paths, and
iterate. Original tournaments stay strong at every larger order, while the
native partners acquire universal sources and retain their respective exact
H counts. This gives infinite strong-source extensions of these numerical
separations. The padded partners are no longer strong, and the enlarged
source need not retain the same three-child substitution shape. No claim is
made that padding preserves a newly enlarged cursor-observation language.
