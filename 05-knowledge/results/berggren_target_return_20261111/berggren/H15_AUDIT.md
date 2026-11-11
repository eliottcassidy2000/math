# Independent audit of the S4 escape and H15 menus

**PASS.** [H15_MENU.md](../reset/H15_MENU.md) correctly separates the
all-height finite-core exclusions, the positive finite return, and the
failure of its particular uniform module lift. The physical root and its
sole incoming neighbour survive both finite maps. No arbitrary-core or
all-future closure claim follows.

## Complete inverse and axes

`audit_H15.py` constructs each literal core/cut/flank tournament without
using the producer's compiler. Nine fresh counts for each of the 272 cuts
recover all coefficients of the inherited exact marked-core formula, for
2448 independent subset-DP probes. The count formula itself remains a
proved inherited dependency; finite interpolation is not being used to
establish the functional form.

The audit then replaces the producer's recursive binary inverse by direct
complete enumeration of the bounded positive-product region and the axes.
For positive A,B, nonnegative terms give `AB <= floor((5993-h)/L)`;
here L=384 for both H15 cores. On each nonuniform axis the positive linear
coefficient c gives `u <= floor((5993-h)/c)`. Every remaining axis has
zero slope and parity term and h=d, so its entire infinite image is
`d*g(u)`. Here d is 15 or 5, neither dividing 5993. Thus no axis is
discarded merely because it has long or noncanonical words.

This independently reproduces empty answer sets for both adjacent H15
cores and all six stated normalized S4 charts. Exterior source/sink
padding changes the presentation size but preserves H; it adds no target
value. The H15 exclusion does not include arbitrary internal singleton
gaps or different seven-vertex cores, exactly as the note states.

## Physical maps and their limits

The explicit B-to-T map retains ten of the twelve source vertices in a
thirteen-vertex target. Fresh induced-matrix equality verifies every
retained arc. The exact trace is

    1997,409,1599,3067,5993.

The old root 0 remains root 0; old in-neighbour 11 becomes target 12 and
is still its sole in-neighbour. Direct counting gives F_end(0)=5998, so
the final positive defect is 5.

The audit also independently counts all four positive ready-return maps
from this same T, including the best displayed trace

    5993,395,1027,2483,5993.

Each preserves both designated physical names, has final defect 2 and
indegree 1, and admits the actual next two births to 17983. The count is
preserved between the endpoints of a return, not at each intermediate
deletion or insertion.

The four target retention optima 10,9,8,8 are explicitly restricted to
this one source, these four fixed strong rooted targets, and the fixed
root/port assignment. Reading `H15_ready.py` and its matching routine
confirms that all source subsets at higher retention levels are tested,
and that the recursive map search prunes only incompatible incidences.
Both fixed endpoints are retained and correctly indexed. I did not
independently rerun this negative embedding census or the inherited
all-order ready inventory. The audit receipt therefore records positive
maps, not a second computational certificate for the optimum. Weak
padded targets and other source charts remain outside that optimum claim.

## Complete single-module obstruction

The audit recomputes all 84 full marked-module polynomials with the
independent backward suffix/run-token method from `audit_reset.py`.
All coefficient lists match. No marked vertex of any of the six S4
templates has second-run coefficient 3402 with all higher coefficients
zero. This criterion is necessary for the proposed target family:
the largest nonzero run coefficient contributes a positive leading
exponential, so a run order above 2 cannot cancel to a dyadic law.

For the explicit mapped placeholder at vertex 1, the complete target
coefficients are `(5993,3684,0,...)`. Together with the source polynomial
they prove, for every m>=1,

    H(T_m)-(3H(B_m)+2)=282(2^m-2).

The retained map expands by matching the whole transitive placeholder in
order, preserving m+9 of the source's m+11 vertices in a target of order
m+12. The fixed-core indegree-one identity gives defects 3 and 5 for
all m. Literal controls m=1,...,4 independently verify the maps and
profiles. The target equality holds only at m=1. Multiple modules,
changing templates with m, or another source family are not excluded.

## Reproduction

Run `python -B .../berggren/audit_H15.py`, or `python -O -B` for the same
checks with Python assertions disabled. `H15_audit.json` contains the
independent inverse rows, all 84 complete polynomials, positive return
controls, and hashes of the frozen producer files. All mathematical
checks use explicit exceptions and pass in both modes.
