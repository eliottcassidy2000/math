# Physical slot covers: Hadwiger–Nelson hinges and the marked fourth star

2026-11-26. **EXACT** contact and overlap tests; the full unrelaxed covering
question is **open here**. The earlier signed-lattice and positive fractional
results are strengthened, not repeated: every individual contact survives
every proper two-center projection, and even adding every clique inequality
does not eliminate the fractional solution.

## 1. The actual board

Keep the graph on 3,250 vertices from [the fourth-star repair](NESTED-CORES-AND-RIGID-FOURTH.md)
and keep the marker `(1,1)--(2,7)`. The source already has the parent and its
three neighbors in center branches `(56,0,28)`. After the marker there
remain 52 exterior branches and 52 unused radius-two slots at each center.
A contact in exterior branch h at sheet y occupies four rows:

    h,  y at center 56,  u_h(y) at center 0,  v_h(y) at center 28.

The pivot coordinate is literally y. The new checker builds the old graph,
checks girth at least five, tries every branch and sheet independently by
distance, and checks **all 1,514,670 pairs** of the 1,741 surviving contacts.
Two contacts conflict if and only if their rows overlap. Indeed any new
short cycle uses either one new edge, ruled out by individual legality, or
two new edges incident at the source. In the latter case the two target
vertices either were adjacent or had a common neighbor; the checker
recomputes that predicate directly from the old graph plus marker. This
equivalence means the physical full-star problem is exactly an independent
set of size 52 in the contact conflict graph, or an exact four-partite
cover. A solution would add one degree-57 vertex; it does not finish the
whole Moore graph. The old 51-contact patch, marker included, still stands.

## 2. Every contact survives deletion of any one center

For each of the 1,741 contacts x and for each choice of one omitted center
family, there exists a choice of 52 actual additional contacts containing x
with the branch and other two center families covered exactly once. This
gives all **5,223** pinned projection witnesses. The omitted center is
allowed to have repeated hits; those witnesses are not patches to the
physical graph. In particular an individual contact cannot be discarded
because it has no extension in some complete two-center projection.

The deterministic sparse exact-cover run produces explicit choices and an
independent checker verifies every witness. No conclusion of unsatisfiability
is inferred from either a search cutoff or a solver timeout.

## 3. A physical hinge in two branches

The inherited [type-M study](../hn_type_m_hinge_rotations_and_forced_radii_20261014.md)
found that just forcing a radius was not enough to glue two coloring pieces.
Even the distinguished monochromatic sets did not capture all the unit
contacts. The later [six–seven–eight boundary study](../hadwiger_six_seven_eight_20261025/SYNTHESIS.md)
gave exact gadgets with full pairwise shadows and different pinned extension
counts. This suggests retaining the actual boundary assignment rather than
its separate feasibility flags. Compare also [Parts' graph-minimization
method](https://arxiv.org/abs/2010.12665), which preserves the chosen predicate
through an explicit checker. The Moore contact graph is not asserted to be
a planar unit-distance graph.

On this board, allow sheets 11 and 28 on branches 3 and 12:

| Contact | Center 56 | Center 0 | Center 28 |
|---|---:|---:|---:|
| (3,11) | 11 | 21 | 30 |
| (3,28) | 28 | 55 | 33 |
| (12,11) | 11 | 55 | 47 |
| (12,28) | 28 | 44 | 30 |

For equal sheets the contacts collide at center 56. For the two crossed
pairs, the collisions are at 28 and 0, respectively. So *each* two-center
projection allows a pair, yet no pair works on all three. The four choices
are legal individually and the four cross-collisions are actual physical
common neighbors. This is a restricted-domain obstruction, not a claim
that those branches have no other choices or that the full board fails.

There are **241** such rectangles on 32 branch pairs, with at most 18 on a
pair. The transport description is useful: for branches h,k and unequal
sheets a,b, one orientation of the obstruction is

    b = v_k^{-1} v_h(a),  a = u_k^{-1} u_h(b),

with all four contacts supported. Thus a mixed return map has a fixed
point even though no single choice of the two sheets satisfies all three
margins. The four contact variables form a clique, so their sum is at most
one in an integral model, a valid inequality involving two different
branch fibers and all three center families.

## 4. Even all clique inequalities leave a positive solution

The board is linear as a hypergraph: two distinct contacts occupy at
most one common row. Within a branch the three charts are permutations.
Between branches two common center rows would have produced an old
4-cycle. For any clique not contained in a single row, at most three of
its contacts can be in one row r: pick a clique contact outside r, which
can intersect those contacts only through its other three parts, at most
once per part by linearity. Fixing one contact in the clique thus bounds
the total size by 1+4(3-1)=9.

On the actual board the bound is stronger. Enumerating mixed triangles
and extending each by intersecting neighbor sets gives the following
numbers of inclusion-maximal mixed cliques:

| Size | 3 | 4 | 5 |
|---|---:|---:|---:|
| Count | 6,948 | 88,248 | 822 |

There are none larger. Every mixed clique contains a triangle with no
common row: take two contacts and a third not on their shared row.
Therefore this enumeration does not need to traverse subsets of the
large ordinary row cliques.

The rounded weights from the [previous exact certificate](OPTIMAL-55-DETOURS-AND-PETERSEN-SHELL.md)
have maximum less than 0.039217 times the common row total. Using the
same exact cofactor correction preserves positivity and keeps every
weight strictly below 1/25; this comparison is checked with integers.
All single-row clique inequalities follow from the row equalities.
Every mixed clique inequality is strictly slack, even using the general
bound of nine. Thus even completing the conflict graph and imposing
**all** its clique inequalities cannot certify that the target has no
exact cover. They can still be useful for propagation in an integer
solver. A future negative result must go beyond this relaxation, for
example combining incompatible partial matchings on a shared boundary
or producing a proof-checked exact-cover refutation.

## Replay

From this directory:

```sh
python nested_core_check.py --board > board.txt
c++ -O3 -std=c++17 slot_dlx.cpp -o /tmp/slot-dlx
/tmp/slot-dlx project 25 71 8000 < board.txt > projected.txt
python slot_interface_check.py board.txt projected.txt
python slot_mixed_cliques.py board.txt
python marked_margin_certificate.py < board.txt
```

The project search emits one witness per forced-contact/omitted-center
pair. `slot-dlx full` is a bounded search for the simultaneous assignment,
and `slot-dlx omit-all` looks for patches missing only one branch. An
empty output from either is not an infeasibility certificate. Bounded
full-board attempts with exact cover and binary optimization did not
produce a physical size-52 cover in this session.
