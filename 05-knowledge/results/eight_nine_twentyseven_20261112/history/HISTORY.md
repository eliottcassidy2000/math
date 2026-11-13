# Why 8 and 9 are actual ports, and why 27 needs a retained coordinate

Date: 2026-11-12. **INHERITED PROVED:** the square-reflection and Collatz
clock mechanisms cited below. **PROVED HERE:** the selective-port factorial
filtration in [PORT_FACTORIAL.md](PORT_FACTORIAL.md). The numeral pattern
guided the choice of objects; the maps below retain their different targets.

The concurrent [Session295 construction atlas](../../session295_construction_atlas_20261113/atlas.md#new-selective-port-filtration),
read at `origin/main` commit `9138a603e2`, independently proves the same
all-r factorial filtration with the same compound-block mechanism. This
is parallel proof credit, not a claim that the present theorem strengthens
that result.

Board: partial reflections / indexed path positions / exceptional excursions /
free block permutations / one lost ternary digit. The productive move was to
retain the exceptional pieces as whole blocks, then permute only the pieces
whose seams remain uniform.

## 1. The important square-sum 8/9 pair has two precise older meanings

The September theorem **THM-4510**,
`01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md`,
contains the actual path

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9.

It visits 1 through 15 and uses only the adjacent-sum targets 9,16,25.
These are the squares of the Pythagorean triple (3,4,5), with 9+16=25.
For a target t, the operation is the partial reflection x↦t−x on the
host [1,15], with fixed points deleted because loops are forbidden.

The reflection at 16 pairs every host vertex except its fixed point
8=16/2. The two outer reflections at 9 and 25 cover complementary parts
of the host and leave the boundary value 9 unpaired. Thus the union of
these three matchings has precisely the endpoints 8 and 9. Two reflections
compose to a translation, and the alternate-vertex orbit is rotation by
9 modulo 16. This is the actual operation behind the endpoint pair.

The full square-sum graph also contains the unused edge 1--3 of sum 4;
the three-target carrier is a specified subgraph, not a claim that the
entire graph is only a path. The old theorem identifies 15,16,17 as the
Hamilton window supplied by this one consecutive-square triple. The
three-target union branches at 18; the full square-sum graph then has
the familiar three-leaf obstruction.

The same theorem has a second exact role for 8/9. In the sumgraph whose
targets are powers of 2 or 3, the pair 8=2^3 and 9=3^2 is consecutive.
Its two reflections therefore generate translation by 1, yielding the
two-target Hamilton paths at host orders 7 and 8. This is an equality of
arithmetic reflection operations. It is not an assertion that a path in
that sumgraph is a Collatz trajectory.

For arbitrary primitive s²+t²=u² with s<t, the same construction gives
the host [1,t²−1]. Its endpoints are s² and t²/2 when t is even, or s²/2
and s² when t is odd. The actual square labels, which endpoint is a
deleted fixed point, and the host interval are the reusable coordinates.

## 2. A graceful sign gauge turns the path into a signed clock

The July source **THM-2761**,
`01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md`,
uses a bipartite vertex-sign gauge to turn edge sums into oriented
differences. Applied to the selected path, write its vertices x_i and
put y_i=(-1)^i x_i. Then

    y_i−y_(i+1)=(-1)^i(x_i+x_(i+1)).

Thus the square-clock word becomes a signed displacement word. Its
reconstruction is exactly

    x_i=(-1)^i x_0+c_i,
    c_0=0,  c_i=s_i−c_(i−1).

The October extension
`05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/SHARED_CLOCKS.md`
retains shared *indexed* vertices between several such paths. One signed
integer root remains per connected component. Even closed clock walks
test a translation; odd ones pin twice the root. Host inequalities and
collision holes must still be intersected. An actual indexed mark fixes
the remaining root; an unindexed assertion that a Hamilton path contains
a given integer supplies no such information.

This gauge applies to the selected bipartite path, not globally to the
full square-sum graph, which already has an odd cycle at order 15. The
mirror factor in the graceful theorem likewise distinguishes oriented
differences from their absolute values; discarding signs changes the
predicate. The operation shared with our new module result is specific:
mark and retain the pieces whose seams carry conditions before permuting
or identifying the other pieces.

## 3. The asymmetric eight-vertex module has permutation arithmetic

S293's
`05-knowledge/results/session293_cycle_returns_20261113/control/MODULE_IDEAL.md`
gives the asymmetric eight-vertex core with rows

    (56,117,129,22,100,76,141,59)

and H=513. The gcd of its ordered cover counts is 9, but every uniform
context has path count divisible by 27. The extra factor comes from
permuting the actual k−1 exterior excursions between k core runs, not
from an automorphism of the core.

The new [selective-port theorem](PORT_FACTORIAL.md) answers what survives
when uniformity is broken at r actual exterior vertices. Contract every
excursion containing a selective vertex into its adjacent core runs.
The resulting compound blocks keep all exceptional endpoint conditions.
Their remaining free orders give the universal k-run divisor

    ((k−r)_+)! ((k−r−1)_+)!.

This is a literal retained-slot construction. At r=0 it recovers the
uniform excursion factor; at r=1 the entire arbitrary-context calculation
modulo 27 needs only the layers k≤6. In an eight-vertex core the k=8
layer cannot break 27 with fewer than three selective exterior vertices,
although lower layers can escape with a single selective vertex.

The exact top-layer family in the new proof supplies a physical positive
control: for k=5 and a one-vertex selective mask of size 1, its top layer
has 576 paths, divisible by 9 but not by 27. This establishes a real
boundary of the factorial filtration rather than predicting a target
from a gcd alone.

The formula also gives h_6=14400 at k=6,s=1, with v_3=2. This is exactly
the sharp one-port example in Session295: seven is the first universally
vanishing run layer modulo27. Our bounded family receipt stops at k=5;
the displayed all-length formula proves the six-run specialization.

Source→target: indexed exceptional seams→compound blocks→free orbit
sizes. Preserved: every directed adjacency and the complete exceptional
incidence. Lost when only the divisor is retained: the positive menu of
actual counts. No equality between the square-clock affine group and
the module permutation groups is needed for this transfer of method.

An additional literal square-sum operation is recovered by
[Session294's bridge](../../session294_context_27_20261113/README.md#square-sum-and-collatz-exact-bridges-with-different-carriers)
from the older [Session53 graph theorem](../../session53_g27_collatz_cooper_20261024/g27-growth/THEOREM.md).
These are inherited results. In the graph G_n on labels 1,...,n, with
edges joining distinct labels whose sum is a square, the new vertex 27
has exactly the old neighbors 9 and 22: 27+9=36 and 27+22=49. Among G26's
1116 spanning two-path covers, exactly 32 expose 9 and 22 on different
components and can be joined through 27. Each compatible cover gives
one unoriented path. Three of G26's twelve Hamilton paths, counted up
to reversal, instead extend at an endpoint, giving the separate split
3+32=35 at G27. Since 9+22=31 is not a square, the joins do not subdivide
an old edge.

This supplies a precise complementary coordinate. Endpoint compatibility
can make a prescribed join sector empty; free excursion permutations can
make a nonempty tournament sector invisible modulo 27. The first retains
actual labels and exposed endpoints, while the second retains exceptional
blocks and counts whole permutation orbits. A minimum cover number or a
context gcd alone loses the distinction between those two mechanisms.

## 4. Collatz 27 is a two-step clock modulus

The older hexagonal-growth packet's
`05-knowledge/results/collatz_hexagonal_growth_20261102/MINIMAL_CHART.md`
proves a precise finite-horizon result for its SF5 inverse branch/clock
observer. Its first branch menu depends on the unit residue modulo 9.
For a chosen exponent a, the update is

    j_a(n)=(2^a n−1)/3.

Predicting the next menu requires one more ternary digit. The two-step
chart therefore needs the unit residue modulo 27; more generally h
steps need modulo 3^(h+1). The integer register is retained when actual
source labels or native meeting times are requested.

A concrete three-lift control uses the odd roots 7,25,43. All are 7
modulo 9 and have first menu (4,6). Choosing exponent 4 gives

| root | residue mod 27 | child | child mod 9 | next menu |
|---:|---:|---:|---:|---|
| 7 | 7 | 37 | 1 | (2,4) |
| 25 | 25 | 133 | 7 | (4,6) |
| 43 | 16 | 229 | 4 | (2,6) |

The three sources were indistinguishable to the first modulo-9 menu but
not to its continuation. Division by 3 transports the previously hidden
ternary digit into the next observable residue. This is an exact finite
observer statement, not a convergence assertion about the orbit of 27.

The square-sum pair, the factorial divisor 27, and the modulus 27 therefore
pose three distinct retained-coordinate problems: the actual matching
endpoint, the exceptional excursion/slot, and the next ternary digit.
The new theorem resolves the second uniformly for every r. The first
and third provide concrete older models of why an apparently small mark
can change the available continuation without changing a coarse current
summary.

## Reproduction

Run `python3 -B check.py` and `python3 -O -B check.py` in this directory.
The receipt contains 72 generated contexts with 0–2 selective vertices,
independent literal path enumeration for total orders through 8, ten
exact top-layer family cases, the actual 15-vertex clock, and the three
Collatz lifts above. These are targeted controls for the maps and theorem;
the old all-order results are inherited from their cited proofs.
