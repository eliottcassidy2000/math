# Independent audit: hierarchical attachment kernels and indexed flags

PASS on the structural arguments described below, with source pins recorded
in the final audit receipt. The construction concerns actual tournaments and
specified named-cell automorphism/path carriers, not every native metagraph.

## Uniform-depth structural checks

The full internal automorphism claim W_k is justified, not merely assumed:
vertices in different top children have exactly (r-1)/2 common outneighbors,
while same-child vertices have at least r. Thus the top three children are
intrinsically recovered, and automorphisms rotate them cyclically and recur
inside each. This proves the iterated wreath group and its order.

The prefix restriction q_l is surjective. Its kernel is the independent product
of3^l subtree groups W_(k-l), giving order3^((3^k-3^l)/2). For l<k all independent
bottom-triple rotations remain, so exactly3^(k-1) actual vertex owners are needed.
The l=k kernel is trivial. The l=0 block is uniform and cannot justify an
informative-attachment statement. Isolation likewise retains the full internal
group, not the nonuniform attachment kernel.

For regular outer orientation and uniform depth, incoming/outgoing block degrees
cancel exactly, giving total degree(m3^k-1)/2. The cell-preserving global group is
a fiber product over one W_l frame, with exact stated order. Its owner base is
m3^(k-1) below terminal visibility and3^(k-1) at full visibility. The proof of
sufficiency is explicitly for named cells; the full uncolored tournament may
have additional automorphisms. All lower bounds still apply to that larger
group because the constructed kernels are actual automorphisms.

The internal lexicographic path is directed because its first changed digit
increments by one. Reversing the outer Hamiltonian cell order gives correct
boundary arcs when every relevant visibility depth is at least1. The same
argument does not automatically apply across a depth0 boundary. Ordered
spanning paths have trivial pointwise stabilizer, so the orbit size is exactly
the group order. Indexed owner and arc observations are correctly interpreted
as pointwise stabilizer constraints on actual path positions.

For uniform depth, the proposed matching of successive bottom triples along
the actual path attains ceil(t/2) indexed arc flags. The lower bound is the
independent bottom-rotation support count. It is orbit-relative, not a claim
that these flags identify the path among all Hamiltonian paths.

## Mixed depths and the necessary regularity guard

A chart can be specified by one rotation a_w at every INPUT prefix w. Restriction
to depth l fixes exactly the parameters at depths less than l. Consequently,
for threshold graph G_d retaining edges of depth at least d+1, each input prefix
of depth d contributes one C3 parameter per component of G_d. The exact group
order is3^(sum_d3^d c(G_d)). The parameters count independently even though their
composition is nonabelian. Isolated cells are counted as graph components.

Let G_full retain depth-k edges. One independent bottom rotation survives per
pair (component of G_full, bottom-prefix address). These are the exact ownership
slots. Marking one owner per slot kills every higher rotation too, so the exact
base is3^(k-1)c(G_full). This includes depth0 and disconnected support correctly.

Arbitrary mixed depths do NOT preserve global regularity merely because the
outer tournament is regular. With block degree d_e=3^(k-l_e), the required
condition at each cell is equality of the outgoing and incoming sums of
2d_e-3^k. The uniform-depth case satisfies it. An odd circulant outer tournament
with depths constant on undirected step classes also satisfies it, since every
step supplies one incoming and one outgoing edge per vertex. Its threshold
component counts are gcd(m,{steps with sufficient depth}), with empty-step
convention gcd(m)=m. Independent tests verify this regular arithmetic family,
including an actual81-vertex degree40 example.

## Exact arc-slot edge cover and the order63 hostile

An indexed arc covers the slots of its two endpoints. Each bottom triple also
supplies a consecutive internal arc, hence a singleton/loop cover for every slot.
A minimum edge cover in this loop-augmented slot graph has size t-nu, where nu
is the maximum matching size in the nonloop graph. Upper bound: take a maximum
matching then one loop per uncovered vertex. Lower bound: an inclusion-minimal
edge cover is a disjoint union of stars and singleton loops; one edge per
nontrivial star gives a matching of size at least t minus the cover size.
Thus merely halving the owner budget is insufficient for arbitrary mixed-depth
slot arrangements.

The order63 example is independently constructed with outer circulant C7,
internal T2, depth2 on the directed triangle{1,3,5}, and depth1 otherwise.
Every vertex has outdegree31. The reverse-cell lexicographic path is verified
arc by arc. There are15 ownership slots. Its nonloop slot graph is bipartite
with parts6 and9 and has an explicit maximum matching of size6; therefore the
exact indexed arc cost is9, not ceil(15/2)=8. Nine actual flag positions are
retained in the independent receipt. The named-cell group order is3^16; neither
that order nor the path observation upper bound is claimed for arbitrary full
uncolored automorphisms or paths outside the specified orbit.

## Independent evidence

- Direct tournament/generator tests through internal depth4, plus actual regular
  three-cell hosts at depths2 and3.
- All81 declared W2 charts and exact local kernel base searches.
- All531441 triples of W2 charts and all27 mixed depth patterns on three cells;
  exact group counts and a separate parameter-owner set-cover dynamic program.
- All120 circulant depth patterns for odd outer orders3,5,7,9 at depth2, with
  threshold-component gcd and degree-balance controls.
- Independent63-vertex hostile, bipartition, maximum matching and actual flags.

Executable tests supplement the all-size proofs. They do not enumerate W3 or
arbitrary tournament isomorphism classes. No source or Git mutation was performed
by the auditor.
