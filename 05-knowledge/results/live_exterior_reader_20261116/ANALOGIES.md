# What the compressed data forget, and how to read it back

## A literal four-vertex two-switch atom

Keep 0→1, 0→3, 1→2, 2→3 fixed, giving minimum outdegrees (2,1,1,0).
Start with 0→2 and 1→3, and use flags x and y to reverse those arcs.
The four settings give the four tournament classes at order four:
transitive, source plus triangle, triangle plus sink, strong. H=1+2(x+y).
The exact ordered-slot polynomial, including empty slots, is
P(t)=t^4+2(x+y)t^2. The four vertices can each choose a slot;
each directed triangle contributes a further weight 2 and forces its
three vertices into one slot, leaving an independent singleton slot.
There are x+y triangles and no other nontrivial odd-cycle pack.

In particular the choices 01 and 10 have all the same unmarked slot
counts. Nevertheless root 1 separates their weighted profiles:
(h,g,R,S)=(3,1,1,2) versus (3,1,3,3), with the first profile for the
source-plus-triangle. Ordinal sums of k such chiral blocks give 2^k
different tournament classes with the very same P(t)=(t⁴+2t²)^k:
slot runs factor for each fixed label, and the ordered strong components
recover each four-vertex block's choice by cumulative vertex sizes.
A uniformly read block cannot distinguish the choices, but a selectively
marked port can. Conversely, this exact example gives no reason for
the marked changes to be small or for independently chosen ports to add.

The four fixed edges are a C4; the two free edges are the complementary
perfect matching of K4. The other complementary matching choices give
three coordinate charts, not three automatically composable tournament
moves. This is the precise shared combinatorial object with the
[old Fibonacci-box note](../../../07-reflections/fibonacci-boxes-farey-and-the-k4-matching-triple-opus-20260728.md).

[Price's very different ternary tree](https://arxiv.org/pdf/0809.4324)
uses four entries (p−q,q,p,p+q), constrained by two linear identities.
The same hypotenuse is a sum on one matching and a difference on another.
Column flips alone do not preserve the box constraints; Price discards
and recomputes the dependent row. That is a better model of a reader
move than an analogy based on having three branches: the operation
*recomputes* the dependent exterior instead of freezing the two before
changing a contact. The arithmetic matching identity does not transport
the path count. The old note's corrected warning about the modular
free-factor slogan is relevant: a common syntax need not be a conjugacy.

## A genuinely composable interface

For a binary word w store length, its forward and reverse addresses
b(w), r(w), the interaction J(w), and its end bits. Concatenation u v
has exact cross term

    J(uv)=J(u)+J(v)+b(u) r(v)+last(u) first(v)

for nonempty inputs (omit the last term for an empty side). The
geometric term factors; the final correction remembers contact
adjacency. Ordinary binary addresses concatenate by shifts. This is
an associative exact compression without storing every pair of sites;
fixed bit-count or a residue alone is not a substitute.

The [Chair44 discussion](../recursive_interfaces_20261016b.md#6-boundary-composition-recursion-modes-and-the-lawvere-boundary)
isolates the difference between indefinite compatibility and termination.
Unique coarsening literally divides a translational period by two
each time; no nonzero integer period survives indefinitely. In our
reader a persistent nonempty choice set at each scale is not that
argument. An exact analogue would need a physically realized move
that divides a well-defined nonzero error by two while retaining a
continuation and terminates at zero. Neither the tangent bound nor
the exact-gap theorem supplies that move yet.

The [Monster/Baby Monster archive](../monster_levels_square_sum_7_21_20261022.md)
explains {4,2,5,3,7,13} as a valuation mnemonic for a factorization of
|M|/|B| after treating 4 as a label for 13, not as a canonical operator.
Contrast the genuinely nested level 1|2|6 rationality data and the
centralizer/incidence data they retain. The reader analogue is to compare
what a given observation can distinguish, not to transplant a numerical
cycle. The four-vertex hidden bank is an exact instance: forgetting the
mark erases all k directions at once, restoring a mark distinguishes
some of them. The next useful question is what smallest selection of
marked observations recovers an exact target without losing those choices.

### How much compression can be exact?

The archive's [consumer-complete quotient](../../../07-reflections/consumer-complete-quotients-and-natural-addresses-codex-20260902.md)
criterion says a summary can discard two states only when the named
consumer cannot distinguish them. Here a single extra physical read
makes the criterion sharp. At fixed synchronized seed and address width n,
with D=2^n, append one further near-sink zero after the usual reading zero.
Its two observed scores satisfy

    H(w0)-H(w)=2ZD+2g b(w).

This follows from b(w0)=2b(w), r(w0)=r(w), J(w0)=J(w).
So the two scores recover **every address bit**, despite possible
collisions in the original scalar score. An exact deterministic summary
sufficient for both observations needs at least n bits on the full
n-bit family. The associative profile compresses the pairwise
calculation, not the information content of the distinguishing address.

This follows the [delayed-consumer discussion](../../../07-reflections/specialization-before-quotient-and-delayed-consumer-synthesis-codex-20260902.md):
choose the contact or hold all relevant phases *before* identifying
states. The [Rule 30 decorated clock](../../../07-reflections/rule30-three-prize-decorated-clocks-and-observer-mass-root-20260824.md)
makes the same distinction through the discarded inverse-boundary ray,
which selects a particular seed. In contrast, the
[A000568 one-tail](../../../07-reflections/a000568-the-second-quotient-is-the-one-tail.md)
is universally aggregable for its stated Burnside consumer: the count
of one-cycles interacts with the rest only through its number of parts.
Here the two exterior words are not such a freely aggregable tail.

There is even an exact shared smooth chart with Price. Let
u=g b/(ZD), v=1-rev(b)/D. The leading product divided by ZD is
(1+u)(1+v). At target ratio N/(ZD)=2 the partner equation is
u+v+uv=1, Price's relation for the two half-angle tangents. At other
targets the constant changes. This is an equality of real level
equations, not an identification of legal bit-reversed addresses with
all rational half-angle boxes. The tangent construction uses a grid of
rational directions to approximate the level without that false step.
