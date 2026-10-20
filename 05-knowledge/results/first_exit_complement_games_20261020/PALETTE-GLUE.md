# Compatible output gluing of two incoming palette gates

**PROVED elementary composition lemma; FINITE-EXACT placement.** This
integrates incoming commit411c1fd35a's
[six-colour palette gate](../kernels_games_detours_20261020f/forcing.md).
The16-point gate and its nine port triples are that packet's result;
they are not independently invented here. The new object is an exact
two-copy composition preserving all first-copy ports.

## Finite avoidance while preserving an actual root

Let F,H be two finite point sets in a complex field containing
`sqrt(-3)`, each with a distinguished point identified with0. Assume
each set has distinct points. For positive integers t use

    rho_t=(1-3t^2+2t sqrt(-3))/(1+3t^2).

These are distinct units. There is a t among
`1,...,4(|F|-1)(|H|-1)+1` such that F and `rho_t H` meet only at0,
and there is no new unit-distance pair between their nonroot vertices.
Thus their induced unit-distance graph is precisely the one-vertex union.

For nonzero a in F and b in H, the condition
`|a-rho_t b|^2=e`, e0 or1, becomes a polynomial in t of degree at most2
after clearing `1+3t^2`. It is not identically zero: as rho varies on
the circle, the real part of `a conjugate(b) conjugate(rho)` is
nonconstant because a,b are nonzero. Each pair therefore forbids at
most four t's for coincidence or unit distance. Count forbidden values
and choose the first surviving integer. Root incidences are preserved
by rotation, and F is untouched. The proof is elementary and does not
claim a new general graph-amalgamation theorem.

## Exact placement and relation

Use the incoming16-point gate with output p0=0, regarded as complex
coordinates in `Q(sqrt(-3),sqrt(-11))`. Its nonzero coordinates and
nine port triples are reproduced with attribution in `palette_glue.py`.
The first successful parameter is t2:

    rho_2=(-11+4 sqrt(-3))/13.

The union has31 distinct vertices and62 induced unit edges. All465
distances are checked in exact rational field arithmetic. The original
gate, its output and all nine input ports remain at their original
coordinates; the second gate introduces fifteen new vertices.

Recall the incoming gate's conditional relation: if its three port
triples have three distinct colours forming A, its output can have
exactly the colours in A, with24 extensions for each colour. The two
copies only share their output. Hence for palettes A and B the new
output multiplicities are

    count(output=c)=576 if c belongs to A intersection B,
                    0 otherwise.                     (1)

All two-copy boundary assignments of this form are proper: each copy's
boundary is proper, their ports are distinct, and there are no cross
edges. The code independently counts assignments on the full31-vertex
graph for fixed A={0,1,2} and every one of the20 possible B's. It obtains
(1), without using the product formula inside the colouring recursion.
Colour permutations cover every A as well.

For B={0,3,4}, the output is forced to0 with576 extensions. For
B={3,4,5}, the proper boundary assignment has no extension. The whole
graph is nonetheless four-colourable: colour each copy in four colours,
align their root colours by a permutation, and combine. A literal
four-colouring is in the receipt. These are externally conditioned
lists, not unconditional equality forcing in every six-colouring.
The incoming18-point gate is a smaller way to force one output under
external colours. Our added operation is independent two-palette input
and exact compatible composition, not size optimality.

## Two complementary kernels

On the20 three-subsets of six colours, put v_A for the incidence vector.
The shared-output kernel is `576 v_A^T v_B`, a positive semidefinite
rank-six Gram matrix. Indeed `V^T V=6I_6+4J_6`; its nonzero eigenvalues
after the factor576 are17280 once and3456 five times.

By comparison, joining two available-colour outputs by an edge uses
`v_M^T(J_6-I_6)v_N`, as in our19-point hinged wheels. Equality and
inequality split the unconstrained count:

    |M intersection N| + [|M||N|-|M intersection N|]
      = |M||N|.

This is a literal composition law with named colours. It explains why
the same small interface can yield a PSD or indefinite kernel depending
on the physical gluing. There is no transfer of either sign property
to the Collatz ROOT operator without its own source-preserving map.

## Remaining obstruction

The three monochromatic input triples are still prescribed from outside.
Nothing here constructs a global geometric synchronizer enforcing them.
The entire construction also remains in the known four-colour carrier.
Further work must change both the available context and the relevant
carrier while retaining the actual input relation. Replacing triples
by abstract shared colour names would discard that obligation.

Reproduce `python3 palette_glue.py`; `python3 -O` gives identical output.
The receipt is `palette-glue.json`, with exact points, all edges, the
four-colouring and20 complete output rows. Review is same-agent; the
incoming gate's separate audit remains separately attributed.
