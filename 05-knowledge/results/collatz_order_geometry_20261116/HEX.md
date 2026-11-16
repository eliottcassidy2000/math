# The six-clock lift: why the Collatz tree cannot be a Euclidean tiling

## Recovering distinctions from the archive

The earlier [coherence note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md)
noted a literal hexagon, the unit cycle mod 9 under multiplication by 2:
(1,2,4,8,7,5), period six, with negation its half-turn. It also
separates the hexagon intertwining identity from a braided-category
analogy. The [pentagon operator archive](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_pentagon_operator_coherence_20260930.md)
recovers Platonic dual pairs under the cycle operator and the
antipodal twist in the icosahedron fixed point. The next
[archive note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_pentagon_fixed_points_argument_styles_20260930.md)
derives square torus duality and identifies a finite honeycomb
quotient whose hexagon graph is triangular. These are exact
identifications of finite carriers or arithmetic phases; they need
not identify the intrinsic metric of an inverse tree.

## Lift one turn of the hexagon to the actual tree

Use ordinary Collatz: positive even (n) goes to (n/2), odd (n)
goes to (3n+1). Consider only the basin of 1, stopping at first hit
(and omitting the outgoing arc of 1), not the unproved assertion that
the basin contains every positive integer. As a rooted tree its inverse
moves are doubling and, at (z\equiv4\pmod6), the odd inverse
(I(z)=(z-1)/3). There is the cycle exception (I(4)=1); it is absent
in the subtree above 8, i.e., the descendants of 8 in inverse time.

Call a vertex live if (3\nmid y). For live (y), walk over
(2^ky), (k=1,\ldots,6). Each is even and prime to 3. Their
residues mod 18, divided by 2 mod 9, traverse the six-cycle exactly
once. Of the six, exactly 4, 10 and 16 mod 18 admit an odd exit,
and their resulting residues mod 3 are 1, 0, 2, respectively.
Thus a six-step doubling segment has two live exits, one sterile
exit, and continuation along the spine. A sterile exit cannot take
any more odd inverses: divisibility by 3 is preserved by doubling,
and blocks the integrality of the odd inverse.

The two live odd exits from the same (y) are incomparable
descendants: for two different exit positions the two paths initially
double, then one takes an odd child while the other still doubles.
Iterating from 8 gives (2^m) distinct vertices at intrinsic distance
at most (7m+3) from 1. No assumption about integers outside the
basin enters the proof.

The square, triangular and honeycomb edge lattices have at most
(O(R^2)) vertices at graph distance (le R) from any specified
vertex. An injective map of the basin vertices into any of these
lattices with edge stretch at most (L) would send the above
radius-(7m+3) ball into radius (L(7m+3)), contradicting exponential
growth. This holds under every relabeling, and for uniformly
bounded-length embedded edge paths as well as literal embedding
onto edges. It does not rule out drawings with unbounded stretch,
noninjective phase quotients or coordinate coverings.

There is a local way to see the three Euclidean tilings beside each
other: their graphs have polynomial growth despite the readily
available six-cycle quotient of the inverse dynamics. For the
basin, six is the clock period of branch *deposits*, not the
dimension or combinatorial face length of an injective Euclidean
carrier.

## What survives arbitrary relabeling: the end space

An end is an infinite root ray. Two rays are close if they agree
up to large depth. Because every level is finite, the end space
is compact, and common-prefix neighborhoods are clopen.

* An end that reaches a vertex divisible by 3 has only the
  doubling continuation. Therefore it is isolated.
* Every end staying prime to 3 is nonisolated, and even after
  the isolated ends are deleted there are two incomparable
  live macro-exits from every unit vertex sufficiently far along
  the ray. Each continues infinitely by doubling.
* The isolated ends are countable and dense: in every
  common-prefix neighborhood either its ray is already sterile
  or there is a sterile exit on a further six-step doubling
  segment in that neighborhood.
* The subspace of live rays is closed, is nonempty (the infinite
  doubling ray from 1 is live), and is perfect and zero-dimensional.
  By the compact Cantor-space characterization it is homeomorphic
  to the Cantor set. The isolated ends are removed in a single
  Cantor–Bendixson derivative, after which the space is unchanged.

Consequently a relabeling can scramble arithmetic, but it cannot
turn the basin into a single path or a flat edge lattice. The perfect
kernel together with the dense countable terminal-ray ends is an
unlabeled invariant. The word "terminal" here refers to no future
branches, not to a finite ray.

## A family of changes from 3 to 5

For the positive (pn+1) version on odd (n), let (p) be an odd
prime for which 2 generates the units mod (p^2). From a unit vertex
one turn along the doubling spine has length (p(p-1)). The even
unit residues mod (2p^2) are traversed exactly once. There are exactly
(p) residues congruent to (p+1) mod (2p), so exactly (p)
admissible odd inverse exits. Their quotients ((z-1)/p\bmod p)
cover all residues mod (p): one sterile, (p-1) live. At (p=3)
the split is (3=2+1) in a six-clock; at (p=5) it is (5=4+1)
in a twenty-clock. The latter clock is exact because
(2^{10}\equiv-1\pmod{25}), (2^4\not\equiv1\pmod{25}).
This is an inverse-tree statement off any root cycle, not a
convergence theorem for either map. In particular different p-values
can have very different forward attractors even when their branch
mechanisms match.

## An infinite tournament on the basin can encode the full finite spectrum

The growth certificate also proves unbounded depth layer sizes, even
in every live inverse subtree. Its (2^m) simultaneous macroaddresses
use at most (7m+1) individual depth levels, so some level contains
at least (2^m/(7m+1)) addresses. No use is made of ordinary natural-number
ordering or of unknown forward orbits.

Make a tournament on the basin by orienting every strict depth edge
toward shallower depth and freely orienting ties. This always
keeps every pruned Collatz edge. If ties come from one linear
ordering, the finished tournament is transitive and every finite
induced Hamilton count is 1. If ties can carry cycles the picture
changes maximally:

**There exists one fixed depth-respecting completion whose finite
ancestor-closed ideals realize precisely the ordinary finite
tournament Hamilton-path spectrum.**

Enumerate the finite tournaments by (U_1,U_2,\ldots).
For the j-th one choose a new live branch anchor deeper than
every earlier distinguished level (a far doubling descendant of 8
will do). Within its subtree take (|U_j|) vertices at a common
depth (d_j), and orient their ties as (U_j). These choices do not
conflict since the distinguished levels are strictly increasing.
Fill all other same-depth ties by a fixed total order.

Now take the ancestor closure of the distinguished set for (U_j).
Its top layer is exactly (U_j). Below the anchor there is just
one vertex per level, so the earlier tie gadgets cannot reappear
as multiple vertices of the same level. From the anchor up to
but excluding the selected depth, all same-depth layers are
transitive; higher gadgets are absent. Thus every lower-layer
factor has count 1, and the Hamilton count of the ideal equals
(H(U_j)) by the ordinal sum law. The reverse inclusion holds
since each finite ideal induces a finite tournament. One may build
any such individual ideal by a parent-before-child sequence or
in batches of available children, without affecting its final count.
Distinct witnessing ideals need not be nested. The equivalence
does not generate a tournament at a prescribed number unless the
ordinary problem supplied it already; it separates what is in the
Collatz skeleton from what has been placed into ties.
