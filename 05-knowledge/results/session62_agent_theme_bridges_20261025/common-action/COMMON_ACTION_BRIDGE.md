# One permutation action, two different fixed-space questions

**Status:** PROVED elementary deduction linking this session's two proof lanes;
independent proof audit PASS. It transports a specified action, not a
coloring predicate or tournament path.

Let A be a finite odd-order permutation group on a nonempty n-element gap
set G. Let r be its number of orbits. Over any field, the common invariant
subspace of the permutation module V=K^G has dimension r, because invariant
functions are exactly those constant on each orbit.

The clone lane proves that the subset-orbit polynomial has vanishing order

    d = min_(g in A) dim V^g = min_(g in A) #cycles(g on G)

at t=-1 (the polynomial proof is over the real/rational coefficient field).
This asks about the smallest fixed space of ONE group element. By contrast,
choose a connected finite gluing graph Gamma with cycle rank beta and
holonomy image A on the SAME permutation module. Rank-nullity gives

    dim H^1(Gamma;V) = n(beta-1)+r.                    (1)

This asks about the simultaneous fixed space of the WHOLE holonomy group.
Equation(1) holds in every coefficient field; its proof uses the dimension
of C^0, C^1 and common flat sections, not semisimplicity.

For the regular A=C3^k action on n=3^k gaps, k>=2, choose a bouquet of k
loops whose holonomies are the standard generators. Equivalently subdivide
those loops into triangles meeting at one vertex to obtain a simple graph.
Then

    r=1, d=3^(k-1), beta=k,
    dim H^1(Gamma;V)=3^k(k-1)+1.

Over F2, n is odd. The augmentation p:V->F2, p(x)=sum x_g, splits by the
constant all-ones vector. Its kernel E has dimension n-1, and E^A=0 for
this transitive odd action. Hence the obstruction forgotten by retaining
only augmentation has dimension

    dim H^1(Gamma;E)=(3^k-1)(k-1).                    (2)

At k=2, the SAME nine-gap action therefore has one common invariant
direction, cancellation order3, full gluing-obstruction dimension10, and
augmentation-forgotten dimension8. These are not four names for one rank.
The clone lane realizes the action by an explicit24-vertex tournament,
so the comparison does not rest only on an abstract permutation example.
The added gluing graph is nevertheless an additional construction, not
an intrinsic graph of tournament moves or a geometric exterior.

The three-pair local-system theorem is the small permutation-module case:
V=F2^3, E=ker p. If holonomy has r orbits on the three pairs, then
E^H has dimension r-1, giving lost dimension2(beta-1)+r-1. This recovers
2beta,2beta-1,2beta-2 for trivial, transposition and transitive holonomy.

Connection contract: start with the full permutation action; send it to
(a) an orbit enumerator and (b) a chosen graph local system. Preserve the
actual matrices, not just n or r. The polynomial forgets which generators
occur along which cycles; the gluing problem additionally requires that
holonomy representation and a specified coefficient field. The regular
C3^k family is the decisive hostile to replacing minimum individual fixed
space by common fixed space. No bound on plane coloring, general tournament
connectivity or counting complexity follows from these dimension identities.
