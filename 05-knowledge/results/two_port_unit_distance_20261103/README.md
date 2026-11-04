# Two-port reconstruction and its unit-distance test

**2026-11-03. Universal two-port reconstruction remains OPEN.** This packet
proves a stronger all-order reconstruction bound, establishes exact two-port
faithfulness through order9, and tests the readout mechanism on fully realized
unit-distance colour boundaries. No six-chromatic graph is claimed.

## Exact progress toward reconstruction

For a named tournament T, the observer is its Hamiltonian-path count after
adding ordered transitive ports with arbitrary named attachment arcs. Labels
are aligned across every query. The key new identity is

    sum_(alpha,beta=+-1) alpha*beta H(T_(alpha,beta))
       = 2 H((T-v)+one collapsed port).

Two adjacent ports have identical external relations; alpha,beta vary their
arcs through old vertex v. The character cancels every path except those
using the consecutive blocks x,v,y or y,v,x. Both collapse to the same
deleted-core path. Thus **four two-port queries recover every one-port
response of a named vertex deletion**. See [R1--R4](reconstruction/PROOF.md).

Iteration gives a proof for every order: if m-port responses collide, every
affected vertex has at least m changed incoming arcs and m changed outgoing
arcs. Positive circulation must survive every deletion of at most m-1 old
vertices. A parity argument additionally excludes odd-order converse pairs:
their ordered two-path-cover matrices cannot be symmetric.

Consequently **floor(n/2) ordered ports suffice for every n>=4**, improving
the recovered n-1 bound. At n=3 two ports are necessary and sufficient;
n=2 needs one, and n<=1 needs none. This is a sufficient bound, not an
optimality claim at all orders. The explicit decoder retains endpoint
tomography and has an exponential finite search; no complexity collapse
is asserted. Independent [algebra proof and controls](algebra/PROOF.md).

Full two-port response has Boolean bidegree at most(2,2), so exactly keeping
the constant, singleton and pair masks in each port gives a sufficient
**[1+n+binom(n,2)]^2-query** battery. This retains the contributions from
ordered three-path covers; the proof does not discard them in favour of a
pairwise score or assume the gradient is globally injective.

## Complete finite result through nine vertices

The recovered two-port battery already determines(H,S,E,C): the Hamilton
count, both named endpoint vectors, and the complete ordered two-path-cover
matrix. Exact canonicalization of this sidecar separates **all191536
order-nine tournament classes**, and every sidecar automorphism preserves
the original tournament. Both tests also pass at every smaller order.

The automorphism check promotes class separation to a statement about
**every labelled tournament through order9**. Thus a counterexample to
universal two-port reconstruction must have at least **ten vertices**.
[Universe, proof of the reduction and receipts](census/RESULT.md).

This extends Session27's finite bound of six, using the complete October26
class bank with its provenance retained. Python and C++ censuses agree
through8; the C++ census covers9. Independent literal permutations and
unrefined canonicalization check every class through6, with six larger
counting controls. No random sampling replaces completeness at order9.

For any remaining candidate pair, retain all of the following: equal H,S,E,C;
equal one-port responses after every single vertex deletion; a changed-arc
graph with at least two incoming and outgoing arcs at each affected vertex
and cycles surviving every one-vertex deletion; and equality of the complete
bidegree-(2,2) battery. The structural conditions alone are not sufficient:
opposite cyclic five-tournaments satisfy the graph conditions but give89/91.

## The geometric test gives complete observation and a forcing obstruction

A private three-edge path u-x-y-v is a genuinely realizable two-vertex
unit-distance probe whenever the existing terminal distance is strictly
between0 and3. It can avoid all unintended contacts with a finite old
placement. Its exact colour-count law is

    Z_(G+ear)(q)=(q^2-3q+3)Z_G(q)-E_G(u,v;q),

where E counts old colourings with equal terminal colours. Joint star-forest
ear queries and integer partition inversion recover the **entire** colour
boundary tensor, not just its pair marginals. All Bell(m) partitions remain
present. See [G1--G5](geometry/PROOF.md) and the independent
[geometry review](algebra/GEOMETRY_AUDIT.md).

The concrete experiment has **37 distinct planar points,57 unit edges and
609 certified nonedges**. Its203 contexts reconstruct all203 six-port
boundary entries of the regular-hexagon wheel, including163 zeros.
The placement, complete tables and a literal three-colouring are saved.

This also proves why this mechanism cannot yield the desired obstruction:
each ear contributes12 or13 extensions at five colours, and20 or21 at six
colours. These factors are positive for every boundary assignment. Any
allowed ear assembly therefore preserves the complete boundary support and
five-/six-colourability of its base. Signed inversion of measurements is
not a physical forcing gadget. [Exact transfer and stopping boundary](TRANSFER.md).

Nor can arbitrary tournament twins simply be drawn as unit probes: two
distinct planar points have at most two common unit neighbours. The lawful
ear construction supplies an actual alternative observer, with its exact
metric scope and its separate positivity limitation.

## Incoming work, audit and reproduction

The packet recovers Session27 and the previous proof-principles packet,
then integrates Session179 and the170807 triangle-return work. The latter's
third-run correction and first symmetry-breaking cost reinforce retaining
actual endpoint/run data; modular spectrum coverage is kept separate from
exact reconstruction. [Incoming comparison](INCOMING.md), [concept board](BOARD.md).

Proof status and finite scopes are kept separately. Mathematical reviews
cover deletion contraction, the odd-converse and triangle-gradient lemmas,
the class-to-labelled census reduction, planar realization and partition
inversion. The actual geometric coordinates are certified with exact
quadratic-field identities and rational isolating intervals, not floating
distance tolerances. Counts and symmetries use exact integers.

From the repository root:

```sh
python3 05-knowledge/results/two_port_unit_distance_20261103/verify.py
python3 05-knowledge/results/two_port_unit_distance_20261103/verify.py --replay
```

The first command checks frozen receipts, artifact/source hashes and the
cross-engine census agreement. The second replays all seven engines normally
and under Python optimization, comparing pre-frozen expectations without
overwriting them. A C++17 compiler is required for the order-nine census;
all other engines use standard-library Python. No software is installed.

The remaining target is explicit: prove two-port injectivity beyond the
finite frontier or supply a fully checked collision at order>=10. For the
six-colour lane, a new geometric interaction with a genuinely zero extension
factor at five colours is required; the complete positive-ear language has
now been ruled out as a route to six-chromatic forcing.
