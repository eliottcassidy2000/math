# Unordered pair transport: the exact obstruction lost by Fano completion

**Status:** PROVED by the elementary arguments below; independent proof audit PASS. This extends
Session61's ordered-pair theorem, not the plane six-color bound.

## 1. The action and coefficient system

At each vertex of a finite graph Gamma, take six named local palette slots
partitioned into three pairs. Unlike Session61, an edge may permute the
pairs as well as reverse their endpoints. Choose a temporary local order
of the three pairs and an endpoint frame. An oriented edge e:u->v then
has transport (sigma_e,omega_e) in F2^3 semidirect S3. The permutation
sigma_e transports pair coordinates; omega_e records endpoint reversals.
Reverse orientation uses the inverse affine map. No assumption that these
edge transports have trivial cycle product is imposed.

Compatible endpoint choices b_v in V=F2^3 satisfy

    b_v + sigma_e b_u = omega_e.                       (1)

The linear parts define a rank-three local system V_sigma on the graph.
There are no two-cells, so any omega is a one-cocycle, and (1) is solvable
iff its class in H^1(Gamma;V_sigma) vanishes. Changing temporary local
frames acts by gauge transformations, preserving this existence predicate.

Each local pairing has two Fano completions, as proved in Session61.
A pair permutation preserves their parity rule; endpoint reversal x changes
the completion bit by p(x)=x_1+x_2+x_3. Hence completion choices eta_v
satisfy exactly

    eta_v + eta_u = p(omega_e).                        (2)

The quotient coefficient system is constant F2, even when V_sigma is twisted.
These statements concern the local pairing/Fano construction. They do not
turn palette slots into graph vertices or conic components.

## 2. Equivariant splitting and exact obstruction dimension

Let E=ker(p), of dimension two. Since three is odd, the invariant vector
j=(1,1,1) has p(j)=1. Thus

    V_sigma = E_sigma direct-sum F2*j                 (3)

as LOCAL SYSTEMS, not merely as vector spaces at a chosen vertex.
The splitting is invariant under all pair permutations. Consequently

    H^1(Gamma;V_sigma)=H^1(Gamma;E_sigma) direct-sum H^1(Gamma;F2),

and observing only the Fano-completion class forgets precisely
H^1(Gamma;E_sigma). Vanishing of the scalar class is necessary, not in
general sufficient, for a full frame lift.

Assume Gamma connected. Write beta=|edges|-|vertices|+1 and let H<=S3
be the holonomy subgroup, obtained after choosing a spanning tree and
transporting the cycle permutations to a root. Its conjugacy class is
independent of the tree and root-frame convention. Put f=dim(E^H).
Euler characteristic of the two-term cochain complex gives

    dim H^1(Gamma;E_sigma) = 2(beta-1)+f.              (4)

Indeed dim C^0=2|vertices|, dim C^1=2|edges|, and H^0 identifies with
holonomy-invariant vectors in the root fiber. Rank-nullity proves (4).

The complete S3 subgroup list yields:

    H trivial:       f=2, lost dimension=2 beta;
    H order two:     f=1, lost dimension=2 beta-1;
    H contains C3:   f=0, lost dimension=2 beta-2.

A transposition fixes {(a,a,0)} in E after relabeling. A three-cycle fixes
only constant triples in V, and the only even constant triple is zero.
The final row covers C3 and S3. The displayed dimensions are nonnegative
for realizable graph holonomy: a tree has trivial holonomy, and a single
cycle cannot generate S3.

This closes the exact local-system qualification left in Session61.
The old value 2 beta is recovered only for trivial action on E; it must
not be carried over unchanged when pair permutations are allowed.

If the full obstruction vanishes, the number of compatible full endpoint
frames is 2^(1+f); the scalar completion equation has two solutions.
For each prescribed compatible completion, exactly 2^f full frames remain.
For disconnected graphs the formulas add over components, rather than
inserting one global holonomy subgroup.

## 3. Hostiles and boundary cases

Take a simple three-cycle with identity transports on two edges, and
permutation sigma on the closing edge. Tree gauge reduces the obstruction
to a class in E/(I+sigma)E.

* sigma=id: two forgotten bits remain. The old gain(1,1,0) gives a false
  positive for the scalar parity observer.
* sigma=(12): one forgotten bit remains. The image of I+sigma on E is
  span(1,1,0); a gain(1,0,1) is not in that image although its parity is0.
* sigma=(123): I+sigma is invertible on E. Every even gain can be repaired,
  so completion compatibility is sufficient on this one-cycle graph.
  This is not a universal sufficiency claim on higher-cycle graphs.

A graph with two independent cycles both carrying the same three-cycle
has f=0 and lost dimension2, not zero. Nontrivial holonomy can reduce a
fixed finite obstruction dimension; it does not eliminate every obstruction.

The zero-gain twisted case also exposes a terminology boundary. Equation(1)
asks for a section of an affine endpoint-choice bundle; it can have a
section even if sigma has nontrivial holonomy. It does NOT assert that all
six palette labels can be globally identified by the prescribed transports,
which would require trivial full permutation holonomy. A section and a
trivialization are different targets.

## 4. Connection contract and scope

Source: marked local six-palettes with paired slots and affine edge transport.
Map: parity V_sigma -> F2, followed by the Fano-completion interpretation.
Preserved: the scalar obstruction and its compatible completion choices.
Lost: E-valued cycle transport, with exact dimension (4).
Sidecars: graph, local pair orders, endpoint frames, every edge permutation
and gain, and the distinction between a section and a trivialization.
Cheap decisive tests: identity/transposition/three-cycle monodromy on the
same triangle, and two cycles with the same three-cycle holonomy.

This is a precise bridge among palette transport, Fano completions and
local-system cohomology. The genus-two matching dictionary and the marked
six-root canonical character remain inherited interpretations of the local
parity map. No claim that arbitrary transports are realized by a planar
unit-distance exterior, and no six-coloring of the plane, follows.

## 5. The eight-state cube is a genuine common carrier here

For a fixed pairing, each bit vector x in F2^3 chooses one endpoint from
each pair, hence one of the eight transversal triples. The two possible
Fano completions select respectively the four even and four odd triples.
Thus their two completion states are exactly the two bipartition classes
of the cube graph on these eight endpoint-choice states. The full action
x -> omega+sigma x preserves cube adjacency and acts on its bipartition
by eta -> eta+p(omega). This is precisely the same quotient as (2).

The concurrent conic-component theorem supplies a second marked cube:
a free hourglass code acts by translations on eight conic labels, with
three independent direction differences d1,d2,d3. Choosing that direction
basis identifies its component LABEL SET with F2^3. Under this marking,
the endpoint-choice cube and the conic-label cube can be equipped with
the same ABSTRACT affine coordinate action. The inherited geometric
translation action agrees with its translation subgroup; coordinate
permutations supplied by the marking are not thereby realized geometric
automorphisms. Their bipartition quotient forgets E exactly as above.

This is a correspondence of marked combinatorial carriers, not a proof
that every affine cube automorphism lifts to a surface automorphism or a
Euclidean isometry. Component intersections, root values, complex geometry,
and actual exterior coloring constraints remain separate sidecars.

In particular this cube is NOT the other agent's eight-state Boolean path
automaton at six colors. That automaton has two globally fixed semantic
states (dead and free) under palette permutations. The cube's translation
action is regular; every nonzero translation has no fixed vertex. The
actions and target predicates must not be identified merely because both
sets have eight elements.

## Sources and reproduction

The local-coefficient formalism is classical: Allen Hatcher, *Algebraic
Topology*, Section3.H, especially the module/bundle definitions and
cellular cohomology discussion (printed pages327–334), primary author PDF:
https://pi.math.cornell.edu/~hatcher/AT/AT-doublepage.pdf .
The relevant definitions were inspected2026-10-25. The graph calculations
here are direct rank-nullity proofs; no new general cohomology theorem is
claimed. The new-to-route application is the explicitly marked unordered
palette/Fano/cube transport and its obstruction table.

Run `python3 check_holonomy.py` and `python3 -O check_holonomy.py`.
The finite universe is every connected labeled simple graph on1–4 vertices,
every assignment of S3 edge permutations, all48 affine one-cycle cases,
and all48 affine cube maps. The program compares independent holonomy
closure with full V/E/scalar coboundary ranks. Finite tests corroborate,
but do not replace, the arbitrary-graph proof.
