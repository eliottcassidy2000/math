# Paths, edge traversals, and the integer winding a small graph forgets

**PROVED elementary constructions and boundaries; CITED braid interpretation.**
No equality between an arbitrary tournament's Hamiltonian-path count and
the parameters of a torus knot is asserted. The useful bridge carries
incidence ports and integer winding explicitly.

## Inheritance and the object map

The closest repository mechanism is **THM-2646**,
[braid-three modular central pullback and full-twist knot fibre](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md).
Its same-modular-image unknot/trefoil pair is the canonical hostile: finite
permutation or modular data alone need not retain central winding. Its
exact B3 pullback retains the full modular class and an integer exponent;
the exponent alone does not classify knots under Markov moves.

The incoming [parity-path matching arithmetic](../spectrum_basins_20261020e/arithmetic.md)
supplies an actual matching/SCC-word bijection and marked Fibonacci/Farey
counts. It explicitly does not assign a knot to each tournament or send
black moves to braid moves. The map here has a different, smaller source:
the fixed torus-braid family `(sigma1 ... sigma_(p-1))^q` only.

## 1. Eulerian edges become Hamiltonian vertices only with the ports

For a finite directed graph D with labelled arcs, form its directed line
graph L(D): its vertices are arcs e of D, with an edge `e->f` precisely
when `head(e)=tail(f)`. A directed Hamiltonian word in L(D) is exactly a
linear Eulerian arc word in D. Both are permutations of the arc set, and
their consecutive-incidence tests are identical. This is a bijection of
words, without dividing by reversal, rotation, or automorphisms. A closed
version additionally tests the last-to-first incidence. An Eulerian word
using m arcs has length m in D, and its Hamiltonian image has m vertices
and m-1 edges in L(D).

For an undirected graph the ordinary line graph loses a necessary bit.
The three edges of K_(1,3) give a line graph K3, which is traceable. But
the star has four odd vertices and no Eulerian trail. Its proposed three
edges all meet at the centre: after traversing the first and second, the
second is exited at a leaf, so the third cannot follow. Pairwise shared
endpoints are insufficient.

The exact repair replaces each edge by its two directed darts and asks
for a compatible word choosing **one dart from each reversal pair**.
This is not a Hamiltonian path on the whole dart graph: that would use
each undirected edge twice. The port and pair-selection constraints are
the data that restore equivalence.

Thus the user's star versus triangle-plus-isolate pair simultaneously
illustrates the complement game and the failure of an incidence quotient.
Their edge counts alone cannot determine an Eulerian traversal. Repeated
edges, on the other hand, give the spanning walks of [game.md](game.md),
where connectivity is sufficient.

## 2. An exact torus-family graph with integer carries

Fix integers p>=2 and q>=1. On residues i=0,...,p-1 define the directed
permutation graph

    i -> j_i=(i+q) mod p,
    c_i=floor((i+q)/p),        i+q=j_i+p*c_i.             (1)

Loops are allowed here, and this is not a tournament. Put d=gcd(p,q).
The graph has d directed cycles of length p/d. Along each cycle C,
telescoping the last identity gives

    sum_(i in C) c_i=q/d,       sum_(all i) c_i=q.        (2)

Consequently each component carries the integer pair `(p/d,q/d)`:
its number of vertices and its winding. In the standard torus-braid
interpretation, the endpoint permutation is a shift by q modulo p
(choose the orientation convention for the cyclic generator accordingly),
and its cycles are the components of the closure. Thus d=1 is precisely
the knot case in this family. The braid convention is stated in de Paiva
and Purcell, [Satellites and Lorenz knots, Definition 2.2](https://academic.oup.com/imrn/article/2023/19/16540/6873560).
The component count itself follows directly by following braid strands.

If d=1, the permutation graph has p directed Hamiltonian words, one for
each start, and a closed Eulerian traversal of p arcs. If d>1, it has
neither a spanning Hamiltonian word nor one Eulerian traversal covering
all arcs. These statements follow from its disjoint-cycle decomposition.
They do not say that q is a path count or an Eulerian length.

Indeed, q and q+kp have the same residue graph. All carries increase by k,
and their total increases by kp. For p=2, q=3 and q=5 have the same two
arcs, the same two Hamiltonian words and the same Eulerian length two,
but their carry vectors are `(1,2)` and `(2,3)`. The braid exponent has
changed. THM-2646 already supplies a stronger same-quotient hostile with
proved distinct knot types, so no new topological distinction theorem is
needed to establish the information loss here.

The source is the pair (p,q); the target is an ordered residue graph with
integer-labelled arcs. Equations (1)-(2) preserve both parameters and
component bidegrees. Forgetting labels destroys the unbounded quotient
floor(q/p), while retaining q as a scalar sidecar suffices within this
one-generator family. It does not suffice for general braids: the full
braid word, crossing order, and embedding information are then missing.

## 3. How to use this in the current research, without a false transfer

There are three precise comparisons on the current concept board.

* **Deficiency-three kernels:** the primitive graph retains every external
  defect partner; its bulk double factorial is a proved completion count.
  Just as the residue graph cannot discard winding, a gap-only kernel
  cannot discard external matching incidences. Here the bounded kernel
  theorem actually proves sufficiency for its owner predicate.
* **Black detours:** a path-edge XOR identity is like an unlabelled
  incidence walk. Each intermediate word still needs the right orientation
  and a blackness certificate. Endpoint-changing detours retain both tests;
  relocation inside one tournament fibre is not itself a black move.
* **Collatz excursions:** the structurally similar candidate is a finite
  residue state plus an unbounded exact affine carry. Equations (1)-(2)
  show that a finite graph can carry unbounded arithmetic exactly. They
  supply no well-founded decreasing rank for Collatz. The missing step is
  a cocycle inequality surviving every precision reset and paying the
  original source threshold, not mere ability to encode the integer.

This last comparison sharpens the earlier “Beyond 27” family direction:
families sharing a residue graph must be distinguished by their affine
carry and reset history. A graph is allowed to be small; its arithmetic
labels need not be bounded. The same principle does not by itself prove
eventual certification of any unresolved fixed integer.

The decisive next construction is a **port-labelled transport relation**:
retain endpoint incidences, exact orientation guards, and an additive
integer correction on each local move; then test whether composition
preserves the desired predicate. A quotient may be taken only after that
test. This is a research target, not an established common category for
tournaments, knots and Collatz.

`game_check.py` checks the star hostile, all 2,900 pairs `2<=p<=30`,
`1<=q<=100`, the telescoping winding on every component, and independent
permutation counts for `p<=7,q<=14`. Those finite controls supplement
the elementary all-order proofs; they are not topological knot tests.
