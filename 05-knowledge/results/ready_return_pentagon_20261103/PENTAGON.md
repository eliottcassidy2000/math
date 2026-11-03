# Five families as an operation diagram, not just a shared cardinality

## P1. A literal pentagon for every prefix length

**PROVED.** Write A=10, B=01 and let E denote the empty tail. At each
k>=1 the five exceptional source types support the cycle

    E -- A -- AB -- BA -- B -- E.                       (P1)

The first two edges append an A and then a B at the end of the chain,
with their specified two-port arcs and the unchanged root. The middle
edge reverses just the arc between these adjacent tail vertices u,v.
All their other incidences stay fixed. Their order becomes v,u, so the
new word is BA. Delete the now-last A to reach B, then delete B to
recover the original source, including every old vertex and old arc.

There is an explicit native implementation of the middle edge. In AB
the Hamilton path

    a, r, u1,...,uk, u, v, b

contains u->v with v nonterminal. The inherited stabilized arc-flip
gadget uses two actual off-spine complements and one temporary vertex,
then deletes that temporary vertex. Peak order is k+6. The script
replays its path guards and complete target relation. This edge changes
H by2^k, from13*2^k+3 to14*2^k+3.

The other edges are literal insertion/deletion operations, so P1 is
not asserted to be a five-cycle in the fixed-order native metagraph.
Its closed tour has identity transport on all persistent vertices.
Prepending11 commutes with every displayed operation: the same diagram
exists at every scale, with the old R growth as its clock. Thus the
five-family observation has become a reusable operation diagram.

As an abstract graph P1 is the four-letter associahedron. One explicit
adjacency-preserving map, in cyclic order, is

    E  -> (((ab)c)d),       A  -> ((a(bc))d),
    AB -> (a((bc)d)),       BA -> (a(b(cd))),
    B  -> ((ab)(cd)).

Each displayed bracket edge is one rotation. This upgrades the old
[September pentagon/coherence analogy](../collatz_pentagon_operator_coherence_20260930.md#4-mac-lanes-coherence-against-terras-and-conway-analogy-the-bijection-typed)
to an actual operation graph on these tournament sources. It does not
establish an associator for general tournament attachment, higher
coherence, or termination of a source machine.

## P2. A genuine involution matching Platonic duality

Interchange the two port names and reverse the port arc to restore the
ordered convention a->b. This sends every letter(A,B) to(B,A). It fixes
E and exchanges A with B and AB with BA. By A1 the joint response C
does not use the port arc, hence is preserved. Strongness, q=s and
the11 prefix are preserved as well. H need not be preserved.

This is a literal involution of the source family, not a tournament
automorphism. As an involution-set it is isomorphic to Platonic duality:

|Source tail|Regular convex polyhedron type|Dual partner|
|---|---|---|
|E|tetrahedron {3,3}|itself|
|A / B|octahedron {3,4} / cube {4,3}|each other|
|AB / BA|icosahedron {3,5} / dodecahedron {5,3}|each other|

This map preserves the `1+2+2` orbit structure; it forgets all tournament
arcs, Hamilton counts and geometric face incidences. There is no claim
that the full symmetry groups agree: our rooted sources are rigid.
The two choices of orientation within each paired row are conventions,
not an arithmetic identification of individual primes with solids.

The geometric five arise from the positive-curvature inequality
`1/p+1/q>1/2` for p,q>=3, which has exactly the five ordered solutions
in the table. The source five arise from the positive joint-response
budget C=2 and the finite-tail law. Both arguments bound exceptional
local configurations before allowing a repetitive part. The shared
method is useful even though the inequalities are different.

## P3. Recovering the Fermat and modular-curve mechanisms

The October [FIVE_LADDER](../fermat_five_pepin_collatz_20261025/FIVE_LADDER.md)
suggested constructibility, Pascal parity and cyclic fivefold symmetry.
The later [GENUS_BRIDGE](../clocks_carriers_genus_20261102/GENUS_BRIDGE.md)
separates the correct geometric objects. Its mechanisms transfer as
follows; none uses a presumed bijection merely because two lists have five
members.

**CITED / inherited.** For odd prime p, all circulant tournaments of order
p have only their translation automorphisms iff p is a Fermat prime:
[THM-4644, Fermat/Wieferich circulant phenomena](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md).
The multiplier result is credited there to Akin, with the prime-degree
automorphism reduction credited to Burnside/Alspach. The mechanism is
the absence of odd multiplier subgroups when p-1 is a power of two.
It supplies a symmetry test for tournament families, not a counting
explanation for our five exceptional tails.

There are five **known** Fermat primes,3,5,17,257,65537; exhaustion is
open. See [Hochster, Fibonacci Numbers and Mathematical Explorations,
chapter12](https://sites.lsa.umich.edu/hochster/wp-content/uploads/sites/1337/2024/11/fib20B.pdf).
The ready-family five, by contrast, is an exhausted classification.

The prime levels with genus-zero X0(p) are exactly2,3,5,7,13. The
classical prime-level genus formula

    g=(p+1)/12-e2/4-e3/3,  e2,e3 in {0,2},

bounds p by13 at g=0; testing the remaining prime residue classes gives
the list. Source: [Snowden, Lecture12](https://websites.umich.edu/asnowden/teaching/2013/679/L12.html).
This is another finite exceptional list produced by a positive growing
term and bounded corrections. Ogg's fifteen concern the different
Fricke quotient X0(p)^+; see [Aricheta's primary exposition of Ogg's
observation](https://celebratio.org/Ogg_AP/article/1157/). The source and
quotient must stay attached to the invariant.

The effective connection to our new return theorem is **precision**.
THM-4644's Wieferich layer and R3 both ask whether a multiplicative clock
lifts at prime-square precision. Our order calculation has valuations
exactly one at7, so it gives the ordinary lift and an additional binary
phase for the mixed2/3 observer. An actual Wieferich failure would alter
that clock; a five-object analogy would not detect it.

## P4. Next definitions and tests

The strongest new organizing object is a **graded source atlas**:
retain the actual source word/core, ordered in-ports, q-s, joint defect,
the exact joint-cover polynomials and the required precision clock.
ATLAS proves that its neutral two-in-port layers have finite boundaries;
RETURN supplies literal edges between different layers.

Three precise open targets follow.

1. **Offset restoration on the graded atlas:** construct operations at
   `11^k00,01` and `10^k00,10` whose output counts are3H+2 and3H+4,
   with a source return inside a specified enlarged grammar. H35 proves
   that demanding a ready return is too restrictive.
2. **The next atom layer:** classify neutral defect6, separating the
   transitive two-port slice from three-in-port atoms and nontransitive
   residuals. Compare actual joint tensors before taking class quotients.
3. **Operational coherence:** determine which larger boundary diagrams
   commute with persistent-port transport and which retain holonomy.
   The local pentagon is an exact base example, not a proof of coherence
   or of the global spectrum.

These tests turn the user's fivefold prompt into constructions, bounds
and operation comparisons rather than a restriction on what may be
explored.
