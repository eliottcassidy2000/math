# An equivariant sign-reversing involution inside each inflated tournament

2026-10-25. **PROVED:** the weighted clone-removal sum cancels separately
over each isomorphism class of inflated tournament. The involution preserves
the literal inflated graph and descends to weighted quotient isomorphism
classes. **FINITE-EXACT:** complete small fibres and symmetric hostiles.
This explains the signed substitution and its clone mergers; it does not
give a faster enumerator or a black-path transport theorem.

## 1. Corrected coefficient meaning and inherited structure

Use the corrected interpretation in
[the independent review](INDEPENDENT_REVIEW.md) and
[EXCESS_EXTENSION.md](EXCESS_EXTENSION.md):

    (-1)^E b_(M,E)
      = number of isomorphism classes of M-vertex tournaments
        with nonnegative integer vertex weights summing to E.       (I1)

Weights are preserved by isomorphisms. This is the integer meaning of the
bivariate coefficients, not just a rational cycle-index expression. In
Burnside language a vertex cycle of length d receives a common weight
and contributes1/(1-z^d). Tournament automorphisms have only odd cycles;
replacing z by -z gives1/(1+z^d), recovering the bivariate specialization.

Inflate a weighted vertex v to a transitive tournament on1+w(v) vertices.
Every cross-block arc has the quotient tournament's direction. The result
has N=M+E vertices. We group the signed objects (I1) by this inflated
tournament T.

The load-bearing inherited theorem is
[THM-4463, homogeneous tournament contacts and pinned scale floor, section 2](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md):
the graph J(T) of homogeneous pairs is a disjoint union of directed paths,
and its vertex components are exactly the maximal transitive modules of T.
A homogeneous pair means that all other vertices see its members alike.
The path orientations are the original tournament directions; the order
inside each component is intrinsic and has no nontrivial automorphism.

The live objects are the fixed inflated graph, its canonical module paths,
their literal gaps, the automorphism action, and the cut partition. The
canonical hostile is a threefold rotational symmetry: a labelled choice
of just one clone block does not descend to unlabelled objects.

## 2. The exact inflation fibre is a cut-set orbit space

Fix a labelled representative T of the inflated isomorphism class. Write
its canonical maximal transitive modules in their directed orders as

    L_1,...,L_q,
    G=the set of adjacent gaps in all these paths.

Every inflated block is a transitive module of T. By THM-4463 it lies in
one L_i, and it must be a contiguous interval there: omitting a vertex
between two members would let that vertex distinguish them. Conversely,
every contiguous interval in an L_i is a transitive module of T. Its
incidence is constant to other components and to all vertices preceding
or following the interval in its own component.

It follows that decompositions of the literal T into inflated blocks are
exactly the subsets C of G: cut at every selected gap. The resulting
interval blocks determine their weighted quotient uniquely, with weight
equal to interval length minus one. In particular,

    M=q+|C|,       E=N-q-|C|,
    sign(C)=(-1)^(N-q-|C|).                          (I2)

Let A=Aut(T), acting on gaps and cut subsets. The weighted isomorphism
classes in this inflation fibre are exactly the A-orbits of cut subsets.
One direction is immediate from an automorphism of T. Conversely, any
weight-preserving quotient isomorphism sends equal-sized interval blocks
to one another. It lifts uniquely to the order isomorphism on each pair
of transitive blocks; uniform cross-block arcs make the resulting map an
automorphism of T. Thus two quotient objects are isomorphic iff their
cut subsets lie in the same A-orbit. Choosing another identification of
an inflated object with T changes the cut subset only by such an orbit.

This step retains clone mergers exactly. Adjacent quotient vertices may
merge into the same maximal transitive module after inflation. They are
represented by different cuts of that fixed module, rather than being
assumed to remain canonical blocks.

## 3. The orbit-sized toggle

Suppose G is nonempty. We choose one A-orbit O of nontrivial maximal
modules, depending only on the uncoloured inflated T, not on C. The
choice can be made invariant under all isomorphisms: encode the marked
object (T,L_i), with L_i distinguished, by the least adjacency-and-marker
code over all vertex labellings, and choose the least code among nontrivial
modules. Two modules have the same code precisely when an automorphism
of T sends one to the other. This selects exactly one orbit O. It is a
finite canonical choice, not an efficient-canonicalization claim.

Every finite tournament automorphism group has odd order. An involution
would exchange a vertex pair and reverse its mutual arc, so cannot be an
automorphism. In any finite group of even order, pairing elements with
their inverses leaves a nonidentity involution. Hence |A| is odd, and
every A-orbit of modules, including O, has odd cardinality.

For every L in O take its first directed gap, between its first and second
vertices. Let F be the set of those gaps. Automorphisms preserve intrinsic
orders and permute the modules of O, so F is A-invariant and

    |F|=|O| is odd.

Define

    iota(C)=C symmetric-difference F.                (I3)

This toggles the first gap simultaneously in every selected module. It
leaves every vertex and arc of T unchanged. On a currently cut first gap
it merges the initial singleton interval with the following interval;
on an uncut gap it splits the initial vertex from its interval. All
resulting block sizes remain positive.

The map is an involution because F depends only on T. It is equivariant
because F is invariant:

    iota(g C)=g iota(C),       g in A.

Moreover

    |iota(C)|-|C|=|F|-2|C intersect F| is odd.        (I4)

Equations (I2) and (I4) reverse the sign. The operation therefore descends
to a sign-reversing involution on weighted quotient isomorphism classes
in the fibre. No cut orbit can be fixed: automorphisms preserve |C|,
whereas the paired cut subset has the opposite cardinality parity.

The requirement that O be selected from T alone is essential. Selecting
it from the changing weighted quotient or its cut pattern would not
automatically produce the same F on the return step.

## 4. Survivors and the identity

If T has a homogeneous pair, G is nonempty and every fibre object is
paired by (I3). Its signed contribution is zero. If T has no homogeneous
pair, every canonical module is a singleton, G is empty, and the sole
decomposition is T itself with all weights zero. It contributes+1.
The empty tournament is also this latter case, with the empty partition.
Thus, separately for every inflated isomorphism class,

    sum_(weighted quotients inflating to T) (-1)^E
       =1 if T is pair-free, otherwise0.            (I5)

Summing (I5) over all N-vertex tournaments gives

    sum_(M+E=N) b_(M,E)
       = number of pair-free N-vertex tournaments.  (I6)

This is a direct combinatorial cancellation proof of the diagonal
clone-removal identity. It also explains why odd-order automorphisms
matter: an invariant orbit of gaps has odd size and hence reverses the
sign. The invariant cut set is the required sidecar. Merely counting
maximal modules or quotient vertices would not specify the operation.

## 5. Two hostile symmetry controls

For T=C3[TT2,TT2,TT2], there are three maximal modules, each of order two.
Its automorphism group C3 rotates them. Their first gaps form one orbit
F of size three. The cut-subset orbits have cut counts0,1,2,3, one orbit
at each count. Their signs from (I2) are respectively minus, plus, minus,
plus. The valid toggle complements all three gaps, pairing0 with3 and1
with2.

Toggling just a labelled first gap g0 is not well-defined on isomorphism
classes. The cuts {g0} and {g1} lie in the same orbit, but toggling g0
sends them to the empty cut set and to {g0,g1}, which have different
cardinalities and cannot lie in the same orbit. This refutes the tempting
single-block rule while preserving the orbit-sized repair.

At T=C3[TT2,TT2,1], there are two nontrivial modules and the automorphism
group is trivial. Toggling the first gaps of *all* nontrivial modules
would toggle two gaps and preserve the sign. Choosing one automorphism
orbit toggles one gap and succeeds. Thus neither a single labelled
module nor the union of every nontrivial module is the general rule.

## 6. Exact finite controls and scope

Run `python check_clone_involution.py` and
`python -O check_clone_involution.py`. Both passed. The checker uses
literal tournament arcs, all vertex permutations, weighted canonical
codes, and gap subsets; it imports no cycle-index or affine-state engine.

- It enumerates every tournament isomorphism class for orders0..5,
  with class counts1,1,1,2,4,12. Across their fibres it checks58 cut-set
  orbits, including empty and pair-free survivors.
- For every fibre it independently compares weighted quotient canonical
  codes with Aut(T)-orbits of cuts, checks the intrinsic module paths,
  uniform block incidences, the exact block count (I2), canonical orbit
  selection, odd toggle size, equivariance, no fixed orbit, and (I5).
- It checks the complete eight-subset fibre of C3[TT2,TT2,TT2], its four
  cut orbits, and the single-labelled-gap failure. The asymmetric
  C3[TT2,TT2,1] control checks why selecting all modules can fail.
- Every gate is an explicit exception and remains active under `-O`.

The proof preserves the inflated tournament even when the source quotient
gains or loses vertices and clone blocks merge. Computing the canonical
module orbit can itself require substantial isomorphism work. No smaller
unrestricted counting state, complexity bound, port transport, or
black-core connectivity result follows merely from this cancellation.
