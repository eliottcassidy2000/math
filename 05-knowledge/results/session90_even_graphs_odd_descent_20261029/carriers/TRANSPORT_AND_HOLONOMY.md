# Square-root transport, centralizer memory, and a nine-vertex obstruction

**Status:** the algebraic statements below are PROVED directly. The minimum
odd tournament order and numerical profiles are FINITE-EXACT with the stated
complete small universe. The Bruck-loop construction is CITED classical theory;
its concrete occurrence here is independently checked. No external novelty
claim is made. This is a new application relative to the scoped sources inspected,
not a claim to have searched the complete repository or literature.

## Inheritance and concept board

The starting point is Session89's
`encoding/ODD_SELF_CONVERSE_LIFT.md` in
`05-knowledge/results/session89_consecutive_tournament_pairs_20261028/`:
aligned pairs of tournaments on n and n+1 vertices correspond to self-converse
tournament classes on 2n+1 vertices. The retained group is the centralizer of
a chosen anti-involution, not the whole automorphism group. The smallest
hostile is the directed triangle: its pair has group1 but its lift has group3.

The incoming `singleton_pair_reconstruction_20261028/CONNECTIONS.md` credits
Eplett's 1979 odd-order representation. The lifted representation and count
identity are inherited, not discoveries here. The new questions concern
explicit transport, its composition, and exact retained symmetry.

Live concepts: odd automorphism groups; anti-involution markings; marked versus
unmarked mass; local versus based transport; twisted subgroups; named ports.
The least-used sidecar in this lane is the *actual conjugating permutation*,
rather than its existence or the centralizer's order.

## 1. A supplied isomorphism can be repaired without searching a group

Permutations compose from right to left. Let W be a self-converse tournament,
A=Aut(W), and let s,t be anti-involutions. A has odd order: an involutory
automorphism would exchange the ends of an arc. Therefore x=ts lies in A and
has odd order. Define

    g(t,s) = sqrt(ts),

where sqrt denotes the unique odd-order permutation square root. Then

    g(t,s) s g(t,s)^(-1) = t,
    g(s,t) = g(t,s)^(-1).

**Proof.** sxs=x^(-1). If x has odd order m, its root is
x^((m+1)/2), so s g s=g^(-1). Hence g s g^(-1)=g²s=xs=t.
Uniqueness of the root follows already in an odd group: a square root commutes
with its square, and an odd-order element is a power of its square. The same
root is unique among odd-order permutations of the ambient vertex set.

For an explicit permutation array, each cycle has odd length l. Replace its
one-step rotation with the rotation by (l+1)/2. This takes O(N) array operations
and O(N) storage on N vertices, including linear validation of the permutation.
There is no need to calculate a potentially huge least-common-multiple order.
The arithmetic on vertex indexes has O(log N)-bit operands. Checking that the
supplied maps really are tournament (anti-)isomorphisms still costs O(N²) in
explicit adjacency storage.

For an unmarked isomorphism phi:W→W' and supplied marks s on W and t on W', put

    s' = phi s phi^(-1),
    corrected_phi = g(t,s') phi.

The result is an isomorphism intertwining the chosen marks. The Session89
decoder then turns it into an aligned-pair isomorphism. Thus **the search for
an automorphism correcting phi is unnecessary** once phi and both marks are
given. This does not find phi, does not find an anti-map in an arbitrary input,
and does not solve tournament isomorphism in polynomial time.

The operation is natural under relabeling:

    g(phi t phi^(-1), phi s phi^(-1)) = phi g(t,s) phi^(-1).

That follows from uniqueness of the odd square root.

## 2. Exact factorization of the symmetry group

Fix a base anti-involution t0 and set

    C = C_A(t0),
    I = {anti-involutions of W},
    K = {x in A : t0 x t0 = x^(-1)} = {s t0 : s in I}.

K contains the identity, inverses, powers, and xyx whenever x,y lie in K.
It is a uniquely 2-divisible *twisted subgroup*, but need not be a subgroup.

Every a in A has a unique factorization

    a = k c,       k in K, c in C.

Indeed set s=a t0 a^(-1) and k=g(s,t0). Then k is in K and k^(-1)a centralizes
t0. Conversely k determines s=k²t0. Since squaring permutes K, distinct k
correspond to distinct anti-involutions, so the factorization is unique.
Consequently

    A = disjoint union over s in I of g(s,t0) C,
    |A| = |I| |C|.

This is an explicit coset transversal, not just a Sylow existence statement.
If all anti-involutions and the full marked group C have already been found,
this formula lists all unmarked automorphisms with no further isomorphism
classification. Enumerating those inputs and outputs can itself be expensive.
The abstract numbers |I| and |C| do not specify multiplication or port action.

## 3. Pair-local square roots can accumulate centralizer error

For three anti-involutions define the triangular composition discrepancy

    h(u,t,s) = g(u,s)^(-1) g(u,t) g(t,s).

Both the direct and two-step maps send s to u, so h lies in C_A(s). Its value
can be nonidentity even though every intermediate map is an exact transport.
This is a concrete reason that the endpoint anti-involution does not record
the entire accumulated named-port permutation.

**Flatness theorem.** The following are equivalent:

1. Every h(u,t,s) is identity.
2. The set K is pairwise commuting.
3. K is an abelian subgroup of A.

When these hold, A=K⋊C, with t0 acting by inversion on K and trivially on C.

**Proof.** Write s=t0, t=x²t0, u=y²t0 with x,y in K. Flatness is

    sqrt(y² x^(-2)) = y x^(-1).

Squaring and cancelling gives xy=yx. Conversely if K commutes, the displayed
identity holds and gives all triangles, using g(t,s)=b(t)b(s)^(-1) with
b(t)=g(t,t0). A commuting K is a subgroup because t0(xy)t0=x^(-1)y^(-1)
=(xy)^(-1). C normalizes K, K∩C=1, and the preceding factorization gives the
semidirect product. The converse implication from an abelian subgroup is
immediate. QED.

The condition is specific to this *natural pair-local square-root rule*.
It is not an obstruction to all possible flat transport conventions.

## 4. A global base mark supplies an exact flat alternative

Retain t0 and let b(s)=g(s,t0). Define

    f(t,s)=b(t)b(s)^(-1).

Then f(u,t)f(t,s)=f(u,s) identically, whether or not K commutes. The difference
between the natural pair-local map and this based map is exactly

    c(t,s)=b(t)^(-1) g(t,s) b(s) in C,
    g(t,s)=b(t)c(t,s)b(s)^(-1).

Thus two safe computational conventions are available:

- keep the reference anti-involution and use flat based transports;
- keep pair-local square-root transports, but retain their centralizer-valued
  corrections when comparing paths or shared ports.

The first convention uses a supplied global reference. It does not canonically
choose an anti-involution in an unmarked object. The second retains more than
an endpoint or an abstract automorphism order. Both are exact. This clarifies
the relation to LCA/alignment work: a common base removes artificial path
dependence; erasing it without retaining corrections can reintroduce ambiguity.
No LRC, Collatz, or arithmetic-conjecture implication is claimed.

## 5. The first odd-order nonflat control is nine vertices

Let W=C3[C3,C3,C3], with vertices (i,j) in F3², encoded as 3i+j.
The arc from (i,j) to (k,l) points forward around the macro C3 if i≠k,
and forward around the micro C3 if i=k. Its full automorphism group is
C3³⋊C3, of order81. The anti-involution (i,j)↦(-i,-j) has centralizer order3
and27 conjugate anti-involutions.

`check_transport.py` independently searches all of S9 to verify that this
explicit wreath group is the entire automorphism group, not a subgroup.
It checks every one of the27³ ordered mark triples. Of these,
**11664 have nonidentity holonomy** and8019 have identity holonomy.

One witness, in permutation-image notation on0,...,8, is

    s = (0,2,1,6,8,7,3,5,4)
    t = (0,2,1,7,6,8,4,3,5)
    u = (3,5,4,0,2,1,6,8,7)
    h = (0,1,2,4,5,3,8,6,7).

The error fixes the anti-involution s but moves named port3 to4. Thus it is
invisible to an unmarked tournament-class observation while remaining visible
to a later operation that refers to these named ports.

All99 odd self-converse classes of orders1,3,5,7 (counts1,2,8,88) are checked,
using the inherited complete class bijection and full permutation searches for
their automorphism groups. All their square-root triangles are flat. Therefore
**9 is the minimum odd tournament order with a nonflat square-root transport**.
This minimum is FINITE-EXACT, supported by the complete lower-order universe
and the explicit order9 witness; it is not a claim about even orders.

Paley7 has |A|=21, |I|=7, |C|=3 and is flat. The nine-vertex example shows why
nontrivial centralizer alone is not enough to cause holonomy: K must fail to
commute. Here the familiar numbers7 and21 describe an exact orbit-centralizer
factorization, whereas27 and81 furnish the first odd-order nonflat control.
The separate twelve-vertex81/27 equal-cyclic-triples pair remains a distinct
example; sharing group orders is not an isomorphism between the two problems.

## 6. A classical nonassociative structure, not a numerical analogy

For K above, define

    x circle y = sqrt(x y² x).

The classical twisted-subgroup construction makes this a left Bruck loop.
See Petr Vojtechovsky, [Three lectures on automorphic loops, Proposition2.3](https://www.cs.du.edu/~petr/data/papers/three_lectures_on_automorphic_loops.pdf),
which gives the construction for uniquely2-divisible twisted subgroups.
We import that theorem, not a new classification of loops.

In the nine-vertex control K has27 elements. The checker verifies closure,
identity, all left and right Latin conditions, the automorphic inverse law,
and every left Bol identity directly. It also finds an explicit failure of
associativity. Therefore this tournament's changing anti-involution frames
carry a concrete **nonassociative Bruck loop of order27**.

This is an exact operation-level bridge: the source is the anti-involution
frame set, the coordinate map is s↦s t0, and the operation is the displayed
square-root formula. Changing the chosen t0 changes the coordinate model.
The entire unmarked tournament is not determined by this loop, and no claim
of loop associativity being equivalent to transport flatness is made. Flatness
has the stronger, explicitly proved ambient-commutativity criterion above.

## 7. Verification and output scope

Run from this directory:

    python3 check_transport.py
    python3 -O check_transport.py

Both return the same report with **73467 active checks**. The script is
standard-library-only and does not use Python assert statements. Its universe
is all natural pair classes n≤3, all their triangle transports, all arbitrary
vertex relabelings and corrected marked isomorphisms through lift order5,
and the complete wreath9 anti-involution set. The lower-order class counts
are not obtained from an OEIS lookup. Full group searches range over S_N
with only a logically necessary degree filter.

The generated JSON is a local audit receipt, not a required checked-in input.
The executable and this note suffice to reproduce all stated finite results.
