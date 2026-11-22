# The geometry of the two absent Hamiltonian counts

Let \(X_n(L)\) be the subcomplex of the \(\binom n2\)-dimensional
**labelled arc-flip cube** consisting of the faces all of whose
orientations have at most \(L\) Hamiltonian paths. A free edge is one
coordinate; whether freeing it stays below \(L\) depends on the other
arcs, not on the unlabelled two-vertex tournament. Quotienting by
relabeling is a different space; no quotient homotopy is claimed.

**PROVED with finite local certificates.** The 7-plateau has a
canonical permutohedral model for every \(n\), with \(X_5(7)\simeq S^3\).
The whole interval \(15\le L\le23\), not just the 21-plateau, has a
common permutohedral homotopy model for every \(n\); natural inclusions
within the interval are homotopy equivalences. At \(n=6\) the model
is \(S^4\). The [nearest literal one-arc crossings](LITERAL-BRIDGES.md)
distinguish the two holes even though the latter is topologically
invisible in this filtration up to homotopy.

## The union-of-completions lemma

For a face \(F\), give fixed pairs their fixed direction and free
pairs both directions. The resulting union digraph has strongly
connected components in an ordered partition \(\pi(F)\), because
its condensation is an acyclic tournament. All intercomponent arcs
are fixed, and the maximum height on \(F\) is the product of the
maxima on the induced partial faces on its blocks.

**Lemma.** A mixed complete orientation on at least three vertices
with strong union has a strong completion. Choose a completion
maximizing the sum of squared SCC sizes. If it is not strong, order
the components. A free arc from component \(i\) to \(j>i\) can be
reversed and merges components if a forward path remains between
its endpoints without using that arc. This is possible if
\(j>i+1\) by an intermediate vertex, or if either endpoint is in
a nonsingleton component by strong connectivity there. No old
component splits. Thus in a nonstrong maximizing completion the
only free cross arcs join consecutive singleton components.
Strong union needs a free arc crossing each cut between consecutive
components, so all blocks must be singletons and all consecutive
arcs free. Reversing two consecutive arcs forms a strong triangle
(the skip arc stays forward), contradicting maximality.

A two-vertex union block is exceptional: its pair is free but
neither orientation is strong; its height is always 1. Denote
by \(C_{\pi(F)}\) the full face with the same block order and all
internal arcs free. If \(F\subseteq G\), the union SCCs of \(F\)
refine those of \(G\), in the same order across distinct blocks;
in particular \(C_{\pi(F)}\subseteq C_{\pi(G)}\).

## The 7-plateau: canonical closure at every order

The strong minima at sizes 3 through 7 are \(3,5,9,15,25\),
and larger orders do not lower the bound. Use [Busch's result](https://doi.org/10.37236/1141)
and the inherited [all-order proof, THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The maxima over all tournaments of sizes 3 and 4 are 3 and 5
(equivalently at most two directed triangles at order 4, with
\(H=1+2c_3\) there). The union lemma and product rule characterize
the retained faces:
- At \(L=1\) all blocks have size at most 2.
- At \(L=3\) additionally one size-3 block may occur.
- At \(L=5\) or \(L=7\) additionally one size-3 **or** size-4 block may occur.

For the last case, any size-5 or larger block admits a strong
completion of height at least 9, and two blocks of size at least 3
already give a joint completion of height at least \(3^2\).
Conversely every listed full block face has height at most 5.
In particular this gives \(X_n(7)=X_n(5)\) directly.

At each cutoff the enlargement \(F\mapsto C_{\pi(F)}\) is a
monotone, extensive, idempotent map on nonempty faces. The
pointwise inequality \(F\le C_{\pi(F)}\) gives the order-complex
homotopy from identity to closure, fixed on its image. The image
poset is exactly the face poset of the matching subcomplex of the
ordinary permutohedron (faces are ordered partitions; refinement is
the boundary relation). A permutohedron face with \(b\) blocks
has dimension \(n-b\), regardless of the cubical dimension.

The \(L=3\) model contains the entire 2-skeleton and is simply
connected; at \(n=4\) it is exactly the boundary \(S^2\).
The \(L=5\) model contains the entire 3-skeleton and is
2-connected; at \(n=5\) it is exactly the boundary \(S^3\).
Connectivity follows by cellular approximation in the full
contractible permutohedron: a sphere of dimension \(i\) and
its filling need only cells through dimension \(i+1\).

## The 21-plateau: the six-vertex decorations collapse

Fix \(15\le L\le23\). The union lemma excludes any block of size
at least 7. Let \(Y_n\subset X_n(L)\) contain the faces with **no
six-vertex union block**. The maximum over tournaments of size 5
is 15 (see the small replay). Thus \(Y_n\) allows arbitrary
singles and pairs, along with any of: up to two size-3 blocks;
one size-4 block with at most one size-3 block; or one size-5
block and no other size 3 or larger. Counts multiply, so this
list follows from the small extrema and is independent of \(L\).
The block-face map is again a closure. Its permutohedron model
contains the entire 4-skeleton, hence is 3-connected. For \(n=6\)
every proper ordered partition is allowed, so \(Y_6\simeq S^4\).

The finite certificate [holes-relative-collapse.cpp](holes-relative-collapse.cpp)
collapses \(X_6(L)\) relative to \(Y_6\). A base-three face digit
0 or 1 fixes an orientation; 2 frees that pair. The face maximum
is the maximum of the two subfaces obtained by fixing the first 2.
Vertex heights come from the subset path recurrence
\(D[S,v]=\sum_{u\to v,\ u\in S-v}D[S-v,u]\), with singleton seed 1.
Mark exactly the retained faces whose union is strong. Only
these are outside \(Y_6\), and any coface of a marked face is also
outside \(Y_6\) since adding arcs preserves strong connectivity.
If an outside face has exactly one immediate coface still present,
delete the pair. The upper face is maximal: otherwise a cube
two dimensions above the free face would provide two different
immediate cofaces. Thus each pair is an elementary collapse
relative to \(Y_6\). All outside cells are removed:

| \(L\) | Cells initially outside \(Y_6\) | Free-face pairs removed | Cells left outside \(Y_6\) |
|---:|---:|---:|---:|
|15|40,320|20,160|0|
|17|158,400|79,200|0|
|19|267,840|133,920|0|
|21|267,840|133,920|0|
|23|587,520|293,760|0|

Oddness covers all intermediate real cutoffs. The \(L=19\) and
\(L=21\) complexes are identical by the known gap; one can also
read this off the six-vertex counts and the product reduction.
There are more than twice as many relative cells at 23 as at
21, but they do not change the homotopy type.

**All-order gluing.** An outside face has a unique size-6 union
block \(S\): any other block of size at least 3 could combine
with a strong completion of \(S\) to reach \(15\cdot3>23\).
Let \(A,B\) be the sets of vertices before and after \(S\).
For each fixed labeled triple \(A|S|B\) the outside faces belong to
\[
 X_A(1)\times X_S(L)\times X_B(1)
\]
with interblock arcs fixed forward. Its intersection with \(Y_n\)
replaces \(X_S(L)\) by \(Y_S\). Products indexed by different
triples meet only inside \(Y_n\): a face outside \(Y_n\)
determines \(A,S,B\) uniquely. The local strong deformation
retraction fixes \(Y_S\); the product of that homotopy with
identities on the outside factors glues to a strong deformation
retraction \(X_n(L)\to Y_n\). This also proves that the inclusions
for increasing \(L\) within the range are homotopy equivalences:
all contain the same homotopy-equivalent subcomplex \(Y_n\).

## Euler series and replay

Let \(a=1-z+z^2/2,\ b=z^3/6,\ c=-z^4/24,\ e=z^5/120\).
A block of size \(j\) contributes signed exponential weight
\((-1)^{j-1}z^j/j!\) in the permutohedron model. A sequence of
singles and pairs gives \(1/a\). Therefore
\[
 \sum_{n\ge0}\chi(X_n(7))\frac{z^n}{n!}
 =\frac1a+\frac{b+c}{a^2},\qquad
 \sum_{n\ge0}\chi(X_n(21))\frac{z^n}{n!}
 =\frac1a+\frac{b+c+e}{a^2}+\frac{b^2+2bc}{a^3}.
\]
The mixed term has two possible block orders. At \(n=6\) the
values are \(-30\) and 2. The raw cubical counts through six,
without using the model, are in [holes-cubical.out](holes-cubical.out).

~~~sh
g++ -O3 -std=c++17 -Wall -Wextra -Werror holes-cubical.cpp -o holes-cubical
./holes-cubical
g++ -O3 -std=c++17 -Wall -Wextra -Werror holes-relative-collapse.cpp -o holes-relative-collapse
./holes-relative-collapse
~~~

## What survives from the planar suggestion?

The identities \(7=2\cdot3+1\), \(21=2\cdot11-1\) correctly put
both holes on the signed odd lattice, but do not specify a
contact, face, or invariant of the flip cube. The two arithmetic
sheets \((p,p)\), \((p,7p)\) can index possible count pairs; the
slope does not decide membership. For example \(49=7\cdot7\)
and \(63=7\cdot9\) do occur even though 7 does not. The real
link is **strong-component multiplication**: at height 21, the
only alternative to a single strong 21-core is a 3-core and a
missing 7-core. The strong 21-core requires a separate exclusion.
A count point on either line is not an arc flip.

The model suggests studying **marked** cells or correspondences
which transform the absolute height, retaining which paths are
gained and lost. Ordinary homotopy of the unmarked sublevel
complex is already too coarse to detect 21, since even 17, 19
and 23 lie in its homotopy range. This neither rules out a finer
geometric invariant nor proves that all other odd counts occur.

See also [sevenfold geometry](SEVENFOLD.md) for a no-go theorem
for the simplest cyclic wrapper and a compressed, weighted-cut
explanation of the strong 21-core exclusion.
