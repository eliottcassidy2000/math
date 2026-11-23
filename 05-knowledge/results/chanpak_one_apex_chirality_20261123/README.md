# The exact threshold for seeing converse orientation

This extends the [continuation observability packet](../chanpak_continuation_observability_20261123/README.md).
All old labels are fixed. An exterior vertex is **uniform** if it
beats the entire old core or loses to the entire core, regardless
of its other exterior contacts. Other exterior vertices can
individualize the core. There is a sharp observer boundary:

- With **one** individualized vertex, even arbitrary uniform hosts
  can be blind to orientation. There are exponentially large
  families of pairwise nonisomorphic strong examples.
- **Two** individualized vertices, with finite extra uniform
  probes, recover every labelled arc of an arbitrary tournament.
  See [TWO-APEX.md](TWO-APEX.md).
- The explicit 21-vertex pair is already distinguished by two
  singleton contacts **with no additional wrap**, by exact ordinary
  Hamilton-path counts.

The zero/one/two statement concerns permitted contacts on the
old core, not how many outside vertices are present.

## Exact one-apex invisibility criterion

Let \(q_m(T)\) count ordered covers of \(T\) by \(m\) directed
paths allowing empty slots. In such a cover distinguish slot 1
and let \(D_m(v)\) be the number in which its path starts at
\(v\), minus the number in which it ends at \(v\). For a vertex
\(x_S\) beating precisely \(S\) of the old core, on fixed labels,
\[
 q_m(T+x_S)-q_m(T^{\rm op}+x_S)
   = m\sum_{v\in S}\bigl(D_m(v)+D_{m+1}(v)\bigr).
\]
To prove this, fix the slot containing \(x_S\); there are \(m\)
identical choices. If \(x_S\) stands alone, there is no change.
If it is an endpoint of a nontrivial path, the two possible
positions contribute \(2\sum_SD_m\): swapping starts with ends
in the core reverses the sign and \(\sum_V D_m=0\).
If \(x_S\) is internal, separate its two nonempty pieces into
two named slots. Let \(K(u,v)\) count covers with the first
piece ending at \(u\) and the other starting at \(v\).
Reversing the core and exchanging the two pieces transposes
this matrix. The difference across the cut from \(S^c\) to
\(S\) is the sum on \(S\) of column totals minus row totals.
Those totals give \(D_{m+1}-D_m\): the distinguished path is
nonempty and the other named slot is now required nonempty, so
subtract the case where that other slot was empty. Combine the
endpoint and internal contributions.

**Necessity.** If all one-individualized hosts agree, apply
[uniform cyclic probes](../chanpak_continuation_observability_20261123/INTERFACE.md)
to the augmented core \(T+x_S\). These recover its whole slot
vector and are uniform towards the original \(T\). For singleton
\(S=\{v\}\), the formula forces \(D_m(v)+D_{m+1}(v)=0\).
Each \(D_m(v)\) is a polynomial in \(m-1\) of degree at most
\(n-2\), for \(n=|T|\ge2\): fix the vertices of the distinguished
slot; a singleton contributes zero, and its complement of size
at most \(n-2\) has a cover polynomial of degree at most its
size. Hence it suffices to test \(m=1,\ldots,n-1\). They
annihilate the polynomial \(D_m+D_{m+1}\), and the leading
coefficient shows that this implies \(D_m=0\). The order-one
case is immediate.

**Sufficiency for arbitrary hosts.** The coefficients of
\(\binom{m-1}{r-1}\) in \(D_m(v)\) are the corresponding
differences for \(r\) ordered nonempty runs; this follows by
choosing which other slots are occupied. They all vanish
when \(D_m=0\). Now fix the outside itinerary and precisely
which gaps between its vertices contain a nonempty old run.
Uniform exterior contacts do not inspect the old run
endpoints. If the sole individualized apex touches no run,
reversing every run preserves the count; if it touches one
run, its test is exactly the vanishing nonempty-run endpoint
difference. If it touches two runs, the change is the cut of
the skew part of the same kernel \(K\); its divergence is
the vanished single-slot difference for \(r\) nonempty runs.
Summing over outside itineraries proves sufficiency. This is
an **if and only if for the entire observer class**, not just
for the unpadded one-apex cut.

Every vertex-transitive tournament has \(D_m(v)=0\):
by symmetry all vertices have the same start count and the
same end count, and their totals agree. Transitivity alone
does not imply isomorphism with the converse.

## A chiral core invisible to one apex

Use vertices \(a+7b\), with \(a\in\mathbb Z/7\), \(b\in\mathbb Z/3\),
and the semidirect product
\((a,b)(c,d)=(a+2^b c,b+d)\).
The steps
\[
 \{1,2,3,14,8,9,10,11,12,13\}
\]
choose one element in each inverse pair. Orient \(g\to h\)
when \(g^{-1}h\) is in this step set. Left multiplication
acts transitively, and steps 1 and 14 generate all vertices
as directed products, so the tournament is strong.

It is not isomorphic to its converse. By transitivity any
such isomorphism can be translated to fix the identity.
In either rooted out-neighborhood the score multiset is
\(3,3,4,4,4,4,5,5,6,7\). But the unique degree-seven
vertex has out-neighbor score multisets respectively
\(3,4,4,4,5,5,6\) and \(3,4,4,4,4,5,6\).
This small local certificate avoids counting any
21-vertex Hamilton paths to establish chirality.

### Exponential ambiguity survives restriction to strong cores

Independently use \(T\) or \(T^{\rm op}\) in each of \(k\) ordinal
blocks, and wrap the ordered join \(B\) in the strong cycle of
blocks \(a\to B\to b\to a\), with singleton anchors \(a,b\).
Every vertex outside a block is uniform toward that block,
except possibly the designated apex testing the entire
core. Replacing the blocks one at a time thus preserves
the whole response under fixed labels, by the criterion.
The resulting \(2^k\) cores have order \(21k+2\), are strong,
and are pairwise nonisomorphic. The anchors are identified by
their unique extreme degrees \(21k\) and \(1\): each old vertex
has ten in- and out-neighbors in its own block. Delete the
anchors and the positions of the strong components are
intrinsic in the ordinal order. Each component's chirality
then recovers its corresponding bit. Equality of responses
is genuine observation loss, not a mystery caused by sorting
two endpoint heights.

## A simple second-apex witness, with a general regular-core read

In the chiral example, let \(x\) beat old vertex 0 only,
\(y\) beat old vertex 1 only, and \(x\to y\). With all other
contacts into the two apices, the absolute heights are
\[
 H(T+x+y)=1134378499638097,\qquad
 H(T^{\rm op}+x+y)=1136398069453171.
\]
This is an exact count at order 23, not just a difference in
a hypothetical upper coefficient. The ordinary subset DP
replay gives these and the two exchanged-contact counts.

A smaller analytic check is useful in its own right. For a
regular tournament of degree \(d\) with adjacency matrix
\(A\), put singleton contacts at distinct old \(u,v\),
respectively, keeping \(x\to y\). The number of directed
5-cycles containing both new vertices is
\[
 d^2+2d-1-(A^2)_{vu}.
\]
Indeed cycles using \(x\to y\) follow a two-arc old path
starting at \(v\), excluding those ending at \(u\). Otherwise
the two possible placements of \(y\) contribute
\(d-[v\to u]\) and \(d-[u\to v]\).
Swapping the two singleton contacts changes this 5-cycle
count by
\[
 (A^2)_{uv}-(A^2)_{vu}=A_{vu}-A_{uv}.
\]
The identity follows by squaring \(A^t=J-I-A\) and
using \(AJ=JA=dJ\). In the Cayley example the mixed counts
are 116 and 117.

In the ordered cover polynomial on \(N\) vertices,
the coefficient of \(m^{N-4}\) is twice the 5-cycle count
plus four times the number of pairs of vertex-disjoint
triangles, by [odd-cycle packing](https://arxiv.org/abs/2307.05569).
Swapping these singleton contacts preserves the latter number:
terms meeting just one new vertex cancel in pairs, terms
with one new vertex in each triangle are symmetric in \(u,v\),
and terms with both apices in one triangle count old
triangles avoiding \(v\) or \(u\). In a regular core the number
of old triangles through a vertex is \(d(d+1)/2\): count arcs
from its out-neighbors to its in-neighbors using the outdegree
sum of the former. Thus this coefficient contrast equals
\(2(A_{vu}-A_{uv})\), another exact two-contact test
when the old core is regular. The [general theorem](TWO-APEX.md)
does not need regularity.

~~~
python3 direction-replay.py
g++ -O3 -std=c++17 -Wall -Wextra -Werror direction-height.cpp -o direction-height
./direction-height
~~~
For the order-23 DP, pair-swapping consecutive positions in a
vertex permutation bounds a Hamilton count by
\(23!/2^{11}<2^{64}\). The same bound at each smaller
subset is no larger, so unsigned 64-bit cells are exact.

[Minimality and regular collision families](MINIMALITY.md): order 21 is the first vertex-transitive chiral example; the prior through-eight census and a new regular order-nine census also exclude smaller examples through ten in the full one-individualized converse class. Nontransitive odd orders 11–19 are not settled. Carousel substitution gives regular strong examples counted by binary necklaces.
