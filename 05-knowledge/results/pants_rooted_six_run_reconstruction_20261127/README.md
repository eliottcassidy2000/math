# Recover a skew conference tournament from one rooted six-run slice

2026-11-27. The [six-run calculation](../pants_skew_conference_quartic_detection_20261127/README.md)
showed that the first nonconstant switch observable for a skew
conference tournament is quartic. It is also enough to recover
every arrow up to reversing all of them. Here is an elementary,
adaptive quadratic-query reconstruction from signed coefficients
through any one fixed root. The only ambiguity, global converse,
is unavoidable.

Let \(A\) be a labelled skew adjacency with \(A^2=-(n-1)I\), \(n\ge8\).
Write \(n=4t+4\). The quartic coefficient at a four-set \(S\) is
\(-3(n-6)!Q_S/4\), where \(Q_S=\sum_x\prod_{i\in S}A_{ix}\).
We take a query to be one exact signed value \(Q_S\), not a count
of direct evaluations of the path function. If the latter function
is given on the labelled switch group, Fourier inversion supplies
these coefficients. An unlabelled histogram is not enough for this
algorithm.

For any chosen root \(0\), the algorithm uses only values
\(Q_{\{0,a,b,c\}}\). In fact, the number of distinct queried
four-sets is at most
\[
\binom{n-2}{2}+2\binom{(n-2)/2}{2}
 +\begin{cases}2,&n\equiv0\pmod8,\\1,&n\equiv4\pmod8.\end{cases}
\]
The residues alone recover the switching class up to converse;
the signs then give the switch bits.

## 1. Residues and a double anchor give an elementary core reconstruction

Switch to make the root a source and remove it. Orthogonality to
its row makes the core regular of outdegree \(2t+1\).
Two core rows agree at \(2t\) other core vertices: their product
at the root is \(+1\), so the sum of products outside their own
labels and the root is \(-1\). Their common out- and in-neighbor
counts agree by regularity and are both \(t\). This is the doubly
regular property. Each oriented core edge lies in exactly \(t+1\)
cyclic triples, counted by taking the out-neighbors of its target
that are not out-neighbors of its source.

For a core triple, the rooted \(Q\) has residue \(n-4\pmod8\) if it
is cyclic, and \(n\pmod8\) if it is transitive. This is the earlier
four-set Pfaffian test with one source. A switch can only negate
\(Q\), which preserves both residues. Thus a query provides the
cyclicity bit in the normalized core, even before knowing its
normalizing switch.

Choose one more anchor \(a\) inside the core. Query every triple
of the form \(a,b,c\). Make the graph on the other core vertices
with edge \(bc\) if that triple is cyclic. It is bipartite between
the out-neighbors \(O\) and in-neighbors \(I\) of \(a\), since only
cross triples can cycle. For \(b\in O,\ c\in I\), the graph edge
\(bc\) means exactly \(b\to c\). Both parts have size \(2t+1\) and
the graph has degree \(t+1\) on both sides, by cyclic completion
counts. Any two vertices in the same part have a common neighbor
because \(2(t+1)>2t+1\). The graph is connected, so its unique
bipartition determines \(O,I\) up to exchange. Pick one ordering;
that choice fixes all arrows involving \(a\) and every crossing
arrow. Exchanging the parts reverses all of these arrows.

For each same-side pair \(b,c\), pick a vertex \(w\) on the other
side whose adjacencies to \(b,c\) are opposite. Such a vertex
exists: for a pair in \(O\), each has \(t+1\) out-neighbors in
\(I\), so their sets cannot coincide by the global common-out
bound \(t\). For a pair in \(I\), use the same argument with
in-neighbors in \(O\). The known two arcs give a two-edge path
through \(w\); one queried cyclicity bit for \(\{b,c,w\}\) then
fixes the missing \(bc\) arc. This proves the quadratic query
bound for the core and gives a constructive, doubly-regular
special case of the general
[Boussaïri–Ille–Lopez–Thomassé cyclic-triple reconstruction theorem](https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf),
without invoking that theorem. The useful feature is that the
second-anchor triplets determine where to look, rather than
requiring every rooted magnitude.

## 2. Signs determine the unique gauge

Reattach root \(0\) to the reconstructed core as a source and
call its matrix \(B\). Under a switch by signs \(\varepsilon\),
\(Q_S\) is multiplied by \(\prod_{i\in S}\varepsilon_i\);
converse leaves \(Q_S\) unchanged. Fix \(\varepsilon_0=1\).
For every nonzero reference coefficient the ratio between the
observed \(Q_S\) and \(Q_S(B)\) gives one parity equation on
the \(n-1\) unknown switch bits. Choosing the other possible
core orientation changes the final result by global converse:
take the converse of \(B\), then switch all nonroot vertices
to make its root a source again.

The residue alone guarantees nonzero coefficients on all
cyclic triples when \(n\equiv0\pmod8\), and on all transitive
triples when \(n\equiv4\pmod8\). These families are 2-designs
on \(v=n-1\) points. Each pair is in \(\lambda=t+1\) cyclic
triples and in \(v-2-\lambda=3t\) transitive ones.

If a nonempty set \(U\) of \(r\) vertices met every mandatory
triple evenly, every such triple with a pair in \(U\) would
contain exactly two vertices of \(U\). Count cross pairs:
\[
 r(v-r)\lambda=2\binom r2\lambda,\qquad r=n/2.
\]
Use the relevant positive design parameter, cyclic or transitive.
For the cyclic design, \(U\) would be a transitive tournament,
whose first two vertices have at least \(n/2-2=2t>t\)
common out-neighbors. For the transitive design, \(U\) would
have no transitive triple, impossible on four or more vertices:
some vertex in any four-vertex tournament has two out- or
two in-neighbors. Thus the mandatory rows span all \(n-1\)
coordinates over \(\mathbb F_2\).

**Only one or two further sign queries are needed.** Many sign
equations are already in the queried star through \(a\). If
the mandatory triples are cyclic, their star graph is the
connected bipartite graph above. The kernel consists exactly
of a constant bit \(u\) on \(O\), a constant bit \(w\) on \(I\),
and bit \(u+w\) on \(a\), so has dimension two. If the mandatory
triples are transitive, their star graph is the complement
on \(O\cup I\): it contains both cliques, and every vertex
has \(t>0\) cross neighbors. It is connected and has an odd
cycle. The equations \(x_b+x_c=x_a\) over its edges force
\(x_a=0\) and a single common bit on \(O\cup I\), a
one-dimensional kernel. In either case use the now-known
core to choose just enough independent additional mandatory
rows. Some might already have been queried in the same-side
step. Solve the resulting binary equations and switch \(B\)
back. Every resulting arrow is certified up to the common
converse ambiguity.

## 3. Replay

The companion Node.js replay starts from separately switched
bordered Paley examples of orders 8, 12, 20, 24, 32 and 44.
For each it checks the conference identity, reconstructs a
core using only the bounded queries, solves the signed binary
system, and compares the returned matrix to the original and
its converse. It also checks the predicted star ranks before
using the last one or two queries. The proof does not require
the Paley construction or a classification of conference
tournaments.
