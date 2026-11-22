# The primitive grammar for literal Hamilton-count seams

**PROVED.** An ordered pair of tournaments on the same vertex set,
differing in one arc and having counts \(N,N+2\), decomposes uniquely
(up to simultaneous isomorphism) as
\[
              \mathrm{TT}_i\ \Rightarrow\ (D,D')\ \Rightarrow\ \mathrm{TT}_j .
\]
The two transitive ends can be empty. The middle pair has no common
proper ordinal cut. If \(N\ge5\), **both** middle tournaments are
strong. There are just two low-count exceptions: the primitive
\(1\to3\) edge on three vertices, and two primitive \(3\to5\) edges
on four vertices distinguished by the orientation of the singleton
component of the lower tournament. An arc reversal that changes
count cannot have two weak middle tournaments.

This is a grammar for literal physical switches, not for pairs made
by choosing representatives of two attainable scalar heights. It
explains why attaching a scalar count witness to another count
witness cannot manufacture the perfect-number seam.

## Common cuts and the gcd constraint

An ordered cut is a nonempty proper subset \(A\) all of whose
vertices beat its complement. These cuts are exactly the
initial segments of the ordinal decomposition into strongly
connected components, so they form a chain. Let \(u\to v\) be
the changed arc from \(T\) to \(T'\). No common cut can separate
\(u,v\). Splitting at *every* common ordered cut therefore leaves
one middle interval with both endpoints. All other intervals
are unchanged strongly connected components; otherwise one
could add a common cut inside them. Conversely, all such
splittings are invariant under simultaneous isomorphism.

A Hamilton path is forced to visit ordinal components in order.
Consequently the product of the unchanged component counts
divides both \(H(T)\) and \(H(T')\). For a unit seam, these are
odd and differ by two, so the product is one. Every strong
tournament with more than one vertex has a directed Hamilton
cycle by Camion, and its rotations give more than one path.
Thus all unchanged components are singletons.

There is no possible nonzero middle switch with both sides
weak. In fact a proper cut \(A\) of the first middle tournament
and a proper cut \(B\) of the other must separate the changed
arc, or be common. Orient names so \(u\in A\setminus B\) and
\(v\in B\setminus A\). If \(a\in A\setminus B\) and
\(b\in B\setminus A\), the two cuts orient \(ab\) oppositely,
so necessarily \(a=u,b=v\). Now \(A\cap B\) and \(A\cup B\)
are common ordered cuts whenever nontrivial: the arcs between
these sections and their complements follow one or other of
the two original cuts and never use \(uv\). With no common
cut their sizes force \(A\cap B=\varnothing\), \(A\cup B=V\),
and hence \(V=\{u,v\}\). The counts would be equal.

## The two primitive weak exceptions

To turn a weak \(T\) into a strong \(T'\) with one reversal, the
arc must go from its first to its last strongly connected component
before reversal, or there is an untouched proper cut. The relevant
[inherited ordinal seam formula, extended in the previous packet](../chanpak_absolute_liouville_leech_20261122/MODULE-OFFSETS.md)
can be seen directly. For a split \(U\Rightarrow V\) across that
arc \(u\to v\), the paths lost are \(es\), where \(e\) counts full
paths of \(U\) ending at \(u\), and \(s\) counts full paths of \(V\)
starting at \(v\). The gained paths have portions
\[
      U_0,\quad V_1\ \text{ending at }v,\quad
      U_1\ \text{starting at }u,\quad V_0 .
\]
Here \(U_0,V_0\) can be empty; \(U_1,V_1\) are nonempty, each
splits its own vertex set with its partner, and all four pieces
are directed paths. Thus the gain is \(PQ\), the product of the
two pointed two-slot counts. Choosing the pointed piece just
the endpoint shows \(P\ge H(U-u)\ge e\) and \(Q\ge H(V-v)\ge s\).
The net jump is \(PQ-es\ge0\); in particular the weak member is
always lower if the difference is nonzero.

For at least three components there is nothing lost: no full
path can skip the middle components to use \(uv\). The two
endpoint components each contribute their pointed count. Every
interior component \(C\) contributes \(\sum_{S\subseteq C}H(S)H(C-S)
\ge 2^{|C|}\) by Rédei, multiplying the other factors. If it has
two or more interior vertices in total the jump is at least four.
If a nontrivial strong endpoint component exists, its factor
is at least three, by choosing the pointed piece of size one,
two, or the full component. The only remaining possibility
is three singleton components, giving \(\mathrm{TT}_3\to C_3\).

For exactly two components, the endpoint excesses \(P-e,Q-s\)
are even by applying Rédei to the same reversal against a singleton
on the opposite side. If both components are nontrivial strong
the excesses are positive, and \(e,s\ge1\) by rotating directed
Hamilton cycles. Therefore the jump is at least eight. If one
component is a singleton, say \(U\), the jump is \(Q-s\).
For a strong \(V\) of order at least four, take a directed
Hamilton cycle cut to a full path ending at \(v\), with suffixes
of sizes \(2,3,|V|\). These are three *different* non-singleton
contributions to \(Q\), on top of \(H(V-v)\ge s\). The even
jump is therefore at least four. The only possible \(V\) left
is a directed triangle, which gives \(Q-s=2\). Reflection yields
the other order of the two components. The triangle automorphisms
act transitively on its vertices, so each order yields one
simultaneous-isomorphism type. These two lower tournaments
are not isomorphic (one has a source and the other a sink).

## A primitive strong seam changes an arc active at both ends

Call an arc dead if no full directed Hamiltonian path uses it.
In a strong tournament, reversing a dead arc raises the count by
**at least four**. Use the [dead-arc separator theorem](https://digital.auraria.edu/files/pdf?fileid=de2a8eac-a8fb-44fa-b93e-e60575695962)
(Busch's dissertation, Theorem 4.12; see [Support and parity](SUPPORT-RIGIDITY.md)).
For a dead \(a\to b\) it gives a vertex \(z\) and the ordered
components \(C_1,\ldots,C_k\) of \(T-z\), with \(k\ge4\),
\(a\in C_1,b\in C_k\), and all \(C_2\to z\to C_{k-1}\).
Choose a full path in each component, starting at \(a\) in \(C_1\)
and ending at \(b\) in \(C_k\) (cut Hamilton cycles if nontrivial).
Write \(M=C_3,\ldots,C_{k-2}\) for the middle list, possibly empty.
After flipping to \(b\to a\), all four of the following are paths:
\[
\begin{array}{l}
z,\ C_{k-1},C_k,\ C_1,C_2,M;\\
M,C_{k-1},C_k,\ C_1,C_2,\ z;\\
C_2,\ z,\ C_{k-1},C_k,\ C_1,M;\\
M,C_k,\ C_1,C_2,\ z,\ C_{k-1}.
\end{array}
\]
The distinct positions (first, last, interior before the switched
arc, interior after it) distinguish them. All four use the reversed
arc at the \(C_k,C_1\) junction. No old paths are lost since the
forward arc was dead. In particular for a gap of two between strong
tournaments, the switched orientation is used on each side: it
cannot be dead on the lower side by this bound, nor on the upper
side by applying the same monotonic bound to the reverse operation.
This narrows the all-height search to genuinely active switches in
the stripped supports rather than arbitrary rewiring of missing arcs.

## The exact size generating function, and an arbitrary-gap extension

Let \(U_{a,b}(z)\) count primitive pairs, with \(a<b\) and no
common cut, up to isomorphism acting on the vertex set
**simultaneously**. Different marked arcs producing the same
unmarked target class can be distinct pair types. Write
\(B_{a,b}(z)\) for the size series without the primitive restriction,
and \(T_d(z)\) for the ordinary size series at count \(d\),
including the empty tournament only at \(d=1\). The unique
common-cut parse gives the exact formula
\[
B_{a,b}(z)
=\sum_{d,e\ge1:\ de\mid\gcd(a,b)}
          T_d(z)\,T_e(z)\, U_{a/(de),\,b/(de)}(z).
\]
Reverse construction is by ordinal concatenation with the
same exterior factors in both tournaments; the insertion does
not make any internal cut of a primitive pair common.
In particular
\[
 B_{N,N+2}(z)=\frac{U_{N,N+2}(z)}{(1-z)^2},\quad
 U_{1,3}(z)=z^3,\quad U_{3,5}(z)=2z^4.
\]
Every primitive carrier has at least one strong endpoint, so
there are only finitely many for any fixed \(a,b\): rotating
a Hamilton cycle bounds its order by \(\max(a,b)\),
with the size bound sharpened by
[Busch's strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
For \(N\ge5\), both strong endpoints give the smaller bound \(N\).
Thus \(B_{N,N+2}(z)\) is rational, and its coefficients are
eventually exactly an affine function with slope \(U_{N,N+2}(1)\).
If the polynomial is nonzero its lowest degree and leading
coefficient are precisely the minimum vertex order and the
multiplicity of the smallest literal joint seam types.
The two cannot be recovered just by taking the minima of
the separate scalar count series.

For \(E_p\pm1\) at \(p\ge3\), this reduces the search to joint
strong cores with no exterior padding. In particular
[the 60-vertex exponent-31 pair](../chanpak_e31_physical_20261122/README.md)
is primitive, though not known to be minimum here. The grammar
does not assert that its primitive polynomial is nonzero at
every exponent. What is needed is an order-growing construction
inside the strong primitive carrier, not more exterior padding
or another independent scalar attainment.

## One missing pair: a recurrence target outside tournaments

There is an equivalent path-cover state on an oriented graph \(D_0\)
with exactly one missing unordered pair \(\{a,b\}\). Write \(C=H(D_0)\).
Let \(X\) count ordered pairs of disjoint directed paths covering the
vertices, with the first ending at \(a\) and the second beginning at
\(b\). Define \(Y\) with \(a,b\) interchanged. Both paths are nonempty.
These are marked two-slot reads of the *incomplete* digraph; for example
\[
 X=\sum_{\substack{W\subseteq V\\a\in W,\ b\notin W}}
      H_{\mathrm{end}=a}(D_0[W])\,
      H_{\mathrm{start}=b}(D_0[V\setminus W]).
\]
Cutting a full path precisely at the missing arc shows that the two
possible completions have counts \(C+X,C+Y\). Thus for a strong unit
seam, after relabeling the marked endpoints, \(X>0\) and \(Y=X+2\);
the exact center is \(C+X+1\), not \(C\) alone.

This suggests an explicit induction contract on incomplete oriented
graphs: carry the exact center \(M=C+X+1\), half-gap
\(g=(Y-X)/2\), and the current power \(B=2^p\). A valid lifting
operation taking \((M,g,B)\) to \((4M+B,g,2B)\), and keeping the
two completions strong, would advance a perfect-center seam by one
exponent since \(E_{p+1}=4E_p+2^p\). This is a **sufficient open
construction target**, not an existing operation. Neutral singleton
padding keeps all three path reads unchanged, while an ordinal block
of count larger than one multiplies the gap, so neither can be the
lifting operation. The exact path-cover state indicates what an
active internal substitution still has to control.

The finite replay checks the parsing identity and the two
weak exceptions over all pairs through order six. The statements
above follow from the arguments, not by extrapolating that census.
