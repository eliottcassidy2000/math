# One new reader vertex cannot close a bounded carry at a large baseline

Let \(T\) be any tournament, \(T'=T+v\) any induced extension by
one vertex, and \(\Delta=H(T')-H(T)\). Write \(U=N_{T'}^+(v)\),
\(W=N_{T'}^-(v)\). Use count \(1\) for the empty tournament.
We prove the order-free bound
\[
       \Delta>0\quad\Longrightarrow\quad
       H(T) \ \le\ (1+1/k)\Delta\ \le\ 2\Delta,
\]
where \(k\) is the number of old \(U\to W\) arcs (necessarily
positive). In particular, for \(k\ge2\) the bound is \(3\Delta/2\).
For every \(\Delta\ge6\) one has the sharper universal bound
\(H(T)\le2\Delta-3\), attained for infinitely many gap values.
This applies even when both the entire old tournament and the
new reader's contact mask may depend on the target. In particular a
bounded nonzero carry does not give a universal **monotone one-vertex**
move, not even with changing cores. The exact two smallest gaps are
\[
\begin{array}{c|c}
 \Delta &(H(T),H(T'))\\ \hline
 2&(1,3),(3,5)\\
 4&(1,5),(5,9).
\end{array}
\]
This is not an obstruction to edge flips, replacements, or joint
multi-vertex moves of zero or unbounded delta.

## The sharp linear bound from two literal decompositions

For each old Hamiltonian path \(P\), insert \(v\) at any
directed seam where the preceding old vertex lies in \(W\)
and the succeeding one lies in \(U\). The initial position
is allowed if \(P\) begins in \(U\), and the final position
if \(P\) ends in \(W\). Thus the number of insertion positions is
exactly \(1+\#(U	o W	ext{ adjacencies in }P)\):
in a binary word the difference between the two directed
transition counts telescopes to its endpoint bits.
The resulting paths are distinct across positions and across
\(P\), since deleting \(v\) recovers \(P\). All old paths
except those with every \(W\) vertex before every \(U\) vertex
have at least one \(U	o W\) adjacency. There are at most
\(h_0=H(T[U])H(T[W])\) exceptional paths. Therefore
\[
                         \Delta\ge H(T)-h_0 .
\]

For the complementary estimate \(\Delta\ge k h_0\) inspect
the odd-cycle packings in the next section. Combining both
estimates gives \(H(T)\le (1+1/k)\Delta\).

## Finite-set readout and an independent selector-surgery bound

Use the tournament specialization of the
[Grinberg–Stanley odd-cycle packing formula](https://arxiv.org/abs/2307.05569):
\(H(G)\) counts collections of vertex-disjoint directed odd cycles
(length at least three), with every cycle freely labeled by a bit.
Vertices not in a cycle remain singleton atoms. This is a concrete
finite-set use of the project's [field-with-one-element heuristic](https://github.com/eliottcassidy2000/math/blob/main/07-reflections/the-fundamental-truths-in-every-language.md):
the cycle-free term is the single empty packing, and choices are
real finite subsets, not coefficients of a formal interpolation.
No assertion about a literal field is needed.

The packings ignoring \(v\) are exactly the old packings; hence
\(\Delta\) counts packings containing \(v\) on a cycle. Let \(C\)
be the number of old labeled packings having a cycle meeting
both \(U\) and \(W\). Then \(H(T)=H(T[U])H(T[W])+C\).
Fix any labeling of vertices to break ties. For each such old
packing select deterministically one crossing cycle, and one
vertex \(x\) on it whose predecessor \(w\) lies in \(W\) and whose
successor \(u\) lies in \(U\). Such \(x\) exists: the length of the
cycle is odd, so moving two positions permutes it cyclically,
and some two-step move must go from \(W\) into \(U\). Replace \(x\)
by \(v\), keeping its cycle label and all the other cycles and labels.
The new cycle is valid because \(w\to v\to u\).

For a fixed image packing, any preimage is determined by a single
unused old vertex \(x\) with \(w\to x\to u\). Write \(m\) for the
number of these eligible vertices. For each even subset \(S\) of
this \(m\)-set, make another packing counted by \(\Delta\): insert
a Hamiltonian path of \(S\cap U\), if nonempty, between \(v\) and \(u\),
and a Hamiltonian path of \(S\cap W\), if nonempty, between \(w\) and \(v\).
These paths exist by Rédei, and their endpoints can be joined across
the two seams because every eligible \(x\) has the two displayed
arcs and has the stipulated contact to \(v\). The modified cycle
is still odd. Distinct \(S\)'s use distinct vertex sets; its
cycle label has two choices. For \(m\ge1\) this gives \(2^m\)
distinct new objects, hence \(m\le\lfloor\log_2\Delta\rfloor\).
No disjointness across different fibers is assumed. Therefore
\[
                           C\le \Delta\lfloor\log_2\Delta\rfloor .
\]

To bound the noncrossing term, note \(\Delta>0\) implies some old
arc \(a\to b\) goes from \(U\) to \(W\), since removing \(v\) from
its cycle leaves a path going from \(U\) to \(W\). For such an arc,
let \(P_a\) be the unsigned weighted subset-path count ending at
\(a\) within \(T[U]\), and \(Q_b\) the same count starting at \(b\)
within \(T[W]\), with weight the Hamiltonian count of the unused
part. Let their signed versions (weight \((-1)^{|\mathrm{subset}|}\))
be \(-s_a,-t_b\); these are nonpositive by the
[signed endpoint identity](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/READY_CLASSIFICATION.md).
Exactly the apex packings whose apex cycle uses only this one forward
\(U\to W\) arc contribute \(P_aQ_b+s_at_b\), if also restricted to
have no other crossing cycles in the packing. Indeed the apex
cycle consists of one run in \(U\) ending at \(a\), the arc
\(a\to b\), and one run in \(W\) starting at \(b\). Its total old
run length must be even, and the factor of two from the cycle bit
cancels the denominator of the parity projector.
By splitting a full side path at its marked vertex,
\(P_a\ge H(T[U])\) and \(Q_b\ge H(T[W])\). Thus, if \(k\) is the
number of forward cut arcs,
\[
    \Delta\ \ge\sum_{a\to b}(P_aQ_b+s_at_b)\
               \ge kH(T[U])H(T[W]).
\]
Also each such arc supports at least one two-labeled apex triangle,
so \(k\le\Delta/2\). In particular \(\Delta\ge k h_0\), proving the promised linear bound.
For \(\Delta=0\), the contact set is a genuine ordinal cut
(no \(U\to W\) arc), so \(T'\) is \(T[W]\to\{v\}\to T[U]\).

## Exact classification of the two smallest carries

When \(k=1\), all crossing cycles (old or new) must use the same
unique old forward arc \(a\to b\), so old crossing cycles have the
same two side runs as above, reclosed by a backward cut arc.
The missing reverse \(b\to a\) could invalidate a two-vertex choice
only; the odd-length selector for an old cycle already excludes it.
Hence, with \(h_0=H(T[U])H(T[W])\), \(D=P_aQ_b\), \(Z=s_at_b\),
\[
                  H(T)=h_0+D-Z,\qquad \Delta=D+Z.
\]
Here \(0\le Z\le h_0\le D\). If \(Z=0\), at least one
marked endpoint is outside the appropriate extreme strongly
connected component: the first component for \(a\) in \(U\),
or last component for \(b\) in \(W\). Indeed each specified
vertex in that component can be a path endpoint by opening a
directed Hamiltonian cycle, with singleton case trivial.
Then its unsigned count is at least twice that side's ordinary
count. An ordinal predecessor (for an end response) or successor
(for a start response) splits freely into two slots, giving
\(F_2(C)\ge2H(C)\) on every nonempty added portion \(C\); the
other ordinal padding factors normally. Thus \(D\ge2h_0\)
and \(H(T)\le3\Delta/2\). If \(Z\ge1\) then
\(H(T)=\Delta+h_0-2Z\le2\Delta-3\). Together with the
\(k\ge2\) bound, this gives \(H(T)\le2\Delta-3\) when \(\Delta\ge6\).

Sharpness is literal. For \(r\ge2\), take \(W\) to be a
transitive chain of length \(r+1\) plus a marked \(b\) beating
only its first vertex; put \(U=\{a\}\) and make \(a\to b\)
the unique exception to \(W\to U\). Then
\(h=H(W)=1+2^r\), \(Q_b=h\), \(t_b=1\). To verify these,
a marked start is either the singleton \(b\) or follows the
first spine vertex and any subset of the \(r\) followers, while
the unique full path ending at \(b\) follows the whole spine.
The ordinary count also follows by counting the odd cycles
through \(b\): the \(r\) followers contribute \(2^{r-1}\)
odd-cardinality subsets. The unique-cut law gives
\((H(T),\Delta)=(2h-1,h+1)=(2\Delta-3,2+2^r)\).

For classification at the two smallest gaps, if \(h_0=1\), both sides are
transitive, \(D=2^j\) where \(j\) counts predecessors of \(a\)
in \(U\) plus successors of \(b\) in \(W\), and \(Z=1\) exactly
when \(j=0\), and \(0\) otherwise. This gives \(\Delta=2\) with
old count 1 or 3; \(\Delta=4\) gives old count 5. If \(h_0=3\)
then one side is a directed triangle with ordinal singleton
padding and the other side transitive: the packing formula allows
just one directed triangle, so no strongly connected component
besides that triangle can be nontrivial. Any root profile on the
triangle side is a positive multiple of 3: in the triangle it is
3 up to dyadic padding; for a singleton root it has a factor 3
or a factor \(F_2(C_3)=12\), according to which side of the
triangle it lies on. Thus \(D\) is a positive multiple of 3.
For \(\Delta=4\), it must be 3, giving \(Z=1\) and old count 5.
No \(h_0\ge5\) occurs at \(\Delta\le4\). Every claimed
unique-cut endpoint pair is realized by transitive sides and
the specified \(j=0,1,2\).

For \(\Delta=2\) there are never two cut arcs. For \(\Delta=4\),
\(k=2\) is possible only when \(h_0=1\) and each one-crossing
arc gives exactly 2, since each already contributes at least 2.
Thus each such arc is one of
\[
  U_1\to W_{\rm last},\quad U_2\to W_{\rm last},\quad
  U_1\to W_{\rm penultimate},
\]
where each side is listed in topological order, omitting
unavailable endpoints. If the second and third arcs both occur,
then \(U_1,W_{\rm penultimate},U_2,W_{\rm last}\) is an additional
four-vertex run, contributing an extra two-labeled apex cycle.
Otherwise the two arcs share head or tail and the old tournament
is transitive: the active triple has transitive orientations
and all earlier \(W\) or later \(U\) vertices are ordinal padding.
These cases give \(H(T)=1, H(T')=5\), realizable using the
three active vertices alone. This completes the classification.

## Corollary: all induced extensions, however many new reader vertices

For an induced extension \(T\subseteq G\), add the missing vertices
in any order. Each step is nonnegative. If \(H(G)>H(T)\), the first
positive step has baseline \(H(T)\), regardless of its position in
the order. Its gap \(\delta\) is at most the total gap \(\Delta\).
Hence
\[
           H(T)\le2\,[H(G)-H(T)] .
\]
The same argument and the least-gap classification give
\(H(T)\le2\Delta-3\) for total \(\Delta\ge6\): if the first
positive step has gap 2 or 4 then \(H(T)\le5\), and otherwise
use its own sharp bound. Total gaps 2 and 4 have the same endpoint
classification as single-vertex gaps (the only way to split 4
into two positive gaps is \(1\to3\to5\)). Sharpness persists
because each of the sharp single-vertex moves is an induced extension.

An equal count is possible **if and only if** every added vertex
is a singleton strongly connected component in \(G\). Indeed
every zero step is insertion at an ordinal cut and leaves all
previous SCCs intact. Conversely such singleton components can
be removed without changing any other SCC or the product of
ordinary counts. In particular an exact neutral return on strong
tournaments of distinct orders cannot preserve the entire earlier
tournament as an induced subgraph. A genuine edit or replacement
of old material is necessary.
