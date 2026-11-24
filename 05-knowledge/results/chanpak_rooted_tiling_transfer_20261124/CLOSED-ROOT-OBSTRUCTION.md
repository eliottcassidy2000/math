# Closed viewpoints can lock: cycles are not a substitute for paths

**PROVED** for the structural, extremal and probabilistic statements;
**FINITE-EXACT** for the small census. Companion to
[the one-vertex path classifier](README.md). All moves below are
changes of viewpoint in one fixed tournament, not reversals of arcs.

## The quotient on the closed page

A closed path word reads a directed Hamilton cycle, with a cut edge
omitted. Moving its first vertex to the end rotates the cut, so all
rotations form one connected necklace (up to automorphisms when using
normalized tilings). Restrict the one-vertex re-rooting graph to
closed words and contract these necklaces. The resulting graph is
exactly the graph of cycle orbits under **cyclic erasure and
reinsertion**: erase one vertex so the remaining circular word is
still directed, then put it in a valid circular gap. Indeed an
interior erasure in a path requires the directed shortcut. An
endpoint-to-interior move from a closed root to another closed root
also requires that shortcut, now as the new closure edge. Endpoint
rotations do not change the cycle. Conversely any cyclic move can
be read with a common retained vertex outside the two changed cut gaps first, giving a legal path
move between closed words. These arguments are equivariant under
automorphisms and hence apply to ordinary tiling fibres.

## Sharp lock threshold for every order

For a directed cyclic word of length \(n\), define a legal
\(r\)-erasure to leave a directed cycle of length at least three,
and permit the erased vertices to be reinserted in any gaps after
that. Here \(n\ge5\). A cycle has **no** erasure of size \(1,\ldots,k\)
iff every forward distance \(2,\ldots,k+1\) along the word has a
backward chord (for \(k\le n-3\)). If a legal erasure exists, take
a nonempty consecutive deleted run: its shortcut must be forward.
Conversely any forward shortcut across a run gives the corresponding
legal erasure by deleting just that run. This characterizes feasibility
of erasure, not necessarily whether a rearrangement will change the
cycle orbit.

Every cycle permits some such erasure of at most \(\lfloor(n-1)/2\rfloor\):
for even length use a diameter, choosing its forward direction; for odd
length use the shorter or longer of the two near-diameter directions
according to the arc. Both arc directions span a retained cycle of
at least three vertices. The bound is sharp, and each sharp example
below has an actual rearrangement escaping its original orbit
using exactly that many erased vertices.

**Odd \(n=2k+3\), \(k\ge1\).** Use the circulant with forward differences
\(1,-2,\ldots,-(k+1)\). The \(+1\) cycle has no legal erasure up
through \(k\). Delete \(k+2,\ldots,2k+2\); the survivors \(0,\ldots,k+1\)
close. For \(k\ge2\), reinsert to get
\[
 (0,1,\ldots,k,\ 2k+2,\ k+2,k+3,\ldots,2k+1,\ k+1).
\]
This is directed, and vertex \(k\) can now be erased: the shortcut
from \(k-1\) to \(2k+2\) is forward. For \(k=1\) use
\((0,3,1,4,2)\), whose first two-step shortcut is forward.
Thus the new orbit is different even after quotienting by
automorphisms. The \(+1\) and \(-2\) cycles also give distinct
**blue** closed roots in the same class, independently of the
displayed rearrangement. Both step sizes are units modulo \(n\)
and position reversal is of the form \(x\mapsto c-x\), an
anti-automorphism of the circulant. Any isomorphism of those two
cycles would be \(x\mapsto c-2x\); multiplication by \(-2\) sends
\(-2\) to the backward difference \(4\) for \(n=5\) or \(n\ge9\),
and sends \(-3\) to backward \(6\) at \(n=7\).

**Even \(n=2m\), \(m\ge3\).** Use the same step-one direction and
send distances \(2,\ldots,m-1\) backward. Choose diameters
\(m\to0\) and \(m-1\to2m-1\); the remaining ones can be set from
their lower to higher endpoint. The \(+1\) cycle has no erasure
up through \(m-2\). Delete \(m+1,\ldots,2m-1\), leaving the
directed cycle \(0,\ldots,m\). For \(m\ge4\), reinsert to get
\[
 (0,1,\ldots,m-1,\ 2m-1,\ m+1,m+2,\ldots,2m-2,\ m).
\]
The four non-unit differences are \(m,-(m-2),-(m-2),m\).
The shortcut from \(m-2\) to \(2m-1\) is forward, so \(m-1\) is
now removable. For \(m=3\) use \((0,1,4,2,5,3)\);
the shortcut \(4\to5\) makes vertex 2 removable. Hence even
orders have the same sharp formula. In fact every maximally locked cycle
escapes at the bound. At odd order the lock fixes the entire
circulant. At even order only the diameters are undetermined;
rotate so \(m\to0\). For \(m\ge4\), the alternate arrangement
\((0,\ldots,m-2,2m-1,m+1,\ldots,2m-2,m-1,m)\)
uses just this diameter and has the forward shortcut
\(m-3\to2m-1\). At order six the replay checks all eight diameter
choices. The universal bound for **arbitrary**, rather than
maximally locked, cycles above is only for a legal erasure.

## The missing open ports, at smallest orders

At order five there is a split even after normalization. Outgoing
masks \((18,12,9,17,6)\) give \(H=15\), automorphism order five,
and roots \(12,18,41\). All three are blue. Closed masks \(18,41\)
are separate necklaces with open mask \(12\) between them.
At order six the rigid tournament with outgoing masks
\((34,4,9,19,39,14)\) gives \(H=29\) and two closed necklaces.
There is a shortest connection
\[
 \underbrace{757}_{\text{blue, closed}}
 \ \to\ \underbrace{250}_{\text{black, open}}
 \ \to\ \underbrace{470}_{\text{blue, closed}},
\]
with path words \(340512,345120,534120\). The replay enumerates
all ordinary classes through order six: none split through order
four; order five has the first split and no rigid split; order
six has three split classes, two rigid.

A labelled cut of any cycle whose distance-two chords all point
backward has at least two distinct one-move open neighbors.
For the first vertex \(v\), read bits along the remaining path
indicating whether that vertex beats \(v\). The first two bits
are \(0,1\), and the last two are \(0,1\). Therefore there is an
interior \(1\)-to-\(0\) transition into which \(v\) fits. The new
root is open because the shortcut across the old cut is backward.
Apply the same argument to the last vertex by rotation. Two
single-vertex moves of different vertices producing the same changed
word from one word would have to be an adjacent swap, which cannot
preserve directed adjacency in both words. The two open ports need
not be distinct after taking automorphism orbits.

## A deterministic factorial family of locked blue necklaces

The obstruction need not be chiral. Fix \(k\ge1\), put \(m=2k+3\),
and use the locked odd circulant on \(m\) blocks. Replace each block
by a transitive tournament on \(s\) ranked vertices, with all
interblock directions uniform. On the \(N=ms\) vertices, read the
block order \(0,1,\ldots,m-1\) for \(s\) rounds, independently choosing
the permutation of the \(s\) ranks in each block. Each such word
is a \(k\)-locked directed Hamilton cycle: every prohibited shortcut
has the same relative block displacement \(2,\ldots,k+1\) as in
the base cycle. There are exactly \((s!)^m/s\) such necklaces,
since exactly \(s\) cyclic cuts begin in block 0.

The automorphism group of the blowup acts diagonally on the blocks,
so has order at most \(m!\). Indeed outdegree distinguishes the
internal ranks (the base circulant is regular), and on each rank
an automorphism must act as an automorphism of the block tournament.
Between two ranks in increasing order the compatibility matrix is
\(B+I\), with \(B\) the base adjacency matrix. This matrix has
distinct columns: between any two columns, the row for the vertex
losing their mutual edge distinguishes them. Simultaneous action
by one rank's automorphism preserves \(B+I\), so comparing it with
the action on the other rank forces those two block actions to agree.
Consequently at least \((s!)^m/(s\,m!)\) cycle orbits are locked.

To count a blue part, use the anti-automorphism
\(\alpha(j,r)=(-j,s+1-r)\), with blocks modulo \(m\) and ranks
\(1,\dots,s\). On positions indexed \(i=mt+j\) require \(\alpha\)
to act as \(i\mapsto ma-i\bmod N\), choosing \(a\) odd if \(s\)
is even. Each of \((m-1)/2\) paired blocks freely chooses its
rank permutation, and its mate is forced. In block 0 the rank
permutation must conjugate the round reflection \(t\mapsto a-t\)
to internal rank reversal. Their fixed-point counts match; the
number of conjugating bijections is
\(2^{\lfloor s/2\rfloor}\lfloor s/2\rfloor!\).
An invariant cut exists: solve \(2i\equiv ma+1\bmod N\)
for its first position, possible by the choice of \(a\).
Thus at least
\[
 \frac{(s!)^{(m-1)/2}\,2^{\lfloor s/2\rfloor}\lfloor s/2\rfloor!}
      {s\,m!}
\]
distinct cycle orbits have a blue closed root and are isolated
under the \(k\)-erasure rule. This gives \(\exp(\Omega(N\log N))\)
blue isolated necklaces along the divisible orders, and shows that
the numerous locks are not a consequence of missing reflection
symmetry.

## Many isolated necklaces, even without symmetry

For a fixed labelled cyclic word, impose the step-one directions
and require distances \(2,\ldots,k+1\) backward, where
\(n\ge2k+3\). These are \((k+1)n\) distinct unordered arcs.
By the erasure characterization, the expected number of these
locked directed cycles (modulo rotation) in a uniform tournament
is exactly
\[
 (n-1)!/2^{(k+1)n}.
\]
For every fixed \(k\), this expected count remains
\((1-o(1))(n-1)!/2^{(k+1)n}\) after discarding all tournaments that
are nonrigid or self-converse. In fact the statement is uniform
for \(k\le(1/6-\epsilon)n\), for any fixed \(\epsilon>0\).

To prove this, condition on the imposed directions for a single
cyclic word. Tournament automorphisms have odd order: an involution
would swap a vertex pair and reverse its arc. A candidate odd-order
permutation moving \(s\) vertices moves at least \(s(n-1)/2\)
unordered pairs and has edge orbits of length at least three on
those pairs. Invariance would therefore impose at least
\(s(n-1)/3\) independent binary constraints on all directions.
At most \(2(k+1)s\) of the imposed edges touch moved vertices,
so conditioning loses at most that many constraints. Summing
\(n^s2^{-s((n-1)/3-2(k+1))}\) over supports \(s\ge1\) tends to zero
in the stated range. If the tournament is rigid, a converse
isomorphism is an involution. It cannot fix more than one vertex,
so at most \(\lfloor n/2\rfloor\) unordered pairs are fixed.
On all other pairs its two-element orbits each impose one equation.
The unconditional exponent is at least
\((\binom n2-\lfloor n/2\rfloor)/2\); conditioning loses at most
\((k+1)n\) and there are at most \(n!\) candidates. This also tends
to zero in the stated range. Consequently linearity of expectation
gives a rigid chiral tournament with at least
\((1-o(1))(n-1)!/2^{(k+1)n}\) isolated necklaces for the \(k\)-erasure
rule, at every sufficiently large order. All its roots are black.
For example the bound is
\(\exp(\Omega(n\log n))\) even with \(k\le(1-\epsilon)\log_2 n\).

The path classifier still connects every fibre. What fails is a
simpler classifier that tries to stay on the closed page, even
with a bounded budget of simultaneous cyclic erasures. For
Chan–Pak this is an obstruction to discarding open ports, not a
new restriction on the ordinary path-count spectrum.

Replay: python3 closed-root-check.py.
