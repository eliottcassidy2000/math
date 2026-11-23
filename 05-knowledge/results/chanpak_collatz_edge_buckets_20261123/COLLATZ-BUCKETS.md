# How many tournament edges can genuinely be overlaid?

**PROVED, finite exact transfer and its limits.** Write \(C(x)\) for the
accelerated positive odd Collatz map. The guard \(j\ge1\) specifies that
\(v_2(3x+1)=j\); on that domain the successor is \((3x+1)/2^j\).
The finite guard coding below is in the spirit of
[Bernstein–Lagarias's parity coding](https://websites.umich.edu/~lagarias/doc/bernstein.pdf),
but counts exact accepting starts for an arbitrary tournament rather than
making a convergence claim.

## 1. Residue buckets without an edge selector: at most three

Fix any odd \(M\); a bucket is a residue modulo \(M\). Join distinct
buckets \(u\to v\) if at least one actual positive odd Collatz edge goes
between them. For fixed \(j\), the source residue alone determines the
target \(2^{-j}(3u+1)\bmod M\). Conversely the guard condition gives a
single residue modulo \(2^{j+1}\), compatible by CRT with any \(u\bmod M\).
Thus the raw adjacency is exactly
\[
 u\to v\quad\Longleftrightarrow\quad [3u+1]=[v],
\]
where brackets are orbits under doubling modulo \(M\).
Each bucket \(u\) is now an arrow of types \([u]\to[3u+1]\), and
adjacency is composition of these arrows.

An induced tournament, ignoring self-loops, has at most three
vertices. If it contains a directed triangle, its three type arrows
form a 3-cycle on distinct types: equality of two types would force
a reciprocal edge. A fourth arrow has only two endpoints and can
compose with at most two triangle arrows. Otherwise the tournament is
transitive. In a transitive four-vertex subset the consecutive arrows
must be \(a_0\to a_1,\ldots,a_3\to a_4\). The two length-two
forward shortcuts force \(a_1=a_2=a_3\), making the middle two
arrows reciprocal, a contradiction.

A second obstruction concerns an actual finite list with one integer
edge per unordered pair, not a residue quotient. If no switching
between representatives is allowed, each starting integer determines
at most one Hamilton traversal, so there can be at most \(\binom n2\)
such traversals in the list. This already fails to transfer the full
count at order five: a strong tournament on five vertices has 15
Hamilton paths (a concrete example is checked by the replay).
One needs families of actual edges, not a single representative
per pair, to preserve all Hamilton paths in general.

## 2. Exact \(\binom n2\)-instruction overlay, with full trajectories

Take \(M=3^a\) and \(d=2\cdot3^{a-1}\ge n\), and assign distinct
unit residues to the vertices of an arbitrary tournament. The element
2 generates the units modulo \(M\). Indeed its order mod 3 is two,
and cubing \(4^{3^b}\) raises the valuation of \(4^{3^b}-1\) by one,
starting from \(v_3(4-1)=1\).

Use two guarded steps for each directed edge \(u\to v\); both steps
together have fixed total guard \(L=d+1\). Their endpoint is
\[
 2^{-L}(9u+3+2^{j_1})\pmod M,\qquad j_2=L-j_1 .
\]
Choose \(1\le j_1\le d\) so that
\(2^{j_1}=2^Lv-9u-3\bmod M\). The right side is a unit since \(v\)
is. Then \(j_2\ge1\), and both guards and the internal residue are
specified by the arc. Only the vertex residues at macro boundaries
are read, using the same two-step clock for every edge. There is
one selected instruction per unordered pair. The opposite directed
instruction is not allowed.

For any finite guard word \(j_1,\ldots,j_t\) of total \(J\), the set
of positive odd starts with precisely that word is one residue class
modulo \(2^{J+1}\). This follows backwards: if the allowed next
value is \(y\equiv b\bmod2^m\) with \(b\) odd, then
\(3x+1\equiv2^j b\bmod2^{m+j}\) holds for a unique odd residue \(x\);
it ensures the guard is exact. Start from the final condition
\(y\equiv1\bmod2\). CRT intersects this one cylinder with any
prescribed starting bucket \(u\bmod M\). This proves actual
coherence of the accepted trajectories, not just separate
existence of each step.

Every Hamilton path in the tournament now specifies one start
class modulo the **same** number \(P=M2^{1+(n-1)L}\).
Distinct paths specify distinct classes because actual accelerated
steps determine the guards and boundary buckets. Conversely an
accepted start (initial bucket in the vertex set, obeying the
registered guards and visiting every boundary vertex exactly once)
traces a Hamilton path. Thus each complete interval of \(P\)
consecutive starts contains exactly \(H(T)\) accepted positive odds
if the interval is positive. This is an exact scalar transfer;
there is no weighting or asymptotic limit.

## 3. One accelerated step per arc, using fixed endpoints instead

With two extra buckets add a universal source \(s\) and sink \(t\),
leaving the old tournament in the middle. Choose \(d\ge n+2\).
For each allowed arc \(u\to v\) choose the single guard
\[
 j_{uv}=d+\log_2(3u+1)-\log_2(v)
\]
with both logarithms chosen between 0 and \(d-1\).
Then \(1\le j_{uv}\le2d-1\), and it lands in \(v\).
All Hamilton paths start at \(s\) and end at \(t\). Their total
guard is the same integer
\[
 J=(n+1)d+\sum_{u\ne t}\log_2(3u+1)-\sum_{v\ne s}\log_2(v).
\]
This is the endpoint correction that makes path mass uniform.
The same CRT proof gives exactly \(H(T)\) accepted starts per
positive period \(M2^{J+1}\).

For a fixed guard per arc, a quadratic guard length is in general
necessary for this type of realization. Outgoing guards from one
bucket must be distinct, or they'd land in the same bucket. On
an arc-transitive regular tournament of order \(n\), each arc is
in the same fraction \(2/n\) of Hamilton paths by counting arcs
per path. The mean total guard is therefore at least
\((n^2-1)/4\): sum the smallest \(d^+(d^++1)/2\) at each
vertex where \(d^+=(n-1)/2\), then multiply by \(2/n\).
There are unbounded arc-transitive regular examples, e.g. the
quadratic-residue orientations on \(\mathbb F_{3^{2k+1}}\).
This bounds guard cost, not the smallest positive lift.

## 4. Information preserved and lost

The unfiltered quotient loses the direction of most pairs.
The guarded quotient restores adjacency by supplying a registry,
and the CRT read preserves the exact number of spanning bucket
traversals, even though the same bucket along different paths often
uses different actual integers. What is not preserved is the
absolute numerical value of an intermediate Collatz integer, an
unmarked continuation state, or any rule to synthesize the registry
from a desired path count. A guard registry is already carrying
the orientation choice for every pair. Thus this construction
gives an exact translation of an existing tournament, not a proof
that a given \(N\) occurs, and in particular does not explain the
exclusions 7 and 21 by Collatz convergence.

A uniform rule intermediate between the raw quotient and a full registry
also has a testable failure. On the primitive modulus \(3^a\) unit
buckets, allow a single consecutive window of guards globally. In
discrete-log order each out-neighbor set is the restriction of a circular
interval, with the vertex itself removed if present. This cannot give
the Paley tournament on seven vertices (differences 1, 2, 4 modulo 7).
If its out-set is not itself a 3-interval, it becomes one 4-interval
after adjoining the vertex itself, and then its in-set is a 3-interval.
Choose one such in/out 3-interval at every vertex. All seven selected
sets are distinct: in–in and out–out at different vertices intersect
in one point, and equality of an in-set and an out-set forces their
vertices equal by summing modulo seven, whereupon they are disjoint.
Thus the seven selected sets exhaust all seven circular 3-intervals.
Two neighboring intervals intersect in two points, so their types
must alternate between in and out around an odd cycle. Impossible.
Local guard policy, not just one cut in a global clock, is needed for
a universal overlay.
