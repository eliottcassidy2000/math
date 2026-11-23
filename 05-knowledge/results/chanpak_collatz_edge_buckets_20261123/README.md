# Collatz edge buckets and the broken 61–63–65 staircase

**PROVED scoped / FINITE-EXACT. Chan–Pak remains open.** Ordinary counts
61, 63, 65 are not perfect-number exponents. The central distinction is
between an arithmetic digit, a coherent trajectory in a guarded quotient,
and a reusable marked move on an actual tournament.

## 1. Centered ternary descent and the two holes

Use the three children \(3h-2,3h,3h+2\) of a positive odd \(h\).
In half-index \(k=(h-1)/2\) these are \(3k,3k+1,3k+2\).
Thus every positive odd \(N>1\) has a unique smaller odd parent,
\(2\lfloor(N-1)/6\rfloor+1\). This is deletion of a ternary digit,
not the Collatz odd-step rule. No convergence conjecture is involved.
For tournaments, only the middle digit is free: take the ordinal join
with a directed triangle, multiplying the count by three.

The two absent nodes 7 and 21 are consecutive on the middle-digit
ray \(7\to21\to63\to189\to\cdots\). Deleting them leaves precisely
six rooted components, at \(1,19,23,61,63,65\); indeed 7's other
children are 19 and 23, and 21's children are the last three. In
half-index the cuts are \(3=(10)_3\) and \(10=(101)_3\).
This is a structural role for the pair, not a proof that there
are no other holes. The roots were already isolated in the
[half-index boundary](../fresh_unit_prime_gluing_20261116/SPINE_AND_HALF_INDEX.md).
Combined with the [SCC-grade ray law](../chanpak_window_endgame_20261123/README.md),
it yields an asymmetric bridge: multiplying 63 repeatedly is legal,
but the grade at \(7\cdot3^r\) stays exactly \(r-1\), not \(r\).
The new 63-core carries the otherwise illegal factors together;
it does not give a physical division-by-three rule.

For comparison the actual accelerated Collatz map is
\[
 C(x)=\frac{3x+1}{2^{v_2(3x+1)}}\quad (x\text{ odd}).
\]
Here \(21\to1\) in one odd step, while \(7\to11\to17\to13\to5\to1\).
Likewise \(61\to23\) in one step, \(63\) reaches 61 after 34 odd
steps, and \(65\to49\to37\to7\). In this Collatz graph neither
7 nor 21 has an exceptional status tied to tournament feasibility:
21 cannot be an image, as is true of *every* multiple of three,
including the attained 63; 7 can be an image, as the 65 word shows.
The full three numerical paths are checked in the Python replay.
The [inherited inverse-Collatz word rays](../chanpak_tight_mod_six_20261114/COLLATZ_RAYS.md)
really transport tournament contact words and have all joint 2–3 phases,
but neither ray hits the absolute target 63.

## 2. What the 61, 63, 65 window physically contains

The [complementary-flank chart](../chanpak_tight_mod_six_20261114/TIGHT_COMPLEMENT.md)
already gives strong witnesses for all three heights, at
\((x,y)=(5,10),(6,9),(7,8)\). That is a window of distinct
graphs, not a promise that successive graphs differ by one arc.

**All-order literal obstruction.** No tournament of count **63**
can be changed to count **65** by reversing a single arc. In contrast,
there is a strong eight-vertex tournament of count 61 with *different*
arc flips to 63 and 65. Its out-neighborhood bitmasks, in vertex order
0 through 7, are
\[
 (158,28,24,16,32,79,159,62).
\]
The flips of vertex pairs \(\{1,5\}\), \(\{5,7\}\), or both give
63, 65, and 69, respectively. Thus the path from 63 to 65 in the
single-arc graph exists, but its distance is two, and it can go down
to 61 before going up. Alternatively it can go above both to 69.
Both routes remain strong. Even this first three-point local arithmetic
window does not give monotone \(+2\) continuation.

**Why eight vertices suffice for the nonexistence proof.** In an ordered
SCC decomposition the path counts multiply, and nonsingleton strong
components have count \(>1\). Every nontrivial proper factorization of
63 uses a factor 7 or 21. Both are excluded, so there is exactly one
nonsingleton SCC. [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
\(H(S)\ge5^{(|S|-1)/3}\) implies any strong component of order at
least nine already has count \(>65\). The exact anchored census finds
no 63 at order at most seven, so the strong core has order eight.

The core-locality theorem in §3 shows that the only possible unit
reversal is internal to that eight-vertex core.

**Finite certificate.** On eight labels fix the path
\(0\to1\to\cdots\to7\) and run over all \(2^{21}\) other orientations.
Every height-63 graph can be labelled along an existing Hamilton path.
There are 126 anchored masks of count 63 (two strong isomorphism classes).
All 28 pairs are then reversed individually, *including* any fixed-path
arc. All 126 masks have a 61 neighbor and none has a 65 neighbor. The
DP computes paths starting at \(i\) on \(S\) by
\[
 D(S,i)=\sum_{j\in S\setminus\{i\}: i\to j}D(S\setminus\{i\},j),
 \qquad D(\{i\},i)=1.
\]
The counted masks are not unlabelled multiplicities. Canonical marking
gives two reduced classes of \(61\leftrightarrow63\), both strong on
both ends. By locality, the ordinary order generating function for
these marked seams is \(2z^8/(1-z)^2\): singleton padding before and
after the unique strong core is the only freedom. The same census
finds two strong marked \(61\leftrightarrow65\) classes at order eight.
Replay instructions are below.

This failure is specific to local continuation of the low three-point
window, not to ternary exits from a tournament of count 63. In the
[archived marked fork](../sixty_three_bridge_20261113/INCOMING.md),
\(63\mapsto187,189,191\) is physically witnessed, with a neutral
rewrite of old arcs to reach the plus branch; no common renewable
interface is known. [Session 302](../session302_ternary_ray_birth_20261113/README.md)
also gives a strong, finite-resource middle chain \(63,189,567,1701\).
Neither claim asserts an edge from count 63 to 65.

## 3. A general no-go for the cheapest ternary lift

**General product obstruction.** A tournament with at least two
nonsingleton SCCs has no literal neighbor of count \(H\pm2\).
If an internal reversal leaves another nonsingleton component
untouched, both counts share its odd factor \(>1\). For a cross
reversal, the same applies if there is an untouched nonsingleton SCC
outside the interval of the two endpoints. Otherwise at least two
strong components lie in that interval. No old directed cycle is
destroyed by a cross reversal. If an interior SCC is strong, each of
three distinct vertices in it supplies a new triangle. If there is
no such interior SCC, both endpoint SCCs are nonsingletons. Use one
out-neighbor of the left contact and one in-neighbor of the right
contact to get two triangles. For a third odd cycle, use a two-edge
path starting at the left contact in its SCC, a cross edge to an
in-neighbor of the right contact, and the last two edges of the cycle.
All three cycles use the reversed edge. Such a two-edge path exists
since the left contact has an out-neighbor and that neighbor also has
an out-neighbor different from the contact.
The [odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
adds at least 2 for each of the three new single-cycle packs;
all old packs remain. So the count increases by at least 6.

**Core locality of unit steps.** If both counts at a literal
single-arc seam are at least 5 and differ by 2, each endpoint has
exactly one nonsingleton SCC. Both strong cores have the *same
vertex set*, and the seam reverses an arc inside it.

Indeed the product obstruction gives at most one strong core,
and the count assumption gives at least one, of order \(\ge4\).
An outside reversal not involving or spanning the core leaves
its odd factor. A reversal spanning it is a cross reversal
and creates a triangle through each vertex of the core.
For a reversal from a contact \(a\) in the core to an outside
singleton \(b\), let \(a\) be the earlier endpoint (take the
converse for the other case). Every out-neighbor of \(a\) in
the core gives a new triangle. If \(a\) has only one such
neighbor \(u\), there is a directed 3-edge path \(a,u,v,w\)
in the core. Otherwise \(u\) could have only one out-neighbor
\(v\) (two would give a path); the remaining vertices would
all beat \(a,u,v\), contradicting strong connectivity. Now
\(a,u,v,w,b,a\) gives a new 5-cycle in addition to the triangle.
Consequently the cross reversal increases the count by at least
4. This excludes all off-core reversals. Apply the result from
both ends of the seam: reversing inside the first core cannot
enlarge it, and reversing inside the second cannot enlarge
that one either, so their vertex sets agree.

In particular even unlimited **ordinal singleton padding between
unit moves** cannot turn a finite seed into a strictly increasing
unbounded \(+2\) staircase above count 3: the order of the unique
strong core would remain fixed. This does not constrain neutral
rewrites or larger jumps that can supply a genuinely new core.

**Sharper estimate for the triangle.** Let \(L\to C_3\to R\) be an ordinal join with a directed triangle
and \(h=H(L)H(R)>1\), using \(H(\emptyset)=1\).
**No single reversal changes \(3h\) to either \(3h-2\) or \(3h+2\).**
This is a uniform restriction on a specific natural neutral host,
not on all \(3h\)-count hosts.

A reversal inside \(L\) or \(R\) leaves a multiple of three; one
inside the triangle gives \(h\). For a reversal \(a\to b\)
between \(A=L\) and the triangle, let \(E_A(a)\) count full paths
ending at \(a\). Let \(U_A(a)\) count full paths starting at \(a\)
plus ordered splits of \(A\) into two nonempty directed paths
whose second path starts at \(a\), with no adjacency demanded
between the two runs. Splitting a full path immediately before \(a\)
gives \(U_A(a)\ge H(A)\). On the triangle side the corresponding
end read is 3: one full path ending at \(b\) and two ordered
two-piece covers with first piece ending there. Paths gained from the
reversed edge consist of an optional initial \(A\)-run, the first
triangle run ending at \(b\), the \(A\)-run starting at \(a\), and
an optional final triangle run. Thus the exact net change is
\[
 [3U_A(a)-E_A(a)]H(R)\ \ge\ 2h .
\]
Reverse all arrows for the analogous triangle–\(R\) case. The lower
bound is sharp when \(a\) is a sink of \(A\) (or, by converse, the
contact in \(R\) is a source). A reversal of an edge skipping the
triangle destroys no old Hamilton path. For each full path of \(L\)
and \(R\), split \(L\) before its contacted vertex and \(R\) after its
contacted vertex. Traverse \(L_{\rm prefix},R_{\rm prefix},
L_{\rm suffix},R_{\rm suffix}\) using the reversed edge; insert a
full triangle path either before \(R_{\rm prefix}\) or after
\(L_{\rm suffix}\). These give \(6h\) distinct gained paths. All
three cases prove the claimed exclusion.

The loophole is the useful part: an exact \(\pm2\) lift must
change the marked core, or first produce a different, non-ordinal
host of count \(3h\). Comparing only the number \(3h\) forgets this
distinction. In particular gluing on a real C3 is a valid product
constructor but not a reusable local additive constructor.

## 4. An unbounded positive rail

The [strong rail templates](STRONG-RAILS.md) evade the product
obstruction. For every \(w\ge1\) they give a strong \(3h\) host
and one literal arc taking it to \(3h+2\), with \(h=3\cdot2^w-1\)
attained by an independent template. A second slope gives
\(h=28\cdot2^w+1\). The arc leaves the weight of two visits to a
uniform module unchanged and adds exactly two to the weight of
one visit. It is the latter coordinate, not the count alone,
that carries the unit correction. Separate strong templates
also realize both signed exits and the neutral one for
\(h=2^{w+1}+1\); we do not identify those with a shared arc
or a renewable state.

## 5. Next constructive task

Keep both a target integer and its marked continuation state.
The six arithmetic seeds are already witnessed, but the 61 fork and
63 obstruction forbid interpreting the edges of their integer graph
as automatic local moves. For completeness one can also bypass the
free middle digit and work only on the prime and single-seven lanes
of the inherited grade reduction. In either approach, the missing
input is a computable non-ordinal port or a carry repair that returns
the new graph to a specified ready class. A readout of a preselected
tournament in the [Collatz buckets](COLLATZ-BUCKETS.md) is exact but
does not synthesize that graph or the ready state.

## Reproduce

    c++ -O3 -std=c++20 seam-census.cpp -o /tmp/seam-census
    /tmp/seam-census 7
    /tmp/seam-census 8 > /tmp/seam-census-8.out
    python3 canonical.py /tmp/seam-census-8.out
    python3 replay.py

The census checks the order-seven and order-eight universes;
smaller orders can pad to seven without changing the count.
The Python replay checks the explicit four-corner fork, both
guarded encodings, and the actual positive Collatz links.
