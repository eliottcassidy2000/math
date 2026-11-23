# One arc is not one complement line

**PROVED scoped.** A Hamilton frame \(P\) is a directed spanning path.
The complement-line operation \(\phi_P\) preserves its \(n-1\) arcs
and reverses every other arc. The ordinary single-arc graph, including
its one-coordinate subcube preserving a fixed frame, is not this graph.

## Low-width tape cannot stay on its rails

For a tournament of order \(n\) put
\[
 w_r(T)=\min_{|X|=r}\max\{|K|:K\text{ is an SCC of }T-X\}.
\]
Let \(\ell=\lceil(n-2r)/(2r+1)\rceil\ge3\). If \(w_r(T)<\ell\),
then **every** frame \(P\) satisfies \(w_r(\phi_P(T))\ge\ell\).

Choose \(X\) realizing the small width, and fix any deletion set
\(Y\) of \(r\) vertices in the output. Removing \(X\cup Y\) from the
ordered frame leaves at most \(2r+1\) contiguous segments with at
least \(n-2r\) total vertices. One has length at least \(\ell\).
It lies in \(T-X\) with all original frame arcs intact. Those SCCs
are ordered, and their orders are below \(\ell\), so the first
and last vertices of the segment lie in different components.
The original arc from first to last is forward and is not a
frame arc (length at least three). The complement reverses it.
The segment therefore closes as a directed cycle of length
at least \(\ell\) in \(\phi_P(T)-Y\). The result follows for all \(Y\).

The property is invariant under ordinary isomorphism and converse.
Thus the low-\(w_r\) classes are independent in the actual
complement-line metagraph, even allowing a re-framing before each
line. With \(r=1\), the explicit six-vertex block tapes in the
[unit switch packet](../chanpak_guarded_ternary_packing_20261123/UNIT-SWITCH.md)
have \(w_1\le6\) after sharing the apex. For order \(n\ge21\),
\(\ell\ge7\), so no one complement line directly connects two
such codewords. The ordinary SCC bit tape is an even simpler
\(r=0\) example once \(n>5\).

## Exact two-line identities

For labelled vertices, represent orientations by bits relative
to a fixed total order. A first complement line at \(P\) and
a second at a Hamilton frame \(Q\) of the intermediate tournament
toggle precisely the unordered edge set
\(E(P)\mathbin{\triangle}E(Q)\): toggling all edges cancels
between the two steps, leaving each frame edge's exemption.
The number of altered edges is even because both frames have
\(n-1\) edges. Hence two lines cannot implement a single labelled
arc reversal. This parity assertion is not passed to the
unlabelled quotient, since a relabelling can change orientation
parity relative to the chosen total order.

There is also an exact score-sidecar bound. If \(d(T)\) is
the vector of outdegrees and \(\partial P=e_{\mathrm{first}(P)}
-e_{\mathrm{last}(P)}\), then
\[
 d(\phi_P(T))=(n-1){\bf1}-d(T)+\partial P,\qquad
 d(\phi_Q\phi_P(T))-d(T)=\partial Q-\partial P.
\]
This follows by first taking full converse and then restoring
the arcs in the frame. Thus a two-line move changes scores
only at the four possible endpoints, and preserves the whole
score vector when the frames have the same endpoints.
These rules are necessary, not a sufficient characterization
of reachable orientations.

These are obstructions to a direct native simulation, not
to longer encodings with ancillary tournaments, changed
deletion width, or a controlled high-width excursion.
The marked Turing tape in the main note is an actual
single-arc simulation but already supplies its own transition
program. An unrestricted complement line supplies neither
that program nor a terminal predicate.
