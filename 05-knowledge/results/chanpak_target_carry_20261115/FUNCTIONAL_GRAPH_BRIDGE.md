# A perfect outgoing rule is not a proof about its terminal basin

## Collatz as a graph, not a candidate Hamiltonian path

On positive integers let T(n)=n/2 when n is even and T(n)=3n+1
otherwise. There is exactly one outgoing edge at every node. The
incoming edges to n are from 2n and, exactly when n=4 mod 6, from
the odd integer (n−1)/3. This is a binary preimage graph, with a
possible missing second branch, not one directed path through the
integers. Deleting the outgoing edge of 1 makes 1 a sink. The
result is a spanning tree directed toward 1 **if and only if** the
Collatz conjecture holds. If all vertices reach 1 the modified graph
is connected and has no undirected cycle: any finite undirected
cycle would force its vertices' unique outgoing edges all to stay
on the cycle, yielding a directed cycle. Conversely, connectivity
with 1 in this functional graph forces every vertex to reach it:
the property of eventually reaching the sink is constant across
each edge, hence across a connected component. No assertion of
that connectivity is made.

For the shortcut map with the odd step (3n+1)/2, the
[Bernstein–Lagarias conjugacy](https://websites.umich.edu/~lagarias/doc/bernstein.pdf)
identifies the 2-adic map with the binary shift through parity
itineraries. At finite precision the natural outgoing edge
reduces precision by one: modulo 2^(r+1) determines the successor
modulo 2^r. It does not induce an ordinary deterministic
self-map modulo 2^r by forgetting that precision cost. The repo's
[ternary residual-state tower](../nextforest_20260926_automata.md)
keeps the affine tail 3^a q+b for exactly this reason: the observed
digit and the true unread-tail state are different data.

## The corresponding carry graph actually has a rank

A carry inverse starts from an absolute target N. Fix an exact sum
candidate, and keep the target's remaining modular debt, the local
digit data, and the remaining number of columns. Backward dynamic
programming marks the transitions that can finish with the correct
sum and residue. At a viable state, choose the first marked outgoing
edge in a fixed lexicographic order. That is its one right edge.
Reducing the remaining length proves termination for the resulting
graph. Before the viability test there is generally no single
universally correct greedy edge.

The planar non-properness comparison suggests keeping a second
observation at every such step. For a chosen sum s and remaining
digit interval, maintain bounds on the **actual** difference x−y
and hence on the concave scalar product term; maintain lower and
upper bounds on the residual cut separately. Eliminate a state
when their sum cannot contain N, even if every modular condition
survives. This is a finite exact integer inverse with interval
pruning; no efficiency claim follows, because the interval test
still need not find any viable state for a missing target.

The open problem is to exhibit an absolute basin for every
admissible N, and to do so without unfolding the entire
carry-state graph. This is structurally parallel to asking about
the basin of 1, not a reduction in either direction. On the
tournament side, lexical choices, residue closure, and the true
integer basin are three distinct results.
