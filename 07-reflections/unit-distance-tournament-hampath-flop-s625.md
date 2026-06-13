# Unit-Distance Tournament Hamiltonian Path Flop S625

The user's mapping is a good way to make the unit-distance problem feel native
to the tiling model: start from a transitive base tournament on the points, then
flip the tile when a pair is exactly unit distance. Or reverse it and leave unit
pairs unflipped while flipping nonunit pairs.

The first thing S625 had to separate was geometry from base order. A unit graph
can have a Hamiltonian path even if the lexicographic tiling order does not
walk along that path. So there are two possible flops:

- a graph-level flop, where the unit-distance graph no longer has any
  Hamiltonian path;
- an order-level flop, where the fixed tiling order no longer makes an all-unit
  directed Hamiltonian path.

The scout finds the second, not the first. In the tested triangular rows through
`n=22` and Moser rows through `n=14`, the unit graph always has a Hamiltonian
path. But under the canonical unit-flip tiling, all-unit directed Hamiltonian
paths vanish at `n=7`.

That `n=7` threshold is not mysterious. The first six points still admit a
canonical order that can be read as a unit snake. At seven points, the compact
hexagon-with-center structure appears. The unit graph still has many paths, but
the lexicographic base order is no longer one of them. The nonunit complement
also shows the hexagon obstruction: it has a Hamiltonian path at `n=6`, loses it
at `n=7` because the center is complement-isolated, and regains it from `n=8`
onward in the tested compact rows.

This suggests a useful modification of the mapping rule. If the goal is to ask
whether the guaranteed Hamiltonian path is genuinely part of the unit-distance
structure, do not fix the base path lexicographically. First choose a
Hamiltonian snake in the unit graph, then use that as the transitive base order.
Under that snake-base tiling, an all-unit directed Hamiltonian path exists
exactly when the unit graph has one.

The canonical tiling is still meaningful, just for a different question. It
measures alignment between a construction and a fixed coordinate/order system.
Its flop at `n=7` says the compact optimal construction has stopped being a
single monotone unit chain. It does not say the unit chain is gone.

The directed Hamiltonian-path histograms are the richer invariant. In the
canonical unit-flip tournament, the profiles do not jump from all-unit to
all-nonunit. The modes sit in the middle: triangular `n=10` has mode `6` unit
steps out of `9`, while Moser `n=10` has mode `5` out of `9`. The mandatory
path is usually mixed. That is probably where Tournament Analysis has the most
signal: not "which side owns the one path?" but "what is the distribution of
unit-step counts over all directed Hamiltonian paths?"

The recursive picture is boundary-snake persistence. Compact triangular
constructions grow by adding boundary pieces, and a unit Hamiltonian snake can
be threaded around the new perimeter. The complement recursion is different:
it becomes path-rich once enough nonunit chords exist, except when a universal
unit-neighbor center isolates itself from the complement.

So the answer to the user's question, in this finite scout, is:

- graph-level unit Hamiltonian path: no flop observed;
- canonical tiling all-unit directed path: first flops at `n=7`;
- nonunit/complement path: appears at `n=6`, fails at `n=7`, reappears at `n=8`;
- best next invariant: Hamiltonian-path unit-step profile, not a binary label.
