# LRC14 Haar Product Tile Discrepancy Reflection

The bridge that now feels live is: discrepancy theory is not only about
averages being small; it is also about locating the first coefficient that
refuses to vanish.  For LRC14, the recent agent work has already made the
state space local: safe intervals, boundary cocircuits, endpoint-owner pairs,
Ramanujan packets, Fejer atom banks, K33/state-lift debt, and smoothing
guardrails.  The missing language was how local pieces multiply.

The two-dimensional Haar product rule gives that multiplication table almost
too cleanly.  In one dimension, two dyadic Haar functions are disjoint, equal,
or nested.  Their product is zero, an indicator, or a signed smaller Haar
function.  In two dimensions, the rule just happens coordinate by coordinate.
That creates exactly the tournament-tiling menu:

```text
zero                  independent tiles
indicator x indicator boundary atom
indicator x haar      owner strip
haar x indicator      owner strip
same-direction nest   recursive refinement
opposite nest         cross handoff / zipper
```

The finite depth-3 census is a good sanity check.  Most rectangle pairs are
orthogonal zero (`86.4%`).  The same-tile atoms are tiny (`0.44%`).  The three
live signed classes are sparse and perfectly balanced.  So a direct positivity
sum is probably the wrong thing to ask for.  A better proof shape is:

```text
if no lonely open interval exists,
then some typed local Haar coefficient must survive;
if none survive,
then the packet is a pure boundary atom;
pure boundary atoms should be AP/GW or a named state-lift obstruction.
```

This also clarifies the role of the tournament tiling model.  A fixed-path
tournament tile is a two-coordinate object; the Haar rectangle is also a
two-coordinate object.  The product rule says interactions are classified by
coordinate equality, nesting, and crossing.  That is what the tiling model has
been doing combinatorially with strips, pins, staircase anti-diagonals, and
edge flips.

The guardrail is just as important as the analogy.  Haar products can tell us
that a local coefficient exists, but they do not automatically remember why it
is legal.  AP and Goddyn-Wong already teach that lesson: the visible boundary
owner skeleton can be the same while the hidden C27 transfer differs.  So any
Haar-discrepancy proof has to carry endpoint owners, exact-period labels,
C27/GW transfer labels, Fejer packet identity, and K33/state-lift debt.

The next concrete computation should attach this table to actual LRC14 packet
fibers.  I would start with the named rows: AP, GW, `12->36`, `10->20`,
`13->26`, `P10+GW`, `12->168`, and the lcm-tail rows.  For each row, make a
rectangular grid whose axes are endpoint walls and proof clocks.  The target is
not a huge coefficient search; it is a typed vanishing test.  If the owner-strip
and cross-handoff coefficients vanish, what boundary skeleton is forced?

That question is much closer to the recent tope/cocircuit and certificate
handoff work than to raw runner enumeration.  It is a discrepancy theorem with
labels, and that feels like the right level of abstraction for the current
LRC14 proof stack.
