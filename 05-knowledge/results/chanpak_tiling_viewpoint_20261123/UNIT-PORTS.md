# Unit gaps in an ordinal SCC ladder

**PROVED** structural exclusion; **FINITE-EXACT** marked seed table and product scan.

Let a single reversal change a Hamilton-path count from \(h\) to
\(h+2\). The [rooted cross-component count](STRONG-CRT.md#one-literal-reversal-two-integer-channels) writes gain
as outside scalar \(R\) times the appropriate marked two-slot product,
also multiplied by \(Q_2\) for every intervening component, with the
old-edge loss subtracted in the adjacent case. That loss is exactly the
product of the two ordinary marked endpoint counts: the old word traverses
its entire first block ending at the first mark and its entire second block
beginning at the other mark. This also gives the formula when there are no
intervening blocks.

An internal reversal scales by the heights of *all* other SCCs, so
no other factor can have odd height greater than one. A cross reversal
likewise has no nonsingleton outside its interval. If any nonsingleton
lies strictly inside, the gain is divisible by four, since in a
nontrivial strong tournament the two-slot count is divisible by four:
complement pairing groups the odd summands of \(Q_2\) into an even
number of pairs.
Thus any nontrivial SCCs must be endpoints. When both endpoints are
nonsingleton, each rooted two-slot factor exceeds its ordinary rooted
endpoint count by at least one: split off the marked vertex for an
injection of old paths, then take the full path starting or ending at
the now-required marked vertex as an additional cover. Strong
tournaments have Hamilton paths with any desired endpoint by rotation
of a Hamilton cycle. Thus with adjacent nonsingleton endpoints, the
gain before the outside multiplier is at least
\((L_A+1)(L_B+1)-L_AL_B\ge3\), and it is even. With an interior
singleton it is already at least eight. Hence h has at most one
nonsingleton SCC.

If h has one strong core and the increasing reversal is exterior,
the core must lie at one endpoint. The other endpoint is a singleton.
They must be adjacent: intervening singletons multiply a rooted count
already at least two. Take the singleton as a source by converse.
Reversing \(s\to u\) creates only new cycles through \(s\), which
must use \(u\to s\). By [odd-cycle packing](https://arxiv.org/abs/2307.05569), every new odd cycle contributes twice the path count of its complement. Indeed no original cycle used a cross-SCC arc, and a packing uses at most one new cycle since all of them contain the reversed arc. A gain two thus requires exactly
one old in-neighbor \(a\) of \(u\) (otherwise there are multiple new
triangles), and forbids any simple path \(x\to y\to a\to u\) in
the core (otherwise it contributes a new five-cycle).
Every old vertex \(y\to a\) is consequently a source of the
tournament with \(\{a,u\}\) deleted, so at most one such \(y\) exists.
Strong connectivity requires one. Now \(a,u,y\) form a cycle and
dominate every other vertex, again by the two unique in-neighbors and
the source property of \(y\). Thus a strong core can only be that
triangle. Its gain two does occur. The exterior unit gaps are solely
\(1\to3\) and \(3\to5\), with outside singleton padding.

In particular whenever \(h>3\), a unit edge reverses an internal
arc of one strong core, and the target's sole nontrivial SCC has
identical support: an internal increase cannot split the SCC, since
the reverse cross-component edit would be nondecreasing.

For the six canonical transitive-precedence seeds (heights
5,17,37,77,129,751), independent directed Hamilton-word enumeration
finds the only absolute-gap-two internal flips as shown below.
Each output fails transitive odd precedence, witnessed by the
specified triangle in its precedence tournament.

| Initial height | Edited height | Edge | Precedence triangle |
|---:|---:|---|---|
|5|3|02|013|
|5|3|13|023|
|17|15|02|012|
|17|19|03|014|
|17|19|25|012|
|17|15|35|345|

Thus no unit-neighbor port of the six-seed ordinal family preserves
transitive precedence, even allowing arbitrary singleton padding.
Numeric adjacency between two separately constructed ordinal products
is a distinct question. Their heights are odd and coprime if they
differ by two, so the two sets of used generators are disjoint.
One side uses at most three types, giving \(O((\log X)^3)\)
possible pairs below \(X\) and \(O((\log X)^2)\) per dyadic shell.
A finite check found no numeric adjacent pair through \(10^{40}\)
among 755771 distinct product values; the global exclusion of numerical
pairs is not claimed.

Replay the marked-arc table with unit-ports-check.py. The separate
arithmetic scan is finite and uses six nonnegative exponents; it
is not needed for the structural result.
