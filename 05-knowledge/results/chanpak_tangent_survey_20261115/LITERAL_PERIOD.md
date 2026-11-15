# A full-alphabet period-length gate for an exact radix ladder

The [preceding syllable sieve](../chanpak_decimal_ports_rank62_20261115/SYLLABLE_BRIDGE.md) excludes a single repeated motif of the 13 floor units as a source of $1000^k+1$. For a motif over **all eight** span-three letters there is a separate finite reduction; it needs no syllable parse.

Fix backward width $b$, literal ends $P,S$, and a nonempty literal period $W$ of length $d$. If $H(PW^kS)=1+B^k$ eventually for $B>1$, then
$$
5^d\le B^3.
$$
In particular every stationary span-three candidate at radix 1000 uses at most twelve literal letters. This is a bound on a fixed period, not a finite reduction for changing modules or arbitrary varying words.

## Two structurally different infinite ribbons

Consider the bi-infinite periodic tournament. If some strong component is infinite, it is unbounded in at least one direction. Two distinct components cannot both be unbounded in that same direction: sufficiently separated vertices give forward arcs between them in both directions. Period translation therefore fixes the infinite component. It is unbounded both ways and contains every vertex, since each vertex can be reached from early members and can reach late members by long forward arcs.

Otherwise all SCCs are finite. Their condensation is linearly ordered by arc direction. Take the initial segment of vertices strictly before one component in that order. This cut differs from a physical half-line only within a finite margin: all vertices sufficiently far before the component beat it, and all sufficiently far after it lose to it. Translation by one period sends the cut to a strictly larger cut. Indeed two initial segments are comparable, and repeated positive shifts eventually contain any given vertex. Consecutive cuts delimit identical finite blocks; the size is $d$ since the net effect of shifting a half-line with finitely many edited sites by $d$ is $d$ extra sites. All arcs between the blocks are ordinal.

In $PW^kS$, take the first and last such cuts a fixed number of periods away from the corresponding physical ends. The first and last truncated blocks stabilize, including the literal end words: all modified end incidences across the chosen cuts are at distance greater than $b$. The resulting count is a constant times a fixed integer base to the power $k$, for every large $k$, by multiplication over the ordered blocks. It cannot be $1+B^k$ with $B>1$.

In the strongly connected case, fix two finite mutual-reachability paths for consecutive vertices in each period phase. The finitely many witnesses have a common maximum distance from their starting indices; translations then work for consecutive vertices outside a fixed physical margin in $PW^kS$. A single strong component contains at least $kd-O(1)$ vertices. The Moon–Busch floor gives $H(PW^kS)\ge5^{(kd-O(1))/3}$. Taking $k$-th roots proves the length gate.

## An exact negative-cycle test and next sieve

The strongly connected case is recognizable on just $d$ residue states. For each residue, include oriented displacements at gaps $1,\ldots,b$ from the word, and all forward displacements $b+1,\ldots,b+d$. Larger forward displacements can be replaced by the least positive one greater than $b$ in the same residue class, decreasing their weights without altering the phase walk.

A negative cycle in this finite graph lifts to a path with negative displacement from one instance of a phase to another. Repeat it to travel arbitrarily far left, using long forward arcs to join any desired start and finish. Conversely strong connectivity gives a path to a smaller vertex in the same phase; reducing large forward steps gives a negative closed walk in the residue graph, which contains a negative cycle. This is an exact finite criterion, not a sample definition of the infinite strong case.

For any candidate and any rotation $V$ of its period, another necessary test is $H(V^t)\le B^t$ for every $t\ge1$: use disjoint literal copies of $V^t$ inside the long periodic bulk and shuffle their Hamilton paths injectively. For span three the period matrix in rank 62 must also have eigenvalues $1,B$, by Cayley–Hamilton applied to the eventual scalar sequence. These give a finite candidate universe at radix 1000, but **the full twelve-letter eigenvalue sieve has not been exhausted**. Passing these tests would still not supply the two physical end readouts.

The [small verifier](verify_literal_period.py) checks the negative-cycle translation against literal reachability for span three through period four, and the exact integer cap $5^{12}\le 1000^3<5^{13}$.
