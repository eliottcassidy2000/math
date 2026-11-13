# Two opposed contacts give absolute +2 corridors, but no ready return

The [one-contact theorem](NEAR_TWIN_GATE.md) rules out a reusable +2
move of that shape. Adding a second contact of the opposite sign changes
the answer. Take the tournament with outgoing rows
$(22,24,18,5,8)$, on vertices $0,\dots,4$, and reverse $0\to1$
to get $(20,25,18,5,8)$. The two end vertices of this switched arc
have the same incidence to vertex 4; their two disagreeing contacts
are $0\to2\to1$ and $1\to3\to0$.

Replace vertex 4 uniformly by any nonempty tournament $U$, using
its old directions on every vertex of $U$. Let $Q_U(t)$ count
ordered covers by $t$ paths with empty slots allowed, and $F_j(U)$
count ordered covers by $j$ nonempty paths. Then the two strong
tournaments have exact counts
$$
 H_-(U)=6Q_U(2)-H(U),\qquad H_+(U)=6Q_U(2)+H(U).
$$
In particular any transitive $U=TT_L$ gives strong realizers of both
$6\cdot2^L-1$ and $6\cdot2^L+1$, differing by exactly two at
arbitrarily large absolute counts. A single literal edge change, not a
choice of congruence class, selects the sign. Both outputs remain
strong because both quotient tournaments are strong and every vertex
in a uniform substituted module can follow external quotient paths.

## Exact interface calculation

Delete the contents of $U$ in a path cover, retaining the placeholder
for each maximal consecutive run. The remaining itinerary uses the four
other vertices once, and the ordered runs are filled by $F_j(U)$
in exactly one way per ordered cover of $U$.
There are no adjacent placeholders in a slot. For a Hamilton path
there can be at most two placeholders, since the root has only one
out-neighbor in the quotient. The exact coefficients of $F_j(U)$,
indexed by the number of ordered slots outside, are:

| Slots | lower chart | upper chart |
|---|---|---|
| 1 | $11,6$ | $13,6$ |
| 2 | $100,98,30$ | $104,100,28$ |

These are finite itinerary counts, not fitted output values.
The standalone recurrence enumerates each choice of subset in the first
slot, keeps its endpoint and number of placeholders, then transfers to
the complementary subset in the next slot. The one-slot row gives the
displayed law using $F_1=H,\ F_2=Q(2)-2H$.

The second row gives the exact **change in the next two-slot response**:
$$
 \Delta Q_{\rm out}(2)=4H(U)+2F_2(U)-2F_3(U)
                    =2\bigl(4Q_U(2)-Q_U(3)-3H(U)\bigr).
$$
It is strictly positive for $|U|\le3$ and strictly negative for
$|U|\ge4$; it never vanishes. For $n=|U|\ge5$, split a path in
an ordered two-cover at any of the $n-2$ internal edges. Each ordered
three-cover has at most two merges of consecutive slots, so
$2F_3\ge(n-2)F_2$, with strict inequality: for any arc $i\to j$,
the three-cover whose first slots are the singleton $j$, then the
singleton $i$, has no first merge; the rest has a Hamilton path.
Also $F_2\ge(n-1)H$ by splitting a Hamilton path. Therefore
$F_3-F_2> (n-4)(n-1)H/2\ge2H$.
At order 4 the odd-cycle packing identity gives
$Q_U(t)=t^4+2c_3t^2$ with $c_3\le2$; the expression inside
parentheses is $-20+8c_3<0$.
At orders 1, 2 and 3 the corresponding expressions are $2,4,2+4c_3$.
This threshold is a discrete loss of interface equality, not a
heuristic about growth. For a transitive module the discrepancy is
$8\cdot2^L-2\cdot3^L-6$.

In particular feeding the outputs into another operation that requires
identical $Q(2)$ cannot treat these signs as exchangeable. A fixed
core with one variable transitive run is also too sparse for absolute
coverage. The move supplies exact adjacent odd targets at unbounded
scale, but does not by itself fill gaps between the centers.
