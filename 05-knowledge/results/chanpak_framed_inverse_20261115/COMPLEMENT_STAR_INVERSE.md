# An inverse with fixed positive undirected Hamilton count

The ordered tile chart is a genuine parameterization, not an identity
of the two Hamiltonian statistics. This construction uses graphs whose
undirected Hamilton path counts are **positive** and exactly known.
One inverse on their labelings gives a new unbounded family of
adjacent physical tournament counts.

## A readable geometry, with binary placements still free

Use the ascending frame. For the rooted remainder W, take a
transitive spine on m vertices 0,...,m−1, then a vertex u=m.
Let u beat a nonempty subset S of the first m−1 positions and
lose to the last spine vertex. Require 0 in S for strongness;
equivalently choose the shortest spine that holds the following
positive address:
```
e = sum_(i in S) 2^(m−2−i).
```
With backward tiles as edges, the smaller graph is a star on
r=|S| leaves and isolated vertices. Use the **complement chart**
instead: the graph is K_m with r edges missing at the distinguished
gap, whose degree is d=m−1−r. This is one ordinary graph geometry whose
order embeddings vary. For m>=3 it has exactly
`d(d+1)(m−2)!/2` undirected Hamilton paths up to reversal.
To count ordered paths, those putting the distinguished vertex
at an endpoint contribute 2d(m−2)!, and internal positions
contribute (m−2)d(d−1)(m−3)!. Divide by two for reversal.
Thus even the entire unlabelled graph, not just its scalar
Hamilton count, is unchanged when choosing the r positions.

The weighted even path count E at u is exactly e. After u, a
path through the spine is increasing; its first spine vertex i
must be in S. The remaining selected spine vertices are an even
subset of the suffix, yielding 2^(m−2−i) choices. Each
complement is transitive and contributes one. The
[one-apex word law](../chanpak_ternary_emptying_20261115/RUN_EMPTYING.md)
gives H(W)=1+p(2e), where p is the inherited binary surcharge.
Consequently a new vertex v beating only u gives the literal
strong count
```
N(e) = H(W)+2E = 1+3e+p(e)+(e mod 2).
```
Strongness follows from u→0 and the forward spine, returning
to u from its last vertex; the new vertex has inbound arcs
from the spine. The new frame is the old spineline followed
by v,u. Its backward tile graph is still the same star, now
with one additional isolated gap. Thus **both assembled
tournaments** also have the same unlabelled complement chart
whenever the m,r values coincide. For an arbitrary frame,
attaching an outlet immediately before the final root has
this same isolated-gap behavior: relocate the last free
column and put one isolated backward gap before it. If
the root is not final, appending at the end instead adds
the [last-labelled leaf](FRAMED_SIMPLE_GRAPH.md).

## The exact load is a pair interaction

Let J(e) sum omega_(j−i) over pairs of set bits i<j of e,
with omega_1=2 and omega_d=2^(d−1) for d>=2. The word p(e)
sums the same weights over the upper-1/lower-0 pairs.
For each fixed upper bit j>0 the sum of all its lower
weights is 2^j. Hence
```
p(e)+(e mod 2)=e−J(e),       N(e)=1+4e−J(e).
```
The graph retains the number of selected positions but
the load includes their **ordered separation**. For example
moving an isolated set bit one position up from j to j+1
increases N by at least 4 when j=0, or 3*2^j+2 if j>=1.
Indeed the unweighted term increases by 4*2^j; the
increase in interactions with lower bits is bounded by
the telescoping sum omega_(j+1)−omega_1, and interaction
with higher bits can only decrease. Thus even a monotone
move inside this geometry is too coarse to give +2 by itself.

## A balanced two-scale placement gives actual adjacent odds

Separate the lower bits x from the upper block 2^k y by at
least one unused bit. Write rho(x)=sum_(i a set bit of x) 2^(-i−1).
Then the cross term of J is exactly 2^k y rho(x), and
```
N(x+2^k y) = 2^k y(4−rho(x)) + 1+4x−J(x)−J(y).
```
Use (x,y)=(6,27) and (5,29). Their coefficients are both
783/8 because 27(4−3/8)=29(4−5/8). The intercepts are +1
and −1, respectively. Each address has six set bits, and
both upper words have length five. For every k>=4 this
gives the strong pair
```
783*2^(k−3)−1,   783*2^(k−3)+1.
```
Each pair has the **same unlabelled complement chart** on
either remainder or assembled tournament, with positive
and equal numbers of graph Hamilton paths. This is an
exact coordinated move missed by a single adjacent swap.
It does not preserve all rooted continuation coordinates
and the pairs are not a single consecutive interval across k.
The first pair is 1565 and 1567, at tournament order 12.

## Prescribed-target inversion within the geometry

Given an odd target N, e lies in [(N−1)/4,(N−1)/3].
Choose a power of two q of order sqrt(N), write e=q a+b,
0<=b<q, and let rev(b) be the q-width bit reversal.
The exact [digit split](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md)
reduces equality to
```
Y = L_b a + p(a) + delta (a mod 2);
Y = N−1−3b−p(b)−(b mod 2);
L_b = 4q−1−rev(b),   delta = 1−top_bit(b).
```
Reject Y<0. Since 0<=p(a)<=a, the only candidates are
between max(0,ceil((Y−1)/(L_b+1))) and floor(Y/L_b).
There are O(1) candidate a per b with this q. Test equality,
and, if a graph shape was specified, enforce its address
width and popcount. Recover the actual remainder by reading
the selected digits. This is an exact partial inverse,
O(sqrt(N) log N) arithmetic operations, not a universal selector. For instance target 63 decodes to
e=19 and (H,E)=(25,19); target 69 decodes to e=18 and (33,18).
Even 15 is absent in this particular channel.

The unrestricted one-row family already has all fixed
odd-compatible residues: (e mod 2M, p(e)/2 mod M) is independently
programmable by the [prefix theorem](../chanpak_absolute_steering_20261115/PREFIX_PROGRAMMING.md).
Nevertheless its presentation count through X is
X log(4/3)+O(sqrt X): apply the [binary row law](../chanpak_viability_shape_20261115/LOCAL_DENSITY.md)
with slope 3 and checksum coefficient 1; the bounded parity
term affects only the cutoff error in the bit-reversed Riemann
sum. Therefore its distinct absolute counts omit at least a
relative 1−2 log(4/3) of odds in lower density. Matching the
geometry and all residue tests is not sufficient for universal
absolute placement. To get beyond it one needs other graph
shapes or another independently selectable response.
