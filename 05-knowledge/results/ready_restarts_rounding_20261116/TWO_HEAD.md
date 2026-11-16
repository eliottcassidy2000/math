# One two-head template hits both restart scales and prescribed primes

This [supplement](README.md) derives a physical family, not a
universal selector. Use a terminal block \(E\) with outgoing masks
\((6,20,8,3,13)\), root \(x=3\), second head \(y=4\).
Prepend a transitive word with ordinary binary value \(X\),
empty for \(X=0\). Every prefix vertex beats every terminal
vertex, except \(y\) beats each prefix vertex whose bit is 1.
Write \(f=f(X)\) for one-apex paths in that prefix word and
\(p=X\bmod2\). The complete profile at \(x\) is
\[
g=f+2X+2,\qquad h=5f+9X+8+p,\qquad
r=3f+4X+3+p,\qquad s=g+r.
\]

**Why this is exact.** In the odd-cycle packing expansion, a cycle
meeting the prefix must use \(y\), and there can be only one such
cycle. No use of the prefix contributes \(H(E)\). Cycles otherwise
using just a prefix segment and \(y\) contribute
\((f-1)H(E-y)\), since the complement splits ordinally.
Mixed cycles match an increasing prefix path starting at a 1-contact
to a nonempty terminal path outside \(y\) ending at a vertex that
beats \(y\). Their path orders must have equal parity. The numbers
of even and odd prefix paths are \((X-p)/2,(X+p)/2\); for each
possible first bit the subsequent choices balance parity unless
that first bit is at the last position. Let \(B_e,B_o\) be the
terminal path sums of \(H(E-y-V(Q))\) in these two parities.
The full mixed contribution is \(X(B_e+B_o)+p(B_o-B_e)\).
All four terminal constants below come directly from subset DP:

| Terminal | \(H\) | \(H\) after deleting \(y\) | \(B_e\) | \(B_o\) |
|---|---:|---:|---:|---:|
| \(E-x\) | 3 | 1 | 1 | 1 |
| \(E\) | 13 | 5 | 4 | 5 |
| \(E\) after one near-sink append at \(x\) | 25 | 11 | 7 | 10 |
| \(E\) after two near-sink appends at \(x\) | 49 | 23 | 13 | 20 |

The first two rows give \(g,h\); the succeeding increments are
\(2r,4r\), giving \(s=g+r\) by the append law. This verifies the
profile for every prefix length.

Consequently after \(k\) appends, with \(d=2^k\),
\[
H_k=(8d+1)X+(6d-1)f(X)+(6d+2)+(2d-1)p.
\]
The [published binary inverse](../prescribed_prime_30000_20261116/README.md)
solves the first, second and fourth terms together as
\(aX+bf(X)+c(X\bmod2)\), stripping the last bit and carrying a
finite change of coefficients. Here the \(k\)-scan is finite
for any target, because the count is at least \(6\cdot2^k+2\).
This is target selection without guessing an approximate prefix.
Full literal subset DP replays the following:

| Prescribed count | \(k\) | A selected \(X\) |
|---:|---:|---:|
| 13 | 0 | 0 |
| 41 | 0 | 2 |
| 125 | 0 | 8 |
| 10009 | 2 | 214 |
| 30011 | 0 | 2673 |
| 30089 | 0 | 2625 |
| 30091 | 0 | 2332 |

Those are individual exact counts, not the much smaller edit-pair
tournaments. The same inverse certifies *no* representation in
this fixed template of 10007, 30013, 30021 or 30027. It is not a
claim about the full tournament spectrum.

**Why this cannot itself give a synchronized +2 pair.** Equality of
the tail slope is equality of \(r\). Since \(f(X)\) is odd, comparing
\(r=3f+4X+3+p\) modulo 2 first forces the same \(p\).
Then \(3\Delta f+4\Delta X=0\) forces
\(\Delta X=6u,\ \Delta f=-8u\) for an integer \(u\), and hence
\(\Delta h=14u\). Thus any two outputs with matched slopes in
this template differ by a multiple of 14, not by 2.
Even when both requested primes have witnesses (30089/30091),
the matched-slope requirement forces a change of physical context.
Run python3 two_head.py for the terminals, the inverse, and the graph
replays.
