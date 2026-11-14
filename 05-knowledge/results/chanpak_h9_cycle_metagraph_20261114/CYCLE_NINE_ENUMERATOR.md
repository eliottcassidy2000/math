# Exact size enumerators for two physical families at Hamilton cycle count nine

These are enumerators of the explicitly constructed classes, not
of the entire HCycle=9 stratum. They combine the
[diagonal cut decoder](h9_bilateral/DIAGONAL_EDGE.md) with the
[cycle-root counts](CYCLE_NECKLACES.md). All classes counted here
are rigid, so an ordinary class has exactly n! labeled copies.

## A diagonal family whose path counts distinguish every class

Let ℓ≥1 and n=2ℓ+6. For each address x in [2^(ℓ−1), 2^ℓ−1]
use the two diagonal masks and the separate three-mask window.
The count formula places all five outputs strictly in the same
interval (84x²,84(x+1)²); the two windows do not meet inside
that interval. Hence the 5·2^(ℓ−1) classes are distinct even
modulo converse and have pairwise distinct path counts.
Each has nine directed Hamilton cycles. Its nine unrooted
cycle decorations each have exactly n rooted closed cuts
because the tournament is rigid.

Let z mark order n. The ordinary class series, the ordinary
unrooted-cycle-marked series, and the rooted closed-cut series
of these selected representatives are, respectively,
$$
 A(z)=\frac{5z^8}{1-2z^2},\qquad 9A(z),\qquad 9zA'(z).
$$
The closed series counts line incidences *before* identifying
converse presentations inside an unmarked representative.

One can even sum all full-path-marked presentations of this
bank explicitly. For any x, sum its two and three path-count
formulas from the diagonal proof. The result is
$$
 420x^2+411x+90g(x)+39p(x)+63,
$$
where p(x)=x mod 2. On the length-ℓ block put M=2^(ℓ−1).
For ℓ≥2 the sums of x, x², and p are
$$
 M(3M-1)/2,\quad M(2M-1)(7M-1)/6,\quad M/2.
$$
The sum G_ℓ of g on the same block is M/2+3M²/4.
Indeed g(2y)+g(2y+1)=2g(y)+y+p(y), and G_2=4.
For the next lengths the parent block has even size, and
summing y+p(y) on a parent block of size m gives 3m²/2.
The recurrence proves the claimed expression.

Consequently the exact full-path-marked count is
980M³+54M²−8M for ℓ≥2, and 1023 for ℓ=1.
It has generating function
$$
 1023z^8+z^{10}\left(
  \frac{7840}{1-8z^2}+\frac{216}{1-4z^2}-\frac{16}{1-2z^2}
  \right).
$$
In particular the ordinary number of presentations does not
grow at the same exponential rate as the class count, and the
closed portion has yet another exact coefficient law.

## A larger independent-tape class census, without distinct-H claims

Use the same rigid five-vertex core and any of the same five
apex masks. Take independent transitive left and right tapes
of lengths ℓ,m≥2. Require that the first two left vertices
are beaten by the apex, and the final two right vertices beat
it. Otherwise paint each tape arbitrarily. Each graph is
strong and has nine Hamilton cycles because the extreme
singleton SCCs close every full path of the puncture.

The apex is the unique vertex whose deletion makes the graph
nonstrong. Deleting it leaves the ordinal puncture. If any
other vertex is deleted, a first left singleton and a last
right singleton with the required contacts survive. Every
Hamilton path of this smaller puncture begins and ends at
the surviving ordered extremes. Closing it through the apex
gives a spanning cycle, proving the smaller tournament strong.
Thus any isomorphism of the original graphs identifies the
apex, then the unique nonsingleton component and both
ordered tapes in its puncture. It determines all contacts:
the core is rigid, and each tape is transitive. This proves
rigidity and distinguishes all configurations, not just
their counts.

At fixed n=ℓ+m+6 there are n−9 allowed width splits and
(5·2^{ℓ-2}2^{m-2}) ordinary classes per split. Hence the
constructed subset has exactly
$$
  5(n-9)2^{n-10}\quad(n\ge10)
$$
ordinary classes and $5z^{10}/(1-2z)^2$ as its class series.
The same exact number survives converse merging. The unique
anti-isomorphism of the rigid core is (0,4,3,2,1) in the row labels;
it maps the complemented five masks (12,13,15,23,14) to
(19,18,2,4,3), none allowed in the selected family. A converse
isomorphism would have to identify the unique apex and the core,
so no two selected classes are converses and none is self-converse.
Their Hamilton path counts need not be distinct.
This is a structural lower bound of order n·2^n within a
single cycle-count fiber, complementary to the smaller
diagonal bank where every path count distinguishes a class.
It does not approach a coverage proof for all odd path counts.
