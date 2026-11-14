# A target-directed unit edge in the single-arc-flip graph with HC=9

The graph operation here is reversing one physical arc. It gives
an edge in the ordinary arc-flip graph of tournament isomorphism
classes, not necessarily a line of the path-complement metagraph.

Take the one strong-five core with outgoing row masks
`(10,28,25,16,1)`. Change an apex's out-mask from 12 to 13,
i.e. reverse just its edge to old vertex zero. The marked
profiles in notation of [PROOF](PROOF.md) are
$$
 (h,f,P,S,E,T;K)=(9,84,34,4,41,5;29)
               \quad\text{and}\quad (9,84,45,5,30,4;31).
$$
For positive integer x, use left and right transitive readers
both with address x, under the ordinary left and
reverse-complement right conventions. The contribution of the
two changes of start/end totals is +11x−11x=0. The contribution
of the two signed changes is +p(x)−p(x)=0. The exact count
pair therefore is
$$
 84x^2+75x+18g(x)+9p(x)+{11,13}.
$$
The two graphs have different path counts, so are not isomorphic.
With normalized words their extreme singleton SCCs close
every path of the puncture; both graphs have exactly nine
Hamiltonian cycles at every x. Thus this is an explicit
unbounded family of actual single-arc edges along which the
Hamiltonian path count changes by two while the cycle count
is constant. Every output is coprime to 3. At x=78 one
obtains 517727=7·73961 and 517729 prime; trial division
certifies both prime numbers, and literal subset DP recounts
both Hamiltonian path and cycle numbers on both 20-vertex graphs.

The inverse requires no height scan and no factorization.
The transitive recurrences give 1≤g(x)≤x+1. To check the
upper bound use induction; for x=2v and v even the bound
is g(v)+v≤2v+1, and for v odd, g(v)=g(floor(v/2))≤v
so g(v)+v+1≤2v+1. Odd x follows immediately. Thus both
counts lie strictly between 84x² and 84(x+1)², since their
excess over 84x² is at most 93x+40. From a prescribed
target N the unique possible address is
x=floor(sqrt(floor(N/84))); evaluate the explicit count
and check both arms. The family yields Θ(sqrt(X))
distinct counts through X, more than a finite bank of
fixed dyadic-height rails, but is still not cofinite.

## A three-odd window under the same diagonal guard

On the same strong-five core the three actual masks 15,23,14
have (P+E,S+T)=(87,7) and K=29,31,33 respectively.
With the same arbitrary x on both sides, they give
$$
 84x^2+87x+18g(x)+7p(x)+{11,13,15}.
$$
The first and last masks differ in just one apex arc; their
path-count change is four, and the middle odd is supplied
by a different mask. All three cycle counts are again nine.
The same square-root inverse works: the excess over 84x²
is bounded above by 105x+40. Choosing a matching offset
constructs the correct physical cut directly from N.

These are protected adjacent profiles in one concrete
cycle-count stratum, alongside a precise source of sparse
absolute count families. The formulas do not prove that
arbitrary odd targets have such a diagonal witness.

## Exact class multiplicities in the diagonal bank

All five masks in the two windows give rigid tournaments. If
the width of x is ℓ≥1 then n=2ℓ+6. A left-tape vertex has
outdegree at least ℓ+5, a right-tape vertex at most ℓ, and a
core vertex is in [ℓ+1,ℓ+4] (its internal score is 1, 2, or 3).
The two tape contact counts of the apex sum to ℓ, so its
score is ℓ+|U| in [ℓ+2,ℓ+4]. The left and right tape
vertex sets are therefore determined by scores, and each is
fixed pointwise because it is transitive. The apex is the
only remaining vertex that beats the first left vertex.
It is fixed too, leaving an automorphism of the strong-five
core. Its internal scores (2,3,3,1,1) have multiplicity at
most two, and a tournament has no automorphism of order two
(a transposed pair would have to reverse its arc). Hence
every score class in the core is fixed pointwise.

Different x have disjoint square-root intervals used in the
inverse. At a fixed x the smallest count in the triple is
larger than the largest in the pair by
12x−2p(x)−2>0. Thus at order 2ℓ+6 the bank gives exactly
5·2^(ℓ−1) different constructed rigid classes, all with HC=9
and all with distinct H. Each has exactly n! labelings.
This is a lower bound on the whole HC=9 stratum, not an
enumeration of all tournaments in that stratum.

## An exponentially larger rigid cycle stratum with independent tapes

In the same construction, let the left and right widths ℓ,m
be independent and at least two. Require both addresses to
begin in binary with 11. This says the first two left vertices
are beaten by the apex, and the last two right vertices beat
it. The contact bits on each tape other than those extreme pairs
may be freely chosen. Keep any one of the five distinct core
masks above. All outputs again have exactly nine Hamiltonian
cycles, regardless of their path counts.

Now the apex is intrinsically recoverable: it is the only
vertex whose deletion yields a non-strong tournament.
Deleting the apex gives the ordinal sum of the two tapes and
the strong-five core, hence is not strong. After deletion
of any other vertex, the remaining puncture has a Hamiltonian
path and still begins and ends in surviving singleton
components with the requisite contacts. Closing that path
through the apex gives a spanning directed cycle, hence the
deletion tournament is strong. The two guard contacts per
side are exactly what make this argument independent of which
vertex is deleted.

After identifying the apex, the invariant order of the
puncture's strongly connected components recovers the
two widths and the ordered singleton vertices, so it also
recovers their contact words. The strong-five core is rigid,
as above, and its core mask is uniquely recovered. Thus at
order n≥10 this bank gives exactly
$$
 5sum_{ell+m=n-6; ell,mge 2}2^{ell-2}2^{m-2}
 = 5(n-9)2^{n-10}
$$
distinct constructed rigid isomorphism classes with HC=9.
The same number are distinct after identifying converse
classes: the unique core anti-isomorphism is
(0,4,3,2,1) in image notation (verify the row masks; uniqueness
is from rigidity), and complementing + relabelling the five
core masks sends them to 19,18,2,4,3, all outside the chosen
set. Each class therefore has a distinct converse missing
from this bank. This class count is a lower bound on the
full HC=9 stratum, not a claim that all its path counts differ.
The diagonal subbank above is where the path-count
separation is also guaranteed.
