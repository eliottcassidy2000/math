# A target-directed unit edge in the single-arc-flip graph with HC=9

The graph operation here is reversing one physical arc. It gives
an edge in the ordinary arc-flip graph of tournament isomorphism
classes, not necessarily a line of the path-complement metagraph.

Take the one strong-five core with outgoing row masks
`(10,28,25,16,1)`. Change an apex's out-mask from 12 to 13,
i.e. reverse just its edge to old vertex zero. The marked
profiles in notation of [PROOF](PROOF.md) are
[
 (h,f,P,S,E,T;K)=(9,84,34,4,41,5;29)
               quadhbox{and}quad (9,84,45,5,30,4;31).
]
For positive integer x, use left and right transitive readers
both with address x, under the ordinary left and
reverse-complement right conventions. The contribution of the
two changes of start/end totals is +11x−11x=0. The contribution
of the two signed changes is +p(x)−p(x)=0. The exact count
pair therefore is
[
 84x^2+75x+18g(x)+9p(x)+{11,13}.
]
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
[
 84x^2+87x+18g(x)+7p(x)+{11,13,15}.
]
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
