# Exact line inversion by deleting singleton components

**Status:** all-order deductions from SCC factorization and Busch's bound; the
small-budget replay checks the implementation, not Chan–Pak coverage.

Write `h=H(T)`, `p=P₂(T)` (ordered pairs of nonempty paths), and
`q=p+2h`. For a linear reader `N=Uh+Vp` with positive integer U,V,
set `W=U−2V`, so that `N=Vq+Wh`. For example the physical
singleton hyperbolic reader is `N=6q−h=11h+6p`.

## 1. An exact vertical fiber, and a cheaper inverse

Delete every singleton SCC, leaving the ordered join S of the nontrivial
strong components. With t deleted singleton components one has
`h(T)=h(S), q(T)=2^t q(S)`: a path of the join is a concatenation,
and the same is true independently in each of the two labelled slots.
The empty core has the formal bookkeeping convention `h=q=1`, and must
have t≥1 to represent a nonempty tournament. Every lift is constructive:
insert t sources or sinks into the component order.

By [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3),
each nontrivial strong component C satisfies `H(C)≥3^(|C|/3)`.
Indeed the exact minima by order modulo three (repo THM-4744) are
`3·5^(k−1), 5^k, 9·5^(k−1)` at orders 3k, 3k+1, 3k+2, and
`5>3`; the smallest initial comparison in each class suffices.
Thus S has at most `floor(3 log_3 h)` vertices, irrespective of |T|.

For fixed odd h, compute the finite set of primitive records
`(q(S), |S|)` by ordered multiplicative convolution of strong profile
records whose h-coordinate divides h. Only component factors of h can be
used. For a prescribed line and h, do **not** enumerate orders up to log₂ N:
require that `D=N−Wh` be positive and divisible by V, and set
`u=D/V`. A primitive record is legal iff `u/q(S)` is an integral power
of two (and the empty-core convention is met). For exactly order n, also
require `n=|S|+log₂(u/q(S))`. Equivalently check the odd part of u
and its two-adic valuation; restoring singleton components cannot correct
a wrong odd part. The reduction replaces a target-sized search by the
prime factorization of h and a bounded component catalogue. No search
bound in h is removed.

An enumeration refinement comes for free: if there are k nontrivial SCCs
and t singleton SCCs, the number of classes obtained from each ordered
nontrivial core S is precisely `binom(t+k,k)`. An isomorphism preserves
the ordered SCC list; the choice of which of k+1 gaps receive the t
singletons is the only freedom. This lifts the repo's existing SCC class
series to a fiber refined by q, including the exact size, before
projection onto a target line. It does not use a converse quotient.

For example, the entire H=3 fiber is `q=12·2^t`, t≥0, and the
H=5 fiber is `q=32·2^t`, t≥0. This follows from the strong floors and
the direct four-vertex classification. For the reader `N=6q−h` the
H=5 fiber therefore realizes exactly `192·2^t−5`. Meeting its dyadic
parity or even all congruences at a fixed finite modulus is not a
substitute for the odd-part test.

## 2. Almost every absolute target requires a growing cyclic core

Let `s` be the number of vertices in the source after deleting its
singleton SCCs. For any fixed finite bank of positive linear readers,
the number of distinct outputs ≤X with `s≤B` is at most
`O((B+1)·2^(B(B−1)/2)·log X)`. Indeed there are at most
`sum_{s≤B} 2^(s(s−1)/2)` reduced sources, even counting labelled
tournaments, and t has only O(log X) possible values:
`Vp≥V(2^|T|−2)`. The factor 2^t in q is a clock, not fresh
independent freedom.

For a fixed finite bank of two-source readers
`N=a h h'+b h p'+c p h'+d p p'` with a,b,c,d positive integers,
the same reasoning with total retained size s gives the bound
`O((B+1)^2·2^(B(B−1)/2)·(log X)^2)`. At a division s=i+j,
the number of reduced profile pairs is at most
`2^(i(i−1)/2+j(j−1)/2)≤2^(s(s−1)/2)`. Both paddings are bounded
by O(log X), using c p and b p' respectively. In particular the bound
holds for every fixed two-visit kernel with a,d>0 and bc>ad. The
constants can be multiplied by the size of an independently chosen bank;
with a subpolynomial bank the following threshold is unchanged.

For every fixed `epsilon>0` less than sqrt(2), set
`B=(sqrt(2)−epsilon) sqrt(log₂ X)`. The outputs of these restricted
families have size `O(X^(1−delta))` for some `delta>0` (absorb the
polynomial logarithmic factors). Hence almost every odd target needs
strictly more than B nontrivial source vertices in any representation
from the bank, even when the total source order is arbitrarily large.
The same power saving is smaller than the number of primes up to X.
Using Busch, the same conclusion holds when the source H, or the product
of the two source counts, is bounded by `3^(B/3)`.

The threshold is a one-sided obstruction, not an estimate of the true
minimum cyclic-core size. It says why high-order checks that get their
order primarily from transitive padding are a poor test of universal
placement. A candidate proof now has a specific obligation: supply
growing cyclic cores on demand, or use a growing/changed reader whose
parameters are not supplied by this fixed-bank projection.

## Source of the analogy

The seeded old-work sample in OLD_SAMPLE.txt includes the selected
Collatz carry-defect certificate, the nested cubic square-sum contact
chain, and the exact no-pole Cauchy transform. Here the useful common
move is to separate the unbounded clock from the finite input: derive
an exact invariant after deleting the clock before doing any search.
The Collatz and algebraic results are inspiration, not dependencies.

## 3. Top order layers freeze: a triangle-budget reduction

There is also an exact regime in which h need not be fixed beforehand.
Write `τ(T)` for the number of directed triangles. The empty packing
and each one-triangle packing in the positive [Rédei–Berge expansion](https://arxiv.org/abs/2307.05569)
give `q/2^n ≥ 1+τ(T)/2`. A direct two-slot proof of the same inequality:
each of the 2^n assignments supports a pair of paths. In the positive
cycle expansion, each directed triangle has two orientations of its
path-count contribution and a common slot for its three vertices, so
its contribution to q is `2·2^(n−2)=2^(n−1)`; distinct triangle
packings are separate positive summands. An arbitrary assignment alone
would not justify adding the triangle contributions.

A nontrivial strong component of size m has at least m−2 directed
triangles (Moon's strong-deletion chain, adding a triangle through the
new vertex at each extension; see repo THM-4744). For the last step,
a new vertex in a strong tournament has nonempty O,I and some arc
O→I, since otherwise no vertex reachable from it could leave O.
All triangles are within SCCs.

Fix V,U>0 and a normalized upper radius K. At order n≥3, if
`N≤KV2^n`, then from `p≥(n−1)h`,
`p≥(n−1)q/(n+1)`. It follows that

    sum_{nontrivial SCCs C} (|C|−2) ≤ τ ≤ floor(2[K(n+1)/(n−1)−1]).

In particular any fixed-width *logarithmic* top shelf is exactly
covered or not covered by a fixed finite catalogue, even as n→∞.
One needs only to enumerate nontrivial SCCs with total cost at most
the displayed bound, then use
`N=V·2^n[q(S)/2^|S|]+W·h(S)` with `|S|≤n`.
This is a finite list of exact affine rays in 2^n, not a dense lattice.
It includes every admissible tournament for that shelf, not just a
selected construction. Outside the shelf it makes no coverage claim.

There is a useful sharper count by cost, rather than total retained
vertices. For cost j=∑(m_i−2), ordered costs are compositions of j.
Putting `g(w)=(w+2)(w+1)/2`, one has
`g(u+v)−g(u)−g(v)=uv−1≥0` for positive integers u,v.
So a cost bound J has at most `2^(J²/2+O(J))` reduced profile
records by counting labelled orientations separately in each component.
If X is large and n≥log₂(X/V)−(1/2)log₂ log₂ X+c, all outputs ≤X
have J≤(2^(1−c)+o(1)) sqrt(log₂ X). For any fixed c>1/2 the
preceding count with O(log X) paddings gives a power-saving upper
bound on the number of outputs ≤X in this high-order regime.
For the target-relative formulation
`n≥log₂(N/V)−(1/2)log₂ log₂ N+c`, the same proof applies in
each dyadic target band or directly by monotonicity of the J bound
for large N. Thus for almost every target N, a fixed linear reader
needs to leave more than roughly half a logarithm-of-a-logarithm
of headroom below its usual largest allowed order. This does not
preclude a witness from a smaller order.
