# Dyadic line obstructions and the entropy cost of an extra visit

**Status:** proved all-order consequences of the inherited odd-cycle
expansion. Finite direct controls in `reproduce.py` and `replay.txt`.
This does not settle line coverage with unrestricted sources or Chan–Pak.
No historical priority claim.

## 1. Which pair lines can keep a small dyadic response?

Let $H=P_1$, $P_j$ be ordered covers by exactly $j$ nonempty directed
paths, and $f(z)=\Pi_z$, so $Q=f(2)=P_2+2H$. Use the inherited
Rédei–Berge slot expansion: sum over disjoint directed odd-cycle packings,
including the empty one. A packing of $k$ cycles leaving $u$ vertices
unused contributes $2^k z^{u+k}$. Cycle rotations are not counted as
distinct, and singleton vertices count as unused, not directed cycles.
The exponent has parity $n$. In particular
$Q\equiv(-1)^n H\pmod 3$.

At $z=2^s$ the packing contributes $2^{su+(s+1)k}$, with
$u+k\equiv n\pmod2$. This gives the following simple all-order tests
for $s\ge1$:

* For odd $n\ge5$,
  $f(2^s)/2^{s+1}\equiv C_n \pmod{2^{2s}}$, where $C_n$ counts
  directed Hamiltonian cycles. A one-cycle packing with $u=0$
  has lowest cost $s+1$; the next possible $u$ for one cycle is
  2. Two cycles require odd $u\ge1$; three have cost at least $3s+3$.
  The empty term's cost $ns$ also exceeds the first band.
  The parity-only statement holds at $n=3$ as well.
* For even $n\ge6$,
  $f(2^s)/2^{2s+1}\equiv C_{n-1}+2D_n\pmod{2^{2s}}$, where
  $D_n$ counts unordered pairs of vertex-disjoint directed odd cycles
  spanning the tournament. Their costs are $2s+1$ and $2s+2$
  respectively, both admissible; one cycle with $u\ge3$ costs
  at least $4s+1$, two with $u\ge2$ at least $4s+2$, and
  three require $u\ge1$ and cost at least $4s+3$. Four or more
  also exceed the band; the empty term has cost at least $6s$.
  The parity-only statement holds at $n=4$ as well.

Taking $s=1$ gives $P_2\equiv-2H\pmod4$ for $n\ge2$, and
$P_2\equiv-2H\pmod8$ for even $n\ge4$. The next bit is the
cycle parity above. For example, in the triangle one-port line
$N=3H+P_2=H+Q$, a prime $N>3$ requires even source order.
If this route seeks a source with $H\equiv N-8\pmod{16}$, it must
have an odd number of $(n-1)$-cycles. Moreover $Q$ multiplies under
ordinal joins: restrict each of the two slot paths to each factor,
then concatenate the corresponding paths. Thus for $n\ge4$ even and
$v_2(Q)=3$, the source is either strong or one odd strong component
plus a singleton; each nontrivial strong component costs at least
two bits, and a singleton costs one, ruling out every other decomposition.
This is a restriction on source choice, not a prime obstruction for
all constructions.

Both lowest valuation bounds are sharp in strong families. In an odd
transitive order $n\ge3$, reverse only $(0,n-1)$. Every cycle uses
that reversal and there is one Hamiltonian cycle. In even order
$n\ge6$, reverse both $(0,n-1)$ and $(1,n-1)$; every cycle
uses one of them. Every single deletion other than vertex $n-1$
allows exactly one spanning cycle on the remainder (it must leave
(n-1) for the earliest remaining vertex, then increase), giving
$n-1$ odd cycles of size $n-1$. For even $n=4$, a sharp
non-strong example is $C_3\Rightarrow TT_1$.
It is not enough to infer the next bit from $H$: in order 5,
orientation codes 40 and 76 both have $H=15$, but $Q=108,120$
respectively, with $Q\bmod8=4,0$. Codes invert lexicographically
listed arcs $(i,j),i<j$ at bit positions set to one.

## 2. The number of outputs of an extra-visit reader is sublinear

A one-apex source $A$ is a tournament obtained by adding one arbitrary
vertex to a transitive spine. Let a **fixed** one-port template $K,x$
read it by homogeneous substitution at $x$. Decomposing every
Hamiltonian path by its $j$ visits to $A$ gives
$F(A)=H(K[x\gets A])=\sum_{j=1}^{r} a_j P_j(A)$,
where $a_j$ are nonnegative integers counting formal traversals with
ordered occurrences of $x$, and $a_r>0$. A cyclic port yields
$r\ge2$ by taking a triangle through $x$ and inserting the remaining
vertices into the resulting path on the graph with $x$ doubled.
Alternatively this highest coordinate is exactly one plus the port
matching number in the previous reader-retention session.
If $r\ge2$, the number of different outputs up to $X$ is
$$
  O_K(X^{\log_{r+1}2}).
$$

**Proof.** If no apex out-neighbor $i$ precedes any in-neighbor $j$
in the spine order, the source is transitive. Then $P_r(A)$ is the
number of surjections from the $n$ vertices to $r$ labelled slots
(one increasing path in each nonempty cell), so the transitive sources
provide at most $O_K(\log X)$ outputs.

Otherwise retain the interval from the first apex out-neighbor to the
last apex in-neighbor, and include the apex. This core is strong:
the apex reaches every vertex in the interval via its first vertex,
and every vertex returns via its last vertex. The prefix and suffix
are ordinal padding, of total size $t$, so each slot polynomial gains
a factor $z^t$, independent of the division into prefix and suffix.
For core order $m\ge3$, there are $2^{m-3}$ choices of the orientation
word after fixing the endpoints. The core contains an odd cycle from
apex to first vertex, then any even subset of the $m-3$ interior
vertices in increasing order, then last vertex and apex. The sum of
their cycle terms alone gives
$$
 \Pi_r(A)\ge r^t r\big[(r+1)^{m-3}+(r-1)^{m-3}\big].
$$
Every directed cycle uses the apex, and there is at most one on any
given vertex set, since on deleting the apex the path is increasing.
For $v\ge1$,
$$
 \Pi_v(\text{core})\le v^m+2v(v+1)^{m-1}\le3(v+1)^m .
$$
Consequently the ratio $\Pi_{r-1}(A)/\Pi_r(A)$ is bounded by a
constant depending on $r$ times
$(r/(r+1))^m ((r-1)/r)^t$. By a union bound over one empty slot,
$P_r(A)\ge\Pi_r(A)-r\Pi_{r-1}(A)\ge\Pi_r(A)/2$ outside
finitely many $m+t$. If $F(A)\le X$, this bounds
$(r+1)^{m-3} \ll_K X/r^t$.
At fixed $t$ summing $2^{m-3}$ through the permitted $m$ gives
$O_K((X/r^t)^{\log_{r+1}2})$; summing over $t$ is geometric. ∎

The triangle case $F=H+Q$ thus has at most $O(X^{\log_3 2})$ outputs
from **all** one-apex sources, not just simple repeated words or dyadic
responses. Finite banks of cyclic one-port readers remain sparse on
that source class. This is a reason to change the reservoir rather
than to search longer in the full one-apex grammar. An extra hypothesis
would be needed to deduce sparsity for two independent apex words:
the present $2^{m-3}$ presentation bound would be replaced by roughly
$4^m$, removing this counting obstruction for $r=2$.

## Provenance

Reproducible sample: sort paths of Markdown results in the tree at
`1b91742ce55a202b73c80d56561ba56d00b503a6` as listed by
`git ls-tree --full-tree -r --name-only`, filter prefix
`05-knowledge/results/`, then Python `random.Random(115203).sample(paths,14)`.
See `sample_paths.txt`. The useful transfers were the saturated Gold
poset barrier's forced prefix/core/suffix and the Collatz cube-theta
weight vs height capacity bound: identify the exact retained width
before charging a cost to its possible words. The cycle expansion is
inherited, not re-proved here. The direct replay tests 1,498 small
labelled orientations, reports exact witnessing core bits, and does not
assume the expansion to compute its data.

### Ordinal products do not remove the entropy obstruction

The same sparsity extends to a substantially larger family: allow
arbitrarily many one-apex strong components in an ordinal source, plus
arbitrary transitive padding. For every fixed cyclic reader of rank
$r\ge2$ and every $\varepsilon>0$ its output count is
$O_{K,\varepsilon}(X^{\log_{r+1}2+\varepsilon})$.

Every component slot polynomial multiplies in an ordinal join, so the
readout depends only on the multiset of strong components, not their
ordering. Use one possible factor per apex presentation (allowing
duplicates only overcounts), including a singleton for padding.
Give factor $B$ weight $W_B=\Pi_r(B)>1$. For
$\sigma>\log_{r+1}2$ the strong component candidate series converges:
the core of order $m$ has at most $2^{m-3}$ presentations, and each
weight is at least $r(r+1)^{m-3}$. Thus the multiset sum
$\prod_B(1-W_B^{-\sigma})^{-1}$ is finite, since its logarithm is
bounded by a constant times $\sum_B W_B^{-\sigma}$. The number of
multisets of product weight at most $Y$ is therefore $O_\sigma(Y^\sigma)$.

To apply this count, $P_r$ is still at least $\Pi_r/2$ outside
finitely many source orders. For each nonempty component all monomials
in its slot polynomial have positive exponent and nonnegative
coefficient, so the ratio $\Pi_{r-1}/\Pi_r$ is at most $(r-1)/r<1$.
A source with many factors thus has small product ratio. One with
boundedly many factors and sufficiently large order has a large strong
one-apex factor, whose ratio tends to zero uniformly by the bound in
Section 2. Hence $F\le X$ again forces $\Pi_r\ll_K X$, completing
the proof.

For the pair line $N=3H+P_2$, the simpler argument needs no comparison:
$N=Q+H\ge Q$, and $Q$ itself is the product weight. In particular,
even arbitrary recombination and repetition of one-apex SCC products
cannot turn the triangle line into a density-one proof. To exploit
ordinal multiplicativity, an actual supply of richer strong components
is needed. This is the tournament counterpart of a product Euler
series failing its dimension budget, not a limitation on mixing cores
inside one SCC.

### Exact order bands do not overlap at distance two (triangle reader)

For a strong normalized one-apex source on $m$ vertices, write its
spine word $b_0\cdots b_L$ ($L=m-2$), with 1 meaning the apex beats
the vertex. Thus $b_0=1,b_L=0$. Put
$J_d=\#\{i:b_i=1,b_{i+d}=0\}$. Choosing first and last spine vertices
of an odd cycle, then its even-sized subset of intervening vertices,
gives the exact gap formula
$$
F=H+Q=1+2^m+\sum_{d=1}^L J_d w_d,\qquad
w_d=2^{\max(1,d-1)}+2^{m-d-1}(3^{d-1}+1).
$$
The two summands are respectively its contributions to $H,Q$.
Both are nondecreasing in $d$: for the second, the increase
from $d$ to $d+1$ is $2^{m-d-2}(3^{d-1}-1)$.

There is an exact minimizing word: choose bit 0 at interior position
$i$ if $i\le L-i$, bit 1 otherwise. Every interior vertex must
contribute an inversion with one endpoint, at cost at least
$\min(w_i,w_{L-i})$; this choice attains all the separate lower
bounds and has no interior inversions. Thus the minimum is
$1+2^m+w_L+\sum_{i=1}^{L-1}\min(w_i,w_{L-i})$.
For the maximum, swapping adjacent 01 to 10 while holding the endpoint
bits fixed increases the sum: it gives one new inversion and pushes
each changed existing inversion farther apart. At fixed first block
size $a$, the extremal word is therefore $1^a0^b$, $a+b=m-1$.
Its gap inventory is $J_d=\min(d,a,b,m-1-d)$, simultaneously
maximized by $|a-b|\le1$.

For these two-block sources the apex cycles can choose independent
nonempty subsets of each of the two blocks, with an even combined size.
Hence
$$
F=2\cdot3^{m-1}+7\cdot2^{m-1}-2(3^a2^b+2^a3^b)-3(2^a+2^b)+5 .
$$
Since $2((3/2)^a+(3/2)^b)\ge7.5$ for $a+b\ge3$, this is
strictly below $2\cdot3^{m-1}$ for $m\ge4$; at $m=3$ the
same inequality reads $15<18$. But at order $m+2$ the endpoint
pair by itself contributes strictly more than $2\cdot3^{m-1}$
to $F$, from the $3^{L-1}$ part of its weight. Thus strong
one-apex sources at orders two apart have disjoint absolute
triangle-reader bands. In particular the even-order bands are
strictly separated, while odd-order sources cannot give primes
above 3. This does not forbid filling the gaps by adding ordinal
padding or changing SCC shape. The accompanying band replay uses
the explicit $J_d$ formula and cross-checks against the independent
subset path DP for all 63 normalized words through order 8; it
enumerates every normalized word through order 21.

### Wider source theorem: bounded feedback and a high-visit reader

There is a more general obstruction when the reader's visit rank is
large relative to the source's allowed feedback. Fix $k\ge1$ and let
every strong component of the ordinal source have a transitive
deletion set of size at most $k$ (meaning at most $k$ vertices need
be deleted to leave a transitive tournament). For a fixed rank-$r$
cyclic reader, with $r\ge 2^k$, the number of outputs at most $X$
is $O_{K,k,\varepsilon}(X^{k\log_{r+1}2+\varepsilon})$ for
every $\varepsilon>0$. In particular rank 4 already suppresses
the entire two-apex-per-component class in density.

**Long return via a short exceptional set.** Let $B$ be one strong
component with an increasing spine of size $m$ and $a\le k$
exceptional vertices. Take a shortest path from the last spine vertex
to the first. A run of three spine vertices in this path has a
shortcut, so the path uses at most $2(a+1)$ spine vertices. Let it
use $b$ spine vertices and $c$ exceptional ones. Join its two endpoints
back along the increasing spine through any subset of the other
$m-b$ vertices. Each choice gives a different cycle; choose the
parity of the subset size so that the cycle is odd. These cycles
alone contribute
$$
r^{a-c+1}\left[(r+1)^{m-b}
                 +(-1)^e(r-1)^{m-b}\right]
$$
to $\Pi_r(B)$, for some $e\in\{0,1\}$. If $m-b\ge1$, this
is bounded below by a positive constant depending only on $r,k$
times $(r+1)^m$, since $b\le 2(k+1)$. All other orders are bounded
in terms of $k$, and the empty packing supplies the same bound
after decreasing the constant. The case where the return has
just one edge is included: if an empty subset would give a
2-cycle, that subset is even and already excluded by odd parity.

**No fix is needed for empty slots.** More generally, for any
tournament of order $n$, a nonempty ordered $j$-cover can be split
at any of its $n-j$ path arcs, placing the suffix in a new
$(j+1)$st slot. A resulting cover has at most $j$ preimages,
each obtained by appending its last path to one of its first $j$
paths. Thus $j P_{j+1}\ge(n-j)P_j$. For fixed $r$, summing the
lower slots with their binomial placement multiplicities gives
$\Pi_r\le(1+O_r(1/n))P_r$ as $n\to\infty$, for all tournaments.
In particular the comparison used above extends beyond one-apex
sources.

There are at most $2^{am+\binom a2}$ spine presentations with
$a$ exceptional vertices; this is at most a constant times
$2^{km}$. Hence the component weight series
$\sum_B\Pi_r(B)^{-\sigma}$ converges at each
$\sigma>k\log_{r+1}2$, by the long-return lower bound, and the
same ordinal multiset product argument applies. Rank $r\ge 2^k$
ensures the exponent can be chosen below 1. This proof identifies
a quantitative reason to try **two** apexes with a **two**-visit
reader: their $4^m$ presentation budget is no longer defeated
by the $3^m$ long-return weight alone. It does not claim such a
reader gives full coverage.
