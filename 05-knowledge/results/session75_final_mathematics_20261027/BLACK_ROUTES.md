# Black routes with named endpoints and a strong complementary core

## Conventions and dependencies

For a directed Hamiltonian path `P`, let `R_P` reverse just its consecutive
arcs. The native move `K_P` instead retains those arcs and reverses every
other pair, so `K_P(T)=(R_P(T))^op`. A path is blue exactly when its positional
reversal is an anti-automorphism of the current labelled tournament; otherwise
it is black. A blue path necessarily has endpoint outdegree sum `n-1`.
The converse implication is not used.

Two dependencies are explicit. Thomassen's
[Hamiltonian-connected tournaments](https://doi.org/10.1016/0095-8956(80)90061-1)
states that every arc of a 3-strong tournament belongs to a Hamiltonian cycle.
Here 3-strong means that deleting at most two vertices preserves strong
connectivity. The primary abstract was checked; its full proof is cited,
not reproduced. The repository's corrected pure-blue classification says
that the strong pure-blue tournaments are `K1,C3,R5`; this is an inherited
classification input, not a new classification proved by the argument below.
See the corrected black-refinement discussion in
[THM-4811](https://github.com/eliottcassidy2000/math/blob/3858d2afca72fb97a46a3601a6ae4b1e6e99a812/01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md).

## A black Hamiltonian path at every prescribed endpoint

In a 3-strong tournament fix `v`. Suppose every Hamiltonian path ending at
`v` is blue. For every out-neighbour `u` of `v`, cut a Hamiltonian cycle through
`v -> u`; its resulting path starts at `u` and ends at `v`. Thus some
anti-automorphism sends `v` to each such `u`.

Let `G=Aut(T)`, `O=Gv`, and choose one anti-automorphism `a`. All anti-images
of `v` form the orbit `O'=Ga(v)=a(O)`, of the same size as `O`, because
`aGa^-1=G`. Hence `N+(v) subset O'`. If `O' != O`, the orbits are disjoint
and `v` has no out-neighbour in `O`. The induced tournament on `O` is
vertex-transitive, so all its internal outdegrees are zero: `|O|=1`.
Then `|O'|=1` and `d+(v)<=1`, contradicting 3-strongness.

Therefore `O'=O`. For any `x=g(v)` in `O`, `N+(x)=g(N+(v)) subset O`.
Strong connectivity forces `O` to be all vertices. Transitivity now transports
the assumption to every endpoint, making every Hamiltonian path blue.
The inherited pure-blue classification contradicts 3-strongness. Thus there
is a black path ending at each `v`; converse gives the initial-endpoint version.

## An exact quarter-semidegree region

Put `delta=min_v min(d+(v),d-(v))`. If deletion of `s` vertices produces a
cut `A -> B`, average indegrees in `A` and outdegrees in `B` to obtain

`delta <= (|A|-1)/2+s`, `delta <= (|B|-1)/2+s`.

Consequently `s >= (4 delta-n+2)/3`. Also every cut with smaller shore of
size `k` has directed capacity at least

`f(k)=k delta-k(k-1)/2`.

Suppose `n>=7` and `delta>=floor(n/4)+2`. Concavity of `f` on
`1<=k<=floor(n/2)` shows every cut has capacity at least `delta`, and every
nonsingleton/non-cosingleton cut has capacity strictly greater than `delta`.
The separator inequality gives 3-strongness.

Fix `v` with indegree `delta`. Choose a black Hamiltonian path ending at `v`.
Work first in the `R` gauge. Every cut changes by

`c_(R_P T)(A)-c_T(A)=1_A(last P)-1_A(first P)`.

Thus all whole-tournament states before move `delta` remain strong.
The deleted tournament `T-v` has minimum semidegree at least `delta-1`,
and `delta-1>floor((n-1)/2)/2`. The same concave-cut bound gives its cut
minimum at least `delta-1`. After the first reversal that minimum remains
at least `delta-2>=1`, so the first deleted core is strong.

Vertex `v` then has `delta-1>=2` in-neighbours. Cut a Hamiltonian cycle of
the deleted core after each such neighbour and append `v`. These paths have
different named starts, so choose the second start different from the first.
After `j` reversals always ending at `v`, if `N_j(u)` counts actual starts,

`d_j(v)=n-1-delta+j`, `d_j(u)=d_0(u)-N_j(u)` for `u != v`.

At the second move its new start has degree at least `delta`, so its endpoint
sum is at least `n`. For `j>=2`, the split first two starts ensures
`N_j(u)<=j-1`, and every later ending-at-`v` path likewise has endpoint sum
at least `n`. All are black. Strong connectivity supplies these paths until
move `delta`.

At the end `v` is a source. Every other vertex has positive outdegree.
If the deleted core were nonstrong, a terminal strong component would form
a zero cut of the whole tournament. A singleton contradicts the positive
outdegree; a larger proper shore contradicts the preserved nonsingleton cut
margin. The complementary core is therefore strong.

Converse gives the sink construction when `d+(v)=delta`. To translate to
native `K`, alternate the tournament's converse and reverse the corresponding
path order whenever needed. After `delta` moves the native endpoint has the
same role for even `delta`, the opposite role for odd `delta`. The lower bound
is exact: in the `R` gauge a named score changes by at most one per move, so
creating an extreme vertex costs at least its smaller initial semidegree.

A regular tournament of order `2m+1>=7` satisfies the region with `delta=m`,
and every vertex is extremal on both sides. Hence **every prescribed vertex
can be either native source or native sink, with strong complementary core,
in exactly `m` black moves**. For nonregular inputs the parity/direction table
must not be replaced by a claim about both roles of every extremal vertex.
The strong pure-blue `R5` is a hostile to dropping the size/connectivity guard.

## A separate explicit family with exact distance two

Let a bottom tournament `S` have marked root `x` and
`S-x=B1 => ... => Bq`, where each block is `K1` or `C3`, `q>=2`, and the last
block is `C3`. The root is uniform on each block, beats the first and loses
to the last; intermediate signs are arbitrary. Write `|S|=h>=5`.
Form `E_t(S,x)`, `t>=1`, from bottom `S`, middle chain
`C=(c1,...,c_(t-1))`, and a top converse copy `sigma(S)`. All interlayer arcs
point down except the upward chain `x -> C -> sigma(x)`.

Name the last bottom cycle `r -> r1 -> r2 -> r`. Let `A=(a1,...,ak)` be a
Hamiltonian path of the earlier bottom blocks, `k=h-4>=1`, and let
`U=(u1,...,uh)` be a top path starting at `sigma(x)`. Concatenating block
paths supplies both choices. The words

`P=(r,r1,r2,x,C,U,A)`

and, in `R_P(T)`,

`Q=(r,x,a1,reverse(U),reverse(C),reverse(a2,...,ak),r2,r1)`

are Hamiltonian. Empty lists are omitted. `P` closes through `ak -> r`;
`Q` closes through `r1 -> r`, so both the initial and first intermediate
are strong. The first endpoint degrees sum to at most `2+5=7<n-1`, and
the second sum is `1+2=3<n-1`. Both words are black.

After both path reversals, the only two original out-arcs of `r` have
been reversed, making `r` a sink. The remaining vertices other than `r1`
lie on the directed cycle

`(x,r2,a2,...,ak,C,U,a1,x)`.

The missing vertex has attachments `a1 -> r1 -> x`, making the entire
complement strong. Each displayed bridge is either an unchanged downward
interlayer edge, a stack edge, or a consecutive edge reversed once/twice;
the supplied checker verifies the literal labelled construction.
The native words are `P` then `reverse(Q)`, with exactly the same final
labelled tournament. Initial minimum semidegree is two (attained at the
terminal C3 and its converse copy), so no single native move can create
an extreme vertex. The strong-core distance is exactly two.

This family theorem does not depend on a universal exceptional classification.
Its lower bound concerns singleton targets, not every possible nonstrong cut.
