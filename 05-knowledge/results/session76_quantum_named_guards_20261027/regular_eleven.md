# Regular eleven: fifteen native moves, a fourteen-black-move hub, and endpoint clocks

## Statement and evidence boundary

Every regular tournament of order eleven has native distance **15 to the
transitive isomorphism class**. More strongly, every such class admits
fourteen black moves to the same strong class `U11`, followed by one blue
move to the transitive class. The fourteen-move part remains strong.
Consequently all regular eleven-vertex classes lie in one component of the
strong black metagraph, with pairwise black distance at most 28.

The first assertion and its discovery belong to
[THM-4815(B8′)](https://github.com/eliottcassidy2000/math/blob/3858d2afca72fb97a46a3601a6ae4b1e6e99a812/01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md).
The present finite verification regenerated a complete named route for each
of the 1,223 ordinary isomorphism classes in
[McKay's regular-tournament catalogue](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html),
then checked paths, native moves, colours, inverse moves, relabellings and
strong connectivity independently. Both ordinary and optimized Python
checks agreed. Catalogue completeness is **cited**, not a new exhaustive
enumeration of all labelled eleven-vertex tournaments. Converse classes
were not identified.

The compact script `check_native_route.py` reproduces one complete
fourteen-black-plus-one-blue certificate and its clock identities. It does
**not** reproduce the entire 1,223-class coverage calculation. That distinction
is important: a sample certificate proves its sample, not a census.

## The move and the lower bound

Let `Q=(q0,...,q(n−1))` be a directed Hamiltonian path. A native move keeps
its consecutive arcs and reverses every other arc. Thus the same named path
remains legal, and applying the move twice is the identity. If `s` is the
outdegree vector, with path endpoints `a,b`, then

    s′ = (n−1)1 − s + e_a − e_b.

For sorted scores `s_(0)≤...≤s_(n−1)`, define

    L(T) = ½ Σ_i |s_(i)−i|.

Converse preserves `L`; changing two endpoint coordinates by opposite units
changes `L` by at most one. This follows either from sorted matching on the
line or from the minimum-over-permutations formula for the L1 distance to
`(0,1,...,n−1)`. A transitive tournament has `L=0`. For regular order
`n=2m+1`,

    L = ½ Σ_(i=0)^(2m) |m−i| = m(m+1)/2 = (n²−1)/8.

At order eleven this is `1+2+3+4+5=15`. Legal fifteen-move certificates give
the matching upper bound. A general theorem producing *some* regular
endpoint from a transitive source does not establish the stronger
quantifier “every regular class”; the finite coverage supplies that part.

## Colour and the common strong hub

A path is blue when its positional reversal is an anti-automorphism of the
current tournament, and black otherwise. This condition is invariant under
the native involution. A transitive tournament has just one Hamiltonian
path, whose reversal is an anti-automorphism. Therefore no black edge is
incident to the transitive class. Its black-only distance from a
nontransitive tournament is infinite, not fifteen.

Starting with the standard transitive order `0,...,10`, make its unique
native move and call the result `U11`. Its sorted scores are

    1,1,2,3,4,5,6,7,8,9,9,

so `L(U11)=1`. It contains the directed cycle
`0→1→...→10→0`, hence is strong. Each regular input has `L=15`, so any route
to `U11`, including any black route, needs at least fourteen moves. The
separately verified one-blue certificate set supplies exactly fourteen
black moves, proving equality. The final blue move is unavoidable, so the
minimum blue budget along an optimal route is one.

This is a class-distance statement. A prescribed labelled transitive target
can impose extra parity restrictions. It is also different from the
five-black-move distance to a singleton plus a strong core.

## A general conditional strong-corridor theorem

The following proof applies to any odd order, **if** a regular/transitive
route attains the score lower bound. It is not an existence theorem for
all regular tournaments at every order.

Put `n=2m+1`, `D=m(m+1)/2`, and

    R(T) = ½ Σ_v |s_v−m|.

The coordinatewise triangle inequality gives `L+R≥D`. Equality holds exactly
when every sorted score `s_(i)` lies between `i` and `m`.

**Lemma.** A nonstrong tournament with neither a source nor a sink has
`L+R>D`.

Indeed, its initial and terminal strong components are disjoint, so one
has size `k≤m`. If this is the terminal component, absence of a sink makes
`k≥3`. Every score there is at most `k−2`; its vertices occupy the first
`k` sorted positions, since all outside vertices beat them. At index
`i=k−1<m`, the inequality `s_(i)≤k−2<i` contradicts equality in `L+R≥D`.
For an initial component use the dual argument: its smallest score is at
least `n−k+1`, at index `n−k>m`, again violating equality.

In a tight `D`-move route from transitive to regular, `R` drops by one and
`L` rises by one at every move. Hence `L+R=D` throughout. An `R`-decreasing
move starts above and ends below the mean; endpoint absolute deviations
fall by one, while other absolute deviations stay unchanged. The first
move removes the transitive source and sink. No later absolute deviation
can restore an extreme score. The lemma therefore forces every state after
the first move to be strong. Reverse the route: it remains strong until
its final transitive target.

## Named endpoint clocks

For any tight regular-to-transitive word, let `s_t` be the score vector after
`t` moves and set

    z_t = (−1)^t (s_t−m1).

Index moves by `t=1,...,D`, with endpoints `a_t,b_t`. Then

    z_t−z_(t−1) = (−1)^t(e_(a_t)−e_(b_t)).

These are roots of `A_(n−1)`, each of L1 norm two. The final centered score vector is
a permutation of `−m,...,m`, so its L1 norm is `2D`, exactly the sum of step
norms. Equality in the triangle inequality forces **no coordinate
cancellation**. A vertex of final score `j` appears as an endpoint exactly
`|j−m|` times; its start/end role is fixed by final rank and integer-clock
parity. In particular, the median-score vertex is interior in every path.
For order eleven the endpoint multiplicities are

    0,1,1,2,2,3,3,4,4,5,5.

These are necessary constraints on a tight word, not a replacement for
checking Hamiltonian-path legality.

## What remains open here

The all-order black-hub formula is not established. It already fails at
order five: the regular five-tournament is a pure-blue atom and cannot
reach the transitive neighbour by black moves. The order-eleven theorem
does not classify the whole strong black metagraph or assert that every
pair of regular classes is distance 28.
