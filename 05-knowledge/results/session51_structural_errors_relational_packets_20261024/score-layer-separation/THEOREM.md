# Equal-multiplicity score layers have disjoint metagraph neighbor bands

24 October 2026. **PROVED**, as an all-size consequence of the inherited actual
Hamiltonian-path transport identity. The finite controls are not the proof.
No primality, cyclic-cell, refinement, or self-converse hypothesis is needed
for the unrestricted path statement.

## 1. Statement and notation

Let n=2p>=2 and let T have exactly two outdegrees a<b, each occurring p times.
The degree sum forces a+b=2p-1. Put k=b-a, a positive odd integer. Write

    c(a) = binom(2p,3) - p[binom(a,2)+binom(b,2)]

for its number of directed cyclic triples. A **score layer** is the set of
all ordinary or converse-merged classes with this degree multiset. Converse
interchanges a and b and preserves their multiplicities, so the layer and
c(a) are well-defined on the merged node set.

For any actual directed Hamiltonian path P, let K_P(T) retain the path arcs
and reverse all other arcs, the native line operation. Its cyclic-triple
count belongs to the necessary-value set

    {c(a)-1, c(a)-k-1, c(a)+k-1}.

For a blue path, only the last two possibilities remain. These are **allowed
necessary values**, not a claim that all values are nonnegative or attained.

**Separation theorem.** Distinct equal-multiplicity two-score layers of the
same even order cannot share a neighbor under the native line operation.
Their full intervals of possible neighbor c3 values are disjoint. No edge
between two vertices in the union of these layers can occur except within
the balanced layer a=p-1, b=p. Within that layer a necessary condition is
that the retained path start at degree b and end at degree a.

This statement covers both blue and black line operations. It gives a wiring
obstruction, not a sufficiency criterion, existence theorem for paths, or
connectivity theorem. The degenerate n=2 operation is the identity because
there are no non-path arcs; the port statement remains true, without imposing
a nondegenerate edge convention at that order.

## 2. Proof from inherited endpoint transport

Session40 `marked-ports/PATH_TRANSPORT_AND_FOREST_RULE.md`, Section1, proves
that for any actual Hamiltonian path with start s and end t,

    d_(K_P(T))(v)=2p-1-d_T(v)+1_(v=s)-1_(v=t).

This source is archived and pinned. Since s and t are distinct when n>=2,
squaring the degrees and using the fixed degree sum gives

    sum_v d_(K_P(T))(v)^2 - sum_v d_T(v)^2
      = 2(d_T(t)-d_T(s)+1).

The cyclic-triple identity

    c3(T)=binom(2p,3)-sum_v binom(d_T(v),2)

therefore yields

    c3(K_P(T))-c3(T)=d_T(s)-d_T(t)-1.

The endpoint score difference is0,+k,or-k, proving the three-value statement.
For a blue path, positional reversal is an anti-automorphism and exchanges
s with t, so d_T(s)+d_T(t)=2p-1=a+b. Since a!=b, the endpoint degrees must
be opposite, giving only the +/-k cases. This does not assert that every path
with opposite endpoint degrees is blue.

The p low-score vertices induce binom(p,2) arcs, so pa>=binom(p,2); hence

    ceil((p-1)/2) <= a <= p-1.

The upper bound follows from a<b and a+b=2p-1. For p=1 or2 this leaves only
the balanced layer a=p-1. For p>2 compare consecutive admissible integer
layer values a and a+1. Their gaps are k and k-2, with k>=3. Direct subtraction
gives

    c(a+1)-c(a)=p(k-1).

The necessary neighbor interval for a is

    I_a=[c(a)-k-1, c(a)+k-1].

The distance from the upper endpoint of I_a to the lower endpoint of I_(a+1)
is

    [c(a+1)-(k-2)-1] - [c(a)+k-1]
      = (p-2)(k-1) > 0.

Thus all layer intervals are strictly separated. Skipping an unrealized or
inadmissible layer cannot weaken the separation. A target class has one
cyclic-triple count, so cannot neighbor two different layers.

Moreover c(a) itself lies in I_a because k>=1. Therefore a vertex of a
different layer cannot be a neighbor: its cyclic-triple count lies in its
own disjoint interval. A neighbor within the same layer would require

    d_T(s)-d_T(t)-1=0.

Among the allowed endpoint differences0,+k,-k, only +k=1 works. Hence k=1,
with the start at b and the end at a. This proves the theorem. ∎

The argument requires no construction of all tournaments in a layer and no
class-count formula. It constrains their actual line neighbors by a retained
path operation, not by an arbitrary feasible graph with matching margins.

## 3. Order-ten application and exact stopping boundary

The positive-error self-converse nodes at order ten have score layers
3^5,6^5 and4^5,5^5. Their c3 values are respectively30 and40. A blue neighbor
therefore has c3 in the disjoint sets

    degree-one crossing layer: {26,32};
    degree-two crossing layer: {38,40}.

For arbitrary paths the corresponding necessary sets include29 and39,
respectively. It would be false to state only two values for black paths.

Consequently the two source families have disjoint complete neighbor sets,
and every positive-to-positive blue edge lies in the balanced4^5,5^5 layer.
This is the all-size structural mechanism behind the independently enumerated
order-ten boundary separation. The separately frozen neighborhood packet
identifies238 zero-error neighbors of the four unbalanced defect nodes and793
zero-error neighbors of the twenty balanced defect nodes. Those cardinalities
and connectivity within their incident-edge graph are finite results, not
consequences of the interval argument alone.

No assertion is made that these two neighborhoods lie in different connected
components of the full metagraph. Other edges among their boundary nodes or
through additional score types can connect them.

## 4. Full-metagraph distance and loop consequences

Because different score layers have neither an edge nor a common neighbor,
any finite distance between a node in one such layer and a node in another
is at least three in the full merged metagraph. This is a distance lower
bound, not infinite separation: longer paths through other score types remain
possible. It applies to the full graph with both line colors, not merely the
positive-error induced slice.

More generally, at any order n>=2, a loop at a converse-merged node requires
an actual Hamiltonian path whose retained-path partner is isomorphic either
to T or to T^op. The cyclic-triple count is unchanged by either isomorphism
or converse. Hence the endpoint identity requires

    d_T(start) = d_T(end)+1.

Therefore a tournament whose outdegree multiset contains no consecutive
values has no native C-loop at its merged node. In particular every unbalanced
equal-multiplicity two-score layer (k>=3) is loopfree in the full merged
metagraph. These are necessary conditions only; adjacent degree values do
not guarantee a Hamiltonian path realizing a loop.

The endpoint leg-law obstruction itself is inherited: THM-4731,
`THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
statement4 proves it for both ordinary self-lines and reversal-neutral
self-lines, citing THM-793 for the former. Its archived copy here is
`sources/THM-4731-self-line-leg-law.md`. The present proof uses the same
inherited endpoint transport to derive the layer-separation and distance
consequences; no originality claim is made for the general loop leg law.

## 5. Controls and the unattained-value hostile

`check.py` checks the interval arithmetic for every possible layer at p<=500,
and exhausts every labelled source of the stated two-score type at n=2,4,6.
For every actual directed Hamiltonian path of each source, it constructs the
partner by literal non-path reversal and counts its cyclic triples directly.
It checks the endpoint jump, the blue restricted list, and the restriction on
partners that themselves lie in a two-score equal-multiplicity layer.

At n=6 the two ordered directed triangles have scores1^3,4^3 and c3=2. One
formal entry of the allowed neighbor list is-2, which is impossible and is not
attained. The checker explicitly preserves this hostile so that the theorem
cannot silently promote a necessary-value list into an attainment claim.

Normal and optimized Python receipts are byte-identical and reproduce after
a real-copy relocation. The all-size proof is algebraic; the finite interval
range and small actual-port universes are independent controls only.
