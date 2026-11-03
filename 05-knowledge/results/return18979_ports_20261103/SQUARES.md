# Minimum connectors require compatible endpoint pairs

**PROVED general deletion/splicing criterion; FINITE-EXACT complete arithmetic
classification.** The original anchors, recovered from
[S173 section4](../session173_modular_clocks_20261103/squares/NOTE.md), are

    {527,2,7,44,5,4,13,3}.

Their three components are527--2--7,44--5--4,13--3. Thus the supplied path

    527--2--7--317--44--5--4--12--13--3

adds **317 and12**. The extension3--1--8 adds1 when8 is another required
anchor. S173 already proves minimal connector cardinality and height for
these paths. We extend that result to a complete catalogue, and credit
its deletion/factor-pair mechanism as inherited.

## SC1. An exact minimum-splice carrier

For any finite simple graph G, deletion of k added vertices from a Hamilton
path leaves at most k+1 nonempty old paths. Hence k>=pc(G)-1, where pc
is the minimum number of vertex-disjoint paths covering G. In particular
k>=c-1 for c components.

When each component is traceable, equality k=c-1 holds iff there exist:

1. an ordering of the components and an oriented Hamilton endpoint pair
   for each, with an actual internal path;
2. one permitted new vertex joining each successive endpoint pair;
3. distinct new vertex names, all outside the old set.

Necessity follows because deletion must give exactly c nonempty pieces,
one per component. Thus no added vertex can be a path endpoint, and no two
added vertices can be adjacent in the path. Each old piece spans its whole
component. Sufficiency is literal concatenation. Extra graph edges do not
invalidate this equivalence.

Endpoint marginals alone lose which endpoints belong to the same spanning
path. The connector lists also require a shared all-different constraint:
two individually legal bridges cannot reuse the same integer. These are
joint-port data, not a scalar connectedness test.

For distinct positive endpoints a>b, every positive common square neighbour
z is obtained from a factor pair uv=a-b, u<=v, u=v mod2, by

    z=((v-u)/2)^2-b.

Filter z>0 and exclude every occupied old vertex. This gives all candidates
without an arbitrary search cutoff. Exhaust component orders, internal
Hamilton paths and these finite lists, then impose distinct connectors.
The construction is therefore a complete search over all positive integers.

## SC2. The complete catalogues

|Required anchors|Components|Minimum additions|Connector sets|Paths modulo reversal|
|---|---:|---:|---:|---:|
|Original eight, including44|3|2|9|9|
|Original eight plus8|4|3|17|17|

Every listed connector set has exactly one path modulo reversal. Full
lists and factor certificates are in `square_splices.json`. The unique
sets minimizing the largest added vertex are {12,317} and {1,12,317}.

There is a small hostile control: adding {212,317} makes the original
graph connected, but yields no Hamilton path. Its edges are the original component edges,
7--317,44--317,317--212,13--212. Its degree-three vertex is317, and its
leaves are527,4,3. A spanning path in this tree is impossible, since it
would have to cover three leaves. Independent subset DP also gives zero.

The attachment-response filtration has its first nonzero Hamilton count
at depth2 for the original set and at depth3 with8: fewer births cannot
join that many old components, and the supplied paths attain the bounds.
This is an existence depth, not a linear response-rank or query-count claim.

## Boundary correction retained explicitly

The initial exploratory code instead treated317 as an original anchor and
44 as a connector. That is a different legitimate three-component problem:
it has13 minimum connector pairs, with unique minimum-height set {12,44}.
The initial commentary's33-triple count likewise referred to that alternate
anchor set plus8. It did not describe the user's inherited anchor problem.
The current receipt keeps the13-pair alternate case separately, and the
maintained source boundary is the original44-anchor set. The correction is
also logged in MISTAKES. The failed implication was reading connector roles
from the final path alone; a path does not determine its original boundary.

The tournament analogy is the deletion/splicing operation and its retained
ordered endpoints. No orientation is assigned to nonedges or ties in this
undirected square-sum graph, and no tournament-spectrum theorem is inferred.
