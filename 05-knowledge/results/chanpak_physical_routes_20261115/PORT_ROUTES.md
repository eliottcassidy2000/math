# Two literal switches in a fixed residue fiber

This note constructs an infinite family of physical paths of length two.
It does not claim that previously found arbitrary below- and above-target
witnesses lie in the same component.

## Port substitution

Encode a tournament on `n` vertices by bits in lexicographic order of pairs
`(i,j), i<j`, least significant first, with a 1 meaning `i→j`. Fix
vertex `v=n−1`; replace it with any nonempty tournament `D`, with every
arc to a vertex outside the module oriented as the original arc to `v`.

Write `P_k(D)` for the count of ordered covers of D by k nonempty directed
paths. For each permutation of the n−1 external vertices, mark an internal
gap mandatory if there is no forward arc across it. Reject the permutation
if any mandatory gap cannot take a visit via the port. A forward internal
gap is optional if both arcs of a visit via the port exist. The two end
gaps are optional precisely when their respective arc to or from the port
exists. If the resulting mandatory and optional counts are b and a, the
response is `z^b(1+z)^a`. Summing gives the response polynomial `F(z)`.
Its coefficient at z^k multiplies `P_k(D)` exactly: assign the k ordered
path runs to the chosen gaps. The constant coefficient is irrelevant
because D is nonempty.

## A route with the previous dyadic clock modulus

The following three templates are all strongly connected; the two edits
are single arc flips external to the port:

| Template (n=7, v=6) | Response F | Next changed bit |
|---|---|---|
| 366909 | 15+27z+9z² | 65536, arc (3,5) |
| 301373 | 15+29z+10z² | 32768, arc (3,4) |
| 268605 | 17+31z+11z² | — |

Let `H=P_1(D)` and `Q=2P_1(D)+P_2(D)`. Then actual physical counts on
this arc-flip path are

    9H+9Q  →  9H+10Q  →  9H+11Q.

All inflated graphs are strong: every path in the quotient lifts to a
path between the corresponding modules, and two vertices within the
same module can connect by following a quotient cycle through it.
The intermediate hits the midpoint, without ever leaving its residue
class modulo `M=Q`. This is a route rather than a pair of separately
constructed congruent witnesses: the only changes are the two listed
arcs, and the entire port D is held fixed.

For all `k≥2`, take D to be the ordinal composition of a directed
triangle and k−2 singletons. The count H is 3 and Q is `12·2^(k−2)`,
since optional-slot counts multiply under ordinal composition: each
ordered slot restricts to an independent ordered slot in each
component, with concatenation in ordinal order. Thus for precisely
the earlier growing powers-of-two modulus `M=3·2^k` we get the
strong fixed-residue route

    9M+27  →  10M+27  →  11M+27.

The attained midpoint here belongs to a prescribed sparse lane: the
construction by itself does not place this route at a general target
of size far larger than M.

## A distinct readout: fixed path-count step, adjustable slot-count carrier

A second all-strong route has n=6, v=5:
`29130→29002→29000`, with external xor masks 128 and 2. Its polynomials
are `9+27z+21z²+4z³`, `11+29z+21z²+4z³`, and
`13+31z+21z²+4z³`. It therefore has counts

    4Q₃+9Q₂−3H  →  4Q₃+9Q₂−H  →  4Q₃+9Q₂+H,

where `Q_j(D)=Σ_k binom(j,k)P_k(D)` counts ordered covers allowing
empty slots. This route has step `2H`, and permits varying Q₂,Q₃
while fixing H, in contrast to the Q-step construction. For a fixed
port D with count h, adding a dominated singleton to D keeps h and
multiplies Q₂ by 2 and Q₃ by 3. Thus the carrier can go to infinity
without changing its residue-step size, although a fixed seed still
has exponential holes between tower levels. In particular this
does not supply even low-scale cofinite coverage.

The companion script verifies the signatures exactly and counts
several inflated examples by direct Hamilton-path enumeration by
subset DP. The identity for arbitrary D comes from the gap proof,
not from those test cases.

## What certifies transport under a second port

The same run enumeration works with two marked modules: record a
bidegree (r,s) for the counts of maximal visits to the two modules;
the coefficient multiplies P_r(D)P_s(E). Each degree is at most
the number of external singletons plus one, since two visits
to a module cannot be separated solely by a single visit to the
other module: the intermodule arcs are one-way. This is finite,
even when both substituted tournaments have arbitrary order.

To retain the exact step Q(D) under every transitive enlargement
E, the difference of adjacent response tables must be precisely
`(2z+z²)w`, where z and w count the respective visits.
Necessity is not merely a conjecture behind a filter: on a transitive
tournament of order t, P_s vanishes when s>t, and P_t=t!. Thus
testing enough transitive sizes in (D,E) gives a triangular
tensor-product system with nonzero diagonal, and determines all
coefficients. Conversely that table gives Q(D)H(E), which is
Q(D) for every transitive E. A useful new strip must also have a
baseline response with some coefficient at w^s for s≥2; otherwise
a larger transitive E just pads the graphs without moving the counts.

The seven-vertex Q-strip above does not have such a second port:
the only vertices external to the changed arcs and outside port D
are 0,1,2. Their first-edge responses contain the following
coefficients of z w², z² w², z w³, z² w³ respectively:
(0,1,0,0), (3,1,0,0), (−5,−2,−2,−2). Each violates the exact
separation test. This rejects the displayed strip as a uniform
two-port construction; it does not reject other strips.

There is a stronger exact failure for these same ports. The second
step minus the first step, after insertion of arbitrary E, evaluates
as follows (h=H(D), p=P₂(D)):

| Second port | Second step minus first |
|---:|---|
| 0 | 2h·P₂(E) |
| 1 | (2h+p)·P₂(E) |
| 2 | (4h+2p)·P₂(E) |

These are the complete second differences of the run coefficients,
not an estimate. Thus no nontrivial inflation at an available second
port can even keep the two steps *equal*, let alone keep them at the
original modulus. Indeed P₂(E)>0 whenever E has at least two
vertices, by splitting a Hamiltonian path. A new reader or a
third coordinated operation must compensate the actual defect.
