# A literal tournament register lift of PRIMEGAME, and its decode wall

Set the physical block for prime 2 equal to a singleton (H=1),
and the block for prime 7 to a strong H=49 tournament, rather than
attempt the forbidden H=7 block. For each of the remaining
participating primes p take a literal strong tournament with H=p:

| register | outgoing masks | H | Q(2) |
|---|---|---:|---:|
| 2 | (0,) | 1 | 2 |
| 3 | (2,4,1) | 3 | 12 |
| 5 | (2,4,9,3) | 5 | 32 |
| 7 | (2,4,25,35,75,23,47) | 49 | 824 |
| 11 | (2,20,9,3,13) | 11 | 100 |
| 13 | (2,12,9,17,7) | 13 | 104 |
| 17 | (2,4,9,19,39,15) | 17 | 224 |
| 19 | (2,4,41,19,7,27) | 19 | 256 |
| 23 | (2,4,25,19,35,15) | 23 | 272 |
| 29 | (2,12,17,21,35,15) | 29 | 320 |

Use v_p(m) copies of each primitive block, in a fixed canonical
prime order, and ordinal-join the blocks. Registers are recoverable
from the graph: these blocks are distinct strong SCC types, and
the singleton lane is a count of singleton SCCs. Every reduced
FRACTRAN denominator is permitted exactly when enough of the
corresponding blocks are present. Delete those blocks and insert
the numerator blocks in their prescribed group positions. Each
such step is a literal tournament construction preserving the
arithmetic priority and register update. Ordinal stacking multiplies
slot polynomials because each occupied slot traverses the SCCs in
their fixed order. Thus this really is an exact physical lift of
the register states, not a request to realise non-tournament
intermediate counts.

What it **is not** is a proof about Hamilton realization for the
output exponent. When an observed integer is a pure power 2^p, its
graph is the transitive TT_p, with H=1 and Q(z)=z^p. The hidden count
p is available to a marked additive observer:

    Q'(1)=p,        sum_v H(TT_p−v)=p.

The second form counts ordered two-slot covers in which the first
slot must contain exactly one vertex. The short uniform observers
cannot impose this length requirement.

Indeed consider **any** uniform exterior context, even chosen
separately for each p. By Matching Sections (Session 331), if the
context's exterior crossing graph has no matching edge, only
the one-visit term survives. Its count is exactly the quotient's
old Hamilton count, so making it p presupposes a p-realizing
quotient. If the crossing graph has a matching edge, the
coefficient of P2(TT_p) is positive, and

    H(output) ≥ P2(TT_p)=2^p−2 > p       (p≥3).

This exhausts the uniform-module decoder options. Moreover no
fixed finite linear combination of such readouts, even allowing
signed rational coefficients, can agree with p for unbounded p.
For fixed auxiliary budget each readout is a finite linear
combination of j^p: P_k(TT_p) counts the surjections from p
vertices to k occupied slots, and inclusion-exclusion gives
Σ_j (−1)^(k−j) binom(k,j) j^p. A fixed finite combination of
these either grows exponentially in magnitude or is eventually
constant. The same obstruction holds if one switches among
finitely many such signed readouts.

This identifies a precise remaining obligation for any useful
PRIMEGAME→Chan–Pak construction: preserve a physically usable
nonuniform mark that can extract a length-one side of a split
(or find a different genuinely additive encoding). Prime
generation and register regeneration alone do not provide this.

## Even nonuniform auxiliary contacts cannot activate the unary register globally

Suppose TT_p is extended by b exterior vertices with arbitrary
directions. The intersection of a strong component with the unary
chain is an interval in the chain order: if i<j belong to it then
i→k→j for each intervening k, and a path j↝i proves k is in the
same component. Every nonsingleton strong component contains at
least one exterior vertex, since the chain alone is acyclic. Thus
the active unary vertices lie in at most b intervals; the rest
are singleton strong components and cannot affect H.

For a prescribed prime count p≥3, there is only one nonsingleton
strong component, because Hamilton counts multiply over components
and a nonsingleton strong component has count at least three.
The same conclusion holds for target 7p with p prime if the target
is realised: every proper factorisation has a forbidden 7 (the
p=7 case gives 7×7). The inherited Busch/core bound then restricts
that single active interval, including its auxiliary vertices,
to O(log p) total order. In particular a decoder whose entire
p-vertex chain must become strong cannot work at large p. A
nonuniform decoder is not ruled out, but it must select a small
active interval or change the physical encoding; simply extracting
the count of marked singletons throughout the whole unary chain
cannot be a strong-output proof.

The same uniform-decoder wall applies to 7p for p≥7: a split
forces at least 2^p−2>7p paths, and a nonsplitting context already
needs a quotient with the output count. This does not identify the
spectrum with PRIMEGAME or rule out all other certificate codings.
