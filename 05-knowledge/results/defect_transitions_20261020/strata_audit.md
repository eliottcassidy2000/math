# Independent deficiency audit and the third polynomial coefficient

**PROVED below:** the exact deficiency-two orbit count and the resulting
third eventual polynomial coefficient, using a complete local collision
classification. **FINITE-EXACT:** the earlier independent normalization
census at budgets one through four agrees with the formulas. No extrapolation
from those four values and no external novelty claim is used.

The audit of [strata.md](strata.md), propositions 1--3, found no correction.
The identity
\(D=2(k-r)+(2r-t)+u\), the endpoint-weight bound, the \(3D\) exceptional
endpoint bound, and the exhaustive deficiency-two menu are valid.
The code retains the displayed weighted cost and recomputes deficiency
after canonical normalization, as required. This note supplies the
previously missing quotient count for that menu.

## 1. The orbit count to be proved

Fix \(k\ge1\), cap \(K=k+1\), and heavy count \(h=2k-1\). Put

\[
 M=(2k-1)!!,\qquad P=(2k-3)!!=\frac{M}{2k-1},
\]

where \((-1)!!=1\). The number \(N_2(k)\) of capped weighted-core orbits
in this stratum is

\[
\boxed{\displaystyle
 N_2(k)=
 M(2k^4+7k^3+2k^2-3k+1)
 +\frac{P}{3}(3k^3-16k^2+13).}
\tag{1}
\]

The values are \(9,264,5365,102750\) for \(k=1,2,3,4\), agreeing with the
existing complete deficiency-two compiler.

Use the canonical transitive order of the heavy vertices. Every light
vertex has an intrinsic gap relative to that order, since no backward
arc can meet a heavy vertex. Isomorphisms preserve the heavy order and
each light gap. This reduces quotient ambiguities to the small clusters
inside exceptional gaps; the ordinary matching background stays marked
by its ordered endpoint positions.

Call an exceptional gap **occupied** if it has a light vertex of weight
\(1,\ldots,k\), and **empty** if its weight is zero.
Every good canonical presentation has at least \(h-1=2k-2\) defect
endpoints, so its core order is at least \(4k-3\). The menu bounds it
above by \(4k+1\). The five core-order bins counted below therefore
exhaust the stratum.

## 2. Two exceptional gaps in a matching

The first menu row consists of a matching on \(2k\) endpoints and two
exceptional gaps among its \(2k+1\) gaps.

### Both gaps occupied: \(k^3(2k+1)M\)

There are

\[
 \binom{2k+1}{2}k^2M=k^3(2k+1)M
 \tag{2}
\]

presentations. Their underlying unweighted cores are unchanged from the
top matching cores, so they remain pair free, of order \(4k+1\).
They are injective even when the two demoted gaps are adjacent.

Here is the local injectivity check. Separated demotions create disjoint
clusters of at most three vertices; a transitive cluster has unique order,
and a directed triangle whose endpoints are internally matched has no
external defect port, so its rotations do not change the template.

Two adjacent internal demotions create a five-vertex cluster, with three
endpoint vertices and two unmarked vertices between them. There is at
most one internal matching edge.

- With no internal matching edge the cluster is transitive.
- An edge between consecutive endpoint positions gives a directed triangle
  followed or preceded by two transitive singleton components. Their
  strong-component order distinguishes those pieces.
- An edge between the first and third endpoint gives
  \(C_3[1,3,1]\). Its unique three-vertex transitive block has a fixed
  internal order; its middle endpoint is additionally distinguished by
  its external matching partner. Both unmarked weights are recovered.

The analogous exterior cluster has at most four vertices: it is transitive
or a singleton ordered beside a directed triangle. These cases recover
the demoted positions and weights. A triangle rotation can change a
chosen vertex within one presentation, but cannot identify different
parameter tuples.

### One occupied gap and one empty gap: \(4k^3M\) new states

The raw count is \(2k^2(2k+1)M\). A merger occurs precisely when the empty
gap lies between two adjacent matched endpoints.
Summed over all matchings there are \(M\) such adjacent matched pairs.
For each pair there are \(2k\) choices of the other, occupied gap and
\(k\) choices of its weight. Hence \(2k^2M\) presentations merge.

If the occupied gap is adjacent to the empty one, the matched pair and
its neighboring unmarked pieces collapse to one heavy chain. The result
is a matching on \(2k-2\) endpoints with all gaps heavy; it is counted in
section 5 below. Otherwise at least two old heavy vertices merge and the
heavy count drops below \(2k-1\), so the output leaves this stratum.

Every remaining template is pair free and injective, giving

\[
 2k^2(2k+1)M-2k^2M=4k^3M
 \tag{3}
\]

new states, of core order \(4k\).

For completeness, the only additional local case compared with (2) is
a four-vertex cluster from adjacent exceptional gaps. If its three
endpoints have an internal matching edge, the cluster is either a
singleton beside a triangle or \(C_3[1,2,1]\). In the latter case the
external matching port identifies which of the two middle-chain vertices
is an endpoint. Its position distinguishes which gap was empty. If the
edge joins the two endpoints bordering the empty gap, it is exactly the
excluded merger. With no internal edge the cluster is transitive.

## 3. Two empty matching gaps and one occupied \(P_3\) gap

These two presentation families have the same core order \(4k-1\), and
have one genuine cross-family collision.

### Pair-free two-empty-gap presentations

For a matching let \(a\) be its number of adjacent matched pairs.
Both chosen empty gaps must avoid these \(a\) positions.
Thus the number of pair-free presentations is

\[
\begin{aligned}
 E&=\sum_{\text{matchings}}\binom{2k+1-a}{2}\\
  &=k(2k-1)M+(k-1)P.
\end{aligned}
\tag{4}
\]

Indeed \(\sum a=M\). To count \(\sum\binom a2\), choose two disjoint
adjacent pairs of endpoint positions and match the remaining vertices:

\[
 \sum\binom a2
 =\frac{(2k-2)(2k-3)}2(2k-5)!!=(k-1)P.
\]

The formula holds separately at \(k=1\), where the second term is zero.

Every excluded presentation either loses a heavy vertex or collapses to
the all-heavy matching family of section 5. Every retained presentation
is injective. Its only nontransitive light cluster is a triangle produced
by two consecutive empty gaps whose outer endpoints are matched.
The middle endpoint has one external matching partner; that external
port fixes the middle endpoint, and the cyclic order then recovers the
presentation.

### Occupied \(P_3\)-gap presentations and the unique collision

There are \((k-1)M\) defect graphs consisting of one \(P_3\) and the
remaining matching edges. Choosing an occupied exceptional gap and its
weight gives

\[
 O=2k^2(k-1)M
 \tag{5}
\]

raw presentations. All remain pair free because every gap is occupied.

The only ambiguity is a three-vertex light gap inducing \(C_3\), with
all three weights equal to one, and exactly one external defect port.
It has exactly three presentations in the union under discussion:

1. one matching presentation with two empty gaps, placing the externally
   connected endpoint in the middle;
2. two \(P_3\) presentations with one occupied gap, placing that endpoint
   at the first or last position in a cyclic ordering of the triangle.

The external endpoint and its external partner are intrinsically marked
by the heavy gaps. The three cyclic positions therefore exhaust the
ambiguity. If the unmarked weight exceeds one, that vertex is identified
by its weight and the two \(P_3\) presentations cannot coincide.

No other overlap is possible. Distinct heavy-gap occupancy patterns
already separate most cases. A transitive light cluster has unique
order; its external defect incidences then distinguish endpoint vertices
from unmarked vertices. A triangle without an external port cannot arise
from the three matching endpoints of the two-empty-gap description.
Clusters of size two are rigid.

The number of these triply presented states is

\[
 J=(2k-2)P:
\tag{6}
\]

choose three consecutive endpoint positions, match the outer two, and
match the remaining \(2k-2\) vertices. The middle vertex's partner is
necessarily external.

Consequently the union count for core order \(4k-1\) is

\[
 E+O-2J=(2k^3-k)M-3(k-1)P.
 \tag{7}
\]

## 4. One empty \(P_3\) gap

The raw count is \(2k(k-1)M\). The omitted internal gap creates a cluster
of two endpoints. There are exactly three merger cases:

| Neighboring endpoints | Number of raw presentations | Result |
|---|---:|---|
| Endpoints of an isolated matching edge | \(2(k-1)(k-2)P\) | Heavy count drops |
| Center and leaf of the \(P_3\) | \(4(k-1)P\) | Leaf merges with its neighboring heavy chain; an all-heavy matching remains |
| The two leaves of the \(P_3\) | \(2(k-1)P\) | Leaves merge to one weight-two endpoint; an all-heavy decorated matching remains |

These counts follow by choosing the adjacent pair of positions, choosing
the remaining \(P_3\) roles, and matching the rest. For example, an adjacent
center/leaf pair has \(2k-2\) positions, two choices of which is the center,
and \(2k-3\) choices of the other leaf, followed by \((2k-5)!!\) matchings.
The unavailable cases at small \(k\) contribute zero.

The list is exhaustive. A nonedge pair can be homogeneous only when its
external defect neighborhoods agree; in this graph that means the two
\(P_3\) leaves. An edge pair either is an isolated matching edge or is
a center/leaf edge. Other pairs are distinguished by their external defect
neighbors, and any vertices in different heavy gaps are distinguished
by an intervening heavy vertex.

After removing these cases, the two-light cluster is transitive and rigid.
The remaining templates are injective and have core order \(4k-2\).
Their count is

\[
 2k(k-1)M-2(k^2-1)P.
 \tag{8}
\]

The merged states that stay in the stratum are already counted below.

## 5. The two all-heavy families

The \(\delta=2\) family has either a four-vertex path or claw, or two
three-vertex paths, with all remaining defect components single edges.
All endpoint weights are one and every gap is heavy. The heavy order
identifies every endpoint separately, so these are injective and pair free.
Their count is

\[
\begin{aligned}
 C
 &=16\binom{2k-2}{4}(2k-7)!!
   +90\binom{2k-2}{6}(2k-9)!!\\
 &=\frac{(k-1)(k-2)(3k-1)}3\,P.
\end{aligned}
\tag{9}
\]

The first term counts the 12 paths and four claws on a chosen four-set.
The second term chooses two unordered triples and one center in each.
Unavailable terms are omitted; the simplified expression is valid also
at \(k=1,2\), where it vanishes.

Finally the \(q=1\) family is a matching on \(2k-2\) endpoints, all gaps
heavy, and either all endpoint weights one or exactly one endpoint
weight two. Its count is

\[
 (1+2k-2)P=M.
 \tag{10}
\]

Both families have core order \(4k-3\). They are disjoint: the fixed heavy
order recovers the defect graph and endpoint weights. The latter family
also absorbs every earlier normalization that stays in the stratum.
At \(k=1\) it is the single heavy transitive core.

An explicit owner rule follows. Assign a unit triangle with one external
defect port to its matching/two-empty presentation, with that port in
the middle, and discard both occupied-\(P_3\) presentations. Assign every
surviving merger to the all-heavy matching family, including its possible
weight-two endpoint. Every remaining pair-free presentation is recovered
uniquely by heavy-gap order and its local cluster structure.

Adding (2), (3), (7), (8), (9), and (10) proves (1). This separates
all raw-presentation collisions by a finite local rule; no generic
symmetry divisor is used.

## 6. The universal third coefficient

The top stratum contributes \(M x^s/(1-x)^{2k+1}\), where
\(s=2k^2+5k+1\). The previous proof gives the entire next stratum:

\[
 \frac M{(1-x)^{2k}}
 \left((2k+1)\sum_{a=1}^k x^{B+a}
             +2k x^B+(k-1)x^{B-1}\right),
 \qquad B=2k^2+4k.
\]

It has \(M(B-1)\) states, with total shift sum

\[
 M\left(B(B-1)+\frac{(2k+1)k(k+1)}2-(k-1)\right).
 \tag{11}
\]

Every state in the deficiency-two stratum contributes leading coefficient
\(1/(2k-2)!\). Lower-heavy strata cannot affect degree \(2k-2\).

For the top term, put \(q=2k\). The second elementary symmetric sum of
\(1-s,\ldots,q-s\), divided by \(q(q-1)\), is

\[
 \frac{s^2}2-\frac{s(q+1)}2+\frac{(q+1)(3q+2)}{24}.
 \tag{12}
\]

Combining (1), (11), and (12), and simplifying, gives

\[
\boxed{\displaystyle
 [n^{2k-2}]\,a_{n,\le k}
 =\frac{k(17-9k)}{2^{k+1}k!}.}
\tag{13}
\]

Here the bracket refers to the eventual counting polynomial, not to a
generating-series coefficient. Together with the two previously proved
coefficients, for \(k\ge2\),

\[
 a_{n,\le k}
 =\frac{n^{2k}-3k n^{2k-1}
          +\tfrac12 k(17-9k)n^{2k-2}}
        {2^k k!}
   +O_k(n^{2k-3}).
\]

For \(k=1\), (13) is the constant 2 in the already exact polynomial
\((n^2-3n+4)/2\).

## Audit scope and controls

The all-budget proof rests on the local cluster classification above.
The incoming checker was independently replayed through \(k=3\), and its
saved \(k=4\) output was inspected. The previously recorded counts match
(1) and the scaled third coefficients are \(4,-1,-15,-38\).

Additional targeted exact normalization probes separated occupied/empty
gap sectors at \(k=1,2,3\). They found no duplicates in the two-occupied
sector; the one-occupied/one-empty full-core counts were \(4,96,1620\);
the empty-\(P_3\)-gap full-core counts were \(6,132\) at \(k=2,3\).
These are controls for the local proof, not its all-budget justification.
The arithmetic simplification of (1), (11), and (12) was independently
checked symbolically.
