# A matrix code without a symmetry hypothesis; the projective trade stars

25 November 2026. Continuation of [the midpoint and Moore parity note](MIDPOINT-MOORE-PARITY.md) and the earlier [LCD projection](../moore_odd_count_lcd_projection_20261124/README.md). The compression assumes only that there is a degree-57 Moore graph. The final geometric result separately assumes that a 400-coclique gives the point-line design of \(\mathrm{PG}(3,7)\). No graph is constructed.

## 1. A principal block decodes the entire graph

For \(A\) a degree-57 Moore adjacency, the projector onto the \(-8\) eigenspace over \(\mathbb Q\) is \((A^2-64A+399I)/975\). Its characteristic polynomial is \(t^{1730}(t-1)^{1520}\), and both entries and identities reduce at 2 because 975 is odd. Over \(\mathbb F_2\) this gives the symmetric idempotent
\[
P=I+J+A,\qquad \operatorname{rank}P=1520.
\]
Choose a maximal independent row set \(S\) of \(P\). Its principal minor \(P_{SS}\) is nonsingular: the row vectors span the nondegenerate image of the orthogonal projector, and their Gram matrix is the principal minor since \(P^2=P=P^{\mathsf T}\). Then
\[
 P=P_{:S}(P_{SS})^{-1}P_{S:}.
\]
The graph data touching \(S\) suffice. The graph can be recovered up to isomorphism by putting \(S\) first and using the resulting labelled representative. Preserving the original site labels adds \(\lceil\log_2\binom{3250}{1520}\rceil=3235\) bits for the original basis subset, with each side ordered by its original labels. The replay recovers the full matrix in Petersen with rank 4 and Hoffman–Singleton with the complementary projector of rank 22.

This is a projection code, not a reading of a single path count. The binary image also preserves the graph's permutation symmetries as a code: a coordinate permutation preserving the image preserves its unique orthogonal projector, and hence \(A\); the converse is immediate. The marked tournament format itself need not preserve those automorphisms.

## 2. Direct coding versus enumerative coding

Let \(i=e(S)\) and \(j\) count edges across to the complement. As \(A+8I\succeq (65/3250)J\), testing the indicator of \(S\) gives
\[
 2i+8\cdot1520\ge (65/3250)\,1520^2,\quad
 i\ge17024,\qquad j=86640-2i.
\]
The number of possible records of the indicated type is at most the sum over \(i\ge17024\) of \(\binom{1154440}{i}\binom{2629600}{j}\). The ratio of consecutive admissible terms is
\[
 \frac{1154440-i}{i+1}\frac{j(j-1)}{(2629601-j)(2629602-j)},
\]
which decreases with \(i\) and is already below \(1/32\) at the first term. The geometric upper bound is strictly below \(2^{499859}\) by exact integer comparison. To rank or unrank the record, first find the offset of \(i\) using the consecutive binomial products, then use two ordinary subset ranks at those weights. The proof needs no orbit action, geometric design or maximal coclique.

There is a smaller **enumerative** variant. Let \(d_v\) denote the degrees inside \(S\), and \(w(d)=\binom{1730}{57-d}\) for \(0\le d\le57\), zero otherwise. Given an internal graph the number of crossing records is \(\prod_v w(d_v)\). Multiply by \(36^{i-17024}\ge1\), and sum over all binary internal edges. Each edge bit appears at just its two endpoints. Summing it by Cauchy–Schwarz replaces both incident functions by their norm in that bit. Iterating eliminates each internal bit and produces
\[
 \frac{\left(\sum_{d=0}^{57}\binom{1519}{d}\,36^d\,w(d)^2\right)^{760}}{36^{17024}} < 2^{493538}.
\]
The relaxation deliberately ignores the degrees of complement vertices and whether a record completes to a Moore graph. Rank by finite enumeration of the relaxed records to get an exact but potentially impractical decoder; the displayed inequality alone is not an efficient ranking algorithm.

| Code | Bits | Marked tournament | With fixed directed spine | Spine and original graph labels |
|---|---:|---:|---:|---:|
| Two subset ranks | 499859 | 1001 | 1002 | 1005 |
| Complete enumeration | 493538 | 995 | 996 | 999 |

Each available arc carries one bit relative to a fixed contact order; the site ordering and reserved arcs are part of the marked format. The larger code permits direct subset ranking; neither code asserts the necessary bit count for actual Moore graphs. The much smaller earlier [249-site order-11 code](../moore_fixed_order_spectra_20261124/PENTAGON-PORT-ATLAS.md) uses a different, additional symmetry assumption.

## 3. A smaller certified Hamilton-count window for the spined codes

The 1002-site payload leaves 641 nonspine contacts. Reserve the closing arcs of 334 disjoint consecutive triples and orient each triple cyclically; reserve a further nonspine arc as an ignored flip marker. The 996-site enumerative version similarly fits 332 triangles and a marker. An arbitrary payload and either marker value preserve the planted triangles. The [odd-cycle packing identity](https://arxiv.org/abs/2307.05569) gives at least \(3^k\) paths from \(k\) disjoint directed triangles, counting each subcollection with weight \(2^{|\text{subcollection}|}\).

For an upper bound, add a site with two-way incidence to all \(n\) sites and take the permanent of the resulting directed adjacency matrix. Each Hamiltonian path extends to a distinct cycle cover, so Brégman's inequality bounds it by \((n!)^{1/n}\prod_i ((d_i+1)!)^{1/(d_i+1)}\), where \(d_i\) are the outdegrees. This is the standard permanent technique in [Alon's Hamilton-path bound](https://www.tau.ac.il/~nogaa/PDFS/hamilton.pdf). Here the maximum of this bound is at balanced integer row sums \(d_i+1\). For completeness, let \(g(r)=\log(r!)/r\), \(D=\log(r+1)-g(r)\). Discrete concavity amounts to \(2D\ge (r+1)\log(1+1/(r+1))\). The right side is less than one and \(D\ge\tfrac12\) by \(\log((r+1)/i)\ge1-i/(r+1)\) averaged over \(1\le i\le r\). For even \(n=2m\) the upper bound is therefore
\[
 (m!)\,((m+1)!)^{m/(m+1)}(n!)^{1/n}.
\]
In these two cases exact power comparisons give \(((m+1)!)^{1/(m+1)}>2^{7.5}\), \((n!)^{1/n}<2^{8.6}\). Multiplying compares the bound to \(2^{1.1}m!(m+1)!\), which the replay checks by tenth powers. The certified windows are
\[
 3^{334}\le H<2^{7563}\quad(n=1002),
 \qquad
 3^{332}\le H<2^{7509}\quad(n=996).
\]
This improves the weaker disjoint-six-block estimate \(2^{7882}\) for the first format. Both marker values decode to the same graph; **no jump or target Hamilton count is promised**.

## 4. A projective center has exactly 400 leaves in the closest-coclique graph

Assume now the separate projective-space branch. Regard the original 400-coclique \(C\) as the points of \(\mathrm{PG}(3,7)\), and the other graph vertices as its 2850 lines. The [previous trade lemma](MIDPOINT-MOORE-PARITY.md#4-reciprocal-spreads-and-the-plane-trade) exchanges the 57 points of each plane \(\Pi\) for its 57 line vertices. Call the new 400-coclique \(C_\Pi\). There are 400 plane choices. For two distinct planes, their old point sets meet in an eight-point line and their exchanged line sets meet in the one line of intersection. Thus
\[
 |C_\Pi\cap C_\Sigma|=400-(57+57-8)+1=295 .
\]

We can show more than a change of design isomorphism type: **the traded design on \(C_\Pi\) has exactly one projective subplane of order seven**.

- Its points consist of 343 old affine points and the 57 new infinity points, which are old plane lines. A block given by an old point of \(\Pi\) consists entirely of new infinity points; those blocks form the dual plane. A block given by an old exterior line consists of its seven old affine points plus one new point \(b\), the unique plane-line in the reciprocal spread of that exterior line.
- For a fixed old direction \(p\in\Pi\), the 49 exterior lines in that direction get the 49 new points \(b\) not incident with \(p\), each exactly once. Indeed the reciprocal spread at \(b\) contains \(b\) and one outside line at each of the 49 points of \(\Pi\) off \(b\). Therefore two old affine parallels remain disjoint in the traded design: they receive different new points.
- If a projective subplane in the traded design contains an old affine point \(x\), its eight lines through \(x\) each have one new point, all distinct. Every other new point in the subplane must also lie on one of those eight lines, so the subplane has exactly eight new points and 49 old affine points. Its affine part is closed under whole old affine lines, by the two-point uniqueness of each block. Over \(\mathbb F_7\) such a set is an affine flat: after translating to contain zero, it is closed under scaling and under addition by taking midpoints and doubling. Having 49 points it is an affine plane, containing old parallel lines. These remain disjoint after the trade, impossible in a projective subplane. A subplane wholly in the new points is the dual plane already exhibited.

In the graph on maximum cocliques with an edge for overlap 343, an edge is exactly an order-seven plane trade. The original projective design has exactly 400 such subplanes: line closure over the prime field \(\mathbb F_7\) forces an ordinary projective flat, and a subplane of that order spans a plane. Explicitly, any three noncollinear points of the subplane span an ordinary plane; closure gives their triangle sides, and any point in that ordinary plane is on a join of a triangle vertex to a point on the opposite side. Every traded coclique has degree one by the preceding proof, with its one reverse trade returning to \(C\). Hence the whole connected component is the **star \(K_{1,400}\)**. This is stronger than merely having 401 maximum cocliques, or proving the 400 traded designs are not projective. In particular no graph automorphism can send the center to any leaf.

The separate [order-11 and free order-13 observations](MIDPOINT-MOORE-PARITY.md#4-reciprocal-spreads-and-the-plane-trade) give a useful conditional multiplier: if either action occurs together with a projective 400-coclique, its orbit consists of 11 or 13 distinct projective centers. Their star components are disjoint. This forces at least \(4411\) or \(5213\) maximum cocliques, respectively. None of these hypotheses is inferred from the degree alone.

Replay: **node projector-compression-check.js**. It verifies the integer bounds, code capacities, and small-graph projector completions. The star is an incidence proof, not a search over unknown completions.
