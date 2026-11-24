# Hamilton counts and unit rails from compressed Moore observers

The compact observer displays retain every graph label as a source, but omit observers for the labels in a hidden closed neighborhood. The undisplayed induced star is decoded from its internal degrees. The surviving order is seven for a pentagon and sixteen for Petersen. This note counts ordinary directed Hamilton paths and literal arc moves; degree 57 remains conditional.

## The seven- and sixteen-vertex readings

In a pentagon, choose any of the five hidden stars and order its three labels and the two observed labels. Keep precisely the orders for which the prescribed compressed blue spine is directed. Of the 60 possible data, 20 survive, yielding exactly two ordinary classes by exhaustive isomorphism testing. With outgoing rows as zero-based bitmasks, the representatives and counts are:

| Rows | \(H\) | Directed Hamilton cycles | Blue open / closed |
|---|---:|---:|---:|
| \((106,100,25,82,35,76,20)\) | 147 | 12 | 5 / 4 |
| \((114,100,89,19,34,76,24)\) | 135 | 12 | 5 / 4 |

Both are strong and rigid. For Petersen take two-subset labels in lexicographic order, hidden tuple \((0,7,8,9)\), and observed tuple \((2,3,4,5,6,1)\). The sixteen outgoing masks are
\[
(43198,28988,57464,21680,36128,7232,923,19254,37421,9279,63943,62286,58069,49656,33649,746).
\]
The chart is rigid, with \(834,157,149\) Hamilton paths, \(27,494,558\) directed Hamilton cycles and blue split \(19,295\) open / \(26,180\) closed. The [standalone replay](moore-compact-rail-final.py) computes paths by subset/end-point DP, cycles by fixing the start and testing the closing arc, and blue paths by enumerating the signed half-spine compatible with the prescribed anti-automorphism.

These displays differ markedly from full observer hosts or sparse/dense spine charts; graph order and degree alone do not fix Hamilton count in its tournament representations.

## Three distinct strong \(+2\) rails

Mark a single vertex \(v\) and replace it by any nonempty tournament \(M\), with its uniform outside incidences. Let \(P_j(M)\) be the number of ordered covers of \(M\) by \(j\) nonempty Hamilton paths on disjoint subsets. For a fixed exterior, the total height is
\[
 H(T[v\gets M])=\sum_j b_jP_j(M).
\]
To compute \(b_j\), order the exterior vertices into a permutation; the \(j\) occupied gaps are the visits to \(M\). An internal gap is forced if its exterior arc is invalid and the arcs through \(v\) allow a visit, optional if both possibilities are valid, and unusable if neither path through \(v\) nor a required exterior arc is legal. At the endpoints use only the one needed boundary arc. Count choices of occupied gaps by binomial coefficients after the forced gaps. This counts every Hamilton path exactly once by deleting its visits and then filling them in order.

For a transitive module of order \(k\), \(P_j=j!\,{k\brace j}\), so \(P_1=1\), \(P_2=2^k-2\), \(P_3=3^k-3\cdot2^k+3\). The finite width-versus-slot table is triangular with nonzero diagonal; hence a fixed exterior flip gives the same \(+2\) at every positive width iff the difference vector is \(2e_1\). This scalar interpolation makes the following search exhaustive **in its fixed single-port architecture**.

For each of the two compressed pentagon classes, choose at most one preparatory arc reversal, a marked port, and a final arc not incident to the port. Retain only strong hosts on both sides of the final reversal, and require growth (a nonzero \(b_j\) for some \(j\ge2\)); the entire inflation then stays strong. There are no direct growing seams without preparation. The complete marked counts with one preparatory reversal are:

| Lower coeffs \((b_1,b_2,b_3)\) | From 147 chart | From 135 chart | Lower height, transitive width \(k\ge1\) |
|---|---:|---:|---|
| \((109,59,0)\) | 4 | 0 | \(59\cdot2^k-9\) |
| \((133,126,32)\) | 12 | 8 | \(32\cdot3^k+30\cdot2^k-23\) |
| \((133,186,68)\) | 6 | 4 | \(68\cdot3^k-18\cdot2^k-35\) |

All unlisted higher coefficients vanish. Each upper height is exactly two greater, on the same \(k+6\) vertices by reversing its specified exterior arc. Counts in the middle columns are labelled edit/port/arc triples, not numbers of isomorphism classes.

Concrete witnesses all begin at the 147 chart. For the first rail reverse \((0,4)\), mark port 0, and use final flip \((3,4)\) after inflation. For the second and third rails reverse \((1,2)\), use final flip \((2,4)\), and choose port 0 versus port 6, respectively. These last two ports are anti-mates of the original chart and have the same original slot profile \((147,163,50)\): the initial inflated height is \(50\cdot3^k+13\cdot2^k-29\) in either case. Both edits are exterior to both ports, so these are two literal two-edit branches from converse-related inflated source presentations of equal height. At width one the branches coincide at \(133\to135\); at width two their lower heights are 385 and 505. Thus the choice of which mirror mate to expand changes the continuation profile, not just the order or initial count.

The first rail has one simpler exponential scale; the other two carry two independent scales. On \(M_0=3\cdot5\cdot7\cdot11\cdot13=15015\), each is periodic in \(k\ge1\) with a common period 60. Exact enumeration of those 60 phases gives 60, 48, and 44 distinct lower residues, with a single overlap between the last two, hence 151 lower residues. Adding two gives 151 other residues disjoint from all lower residues. This is a literal modular selection gain, but a fixed modulus and three sparse families do not give absolute selection of all sufficiently large odds.

## What a degree-57 pattern would need to add

Every diameter-two Moore graph of degree at least two contains a pentagon. Consequently these local pentagon-derived rails are already available without a degree-57 solution. To make degree57 integral to a Chan–Pak proof would require a genuinely global mechanism, not another instance of this same local motif.

There are two exact barriers. First, a fixed finite bank of bounded-order cores with independent transitive port widths has only polylogarithmically many distinct counts at most \(X\). After ignoring transitive ordinal padding, each relevant core vertex belongs to a directed triangle; a module of width \(w\) at such a port contains an induced triangle inflation of count \(2^w+1\), which injects into paths of the ambient tournament by vertex insertion. Thus each relevant width is at most \(\log_2X\), so finitely many ports cannot cover positive density.

Second, the present difference for an arbitrary filler \(M\) is really \(2P_1(M)=2H(M)\), not a unit after arbitrary nontransitive substitution. More generally, for any fixed set of simultaneously marked ports the same exterior-permutation proof gives a polynomial in the products \(\prod_i P_{j_i}(M_i)\). If a fixed exterior edit gives \(+2\) for **all independent transitive widths**, tensoring the triangular interpolation shows its difference on arbitrary fillers is exactly \(2\prod_i H(M_i)\). This leaves a precise constructive question: can a global Moore57 transport system produce a growing *family of compatible marked states* with coupled continuation profiles or adaptive port assignments that hit prescribed absolute integers, rather than independent width rails? The anti-mate branching above is a small concrete example of distinct continuation profiles, not such a selector.
