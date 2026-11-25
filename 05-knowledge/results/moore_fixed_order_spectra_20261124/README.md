# Moore representations versus their full fixed-order Hamilton spectra

[Follow-up: exact 14- and 50-site Hoffman–Singleton phase codes and longer fixed-order prefixes](EXACT-PHASE-CHARTS.md).

\(H(T)\) means the number of directed Hamiltonian paths in a tournament, not the number of undirected paths in the graph being encoded. Write \(S_n=\{H(T):|V(T)|=n\}\). Counts of classes below are **ordinary isomorphism classes**, not labelled tournaments, not classes identified with converse, and not marked paths. Adopt the nontrivial diameter-two convention for the main comparison, and include \(K_2\) as the diameter-one boundary check.

## Representation and placement

| Moore graph | Tournament order and encoding | Hamilton-path placement |
|---|---|---|
| \(K_2\), boundary case | 2, the unique arc, marked spine is the graph edge | \(H=1\); entire \(S_2\) |
| \(C_5\) | 5, reverse only the nonspine graph chords along a graph spine, or reverse its tile mate | \(H=9,11\); the fourth and fifth of seven heights. A different compact circulant encoding has \(H=15\), the maximum |
| Petersen | 10, graph spine or compact affine encoding | Graph-spine representatives give \(5213,6961\), with 2127 classes at **each** height. The compact encoding gives \(8265\), with 2108 classes at that height |
| Hoffman–Singleton | 50, compact affine ordinary class, or native marked graph spine | Exact chart and certified placement below: the compact count is below the mean of uniform labelled tournaments, and both native counts lie below 53% of that mean. Neither the full counts nor \(S_{50}\) were enumerated |
| degree 57 | 3250, **conditional** native marked complement spine | Conditionally below \(2^{-970}\) times the order-3250 labelled mean; graph, height, and whole \(S_{3250}\) remain unresolved |

To form a native chart from an order along a graph Hamiltonian path, orient the displayed spine forward, orient all other pairs forward except for graph chords, and retain the spine as the decoder tag. To decode, add the spine edges to the backward pairs. The corresponding rule on a Hamiltonian path of the *complement* uses all graph edges as backward pairs and does not add the spine. A tournament isomorphism class together with its marked spanning path specifies a layout; the bare counts do not identify the intended decoder. The compact affine Petersen and Hoffman–Singleton classes have a decoder on the **unmarked** tournament (see [compact rule and decoder](../moore57_compact_blue_routes_20261124/COMPACT.md) and [Python](../moore57_compact_blue_routes_20261124/compact.py)). No padding is needed for these known graphs.

In the row-by-row lexicographic pair convention (bit 1 when the arc goes backwards), mask 8 represents the native \(C_5\) graph-spine class; optional-tile complementation gives mask 358 on the same spine. The circulant \(i\to j\) iff \(j-i\in\{1,2\}\pmod 5\) represents the compact one. For native Petersen, take pentagon vertices 0–4, pentagram vertices 5–9, and join mates \(i,i+5\). The all-spine classification in [the native ledger](../moore_merged_count_addresses_20261124/README.md) applies. For native Hoffman–Singleton, use [the explicit affine graph spine](../moore_chanpak_counts_20261124/counts.py), and for the compact class use the side-swapping, self-converse tournament from [field_host(5)](../moore57_compact_blue_routes_20261124/compact.py). The existence of a native degree-57 complement spine is conditional only on the graph itself: the complement has minimum degree \(3192>3250/2\).

The displayed small masks can also be taken from the exhaustive spectrum CSV using their count and the old canonical ledger, not inferred from height alone. The table does not assert uniqueness of any representation over every possible graph decoder.

## Exact spectra through order ten

| order | ordinary classes | distinct heights | highest \(H\) | first omitted odd other than 7 or 21 *below* the maximum |
|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | — |
| 3 | 2 | 2 | 3 | — |
| 4 | 4 | 3 | 5 | — |
| 5 | 12 | 7 | 15 | — |
| 6 | 56 | 19 | 45 | 35 |
| 7 | 456 | 77 | 189 | 63 |
| 8 | 6,880 | 320 | 661 | 611 |
| 9 | 191,536 | 1,520 | 3,357 | 2,883 |
| 10 | 9,733,056 | 7,663 | 15,745 | 14,685 |

In particular
\[
S_5=\{1,3,5,9,11,13,15\},
\quad
\text{with multiplicities }(1,3,2,2,1,1,2).
\]
The compact \(C_5\) class is one of two maximizers: the other has sorted scores \(1,2,2,2,3\), not \(2,2,2,2,2\). Neither 9 nor 11 first appears at a smaller order.

The full [order-ten ledger](order-10.csv) consists of every odd from 1 through **14,683** other than 7 and 21, and 323 additional odd heights from 14,685 through 15,745; there is nothing above 15,745. Its rows give multiplicity and a certifying mask. The nine-vertex maximum 3357 means that every native Petersen height in the previous source and mate families is first attained at order ten, even when other unrelated tournaments share that height. In the following table rank counts distinct heights; the cumulative class share counts all ordinary classes at or below \(H\).

| Petersen example | \(H\) | class multiplicity | rank in 7663 heights | cumulative class share |
|---|---:|---:|---:|---:|
| native graph-spine | 5213 | 2127 | 2605 | 29.55% |
| marked +2 seam, source | 5625 | 2139 | 2811 | 33.91% |
| marked +2 seam, mate | 5627 | 2176 | 2812 | 33.94% |
| graph-spine tile mate | 6961 | 2127 | 3479 | 48.78% |
| compact affine | 8265 | 2108 | 4131 | 63.45% |
| exceptional self-converse native mate | 15549 | 11 | 7636 | 99.9991% |
| highest native complement-spine mate | 15571 | 4 | 7644 | 99.9995% |

Thus the same graph can sit in an ordinary part of the spectrum or near its extreme tail, with only four classes at its highest native mate height. The 333 complement-spine source heights and 362 mate heights overlap in 29 **scalar** values, even though the merged source and mate nodes are disjoint. No source height enters the incomplete fringe above 14,683; 29 mate heights do. The previously identified 5625/5627 move is a tile complement, not a single-arc reversal.

### The 208 extra holes are all fixed at order eleven

For each of the 208 odd heights missing from \(S_{10}\) below its maximum other than 7 and 21, [the fill ledger](order-11-fill.csv) gives an explicit 11-site tournament. All masks are independently checked by subset DP. Thus every such hole is **exactly** an order-ten phenomenon, and every odd value through 15,745 other than 7 and 21 belongs to all \(S_n\) for \(n\ge11\): append a new sink, which must be last, to raise the order without altering \(H\). This assertion is a constructive band, not an exact enumeration of \(S_{11}\) or a claim about larger targets. The standing all-orders question is [Chan–Pak's Hamilton-path spectrum conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).

### Enumeration mechanism and independent replay

[McKay's nauty generator](https://users.cecs.anu.edu.au/~bdm/data/generators.html) supplies one tournament per ordinary class. Use its default (noncanonical, upper-triangle text) output; it keeps the same last-vertex deletion core together. For a fixed \((n-1)\)-site core let \(P(S,i)\) count paths on \(S\) ending at \(i\), and \(R(S,j)\) count those starting at \(j\). Compute both by endpoint subset DP, then cache
\[
Q_{ij}=\sum_{S\ni i,\ S\not\ni j}P(S,i)R(\bar S,j).
\]
With apex outcontact bits \(x_i\), paths either put the apex first, last, or between two nonempty parts, giving
\[
H(x)=\sum_i x_i R(V,i)+\sum_i(1-x_i)P(V,i)+\sum_{ij}(1-x_i)x_j Q_{ij}.
\]
All apex words for a core can be evaluated by quadratic subset increments. Order ten used 191,536 cached cores and counted 9,733,056 classes. A separate implementation checked one labelled witness for every height listed at orders 2–10 and all 208 order-eleven witnesses using the unseparated endpoint recurrence. This checks their counts independently, not just the apex formula. The complete census relies on the isomorph-free generation.

## Where the large native charts sit

For the concrete Hoffman–Singleton charts, a whole-cohort DP gives lower bounds by keeping only paths that visit each displayed 5-site cohort contiguously. These are not exact \(H\)'s. I combined them with rigorous upper bounds by adding a universal bidirectional apex, making each Hamilton path a single-cycle permutation. An upper bound on the permanent then bounds the number of single cycles.

A useful sharpening over plain Brégman: expand a row with few possible columns, remove the row and chosen column, and track the directed paths of prescribed entries. Never allow an arc that closes a component before the cycle is spanning, and delete the arc closing the newly combined path from the remaining matrix. At each leaf upper-bound the remaining count by the product of \((r!)^{1/r}\) for its row sums. This tree is an upper bound even though the leaf counts still include disconnected cycle covers. The replay uses exact dyadic upper multiplicative factors obtained by integer-power comparisons and a conservative final IEEE rounding allowance.

| 50-site chart | verified lower bound | upper bound, as a strict power of two | upper as share of labelled mean |
|---|---:|---:|---:|
| compact affine \( \texttt{field\_host(5)} \) | 887,480,820,590,625 | \(2^{164.878}\) | \(<0.80\) |
| native graph spine | 2,496,220,189,141 | \(2^{164.221}\) | \(<0.51\) |
| its tile mate | 18,811,695,423,864 | \(2^{164.285}\) | \(<0.53\) |

The replay evaluates the **self-converse field host** used in the table directly. Its uniform cyclic translation acts freely on ordered paths, and reflection gives \(H\equiv1\pmod4\) for even self-converse tournaments by the repo's odd-cycle argument, so its count is \(5\bmod20\).

### An exact variance formula makes this a percentile bound

For uniformly random **labelled** tournaments put \(\mu_n=n!/2^{n-1}\). Fix a first ordered path and choose a uniformly random second permutation. For an edge from the first path, the contribution to the joint probability ratio across the two paths is 2 if reused forward, 0 if reused backward, 1 otherwise. Write it as \(1+I_+-I_-\), and expand the product along the first path. For a chosen set of edges, split it into runs of consecutive first-path edges. A run of length \(m\) can occur in the second permutation only as an oriented block, forward or reversed, and its aggregate weight is \(1+(-1)^m\). Only runs of even length survive. Contracting all \(2j\) selected edges gives \((n-2j)!\) placements for each direction assignment. Therefore
\[
\frac{\operatorname{Var}(H)}{\mu_n^2}=
\sum_{j=1}^{\lfloor(n-1)/2\rfloor}
\frac{(n-2j)!}{n!}
\sum_{r=1}^{\min(j,n-2j)}
2^r\binom{j-1}{r-1}\binom{n-2j}{r}.
\]
The binomial factors split \(j\) into \(r\) even-length runs and place them on the first path with at least one missing edge between successive runs. At 50 the exact sum is about 0.03995937. The one-sided variance inequality now shows that more than **50%**, **85%**, and **84%** of labelled 50-site tournaments, respectively, have greater \(H\) than the three displayed charts. This is not a statement about percentiles over *distinct values* or ordinary classes.

For the conditional degree-57 sparse chart, \(r_i\le\min(n,n-i+57)\) at site \(i=0,\dots,n-1\) after adding the apex, whose own row sum is \(n\). The orientation starts transitive and reverses at most 57 arcs incident to each site. Brégman's product for these upper row sums is about \(2^{28988.73}\) compared with \(\mu_n\approx2^{29984.63}\). A coarse rigorous gap exceeds 970 bits: use
\(\log(r!)\le (r+\tfrac12)\log r-r+1\), from trapezoids for the concave log, and \(\log(57!)\ge57\log57-56\). With \(n=3250\),
\[
\log(\mu_n/U)\ge
 n+1-(n-1)\log 2+\log(57!)-58\log n
 -(1+\tfrac12\log n)\bigl(59/n+\log(n/57)\bigr)>970\log 2.
\]
For the last bound one can use \(\log2<.694,\log n<8.1,\log57>4\) and \(\log(n/57)<4.1\). Thus unlike the upper-tail Petersen mate, the native **sparse** degree-57 representation would be far below the ambient mean on this scale, independently of the graph or complement spine. This does not bound its dense tile mate. See the [earlier exact count-address construction](../moore_merged_count_addresses_20261124/README.md) for an alternative, larger HS encoding with a fully evaluated count and a conditional degree-57 analogue.

## Replay

Build nauty 2.9.3's gentourng. From this directory, compile spectrum.cpp, fill.cpp, verify.cpp and hs-bracket.cpp with C++17 and optimization. For n=2,..,10, use gentourng -q n piped to spectrum n, and compare order-n.csv. Run fill order-10.csv; verify n < order-n.csv (or verify 11 fill < order-11-fill.csv). The last column of the fill ledger gives the zero-based PRNG attempt at which the mask was found, not a class multiplicity. Run hs-replay.py after compiling hs-bracket; it loads the older chart definitions and computes the integer dyadic factors, count bounds and variance thresholds. Compile petersen-spine-heights.cpp for the extra source/mate height lists. The boundedness tests are not a substitute for the unresolved full spectra at 50 and 3250.
