# How many actual arc edits does the Kaprekar cycle need?

**PROVED scoped / FINITE-EXACT.** The [decimal Kaprekar quotient](KAPREKAR-IMAGES.md) cycles through \(09,81,63,27,45\). All five integers occur as ordinary Hamilton counts. Requiring the steps to be actual arc reversals is stricter than requiring the vertices of the numerical cycle to be attainable.

## The sharp physical answer

The shortest arc walk visiting the five heights in the displayed order and returning **up to isomorphism** takes **seven edits**. It takes **eight edits** to return to the same labelled orientation. The following witnesses meet these bounds with only genuine off-cycle detour heights, not by waiting at the same numerical phase.

| Closure | Vertices | Actual height walk |
|---|---:|---|
| Up to relabelling | 9 | \(9\to81\to115\to63\to37\to27\to45\to9\) |
| Exactly labelled | 8 | \(9\to81\to83\to109\to63\to73\to27\to45\to9\) |

For the first row start with outgoing masks \((486,484,488,499,487,448,384,256,0)\) and reverse
\[
(2,8),(2,4),(7,8),(6,8),(2,4),(6,8),(2,8)
\]
in order, using zero-based vertices. The final orientation differs only on edge \((7,8)\); those two vertices were terminal singleton pads, so swapping their labels is an isomorphism. For the second row start with masks \((250,252,249,208,224,200,128,0)\) and reverse
\[
(0,7),(0,3),(0,6),(3,5),(0,3),(0,7),(3,5),(0,6).
\]
Each arc occurs twice. Both walks are replayed by subset DP. Eight is the minimum vertex order for a full walk (height 63 already requires eight). The separate order-eight seven-edit finite test rules out realizing the seven-edit optimum there, so the shortest relabelled walk requires nine vertices.

## Why there must be two detours, at every order

We use the [strong-order bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), the SCC product law and the previously established exclusion of counts \(7,21\). Any tournament of height \(63=3^2\cdot7\) has exactly one nontrivial SCC, because splitting off any nontrivial factor would require an SCC of count 7 or 21. The strong-order bound places this core at order at most eight. An exhaustive order-seven Hamilton-framed scan gives no 63 (and ordinal singleton padding handles smaller orders), so the core has exactly eight vertices. In every framed scan all orientations of the nonconsecutive pairs are tested, and every tournament admits a Hamilton path to use as a frame.

A reversal taking a 63 tournament to 81 either edits this core internally or merges it with exactly one singleton:
- Internally, a split of height 81 would need factors \(3\cdot27\), \(9\cdot9\), \(3\cdot3\cdot9\), or \(3^4\). Their minimum vertex orders, using the maxima at orders four and five, give at least nine vertices in every case. Thus an internal result remains strong.
- Across components, a reversal merges the interval between its endpoint SCCs when there is a retained forward path between them. If that interval misses the original core, 63 persists as a factor and cannot give 81. If it includes the core, it merges at least one singleton; height 81 cannot have a strong core larger than nine by the same bound.

Thus the 81 state has a single strong core on eight or nine vertices. Any reversal to 9 would have to be internal to that core: otherwise the old strong core stays a component or remains an induced subtournament of the merged component. Induced-subtournament counts are lower bounds on ambient counts, by inserting missing vertices into each directed path while preserving its old vertex order. Consequently \(9\to81\to63\) at arbitrary order would restrict to the same two steps on the 81 core.

The finite test starts at height 9: either ordinal pads around one of the two strong order-five height-9 cores, or two directed triangles with pads. For the two orders those are 14 and 20 representative shapes. Reversing every arc gives 12 and 13 strong 81 arrivals, respectively; checking every next arc gives **no** 63. Hence at least one detour is needed across the first two numerical arrows. The short path \(9\to81\to115\to63\) from the first witness shows that one extra edit suffices for that segment at arbitrary order.

Independently there is **no** \(63\leftrightarrow27\) arc. Any downward flip from 63 is internal by the same monotonicity argument. Enumerate height 27 on eight vertices: padded strong cores of order six or seven, or a 3-core and a 9-core in either ordinal order. Other decompositions do not fit, and a larger strong 27 core is precluded by the strong bound. The exact replay checks all 220 such presentations (may include repeated isomorphism classes), with zero reversals to 63. Thus the second detour is necessary in a disjoint segment.

Five prescribed steps plus these two independent detours give the global lower bound seven. A labelled closed walk in the orientation cube has even length because every edit changes exactly one orientation bit, so a labelled loop needs at least eight. The argument does not depend on any closure convention to rule out the direct five-edge schedule.

Replays: `python3 kaprekar-metagraph-audit.py` for the two all-order reductions; `python3 kaprekar-cycle-check.py` for the literal loops; `c++ -O3 -std=c++20 kaprekar-cycle-exact.cpp -o /tmp/kap-cycle && /tmp/kap-cycle 8 2` for the order-eight, seven-edit nonexistence test (second argument is the number of extra edits beyond five). It allocates the extra steps in every possible gap, starts from the complete set of height-9 shapes, and checks closure up to relabelling; the result is zero. The [base-16 comparison](BASE16-OBSTRUCTION.md) fails already at its first individual jump.
