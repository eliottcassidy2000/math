# Count-fiber class series through 73 and the three-bit Collatz lens

This extends the [series through 25](../chan_pak_count_fiber_grammar_and_hidden_cycle_trades_20261126/README.md), the [9/17 rooted and uniform readouts](../chanpak_uniform_and_rooted_readout_equivalences_20261126/README.md), and earlier [Collatz parity-word work](../collatz_weighted_passages_first_crossing_20261126/README.md). It gives an exact local/global obstruction for labelled cyclic triples, not a claim about the long-run behavior of Collatz or a proof of global Chan–Pak.

## 1. Complete class series for every odd height 27 through 73

Let \(t\) mark vertices and \(R=1-t\). In the table a row denotes \(A/R^2+B/R^3+C/R^4\), an ordinary generating function for tournament **isomorphism classes**, not labelled orientations or the quotient by converse.

| \(h\) | \(A\): one strong component | \(B\): two | \(C\): three |
|---:|---|---|---|
|27| \(2t^6+3t^7\) | \(4t^8\) | \(t^9\) |
|29| \(4t^6+2t^7\) | 0 | 0 |
|31| \(2t^6+2t^7\) | 0 | 0 |
|33| \(5t^6+3t^7\) | \(2t^8\) | 0 |
|35| \(2t^7\) | 0 | 0 |
|37| \(5t^6+2t^7\) | 0 | 0 |
|39| \(2t^7\) | \(2t^8\) | 0 |
|41| \(t^6+5t^7\) | 0 | 0 |
|43| \(2t^6+2t^7\) | 0 | 0 |
|45| \(2t^6+5t^7+5t^8\) | \(4t^8+8t^9\) | \(3t^{10}\) |
|47| \(2t^7\) | 0 | 0 |
|49| \(7t^7+t^8\) | 0 | 0 |
|51| \(8t^7+4t^8\) | \(4t^9\) | 0 |
|53| \(4t^7+6t^8\) | 0 | 0 |
|55| \(4t^7\) | \(2t^9\) | 0 |
|57| \(4t^7+3t^8\) | \(4t^9\) | 0 |
|59| \(8t^7+4t^8\) | 0 | 0 |
|61| \(t^7+4t^8\) | 0 | 0 |
|63| \(2t^8\) | 0 | 0 |
|65| \(11t^7+3t^8\) | \(2t^9\) | 0 |
|67| \(6t^7\) | 0 | 0 |
|69| \(12t^7+6t^8\) | \(8t^9\) | 0 |
|71| \(7t^7+4t^8\) | 0 | 0 |
|73| \(8t^7+11t^8\) | 0 | 0 |

**Certification and recursion.** [Busch's exact strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) is 75 at order nine, and the lower bound \(5^{(n-1)/3}\) exceeds 73 at every larger order. Camion's theorem allows fixing a directed spanning cycle to enumerate strong types through order eight by its free chords. The replay computes \(H\) independently by endpoint subset DP and odd-cycle packs for every retained canonical class. It retains 207 strong classes of \(H\le73\) at orders 3 through 8: respectively 1, 1, 6, 35, 111, 53. The five entries of each pack vector are \(a_3,a_5,a_7,d_{33},d_{35}\), counts of cycles and the indicated disjoint pairs; no larger pack fits. The raw rows give order, height, chord reversal mask, automorphism order, self-converse flag and that vector.

If \(G_h(t)\) is the enumerator of strong cores of height \(h>1\), use \(W_{1,0}=1,\ W_{h,0}=0\) for \(h>1\), and
\[
W_{h,j}=\sum_{d\mid h,\ d>1}G_d\,W_{h/d,j-1}.
\]
The denominator \(R^{j+1}\) inserts singleton chains into the \(j+1\) gaps of each ordered core word. The \(t^9/R^4\) term is literally the three-component word of three directed triangles, while \(2t^6/R^2\) comes from two strong converse classes whose automorphism groups have order 3. Same count, different roles for the prime 3. At 45, the two-component terms come separately from \(3\cdot15\) and \(5\cdot9\), while the three-component term comes from the three orders of \(3,3,5\).

## 2. New ways in which equal counts hide geometry

Masks below reverse chords in lexicographic order from a directed spanning cycle \(0\to1\to\cdots\to n-1\to0\), other chords forward. Equalities of the full path-cover polynomial use the normalization \(D_T(x)=1+2a_3x+(2a_5+4d_{33})x^2+(2a_7+4d_{35})x^3\) for \(n\le8\), from the [Rédei–Berge expansion](https://arxiv.org/abs/2307.05569).

- **27, order seven.** Masks 1 and 23 exchange two 5-cycles for a pair of disjoint triangles, preserving \(D=1+10x+14x^2+2x^3\), extending the previous trade at 17. One is not self-converse and the other is.
- **29, order six.** Masks 21 and 40 have the same full unlabelled pack profile \((6,6,0,1,0)\). The first is self-converse and the other is not. Even retaining the pack species does not recover whether a class is paired in the metagraph.
- **33, order six.** Masks 12 and 22 have \(D=1+14x+18x^2\), with 9 versus 7 five-cycles and respectively 0 versus 1 disjoint triangle pair. They have different automorphism orders, 3 versus 1, and hence 11 versus 33 tilings in their classes. By the earlier uniform-host theorem, no homogeneous one-port observation can distinguish them. This upgrades the earlier 15 example, where an automorphism jump used a source-padded witness: here both tournaments are strong.
- **61, 63, 65.** The first and last have witnesses already at order seven; 63 first appears at eight. Although \(63=3^2\cdot7\) is composite, every witness must have one nontrivial strong component: every nontrivial product factorization using only smaller odd numbers would involve 7 or 21. Its two order-eight witnesses (masks 6 and 19) are a rigid converse pair with identical profile \((8,17,6,0,0)\), so the two classes collapse to one merged node with 126 tilings. By contrast, \(65=5\cdot13\) has two ordered-join cores of size nine in addition to its single-component cores. This is a concrete way the gap at 7 still affects later heights without excluding all its multiples.

## 3. The labelled triangle predicate in Collatz parity coordinates

Use the halved Collatz map \(T(x)=x/2\) for even \(x\), and \(T(x)=(3x+1)/2\) for odd \(x\). Each \(k\)-bit itinerary occurs on exactly one residue class of \(x\bmod 2^k\). To see this inductively, two inputs differing by \(2^k\) agree on the first \(k\) bits. Along the common branch with \(a\) odd steps their \(k\)-th iterates differ by \(3^a\), which is odd, so the next itinerary bit distinguishes the two lifts.

Assign the \(\binom n2\) parity bits to the unordered edges in any fixed order, with \(b_{uv}\) for \(u<v\) and \(b_{vu}=1-b_{uv}\). This bijects the \(n\)-vertex **labelled orientation cube** with the residue classes modulo \(2^{\binom n2}\). Ordering contacts by successive largest vertex gives a triangular filtration: extending the tournament by one vertex is exactly lifting the parity cylinder by another \(n-1\) bits. Fixing a Hamiltonian spine is just prescribing its corresponding bits. The orientation status of a labelled triple \(u<v<w\) is read by
\[
\mathrm{cyclic}\iff b_{uv}=b_{vw}=1-b_{uw}.
\]
If written in cyclic order, the three bits are monochromatic. For a consecutive three-step parity word this is precisely the cylinder with local input \(0\) or \(7\bmod8\) (words 000 and 111). A triple in the full tournament generally reads three *nonconsecutive* bit positions with one sign reversed. The local \(0,7\) rule cannot be substituted for that signed contact routing. Reversal of all arrows preserves the monochromatic test, just as complementing all three queried bits does.

**Exact local/global test, including counts.** Treat the edge bits as Boolean variables. Every required cyclic triple contributes two signed equalities over \(\mathbf F_2\). Forbidden cyclic triples are signed not-all-equal conditions. For a required set \(A\), let \(L(A)\) be 0 if the signed equations are inconsistent and \(2^c\) otherwise, where \(c\) counts free connected components in the equality system, including isolated edge variables. Fixed spine contacts pin their components and may create a contradiction. For required triples \(A\) and forbidden triples \(B\), the number of compatible labelled tournaments, hence compatible Collatz cylinders, is
\[
\sum_{J\subseteq B}(-1)^{|J|}L(A\cup J).
\]
This is signed union-find plus inclusion–exclusion. It counts realizations of marked triple data even when the whole unordered Hamiltonian fiber is much larger. It is not a count of isomorphism classes; to get those one must also quotient by relabeling.

## 4. A five-vertex joint obstruction, and what cannot be imported from Collatz

Require exactly two cyclic triples \(0ab\) and \(0cd\), intersecting only at \(0\), and no others. Every four-vertex restriction is realizable; the two positive constraints are consistent, and no individual forbidden triple is forced by the positive constraints alone. Nevertheless the whole pattern is impossible. After naming the vertices so \(0\to a\to b\to0\) and \(0\to c\to d\to0\), forbidding \(0ad\) and \(0cb\) forces \(d\to a\) and \(b\to c\). The outside four vertices now span the cycle \(a\to b\to c\to d\to a\). Its two diagonals independently make exactly two cyclic triples on those four vertices. Thus allowing just one additional cyclic triple outside 0 also cannot work. This is a small instance of the same proof strategy as [the Moore route-color tests](../moore_route_colors_tournament_ports_20261126/README.md) and [the joint Pollock moment tests](../pollock_joint_local_moments_and_selector_20261126/README.md): a shared contact must satisfy all the constraints at once, even when the projected tests pass. The actual graph or Diophantine constraints in those problems are different.

Exhausting labelled quintuple patterns gives 253 passing all four-vertex restrictions and 208 globally realizable. The 45 false positives consist of exactly two relabelling types: the 15 two-triple patterns just described and their 30 variants with one additional outside triple. So this joint obstruction, together with the four-vertex checks, suffices at order five. No completeness claim is made at higher order.

By contrast, attempting to identify tournament **count descent** with the accelerated odd map \(x\mapsto(3x+1)/2^{v_2(3x+1)}\) fails in a provable way. Realizable counts 9 and 37 both go directly to the impossible count 7; \(65\to49\to37\to7\) gives another instance. The other gap 21 goes directly to 1 and has **no** incoming odd step at all, since \(3x+1\) is never divisible by 3 and stripping powers of 2 cannot change that. Thus the two gaps have genuinely different dynamics despite sharing an obstruction by strong-component factorization. The useful analogy is the parity code of marked contacts and its compatibility constraints, not closure of Hamiltonian heights under ordinary Collatz.

## Reproduction

Run the following from this directory (the output order is pinned by Gray-code iteration):
```sh
g++ -O3 -std=c++20 census.cpp -o census
./census > census.replay
diff -u census.out census.replay
python3 analyze.py
python3 obstruction.py
```
