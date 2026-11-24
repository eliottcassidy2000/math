# How much of the staircase lies on doubled nodes?

**PROVED.** Fix the directed path on the ordered vertices \(1,\ldots,n\), \(n\ge 2\). Write \(Q=\binom n2\), \(h=\lfloor n/2\rfloor\), \(f=n-2h\), \(o=(Q+h)/2\), and \(I_n=n!/(2^h h!)\), with \(I_0=1\). Among the \(2^{Q-n+1}\) boards let \(M_n\) count those whose ordinary tournament is self-converse. The complement counts exactly the boards on merged nodes of even mass. We obtain an exact union formula for \(M_n\), a positive finite sum for its marked version, and
\[
\boxed{M_n=I_n2^{o-n+1}\left(1+\frac{1+f}{n}+O(n^{-2})\right).}
\]
In particular the fraction of boards on doubled nodes is
\[
1-I_n2^{-(Q-h)/2}\left(1+\frac{1+f}{n}+O(n^{-2})\right).
\]
The unmerged fraction has base-two logarithm \(-n^2/4+(n/2)\log_2n+O(n)\). These are board-weighted probabilities. The mass of an isomorphism class weights its presentations by \(H/|\mathrm{Aut}|\). The rooted-path distribution adds a first-order asymmetry depending on the parity of \(n\).

The board convention and reflection come from [THM-409, self-converse perspective flip](../../../01-canon/theorems/THM-409-self-converse-perspective-flip-involution.md) and [THM-643, blue and black line parity](../../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md). The [earlier marked mass](../session90_even_graphs_odd_descent_20261029/carriers/MARKED_MASS.md) retains a symmetry mark on an unpinned object. Eplett describes the structure of self-converse tournaments. McKay gives the unweighted class census through order 13. The cancellation below comes from two orientations of an interval in the pinned path. [Eplett, *Self-Converse Tournaments*](https://www.cambridge.org/core/journals/canadian-mathematical-bulletin/article/selfconverse-tournaments/ED7BC2A1F6FF2D7F31547D771CE69B80) · [McKay, combinatorial data](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html)

## Exact union of affine families

A nonidentity involutive automorphism would exchange the endpoints of an arc and reverse that arc. Thus the automorphism group has odd order. For a self-converse tournament the union of the automorphisms and all reversing permutations is a group of twice that order. It has an involution, necessarily reversing. Such a permutation fixes at most one vertex, since the arc between any two fixed vertices would be reversed. There are \(I_n\) candidates with \(h\) transpositions and \(f\) fixed vertices.

For each candidate \(\pi\) let \(E_\pi\) be its event on pinned boards. For \(i<j\), let \(x_{ij}=0\) mean \(i\to j\). The pins have \(x_{i,i+1}=0\). For every \(i<j\) impose
\[
x_{\{\pi i,\pi j\}}=x_{ij}\mathbin\oplus [\pi i<\pi j],
\]
with the first subscript sorted. For a nonempty subset \(A\) of candidate involutions form a graph on the \(Q\) pair variables and one zero anchor. Join the anchor to every pin by an edge marked zero. Add an edge for each displayed equation for each \(\pi\in A\), marked by its XOR offset. Let \(c(A)\) be the component count including the anchor and \(b(A)=1\) when every cycle has even total offset, zero otherwise. Then
\[
 M_n=\sum_{A\ne\varnothing}(-1)^{|A|+1}b(A)2^{c(A)-1}.
\]
A consistent component has one free bit except the anchored component, whose value is forced. Inclusion-exclusion gives the formula.

For one involution the \(Q\) variables have \(o\) orbits, since its transposed vertex pairs supply the \(h\) fixed variables. Call a two-element orbit of **distinct** pinned edges parallel if their endpoints are exchanged in increasing order, reversed otherwise. Denote their counts by \(a(\pi),d(\pi)\). A fixed pinned variable is necessarily reversed in vertex order and always consistent. A parallel orbit forces contradictory pins. For \(a(\pi)=0\), orienting the pins costs \(n-1-d(\pi)\) independent bits. Hence
\[
 |E_\pi|=2^{o-n+1}\,\mathbf 1_{a(\pi)=0}\,2^{d(\pi)}.
\]
The largest possible \(d\) is \(\lfloor(n-1)/2\rfloor\), achieved by staircase reflection \(i\mapsto n+1-i\). This gives exactly \(B_n=2^{o-n+1+\lfloor(n-1)/2\rfloor}\) blue boards. Other reversing involutions can put black boards on the same kind of self-converse node. In fact \(\log_2(B_n/M_n)=-(n/2)\log_2 n+O(n)\).

## Positive formula by cancellation on path intervals

Write \(W_n=\sum_\pi|E_\pi|\), counting a board once for each of its reversing involutions. For fixed \(\pi\), expand
\[
 \mathbf 1_{a(\pi)=0}2^{d(\pi)}=(1-1)^{a(\pi)}(1+1)^{d(\pi)}.
\]
A term marks any subset of the two-edge pin orbits. Its sign is \((-1)^{\text{number of marked parallel orbits}}\). On the base path, the union of the marked edges splits into maximal runs. A run of length \(L\) occupies \(L+1\) consecutive vertices. The permutation of the marked subgraph takes a run onto one of equal length. On two mutually exchanged runs the involution is either increasing, giving sign \((-1)^L\), or decreasing, giving sign \(+1\). Both specify \(L+1\) disjoint transpositions and have the same number of completions. Summing these two orientations cancels every odd-length pair of runs. Each even-length pair contributes weight two.

A run mapped to itself must reverse. The identity would individually fix its edges, which cannot be marked. An odd-length reversal also individually fixes the middle edge. Thus a self-run has even length, uses the unique fixed midpoint, and contributes weight one. It is possible only for odd \(n\), at most once. These choices exhaust the terms, since different runs have disjoint vertex sets.

Explicitly, choose \(j\ge0\) paired runs of edge lengths \(2a_1,\ldots,2a_j\) with each \(a_i\ge1\). Each listed length specifies **two** runs. Optionally choose a self-run of length \(2d\) when \(f=1\), with \(d=0\) indicating its absence, and put \(t=[d>0]\). The number of intervals and their total number of vertices are
\[
 R=2j+t,\qquad V=4\sum_i a_i+2j+2d+t.
\]
Then the finite, positive marked sum is
\[
 \boxed{\frac{W_n}{2^{o-n+1}}=
 \sum_{j,(a_1,\ldots,a_j),d\, :\, V\le n}
 \frac{R!}{j!}\binom{n-V+R}{R} I_{n-V}.}
\]
Indeed \(R\) labeled intervals can be arranged in \(R!\) orders, and the \(n-V\) unused vertices fill \(R+1\) gaps of nonnegative size. Dividing by \(2^j j!\) forgets the orders of the two intervals in each pair and the labels on pairs, when summed over ordered lengths. The two orientations of every surviving pair cancel the \(2^j\) denominator. There are \(I_{n-V}\) completions on the remaining vertices.

## The first correction and the cost of forgetting the mark

Divide the last sum by \(I_n\). The empty choice contributes one. The choice of one pair of length-two runs contributes \((n-4)(n-5)I_{n-6}/I_n=1/n+O(n^{-2})\) for \(n\ge6\). At odd \(n\), the choice of one self-run of length two contributes \((n-2)I_{n-3}/I_n=1/n\) exactly.

Here is a uniform estimate for all the remaining terms. A prescribed collection of \(r\) transpositions in a uniform candidate has probability at most \((C/n)^r\) for an absolute \(C\). For even \(n=2h\) the denominator is a partial product of \(n-1,n-3,\ldots,1\), whose geometric mean is at least that of the whole list, at least \((h!)^{1/h}\ge h/e\). Odd \(n\) has larger odd factors. The same estimate with a larger \(C\) holds conditionally on a specified fixed vertex, whose probability is \(1/n\). For \(j\) pairs the placements of interval starts are at most \(n^{2j}/j!\) when summed over ordered pair lengths, allowing a factor 2 per pair for orientation. One pair of even length \(L\) specifies \(L+1\) transpositions. Summing its bounds over \(L=2,4,\ldots\) gives \(u_n=2n^2\sum_{L=2,4,\ldots}(C/n)^{L+1}=O(1/n)\), with the \(L\ge4\) part \(O(n^{-3})\). The bound for \(j\) pairs is \(u_n^j/j!\). The single optional self-run has at most \(n\) positions, requires a fixed vertex, and specifies \(d\) additional transpositions. Its sum contributes a factor at most \(\sum_{d\ge1}(C/n)^d=O(n^{-1})\), with \(d\ge2\) contributing \(O(n^{-2})\). All terms other than the three specified choices are thus \(O(n^{-2})\) in total.

It remains to pass from \(W_n\) to \(M_n\). Fix a compatible \(\pi\), and take a different \(\rho\) with nonempty intersection. The product \(g=\pi\rho\) is a nonidentity automorphism of any board in the intersection. Its order is odd and its support \(S\) has size \(s\ge3\). Since \(\pi g\pi=g^{-1}\), the set \(S\) is invariant under both reversers, which agree off \(S\). The \(q_s=s(2n-s-1)/2\) pairs incident to \(S\) have \(g\)-orbits of length at least three. A singleton would fix both endpoints or exchange them, and an automorphism cannot exchange endpoints of an arc. Under \(\pi\) alone, these variables have at least \(q_s/2\) orbits, with at most \(2s\) spine edges pinned. Thus at least \(q_s/2-2s\) free bits lie here. Under both involutions there are at most \(q_s/3\) bits on the same variables, since a bit on each \(g\)-orbit fixes the rest. All complementary variables are independent of these and subject to identical restrictions for the two counts. Consequently
\[
 |E_\pi\cap E_\rho|/|E_\pi|\le 2^{-q_s/6+2s}.
\]
There are at most \(n^s\) possible \(g\) of support size \(s\), each determining \(\rho\). Since \(q_s\ge s(n-1)/2\), the sum of the ratios over \(\rho\ne\pi\) is at most \(z^3/(1-z)\) whenever \(z=4n2^{-(n-1)/12}<1\). The union bound summed over marks gives \(0\le W_n-M_n\le W_n O(n^3 2^{-n/4})\). This proves the displayed asymptotic and shows that a uniformly sampled self-converse board has a unique reversing involution with probability tending to one.

## Finite census

**FINITE-EXACT.** Enumerating compatible affine families gives the following counts. The column headed “multiple” counts boards admitting at least two reversing involutions. The positive interval sum gives the marked column independently. The census enumerates board presentations with the specified Hamilton path, including every choice of off-path bits.

| \(n\) | All boards | Compatible involutions | Marked \(W_n\) | Self-converse \(M_n\) | Multiple | Blue \(B_n\) | On doubled nodes |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 | 0 | 1 | 0 |
| 3 | 2 | 3 | 4 | 2 | 1 | 2 | 0 |
| 4 | 8 | 2 | 6 | 6 | 0 | 4 | 2 |
| 5 | 64 | 12 | 76 | 52 | 9 | 16 | 12 |
| 6 | 1,024 | 10 | 272 | 240 | 16 | 64 | 784 |
| 7 | 32,768 | 78 | 8,320 | 7,308 | 408 | 512 | 25,460 |
| 8 | 2,097,152 | 68 | 59,904 | 58,712 | 596 | 4,096 | 2,038,440 |
| 9 | 268,435,456 | 672 | 4,648,960 | 4,507,860 | 63,741 | 65,536 | 263,927,596 |

Run `python3 -B mass.py` in this directory. The size-nine union occupies two 32 MiB bitmaps. The geometric argument gives unconditional frequencies. Selecting a board at a prescribed ordinary Hamilton count, the condition needed in the Chan–Pak question, remains open.
