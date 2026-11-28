# Finite zeta nonvanishing, reciprocal tournament readouts, and Moore pentagon matchings

2026-11-28. Three distinct results: the announced arithmetic nonvanishing transfers exactly to two **marked** readings of a transitive tournament; ordinary short-cycle totals miss a Hamiltonian-path change even on a single arc; and length five yields 28,467 new simultaneous matching boards for the current degree-57 construction.

## 1. Arithmetic: two independent readings of one transitive class

[Matsusaka's Nov. 27 announcement](https://sites.google.com/site/tmatsusaka556/) and [linked manuscript](https://drive.google.com/file/d/1bRca-tEwQW5aieYNZTHlYkLCfMHmjEBb/view?usp=share_link) state that \(Z_{\mathcal A}(3)\) and \(Z_{\mathcal A}(5)\) are not both zero. Here we use the stated conclusion, not a verification of its proof. In [Kaneko–Zagier's notation](https://www2.math.kyushu-u.ac.jp/~mkaneko/papers/FMZV.pdf), equation (11), \(Z_{\mathcal A}(k)\) has components \(B_{p-k}/k\bmod p\), for sufficiently large primes. In \(\mathcal A=(\prod_p\mathbf F_p)/(\bigoplus_p\mathbf F_p)\), nonzero means that infinitely many prime components are nonzero. The pair conclusion does not specify which component works, or require both to work at the same prime. These values are not the ordinary Riemann zeta values and not the depth-one sums \(\sum_{x\ne0}x^{-k}\), which vanish for sufficiently large \(p\).

On \(\mathbf F_p^\times\), put the field labels in the order \(c,2c,\ldots,(p-1)c\) and orient forward. For \(k=3,5\) give an arc \(x\to y\) the marked reading \(x^{-1}y^{1-k}\), and let \(W_k(c)\) be the sum. For all primes \(p>5\),
\[
W_k(c)=c^{-k} B_{p-k}\pmod p.
\]
Indeed, at \(c=1\) sum \(i^{p-2}\) over \(i<j\) by the Bernoulli polynomial formula and multiply by \(j^{p-k}\). In the resulting sum over \(j\ne0\) the only surviving Bernoulli index is \(m=p-k\), since the exponent \(p-1-m+p-k\) must be divisible by \(p-1\). Its coefficient \(\binom{p-1}{m}/(p-1)\) is \(-1\) since \(m\) is even, and summing the power of \(j\) supplies the second minus sign. Dilation contributes \(c^{-k}\). Reversing the order at fixed labels negates each reading: the sum against the two complementary orientations is \((\sum_x x^{-1})(\sum_y y^{1-k})-\sum_x x^{-k}=0\).

One need not know which weight survives. Use the **same** mixed arc reading \(x^{-1}(y^{-2}+y^{-4})\) at \(c=1\) and \(c=2\). The measurement matrix on \((3Z_{\mathcal A}(3)(p),5Z_{\mathcal A}(5)(p))\) is
\[
\begin{pmatrix}1&1\\2^{-3}&2^{-5}\end{pmatrix},\qquad \det=-3/32.
\]
Thus the two finite measurements recover both weights at every \(p>5\), and on infinitely many primes at least one measurement is nonzero by Matsusaka's stated result. Every tournament in the construction is transitive, with Hamiltonian-path count 1. This is an exact example of arithmetic information erased by passing to the unmarked tournament class; it is not a statement about new integer heights. More generally for distinct fixed weights \(k_j\), dilation by \(2^i\), \(0\le i<r\), gives the Vandermonde matrix on the nodes \(2^{-k_j}\). Only finitely many primes divide the resulting determinant. The extension is an inverse for marked weight profiles, not a claim that new weights are arithmetically nonzero.

A related caution for the [Cullen–Woodall split-norm work](https://github.com/eliottcassidy2000/math/commit/f018aa79e5a2c13de34fe8c4750e5040607a61f2): the product of two nonzero elements of \(\mathcal A\) can vanish on disjoint supports. A coordinatewise noncancelling scalar detector does exist: choose a nonsquare \(\nu_p\) in each odd \(\mathbf F_p\). Then \(a_p^2-\nu_p b_p^2\) is zero exactly when both coordinates are zero. This is a norm from the nonsplit quadratic extension at each prime, not a fixed rational combination or an assertion that its class is itself an FMZV.

## 2. Chan–Pak: a single arc that 3, 4, and 5 cannot see

Take the order \(0,\ldots,5\) and reverse exactly the pairs (0,1), (0,5), (1,4). Reverse the additional arc at (1,2) for a second tournament. Both have directed cycle counts
\[
(C_3,C_4,C_5)=(6,7,6).
\]
The first has two pairs of disjoint directed triangles and 33 Hamiltonian paths; the second has one pair and 29 paths. The replay checks the path counts directly over all 720 orders of the six sites.

This difference is also explained by the [Grinberg–Stanley Rédei–Berge expansion](https://arxiv.org/html/2307.05569): each collection of vertex-disjoint nontrivial odd cycles contributes \(2^{\#\text{cycles}}\). On six vertices \(H=1+2(C_3+C_5)+4D_{3,3}\). Therefore six is the smallest order for which the pair \((C_3,C_5)\) can hide a path count at fixed order. The first missing datum is **disjointness**, not a higher odd cycle. The example neither explains the gaps 7 and 21 nor assumes independence of disjoint-cycle collections.

## 3. Moore: length five as a local matching constraint

For a simple \(d\)-regular, square-free graph (\(d\ge2\)) on \(d^2+1\) vertices, either a triangle gives a nonzero third trace, or its radius-two neighborhoods attain the Moore bound. In the latter case the graph is Moore and every edge is on \((d-1)^2\) pentagons. For an edge \(uv\), choose the neighbor \(x\ne v\) of \(u\) and the neighbor \(y\ne u\) of \(v\). Girth at least five gives distinct nonadjacent \(x,y\); their unique common neighbor closes the pentagon. Thus \(C_5=(d^2+1)d(d-1)^2/10\). This is a graph-theoretic 3-or-5 statement independent of the arithmetic theorem.

Use this more strongly on a partial Moore chart. Suppose two saturated sites \(u,v\) are nonadjacent, with their unique common neighbor \(w\). Their neighborhoods cannot change. In any completion, the edges between
\[
L=N(u)\setminus\{w\},\qquad R=N(v)\setminus\{w\}
\]
form a **perfect matching**, not simply a sufficient number of edges. For \(x\in L\), \(x\) is not adjacent to \(v\), and its unique two-step route to \(v\) cannot use \(w\) (that would make a triangle through \(u\)); its middle site is in \(R\). Two sites of \(L\) cannot use the same middle site, since they already share \(u\). Each match is one pentagon through the literal wedge \(u-w-v\). A board built from known edges and physically legal optional edges must therefore pass Hall; alternating cycles give exact per-board supports. The same optional edge can occur in multiple boards, so independent matchings are not a joint witness.

For a compact algebraic view, let \(A_0\) be the partial adjacency, \(E\) be the rectangular selector of all saturated vertices, \(L_0=A_0E\), and \(X\) the unknown edges of a completion. The cubic identity \(A^3=dA-(d-1)I+(d-1)J\) gives the shared linear system
\[
L_0^\top X L_0=E^\top\{dA_0-(d-1)I+(d-1)J\}E-L_0^\top A_0 L_0.
\]
A matrix entry gives the quota; the matching condition allocates it to distinct slots. This links back to the repo's [Wilson route defect](../moore57_quadratic_tensors_20261125/WILSON-ACTION-AND-RIGID-MOORE-CHARTS.md), where length three and four only penalize collisions. Here length five imposes simultaneous positive placement obligations.

The [15-star witness](../moore57_fifteenth_complete_star_20261128/README.md) has **60,756 pentagons**; any completion would have **58,094,400**. The saturated wedges give **28,467** boards, with 154,070 already realized matching incidences and 1,440,082 still required. All 28,467 pass the perfect-matching and alternating-support checks; the checks delete or force no candidate. The least number of optional edges on an unfinished board is 1,286. These are counts of board incidences, not counts of independently placeable edges. No new neighborhood or completion is claimed. The clean next use is to enforce shared edge variables across several boards, rather than repeating their individual projections.

## Replay

~~~sh
g++ -O3 -std=c++17 wedge.cpp -o /tmp/moore-pentagon
for n in 5 10 50; do python3 controls.py "$n" | /tmp/moore-pentagon; done
python3 ../moore57_fifteenth_complete_star_20261128/replay.py --adjacency | /tmp/moore-pentagon
python3 finite_readouts.py
~~~

The three complete controls give 1, 12 and 1,260 pentagons. The arithmetic replay tests the identities at 50 primes; it does not prove nonvanishing infinitely often.
