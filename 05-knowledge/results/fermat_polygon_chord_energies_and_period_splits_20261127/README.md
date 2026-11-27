# Squares, nonsquares, and how to draw the Fermat-prime polygon

This is a geometric reading of Gauss's classical period method, not a new constructibility theorem. It builds on the [finite-field half-cycles](../abc_pythagorean_fermat_half_cycles_20261127/README.md) and the [Fermat-neighbor experiments](../clocks_carriers_genus_20261102/FERMAT_NEIGHBOURS.md). Gauss's [Disquisitiones, articles 353–365](https://la.wikisource.org/wiki/Pagina:Gauss,_Carl_Friedrich_-_Werke_%281870%29.djvu/484) give the original period-to-quadratic construction.

## 1. The two half-cycles are two unequal chord-energy classes

Place the regular \(p\)-gon on the complex unit circle with vertex \(\zeta^j\), \(\zeta=e^{2\pi i/p}\). Use the **exponent** labels \(j\in\mathbb F_p\). For \(p\equiv1\bmod4\) both the nonzero squares \(S\) and nonsquares \(N\) are closed under \(j\mapsto-j\), so their vector sums
\[
s=\sum_{j\in S}\zeta^j,\qquad n=\sum_{j\in N}\zeta^j
\]
are horizontal, i.e. real. Under the previous odd-pair readout, for Fermat prime \(p\ge17\) these sets are exactly the two cycles of multiplication by \(\lambda=9/4\) on the exponent labels. This multiplication is **not** an ordinary Euclidean rotation of the p-gon: an ordinary rotation adds a fixed exponent. The map is a symmetry of the labelled modular data and provides an orbit along which to sum Euclidean vectors.

Writing \(E_S,E_N\) for the sums of squared chord lengths from vertex 1 to each class of vertices, the Gauss periods give
\[
s=\tfrac12(-1+\sqrt p),\quad n=\tfrac12(-1-\sqrt p),\qquad
E_S=p-\sqrt p,\quad E_N=p+\sqrt p.
\]
In fact each squared chord contributes \(2-2\operatorname{Re}\zeta^j\), so the first radical is literally **half the difference of two chord energies**, not the square root of an individual polygon side. This is why the additive phase on a half-cycle does useful geometric work while a generic multiplicative character on the same cycle sums to zero.

More generally for two symmetric label sets \(U,V\), the total cross-chord energy is
\[
\sum_{u\in U,v\in V}|\zeta^u-\zeta^v|^2
=2|U||V|-2\eta_U\eta_V,\qquad
\eta_U=\sum_{u\in U}\zeta^u .
\]
Thus the usual sum and product of two daughter periods amount geometrically to a total horizontal vector and a cross-chord energy.

## 2. Why the geometry can recurse by quadratic constructions

If \(g\) generates \(\mathbb F_p^\times\), let \(H_j=\langle g^{2^j}\rangle\). A parent period on coset \(aH_j\) splits into two daughter periods on \(aH_{j+1}\) and \(ag^{2^j}H_{j+1}\). Their sum is the parent. Their product expands as a sum over \(x+y\) from the two daughter sets. Multiplying both inputs by an element of \(H_j\) preserves the count of such pairs (possibly swapping their roles), so the coefficient at \(x+y\) is constant on each **parent-level** coset, with a separate constant for zero. Their product is therefore an integer linear combination of the already constructed parent periods. Counting in the residue field computes the coefficients; the actual polygon is not needed as an input.

If the daughter periods have sum \(u\) and product \(v\), their horizontal positions are \((u\pm\sqrt{u^2-4v})/2\). Signed lengths are constructed by parallels and similar triangles; the square root of a known nonnegative number is a mean proportional obtained from a semicircle. Thus the geometric split is a midpoint plus or minus half a constructible difference. For \(p-1=2^r\), stop at \(\{\pm a\}\), giving \(2\cos(2\pi a/p)\). Each construction up to this point lives on the real projection. The unit-circle intersection finishes the polygon.

For an odd prime \(p\), the degree of this cosine over \(\mathbb Q\) is \((p-1)/2\), since the cyclotomic group identifies just \(a\) and \(-a\) in the real projection. Constructibility forces that degree to be a power of two. Conversely the dyadic periods construct it whenever the degree is a power of two. If \(2^r+1\) is prime then \(r\) is a power of two, since an odd factor of \(r\) would factor the sum of powers. The inability to keep halving the group at an odd factor is the real obstruction, not the mere presence of square and nonsquare labels (which exist at every odd prime).

**The projection and the full vertex have different costs.** The cosine needs \(r-1\) layers of quadratic extension. But \(\sin(2\pi/p)=\cos(2\pi(p-4)/(4p))\) is a primitive real \(4p\)-cyclotomic coordinate since \((p-4,4p)=1\). Its degree is \(p-1=2^r\), requiring a further layer. The final circle intersection has a precise algebraic role. At 17 the projection has degree 8 and the full oriented coordinate has degree 16.

## 3. A short geometric recipe for 17, using the same \(9/4\) labels

At \(p=17\) use \(g=3/2=10\), with \(g^2=15\). Set \(P_{j,k}=\sum_{a\in g^k H_j}\zeta^a\) for \(0\le k<2^j\). Direct modular pair counts give:

| Daughters multiplied | Product |
|---|---|
| \(P_{1,0}P_{1,1}\) | \(-4\) |
| \(P_{2,0}P_{2,2}\) | \(-1\) |
| \(P_{2,1}P_{2,3}\) | \(-1\) |
| \(P_{3,0}P_{3,4}\) | \(P_{2,3}\) |

For example \(P_{2,0}\) has labels \(\{\pm1,\pm4\}\); at the last split, adding the two signed pairs \(\{\pm1\}+\{\pm4\}\) gives \(\{\pm3,\pm5\}\), exactly the labels of \(P_{2,3}\). At level one the two parent classes have chord energies \(17-\sqrt{17}\) and \(17+\sqrt{17}\); the complete cross-chord energy of those two eight-vertex sets is 136. The level-two sibling products \(-1\) correspond to cross energies 34.

To draw the next vertex on the unit circle starting at \((1,0)\):

1. Construct \(d_0=\sqrt{17}\), \(s=(d_0-1)/2\), \(n=(-d_0-1)/2\).
2. Construct \(d_1=\sqrt{s^2+4}\). Form \(A=(s+d_1)/2\) and \(B=(n+2d_0/d_1)/2\).
3. Construct \(d_2=\sqrt{A^2-4B}\) and draw the vertical at \(x=(A+d_2)/4\). Its upper circle intersection is the next vertex. Carry the resulting chord around the circle.

The second discriminant \(\sqrt{n^2+4}\) needs no independent square-root construction, since \((s^2+4)(n^2+4)=68\) and both radicals are positive. The sign choices are exact: \(A=P_{2,0}>0\) because its two cosine angles lie inside \((0,\pi/2)\), and \(B=P_{2,3}=2\cos(6\pi/17)+2\cos(10\pi/17)>0\) by the cosine addition identity. Finally \(2\cos(2\pi/17)\) is the larger daughter of \(A\), because the other is \(2\cos(8\pi/17)\). The construction uses three square-root extractions for the real coordinate and one circle intersection, attaining the degree lower bound. The unit-circle side length is approximately \(0.3674990356\); rounding is only a check, never part of the construction.

## 4. A matrix check and the boundary of the analogy

Color the unordered chords by the square class of the label difference; at 17 this is unambiguous because \(-1\) is a square. The adjacency \(A_S\) formed by square chords satisfies
\[
A_S^2+A_S-4I=4J .
\]
Each Fourier wave of frequency \(k\ne0\) has eigenvalue \(s\) or \(n\) according as \(k\) is square or nonsquare. The constant frequency has eigenvalue 8. Substituting the period equation \(x^2+x-4=0\) gives the matrix identity. This is structurally close to the [Moore Gram equation](../odd_pair_matrices_moore_gram_and_pollock_fibers_20261126/README.md), but a pair here has three or four common neighbors, not the Moore unique-route law. The first matrix split detects two eight-dimensional classes; the deeper period refinements are what actually locate an individual vertex projection.

Run `node replay.js` to reproduce the integer period-product rows, the 289 entries of the matrix identity, and numerical checks of the constructed chord.
