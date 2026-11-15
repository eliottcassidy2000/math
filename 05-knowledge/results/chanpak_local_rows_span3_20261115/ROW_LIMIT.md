# A mesoscopic local limit for a physical one-tape row

Let p be the carrier checksum
```
p(0)=0,   p(2u+1)=p(u),   p(2u)=p(u)+u+(u mod 2).
```
Consider G(n)=c n+k p(n)+φ(n mod 2)+b₀, with fixed integers c,k>0, k odd, b₀ odd, c+φ even. These conditions include the actual physical mask-zero row F(n,1)=189n+11p(n)+11(n mod 2)+11 and the 3-cycle row 18n+3p(n)+2(n mod 2)+9 of [HYP-9752](../../hypotheses/HYP-9752-row-collision-lemma-a-carrier-row-cu-plus-kp-u-has-linear-energy-which-would-give-the-chan-pak-spectrum-positive-lower-density.md). Restrict n≥1 for physical values.

**Result.** For N=2^L, the row has the following explicit density at macroscopic scale:
```
f(t) = (1/k) log((c+k)/max(c,t))    if 0<t<c+k,
       0                          otherwise.
```
The interval counting error, without residue conditions, is O(√N) uniformly in the interval. Moreover:
- Let 0<η<1/5 and α>(1+5η)/2, with α<1. For every m≤ηL the physical counts in each odd output residue modulo 2^(m+1), and every interval of length at least N^α lying in a fixed compact subinterval of the interior of the scaled support, obey the same local limit with a relative power-saving error.
- In relative RMS over the residue classes, the same result holds for 0<η<1/3 and α>(1+3η)/2.

For an interval I here the predicted count in one residue is (2^{-m}Nint_{I/N}f(t),dt). This is a mesoscopic limit theorem for **presentations**, not a singleton-resolution result. It does not prove the row collision lemma.

## 1. Exact digit geometry and the ordinary limit

Write the L binary digits of n in little-endian order as b_i, with b_L=0, and let R(n)=rev_L(n). Iterating the recursion gives
```
p(n) = sum_{i=0}^{L−1}(1−b_i)[floor(n/2^(i+1))+b_(i+1)]
     = n [1−R(n)/N] + O(L+1).
```
For the second expression, the unfloored digit sum has coefficient 1−1/N−R(n)/N, while each remaining term has absolute value at most one. Thus G(n)/N differs by O((L+1)/N) from the two-coordinate surface
```
ψ(x,y)=x[c+k(1−y)],     x=n/N, y=R(n)/N.
```
Under uniform x,y on the unit square, c+k(1−y) is uniform on [c,c+k]. Conditional on its value q, ψ is uniform on [0,q]; integrate (1/k)·(1/q), over q≥t, to obtain f.

Fixing the first s digits of x and y fixes respectively the highest and lowest s digits of n. For 2s≤L, every square in the resulting grid has exactly 2^(L−2s) addresses, equal to N times its area. This is the digit-reversal net property, with no probabilistic assumption on n. Take s=floor(L/2). In any cell the variation of ψ is O(2^(-s)) by bounded derivatives. Comparing wholly contained cells with all cells meeting the interval cut leaves only cells in two strips whose ψ-values are within O(2^(-s)+(L+1)/N) of the interval's ends. Their total area is of that order because f is bounded. This gives the uniform O(√N) counting error (including the exclusion of n=0). No Lipschitz condition on p itself is needed.

## 2. Keep free middle digits for the residue refinement

Fix instead s=floor(σL), leaving R=L−2s consecutive free digits inside each grid cell. The symbol R here is a length, not the reversal in §1. We use the triangular orthogonality in the [dyadic-rank proof](../chanpak_dyadic_collision_ports_20261115/DYADIC.md). For a character of denominator 2^d in q(n)=p(n)/2, its mixed digit interaction vanishes beyond distance d+1, and is 2^(d−1) times an odd unit at precisely that distance. Alternating digit blocks separated into two assignment groups by pairing at distance d+1 gives at least (R−d−1)/2 pairs and exactly orthogonal rows. The last differing lower assignment hits its paired upper digit with coefficient 2^(d−1), while all earlier differences vanish at that upper digit. This proof admits any linear twist in the input, absorbed into the two side weights.

Each cell fixes input parity (take s≥1). Thus (G(n)−1)/2 is a constant plus k q(n) plus a linear form in the free digits. An output character of reduced denominator 2^d is of exactly the preceding kind, because k is odd. There is no need to fix the whole input residue to extract this estimate. The normalized character bound is
```
|mean χ| ≤ min(1, 2^(-(R−d−1)/4)).
```
Fourier inversion over the 2^m odd output classes bounds relative sup error on a whole cell by O(2^((5m+2s−L+1)/4)). Orthogonality (Parseval) bounds relative RMS error by O(2^((3m+2s−L+1)/4)). These estimates are independent of which high and low digits are fixed.

Choose σ strictly between 1−α and (1−5η)/2 for the sup result, replacing 5 by 3 for RMS. On the specified compact part of the support the interval's predicted mass is comparable to |I|/N. The strip error O(2^(-s)+(L+1)/N) is smaller by a fixed power of N. Add full-cell estimates for the inner approximation, bound the boundary contributions by the whole straddling cells, and compare with the square area formula. For RMS use the triangle inequality in the normalized ℓ² norm on residue classes; the boundary vector is bounded componentwise by a sum of whole-cell vectors of total mass O(2^(-s)+(L+1)/N). This gives the claimed relative power errors. The possible n=0 correction is negligible even relative to the residue main term because α>η.

The remaining obstacle to density one or the local limit at integer resolution is exactly what the geometry does not retain: the high quotient simultaneous with the prescribed lower residue bits. The phase count can be uniform in each short *macroscopic* window while all exact collisions within that window remain uncontrolled.

## 3. A deterministic approximate-target corollary

Fix a compact J inside (0,c+k). Given an absolute scaled target tN with t in J and any prescribed odd residue modulo 2^(m+1), choose a grid cell containing a point of ψ=t. For instance take y=0 and x=t/(c+k); the endpoints are interpreted as closed limits of the chosen cell. Take s=floor((L−5m−C)/2) with a sufficiently large absolute C and s≥1. The sup estimate in §2 is then <1, so **every** odd residue occurs among the addresses in that cell. The cell and the digit approximation give a physical count at distance
```
O(√N · 2^(5m/2) + L)
```
from the prescribed target, in its exact residue. The size of the free block is 5m+O(1); enumeration of its middle digits is a deterministic decoder using O(2^(5m)·L) elementary recurrence steps. This trades arithmetic precision for quantitative absolute error, without claiming that the error falls below the modulus.
