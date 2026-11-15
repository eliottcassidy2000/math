# Stronger dyadic local limits using both complete digit tapes

Use the physical mask-zero carrier
```
F(x,y)=100xy+89x+11[p(x)+p(y)+1+(x mod 2)]
p(0)=0,   p(2u+1)=p(u),   p(2u)=p(u)+u+(u mod 2).
```
The [earlier dyadic proof](../chanpak_dyadic_collision_ports_20261115/DYADIC.md) conditions on low input residues so that the bilinear term is constant. This conditioning is not needed. We reach pointwise relative accuracy up to modulus size T^(1/6−ε) and relative RMS accuracy up to T^(1/4−ε) for mesoscopic windows on the true sharp shell. These are local laws for presentations, not estimates of exact collision energy.

## 1. Orthogonality with the physical bilinear term present

Put q(z)=p(z)/2, G=(F−1)/2. Then
```
G=50xy+89 floor(x/2)+50(x mod 2)+11[q(x)+q(y)]+5.
```
In q, the interaction between a lower bit and an upper bit has coefficient 1 at distance one, and 2^(distance−2) at greater distance, as in the earlier digit proof. For an odd character of denominator 2^d, the interactions vanish beyond distance w=d+1, and at w the coefficient has order two. In each consecutive free-digit range of length L pair bits separated by w, using the lower and upper halves of successive 2w blocks. The number paired is
```
k(L,w)=w floor(L/(2w))+max(0,(L mod 2w)−w) ≥(L−w)/2.
```
Label all lower paired positions of both tapes A, all upper B, and fix the unpaired digits.

Every B position is at least w in its tape, so 50xy modulo 2^d is independent of **all** B digits even if the rest of both tapes varies. It enters only a weight on A. The linear and parity terms separate as well. The phase matrix indexed by the two joint assignment groups A,B thus has orthogonal rows: choose a tape in which two A assignments differ, take its highest differing A bit i, and sum on the B bit i+w in the row inner product. Lower differences vanish there, higher positions agree, and the matching difference contributes an odd multiple of 2^(d−1). The two B choices cancel exactly. For k pairs its spectral norm is 2^(k/2), yielding a normalized average at most 2^(−k/2) by Cauchy–Schwarz. Unfixing unused digits preserves the bound.

Consequently in any aligned rectangle with side lengths 2^r,2^s, uniformly in a reduced character h/2^d with d≥1,
```
|mean exp(2πi h G/2^d)| ≪ min(1,2^(-(r+s−2d)/4)).
```
The bound remains valid if x is restricted to odd values, losing at most one free digit. It is uniform even when one tape is too short to furnish a pair. For a single tape of length 2^r with the other address fixed, the bound is instead O(min(1,2^(-(r−d)/4))).

## 2. Sharp-cut local limits

In an aligned rectangle let M be its presentation count and A_b its count with G≡b mod 2^m. Fourier inversion by reduced conductor, using 2^(d−1) characters of denominator 2^d, gives
```
sup_b |2^m A_b/M−1| ≪ 2^((6m−r−s)/4).
```
Parseval with the same grouping gives the relative RMS error
```
(2^(-m) Σ_b |2^m A_b/M−1|²)^(1/2) ≪ 2^(m−(r+s)/4).
```
For one aligned row of length 2^r the exponents are (5m−r)/4 and (3m−r)/4, respectively. Here b corresponds to F≡2b+1 mod 2^(m+1).

Fix 0<δ<1/2, D=ceil(T^δ) and ε>0. In the domain x,y≥D and U<F≤V with T/2≤U<V≤T, there exists κ>0 depending on δ,ε so that, for V−U≥T^(1−κ), the pointwise relative error is O(T^(-κ)) for 2^m≤T^(1/6−ε), and the relative RMS error is O(T^(-κ)) for 2^m≤T^(1/4−ε). Normalize by the ordinary presentation mass of the chosen domain, or by its odd-x subfamily mass if that restriction is used. The RMS error implies that all but a vanishing fraction of classes have the expected order of mass in each such window, not necessarily that the exceptional classes are common to all windows.

For details, the bound 0≤p(z)≤z gives F=100xy+O(x+y+1). In every dyadic box meeting the shell, xy is comparable to T and the error is O(T/D), even on the full box. Divide each side into 2^j aligned pieces with j=floor(η log₂T), where 0<η<min(δ,ε) is fixed and sufficiently small. Each cell's free lengths sum to log₂T−2j+O(1), so the whole-cell errors in the claimed ranges are negative powers of T. Approximate from inside and outside by complete cells. A cell straddling the domain must lie in a product strip of width O(T(2^(-j)+D^(-1))) at one count boundary, or an address strip of width O(D2^(-j)+1) at an address cutoff. Summing row lengths plus rounding gives ordinary mass difference O(T(2^(-j)+D^(-1)) log T). The total mass is at least a constant times (V−U) log T: sum the inward-perturbed hyperbola intervals over D≤x≤T/(400D); the sum of 1/x remains comparable to log T on odd x as well. Choose κ smaller to make the boundary mass a relative negative power.

Adding centered residue vectors of disjoint complete cells is valid in sup or normalized RMS norm by the triangle inequality. For the actual domain its vector lies coordinatewise between the complete-cell vectors from inside and outside. Its error from uniformity is therefore bounded by both envelope errors plus twice the ordinary boundary mass, with the same normalization. This is why the RMS law survives the exact cutoff; replacing the cutoff smoothly is unnecessary.

## 3. An exact collision outside the known valuation-transfer orbit

In the [dyadic transport quotient](../chanpak_dyadic_collision_ports_20261115/MASK_ZERO_ORBITS.md), put Y=y+1 and use odd cores a,b and total valuation j. Then
```
F(x,Y−1)=100·2^j ab+11+22[g(b)−q₀(a)]
g(b)=p(b−1)/2,     q₀(a)=(a−1−p(a))/2.
```
The cores (3,15) and (5,9) are distinct but have the same odd product 45 and checksum difference 3. Both yield exactly 4500·2^j+77 for all j≥0, even after the single-orbit quotient. Split j approximately in half between x and Y to realize these collisions with both addresses comparable to the square root of F. This single ray still contributes only O(log T) different outputs up to T.

Indeed any two fixed core rays collide infinitely often if and only if both their products and checksums agree. Otherwise a nonzero checksum difference fixes a nonzero right-hand side in the equality of the two scaled products, bounding the smaller exponent by its valuation; each bounded exponent fixes at most one opposite exponent. A transpose alone cannot supply a distinct equal-product, equal-checksum core: for a=2u+1, the sum g(a)+q₀(a)=(3u+(u mod 2))/2 strictly increases. This exact nontransposed collision is a control for stronger arguments which might mistakenly assume that the transport orbit is the only ambiguity.

Run `python verify_dyadic_local.py` for finite checks of the paired phase bound and the exact alias. The missing work is still the high quotient simultaneous with the dyadic residue, not another local factor calculation.
