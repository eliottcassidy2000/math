# Polynomial-size dyadic mixing for the actual physical carrier

For the order-five mask-zero physical carrier, write
```
F(x,y)=100xy+89x+11[p(x)+p(y)+1+(x mod 2)]
p(0)=0,   p(2u+1)=p(u),   p(2u)=p(u)+u+(u mod 2),
```
and set q(x)=p(x)/2. The [previous all-modulus mixing argument](../chanpak_voltage_profile_20261115/README.md) allows increasing moduli of logarithmic size. This note uses the exact digit interaction to give a power-size result for the dyadic part, inside the actual sharp-cut domain. It does not assert an unrestricted modulus bound or distinct-value density one.

**Main result.** Fix 0<δ<1/2 and 0<η<min(δ,1/8). Put D=ceil(T^δ). Uniformly in integers m with 0≤m≤(min(δ,1/8)−η)log₂ T and in odd classes b mod 2^(m+1),
```
#{x,y≥D : T/2<F(x,y)≤T, F(x,y)≡b mod 2^(m+1)}
   = 2^(−m) P_T [1+O(T^(−κ))],
```
with κ>0 depending only on δ,η and the form; P_T counts all pairs in that shell. The same relative estimate with a possibly smaller κ holds on every subinterval (U,V] of [T/2,T] of length at least T^(1−κ). It also holds after restricting to x odd, with the corresponding restricted presentation mass.

## 1. A triangular matrix on the digit interaction

If z has binary digits b_i counted from the least significant end, the defining recurrence gives
```
q(z)=sum_i (1−b_i)b_(i+1)
       + sum_(j≥i+2) (1−b_i)b_j 2^(j−i−2).
```
Indeed removing a trailing zero contributes floor(z/2)+(floor(z/2) mod 2) to p. The bit adjacent to that zero therefore occurs twice; after division by 2 its coefficient is 1. Trailing one contributes zero. This proves the formula, and integrality of q, by iteration.

Consider the phase exp(2πi h q(z)/2^d) with h odd and d≥1. Mixed coefficients beyond distance w=d+1 vanish modulo 2^d. At distance w, the coefficient is 2^(d−1), so switching both endpoints gives a sign. Suppose L consecutive digits are free, with every other digit fixed. Partition those digits into successive blocks of length 2w, possibly a partial last block. Within each block pair its i-th bit with the bit at distance w above it whenever both exist. Put the lower positions into A, the upper into B, and condition on the unused bits. There are
```
k = w floor(L/(2w)) + max(0, (L mod 2w)−w)
  >= (L−w)/2
```
pairs.

The rectangular phase matrix on the 2^k assignments to A and the 2^k assignments to B has orthogonal rows. Given two different A assignments, choose their **highest** differing position i. In the column bit paired with i, at i+w, the term from i changes the exponent by 2^(d−1) modulo 2^d. Every lower differing A position is at distance exceeding w and contributes zero; every higher A position agrees by choice. Thus the sum of the two B weights for this one bit in the row inner product is exactly zero, with all other B bits fixed. One-sided terms in the digit formula can be absorbed into the row and column weights. The same proof permits an arbitrary linear phase in z or any additional unit weight that factors between the two assignment sets.

The operator norm of the phase matrix is 2^(k/2), while the ℓ² norms of the two unit-weight vectors are 2^(k/2) each. Therefore the absolute normalized character average is at most 2^(−k/2). This is cancellation across all pairs at once: there is no need to insert a buffer between consecutive blocks.

## 2. Conditional tape mixing and a two-tape gain

On an aligned interval of length 2^r, fix z modulo 2^(m+1), with r≥m+1. The remaining free digits form a consecutive block of length L=r−m−1. For a nonzero character of q modulo 2^m, reduce to h/2^d with d≤m, h odd. Section 1 gives
```
|mean e(h q(z)/2^d)| <= min(1, 2^(−(r−m−d−2)/4)).
```
The constants do not depend on the fixed high prefix or the fixed low digits. In particular a single tape has relative pointwise mixing at powers 2^m below its length to a power strictly less than 1/6, by Fourier inversion.

On a product of two aligned intervals of lengths 2^r and 2^s, fix both input residues modulo 2^(m+1). Then F mod 2^(m+1) depends only on 11(q(x)+q(y)) mod 2^m and on the fixed polynomial and parity terms. Multiplication by 11 is invertible modulo 2^m. Orthogonality on that group, independence of the two tapes, and the product of their character bounds show the **relative** error from uniformity on the odd classes is at most
```
2^m * 2^(−(r+s−4m−4)/4) = 2^(2m−(r+s)/4+1).
```
This holds separately for every pair of input residue classes, so also after any restriction on those classes, such as x odd. The two-tape bound uses the sum of the interval logarithms, not their minimum, except that both tapes must be long enough to fix the requested low digits.

## 3. The actual sharp shell is not a smooth replacement

Take all dyadic boxes [2^r,2^(r+1)) × [2^s,2^(s+1)) meeting the central shell x,y≥D, T/2<F≤T. The elementary bound 0≤p(z)≤z from the [sharp-cut proof](../chanpak_true_cut_fourier_20261115/README.md) shows that on such a box x,y≥D/2, 100xy is comparable to T, and F−100xy=O(T/D), uniformly. In particular r+s=log₂ T+O(1).

Fix a small positive ε<η/2 and subdivide each box into a grid of aligned cells whose side lengths are divided by 2^j, with j=floor(ε log₂ T). Every cell has both logarithmic side lengths at least δ log₂ T−j−O(1), and its sum of logarithmic lengths is log₂ T−2j+O(1). Thus §2 applies with a relative power saving uniformly to every cell. Entire cells in the central shell give a lower bound. All cells that meet the shell give an upper bound. The difference in ordinary point count is
```
O(T(2^(−j)+D^(−1)) log T).
```
Indeed a straddling cell lies in an O(T(2^(−j)+D^(−1))) product strip at either cutoff, or within O(D 2^(−j)+1) of x=D or y=D. In a product strip sum row intervals of length O(T(2^(−j)+D^(−1))/x+1) over D/2≤x≪T/D. In a cutoff strip the other coordinate has O(T/D) choices. This proves the asserted difference, without requiring monotonicity of the checksum or enumerating its fluctuating boundary. Both bounds include the same factor 2^(−m) because the inner and outer approximations are unions of **whole** cells to which §2 applies. That avoids paying the modulus in the boundary error.

On a shell subinterval (U,V] of length Δ≥T^(1−κ), the ordinary mass is asymptotic to its hyperbolic area, of order Δ log T: the area derivative for 100xy at shell level w is (1/100)log(w/(100D²)), and a sum of row errors and checksum perturbation costs O(T log T/D). The same mass order after x odd follows by summing every other x row (or by rectangle approximation). Shrink κ below the cell saving and ε/2. The same two-sided count now gives a relative power error on every such subinterval. All estimates apply to fixed physical addresses, not a randomized surrogate.

## 4. These stronger local laws still do not force density one

This limitation can be shown on actual presentations. Fix an exponent η'<min(δ,1/8) and a small ε>0 with η'+ε<δ, taking ε small enough for the preceding local theorem. On each dyadic shell independently retain whole output fibers with probability (log T)^(-β), where 0<β<1. The selected support has density zero almost surely, by concentration on the number of selected integers per shell and geometric summation. The mean selected presentation count per target across the full shell still diverges.

More strongly, simultaneously for every subinterval of length at least T^(1−ε) in that shell and every odd class modulo 2^(m+1) with 2^m≤T^η', the retained presentation mass is its original mass times (log T)^(-β), with relative error tending to zero. To prove this, write M for the original mass in the query. The theorem gives M≫T^(1−ε−η')log T. The maximum size of an output fiber is O(T/D): fix the larger coordinate and use the strict shorter-coordinate monotonicity in the [sharp-cut proof](../chanpak_true_cut_fourier_20261115/README.md), with O(T/D) choices for the larger coordinate. Thus the sum of the squared ranges in Hoeffding's inequality is O((T/D)M). Asking for relative error T^(−c) on the retained mass gives an exponent at least a constant times
```
T^(δ−ε−η'−2c) (log T)^(1−2β).
```
Choose c>0 with δ−ε−η'−2c>0. There are only polynomially many simultaneous queries (integer cutoffs and residue classes), so the exception probabilities are summable over the dyadic shells. The same random thinning can preserve the [relative full-circle Fourier profile](../chanpak_near_floor_crossing_20261115/DENSITY.md). This is a foil for local information, **not evidence of sparse support of the unthinned carrier**.

There is consequently no dyadic residue obstruction even for polynomially many classes, locally in mesoscopic absolute windows and within the odd-x subfamily. The missing bound concerns distinct presentations of the **same integer**, not just congruences. To exploit the new local precision for coverage, one needs a construction or collision estimate that retains the coarser quotient as well.

[Replay](dyadic.py) checks the integer digit formula, the triangular rank argument on 192 templates, and a literal carrier-box residue histogram. The theorem follows from the preceding inequalities, not from those diagnostics.
