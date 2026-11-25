# Exponent-61 dyadic pruning and the first order-three sign case

## Dyadic restrictions on the exact height

Write the [two-tail heights](VARIABLE-TAIL-SELECTOR.md) as \(H=2(avu+bv+cu-h(q(v)+q(u)))+k\), with \((h,a,k)=(29,640,111)\) or \((651,17616,2313)\). Here \(q(2n)=q(n)\), \(q(n)\le n/2\), and both linear coefficients exceed \(h/2\). For \(H=E_p+\delta\), \(\delta=-3,-1,+1\), suppose both positive flanks are divisible by \(T=2^r\), \(r\le p-2\), and write \(v=Tx,u=Ty\). Reducing the normalized height modulo \(T\) forces \(q(x)+q(y)\equiv Q_0=h^{-1}((k-\delta)/2)\bmod T\), using the least nonnegative residue. Therefore \(x+y\ge2Q_0\), \(xy\ge2Q_0-1\), and
\[
 2aT^2(2Q_0-1)\le H.
\]
With flank ratio at most two, \((x+y)^2\le(9/2)xy\) strengthens this to \(16aT^2Q_0^2\le9H\). These are absolute size restrictions, not statements about a fixed local modulus.

For \(p=61\), the first rejected common valuations in the six-site chart for the three queries are \(37,37,37\), or \(29,30,29\) in the balanced case. In the eight-site chart they are \(36,36,36\), or \(27,27,27\) if balanced; the central query \(\delta=-1\) is the one not already excluded modulo 3. The failure at a given valuation applies to every higher valuation by divisibility. This does not cover the axes or lower valuations.

## Order three in characteristic three

Apply the [Gewirtz branch-sign constraints](NO-ORDER13-SIGN-SYMMETRY.md), pre-gauging the odd-order symmetry. A row pair within a moving 3-cycle agrees on every fixed site. Each moving block contributes 0 (constant signs) or 2 (nonconstant). Their distance is a multiple of three, and the route lower bound forces it to be at least 12. Thus there are at least six moved cycles. At exactly six, every moving block is nonconstant and the internal blocks are \((-1,+1,+1)\), so all intra-cycle route counts are 10.

The symmetric moving block sum matrix \(S\) over \(\mathbf F_3\) has diagonal \(+1\), every entry \(\pm1\), and \(S^2=0\) (the fixed-site contribution to that constant block is multiplied by three). Its self-orthogonal row space has rank at most two: rank three would give a systematic generator \([I_3\mid B]\) with \(BB^{\mathsf T}=-I_3\), impossible by determinants. Gauge the first row all positive. If a row \(r\) is nonconstant, the only sign rows in the span of \(1,r\) with first sign positive are \(1,r\). At any coordinate where \(r\) is negative, the row there cannot be \(r\) by the positive diagonal; but it cannot be \(1\) by symmetry with \(r\). So \(S=J\), the all-ones matrix.

Each block now has one minus sign, at a lag \(\ell_{ij}=-\ell_{ji}\), with \(\ell_{ii}=0\). Write \(L=(\ell_{ij})\) and \(K_{ij}=1+\binom{\ell_{ij}}2\) in \(\mathbf F_3\), using representatives \(0,1,2\). Let \(A\) be the fixed block and \(P\) the 39-by-6 profiles. Expanding the moving polynomial \((z-1)^2+z^\ell\) modulo \((z-1)^3\) gives
\[
 A^2=0,\quad AP=-PJ,\quad LJ+JL=0,\quad L^2+JK+KJ+P^{\mathsf T}P=0.
\]
The columns of \(P\) are distinct even up to sign. Indeed, for two moving orbits the moving-coordinate distances at the three relative shifts average 8. Equal profiles give some total distance at most 8; opposite profiles add 39, forcing some odd multiple of three at least 51. Both violate the route bounds.

More precisely, for two different moving orbits \(i,j\), the multiplicities \(n_s\) of the three differences \(s=\ell_{ik}-\ell_{jk}\) as \(k\) runs over the six cycles obey \(n_0+n_1+n_2=6\). The moving-coordinate distance at shift \(s\) is \(12-2n_s\) (orient the shift to agree with this difference). The fixed part is constant, so the three \(n_s\) are congruent modulo 3. Up to permutation they are \((0,0,6),(0,3,3),(1,1,4),(2,2,2)\). The direct negative lag has multiplicity at least two, from \(k=i,j\). These restrictions reduce the 15-phase search; they do not exclude six or more moving cycles.

Replay the integer thresholds and multiplicity list with the adjacent JavaScript file.
