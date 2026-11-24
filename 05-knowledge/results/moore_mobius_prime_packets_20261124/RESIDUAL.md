# A second Möbius readout from the residual graph

This extends [the pentagon-array calculation](PENTAGON.md) from minimum exterior height to **every** cycle type. All degree-57 graph statements are conditional. A contact deck consists of the exterior pentagram lift on $X$ and the five-column array on $Y$ with the prescribed missing matching for every pair of columns. Write $A$ for the adjacency on $X$, $N$ for the $Y$-by-$X$ incidence matrix, $B$ for the proposed adjacency on $Y$. Sizes are 275 and 2,970. In this note $C_j(B)$ counts unoriented simple cycles, not primitive walks.

## Recover the prime packets from residual walks

Orient the exterior cycles consistently with the order of the five pentagram columns and let $Q$ advance one step, so $A=Q+Q^{-1}$. On the five column-centered spaces,
$$
 N^{\mathsf T}N=55I-(I+Q+Q^{-1}+Q^2+Q^{-2})=56I-A-A^2.
$$
This follows from the pair-of-columns rule: the diagonal blocks are $54I$, and an off-diagonal block is the all-one block minus the matching for the unique one- or two-step exterior route. The five-term sum has eigenvalues of absolute value at most five since $Q$ is orthogonal. The Gram is therefore positive definite on 270 column-centered dimensions and the column-constant sector contributes rank one. Thus $\operatorname{rank}N=271$ for **every** exterior type.

For a completed residual graph,
$$
 BN=J-N(I+A),\qquad
 B^2+B=56I+J-NN^{\mathsf T}.
$$
On the 270 nonconstant contact-image dimensions, $B$ acts as $-I-A$ after transport by $N$; the five discarded column-constant directions of $A$ carry a single $C_5$ spectrum. The constant direction of $B$ has eigenvalue 52. The image trace is $52-270=-218$ since the exterior adjacency, both before and after removal of the column-constant sector, has trace zero. On the remaining 2,699 dimensions $B^2+B=56I$, and trace fixes the multiplicities of 7 and $-8$ at 1,454 and 1,245. In particular for $P_j(s)=\det((s+1)I+A(C_j))$ and exterior cycle census $c_\ell$,
$$
 \chi_B(s)=(s-52)(s-7)^{1454}(s+8)^{1245}
       \frac{\prod_\ell P_{5\ell}(s)^{c_\ell}}{P_5(s)}.
$$
This recovers the previously calculated minimum spectrum when $c_1=55$.

One can decode from traces without numerical eigenvalues. Use the integer polynomials $D_0=2,D_1=x,D_{k+1}=xD_k-D_{k-1}$. Since $D_k(A)=Q^k+Q^{-k}$, define
$$
 T_m=\frac{
 \operatorname{tr} D_{5m}(-I-B)
  -D_{5m}(-53)-1454D_{5m}(-8)-1245D_{5m}(7)+10}{10}.
$$
Then $T_m=\sum_{\ell\mid m}\ell c_\ell$; the final $+10$ restores the removed $C_5$ quotient. Apply ordinary Möbius inversion to get $c_\ell$, then form the actual prime packets of the exterior tournament height. Thus one can cross-check exact integer divisibility **against independent residual walk data**.

There is a simple first-detection rule even for ordinary traces. Up to step $5m-1$, a $5m$-cycle has the same per-vertex return counts as the infinite line (for $I+A$ as well as $A$). At step $5m$, each of its vertices has two extra paths winding once around. Consequently, for two types whose first differing count is $c_m$, the first differing moment is $5m$, with difference $(-1)^{5m}10m\,\Delta c_m$. The fixed exterior size means their line contributions are identical. For example $c_2$ first appears at moment 10:
$$
 \operatorname{tr}B^{10}=144556853479042730+14530c_1+20c_2.
$$
The coefficients count displacements 5 and 10 in ten steps with increments $-1,0,1$, with the $C_5$ quotient subtracted exactly once.

## A completion certificate and its spectral tolerance

Assume only the contact-deck rules, a simple symmetric 52-regular $B$, and $BN=J-N(I+A)$, not Moore completion. Set
$$
 E=B^2+B-56I-J+NN^{\mathsf T}.
$$
From the Gram identity and $NN^{\mathsf T}\mathbf 1=270\mathbf1$, $E$ is zero on $\operatorname{im}N$; on its orthogonal complement $E=q(B)$ for $q(s)=(s-7)(s+8)$. The other blocks of the complete Moore identity already follow from the deck: central-to-exterior uses the prescribed one- or two-step routes, central-to-$Y$ uses the one contact per column, and exterior-to-exterior uses the forbidden pairs plus the matching contact rule. Thus $E=0$ is the exact outstanding condition.

Take the formal target spectrum above, even if $B$ is not yet a solution. The degree, zero diagonal and order ensure agreement with its moments of order zero through two; they also agree fully on the contact image. The formal moments through order four are $2970,0,154440,0,15907320$, independent of the lift, because a cycle of length at least five has no winding walks before order five. Expanding $q^2$ on the complement and using the regular fourth-moment identity yields
$$
 \operatorname{tr}E^2
 =\operatorname{tr}B^4-15907320+2\operatorname{tr}B^3
 =8C_4(B)+12C_3(B).
$$
It follows that eliminating triangles and squares **within** $Y$ certifies the whole completion, once the contact-deck and linear compatibility checks pass.

A second, robust certificate uses only a spectral interval. Here $E$ is symmetric integral, diagonal zero (using the degree and the five contacts of each row), with $E\mathbf1=0$. Any nonzero matrix with these properties has largest eigenvalue at least $\sqrt2$: an entry of magnitude at least two gives a two-site witness; an induced three-site path in the support gives a three-site witness. If neither occurs, each nontrivial support component is a signed clique. A positive-product triangle has eigenvalue two on its principal minor, up to sign switching. If all triangle products are negative, the block is of the form $I-ss^{\mathsf T}$ for a sign vector $s$, after choosing signs from a pivot (also for a two-site component). Such a block of size greater than one cannot have zero row sums, since $s_i\sum_j s_j=1$ for all $i$ would force all signs equal and block size one. This proves the lower bound by interlacing. Consequently it suffices to put the spectrum of $B$ on $(\operatorname{im}N)^\perp$ inside
$$
 [-8-\epsilon,7+\epsilon]\qquad\text{with }\epsilon(15+\epsilon)<\sqrt2.
$$
Indeed the maximum of $q$ on that interval is $\epsilon(15+\epsilon)$, so the remaining integer defect must vanish. For example the exact tolerance $\epsilon=9/100$ works. This is a certificate if the interval is verified rigorously, not a claim that a spectral relaxation finds a candidate.

## Cross-checks at cycle lengths five and six

For any completed $B$ the characteristic formula gives
$$
 C_5(B)=36383149-c_1,\qquad
 C_6(B)=1518993850+5c_1 .
$$
For exponents 0 through 6, the return counts of $I+A$ on the exterior line are $1,1,3,7,19,51,141$. Each exterior pentagon adds 10 to the fifth trace and 60 to the sixth, and the removed quotient also has these contributions. Multiply by $(-1)^k$ for the shifted eigenvalues and add the 52, 7 and $-8$ sectors. Then use $\operatorname{tr}B^5=10C_5(B)$ and $\operatorname{tr}B^6=2970\cdot52(5\cdot52^2-6\cdot52+2)+12C_6(B)$. For the second identity a closed length-six walk on a girth-five graph is either entirely backtracking along a tree or a simple hexagon, counted by twelve rooted traversals. Counting tree excursions at length six gives the regular tree term (the five Dyck path shapes).

The Hoffman–Singleton case forces the minimum independently of the finite census. There $t=5$. A longer exterior cycle $5\ell$ has a phase $2\pi/(5\ell)$ that is not column-constant; the Gram $5I-(I+Q+Q^{-1}+Q^2+Q^{-2})$ is positive on that phase. Its contact image would give a residual eigenvalue $-1-2\cos(2\pi/(5\ell))<-2$, impossible for a 2-regular graph on 20 sites. The exterior therefore is five pentagons, the contact rank is 17, and the three complementary eigenvalues must all be 2 (the other root $-3$ is too small). Including the constant gives four copies of eigenvalue 2, hence four residual cycles; since their lengths are at least five and total 20, all four are pentagons. The degree-57 residual is 52-regular, so this easy spectral exclusion of long exterior cycles does not carry over.

The $\operatorname{tr}E^2$ and interval criteria are finite analogues of using a definite-sign global bound to force local multiplicity constraints, as discussed in [the previous packet](README.md). They do not use RH. See **residual-check.py** for the constants and the known-graph calibration.
