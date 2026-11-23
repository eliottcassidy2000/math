# The six branches become a genus-two curve

**PROVED scoped.** Outside inspiration: the branch-label suggestion
and its genus-two model came from the parallel background pass.
Here the arithmetic scale is derived from the tree, and the
Frobenius claim is strengthened to an equality of operators.

## Arithmetic rather than an arbitrary relabelling

For an exit \(H\) at relative depth \(d\) below 7, use the
centered height
\[
 x=\frac32\left(\frac{H}{3^d}-7\right).
\]
In boundary order this is
\[
\begin{array}{c|rrrrrr}
H&1&19&61&63&65&23\\ \hline
x&\infty&-1&-1/3&0&1/3&1\\
x\bmod5&\infty&4&3&0&2&1
\end{array}
\]
The exterior is deliberately marked at infinity. The values
distinctly reduce at every prime other than 2 and 3, as is immediate
on comparing the five finite rationals by cross multiplication.
Only at \(p=5\) do they exhaust the projective line. Ordinary
reduction of the *path counts* modulo 5 does not give this map.

The rooted mirror is \(x\mapsto-x\). The hexagon half-turn becomes
\(x\mapsto3/x\) in these finite labels: it pairs \(1,63\),
\(19,65\), and \(61,23\). These are statements about two
specified symmetries; no claim that the numerical cyclic order
is itself an orbit of a Möbius rotation is needed.
The translation \(x\mapsto x+1\) fixes the exterior and
cycles the exits \((63,23,65,61,19)\). This is not the cyclic
order of the hexagon. The two outer exits are quadratic residues;
the two inner noncentral exits are nonresidues.

## An exact Paley–Frobenius identity

The six rationals are the ramification points of the smooth
characteristic-zero double cover
\[
 C:\ y^2=x(x^2-1)(9x^2-1).
\]
At 5 its smooth reduction is \(C_5:y^2=x-x^5\).
Let \(\tau(x,y)=(x+1,y)\) and let \(F(x,y)=(x^5,y^5)\) on
\(C_5\). Write \(T=\tau^*\), \(\Phi=F^*\) on
\(H^1_{\rm et}(C_{5,\bar{\mathbf F}_5},\mathbf Q_\ell)\), \(\ell\ne5\).
Then on this four-dimensional vector space
\[
 \boxed{\ \Phi=-\bigl(T+T^{-1}-T^2-T^{-2}\bigr)\ }.
\]
This is not just a coincidence of eigenvalues.

For \(a\in\mathbf F_5\), a fixed affine point of \(\tau^a F\)
has \(x-x^5=a\) and \(y\in\mathbf F_5\), with \(y^2=a\).
There are five distinct solutions for \(x\). Including infinity
gives \(6+5\chi(a)\) fixed points, with \(\chi\) the quadratic
character and \(\chi(0)=0\). All are counted with multiplicity
one because the differential of the composite vanishes. The
curve trace formula then gives
\(\operatorname{tr}(T^a\Phi)=-5\chi(a)\).
After adjoining a primitive fifth root \(\zeta\), take the
Fourier transform:
\[
 \operatorname{tr}(\Phi\mid T=\zeta^k)
 =-\sum_a\chi(a)\zeta^{-ka}.
\]
For each of the four nonzero \(k\), this is the nonzero
quadratic Gauss sum with a sign. Each such eigenspace must
exist, so all four are one-dimensional and they exhaust
\(H^1\). Since \(\chi(-a)=\chi(a)\), the displayed polynomial
in \(T\) has exactly those scalars, proving the identity.
The Gauss sum is \(\pm\sqrt5\) (square 5), so each sign in
\(\{\sqrt5,-\sqrt5\}\) occurs twice. This is precisely the
signed difference of the two distance orbits of a five-cycle
on its nonconstant modes. Ordinary adjacency has eigenvalues
\((\sqrt5-1)/2\) and \((-\sqrt5-1)/2\); the signed difference,
not ordinary adjacency, gives \(\pm\sqrt5\). This uses the
same distance-mode split as the repository's
[golden fold atlas](../collatz_coincidence_atlas_20260927.md),
but now the signed five-cycle is the actual operator, not just
a comparison of characteristic polynomials.

As consequences, the characteristic polynomial is
\(t^4-10t^2+25\) and \(J(C_5)(\mathbf F_5)=J(C_5)[2]\) has
order 16: all six rational branch points give the full
two-torsion, and evaluating the polynomial at 1 gives 16.
Independently \(C_5(\mathbf F_5)\) and \(C_5(\mathbf F_{25})\)
both have six points: for \(x\notin\mathbf F_5\),
\(s=x-x^5\) is nonzero and \(s^5=-s\), hence \(s^{12}=-1\)
and it is nonsquare. The Cartier–Manin matrix is zero:
\((x-x^5)^2=x^2-2x^6+x^{10}\) has no terms at indices
\(4,3,9,8\). The resulting superspecial reduction is an instance
of the [known family](https://www.degruyterbrill.com/document/doi/10.1515/jmc-2019-0021/html).
The replay also counts \(C(\mathbf F_7)\) and \(C(\mathbf F_{49})\)
as 8 and 46, giving \(|J(C)(\mathbf F_7)|=48\).
Prime-to-5 torsion over \(\mathbf Q\) injects at 5 and the
5-primary part injects at 7, so **the entire rational torsion
of \(J(C)\) is \((\mathbf Z/2)^4\)**.

Over characteristic zero, put \(X=\sqrt3\,x\),
\(Y=3^{1/4}y\). Then \(Y^2=X(X^4-\frac{10}3X^2+1)\).
The reciprocal involution and its product with the deck
involution have quotients
\[
 v^2=(u\pm2)(u^2-16/3),\qquad u=X+X^{-1},
 \quad v=Y(X\pm1)/X^2.
\]
They both have \(j=54000\). For the plus sign translate the root
to write \(v^2=z(z^2-4z-4/3)\). Quotienting by the point
\((0,0)\) gives the standard 2-isogenous model with quadratic
coefficients \(8,64/3\) and thus \(j=0\).
This is a separate elliptic-curve bridge to the earlier
supersingular work, not an equivalence with its Ogg-prime test.

## Mark infinity: six points, five labels, fifteen pairings

The [standard two-torsion model](https://people.maths.ox.ac.uk/flynn/arts/art31.pdf)
is the vector space of even subsets of the six branch points modulo
complement. It has intersection parity as symplectic form, and its
15 nonzero classes are represented by the 15 two-subsets. Five
two-subsets contain \(\infty\); the other ten are finite pairs.
The graph on the ten pairs whose distinct vertices are adjacent
when disjoint is precisely the Petersen graph, the Kneser graph
\(KG(5,2)\). The full unmarked 15-pair graph does not have only
these ten vertices.

Every partition of the six points into three pairs gives a
maximal isotropic subgroup with three nonzero points: two
disjoint pairs add to the third by complement. Conversely
any two independent orthogonal pair labels are disjoint.
Dropping the pair through \(\infty\) gives a **bijection**
from the fifteen partitions (Richelot kernels) to the fifteen
Petersen edges. These split into five triples according to
the finite partner of \(\infty\): choosing that partner leaves
three ways to pair four points. Thus the five labels carry a
specific five-colouring of the Petersen edges into matchings
(each class leaves the four pairs containing that label uncovered).
The sixth label is what makes each edge determine a full partition.

This brings the actual seeds into the prior
[Richelot–Frobenius and outer six-set](../session97_overnight_arithmetic_carriers_20261029/supersingular/RICHELOT_FROBENIUS_AND_FANO.md).
Forgetting the mark yields the older \(K_6\) description:
15 matchings, and six pentads of five matchings which partition
the edges of \(K_6\), with every matching in two pentads.
Projective \(\mathbf P^1(\mathbf F_5)\) also gives the six
Sylow-5 subgroups of \(A_5\): under \(PSL_2(\mathbf F_5)\cong A_5\),
each projective point is the fixed point of exactly one unipotent
subgroup of order 5. Fixing infinity leaves the cyclic five-point
translation action. These are different six- and five-set actions,
not correspondences with the fifteen Ogg primes or with physical
tournament flips.


On the five finite labels there are \((5-1)!/2=12\) undirected
pentagons in \(K_5\), paired by edge complementation into six.
Each pentagon's edges are a 5-cycle in the Petersen disjointness
graph; the edges crossing to the complementary pentagon form
a perfect matching of Petersen. The crossing matching recovers
the unordered pair of pentagons by deletion. This recovers six
objects from the five marked exits without choosing a projective
coordinate, and realizes the corresponding Sylow-5 subgroup as
the cyclic rotation group of either pentagon. It relates directly
to the [earlier Petersen and icosahedral models](../planar_shadows_petersen_switching_20261116/README.md),
where an antipodal fold gives six axes and a fivefold rotation
fixes one of them. By contrast neither the Petersen edge relation
nor icosahedral adjacency is the tournament single-arc relation.
