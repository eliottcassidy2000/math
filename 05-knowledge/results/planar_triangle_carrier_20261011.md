# Three branch points obstruct rational times of monomial source repairs

2026-10-11. **PROVED, using CITED curve theory:** the generic-curve genus
formula and nonrationality at every nonzero time for the family below.
**FINITE-EXACT:** independent chart, ramification and hostile controls.
Section J5 rederives the weighted-homogeneous part of an already stronger
repository theorem. **SCOPE CORRECTION, 2026-10-11:** arbitrary polynomial
universal carriers were already settled by the independently audited
[September7 theorem](planar_jc48_sep07_universal_genus.md); the
[September8 cusp-ideal theorem](planar_jc48_sep08_cusp_ideal.md) extends the
nonrational-time obstruction to every nonconstant member of K+D K[p,y].
The first checkpoint of this note incorrectly called mixed weighted
degrees open. Its formulas survive as a rederivation and quantitative
illustration, not a new enlargement of the known nonrational family.
This packet has a self-audit, not an independent peer audit or canon
promotion. **OPEN:** compositions of different flows, useful rational
replacements outside the settled carrier, and planar JC.

## The precise connection to the shape proposal

A triangle and a circle can represent the same topological boundary while
retaining different face data. For the current planar Jacobian target,
there is a more useful exact triangle: a sphere with three marked branch
points. The branch orders determine whether its orbifold geometry is
spherical, Euclidean or hyperbolic. Here it is always hyperbolic, and the
covering curve has genus at least two. This turns the shape suggestion
into an obstruction to a particular family of repair operations.

The source is the completed source-preserving Hamiltonian carrier; the
target is its primitive generic invariant curve; the map is restriction
to the invariant function field. It preserves the invariant and the
flow's nonzero displacement, but forgets affine boundary and polynomial
presentation. The necessary sidecars are primitive power, branch
multiplicities and the completion-to-source comparison. The cheapest
hostile is a nonprimitive invariant: a larger Newton triangle can describe
the same primitive curve, not a larger genus.

The closest proved mechanism is the
[September7 universal exceptional-cover proof](planar_jc48_sep07_universal_genus.md),
read together with its [independent audit](planar_jc48_sep07_universal_genus_audit.md)
after the overly narrow initial inheritance pass was corrected. Its
antecedent is the [September6 family](planar_jc_long_20260906_nonrational.md).
The genus-zero rational Hamiltonian flow x^2t is the canonical hostile. The corrected near miss is
non-local-nilpotence=>nonrationality. The underused sidecar is the
three-point Kummer cover. Board: **primitive invariant / ramification /
completion / displacement / finite specialization / source**.
META-PATTERNS used: **Inventory retained power classes before classifying
base-curve genus** and **Correct the object before sharpening the technique**.

The incoming [Collatz integration](collatz_fundamental_integration_20261009.md)
provides a useful comparison: a faithful formal carrier still needs a
marked ordinary source and terminal. Here an all-order formal repair still
needs rational or polynomial specialization. The barriers are different;
neither theorem is transported to the other problem merely by analogy.

## J1. Statement

Over a characteristic-zero field K use the inherited coordinates

\[
 u=x^2t,\quad p=t(1+x^2t),\quad y=xtp,\quad D=p^3-y^2,
 \qquad\{F,G\}=F_xG_t-F_tG_x.
\]

Let a>=2, b>=0, c>=1 be integers, put

\[
                  I=p^a y^b D^c,
\]

and let f be any nonconstant polynomial over K. The scalar-time completed
Hamiltonian flow E_lambda=exp(lambda{-,f(I)}) exists by the inherited
[completed Hamiltonian theorem](planar_jc_long_20260906_hamiltonian.md),
because f(I)-f(0) belongs to p^2D K[p,y].

**Theorem.** For every nonzero lambda in K, the two images E_lambda(p)
and E_lambda(y) cannot both belong to K(p,y). Consequently this operation
cannot specialize to a polynomial source map at that time.

The parameters b and c are unrestricted in the stated ranges. The
September6 result is the subfamily b=0,c=1; the September7 theorem already
contains the entire present family and arbitrary polynomial multipliers.
No claim of novelty or literature priority is made.

## J2. Strip the power before computing the curve

Set d=gcd(a,b,c), (a0,b0,c0)=(a,b,c)/d and
J=p^{a0}y^{b0}D^{c0}; thus I=J^d. The flow fixes J. It is this primitive
invariant, not a generic level of I after splitting, that defines the
connected geometric curve.

For the next calculation rename (a0,b0,c0) as (a,b,c). Thus a,c>=1,
b>=0 and gcd(a,b,c)=1, even though the original source carrier required
the original a>=2. On the dense torus set

\[
 r=y/p,\qquad z=y^2/p^3,\qquad p=r^2/z,\quad y=r^3/z.
\]

Writing w=2a+3b+6c and k=a+b+3c gives the exact identity

\[
 J=r^w z^{-k}(1-z)^c.
\]

Over F=K(J), the generic curve is therefore the Kummer cover

\[
                  r^w=J\,z^k/(1-z)^c.                 \tag{1}
\]

The exponents of the divisor on the z-sphere at0,1,infinity are
k,-c,c-k. The cover is geometrically connected: their gcd with w is
gcd(w,k,c)=gcd(a,b,c)=1, so the right side is not a proper power of any
prime divisor of w over the geometric rational-function field. This also
justifies its degree w. No generic irreducibility is inferred solely from
the appearance of a trinomial.

In characteristic zero, a local equation v^w=z^m times a unit has
gcd(w,m) branches, each of ramification index w/gcd(w,m). One can see this
after extracting the unit's formal w-th root and using the coprime
parametrization of v^{w/g}=z^{m/g}. There are no other ramification points.
Put

\[
 g_0=\gcd(w,k),\quad g_1=\gcd(w,c),\quad
 g_\infty=\gcd(w,k-c).
\]

[Riemann–Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) then gives

\[
 \boxed{\displaystyle g=\frac{w+2-g_0-g_1-g_\infty}{2}.}            \tag{2}
\]

This use of the cited theorem is for a separable map of smooth proper
geometric curves; the displayed affine model need not itself be smooth
at its missing boundary.

**Uniform lower bound.** If b=0 then g0=k=a+3c, g1<=c and
g_infinity<=2c, so g>=(a+2)/2. Since genus is integral and a>=1, g>=2.
If b>0 then

\[
 g_0=\gcd(w,k)\le b,\qquad g_1\le c,\qquad
 g_\infty=\gcd(w,k-c)\le b+2c.
\]

Thus g>=(2a+b+3c+2)/2>=4. In particular every permitted primitive curve
has genus at least two.

The three branch orders are e_i=w/g_i, with

\[
                  2-2g=w(1/e_0+1/e_1+1/e_\infty-1).     \tag{3}
\]

Their reciprocal sum is strictly below1: the triangle signature is
hyperbolic. This is a precise sphere/triangle/curvature correspondence.
The familiar5/6/7 transition is the special family(2,3,m): its reciprocal
sum minus1 is1/30,0,-1/42 at m=5,6,7. The three angles pi/e_i thus have
sum greater than, equal to, or less than pi. Here the branch orders are
derived from the actual polynomial rather than assigned by resemblance.

| Original (a,b,c) | Primitive triple | Branch orders at0,1,infinity | Genus |
|---|---|---|---:|
| (2,0,1) | (2,0,1) | (2,10,5) | 2 |
| (3,0,1) | (3,0,1) | (2,12,12) | 3 |
| (2,1,1) | (2,1,1) | (13,13,13) | 6 |
| (2,0,2) | (1,0,1) | (2,8,8) | 2 |

At c=1, formula(2) becomes

\[
 g=\frac{2a+3b+7-\gcd(a+3,b)-\gcd(a,b+2)}2.
\]

It equals the number of interior lattice points in the Newton triangle
with vertices(0,0),(a+3,b),(a,b+2). For b=0 it recovers
ceil(a/2)+1 from the inherited hyperelliptic proof. For c>1 the repeated
edge factor D^c and possible invariant power must be retained; a raw
interior-point count need not be the geometric genus. Rounding the
triangle into an unlabelled circle loses exactly this arithmetic boundary.

## J3. A nonzero time has infinite order

Use the inherited logarithmic chart

\[
 p=s^2+\tau,\quad y=sp,\quad D=\tau p^2,
 \qquad\{F,G\}=\tau(F_sG_\tau-F_\tau G_s).
\]

Return to the **original**, unnormalized a,b,c and w=2a+3b+6c. Then

\[
                  I=\tau^c s^b(s^2+\tau)^{a+b+2c}.
\]

Write f(z)=f(0)+f_h z^h+... with h>=1 and f_h nonzero. The derivation
delta={-,f(I)} raises tau valuation by at least ch on K(s)((tau)).
Its s coefficient has valuation at least ch and its tau coefficient
at least ch+1. Thus its exponential converges for every scalar lambda,
has inverse E_{-lambda}, and obeys the additive group law. Direct
differentiation gives

\[
 \delta p=2ch f_h\,s^{h(w-2c)+1}\tau^{ch}
                 +O(\tau^{ch+1}).
\]

All higher iterates have strictly larger tau order, so

\[
 E_\lambda(p)-p
   =2\lambda ch f_h s^{h(w-2c)+1}\tau^{ch}
                   +O(\tau^{ch+1}).                    \tag{4}
\]

For every positive integer N, replacing lambda by N lambda leaves this
coefficient nonzero. Consequently E_lambda has infinite order already on
p. The generic variable s is not replaced by a special numerical value.

If both images were rational, E_lambda would restrict to an injective
F=K(J) field map of K(p,y) to itself. The
[function-field/curve correspondence](https://stacks.math.columbia.edu/tag/0BY1)
turns this into a nonconstant selfmap of the smooth proper generic curve.
Since its genus is at least two, Riemann–Hurwitz forces degree1. Its
[automorphism group is finite](https://stacks.math.columbia.edu/tag/0DST):
in characteristic zero the tangent sheaf has negative degree and no
sections, yielding a zero-dimensional finite-type automorphism scheme.
Thus some positive iterate fixes p, contradicting(4). This proves J1.
Rationality of the inverse was not assumed in making this argument.

The source specialization uses the same injective comparison as the
inherited proof: K[[s,tau]] -> K[x][[t]], s->xt, tau->t. The monomial
s^i tau^j maps to x^i t^{i+j}, so no coefficients collide. Images of p,y
belong to this domain and the derivations intertwine. Clearing a proposed
rational denominator therefore pulls source rationality back to Laurent
rationality, already excluded. This does not identify the entire actual
D-adic completion with the Laurent completion.

## J4. Recovered classification: every multiplier depending only on p

This formula is equivalent to Theorem A' of
[the September7 carrier classification](continuing9_20260907_flow_carrier_genus.md).
If H=p^2 R, then deg H=deg R+2 and o(H) is the degree of the squarefree
part of R. This exactly identifies the two formulas, rather than treating
different notation as a new result.

Let H(p) be a nonzero polynomial divisible by p^2 and let
I=H(p)D. The same nonrationality conclusion holds for f(I) at every
nonzero time. This covers all nonzero univariate multipliers R(p) in the
universal carrier p^2D R(p), not just monomials.

Put h=deg H>=2 and let o(H) be the number of distinct roots of H having
odd multiplicity over an algebraic closure. Over K(I) the invariant
curve has equation

\[
                     y^2=p^3-I/H(p).
\]

The numerator p^3H(p)-I has h+3 simple roots and no common root with H.
Simplicity follows because I is transcendental, whereas every critical
value of the fixed polynomial p^3H lies in the algebraic closure of K.
The finite poles at the roots of H branch exactly at odd multiplicity;
infinity has pole order3 and branches. Thus

\[
                      g=(h+o(H)+2)/2\ge2.               \tag{5}
\]

The squarefree generic numerator also proves geometric integrality. The
leading displacement for the basic Hamiltonian is
{p,I}=2s tau p^2 H(p), nonzero. For a nonconstant f the same valuation
argument as J3 works, so the curve proof applies verbatim. This genus
calculation allows more than three branch points; it shows why a
nonmonomial multiplier can increase complexity without offering the
genus-zero/one escape sought in the previous board.

## J5. The whole weighted-homogeneous carrier has the same obstruction

The preceding branch calculation extends to **every nonzero polynomial I
in p^2D K[p,y] homogeneous for weights wt(p)=2, wt(y)=3**, and then to
f(I) for any nonconstant f. Thus this lane is not restricted to generators
that are monomials in p,y,D.

Over an algebraic closure, weighted factorization writes

\[
 I=C p^a y^b\prod_{i=1}^s(p^3-\beta_i y^2)^{c_i},       \tag{6}
\]

where a>=2,b>=0, all c_i>=1, the nonzero beta_i are distinct, and one
beta_i=1. To prove this factorization, remove the minimum p and y powers.
Equal weighted degrees then make all p-exponent differences multiples of3
and all y-exponent differences multiples of2. The remaining polynomial is
an ordinary binary homogeneous polynomial in p^3,y^2 with both endpoints
nonzero, hence factors as claimed. Work over the algebraic closure only
strengthens the rationality obstruction for the original field.

Strip d=gcd(a,b,c_1,...,c_s) and use the primitive invariant J. Rename its
exponents a,b,c_i and put Csum=sum c_i, w=2a+3b+6Csum, k=a+b+3Csum.
Now a>=1 and the Kummer model is

\[
           r^w=J z^k\Big/\prod_i(1-\beta_i z)^{c_i}.
\]

The same gcd argument gives a connected geometric cover. Its branch
divisor orders are k,-c_1,...,-c_s,Csum-k, so

\[
 2g-2=s w-\gcd(w,k)-\gcd(w,k-Csum)-\sum_i\gcd(w,c_i).    \tag{7}
\]

For b=0, the subtracted terms are bounded by a+3Csum,2Csum,Csum.
Consequently 2g-2 >= (2s-1)a+6(s-1)Csum>0.
For b>0 they are bounded by b,b+2Csum,Csum, giving
2g-2 >= s(2a+3b+6Csum)-2b-3Csum>0. Thus g>=2 in every case.

In the logarithmic chart, the factor with beta_i=1 has positive tau
order, while every other factor is a unit over K(s). Therefore the
leading-displacement proof J3 applies with a nonzero scalar multiplier:
every nonzero time has infinite order. Its fixed primitive invariant has
finite automorphism group. Every nonzero time is again nonrational.

**A constructive hostile for the next route.** Let I1=p^2D and I2=p^2yD.
Their Hamiltonian bracket is

\[
             \{I_1,I_2\}=-p^2D^2(5p^3-2y^2).           \tag{8}
\]

Its primitive data have s=2, w=22,k=11, Csum=3, and branch-cycle gcds
11,2,2,1. Formula(7) gives genus15. Thus taking the first Lie bracket
does not escape the high-genus family. This is a statement about the
Hamiltonian(8); the group commutator of the two original flows has higher
terms and need not fix(8). No obstruction to that entire composition is
inferred from its first Baker–Campbell–Hausdorff term.

## Scope, controls, and next constructive direction

The identity time and constant f are equality boundaries. The Hamiltonian
x^2t has the rational non-LND flow

\[
 x\mapsto x/(1-\lambda x),\qquad
 t\mapsto t(1-\lambda x)^2,
\]

so non-LND alone remains insufficient. The finite controls retain this
hostile and the nonprimitive example(2,0,2). The earlier universal and
cusp-ideal theorems already cover every single polynomial generator in
their stated ideals, including arbitrary sums. Compositions that change
invariants and cancellations between operations require their own audit.

The constructive next target is consequently outside this already closed
lane: a rational Hamiltonian outside K+D K[p,y], or a composition changing
invariants. It must still pay actual source preservation and later finite
compatibility. The September8 theorem shows why genus alone is insufficient
at the boundary: D has elliptic generic fibres, but its vector field has
poles at p=0; a commuting selfmap must retain that pole divisor.

The present z=y^2/p^3 chart and the earlier w=p^3/y^2 exceptional chart
are related by w=1/z, while both use r=v=y/p. Thus the branch data are
the same sphere with zero and infinity exchanged. The earlier proof
additionally embeds a compact bordered cover into an actual nearby fibre,
paying the local-to-global step for arbitrary polynomial tails. Replacing
that embedded piece by an unlabelled triangle or circle would discard the
component gcd and branch multiplicities that make the obstruction work.

The [rational-carrier continuation](planar_rational_carrier_saturation_20261011.md)
tests the suggested change of domain directly: three polynomial source
responses force any rational Hamiltonian back into K+p^2D K[p,y]. A
rational replacement must therefore weaken the universal source condition
or change the source algebra, retaining a proof of the required final return.

The [checker](../../04-computation/experiments/planar_triangle_carrier_20261011.py)
and [output](planar_triangle_carrier_20261011.json) compare two derivations
of the chart, verify exact ramification counts for the stated bounded
universe, and check the inherited subfamily and rational hostile. They
are supporting controls, not a proof by finite enumeration.
Normal and optimized runs agree byte for byte on31,404 always-active exact
gates. The general formulas are proved above; the bounded controls do not
stand in for a classification of all carriers.

```text
python3 -B 04-computation/experiments/planar_triangle_carrier_20261011.py
python3 -B -O 04-computation/experiments/planar_triangle_carrier_20261011.py
```
