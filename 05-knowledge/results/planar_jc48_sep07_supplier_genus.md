# The actual completed-response supplier has generic genus twenty-seven

**Status: PROVED + FINITE-EXACT; independent analytic/source audit PASS.**
This proves nonrationality of every nonzero scalar time for the named
completed-response Hamiltonian, by a degree-seven generic-curve argument.
It does not revoke its exact finite-row action, assert termination, or
classify all polynomial universal carriers.

## 1. Actual object, inheritance, and precise class

Over any characteristic-zero field K, set D=p³-y² and consider

```
I=p² y³ D(A(p)+b y²),  b in K*,  deg A=m>=4.            (1)
```

The actual supplier in [the completed-response theorem](planar_jc48_sep07_completed_response.md)
has b=532468/6235325 and

```
A=12015/15962432+(28467/3990608)p+(2499/498826)p²
  -(4806/1247065)p³+(639368/11223585)p⁴.                 (2)
```

Thus its Sstar=p²D Rstar is exactly (1), with m=4. No coefficient is
replaced by a generic surrogate in this identification.

The inherited mechanism is the completed universal source carrier and
its fixed-input rational-flow comparison. The most recent y-linear
genus results, including [the variable-coefficient extension](planar_jc48_sep07_carrier_frontier.md),
do not include Rstar: it has y-degree five. The canonical hostile to a
bare non-local-nilpotence argument remains the rational flow of x²t.
The corrected near miss here would be to quotient by y->-y while fixing
I: this involution sends I to -I and does not act on the generic fibre
I=c. We instead retain the *squared critical values* of the p-projection.
That sidecar determines actual branch points without discarding c.

The live concepts are the actual six-term supplier, its retained
invariant, critical-value separation, ramification at both boundaries,
and scalar-time iteration. Root independently identified the degree-four
Newton polygon with27 interior lattice points; the proof below uses
Riemann--Hurwitz and does not require a Newton-genus theorem.

**Theorem.** The generic curve of (1) is geometrically integral and has

```
g=(9m+21+(m mod2)-gcd(3,m+5))/2.                        (3)
```

In particular the actual supplier has genus27. Every nonconstant
polynomial outer Hamiltonian f(I), at every scalar time lambda!=0,
fails to have both p and y images rational in K(p,y). This includes
f(I)=I and the nonzero times used by the actual finite supplier.

## 2. Exactly two squared critical values

Let c be transcendental over K. Direct differentiation gives

```
I_y=p² y² Q(y²),
Q(z)=3p³A+5(bp³-A)z-7b z².                             (4)
```

At a nonzero finite p on I=c, y is nonzero. Thus every ramification
point comes from one of the two roots z1,z2 of Q. Define

```
Vi=p⁴ zi³(p³-zi)²(A+b zi)²,  i=1,2.
```

These are the squared critical values of I at y=±sqrt(zi).
Since Q has constant nonzero leading coefficient, its roots are integral
over K[p], and the symmetric expressions V1+V2,V1V2 belong to K[p].
For explicit verification put X=bp³ and

```
T=3125A⁶+11875A⁵X+18700A⁴X²+19036A³X³
  +18700A²X⁴+11875AX⁵+3125X⁶.
```

Reduction modulo Q gives

```
V1+V2=4p⁴(X-A)T/(7⁷b⁵),
V1V2=-432A⁵p²³(A+X)⁴/(7⁷b⁵).                          (5)
```

Hence the residual branch polynomial is the monic polynomial in c

```
N=(c²-V1)(c²-V2)=c⁴-(V1+V2)c²+V1V2 in K[p,c].         (6)
```

Its p-degree is23+9m: the product in (5) has that degree and nonzero
leading coefficient -432 a_m^9/(7⁷b⁵), while the trace has degree4+7m.
Also N(0,c)=c⁴. As an independent check of every constant and sign,
literal degree-seven discriminant expansion gives

```
disc_y(I-c)=-7⁷ b⁶ p¹² c² N.                           (7)
```

The four nonzero critical values are genuinely distinct at generic p.
The discriminant of Q and the linear coefficient of Vi reduced modulo Q
are, respectively,

```
DQ=25A²+34AX+25X²,
L=4p⁴ DQ(25A⁴+19A³X-39A²X²+19AX³+25X⁴)/(7⁶b⁴).
```

Both are nonzero polynomials because m>3. The same is true of V1V2.
Moreover V1-V2=L(z1-z2), so its square is a nonzero polynomial in p.
Consequently collisions, zero critical values, and multiple roots of Q
occur only at finitely many p-values algebraic over the coefficient
field. At any such finite value, the monic quartic (6) cannot vanish at
the transcendental c. None is a generic branch value in (6).

Every zero of N is simple. At such a zero exactly one Vi equals c²,
and exactly one of y=±sqrt(zi) has I=c. This is a nontriple critical
point of the y-polynomial. On its local critical branch yi(p),

```
dVi/dp=2I(p,yi) dI(p,yi)/dp=2c I_p(p,yi),
```

because I_y=0 there. A multiple zero of N would therefore give
I_p=I_y=0 at I=c. This is impossible by the inherited finite-critical-
*values* lemma: the algebraic critical locus of a polynomial has finitely
many components, and its restriction is constant on every curve component
in characteristic zero. It thus has finitely many values algebraic over
K, none equal to c. This argument tolerates critical curves and repeated
factors of A. We obtain exactly23+9m simple finite nonzero branch points,
each with one index-two point. All other finite nonzero points are
unramified, since the leading y coefficient -bp² is a unit there.

## 3. Complete boundary count and geometric integrality

At zero, set p=r⁷ and y=r^-2Y. The initial equation of I-c is
-bY⁷-c=0: all terms involving A have higher r-order, as do the lower
y powers. Its seven roots are simple. The deck action on these roots
is transitive, since gcd(7,2)=1. They descend to one point of index seven,
contributing six to ramification. This also proves geometric
integrality: after extending K(c) to its algebraic closure, the complete
degree-seven p-cover has a component already containing a point of degree
seven over the local p-parameter, leaving no other component. The generic
equation is irreducible before constant extension as well, since I-c is
irreducible and linear in the variable c.

At infinity there are three separate groups, since m>3:

* Two roots have y~±p^(3/2), from D. The chart p=r^-2,y=r^-3Y,
  multiplied by r^(19+2m), has initial polynomial
  a_m Y³(1-Y²). The simple nonzero pair is exchanged by the quadratic
  deck action. Its contribution is one.
* Two roots have y~±sqrt(-a_m/b)p^(m/2), from A+b y².
  With p=r^-2,y=r^-mY, multiplication by r^(4+7m) gives initial
  polynomial -Y⁵(a_m+bY²). Its two simple nonzero roots give two
  unramified points for even m, and one index-two point for odd m.
  The contribution is m mod2.
* Three roots are small: y~(c/a_m)^(1/3)p^(-(m+5)/3).
  With p=r^-3,y=r^(m+5)Y the initial polynomial is a_mY³-c.
  It has simple nonzero roots. They descend to three unramified points
  when3 divides m+5, and one index-three point otherwise. The contribution
  is3-gcd(3,m+5).

The groups have respective degrees2,2,3, so these lists exhaust the
degree-seven cover. Simple initial roots justify the Hensel lifts and
exclude hidden local ramification. The total ramification is

```
(23+9m)+6+1+(m mod2)+3-gcd(3,m+5).
```

Applying the characteristic-zero
[Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
to the p-map gives (3). For m=4 the entries are59 finite simple points,
zero contribution6, and infinity contribution1; hence2g-2=-14+66=52
and g=27. The source quartic coefficients in (2) satisfy every hypothesis.

## 4. Scalar-time consequence and a genuine failed extension

Every f(I)-f(0) is divisible by p²D, so the completed source automorphism
and exact scalar-time group law are supplied by
[the audited all-order carrier](planar_jc_long_20260906_hamiltonian.md).
Use the original logarithmic source coordinates p=s²+tau,y=sp. Then

```
I=tau W(s)+O(tau²),
W=s¹⁷(A(s²)+b s⁶) != 0.
```

The highest term of A(s²) has degree2m>=8, so it cannot cancel b s⁶.
If the first nonconstant outer term is f_k I^k, the displacement of p is
2lambda k f_k s W^k tau^k+O(tau^(k+1)). Every positive iterate remains
nontrivial in characteristic zero. The exact fixed-input comparison with
the actual source completion is the one audited in
[the earlier nonrationality theorem](planar_jc_long_20260906_nonrational.md).
It is not an identification of all completed elements with Laurent
rational functions.

If both p and y images were rational at lambda!=0, they would induce
a nonconstant selfmap of the retained smooth proper generic curve over
K(I). Its genus is greater than one, so Riemann--Hurwitz makes it an
automorphism, and its geometric automorphism group is finite. This
contradicts the explicit displacement at every positive time iterate.
The primary [function-field correspondence](https://stacks.math.columbia.edu/tag/0BY1)
and [finite-automorphism input](https://stacks.math.columbia.edu/tag/0DST)
are precisely those already read and audited for the incoming theorem.

The requirement m>=4 cannot be removed by formal extrapolation of (3).
At degree three take A=-bp³. Then

```
I=-b p² y³ D².
```

Two infinity factors have merged. To compute its actual genus, put
w=y²/p³ and h=y/p; inversely p=h²/w and y=h³/w. Its generic curve is

```
h²⁵=(-c/b) w¹¹/(1-w)².
```

The degree25 cyclic cover has precisely the three branch points0,1,infinity
of orders11,-2,-9. Each is totally ramified, since those orders are
coprime to25. Riemann--Hurwitz gives genus12. Extrapolating (3) at m=3
would give24, so this is an exact hostile to that extension, not a
low-genus escape inside the theorem. No conclusion about all other
degree-three specializations is needed.

For the actual completed-response supplier, write its external boundary
parameter as sigma, unrelated to the logarithmic source coordinate s.
The selected time j(sigma) is a scalar in the coefficient field at each
specialization; equivalently one first extends that field to K(sigma).
No source-dependent time is used in the scalar group law. Every nonzero
specialized j(sigma) is excluded as a rational source map. At its zeros
the selected time is zero and the map is the identity; those are
explicitly exempt.
Its finite row-fifteen response remains exactly as proved. This result
does not imply that each separate image is transcendental, does not
classify every lift in the five-dimensional finite response image, and
does not exclude compositions retaining different invariants.

## 5. Exact reproduction

[The standalone source](../../04-computation/planar_jc48_sep07_supplier_genus.py)
imports no research producer. It reconstructs the derivative, the
quadratic quotient trace/norm and critical separation, and the full
degree-seven discriminant. Ten named controls include the exact supplier,
repeated or vanishing A factors, degrees4..10, and nonunit b. A
degree-preserving squarefree finite-field specialization checks the
entire residual branch polynomial for each named row. Independent local
charts cover degrees4..12. The degree-three hostile is checked by its
literal birational cyclic-cover equation. These finite controls support
the preceding all-parameter proof; they do not replace it.

```
python3 -B 04-computation/planar_jc48_sep07_supplier_genus.py
python3 -B -O 04-computation/planar_jc48_sep07_supplier_genus.py
```

Normal and optimized replays pass130 always-active gates with identical
401-byte output. The frozen pins are

```
source a23ac00b69e027efca1e720b9af4ffaf285203551ebd287f05cc28a6c534c83b
output 1d70b64d691d5aedd69246ed3d492fd0fe66060fa4dbaeeba355f439e0680b62
semantic 36b018dd04b607e83171a2086f952b4932fecf5f0c37b0996ea1a89558e2e51c
```

The source and output are frozen. Independent analytic audit is pending.

The [independent analytic/source audit](planar_jc48_sep07_supplier_genus_audit.md) passes, including separate normal and optimized replays. The actual source, output and all-parameter proof are frozen.
