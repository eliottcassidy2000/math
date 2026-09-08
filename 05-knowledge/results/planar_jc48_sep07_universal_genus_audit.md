# Independent audit: the entire universal carrier has no nonzero rational time

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS.**
The complete theorem in [the primary proof](planar_jc48_sep07_universal_genus.md)
is accepted. Every geometric generic component of every nonzero polynomial
divisible by p(p³-y²) has genus at least two. For every nonconstant
Hamiltonian in K+p²(p³-y²)K[p,y], this excludes simultaneous rationality
of its two source images at any nonzero scalar time. The proof is valid
over every characteristic-zero coefficient field, with the finite-field-
of-definition paragraph added at the parent's request. No primitive-fibre
assumption or generic-coefficient restriction remains.

This does not exclude compositions retaining different invariants,
assert polynomial termination of a completed map, solve later finite
Keller equations, or solve JC(2).

There is a precise incoming antecedent: the final scope discussion in
[continuing9 flow-carrier genus](continuing9_20260907_flow_carrier_genus.md)
explicitly left the local-cusp high-genus claim unproved for arbitrary
multiplicities and extra factors. The present proof pays that stated
local/global and multiplicity obligation. The underlying research question
is inherited, rather than newly proposed in this bundle.

## 1. The central cover is an actual embedded surface

I independently rederived the weighted factorization. The least
(2,3)-weight form of I=p(p³-y²)R contains both mandatory factors.
Over C it is a nonzero constant times

```
p^A y^B product_i(p³-lambda_i y²)^e_i,
A>=1, B>=0, e_i>=1,
```

with distinct nonzero lambda_i, one equal to1. Removing common p,y
powers reduces any weighted homogeneous polynomial to a binary
homogeneous polynomial in p³,y², which proves the asserted complete
factorization. Coincident binomials are combined into their multiplicity;
they are not counted as distinct branch points.

The substitution p=v²w,y=v³w has inverse v=y/p,w=p³/y² off the two
coordinate divisors. It gives exactly

```
I=v^M(Phi(w)+v Psi(v,w)),
M=2A+3B+6E, Phi=C w^(A+B+2E) product_i(w-lambda_i)^e_i,
E=sum e_i.
```

Higher weight is strictly higher v-order, so Psi is polynomial here.
No completed embedded resolution of all factors of R is being assumed.

The proposed embedded-piece argument is stronger than a heuristic
identification of a Milnor fibre with its initial form. On a neighborhood
of the compact punctured sphere K, Phi is nonzero. Choose one uniform
v-disk where |v Psi/Phi| is small. The convergent binomial root near1
is single-valued throughout this neighborhood, even though Phi itself
need not have a single-valued M-th root. The coordinate

```
vtilde=v(1+v Psi/Phi)^(1/M)
```

satisfies I=Phi(w) vtilde^M exactly. Its derivative is uniformly close
to1. Shrinking the disk gives injectivity and a uniform inverse disk
in each fibre, by the derivative estimate and Rouché/implicit inversion.
For every sufficiently small c!=0, all M roots of c/Phi lie in that
inverse disk. Thus the entire actual cover over K is identified with
vtilde^M Phi(w)=c, not merely a selection of its sheets.

The piece has v,w nonzero; blowdown is injective there with the displayed
rational inverse. It is therefore a compact bordered surface embedded
in the actual affine fibre I=c. Capping its boundary later computes its
genus only: no capped compact complex curve is claimed to embed in an
affine surface or in this fibre.

## 2. Component gcd, primitive normalization, and the sharp inequality

The actual deck increments are -n0,-e_i,ninf, with
n0=A+B+2E and ninf=A+B+3E. Their sum is zero. The subgroup they generate
in Z/M has index

```
d=gcd(M,n0,e_1,...,e_r)=gcd(A,B,e_1,...,e_r).
```

The omitted ninf adds no condition because ninf=n0+E. Conversely,
M-2n0=B+2E recovers divisibility of B, then n0-B-2E recovers A.
This checks both directions of the component count. There are d
components, each of degree M/d; capping their lifted boundary circles
gives the usual cyclic covers with precisely these branch exponents.

The total ramification at an exponent n is M-gcd(M,n). Summing all
r+2 punctures and applying Riemann--Hurwitz gives exactly

```
2d(g-1)=rM-gcd(M,n0)-gcd(M,ninf)-sum_i gcd(M,e_i).
```

All components have the same genus. Dividing A,B,e_i by d divides the
right side by d, rather than leaving it unchanged, and preserves g.
The primary proof and source normalize this correctly. The normalized
A remains positive and the normalized E is positive.

If B>0, the bounds on the three subtractions are B+2E, B, E;
the right side is at least rM-2B-3E>=2A+B+3E>=6.
If B=0, ninf=M/2 exactly; the remaining bounds give
(2r-1)A+6(r-1)E>=1. This integer is bounded above by the actual
even integer2(g-1), so that actual integer is at least two. Thus g>=2
in every case. In particular the proof requires neither A>=2 nor
squarefreeness of the initial or full polynomial.

The equality and failure boundaries can also be checked on actual
generic curves. For pD=c, putting v=py gives v²=p⁵-cp, a squarefree
degree-five hyperelliptic polynomial. For p²D=c it gives v²=p⁵-c.
Both have genus two. The power (pD)^n has n geometric generic components,
each of genus two. Conversely D=c gives the squarefree cubic
y²=p³-c, with three finite double-cover branch points and the point
at infinity, hence genus one. The fibre p=c has free coordinate y
and genus zero. These verify actual global sharpness and the two
missing-factor hostiles, beyond the source's arithmetic sanity gates.

## 3. Generic components and field descent are paid

Resolve the rational polynomial pencil on P² to a morphism from a
smooth projective surface and take its finite Stein base. The base is
a connected curve, and after deleting finitely many points its connected
fibres are smooth proper curves of constant genus. Vertical exceptional
curves occur only at finitely many values; horizontal exceptional curves
remove only finitely many points on a general fibre. Thus these smooth
proper fibres compactify the corresponding components of the original
affine pencil. One can choose c both sufficiently small for the embedded
cover and outside the finitely many excluded values.

Each embedded bordered genus-g surface supplies g handle pairs with a
nondegenerate intersection form in its containing compactification.
Their intersection numbers are unchanged by inclusion, so the ambient
genus is at least g. At least one component of a general fibre therefore
has genus at least two. All geometric generic components have the same
genus: I-c is irreducible before extending constants, and the embeddings
of the finite relative algebraic constant extension permute these
components transitively. Equivalently their genus is the constant genus
on the connected Stein base's smooth open. This proves the required
statement on *every* geometric generic component.

The added arbitrary-field paragraph is valid. The coefficients of a
fixed polynomial descend to a field finitely generated over Q, which
embeds into C. Genus of the geometric components is unchanged by passing
between algebraically closed extensions of that field. For the rational-
time argument include the finitely many coefficients of the proposed
rational images and lambda in the same field. Clearing denominators
reduces the asserted equality to coefficientwise formal identities;
the derivation's iterates use only field operations and factorial
denominators. A field embedding preserves those identities and the
nonzero time. It is unnecessary, and in general impossible, to embed
an entire arbitrary characteristic-zero coefficient field into C.

## 4. The completed-flow consumer also handles composite invariants

For S=c0+I in the universal carrier, write I=D^e T with e>=1 and D not
dividing T. In the actual logarithmic chart,

```
I=tau^e W(s)+O(tau^(e+1)),
W=s^(4e)T(s²,s³) != 0.
```

The kernel of the cusp parametrization is exactly (D), proving this
nonvanishing. With the original tau-Poisson bracket the first displacement
of p is2lambda e sW tau^e. The derivation raises tau order by at least e,
so higher iterates cannot cancel it. This remains nonzero after replacing
lambda by any positive integer multiple. This check uses the actual
source p=s²+tau, not an abstract coordinate on a quotient curve.

If both images were rational, the faithful fixed-input comparison and
formal invertibility give an injective field endomorphism of K(p,y)
fixing K(I). Let E0 be the finite relative algebraic closure of K(I).
Its image is contained in E0 and, by finite-dimensionality, is an
automorphism of E0 over K(I). This finite group has finite order, so
a positive power fixes E0 pointwise. On the smooth proper curve with
function field K(p,y)/E0, the resulting nonconstant selfmap has genus
at least two. Riemann--Hurwitz forces degree one; finiteness of its
geometric automorphism group then makes a further positive power the
identity. This contradicts the scalar-time group law and the explicit
nonzero displacement. An explicit polynomial primitive is not needed.

The argument uses the existing audited completed universal-carrier
theorem and the fixed-input comparison. It does not claim that the entire
completion is rational. Scalar time may be an external coefficient-field
parameter or its specialization; it is not source-dependent time.
The field inverse x=yp/D,t=D/p² also correctly transfers simultaneous
rationality between p,y and the original x,t coordinates.

The conclusion covers every nonconstant polynomial lift in the full
five-dimensional finite response image, including all six terminal
carrier directions and arbitrary higher polynomial tails. It concerns
each single invariant-preserving flow. Inverse-time cancellation and
compositions with different invariants remain outside this exclusion.
The earlier positive finite-row supplier remains valid.

## 5. Independent source read, full replay, and pins

I read the entire standalone source and replayed it normally and under
python -O. Both outputs are byte-identical to the frozen395-byte report.
They pass105404 always-active gates:13104 declared multiplicity tuples,
180 separate literal sheet-permutation cover cases, actual polynomial
chart substitutions, sharp powers and boundary sanity checks.

The literal-cover engine globally reverses the signs of the actual
increments; this reverses each puncture permutation and preserves its
orbits, ramification, component count and genus. The analytic proof keeps
the actual signs. No marked monodromy inference depends on that reversal.
The arithmetic universe tests the all-size formulas; the embedding and
constant-field arguments are supplied by the preceding independent proof
review, not extrapolated from that universe.

```
python3 -B 04-computation/planar_jc48_sep07_universal_genus.py
python3 -B -O 04-computation/planar_jc48_sep07_universal_genus.py
```

```
source ffed4b59e5a5832e43aa74fd325e8589c6973cc008066276be702f8f0b00f872
output 910002d880a2849c8509c00cddb577a87155e719c48a88c66a5b79747abe31d6
semantic 629e17482bfce870f7145297a03c205fe02d8b8f7661c8b8bdf2a9187a627d3a
```

Source size3864 bytes; output size395 bytes. No source changes were
made in this audit. The authorized field-descent paragraph is the only
primary-proof edit made by the auditor. This audit is frozen for parent
promotion and integration.
