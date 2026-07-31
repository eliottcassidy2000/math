---
id: THM-2811
title: "Special-image linear intertwiner rigidity and JC degree wall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every
  output-compatible affine-linear substitution intertwining the
  special-image contractions has matrices L,S satisfying L S^T=I.  Hence it
  can only embed, never lower, the number of variable pairs and it preserves
  the separate derivative/coordinate top degrees.  The two THM-2801
  witnesses therefore cannot enter Zhao's xi-linear cubic or separable
  quartic sectors through any such substitution.  Independently, every
  admissible two-pole e=2 response has no totally ramified fibre, its
  potential has degree 2N and at least two projective roots.  Across the
  complete e=2 response layer, h=1 and h=3 are polynomializable while h=2
  is not, but both polynomializable chambers have multiple projective roots
  and hence miss the binary Hessian-nilpotent pure-power locus.  The two
  degree-four boundary polynomials fail Hessian nilpotence explicitly.
  Response and SIC cancellations do share an exact barycentric-stencil
  mechanism, but their residue alphabets and, crucially, their
  fixed-versus-growing grading are different.
source: root/jc-sic-intertwiner-2026-07-28
audit: thm2811-independent-hostile-audit-2026-07-28
depends_on:
  - THM-2801-sharp-special-image-boundary-and-beta-shift-witness
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2799-one-pole-two-double-zero-chebyshev-response-classification
  - THM-2800-two-pole-two-double-zero-stieltjes-recurrence-and-first-nielsen-pair
  - THM-2805-general-two-pole-e2-maxwell-eliminant-and-nielsen-classification
  - THM-2808-three-pole-e2-maxwell-polynomial-and-finite-accessory-classification
related:
  - MISTAKE-237
script: 04-computation/sic_jc_intertwiner_response_stencil_thm2811.py
output: 05-knowledge/results/sic_jc_intertwiner_response_stencil_thm2811.out
script_sha256: 534fc539da0efeccd446f4f22ab84b425b4dc771c6012d76128b6d7b23f37ee0
output_sha256: e8d61aab5ca40bc270afcf1e2ad818a917a2405a1627d31d8fd3adb79c3f5270
independent_script: 04-computation/sic_jc_intertwiner_response_stencil_independent_audit_thm2811.py
independent_output: 05-knowledge/results/sic_jc_intertwiner_response_stencil_independent_audit_thm2811.out
independent_script_sha256: 8c6c586694d0733caab5acd789beac699d3e9cef78e75f6fbf233fae079f41a1
independent_output_sha256: 96b2e568583f5761f7379c0e226a0f047d6bf3c202371bd6c35d1a6d6b828fc7
hash_basis: LF-normalized bytes
---

# THM-2811 -- special-image linear intertwiner rigidity and JC degree wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

There is a real common mechanism behind the newest special-image and planar
response calculations: both use a maximal-order Lagrange annihilator.  The
common mechanism does not supply a transport.  The operator-compatible
linear maps are rigid, the two-pole response maps are not polynomializable,
and the one- and three-pole chambers that are polynomializable miss
Hessian nilpotence.

This theorem makes all three statements exact.

## 1. Affine-linear intertwiners of the special-image operator

Put

```text
A_n=C[xi_1,...,xi_n,z_1,...,z_n],
E_n(P(xi)Q(z))=P(partial_z)Q(z).
```

Let

```text
phi_0:C[z_1,...,z_n] -> C[y_1,...,y_k]
```

be an affine-linear substitution

```text
z=S y+a,                    S in Mat_(n x k)(C).       (1)
```

Consider an affine-linear algebra homomorphism

```text
Phi:A_n -> A_k=C[eta_1,...,eta_k,y_1,...,y_k]
```

which extends `phi_0` and is output-compatible, so initially

```text
xi=K y+L eta+c,
K,L in Mat_(n x k)(C),       c in C^n.                (2)
```

Then

```text
E_k Phi=phi_0 E_n                                      (3)
```

on every polynomial if and only if

```text
K=0,                    c=0,                    L S^T=I_n. (4)
```

### Proof

Apply `(3)` first to `xi_i`.  Its right side is zero, while its left side is

```text
sum_alpha K_(i,alpha)y_alpha+c_i.
```

Thus `K=0` and `c=0`.  Apply `(3)` to `xi_i z_j`.  The right side is
`delta_(ij)`, and the left side is

```text
sum_alpha L_(i,alpha)S_(j,alpha)=(L S^T)_(ij).
```

This proves necessity.

Conversely, by linearity it is enough to take `P(xi)Q(z)`.  The chain rule
and `(4)` give

```text
E_k Phi(PQ)
 =P(L partial_y) Q(Sy+a)
 =P(L S^T partial_z)Q(z) evaluated at z=Sy+a
 =phi_0(E_n(PQ)).
```

This proves sufficiency and the exact classification.

### Consequences

Equation `(4)` forces

```text
rank(S)=rank(L)=n,                    hence k>=n.       (5)
```

Both the `z` and `xi` affine-linear forms are algebraically independent, so
`Phi` is injective.  In particular:

1. an exact affine-linear intertwiner cannot reduce the number of pairs;
2. the `xi` degree of every nonzero bihomogeneous polynomial is preserved;
3. the top `z` degree is preserved, even when `a!=0`; and
4. when `a=0`, the full bidegree is preserved.

The familiar cotangent gauge is the square case

```text
k=n,             z=S y+a,             xi=S^(-T)eta.   (6)
```

For `k>n`, `(4)` describes an injective skew padding.  Pair deletion is not
an operator intertwiner; this explains algebraically why setting a pair to
zero need not preserve special-image kernel membership.

## 2. The exact degree wall for the THM-2801 witnesses

The two counterexample bases of THM-2801 have bidegrees

```text
rank-one two-pair base F:       (4,4),
four-term three-pair base f_1:  (2,2).                 (7)
```

Zhao's two JC-bearing bases have the forms

```text
xi-linear cubic:       sum_i eta_i H_i(y),             (1,3),

separable Hessian:     (sum_i eta_i^2)P(y),             (2,4),
                       P homogeneous quartic.           (8)
```

By Section 1, no output-compatible affine-linear intertwiner sends either
nonzero polynomial in `(7)` to either sector in `(8)`.  Translation may add
lower `y`-degree terms, but it cannot remove the nonzero top component.
Injectivity excludes collapse to zero.

This strengthens the linear-gauge boundary in THM-2801 in two ways: it
includes affine coordinate changes and it covers the bidegree-minimal
three-pair witness as well as the two-pair witness.  The coordinate observer
still has the right superficial type: a `z_i` becomes an affine-linear
coordinate observable.  The obstruction is the base polynomial.

An arbitrary linear change mixing derivative and coordinate variables is not
covered by `(3)` and need not preserve `ker(E)`.  An arbitrary
dehomogenization, polarization, or differential projection is likewise not
multiplicative on powers.  None is a transport of the SIC hypothesis unless
a separate intertwining theorem is supplied.

## 3. Two-pole responses are intrinsically nonpolynomial

The next obstruction is independent of Section 1.  Let `N>=4`, let
`d+b=N` with `d,b>=1`, and consider any admissible balanced
`e=2,h=2` two-pole response datum in the normalized form

```text
D=x^d(x-1)^b,
E=(x-gamma)(x-delta),
B=S E^2,
deg(S)=N-4,
B-D=A(x-z_0),                    A!=0,                 (9)
```

with the squarefree/disjointness gates of THM-2796.  Put

```text
R=B/D.
```

It is a degree-`N` rational map.  Its three distinguished fibre partitions
are

```text
R^(-1)(0):          2^2 1^(N-4),
R^(-1)(infinity):   (d,b),
R^(-1)(1):          (N-1,1).                          (10)
```

For the last line, `B-D` has one finite simple root, while
`R-1=A(x-z_0)/D` has a zero of order `N-1` at infinity.  The ramification
defects in `(10)` sum to

```text
2+(N-2)+(N-2)=2N-2.                                  (11)
```

Thus Riemann--Hurwitz leaves no further branch value.

A degree-`N` rational map is source/target-Mobius equivalent to a polynomial
only if it has a totally ramified fibre `(N)`: send that fibre and its unique
point to the two infinities.  Conversely that condition plainly suffices.
No partition in `(10)` is `(N)`.  Therefore

```text
no admissible two-pole response is Mobius-equivalent to a polynomial. (12)
```

This applies to the proved edge atlas `d=1` of THM-2800 and to every proved
middle-partition object of THM-2805.

Two immediate JC walls follow.

1. The projectivization `[H_1:H_2]` of a planar homogeneous cubic vector
   field has rational degree at most three after its common factor is
   removed.  It cannot equal the degree-`N>=4` response map after Mobius
   changes.
2. A dehomogenized binary form is polynomial in some projective chart.
   Equation `(12)` forbids identifying the response map itself with such a
   form, including a binary quartic.

These are no-go statements for the direct identifications.  The response is
a downstream invariant of a special planar source-fibre chart; they do not
exclude a nonlinear reconstruction of a Keller pair followed by a
dimension-raising symmetric reduction.

### The complete `e=2` polynomializability trichotomy

THM-2796 forces `h<=e+1=3`, and THM-2799/2800/2805/2808 now classify all
three chambers.  Their total-fibre ledger is:

```text
h=1:  pole fibre (N)                         polynomializable;
h=2:  fibres 2^2 1^(N-4), (d,N-d), (N-1,1) not polynomializable;
h=3:  third fibre (N)                        polynomializable. (13)
```

For `h=1`, send the unique pole to infinity; in THM-2799's normalization
this is the source inversion of `B/x^N`.  For `h=3`, THM-2808 has

```text
F=(D-v)/D,
1/(1-F)=D/v,                                           (14)
```

so the displayed target Mobius transformation is already a polynomial.
Section 3 proves the middle line.  Thus polynomializability is not the
missing `e=2` accessory coordinate: it alternates `yes/no/yes` with the
number of poles.

There are no hidden total fibres.  For `h=1`, the defects of

```text
2^2 1^(N-4),                    (N),                    (N-2,1,1)
```

are `2,N-1,N-3`, summing to `2N-2`.  For `h=3`, the defects of

```text
2^2 1^(N-4),                    (a,b,c),                (N)
```

are `2,N-3,N-1`, again summing to `2N-2`.  Riemann--Hurwitz therefore
exhausts the ramification in both polynomializable chambers, just as
`(11)` did for two poles.

## 4. The response potential misses the Hessian sector too

The two-pole response potential is

```text
V_resp=c S D [x(x-1)]^2
      =c S x^(d+2)(x-1)^(b+2),              c!=0.     (15)
```

It has degree

```text
deg(V_resp)=2N>=8                                   (16)
```

and has at least the two distinct projective roots `0` and `1`.  Degree
alone excludes Zhao's quartic `P`.  The root pattern gives a stronger
obstruction.

### Binary Hessian-nilpotent lemma

Let `P(X,Y)` be a nonzero binary homogeneous form of degree `q>=2`.  Its
Hessian matrix is nilpotent if and only if

```text
P=c(X+iY)^q              or              P=c(X-iY)^q (17)
```

after a complex orthogonal change of coordinates.

Indeed, nilpotence forces both trace and determinant to vanish.  With

```text
u=X+iY,                    v=X-iY,
Delta=4 partial_u partial_v,
```

the trace condition makes a homogeneous `P` harmonic, hence

```text
P=A u^q+B v^q.
```

Its Hessian determinant is

```text
-4q^2(q-1)^2 A B u^(q-2)v^(q-2).                    (18)
```

Thus `AB=0`.  Conversely the Hessian of a pure power of an isotropic linear
form has square zero.  In particular every nonzero binary
Hessian-nilpotent form has one projective root, with full multiplicity.

The homogenization of the two-pole potential has at least two roots, so it
is not Hessian-nilpotent in any degree.  In either polynomializable
chamber, the unique total fibre from Section 3 must be sent to polynomial
infinity.  Every finite zero fibre after an arbitrary target-affine
normalization is therefore one of the listed non-total branch fibres or an
unramified fibre.  It can never be `(N)`, so the resulting binary form can
never be a pure `N`th power.  In the displayed representatives this
obstruction is visible as

```text
h=1: the polynomialized response has the two distinct double-zero roots;
h=3: D/v has the three distinct pole roots 0,1,lambda. (19)
```

Hence no member of the complete abstract `e=2` response layer gives Zhao's
binary Hessian-nilpotent form by direct source/target-Mobius
polynomialization and binary homogenization.  Taking the squareclass,
squarefree radical, numerator, or a derivative is not a coordinate change
and does not preserve the response equation or the Hessian predicate.

### Both polynomializable degree-four boundaries are hostiles

An invertible source change cannot coalesce the two poles.  Even if one
leaves the two-pole chamber and passes to the adjacent one-pole partition
`(4)`, the unique degree-four response from THM-2799 is

```text
R(x)=(x^2-1)^2/x^4.
```

After the source inversion `x=1/y`, it becomes the polynomial

```text
p(y)=(1-y^2)^2,
P_1(X,Y)=(X^2-Y^2)^2.                                (20)
```

Exact differentiation gives

```text
Delta P_1=8(X^2+Y^2),
det(Hess P_1)=-48(X^2-Y^2)^2.                        (21)
```

Thus this closest quartic boundary is not Hessian-nilpotent.  This is not an
artifact of the chosen polynomial normalization.  The response passport is

```text
(2,2),                    (4),                    (2,1,1).
```

In every polynomial presentation the unique `(4)` fibre must be the pole
at infinity.  Every finite zero fibre is therefore `(2,2)`, `(2,1,1)`, or
unramified, never `(4)`.  The binary lemma rules out every such polynomial
representative, because a Hessian-nilpotent quartic is a fourth power and
has a `(4)` zero fibre.

More directly, for the displayed normalization and the separable
special-image base

```text
a=(eta_1^2+eta_2^2)P_1,
```

one already has

```text
E_2(a)=Delta P_1!=0.                                 (22)
```

The three-pole degree-four chamber has, up to reordering,

```text
(a,b,c)=(1,1,2),                    lambda=1/2.
```

Equation `(14)` homogenizes, up to a nonzero scalar, to

```text
P_3(X,Y)=X(X-Y)(2X-Y)^2.                              (23)
```

It has three projective roots, as the general argument predicts, and the
exact first tests are

```text
Delta P_3=2(29X^2-27XY+5Y^2),
det(Hess P_3)
 =-3(2X-Y)^2(8X^2-8XY+3Y^2).                        (24)
```

Thus the associated separable base also has first contraction
`Delta P_3!=0`.  Both polynomializable `N=4` endpoints fail before any high
power is examined.  The `h=2,N=4` response remains intrinsically rational
by Section 3.

The other tempting singular operation also leaves the chamber.  Coalescing
the two double zeros lands on a removed diagonal inflection.  There
`B=B'=B''=0`, but the non-pole root of `D''` is simple, so `B'''!=0`:
the contact is exactly triple, not the quadruple contact required by
`E=(x-u)^2` and `E^2|B`.  For the middle partitions the discriminant of the
inflection quadratic is

```text
16 d b (N-1)!=0,
```

and for the edge partition the non-pole factor of `D''` is linear.

## 5. The exact bridge: both cancellations are barycentric stencils

The negative results above do not mean the resemblance is cosmetic.  It
has one exact common core.

For any balanced response in the normal form of THM-2796, put

```text
M=S E T,                    deg(M)=r+1.
```

THM-2796 gives

```text
F_resp'/F_resp=C/M,             C!=0.                 (25)
```

At a root `a` of `M`, let `w_a=ord_a(F_resp)`.  Taking residues in `(25)`
gives

```text
w_a=C/M'(a),                                         (26)
```

with the exact residue alphabet

```text
{w_a}={1^s,2^e,-p_1,...,-p_h}.                       (27)
```

Lagrange interpolation now gives, for every polynomial `H`,

```text
sum_(M(a)=0) w_a H(a)
 =C [x^r](H mod M).                                  (28)
```

Consequently the response moments vanish through degree `r-1` and the next
one is forced:

```text
sum_a w_a a^j=0          (0<=j<=r-1),
sum_a w_a a^r=C.                                      (29)
```

This is the exact mechanism behind the Padé moment prefix, not another
analogy.  In the two-pole `e=2` chamber, `r=N-1`,
`M=SEx(x-1)` has degree `N`, and `(27)` specializes to

```text
{1^(N-4),2,2,-d,-b}.
```

Now take the four-term SIC base `f_1` of THM-2801.  At power `m`, its pure
contraction is a sum over the channel count `j=0,...,m` with weights

```text
lambda_(m,j)=(-1)^(m-j) binom(m,j).                   (30)
```

If

```text
M_m(X)=product_(j=0)^m(X-j),
```

then

```text
lambda_(m,j)=m!/M_m'(j).                              (31)
```

Thus the binomial cancellation is the equispaced instance of the same
Lagrange identity.  Multiplication by the coordinate observer shifts and
truncates this stencil; the missing endpoint leaves the nonzero coefficient
recorded in THM-2801.

The typing boundary is exact:

```text
response: one fixed (r+1)-node geometric stencil, with next defect C!=0;
SIC:      a coherent new (m+1)-channel stencil generated by every power m.
                                                                  (32)
```

No raw affine relabelling identifies the two weight systems.  For `N>=6`,
the response signs split as `(N-2,2)`, whereas the binomial signs are as
balanced as possible.  At `N=4,5`, the respective response multisets

```text
{2,2,-3,-1},
{1,2,2,-4,-1} or {1,2,2,-3,-2}
```

are not proportional to the binomial multisets for `m=3,4`.

Therefore a transfer would need, at minimum:

1. an inter-degree tower of response node sets;
2. a multiplicative law realizing that tower as powers of one base;
3. marked-node observers compatible with the Faber flux equations; and
4. actual Keller-chart entry.

THM-2799/2800/2805/2808 provide fixed-degree response/Nielsen atlases, not
these sidecars.  The useful positive lesson is computational: any proposed
flux functional on a response should first be reduced modulo `M`; equation
`(28)` shows that only its `x^r` remainder coefficient can survive
(`r=N-1` in the two-pole chamber).

## 6. Exact audits and scope

The primary lightweight companion checks:

1. the mixed contraction matrix `L S^T` on the generators;
2. the two-pole fibre defects/potential degree and the complete
   `h=1,2,3` polynomializability trichotomy through `N=40`;
3. the exact `N=4` response residues and barycentric moments;
4. deterministic rational-node Lagrange identities;
5. the equispaced derivative formula `(31)` through `m=40`;
6. the response/binomial weight mismatch through `N=40`; and
7. both polynomializable quartic Laplacian/Hessian hostiles `(21),(24)`.

Run

```bash
python 04-computation/sic_jc_intertwiner_response_stencil_thm2811.py
python -O 04-computation/sic_jc_intertwiner_response_stencil_thm2811.py
```

The independent companion does not import the primary.  It uses a direct
monomial implementation of `E`, a separate rational-polynomial remainder
engine, a different signed node family for `(28)`, and scalar symbolic
Hessian calculations.  It checks the affine identity on a full small
monomial box, all three passport ledgers, both quartic hostiles, and the
stencil boundary in `20,163` exact gates.  Run

```bash
python 04-computation/sic_jc_intertwiner_response_stencil_independent_audit_thm2811.py
python -O 04-computation/sic_jc_intertwiner_response_stencil_independent_audit_thm2811.py
```

Normal, optimized, and stored transcripts agree for both companions, and
neither script contains an `assert` AST node.

The finite controls support but do not replace the all-degree chain-rule,
Riemann--Hurwitz, binary-harmonic, and Lagrange proofs.

This theorem proves no restriction on an unspecified nonlinear transform.
It does not prove that an abstract response enters a Keller chart, prove a
restricted image conjecture, or settle `JC(2)`/`DC(2)`.
