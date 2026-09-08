# Independent audit: genus and extremal rational mates on every W_m

**Status: PROVED / INDEPENDENT ANALYTIC, SOURCE AND FULL REPLAY AUDIT PASS.**
This audit accepts the complete all-m filtration, the sharp geometric generic
genus bound for rational-mated quadratic sections, the full extremal iff for
m>=2, the actual sharp examples and their field degrees, and the recovered
no-polynomial-mate corollary in the declared genuine quadratic class. It does
not classify all rational mates below maximal genus, all functions on W_m,
or arbitrary plane quartics, and it makes no planar Jacobian claim.

Auditor: `three_ray_geometry`, independently of `certificate_audit`, the producer.
The [dg_genus primary](planar_jc48_sep08_dg_genus.md) and its complete standalone
source were read in full. Both independent normal and optimized replays pass
**1,550 always-active gates** and are byte-identical to the frozen 856-byte
output. I independently reconstructed the shifted all-m mate bracket and the
finite universe counts. No mathematical, source or prose correction was
required. Producer files were not edited by this auditor.

## 1. Actual surface charts and complete section space

The surface is exactly

    W_m=(P1_x x P1_z) minus {z=x^m},      m>=1.

Over finite x, the coordinate `t=1/(z-x^m)` parametrizes the entire fibre
P1 minus its graph point; z=infinity corresponds to t=0 and is included.
Over the other x chart, `r=1/x` and `b=z/(1-r^m z)` similarly give the
entire punctured fibre, including the projective z=infinity point when it
belongs to W_m. The formulas

    x=1/r, t=-r^m-r^(2m)b,
    b=-x^m-x^(2m)t, z=b/(1+r^m b)

are mutually consistent on the overlap. The denominator-zero locus in the
last expression is a projective-coordinate value, not an omitted locus in
the affine (r,b) chart. At r=0, the excluded graph point is z=infinity,
so D=W_m minus U0 is exactly the affine line r=0 with coordinate b.

The determinant of the actual source transition is `r^(2m-2)`. Therefore
`dx wedge dt=r^(2m-2)dr wedge db`. The proof correctly includes m=1,
where this form has no zero along D. It does not copy the order-two vanishing
of the earlier m=2 surface into every m.

For fixed n>=0 any source polynomial of t-degree at most n has a unique
numerator representation `R(x,z)/(z-x^m)^n`, with z-degree at most n.
At the generic point of D this becomes

    r^(mn) R(1/r,z)/(r^m z-1)^n.

The denominator is a unit. The top x-coefficient of R is a nonzero polynomial
in z, so it cannot vanish at the generic point of D. Hence regularity forces
`deg_x R<=mn`. Conversely each monomial of this whole degree box yields

    E_(i,j)=x^i t^(n-j)(1+x^m t)^j,
    0<=i<=mn, 0<=j<=n,

whose second-chart expression is
`(-1)^n r^(mn-i)b^j(1+r^m b)^(n-j)` and is everywhere regular.
Independence is the independence of numerator monomials with one fixed
nonzero denominator. This proves completeness, including lower t-degree
members in the same space and the n=0 constants. Every global function
restricts to a polynomial on the source A2, so the union of these spaces
is the full global ring.

For n=2 the three numerator coefficients A,B,C each have degree at most
2m and vary independently. Direct expansion gives

    N=A*x^(2m)+B*x^m+C,
    P=2A*x^m+B, Q=A,
    P^2-4NQ=B^2-4AC.

Thus both the genuine leading coefficient N and the constant discriminant
D0 have degree at most 4m. This retains every global coefficient; it is not a
special sparse chart or a leading-degree truncation. The hypothesis N!=0
is explicitly separate and essential.

## 2. Generic components, exactness and the genus bound

The field identity for the fibre H=lambda is

    V=2Nt+P, V^2=D0+4lambda*N, t=(V-P)/(2N).

The possible vertical-line components over roots of N cause only finitely
many fibre values: such a component also requires P=0 and fixes lambda=Q
at that x-coordinate. Thus the generic radical presentation accounts for
all geometric generic components. No entire boundary branch is discarded
merely because this inverse formula divides by N.

With lambda held fixed, a rational Jacobian mate restricts to
`dG=-dx/V`. A fixed denominator can vanish identically on fibre components
at only finitely many values, since it has finitely many irreducible factors
and any factor lying in a fibre determines that one fibre value. The argument
therefore applies to each generic component on which the rational function
is defined. A square discriminant gives two rational geometric components;
they are not counted as a connected double cover. A nonsquare discriminant
gives the actual connected normalized radical curve.

The inherited [genus_capacity theorem](planar_jc48_sep08_genus_capacity.md)
was independently audited before this bundle. I checked its use and the
included proof. Multiplicity two gives nonzero logarithmic residues. A high
root of multiplicity j contributes at most j-2 to the total primitive pole
degree, with the two unramified points counted separately for even j. For odd
polynomial degree n>=3, the unique infinity point forces primitive local
degree n-2, so only a pure power survives and its genus is zero. For even
n>=4 either infinity point forces local degree n/2-1; their two degrees are
not added. With s simple roots and h_o,h_e odd/even high roots,

    B=s+h_o=2g+2,
    P_fin=n-s-2h_o-2h_e>=n/2-1,
    B+h_o+2h_e<=n/2+1.

Some high root is required, yielding `g<=floor((n-4)/4)`. Small degrees
and square fields are handled separately. Since the actual discriminant
degree is at most 4m, every geometric generic component has genus at most
m-1. This includes m=1 and makes no connectedness assumption on a split
geometric generic fibre.

## 3. Equality fixes a root across the whole pencil

At m>=2 the maximal genus m-1 is positive. A square or odd-degree generic
discriminant cannot have that genus, and any smaller even degree has genus
at most m-2. Thus equality requires degree exactly 4m. Equality in the
capacity arithmetic forces one high root of odd multiplicity 2m+1 and
2m-1 simple roots.

The inherited extremal coefficient theorem, also rederived in the primary,
uses `u=x-p`, `X=1/u`, and `Y=V/u^(2m)`. The model has a squarefree
polynomial of degree 2m-1 and a unique infinity point. The differential's
only pole has order 2m, so every primitive lies in the complete pole space

    span{1,X,...,X^(m-1),Y}.

This follows from the normal affine coordinate ring and the distinct even
and odd pole orders. Odd averaging leaves only a scalar multiple of Y.
Its derivative forces the entire residual polynomial to have only constant
and top coefficients. Therefore each generic discriminant has the full form

    (x-p)^(2m+1) [a*(x-p)^(2m-1)+b],  ab!=0.

This is not a multiplicity-only inference, and no unproved period vanishing
is used.

I checked the separate proof that p is fixed, which is necessary for the
claimed pencil and mate formula. If N,D0 are independent, factor their gcd
G0 and write the coprime residuals N1,D1. A repeated root of
`D1+4lambda*N1` has N1 nonzero there, is a zero of the fixed nonzero
Wronskian `D1'*N1-D1*N1'`, and determines at most one lambda. Generic
residuals also avoid the finite roots of G0. Hence every generic repeated
root is a fixed gcd root with fixed multiplicity. The unique high root is
one fixed p. When N,D0 are proportional, the pencil is a scalar multiple
of N and the same statement holds immediately, including D0=0. The proof
does not divide by a vanishing Wronskian in that branch.

After p is fixed, two distinct generic fibre values imply both N and D0
belong to `span{(x-p)^(4m),(x-p)^(2m+1)}`. This proves the whole-pencil
condition, with affine coefficient functions a(lambda),b(lambda), neither
identically zero. Generic multiplicities, residual separation and degree
are stable outside a finite exceptional set. There is no assumption that
separately chosen fibre roots or primitive constants descend automatically.

## 4. Rational parameter descent and the shifted source identity

Conversely, the stated whole pencil has a connected generic radical curve
of genus m-1. The explicit primitive for the relative form is

    G=2V/[(2m-1)b(H)(x-p)^(2m)].

The coefficients are rational in the original parameter and radical variable;
no square root of a parameter is adjoined. The denominator is a nonzero
rational function because b is a nonzero affine polynomial, H is nonconstant,
and x-p is not the zero field element. Poles on exceptional fibres are
allowed. Derivatives of b(H) contribute nothing in a tangent direction of H,
so substitution of the fibre parameter by H is legitimate.

I independently derived a direct all-m source verification, retaining a
completely arbitrary fixed shift p and all four pencil coefficients. Put
`u=x-p`, `e=2m-1`,

    N=n0*u^(4m)+n1*u^(2m+1),
    D0=d0*u^(4m)+d1*u^(2m+1),
    b(H)=d1+4*n1*H.

In coordinates (x,V), `H=(V^2-D0)/(4N)` and the source Jacobian has
transition factor 2N. A direct differentiation gives

    J_(x,t)(H,V*u^(-2m))
      = u^(-2m)/2 * [(4mD0/u-D0')
                    +4H(4mN/u-N')].

The bracketed numerator is exactly

    e*(d1+4H*n1)*u^(2m).

Therefore multiplying by `2/[e*b(H)]` gives Jacobian one. I checked this
numerator independently for m=2..9 with all four coefficients left symbolic;
the displayed derivation itself proves every m. Since du=dx, no assumption
that the root translation is a global automorphism of W_m is involved.
The actual H is global by its own numerator box. Formula (2) is a rational
function on the source, and sufficiency only asserts that rational property.

This establishes both directions of the maximal-genus iff, including the
proportional-pencil case and the original, possibly nonzero linear source
coefficient P. The coefficient criterion is imposed on a genuine global H;
the proof does not claim arbitrary sparse N,D0 automatically have a global
polynomial realization.

## 5. Actual sharp family and the field-degree sidecar

For every m>=1 let e=2m-1 and delta!=0. The producer's family is

    h=x^m+x^(2m)t,
    H=(x^(2m)+delta*x)(1+x^m*t)^2+beta*h+q.

Its full other-chart values are `h=-b` and
`H=(1+delta*r^e)b^2-beta*b+q`, regular at every point. The three numerator
polynomials explicitly lie in the complete degree-2m box. Thus globality is
paid directly, including all special values beta=0 or q=0; no torus is
silently substituted for the surface.

The rational coordinates `u=x^(-e),h` satisfy `J(u,h)=-e` and
`H=(1+delta*u)h^2+beta*h+q`. This proves directly that
`G=1/(e*delta*h)` has the correct source Jacobian. Its pole at h=0 is
retained. The field map has degree exactly e: `C(x,t)=C(x,h)` and v=1/x
satisfies `v^e=u`, Eisenstein at u over `C(h)[u]`. This argument includes
composite e, and e=1 is the explicitly birational m=1 case. For larger m
the map is not mislabeled as a source automorphism.

The actual full discriminant pencil is

    (4lambda+beta^2-4q)x^(4m)
       +4delta*(lambda-q)x^(2m+1).

Both coefficients are generically nonzero. Its 2m branch points give genus
m-1 for every m, proving sharpness by actual global functions and actual
rational mates. At beta=0 the pencil is proportional, but still has the
same generic genus and mate. In particular m=3 gives genus two, so the
surface parameter and the discriminant degree cannot be dropped from an
elliptic bound inherited from m=2.

## 6. Recovered polynomial obstruction and its boundaries

I read the exact current
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md),
especially Sections3–4. It applies to any polynomial mate, with no bound on
its fibre or total degree. In the present notation it forces N=N0 in C* and
`Q-P^2/(4N0)` affine in x with nonzero slope. A targeted correction search
found no applicable retraction of this supplier.

The remaining incompatibility with globality is independently correct.
From `A*x^(2m)+B*x^m+C=N0`, with all three degrees at most 2m, the top
term forces deg A<=m. Writing E=A*x^m+B, the equation for C forces
`deg E<=m`. Conversely every such A,E gives exactly the full surviving box:

    B=-A*x^m+E, C=N0-x^m*E,
    P=A*x^m+E, Q=A.

If deg A=k>=1, then P has degree m+k without cancellation by E, so its
square dominates Q and the centered degree is 2(m+k)>1. If A is constant,
a nonconstant P gives even centered degree at least two, and a constant P
gives centered degree zero. Zero polynomials are included in these constant
cases. None permits a nonzero affine slope. The reasoning includes m=1.

Thus no genuine global quadratic H has a polynomial mate even if G need not
be globally regular on W_m. This is correctly attributed as a consequence of
the proved plane theorem plus the new global coefficient calculation, not a
new theorem about unrestricted planar Keller pairs. Both boundaries are real:
t is global of degree one and has polynomial mate -x; t^2 is global and
genuinely quadratic with rational mate -x/(2t), but has no polynomial mate.
The sharp family has rational mates and higher genus. These positive controls
prevent either deletion of the genuine-degree hypothesis or replacement of
“polynomial” by “rational.”

## 7. Source universe, independent controls and frozen acceptance

The entire standalone source uses exact symbolic cancellation and explicit
exceptions, preserving all gates under optimization. No inherited producer
implementation is imported. I checked each declared finite universe:

* Every monomial of each numerator box m=1..6, n=0..3 is tested in the full
  actual chart. Independently summing `(mn+1)(n+1)` gives exactly 480 rows.
* For m=1..4, all source monomials `x^i t^j`, i<=4m, j<=2 are retained in
  a matrix of every negative r-power coefficient, including its b-degree.
  Its kernel dimension, the exhibited basis rank and full kernel membership
  are checked. The ambient bounds follow from the complete filtration, not
  an unproved search cutoff. The dimensions are 9,15,21,27.
* The sharp families m=1..8 retain free delta,beta,q, their entire numerator
  box, discriminant, other chart, rational bracket and proportional branch.
  Symbolic identities with beta,q initially nonzero extend to zero because
  the cleared identities have no beta or q denominator; beta=0 is also
  checked explicitly.
* The four-coefficient sparse-pencil formula is checked directly in (x,V)
  for m=2..6, with primitive and residual-root separation identities. The
  independent shifted numerator calculation in Section4 supplements it.
* Every post-cancellation nonconstant-A degree pattern for m=1..12 is
  retained. The independent count `sum m(m+2)` is 806. Constant-A patterns
  and complete symbolic cancellations are checked separately.

The raw kernel and finite genus controls are supporting computations. The
complete filtration, all-m pole inequality, fixed-root descent, parameter
rationality, field degree and polynomial corollary are analytic arguments;
none is inferred from sampling m or coefficient values. The named special
pencil fibre explicitly distinguishes one repeated exceptional root from the
generic squarefree residual. All proof claims retain the scope needed by
the positive and hostile controls.

Independent replays from the worktree root:

    python3 04-computation/planar_jc48_sep08_dg_genus.py
    python3 -O 04-computation/planar_jc48_sep08_dg_genus.py

Both completed successfully and equal the frozen 856-byte output byte for
byte, including 1,550 gates and semantic hash
`c0e665ed29f7969341556b17978f101f4db88223d0ddc9bb5c4ba3691a0d5f9b`.
Source, output and primary were not edited by this auditor.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `04-computation/planar_jc48_sep08_dg_genus.py` | 9148 | `45f25cbe2dea8486f8988136a23d11d386672ac9a6445d9110c929711fd8cff8` |
| `05-knowledge/results/planar_jc48_sep08_dg_genus.out` | 856 | `386c8905a444b7184dffdff3dfc1694f841a8e7b27f3504026dc02875fd2d9f9` |
| Primary proof at acceptance, before status promotion | 19836 | `ccabd96c5832fd5bf8a98fe3c06fe27a6fb3da75f15a533e20aad5aa50b083fa` |

No repair remains pending. Root owns subsequent primary promotion and
integration; those status edits do not alter this accepted proof or the
frozen source/output pins.
