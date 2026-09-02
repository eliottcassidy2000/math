---
id: THM-4328
title: "Seam-covariant Student--Stein face visibility"
status: >
  PROVED RELATIVE TO THM-3992/3997/4298/4308/4315 + VERIFIED-EXACT.  In the
  normalized reduced (2,3) seam, every applicable row tangent cokernel has a canonical
  algebraic Student--Stein functional depending only on the fifth-root
  invariant kappa=4gamma^2/a.  Composed with THM-4298's labelled weighted-
  face flag, the functionals retain every channel at positive even weight
  and kill every channel at odd weight.  At M=12 they give an invertible three-
  coordinate system and exact equations for the U, Z, Lambda, and D walls.
  Explicit U=0 and Z=0 points survive the finite row-eight projection and
  row-nine bracket gate, while the direct contribution of an odd weight-19
  perturbation is invisible to every scalar Stein functional.  This is not seam entry, an all-row lift,
  endpoint-wall extinction, JC(2), or DC(2).
source: root + seam_stein_transport / planar-Jacobian endpoint-wall session, 2026-09-01
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4298-weighted-face-source-normal-unimodular-visibility-transform
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-539
primary_script: 04-computation/jc2_seam_covariant_student_stein_face_visibility_thm4328.py
primary_output: 05-knowledge/results/jc2_seam_covariant_student_stein_face_visibility_thm4328.out
primary_script_sha256: 083a624585750ac79c3462a39c920d88b413d36baff925558d7af0873d622ddd
primary_output_sha256: e29315d178b7fbc59420e8a06bc99283f04cec5f26f56a4efa0b2c0d48b46e5b
semantic_sha256: 95527d31e118108182eb8c6f367e8465906803d40ea788f9dc3726d4b3e66c20
hash_basis: raw LF bytes; semantic digest over the certificate ledger
audit: >
  PASS.  A standard-library Fraction certificate executes 71,016
  optimization-stable gates.  It checks the row tangent and cokernel through
  m=24, moment recurrence, fifth-root action, stationary density and filter,
  all even/odd face channels through weight 72, the M=12 determinant and
  four wall equations, exact evaluation of both imported finite-gate
  formulas, transversality, and the direct weight-19 parity hostile.
  Normal and optimized streams byte-match
  the frozen transcript.
---

# THM-4328 -- Seam-covariant Student--Stein face visibility

**PROVED RELATIVE TO THM-3992/3997/4298/4308/4315 + VERIFIED-EXACT.  THE STOCHASTIC
INTERPRETATION IS REAL-ANALYTIC ONLY FOR `kappa>0`; THE COKERNEL FUNCTIONAL
AND FACE-VISIBILITY THEOREMS ARE ALGEBRAIC IN CHARACTERISTIC ZERO.  THIS
DOES NOT PROVE SEAM ENTRY, AN ALL-ROW LIFT, ENDPOINT-WALL EXTINCTION,
`JC(2)`, OR `DC(2)`.**

## 1. Statement and inheritance

Work in THM-3992/3997's normalized reduced `(2,3)` cell.  Before imposing
the live seam, write

```text
q=gamma*x^2+3a/(2gamma),                 a*gamma != 0,
kappa=4gamma^2/a.                                      (1)
```

The general-gauge operator is a formal transport of THM-4308's normalized
one.  Put `q0(y)=-(y^2+6)/2` and `y=sqrt(kappa)*x`.  Then
`q=-(a/(2gamma))q0(y)`.  Under `y=lambda*x` and `q=c*q0`, the chain rule
sends

```text
q' theta-(q/m)theta'
  -> c*lambda [q0' thetahat-(q0/m)thetahat'].           (1a)
```

Thus THM-4308's identity transports after adjoining `sqrt(kappa)`, and the
rational formula below descends to the original characteristic-zero field.

For every algebraic index `m>=1`, define the operator below.  In the bracket
continuation where the tangent is free (`m>=2`, and `m>=5` in the current
application), changing `v_(m-1)` by `theta v_0'`, with
`deg(theta)<=m-1`, changes the row-`m` compatibility by

```text
D_m theta=q' theta-(q/m)theta'
         =-(a/(2gamma)) S_(m,kappa) theta,

S_(m,kappa)theta=
 [(kappa*x^2+6)theta'-2m*kappa*x*theta]/(2m).          (2)
```

The theorem has four parts.

1. There is a unique normalized algebraic functional `ell_(m,kappa)` on
   `k[x]_(<=m)` whose kernel is `im(D_m)`, with explicit rational moments.
2. It is invariant under the residual fifth-root gauge, and for real
   `kappa>0` it is expectation under an explicit Student/Pearson stationary
   law with an exact row-to-row retention filter.
3. Applied row by row to any nonempty positive residual weight-`M` face flag
   from THM-4298, it is lossless for even `M` and identically zero for odd `M`.
4. At `M=12` it supplies exact endpoint-wall coordinates, but explicit
   finite-gate survivors and an odd perturbation show why the scalar
   stochastic observer cannot replace the labelled face or Newton geometry.

The inheritance pass was:

- closest proved mechanism: THM-4315's normalized `kappa=1` Student--Stein
  cokernel at one row;
- canonical hostile: THM-4315's characteristic-19 finite-row survivor;
- corrected near miss: equal finite-dimensional row fibres do not transport
  source or depth labels, and MISTAKE-539 forbids assigning a double section
  the ambient cubic's two-dimensional deformation algebra;
- least-used sidecar: THM-4298's labelled face channel before scalar
  expectation.

The live concept board was

```text
{row tangent, stationary law, fifth-root quotient, labelled face,
 parity firewall, endpoint walls, Newton component}.                 (3)
```

## 2. The covariant algebraic cokernel

Define `ell=ell_(m,kappa)` on monomials of degree at most `m` by

```text
ell(1)=1,                         ell(x^(2r+1))=0,

ell(x^(2r))=
 6^r (2r-1)!! /
 [kappa^r product_(j=1)^r(2m-2j+1)]                  (4)
```

whenever the displayed monomial has degree at most `m`.  For
`theta=x^d`, equation `(2)` gives

```text
2m S_(m,kappa)(x^d)
 =6d*x^(d-1)+(d-2m)kappa*x^(d+1).                     (5)
```

Odd `d=2r-1` and the recurrence in `(4)` show
`ell(S_(m,kappa)(x^d))=0`; even `d` gives two odd monomials and vanishes
termwise.  Thus `im(D_m)` is contained in `ker(ell)`.

The coefficient of `x^(d+1)` in `(5)` is nonzero for
`0<=d<=m-1`.  Hence `D_m:k[x]_(<=m-1)->k[x]_(<=m)` is injective.  Its image
and `ker(ell)` both have dimension `m`, so

```text
f in im(D_m)                 iff ell_(m,kappa)(f)=0.     (6)
```

This proves the algebraic cokernel statement for every `m>=1` over every
characteristic-zero field; no positivity or integration is used.  Its
bracket-tangent reading starts only at rows where `v_(m-1)` is a free
tangent, as specified above.

The residual fifth-root action of THM-3997 is

```text
(gamma,a)->(zeta*gamma,zeta^2*a),          zeta^5=1.    (7)
```

It fixes `kappa`, while the scalar `-a/(2gamma)` in `(2)` is multiplied by
`zeta`.  Therefore `D_m` changes by a unit and its image and normalized
cokernel functional are invariant.  On the live seam

```text
gamma=-a^3/2,                         kappa=a^5=A5,     (8)
```

so the functional descends literally to THM-3997's fifth-root quotient.

## 3. Exact Student/Pearson interpretation

Assume now that the coefficients are real and `kappa>0`.  Let

```text
mu_(m,kappa)(dx)=c_(m,kappa)
                 (kappa*x^2+6)^(-(m+1)) dx.            (9)
```

Its moments are exactly `(4)`.  It is stationary for

```text
L_(m,kappa)h=(kappa*x^2+6)h''-2m*kappa*xh',            (10)
```

because, for its density `rho`,

```text
d[(kappa*x^2+6)rho]/dx=-2m*kappa*x*rho.                (11)
```

Thus `(6)` is precisely the polynomial Stein identity for this diffusion.
If `X` has law `mu_(m,kappa)`, retain it with probability

```text
6/(kappa*X^2+6).                                        (12)
```

The survival probability is `(2m+1)/(2m+2)`, and the conditional survivor
has law `mu_(m+1,kappa)`.  For real `kappa<0`, `(9)` has real poles and is
not a probability law; only the algebraic functional `(4)` survives.  This
sign boundary is part of the theorem, not an analytic continuation claim.

## 4. Parity decides weighted-face visibility

Let a nonempty positive residual weight-`M` face (`M>=2`) have THM-4298's ordered monomials

```text
(i_s,j_s),                    2i_s+3j_s=M,
```

and minimal source-normal flag `(h_s)`.  Its labelled contribution at row

```text
m_s=(M+j_s)/2                                          (13)
```

is `h_s x^(j_s)`.  Apply the row-specific functional
`ell_(m_s,kappa)` and put

```text
e_s=ell_(m_s,kappa)(h_s x^(j_s)).                       (14)
```

Since `j_s` has the same parity as `M`, equations `(4)` and `(14)` give the
sharp dichotomy

```text
M even:  e_s=c_s(kappa) h_s with every c_s(kappa)!=0;
M odd:   e_s=0 for every s.                             (15)
```

Consequently the combined Student observer is a diagonal rescaling of the
entire labelled face flag at even weight, but annihilates the entire face at
odd weight.  The result uses a different row functional for each labelled
channel.  Summing the channels first would destroy precisely the label that
makes the even-weight statement lossless.

## 5. Exact `M=12` stochastic wall coordinates

For

```text
R_12=U p^6+W p^3y^2+Z y^4,                             (16)
```

THM-4298 gives

```text
(h0,h1,h2)=(U,6U+W,15U+5W+Z).                          (17)
```

The three channels have `(m,j)=(6,0),(7,2),(8,4)`.  Hence

```text
e0=h0,
e1=6h1/(13kappa),
e2=36h2/(65kappa^2),                                   (18)

det(d(e0,e1,e2)/d(h0,h1,h2))=216/(845kappa^3)!=0.      (19)
```

Substituting `(18)` in THM-4298's inverse wall equations gives

```text
U=0       iff e0=0,

Z=0       iff 13kappa^2 e2-78kappa e1+108e0=0,

Lambda=0  iff 65kappa^2 e2-312kappa e1+360e0=0,

D=0       iff 169kappa^2 e1^2+624kappa e0e1-864e0^2
                  -260kappa^2 e0e2=0.                  (20)
```

Thus stochastic expectations give honest coordinates on the even
weight-twelve face.  They do not say that a point satisfying one equation in
`(20)` lifts through later source rows or that its Newton components admit a
nonconstant Keller map.

## 6. Two finite-gate survivors separate rows from geometry

On THM-4308/4315's normalized `kappa=1` slice, the exact certificate
evaluates the inherited response and bracket formulas and intersects the wall equations `(20)` with the
THM-4308 projected row-eight gate and THM-4315's row-nine bracket equation.
Both intersections are nonempty and transverse in the inherited
`(Phi,eta,xi_10,alpha_11)` response chart.

One `U=0` control has

```text
Phi=1, eta=0, xi_10=237757952/54675,
upsilon_5=-731648/2025,
alpha_11=-189841784364917/5646560625,
W=-169749098233/9841500,
Z=5200771070534/66430125,                              (21)
```

with

```text
E5=alpha_11^2-4W upsilon_5
  =35245115050720811582989632889/31883646891800390625 !=0
```

and `Z*D*Lambda!=0`.  One `Z=0` control has

```text
Phi=1, eta=0,
xi_10=1563243041171/115860375,
alpha_11=207913052665031843/2393096045625,
beta_11=1, zeta_3=-3/2,
U=-5200771070534/1042743375,
W=10155617023591/463441500,
U+W=70597468930183/4170973500 !=0.                     (22)
```

The transcript records the remaining inherited coordinates and checks both
finite gates exactly.  At the controls, the decisive derivatives include

```text
dU/dxi_10=-6/11,
dZ/dxi_10=-22886/2673,
dE9/dalpha_11=-511211250 Phi !=0.                       (23)
```

These points are controls, not Keller candidates: they certify only
projected depth through row eight and the scalar row-nine bracket condition.
They do not certify full row-nine depth, all later rows, or a regular-model
map.  THM-4327 in fact excludes the displayed endpoint strata by Newton
component geometry.

## 7. The odd-weight entry firewall

Add the nonzero residual perturbation

```text
epsilon p^8y,                         epsilon!=0.        (24)
```

It has weight `19`.  In source-normal coordinates,

```text
p^8y=x*t^10(1+x^2t)^9.                                  (25)
```

It changes no row through nine.  Every direct later term in `(25)` has odd
`x`-degree, so every row-specific Student functional kills that isolated
face contribution by `(4)`.  But
THM-4298's labelled weight-19 face flag is nonzero and lossless before the
scalar expectation.  This is a literal face-observer hostile:

```text
all row-wise Stein observations of the isolated direct face vanish
        does not imply that labelled residual face vanishes.             (26)
```

Nonlinear later bracket predictions can mix this perturbation with changed
tangents; no claim is made that the full all-row obstruction sequence is
unchanged.  Therefore the stochastic process supplies exact cokernel
coordinates and a useful row transition, but cannot prove seam entry without
the labelled face/Hasse sidecar.  Endpoint extinction additionally needs the Newton
component and good-differential sidecar of THM-4327.

## 8. Typed connection ledger and scope

```text
source:     row-m bracket defect polynomial
target:     algebraic expectation coordinate ell_(m,kappa)
map:        quotient by im(D_m), represented by moments (4)
preserved:  solvability of the degree-capped row tangent equation
destroyed:  tangent primitive, source/depth label, Newton owner, component
sidecar:    degree cap plus THM-4298's ordered face flag and later depth rows
cheap test: the moment recurrence (5), then the odd perturbation (25).
```

What is proved is the algebraic cokernel formula, fifth-root covariance,
real positive-`kappa` diffusion/filter, parity visibility theorem, `M=12`
wall dictionary, and the two finite-gate controls.  What remains open is
entry into the reduced seam, any all-row consequence from the scalar
process, the endpoint subwalls not covered by THM-4327, and `JC(2)` itself.

## 9. Reproduction

The exact script and frozen transcript named in the frontmatter reconstruct
`(2)--(26)` over `Fraction` arithmetic, including the moment recurrence,
dimension count, gauge action, retention ratio, parity diagonal, wall
dictionary, controls, transversality, and odd perturbation.  Normal and
optimized runs are required to byte-match the frozen transcript.
