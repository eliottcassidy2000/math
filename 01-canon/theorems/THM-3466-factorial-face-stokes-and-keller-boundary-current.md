---
id: THM-3466
title: "Factorial face Stokes, Keller flux, and the boundary-current Krylov block"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The n-variable
  factorial functional has an exact Boolean face descent:
  B_I(h)=L_n(prod_(i in I)(1-partial_i)h).  Thus at a factorial-null power,
  every gradient-multiplier observation is exactly a codimension-one face
  debt; on the coordinate interior ideal, every fixed pure-derivative bank
  is eventually blind because the first accessible boundary jet grows
  linearly with the power.  For a planar Keller pair, weighted Hamiltonian integration by parts
  makes the corresponding multiplier/face recurrence endogenous.  On a
  triangle, a constant-J polynomial g=p+iq obeys
  integral_boundary g^n d(conjugate g)=-2nic integral_Delta g^(n-1).
  Under HFC nullity this boundary response is a single nonzero dipole and its
  operation-kernel quotient is exactly one length-two nilpotent Krylov block.
  The origin lies on the boundary image or in a bounded complementary cell,
  and g^2 d(conjugate g) has a pointed polynomial primitive on the source
  boundary.  The scalar packet is only a reindexed moment sequence and the
  primitive need not separate sheets; no FC(3), HFC(3), or JC(2) conclusion
  is asserted.
source: root/factorial-jacobian-alternation/2026-08-15
audit: >
  Independent derivations checked all signs in the orthant, quadrant, and
  triangle formulas; the current response, Cauchy-transform consequence,
  Krylov independence, polynomial degree bound, anti-tangent first-jet
  boundary, and slit-annulus hostile.  The exact companion independently
  integrates all three triangle edges and checks identity, triangular, and
  equilateral controls through response order four.
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
  - THM-3328-boundary-cone-overlap-and-anti-tangent-keller-passport
related:
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
  - THM-2856-sparse-factorial-lowering-decoder-and-laguerre-rank-one-defect
  - THM-3383-terminal-monomial-cone-polynomiality-fork
script: 04-computation/factorial_face_stokes_keller_boundary_current_thm3466.py
output: 05-knowledge/results/factorial_face_stokes_keller_boundary_current_thm3466.out
script_sha256: c964dd0a33264bf182087180c05c52c459aa6942e41c019374e3bfe0d72e9716
output_sha256: 0801e63667b1acfcd739bd318a6f1174932bc35cf53869dccaeb6d8644cf0d7a
semantic_sha256: e925107c9d9a14f7434ada9c2195e761cae97861f76930402aa37ac81afb9527
hash_basis: raw bytes
---

# THM-3466 -- factorial face Stokes and the Keller boundary-current block

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Factorial integration by parts is a face map

For `h in C[x_1,...,x_n]` put

```text
L_n(h)=integral_(R_+^n) h(x) exp(-sum_i x_i) dx,
```

so `L_n(x^alpha)=alpha!`.  For a coordinate set `I`, let `B_I(h)` be
the lower-dimensional factorial functional after setting `x_i=0` for every
`i in I`.  In particular, integration by parts in one coordinate gives

```text
L_n(partial_i h)=L_n(h)-B_{i}(h).                                  (1)
```

The boundary at infinity vanishes because a polynomial is multiplied by an
exponential.  The coordinate operators commute, so iterating `(1)` proves the
Boolean face formula

```text
B_I(h)=L_n(prod_(i in I)(1-partial_i)h).                            (2)
```

This includes every codimension, including `B_empty=L_n` and evaluation at
the origin when `I` is the full coordinate set.

For a power `f^m`, `(1)` becomes

```text
B_i(f^m)=L_n(f^m)-m L_n(f^(m-1) partial_i f).                       (3)
```

Hence at an FC-null orbit,

```text
L_n(f^m)=0  =>  m L_n(f^(m-1) partial_i f)=-B_i(f^m).               (4)
```

Equation `(4)` identifies the missing multiplier observation exactly: it is
the negative face moment, not a new scalar consequence of nullity.

At `n=1`, `(3)` is THM-2848/2856's endpoint identity.  THM-2856 proves
that bare chain lowerings contain only that endpoint modulo scalar moments.
The multivariable formula does not evade the access boundary.  If `x_i`
divides `f`, then the face term is identically zero; on any finite scalar-null
window the corresponding lowerings vanish too.  Thus face descent alone is
not a detector for the interior ideal `(x_1...x_n)`.

There is a stronger jet-level failure.  Suppose

```text
ord_(x_i)(f)=e_i>=1                         (1<=i<=n).       (4a)
```

For every multiindex `beta` satisfying `beta_i<m e_i` for all `i`, repeated
use of `(1)` has no boundary term: every boundary jet encountered before
`partial^beta(f^m)` still contains a positive power of the relevant
coordinate.  Hence

```text
L_n(partial^beta(f^m))=L_n(f^m).                            (4b)
```

On an FC-null orbit, every fixed finite bank of pure derivatives is therefore
eventually identically zero; its order must grow at least linearly with `m`
to reach the first possible boundary jet.  This is a response-depth
obstruction, not only a first-face obstruction.

The phenomenon already has a finite exact-prefix control.  Let `H` be the
nonzero support-`{1,2,3,4}` polynomial from THM-2846/2856 with

```text
H(0)=0,             L_1(H)=L_1(H^2)=L_1(H^3)=0,
```

and put `f(x)=H(x_1)x_2...x_n`.  Then, for `m=1,2,3`,

```text
L_n(f^m)=L_1(H^m)(m!)^(n-1)=0,                             (4c)
```

and every derivative in the region `beta_i<m e_i` is blind by `(4b)`.
Polynomial multipliers remain external: when the `i`-face of `h` vanishes,

```text
L_n(g partial_i h)=L_n((g-partial_i g)h),                  (4d)
```

which does not reduce an arbitrary multiplier to scalar moments.

## 2. A Keller bracket endogenizes a weighted multiplier recurrence

In two variables write

```text
{P,H}=P_x H_y-P_y H_x,
E_x(H)=integral_0^infinity H(0,y)e^(-y)dy,
E_y(H)=integral_0^infinity H(x,0)e^(-x)dx.                         (5)
```

The Hamiltonian vector field `X_P=(-P_y,P_x)` is divergence-free and
`X_P(x+y)=P_x-P_y`.  Apply the divergence theorem to
`e^(-x-y)H X_P` on the positive quadrant.  The two coordinate faces have
outward fluxes `H P_y` and `-H P_x`, respectively, and the flux at
infinity is zero.  Therefore

```text
L_2({P,H})
 =L_2(H(P_x-P_y))+E_x(H P_y)-E_y(H P_x).                           (6)
```

If `{P,Q}=c in C^*` and `H=P^aQ^b`, this gives the exact factorial
Keller recurrences

```text
bc L_2(P^aQ^(b-1))
 =L_2(P^aQ^b(P_x-P_y))
  +E_x(P^aQ^bP_y)-E_y(P^aQ^bP_x),                                 (7)

-ac L_2(P^(a-1)Q^b)
 =L_2(P^aQ^b(Q_x-Q_y))
  +E_x(P^aQ^bQ_y)-E_y(P^aQ^bQ_x).                                 (8)
```

Equations `(7)--(8)` are the direct full-FC/JC interface.  They make a
Hamiltonian derivative observable, but they retain both the exponential
drift and the labelled faces.  Dropping either term gives a false closure.
For the triangular automorphism `P=x,Q=y+x^2`, the `y=0` face term is
already nonzero and is load-bearing.

## 3. Triangle moments are a boundary current

Let `Delta` be a positively oriented triangle, let

```text
g=p+iq,                     J={p,q},
I_m=integral_Delta g^m dA.                                         (9)
```

For every `n>=1`,

```text
d(g^n d(conjugate g))
 =n g^(n-1) dg wedge d(conjugate g)
 =-2ni g^(n-1)J dx wedge dy.                                      (10)
```

Stokes therefore gives, without assuming `J` constant,

```text
integral_boundary g^n d(conjugate g)
 =-2ni integral_Delta g^(n-1)J dA.                                (11)
```

Also `g^n dg=d(g^(n+1)/(n+1))` has zero boundary integral.  Resolving
`dg=dp+i dq` and `d(conjugate g)=dp-i dq` gives

```text
integral_boundary g^n dp=-ni integral_Delta g^(n-1)J dA,
integral_boundary g^n dq= n  integral_Delta g^(n-1)J dA.            (12)
```

Now impose the constant-Jacobian condition `J=c!=0`.  Then

```text
integral_boundary g^n d(conjugate g)=-2nic I_(n-1).                 (13)
```

If `g` is an HFC null candidate, `I_m=0` for every `m>=1` while
`I_0=Area(Delta)`.  Its intrinsic boundary response is exactly

```text
n=0:   0,
n=1:  -2ic Area(Delta),
n>=2:  0.                                                          (14)
```

This is a first-nonzero current packet.  It keeps the oriented edge source
that scalar area moments forget.

## 4. The operation-twisted kernel is exactly one Jordan string

Let `V` be the complex vector space of piecewise-polynomial one-forms on the
labelled source boundary.  Define

```text
C(alpha)=integral_boundary alpha,       T(alpha)=g alpha,
K_infinity=intersection_(j>=0) T^(-j)(ker C).                        (15)
```

Take `alpha=d(conjugate g)`.  Equation `(14)` says

```text
C(alpha)=0,              C(T alpha)=R=-2ic Area(Delta)!=0,
C(T^j alpha)=0 for every j>=2.                                    (16)
```

Thus `T^2 alpha in K_infinity`.  Moreover `alpha,T alpha` are linearly
independent modulo `K_infinity`: if
`a alpha+bT alpha in K_infinity`, applying `C` gives `bR=0` and
applying `C T` then gives `aR=0`.  Hence their pointed Krylov
quotient is exactly

```text
[alpha]  --T-->  [T alpha]  --T-->  0,                             (17)
```

one nilpotent Jordan string of length two.  This is an exact
operation-response statement, not a metaphor and not a finite-prefix
inference: terminal absorption uses every HFC moment.

## 5. The origin-cell constraint

For large `zeta`, expand `1/(zeta-g)` and use `(14)`:

```text
integral_boundary d(conjugate g)/(zeta-g)
 =-2ic Area(Delta)/zeta^2.                                        (18)
```

Both sides are analytic on the unbounded component of
`C minus g(boundary Delta)`, so `(18)` continues throughout that component.
If zero lay there, it would be off the compact boundary image and the left
side would be analytic near zero, while the right side has a nonzero
second-order pole.  Therefore

```text
0 lies on g(boundary Delta)
or in a bounded complementary cell of g(boundary Delta).           (19)
```

This augments THM-3328's collision/cone passport with a distinguished target
cell.  It does not say which edge owns the cell or how many interior sheets
cover it.

## 6. The pointed primitive and its exact limitation

Since `C(T^2 alpha)=0`, choose a boundary basepoint and define

```text
kappa(s)=integral_(basepoint)^s g^2 d(conjugate g).                  (20)
```

The vanishing total integral makes `kappa` periodic.  On an affine edge,
`G(t)=g(gamma(t))` has degree at most `d=deg g` and

```text
kappa'(t)=G(t)^2 conjugate(G'(t)),                                  (21)
```

so `kappa` is polynomial of degree at most `3d` on each labelled edge.
At a collision, `kappa(a)-kappa(b)` is a pointed clutch that the target value
`g(a)=g(b)` does not record.

No universal separation is claimed.  Where `G'!=0`,

```text
d kappa/dG=G^2 conjugate(G')/G'.                                    (22)
```

At THM-3328's no-overlap anti-tangent collision
`G_2'=-lambda G_1'` with `lambda>0`, the two first jets in `(22)`
agree automatically.  Differentiating once more along the target gives

```text
d^2 kappa/dG^2
 =2G conjugate(G')/G'
  +G^2(conjugate(G'')G'-conjugate(G')G'')/(G')^3.          (22a)
```

The first term also agrees across an anti-tangent pair, so only the signed
curvature term can split its second lifted jet.  At `G=0`, both `(22)` and
`(22a)` vanish.  The next lawful observer is therefore a curvature/clutch
jet, not another tangent scalar.

The slit-annulus control is sharp again.  On its two identified radial seams,
`d(conjugate g)=dg` and

```text
kappa=(g^3-rho^3)/3                                                   (23)
```

on both sides, so the primitive descends and does not separate the sheets.

## 7. Hostile controls and scope

- **Nonconstant Jacobian:** `(11)` contains the multiplier `J`.  Scalar HFC
  moments cannot replace it.
- **Finite response prefix:** the affine equilateral coordinate has
  `I_1=I_2=0` but `I_3=1/20`.  Its responses at `n=2,3` vanish and
  the response at `n=4` is nonzero.
- **Topology:** the slit annulus realizes `(19)` with zero in a bounded hole
  and realizes the primitive descent `(23)`.
- **Access:** `(2)--(4)` expose face debts; they do not supply multiplier
  observations from scalar moments.

The theorem proves the displayed identities, the Krylov block, and the
origin-cell alternative.  It proves no injectivity, inverse polynomiality,
HFC(3), full FC, JC(2), or universal clutch separation.

## 8. Information contract

```text
source:      factorial powers / a labelled triangle boundary
target:      coordinate faces / a boundary-current Krylov module
maps:        prod(1-partial_i), Hamiltonian flux, multiplication by g
preserved:   exact bulk moment and constant-J bracket response
destroyed:   interior support after face restriction; global sheet gluing
sidecars:    face labels, exponential drift, boundary orientation, basepoint
tests:       triangular quadrant map, equilateral finite prefix, slit annulus
```

## 9. Exact verification

Run

```text
python3 04-computation/factorial_face_stokes_keller_boundary_current_thm3466.py
python3 -O 04-computation/factorial_face_stokes_keller_boundary_current_thm3466.py
```

and compare raw bytes with the declared output.  The companion performs exact
factorial sums and exact polynomial edge integration.  It checks all faces of
a three-variable control, both Keller flux recurrences on identity and
triangular maps, and all four current identities through order four on three
triangle maps.  Normal and optimized runs are byte-identical.

**QED.**
