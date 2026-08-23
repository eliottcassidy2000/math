---
id: THM-3890
title: "Universal quintic common-zero resolvent class-group dichotomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let H4 be any
  nonzero homogeneous binary quartic and let L be a
  linear form.  If L divides H4, then H4+lambda L^5 is reducible.  Otherwise
  the normal quadratic surface W^2=H4+lambda L^5 has class group Z when H4
  is nonsquare and Z/5 when H4 is a square.  Its units are constant in both
  cases, so it has no three-torsion and supports no connected normal
  finite-flat cubic algebra with that discriminant.  Repeated-root quartics
  are included.  Consequently the irreducible one-projective-point
  degree-five common-zero route is empty for normal cubic covers.  Degree at
  least six,
  reducible or nonnormal cubic orders, unit-ideal nonmonogenic forms, Keller
  realization, and JC(2) remain open.
source: root / post-THM-3887 universal quintic closure, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED on 2026-08-23 by two disjoint audits.
  The proof checks irreducibility, branch normalization and normality without
  assuming that the quartic is squarefree.  Nagata localization is computed
  from the full unit group of k[t,y,(y^2-D)^-1].  In the square case the two
  localization units have boundary valuation vectors (3,-2) and (-2,3),
  whose Smith form is diag(1,5); in the nonsquare case the sole vector is
  (1,1).  The codimension-one tame-inertia/Kummer argument is inherited
  with all hypotheses exposed from THM-3865.  Both audits replayed the 35
  exact gates in normal and optimized modes and independently rederived the
  normalization, class groups, unit groups, genus budget, and Kummer
  exclusion.  They required the explicit non-tangency hypothesis and the
  irreducible-route scope.  Hostile tangent control
  H4=L^2M^2 shows that dropping non-tangency can instead expose Z/3 after
  normalization; that case is reducible and the raw quadratic ring is
  nonnormal.  The assertion-free companion checks the universal charts and
  tangent factor, both Smith forms, two repeated-root binary-cubic controls,
  and the projective parametrization.  Normal and optimized runs byte-match
  the frozen 35-gate transcript.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3865-one-place-inverse-discriminant-resolvent-class-group
  - THM-3887-binary-cubic-common-zero-quintic-one-tangent-obstruction
related:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3889-maximally-confluent-quadratic-binary-cubic-two-place-obstruction
script: 04-computation/jc2_universal_quintic_common_zero_resolvent_class_group_thm3890.py
output: 05-knowledge/results/jc2_universal_quintic_common_zero_resolvent_class_group_thm3890.out
script_sha256: 53a93831008db8ff6c144fca2994bac05f8b483de8758e1df269e38a2730756d
output_sha256: 2a685cadfa86c793337f9e63e872cb8095e61dcf452d9ac9d13805159db04ef1
semantic_sha256: b07d58be28a8d1ab2a06b38310aed8e5228c5744cb64a47c8f9b4003cc30025d
hash_basis: raw LF bytes
---

# THM-3890 -- the irreducible one-point common-zero quintic misses a normal cubic cover

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `k` be an algebraically closed field of characteristic zero, put

```text
R=k[M,L],
```

let `H4` be a nonzero homogeneous form of degree four, and let
`lambda in k*`.  Define

```text
delta=H4(M,L)+lambda L^5,
S=R[W]/(W^2-delta).                                      (1)
```

There are two exhaustive cases.

1. If `L|H4`, then `delta` is reducible.
2. If `L` does not divide `H4`, then `delta` is irreducible, `S` is a
   normal domain, `S*=k*`, and

```text
                 / Z,     if H4 is not a square,
Cl(S) isomorphic <
                 \ Z/5,   if H4 is a square.              (2)
```

In particular

```text
Cl(S)[3]=0,                                                (3)
```

and there is no connected normal finite-flat rank-three `R`-algebra whose
discriminant ideal is `(delta)`.

The conclusion includes every repeated-root factorization type of `H4`.
Squarefreeness of the quartic tangent cone is nowhere assumed.

## 1. The tangent case is a factor, not an escape

If `H4=L H3`, then

```text
delta=L(H3+lambda L^4).                                   (4)
```

Thus a direction that is simultaneously the sole degree-five direction and
a quartic tangent direction makes the quintic reducible.  The rest of the
proof assumes

```text
L does not divide H4.                                     (5)
```

Write

```text
H4(M,L)=L^4 D(M/L),
D(t)=kappa t^4+d3 t^3+d2 t^2+d1 t+d0,
kappa=H4(1,0) in k*.                                      (6)
```

Repeated roots of `D` are allowed.

The hypothesis `(5)` is load-bearing for the class-group formula, not just
for irreducibility.  If it is dropped, take `H4=L^2M^2`.  Then

```text
delta=L^2(M^2+lambda L^3),
```

the raw ring `W^2=delta` is nonnormal, and its quadratic-field normalization
has `z=W/L` and

```text
(z-M)(z+M)=lambda L^3.
```

Its boundary valuation determinant is three, not five.  This is the exact
place where tangent nonreduced geometry can carry three-torsion.  It does
not rescue the desired target because `delta` is already reducible.

## 2. The branch is irreducible and has an explicit normalization

After inverting `L`, the equation `delta=0` is, up to the unit `L^4`,

```text
D(t)+lambda L=0,                 t=M/L.                    (7)
```

The right side is linear in `L`, hence irreducible in `k[t,L]`.  Any
factorization of `delta` in `k[M,L]` would therefore acquire a factor
`cL^n` after localization.  Since `(5)` says that `L` does not divide
`delta`, `n=0`; hence `delta` is irreducible.

The polynomial map

```text
nu:A1_t -> V(delta),
L=-D(t)/lambda,                  M=-tD(t)/lambda             (8)
```

is birational because `t=M/L` off the origin.  It is finite: `t` satisfies
the monic scalar multiple of

```text
D(T)+lambda L=0                                                (9)
```

over the coordinate ring of `V(delta)`.  Since `k[t]` is normal, `(8)` is
the normalization.

The finite normalization addresses over `(M,L)=(0,0)` are exactly the
distinct roots of `D`.  Thus a repeated-root pencil merely changes the
finite address packet; it does not change the global argument.  If `D` is
squarefree, these are four smooth branches with distinct tangent lines, so
the origin is an ordinary quadruple point.  Its local delta invariant is

```text
binomial(4,2)=6,
```

exactly the full arithmetic genus `binomial(4,2)=6` of a plane quintic.
More generally, `(8)` shows that the projective normalization is `P1`, while
the branch is smooth away from the origin and at infinity; hence every root
partition of `D` still concentrates the full genus budget six at the common
zero.  The projective extension of `(8)` is

```text
[T:Z] |-> [-T D_h(T,Z) : -Z D_h(T,Z) : lambda Z^5],        (10)
```

where `D_h=H4(T,Z)`.  The sole parameter point `Z=0` maps to
`[1:0:0]`.  Consequently this family has both one projective point and one
normalization place at infinity.  The implication is proved here rather
than assumed; outside this family, one projective point need not mean one
place.

Equation `(7)` also shows that the branch is smooth wherever `L!=0`, since
its `L`-derivative is `lambda`.  On `L=0`, equation `(6)` leaves only the
origin.  Thus the branch curve has no singular point away from the origin.

## 3. The quadratic surface is normal and has a factorial localization

The hypersurface `S` is Cohen--Macaulay and hence `S2`.  Its singular locus
lies over the singular locus of the branch curve and is therefore contained
in the single closed point `(M,L,W)=(0,0,0)`.  This has codimension two in
the surface, so `S` is `R1` and therefore normal.  Since the irreducible
factor `delta` has odd valuation in `Frac(R)`, it is not a square; hence `S`
is a domain.

In `S_L`, put

```text
t=M/L,                  y=W/L^2,
F(t,y)=y^2-D(t).                                           (11)
```

The surface equation gives `F=lambda L`.  Conversely,

```text
L=F/lambda,             M=tF/lambda,
W=yF^2/lambda^2.                                           (12)
```

These maps are inverse and give

```text
S_L isomorphic to k[t,y,F^(-1)].                           (13)
```

The ring on the right is a localization of the UFD `k[t,y]`; therefore

```text
Cl(S_L)=0.                                                 (14)
```

At `L=0`, equation `(1)` becomes

```text
W^2=kappa M^4.
```

Choose `mu^2=kappa`.  The only height-one primes above `L` are

```text
P_+=(L,W-mu M^2),                P_-=(L,W+mu M^2),          (15)
```

and both have multiplicity one because `S/(L)` is reduced at their generic
points.  Hence

```text
div_S(L)=P_++P_-.                                         (16)
```

## 4. Nagata localization gives exactly `Z` or `Z/5`

Nagata's sequence for `(13)` is

```text
(S_L)*/S* -> Z[P_+] direct_sum Z[P_-]
           -> Cl(S) -> Cl(S_L) -> 0.                       (17)
```

It remains only to compute the first valuation map.

### 4.1. The nonsquare quartic

If `D` is not a square in `k[t]`, then `F=y^2-D` is irreducible and

```text
(S_L)*=k* F^Z.                                             (18)
```

Since `F=lambda L`, its valuation vector at `(P_+,P_-)` is `(1,1)`.
The valuation map has trivial kernel modulo `k*`, so `S*=k*`, and `(17)`
gives

```text
Cl(S)=Z^2/<(1,1)> isomorphic to Z.                         (19)
```

This includes nonsquare quartics with double or triple roots.

### 4.2. The square quartic

If `D=q^2`, then `q` has degree two by `(6)`, and

```text
H4=Q^2,                  Q(M,L)=L^2q(M/L).                 (20)
```

Now

```text
F=(y-q)(y+q),
(S_L)*=k*(y-q)^Z*(y+q)^Z.                                 (21)
```

On `S`, put `X=W-Q` and `Zeta=W+Q`.  Then

```text
X Zeta=lambda L^5.                                        (22)
```

At the boundary prime where `X` vanishes, `Zeta` is a generic unit, so
`v(X)=5`; at the other prime, `v(X)=0`.  Dividing by `L^2` to recover
`y-q=X/L^2` gives the valuation vector `(3,-2)`.  Similarly `y+q` gives
`(-2,3)`.  Thus the image lattice in `(17)` has matrix

```text
V = [ 3 -2 ]
    [ -2 3 ],                 Smith(V)=diag(1,5).           (23)
```

The matrix has trivial kernel, so again `S*=k*`, while

```text
Cl(S)=coker(V) isomorphic to Z/5.                           (24)
```

The determinant five is the exact split-localization/fifth-contact effect
anticipated by the tangent-family scout.  The boundary here is reduced; it
produces five-torsion, not three-torsion.  The genuinely nonreduced tangent
control is the normalization example in Section 1.
Equations `(19),(24)` prove `(2),(3)` for every quartic tangent cone subject
to the load-bearing non-tangency `(5)`.

Here `Cl(S)` is the Weil divisor class group.  Equivalently it is the
Picard group of the regular locus of `Spec S`; no equality with `Pic(S)` is
claimed.  The deleted chart `(13)` itself is factorial.

## 5. The missing three-torsion forbids a normal cubic

Suppose that `T` were a connected normal finite-flat rank-three `R`-algebra
with discriminant ideal `(delta)`.  Its fraction ring is a separable cubic
field `E/K`, `K=k(M,L)`.  Because `delta` is nonsquare, the Galois closure
`N/K` has group `S3`, with quadratic subfield

```text
K(sqrt(delta))=Frac(S).                                    (25)
```

The cyclic extension `N/Frac(S)` has degree three.  It is unramified over
every height-one prime of `S`.  Away from `delta`, this follows from the
unit discriminant.  At the generic point of `delta`, the discriminant
valuation is one.  Tame cubic ramification is therefore `(2,1)`, so the
inertia in `S3` is a transposition and intersects `A3` trivially.  Passing
to the quadratic subfield removes the only inertia.

Since `k` contains the cube roots of unity, Kummer theory writes

```text
N=Frac(S)(cuberoot(h))                                     (26)
```

for some `h in Frac(S)*`.  Codimension-one unramifiedness gives

```text
v_P(h)=0 mod 3                    for every P in S^(1).     (27)
```

Hence `div(h)=3D0`.  Equation `(3)` implies that `D0` is principal.
Dividing `h` by a cube leaves a divisor-zero rational function.  Normality
makes that function and its inverse regular on `Spec S`, so it is a unit.
The unit computation above makes it an element of `k*`, and algebraic
closedness makes that scalar a cube.  The extension `(26)` is trivial, a
contradiction.

## 6. Degree-five common-zero consequence

Let

```text
I=aX^3+bX^2Y+cXY^2+dY^3,                 a,b,c,d in k[u,v]
```

be a binary cubic whose four coefficients vanish at a point `p`, and let
`Delta=Disc(I)` be nonzero.  Assume that `Delta` has total degree at most
five, is irreducible, and that its projective closure has only one point on
the line at infinity.

Translate `p` to the origin.  Homogeneity of the cubic discriminant gives
`Delta in (u,v)^4`.  A degree-four discriminant, or one of order five, would
be a homogeneous binary form and hence reducible over `k`.  Therefore
`deg Delta=5`, its order is four, and

```text
Delta=H4+lambda L^5.                                      (28)
```

The unique projective infinity point is exactly the statement that the
degree-five leading form is one linear fifth power.  If `L|H4`, equation
`(4)` contradicts irreducibility.  Otherwise Sections 2--5 show that no
connected normal finite-flat cubic algebra can have discriminant `(Delta)`.

Thus

```text
common coefficient zero + irreducible one-point discriminant
  + connected normal finite-flat cubic
       ==> discriminant degree at least six.               (29)
```

For a degree-five common-zero index form considered only as a curve,
THM-3887 further says that the finite packet in `(8)` has at least two
distinct addresses; generically it has four.  The present theorem closes
all such packets for a normal cubic cover, including every repeated-root
quartic pencil.

## 7. Exact boundary and replay

This theorem does not exclude:

- degree-six or higher common-zero discriminants;
- reducible discriminants or nonnormal cubic orders;
- an intrinsically nonmonogenic binary index form whose coefficients
  generate the unit ideal and whose form represents no scalar unit;
- a different discriminant geometry without the one-point degree-five
  normal form; or
- a planar Keller map or `JC(2)` counterexample.

The distinction between one projective point and one normalization place is
load-bearing in general.  It is harmless only here: the explicit
normalization `(8)--(10)` proves the one-place property.  Equivalently, in
the projective equation `Z H4+lambda L^5=0`, the unique infinity point
`[1:0:0]` satisfies

```text
partial_Z(Z H4+lambda L^5)=H4(1,0)=kappa!=0,
```

so it is smooth and has exactly one normalization place.

Reproduce the exact algebra with

```bash
python3 04-computation/jc2_universal_quintic_common_zero_resolvent_class_group_thm3890.py
python3 -O 04-computation/jc2_universal_quintic_common_zero_resolvent_class_group_thm3890.py
```

Both modes must byte-match the frozen output named in the metadata. **QED.**
