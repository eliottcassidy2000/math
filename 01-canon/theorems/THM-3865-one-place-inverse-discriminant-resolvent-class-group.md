---
id: THM-3865
title: "Every non-tangent fifth-degree inverse discriminant has torsion-free quadratic-resolvent class group"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  nonzero lambda and every
  linear L not proportional to one of the four tangent factors of Delta_0,
  the normal affine surface W^2=Delta_0+lambda L^5 has divisor class group
  Z.  Consequently it has no three-torsion and supports no connected normal
  finite-flat cubic algebra over k[A,C] with that discriminant.  This closes
  both THM-3853 targets and their full non-tangent family; it is not a general
  inverse-discriminant or planar-Jacobian obstruction.
source: root / inverse-discriminant global algebraization lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).  The
  audit independently checked normality from the immersive four-address
  normalization, both directions of the factorial localization chart, the
  irreducible localization divisor and its full unit group, and the two
  reduced multiplicity-one primes over C.  It rederived the exact Nagata
  sequence, with no hidden relation, and checked the decisive tame-local
  step: discriminant exponent one gives transposition inertia, whose
  intersection with A3 is trivial over the quadratic subfield.  It then
  verified the Kummer/class-group/unit argument, including possible
  unramified residue degree and codimension-two leakage.  The proof uses the THM-3853 normalization
  to establish normality, an explicit UFD chart after inverting C, the two
  reduced height-one primes over C, Nagata's localization sequence, and a
  codimension-one tame-inertia/Kummer descent.  The assertion-free companion
  freezes the chart and inverse chart, irreducible localization divisor,
  special fibre, four-address normalization, and binary-cubic control.
  Normal and optimized replays byte-match the frozen 30-gate transcript and
  both recorded hashes.  A second hostile audit of the coordinate-free
  strengthening (jc_quartic_c3_construct, 2026-08-23) checked the tangent-line
  criterion, degree and squarefreeness of D_L, normalization and infinity
  place, both localization maps, the kappa sign, the two reduced primes, and
  transfer of the Kummer argument.  The companion includes L=A+C as an
  independent orientation.
depends_on:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
  - THM-3853-quadratic-depth-inverse-discriminant-one-place-gluing-obstruction
  - THM-3855-formal-inverse-discriminant-lift-and-algebraization-gate
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff
script: 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
output: 05-knowledge/results/jc2_one_place_inverse_discriminant_class_group_thm3865.out
script_sha256: ce1dbc2de8815f95e4f822824e91e9229f2161676b73e4fe4bfb627ab1e317d3
output_sha256: e8262877989c12ee98298bf67270ef085e1cf802293734ca83b15f7522c53f58
semantic_sha256: 7c662f08ea830a5690d2b47f1e8057279521611080700778941e27b07d6e5648
hash_basis: raw LF bytes
---

# THM-3865 -- the formal one-place lifts have a global class-group obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be an
algebraically closed field of characteristic zero, let `lambda in k*`, and
put

```text
Delta_0=A(C+5A)(4C+19A)(3C-17A),
delta=Delta_0+lambda C^5,
R=k[A,C],
S=R[W]/(W^2-delta).                                           (1)
```

More generally, let `L` be any nonzero linear form not proportional to one
of

```text
A,                 C+5A,                 4C+19A,
3C-17A.                                                     (1a)
```

Put `delta_(L,lambda)=Delta_0+lambda L^5` and
`S_(L,lambda)=R[W]/(W^2-delta_(L,lambda))`.  Then every
`S_(L,lambda)` is a normal domain and

```text
Cl(S_(L,lambda)) isomorphic to Z.                            (2)
```

Choose a complementary linear coordinate `M` and define `kappa in k*` by
`Delta_0(M,0)=kappa M^4`.  The two height-one primes

```text
P_+=(L,W-mu M^2),                 P_-=(L,W+mu M^2),
mu^2=kappa,                                                   (3)
```

generate the class group, with the sole relation

```text
[P_+]+[P_-]=0.                                                (4)
```

In particular `Cl(S_(L,lambda))[3]=0` and `S_(L,lambda)*=k*`.

Consequently there is no connected normal finite-flat rank-three
`R`-algebra whose discriminant ideal is `(delta_(L,lambda))`.  Equivalently,
no polynomial binary-cubic row with discriminant a scalar multiple of any
member of this family can have a normal domain as its Delone--Faddeev
algebra.

This consequence is deliberately exact-target.  THM-3855 constructs a
connected normal nonmonogenic `S3` cubic algebra after completing at
`(A,C)`, and even shows that every fixed-linear formal coefficient row is
locally gauge-equivalent to the THM-3808 packet.  The present theorem proves
that this local cubic cannot algebraize over the whole affine target with
the same discriminant.  Other one-place curves, nonnormal cubic orders,
different tangent cones, and planar Keller maps remain open.

Sections 1--4 give every detail for `L=C,M=A`, where `kappa=-1615`.
Section 5 records why the proof is unchanged for arbitrary non-tangent `L`.

## 1. The branch and the quadratic surface are normal

Put

```text
D(t)=t(5t+1)(19t+4)(3-17t).                                  (5)
```

Then `Delta_0(tC,C)=C^4D(t)`.  THM-3853 gives the normalization of
`V(delta)`:

```text
C=-D(t)/lambda,                     A=tC.                     (6)
```

The four roots of `D` are distinct.  The map `(6)` is immersive: if both
coordinate derivatives vanished, then `D(t)=D'(t)=0`.  It is injective
away from `C=0`, because there `t=A/C`, and exactly the four roots of `D`
map to `(A,C)=(0,0)`.  Thus the irreducible branch curve has only one affine
singular point, the ordinary four-branch origin.

The hypersurface `S` is Cohen--Macaulay, hence satisfies Serre `S2`.  Its
singular locus is cut out by

```text
W=delta=partial_A(delta)=partial_C(delta)=0,                  (7)
```

so the preceding branch audit makes it zero-dimensional.  Therefore `S`
satisfies `R1` and is normal.  Since `delta` is irreducible and not a square,
`S` is a domain.

## 2. Inverting `C` exposes a factorial chart

In `S_C` define

```text
t=A/C,                         y=W/C^2,
F(t,y)=y^2-D(t).                                             (8)
```

Equation `(1)` becomes

```text
F=lambda C.                                                   (9)
```

Conversely

```text
C=F/lambda,               A=tF/lambda,
W=yF^2/lambda^2.                                             (10)
```

Equations `(8)--(10)` are inverse ring maps and give

```text
S_C isomorphic to k[t,y,F^(-1)].                              (11)
```

The polynomial `D` has four simple roots, so it is not a square in `k[t]`.
Hence `F=y^2-D(t)` is irreducible.  The right side of `(11)` is a
localization of the UFD `k[t,y]`, and therefore

```text
Cl(S_C)=0,
(S_C)*=k* F^Z=k* C^Z.                                        (12)
```

This chart also explains the geometry hidden by the completed calculation.
The rational map is a plane chart, while the affine elliptic curve
`y^2=D(t)` is contracted to the four-branch origin.  Its three-torsion is
available formally, but need not extend to a global Weil class on `S`.

## 3. Nagata's sequence computes the global class group

At `C=0`,

```text
delta(A,0)=-1615 A^4,
S/(C)=k[A,W]/((W-mu A^2)(W+mu A^2)).                         (13)
```

Thus `(3)` are exactly the height-one primes containing `C`.  They are
distinct and occur with multiplicity one.  In particular

```text
div_S(C)=P_++P_-.                                             (14)
```

Any unit of `S` becomes `cF^n` in `(12)`.  Its valuations at `P_+` and
`P_-` are both `n`, whereas a unit has valuation zero.  Hence `n=0` and

```text
S*=k*.                                                        (15)
```

Nagata localization for the normal domain `S` now reads

```text
(S_C)*/S* -> Z[P_+] direct_sum Z[P_-] -> Cl(S) -> Cl(S_C) -> 0.
                                                                    (16)
```

By `(12),(14)`, the image of the first arrow is exactly the diagonal
subgroup generated by `(1,1)`, and the last group vanishes.  Therefore

```text
Cl(S)=Z^2/<(1,1)> isomorphic to Z,                             (17)
```

which proves `(2)--(4)`.

## 4. Why torsion-free `Cl(S)` forbids the global cubic

Suppose, for contradiction, that `T` were a connected normal finite-flat
rank-three `R`-algebra with discriminant ideal `(delta)`.  Its fraction
field is a separable cubic field.  Since `delta` is not a square, the
generic Galois closure `L/K`, for `K=k(A,C)`, has group `S3`, and its
quadratic subfield is

```text
M=K(sqrt(delta))=Frac(S).                                    (18)
```

The cyclic extension `L/M` has degree three.  It is unramified at every
height-one prime of `S`.  Away from `(delta)`, this follows from the unit
discriminant of the finite-flat cubic algebra.  At the generic point of
`(delta)`, the discriminant valuation is one; tame cubic ramification is
therefore of type `(2,1)`, so inertia in the `S3` closure is a transposition.
Its intersection with `A3=Gal(L/M)` is trivial.  Passing to the quadratic
subfield consequently removes the only inertia.

Because `k` contains the cube roots of unity, Kummer theory writes

```text
L=M(cuberoot(h))                         for some h in M*.     (19)
```

Unramifiedness at every height-one prime forces

```text
v_P(h)=0 mod 3                         for every P in S^(1).   (20)
```

Thus `div_S(h)=3D_0` for a Weil divisor `D_0`, and `(17)` says `D_0` is
principal.  After dividing `h` by a cube, one obtains a divisor-zero
rational function, hence a unit of the normal affine domain `S`.  By `(15)`
that unit lies in `k*`, and algebraic closedness makes it a cube.  Equation
`(19)` is then trivial, a contradiction.

This proves the global cubic nonexistence and closes the polynomial
algebraization question for the targets `Delta_0+lambda C^5` in this
orientation.  It does not turn the local formal obstruction into a general
theorem about all one-place discriminants.  **QED for `L=C`.**

## 5. Every non-tangent direction has the same two-prime chart

Now choose any `L` satisfying `(1a)` and any complementary linear coordinate
`M`.  Homogeneity gives

```text
Delta_0(M,L)=L^4 D_L(M/L),
D_L(t)=Delta_0(t,1) in the (M,L) coordinates.                 (21)
```

The four distinct tangent factors of `Delta_0` become four distinct roots
of `D_L`.  The non-tangency assumption is exactly what makes all four roots
finite, or equivalently

```text
deg D_L=4,                       lc(D_L)=kappa!=0.             (22)
```

Thus every branch argument in Section 1 repeats with

```text
L=-D_L(t)/lambda,                       M=tL.                  (23)
```

In particular the curve is irreducible, has four distinct addresses over
the origin and one place at infinity, and its quadratic double surface is a
normal domain.

After inverting `L`, put

```text
t=M/L,                       y=W/L^2,
F_L=y^2-D_L(t).                                                (24)
```

The inverse chart is

```text
L=F_L/lambda,             M=tF_L/lambda,
W=yF_L^2/lambda^2,                                           (25)
```

so `S_(L,lambda)[L^(-1)]=k[t,y,F_L^(-1)]` is again factorial
with units `k*F_L^Z`.  At `L=0`,

```text
delta_(L,lambda)(M,0)=kappa M^4,
W^2-kappa M^4=(W-mu M^2)(W+mu M^2),       mu^2=kappa.         (26)
```

These are the only two height-one primes over `L`, both reduced and of
multiplicity one.  Nagata's sequence is therefore again
`Z^2/<(1,1)>`, and the unit group is again `k*`.  Section 4 uses only these
two conclusions and the fact that the branch divisor is irreducible and
generically reduced, so its tame-inertia/Kummer proof applies unchanged.
This proves `(2)--(4)` and the cubic nonexistence for every non-tangent
`L`, including both `L=C` and `L=A+C` from THM-3853.  **QED.**

## 6. Exact replay

```bash
python3 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
python3 -O 04-computation/jc2_one_place_inverse_discriminant_class_group_thm3865.py
```

Both modes must byte-match the frozen transcript in the metadata.
