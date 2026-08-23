---
id: THM-3850
title: "Every nonconstant cubic profile has a nonpolynomial branch component"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  nonconstant polynomial profile b(C),
  at least one irreducible component of the depressed-cubic branch has
  affine normalization missing at least two projective points.  In the
  irreducible case the exact count is deg(rad b)+nu_infinity.  Primitive
  reducible branches split into two explicit rational denominator graphs;
  origin-vanishing profiles force a residual component through the finite-C,
  A-infinity corner.  Individual affine-line components can occur, but the
  whole branch packet never consists only of polynomial curves.
source: jc_sparse_direct_search / nonlinear cubic profile and branch-normalization lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct and root /
  jc-cohn3709, 2026-08-23).  The audits checked the primitive/irreducible
  equivalence, exact finite punctures and infinity parity, the minimal
  three-puncture residual, both primitive denominator graphs, and every
  origin-vanishing valuation row.  In the last case they separately verified
  that the component through `(C=0,A=infinity)` is nonvertical and acquires a
  distinct normalization point over `C=infinity`, even if projective branches
  collide.  Normal and optimized runs byte-match the frozen 47-gate
  transcript and both hashes.
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3842-nonlinear-cubic-tower-trace-shift-eightfold-base-change
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
script: 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
output: 05-knowledge/results/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.out
script_sha256: b5b673c26ba0fc7c32faa0ad6e4936feaec2a1680e04a8531953ba5f1eb21cde
output_sha256: 767e4102e8b45d3a9f952ad3beddb28074d1f781e78cca0217d5a2572b932be3
semantic_sha256: 951ceb8509eff3102902aa0e9a077a7a95b415fb55c08018c0936ddfed1aa4aa
hash_basis: raw LF bytes
---

# THM-3850 -- every nonconstant profile has a bad branch component

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of
characteristic zero.  Let `b=b(C) in k[C]` be nonzero and define

```text
p=3/2+AC,                         u=1+AC+A^2b(C),                 (1)

Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                         (2)
```

Then

```text
8p^3-27u^2=A^2 Delta_b,
disc_A(Delta_b)=-8L^3,                    L=9b-2C^2.              (3)
```

For every nonconstant `b`, the reduced plane curve `V(Delta_b)_red` has at
least one irreducible component whose affine normalization is obtained from
its smooth projective normalization by deleting at least two distinct
points.  Consequently the full branch packet never consists solely of
polynomial curves.  This remains true even when other components are affine
lines.

Assume that the affine branch curve

```text
Gamma_b=V(Delta_b) subset A2_(A,C)                               (4)
```

is irreducible.  Write

```text
L=h^2 ell,                    ell squarefree,                    (5)
r=deg ell,                    s=deg rad(b),                      (6)
```

where scalar square factors are absorbed into `h` and
`rad(b)=b/gcd(b,b')` up to a nonzero scalar.  Then `r>=1`, the smooth
projective normalization of `Gamma_b` is the smooth projective completion of
the hyperelliptic affine model

```text
H_b: W^2=ell(C),                                                 (7)
```

and the **affine normalization** of `(4)` is exactly

```text
Gamma_b^nu=H_b minus D_b,                                       (8)

|D_b|=s+nu_infinity,

nu_infinity = 1  if r is odd,
              2  if r is even.                                 (9)
```

More precisely, `D_b` contains one point over each distinct root of `b` and
all points of `H_b` over `C=infinity`; there are no other deleted points.
Its projective genus is

```text
g(H_b)=floor((r-1)/2).                                          (10)
```

Because `b` is nonconstant and `(4)` is irreducible,

```text
s>=1,                         nu_infinity>=1,                    (11)
```

so `|D_b|>=2`.  In particular no irreducible nonconstant-profile branch has
affine normalization `A1`; it cannot itself be a polynomially parametrized
Jelonek component.

At the first reducible boundary, `b=kappa C` with `kappa in k*`, the reduced
branch is the disjoint union of

```text
C=0 ~= A1_A                                                     (12)
```

and an irreducible residual rational curve whose affine normalization is
`P1` minus three points.  Thus the vertical line is a genuine one-place
component, but it is accompanied by a three-place component.  This boundary
is the smallest positive control showing that an individual one-place
component may occur without repairing the full packet.

The theorem distinguishes three different objects that must not be merged:
the singular affine branch `(4)`, its smooth affine normalization `(8)`, and
its smooth projective normalization `(7)`.  The integer in `(9)` counts the
projective points missing from `(8)`.

## 1. Irreducibility forces a nonzero profile root

If `b(0)=0`, every term in `(2)` is divisible by `C`.  Since `(2)` has
positive `A`-degree, this makes `Delta_b` reducible.  Under the irreducibility
hypothesis,

```text
b(0)!=0.                                                         (13)
```

In particular `gcd(b,C)=1`.  The quadratic `(2)` is primitive as a
polynomial in `A` over `k[C]`: if an irreducible polynomial divided all three
`A`-coefficients, the leading coefficient would force it to divide `b`;
modulo that divisor the other two coefficients become `8C^3` and `9C^2`,
forcing it to divide `C`, contrary to `(13)`.

Therefore Gauss's lemma and `(3)` give the exact criterion

```text
Delta_b irreducible  iff  L is not a square in k(C).              (14)
```

Indeed `-8` is a square in `k`, and `L^3` is a square precisely when `L` is
a square.  Thus the squarefree part `ell` in `(5)` is nonconstant, proving
`r>=1`.  Since `k` is algebraically closed and `b` is nonconstant, `(13)`
also gives at least one distinct nonzero root of `b`, proving `s>=1`.

## 2. Projective normalization and its infinity fibre

Choose `epsilon in k` with `epsilon^2=-8`, and first use the possibly
singular quadratic model

```text
w^2=L.                                                           (15)
```

The quadratic formula for `(2)` is

```text
A={8C^3-54Cb-epsilon Lw}/{54b^2}.                                (16)
```

Since `w=hW`, equations `(5)`, `(7)`, and `(16)` identify the function field
of `Gamma_b` with `k(C)(W)`, `W^2=ell`.  Taking the smooth projective model
gives `(7)` and the standard hyperelliptic genus formula `(10)`.

The fibre of `C:H_b->P1` over infinity has one point when `r` is odd and two
points when `r` is even.  In the even case the two leading square roots are
distinct because `k` is algebraically closed and the leading coefficient of
`ell` is nonzero.  Every such point is absent from the affine normalization:
the coordinate `C` itself has a pole there.  This proves the infinity term
in `(9)`.

## 3. Exactly one finite puncture per distinct root of b

Let `c` be any root of `b`.  By `(13)`, `c!=0`, and

```text
L(c)=-2c^2!=0.                                                    (17)
```

Hence the double cover `(7)` has two distinct, unramified points over `c`.
Choose `sigma in k` with `sigma^2=-2` and take `epsilon=2sigma`.  In `(16)`
the two numerator values are

```text
8c^3-epsilon(-2c^2)( sigma c)=0,
8c^3-epsilon(-2c^2)(-sigma c)=16c^3.                             (18)
```

Thus `A` has a pole at one of the two points.  The other point is present in
the affine branch: direct specialization gives

```text
Delta_b(A,c)=c^2(8Ac+9),                                         (19)
```

so the unique finite point is

```text
A=-9/(8c),                     partial_A Delta_b=8c^3!=0.          (20)
```

It is already smooth and has exactly one normalization preimage.  Therefore
each **distinct** root of `b`, independent of its multiplicity, deletes
exactly one point of `H_b`.

Formula `(16)` has no finite pole away from `b=0`.  Roots of `h` or `ell`
may be finite ramification or singular-resolution points of the quadratic
model `(15)`, but their denominator in `(16)` is a unit, so they remain in
the affine normalization.  Consequently the points in Sections 2 and 3 are
all the deleted points, proving the exact formula `(8)`--`(9)`.

## 4. The minimal reducible profile

Set `b=kappa C`, `kappa in k*`.  Exact factorization gives

```text
Delta_(kappa C)=C H_kappa,                                       (21)

H_kappa=-27kappa^2A^2C+8AC^2-54kappa AC+9C-54kappa.              (22)
```

Since

```text
H_kappa mod C=-54kappa,                                          (23)
```

the vertical line and residual curve are comaximal and disjoint in the
affine branch.  The residual quadratic has

```text
disc_A(H_kappa)=-8C(9kappa-2C)^3.                                (24)
```

This is nonsquare in `k(C)`, so `H_kappa` is irreducible.  Its smooth
projective normalization is the conic

```text
V^2=C(9kappa-2C).                                                 (25)
```

The conic has two points over `C=infinity`.  It has one ramified point
`C=V=0`, but `(23)` says the residual affine curve has no point over `C=0`.
Equivalently, the quadratic formula has numerator of order one in the local
parameter `V` and denominator of order two, so `A` has a pole there.  These
are its only missing points: one over zero and two over infinity.  Hence its
affine normalization is `P1` minus three points, while the separate vertical
component `(12)` is `A1`.

This is the smallest-degree demonstration that producing one acceptable
branch component does not repair the whole branch packet.

## 5. Primitive reducible profiles are two denominator graphs

Suppose now that `b(0)!=0` but `Delta_b` is reducible.  The primitivity
argument in Section 1 still applies, so `(3)` forces

```text
L=h^2.                                                           (26)
```

Choose `sigma in k` with `sigma^2=-2` and put

```text
u=h-sigma C,                 v=h+sigma C,
uv=h^2+2C^2=9b.                                                  (27)
```

The two quadratic numerators factor without a remainder, giving

```text
Delta_b=-(1/3) E_u E_v,                                          (28)

E_v=v^2A-3(C-sigma h),
E_u=u^2A-3(C+sigma h).                                           (29)
```

These are distinct irreducible rational components.  Indeed `u(0)v(0)` is
nonzero.  The coefficient and constant term of each linear polynomial in
`A` are coprime by `(30)`, so each is irreducible.  If the two factors were
associates, `u/v` would be constant; then `(27)` would make `h` proportional
to `C`, contradicting `h(0)^2=9b(0)!=0`.  At a root `c` of `v`, one has
`h(c)=-sigma c` and

```text
C-sigma h=-c!=0;                                                 (30)
```

dually, at a root of `u`, `C+sigma h=-c!=0`.  Thus the two affine
normalizations are exactly

```text
(E_v)^nu = P1 minus ({roots of v} union {infinity}),
(E_u)^nu = P1 minus ({roots of u} union {infinity}),             (31)

number of punctures = 1+deg rad(v),  1+deg rad(u), respectively. (32)
```

Since `uv=9b` is nonconstant, at least one of `u,v` is nonconstant.  The
corresponding component has at least two punctures.  If one factor is a
scalar, its corresponding component is genuinely an affine line; for example

```text
h=sigma C+a,       a in k*,       u=a,       v=2sigma C+a,       (33)
```

and `E_u=0` is the polynomial graph

```text
a^2A+3C-3sigma a=0,                                               (34)
```

while `E_v=0` has one finite denominator puncture and the point at infinity.
This is the primitive counterpart of the vertical-line control `(21)`.

## 6. Origin-vanishing profiles force a projective corner

It remains to allow `b(0)=0`.  Write

```text
b=C^m a(C),                    m>=1,        a(0)!=0.              (35)
```

Let `nu=ord_C(Delta_b)`.  Reading the three `A`-coefficients in `(2)` gives
the complete valuation table

```text
nu=1,                              m=1;
nu=2,                              m=2 and a(0)!=1/6;
nu=3,                              m=2 and a(0)=1/6;
nu=2,                              m>=3.                          (36)
```

The exceptional value in the third row is harmless: after the constant term
cancels to order at least three, the coefficient of `A` has leading term
`-AC^3`, so the order is exactly three.

Put `H=Delta_b/C^nu` and homogenize only the `A`-coordinate:

```text
Hbar(A,Z;C)=d_2(C)A^2+d_1(C)AZ+d_0(C)Z^2.                       (37)
```

In every row of `(36)`,

```text
ord_C(d_2)=2m-nu>0.                                              (38)
```

Therefore the residual projective closure contains

```text
P_0=(C=0,[A:Z]=[1:0]).                                          (39)
```

The order `nu` is exact, so `C` does not divide `H`.  Hence the vertical
fibre `C=0` is not a component of the residual divisor.  A hypersurface in
the smooth surface `P1_A x A1_C` has no isolated component, so choose any
reduced irreducible curve component through `P_0`; it is nonvertical and
therefore dominates the `C`-line.

Complete that component in `P1_A x P1_C` and normalize it.  It has a point
above `P_0`, missing from the original affine branch because `A=infinity`.
The nonconstant projective map to `P1_C` is surjective, so it also has at
least one point above `C=infinity`.  These normalization points are distinct:
their `C`-values are `0` and `infinity`.  This remains true if several
projective branches collide at either corner, because normalization separates
the branches while preserving their two different `C`-values.  Thus this
residual component has at least two punctures.

Sections 1--3, 5, and 6 exhaust respectively the irreducible,
primitive-reducible, and origin-vanishing cases.  They prove the complete
all-profile assertion at the start of the theorem.

## 7. Hostile controls and scope

The constant profile of THM-3847 has `s=0` and squarefree `L` of even degree
two, giving `0+2=2` punctures.  It lies just outside the nonconstant theorem
and is the sharp two-place control.

For nonconstant controls:

```text
b=C+1:          s=1, deg ell=2,       |D_b|=1+2=3;
b=(C-1)^2:      s=1, deg ell=2,       |D_b|=1+2=3;
b=C^3+1:        s=3, deg ell=3,       |D_b|=3+1=4.                 (40)
```

The middle row checks that multiplicity of a profile root does not inflate
the puncture count; only `rad(b)` appears.

The assertion-free exact companion verifies the universal discriminants,
the two local numerator values `(18)`, the smooth finite point `(19)`--`(20)`,
all three controls `(40)`, and the full minimal factorization and conic
normalization `(21)`--`(25)`, the factorization `(28)`, the affine-line
positive control `(33)`--`(34)`, and every valuation row in `(36)`--`(38)`.
It has 47 active gates; normal and optimized replay must byte-match the
frozen transcript.

This theorem does not assert that every branch component is nonpolynomial:
`(12)` and `(34)` show the contrary.  It proves that at least one component
is nonpolynomial for every nonconstant profile, so the **whole** branch
packet can never be a union solely of polynomially parametrized curves.  It
does not classify arbitrary target deformations beyond `(1)`, and it does
not prove `JC(2)`.

Reproduction:

```bash
python3 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
python3 -O 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
```
