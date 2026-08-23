---
id: THM-3850
title: "Nonconstant cubic profiles force punctures on every irreducible branch"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
  INDEPENDENT HOSTILE AUDIT.  For every nonconstant polynomial profile b(C)
  whose depressed-cubic branch Delta_b is irreducible, the affine branch
  normalization is its smooth projective hyperelliptic model with exactly
  deg(rad b)+nu_infinity points deleted.  Here nu_infinity is one or two
  according as the squarefree part of 9b-2C^2 has odd or even degree.  The
  total is at least two, so no such branch is a polynomial curve.  The
  minimal reducible profile b=kappa C splits off an affine line but leaves a
  residual rational component with three punctures.  Arbitrary reducible
  profiles remain open.
source: jc_sparse_direct_search / nonlinear cubic profile and branch-normalization lane, 2026-08-23
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3842-nonlinear-cubic-tower-trace-shift-eightfold-base-change
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
script: 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
output: 05-knowledge/results/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.out
script_sha256: 0d7212f699c1d81457bf3167c454d3d85cf1a6b9ccfe35045643498b81e2a726
output_sha256: 29da8c7562e9a49c2a691e986176183c18fc1e4145f3b18f703685d5df50cb38
semantic_sha256: 6642655c1b115fb099da9d9e02d89e7c7ba19311428813b2c40f60e1ae082498
hash_basis: raw LF bytes
---

# THM-3850 -- every irreducible nonconstant profile pays at least two places

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING
INDEPENDENT HOSTILE AUDIT.**  Work over an algebraically closed field `k` of
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
is a hostile control, not a classification of every reducible profile.

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

## 5. Hostile controls and scope

The constant profile of THM-3847 has `s=0` and squarefree `L` of even degree
two, giving `0+2=2` punctures.  It lies just outside the nonconstant theorem
and is the sharp two-place control.

For nonconstant controls:

```text
b=C+1:          s=1, deg ell=2,       |D_b|=1+2=3;
b=(C-1)^2:      s=1, deg ell=2,       |D_b|=1+2=3;
b=C^3+1:        s=3, deg ell=3,       |D_b|=3+1=4.                 (26)
```

The middle row checks that multiplicity of a profile root does not inflate
the puncture count; only `rad(b)` appears.

The assertion-free exact companion verifies the universal discriminants,
the two local numerator values `(18)`, the smooth finite point `(19)`--`(20)`,
all three controls `(26)`, and the full minimal factorization and conic
normalization `(21)`--`(25)`.  It has 33 active gates; normal and optimized
replay must byte-match the frozen transcript.

Arbitrary **reducible** profiles beyond `(21)` remain open.  In particular
this theorem does not classify how several branch components can share or
trade punctures, and it does not prove `JC(2)`.  Its exact consequence is
that a successful nonconstant-profile cubic completion cannot keep an
irreducible branch and hope to turn it into a polynomial nonproperness curve;
it must exploit a reducible branch packet and solve every component.

Reproduction:

```bash
python3 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
python3 -O 04-computation/jc2_nonconstant_cubic_profile_branch_punctures_thm3850.py
```
