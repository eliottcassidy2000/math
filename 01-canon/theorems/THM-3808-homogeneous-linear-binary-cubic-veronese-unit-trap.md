---
id: THM-3808
title: "Homogeneous linear binary cubic Veronese unit trap"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  A generic
  Delone--Faddeev cubic over k[A,C] whose four binary-cubic coefficients are
  homogeneous linear forms and whose discriminant is squarefree is forced
  to be the third-Veronese cone.  Its maximal etale open has a rank-four
  group of nonconstant units.  An explicit rational four-branch packet is a
  normal nonmonogenic S3 cubic with exactly one visible principal companion
  per branch, so it meets every algebraic gate in THM-3801 except the
  constant-unit hypothesis.  This is a construction boundary, not a planar
  Jacobian counterexample.
source: jc_quartic_c3_construct / binary-index design lane, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The proof checks the Delone--Faddeev sign
  convention, squarefree-discriminant normality, the four-point
  Riemann--Hurwitz calculation, reconstruction of the full O_P1(3) section
  ring, the exact arrangement-complement invariant ring and its unit
  lattice, and every branch/companion factor in the rational example.  The
  exact companion verifies all polynomial identities, trace and index
  determinants, parametrization, ramification, critical values, unit-lattice
  basis, and boundary controls.  Independent hostile audit remains due.
depends_on: []
related:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
script: 04-computation/jc2_homogeneous_binary_cubic_veronese_unit_trap_thm3808.py
output: 05-knowledge/results/jc2_homogeneous_binary_cubic_veronese_unit_trap_thm3808.out
script_sha256: c9a607f905a96b049eb3d615f4696c4d25aa4de0fd33eae5a59a1dca8b383fc5
output_sha256: 26bc0b5800258e70510b68564495957485026092c5ceb28ada68b0f5cd8ad745
semantic_sha256: e2b8d32b332c278b697c28ce32f36babad76061cde9a587d2957fbb0ab140a69
hash_basis: raw LF bytes
---

# THM-3808 -- the homogeneous-linear cubic construction is a Veronese unit trap

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  Let `k` be
an algebraically closed field of characteristic zero and put

```text
R=k[A,C],                         deg A=deg C=1.                     (1)
```

For four linear forms `a,b,c,d in R_1`, put on the free module

```text
S=R*1 direct_sum R*omega direct_sum R*theta                         (2)
```

the Delone--Faddeev multiplication law

```text
omega^2=-ac+b omega-a theta,
omega theta=-ad,
theta^2=-bd+d omega-c theta.                                       (3)
```

Suppose the generic algebra `S tensor_R k(A,C)` is a field and the binary
cubic discriminant

```text
Delta=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd                          (4)
```

is nonzero and squarefree.  Then:

1. `S` is a normal standard-graded domain and

```text
S ~= k[s,t]^(3)=direct_sum_(n>=0) k[s,t]_(3n);                     (5)
```

2. its binary index form is

```text
I(x,y)=-(a x^3+b x^2y+c xy^2+d y^3),                              (6)
```

   so `I(R^2)` never meets `k*`; the only nonmonogenic fibre is the
   square-zero fibre over `(A,C)=(0,0)`;
3. the projective cubic cover has four distinct simple branch points and
   Galois closure `S3`; over each branch there is one double ramification
   ray and one simple companion ray; but
4. after deleting the four ramification rays, the maximal etale open
   `U=Spec B` has

```text
B*/k* ~= { (n_0,n_1,n_2,n_3) in Z^4 :
            n_0+n_1+n_2+n_3 = 0 mod 3 },                           (7)
```

   a free abelian group of rank four.

Thus **every** homogeneous-linear, squarefree binary-index attempt fails the
constant-unit premise `B*=k*` of THM-3801.  The failure is not the
discriminant, monodromy, nonmonogenic locus, or companion-sheet grammar; it
is the radial unit created when a torsion ray is deleted.

## 1. Index, square-zero fibre, and normality

In the basis `(1,omega,theta)`, multiplication by `omega` and `theta` is

```text
M_omega = [[0,-ac,-ad],       M_theta = [[0,-ad,-bd],
           [1, b,  0],                   [0,  0,  d],
           [0,-a,  0]],                  [1,  0, -c]].              (8)
```

The determinant of the three columns `(1,T,T^2)`, for
`T=x omega+y theta`, is exactly `(6)`, and the determinant of the trace
pairing is `(4)`.  Since all four coefficients vanish at the homogeneous
vertex, every value of `(6)` belongs to `(A,C)` and hence cannot be a unit.
Modulo `(A,C)`, law `(3)` becomes

```text
S/(A,C)S = k[omega,theta]/(omega,theta)^2.                          (9)
```

Conversely, squarefreeness of `(4)` forces the linear forms `a,b,c,d` to
span `R_1`: if they had a common linear factor, `Delta` would have its fourth
power.  Away from the vertex the fibre index cubic is not identically zero,
so it takes a nonzero value over the algebraically closed residue field.
Lifting that pair `(x,y)` makes the index a local unit.  Thus `(9)` is the
exact local nonmonogenic locus.

The generic-field hypothesis makes the torsion-free ring `(2)` a domain.
It is Cohen--Macaulay because it is free over the regular ring `R`.  At a
height-one prime away from `Delta`, it is etale and regular.  At a factor of
`Delta`, the DVR discriminant has valuation one.  The discriminant-index
formula between this order and its integral closure changes valuation by
twice the index length, so the order is already maximal there.  Hence `S`
satisfies `R1` and `S2`, and is normal.

## 2. Riemann--Hurwitz forces the third-Veronese cone

Let

```text
P=Proj S -> Proj R=P1.                                             (10)
```

Normality makes `P` a smooth projective curve.  The map has degree three.
Because `(4)` is a squarefree homogeneous quartic, it has four distinct
zeros on the base; valuation one gives one tame transposition over each.
Riemann--Hurwitz therefore gives

```text
2g(P)-2=3*(-2)+4=-2,             so P ~= P1.                       (11)
```

Moreover `O_P(1)` is the pullback of `O_P1(1)`, hence has degree three and
is `O_P1(3)`.  From the graded free module `(2)`, with
`deg omega=deg theta=1`,

```text
dim_k S_n=(n+1)+2n=3n+1       for n>=1.                            (12)
```

This equals `h^0(P1,O(3n))`.  The natural inclusion of `S_n` into these
sections is therefore equality in every degree.  Multiplication is
preserved, proving `(5)`.  Equivalently, every packet in the theorem is a
third-Veronese cone equipped with a different base-point-free pencil of
cubics.

The four simple inertia elements are transpositions.  Since the generic
cubic field is connected, its transitive monodromy group is `S3`.

## 3. Exact unit lattice on the maximal etale open

Use `(5)` and let the scalar group `mu_3` act on `A2_(s,t)`.  The four
ramification points of `(10)` lift to distinct lines

```text
ell_0=0, ell_1=0, ell_2=0, ell_3=0,       H=ell_0ell_1ell_2ell_3.  (13)
```

The vertex lies on every ray.  Away from it, the quotient by `mu_3` is
etale, so the maximal etale locus of `Spec S -> A2` is exactly the cone with
these four rays deleted.  It is the affine principal open

```text
U=D(H^3),
B=S[1/H^3]=k[s,t,1/H]^(mu_3).                                     (14)
```

The localization before taking invariants is a UFD with

```text
k[s,t,1/H]*=k* ell_0^Z ell_1^Z ell_2^Z ell_3^Z.                   (15)
```

Scalar `zeta in mu_3` multiplies every `ell_i` by `zeta`.  A monomial in
`(15)` is invariant exactly when the exponent sum is divisible by three.
Finite invariants preserve units in both directions, proving `(7)`.

There are two pieces in this rank-four obstruction.  The three ratios
`ell_i/ell_0` give the usual unit lattice of `P1` minus four points.  The
fourth unit is radial: `ell_0^3`.  In divisor language

```text
Cl(S)=Z/3,                   div(ell_i^3)=3E_i,                     (16)
```

so deleting even one ramification ray turns its torsion multiple into a
nonconstant unit.

For a branch line `g_i=0`, the homogeneous cubic pullback has the form

```text
g_i(A,C)=ell_i^2 m_i,                                              (17)
```

where `m_i` is a different line.  Distinct branch values ensure that no
`m_i` is one of the deleted ramification lines.  On `U`, `ell_i` is a unit,
so `(17)` becomes precisely the one reduced principal visible companion
required by THM-3801.  Thus companions survive; constant units do not.

## 4. A rational four-branch packet

The following integral packet realizes every positive condition above with
small exact coefficients:

```text
(a,b,c,d)=(A,C,7A,-3A),                                           (18)

omega^2=-7A^2+C omega-A theta,
omega theta=3A^2,
theta^2=3AC-3A omega-7A theta.                                    (19)
```

Its index and discriminant are

```text
I=-(A x^3+C x^2y+7Axy^2-3Ay^3),                                  (20)

Delta=A(C+5A)(4C+19A)(3C-17A).                                   (21)
```

There is an explicit isomorphism to the third Veronese:

```text
A=s^2t,                         omega=s^3,
C=s^3+7st^2+3t^3,              theta=3st^2.                       (22)
```

The four cubics in `(22)` span `k[s,t]_3`, and direct substitution verifies
all of `(19)`.  On the projective coordinate `z=s/t`,

```text
C/A=(z^3+7z+3)/z^2,                                                (23)
```

whose critical points and values are

```text
z=0,-1,-2,3,                 C/A=infinity,-5,-19/4,17/3.          (24)
```

The complete ramification/companion packet is

```text
A             =s^2 t,
C+5A          =(s+t)^2(s+3t),
4C+19A        =(s+2t)^2(4s+3t),
3C-17A        =(s-3t)^2(3s+t).                                   (25)
```

Thus the ramification arrangement is

```text
H=s(s+t)(s+2t)(s-3t).                                             (26)
```

On `U=D(H^3)`, a concrete basis of the lattice `(7)` is

```text
omega=s^3,
(omega+A)/omega=(s+t)/s,
(omega+2A)/omega=(s+2t)/s,
(omega-3A)/omega=(s-3t)/s.                                      (27)
```

Every item in `(27)` is a nonconstant unit.  The first is the radial unit;
the last three are the punctured-projective ratios.  Hence `(18)` is the
first clean positive cubic-normalization grammar isolated by THM-3801, but
it is also a decisive hostile to using homogeneous linear coefficients for
a planar Jacobian counterexample.  A survivor must break the cone grading
(or otherwise prevent deleted ramification classes from producing units).

The exact companion named in the metadata verifies `(3)--(9)`,
`(18)--(27)`, the degree-three field law, and the four-dimensional lattice
index.  Normal and optimized executions byte-match the frozen transcript.
**QED, conditional only on independent hostile audit.**
