---
id: THM-3855
title: "Formal inverse-discriminant lift and algebraization gate"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  At the THM-3808 rational four-ray packet, the four coefficient
  gradients of the binary-cubic discriminant form a basis of all homogeneous
  cubics.  Hence every discriminant deformation beginning in total degree
  five has an all-orders formal coefficient lift with fixed linear part.
  The one-place THM-3853 deformation therefore has a connected normal
  nonmonogenic formal S3 cubic completion.  Polynomial algebraization,
  constant units on a global etale open, and a Keller atlas remain open.
source: root / inverse binary-cubic discriminant lane, 2026-08-23
audit: >
  SELF-HOSTILE EXACT CANDIDATE.  The proof computes the complete coefficient
  gradient matrix, whose determinant is 640000, identifies its ideal with
  `(A,C)^3`, and gives a homogeneous right-inverse recursion.  The companion
  replays the one-place target through total degree twelve, verifies the
  explicit first quadratic correction and a nonzero degree-thirteen residual
  for that truncation, and contains no Python asserts.  Normal and optimized
  runs byte-match the frozen transcript.  Independent hostile audit remains.
depends_on:
  - THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap
related:
  - THM-3853-quadratic-depth-inverse-discriminant-one-place-gluing-obstruction
script: 04-computation/jc2_formal_inverse_discriminant_lift_thm3855.py
output: 05-knowledge/results/jc2_formal_inverse_discriminant_lift_thm3855.out
script_sha256: 1634b5a16fc5addcde2245f378ea7b7c7e65ea2e297a8b70aef9989a57c28b65
output_sha256: 28b97d014ee1ca4bb4df873c0dd064b37c3dbb7b681184ef7af26d6ba4dd585a
semantic_sha256: 264a006df84eaac70a48a6a6836809ad52a29d772be4e312703e910bb370eb46
hash_basis: raw LF bytes
---

# THM-3855 -- the one-place inverse discriminant has no formal obstruction

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**  Work over an algebraically closed field `k` of
characteristic zero.  Put

```text
R=k[A,C],               m=(A,C),               Rhat=k[[A,C]].   (1)
```

For a binary cubic coefficient row `f=(a,b,c,d)`, write

```text
Disc(f)=b^2c^2-4ac^3-4b^3d-27a^2d^2+18abcd.                    (2)
```

At the THM-3808 row

```text
f_1=(A,C,7A,-3A),                                                (3)
```

one has

```text
Delta_0=Disc(f_1)
       =A(C+5A)(4C+19A)(3C-17A).                                (4)
```

The inverse-discriminant map at `(3)` is formally surjective with a
three-degree shift:

> For every `Phi in m^5 Rhat`, there is an
> `eta=(eta_a,eta_b,eta_c,eta_d) in (m^2 Rhat)^4` such that
>
> ```text
> Disc(f_1+eta)=Delta_0+Phi.                                    (5)
> ```

In particular, for every `lambda in k*`, the irreducible one-place target

```text
delta_lambda=Delta_0+lambda C^5                                 (6)
```

from THM-3853 has an exact formal binary-cubic lift with fixed linear part.
The corresponding Delone--Faddeev algebra over `Rhat` is a connected normal
rank-three domain, is globally nonmonogenic over `Rhat`, and has generic
Galois closure `S3`.

For completeness, the one-place assertion in `(6)` is visible without using
the inverse-discriminant calculation.  On `C!=0`, put `t=A/C` and

```text
F(t)=t(1+5t)(4+19t)(3-17t).                                    (6a)
```

Then `(6)` is `C^4(F(t)+lambda C)=0` on this chart, and the curve has the
polynomial normalization

```text
C=-F(t)/lambda,                    A=-tF(t)/lambda.             (6b)
```

The ratio `A/C=t` is a rational inverse.  The four distinct roots of `F`
map to the origin, while every finite `t` remains affine and the polynomial
parametrization has only `t=infinity` on its projective completion.  Thus
`delta_lambda` is irreducible with normalization `A1` and one place at
infinity.  Its degree-four tangent cone is exactly the four distinct lines
in `(4)`.

This does **not** produce a polynomial cubic over `k[A,C]`.  THM-3853 proves
that a lift consisting only of `(3)` plus homogeneous quadratic corrections
does not exist.  Corrections of degree at least three, finite termination,
the unit group of a global deleted-ramification open, and any planar Keller
map remain open.

## 1. The coefficient gradients span every cubic

Differentiate `(2)` with respect to `(a,b,c,d)` and evaluate at `(3)`.  In
the ordered monomial basis

```text
(A^3,A^2C,AC^2,C^3),                                            (7)
```

the four resulting homogeneous cubics are the columns of

```text
M=[[-1858, -378, -588,  162],
   [ -378,   98,  -54,  126],
   [    0,   36,   14,    0],
   [    0,    0,    0,   -4]].                                  (8)
```

Directly,

```text
det(M)=640000 !=0.                                               (9)
```

Thus the coefficient gradients form a basis of `R_3`, and their homogeneous
ideal is exactly

```text
(Disc_a(f_1),Disc_b(f_1),Disc_c(f_1),Disc_d(f_1))=m^3.          (10)
```

This is stronger than a single favorable tangent direction: multiplication
by arbitrary forms gives a surjection

```text
(R_(n-3))^4 -> R_n,       u |-> sum_i Disc_i(f_1)u_i            (11)
```

for every `n>=3`.

## 2. Homogeneous recursion proves the formal lift

Write `Phi=sum_(n>=5) Phi_n` by total degree.  Suppose corrections have been
chosen so that `(5)` is true through degree `n-1`, and let `E_n in R_n` be
the remaining degree-`n` error.  By `(11)`, choose a homogeneous coefficient
correction

```text
u_(n-3) in (R_(n-3))^4,
sum_i Disc_i(f_1)(u_(n-3))_i=-E_n.                              (12)
```

The degree-`n` change in the discriminant is exactly the left side of
`(12)`.  Indeed, all previous corrections have degree at least two, so
replacing the derivative at `f_1` by the derivative at the current row costs
at least one extra total degree; terms quadratic in `u_(n-3)` do likewise.
The new row therefore solves `(5)` through degree `n` without disturbing
lower degrees.

Starting at `n=5` and using completeness of `Rhat` gives a convergent
`m`-adic coefficient row `f_1+eta` satisfying `(5)` exactly.  No analytic or
finite-degree convergence is asserted.

For the target `(6)`, the inverse of `(8)` gives the particularly concrete
first correction

```text
eta_2=lambda C^2(
  91449/40000,
 151263/40000,
-194481/20000,
          -1/4).                                                (13)
```

It creates exactly `lambda C^5` at degree five.  The deterministic companion
uses a fixed monomial right inverse of `(8)` to continue through degree
twelve.  The degree-thirteen error of that finite truncation is nonzero, as
it should be: only the infinite recursion has been proved exact.

## 3. The formal cubic passes the index, normality, and S3 gates

The binary index form of the Delone--Faddeev algebra is

```text
I(x,y)=-(a x^3+b x^2y+c xy^2+d y^3).                            (14)
```

Every coefficient of `f_1+eta` belongs to `m`, because `(3)` is linear and
`eta` starts in `m^2`.  Hence

```text
I(Rhat^2) subset m,                                             (15)
```

so the index form represents no unit and the formal cubic algebra is not
monogenic.

For `(6)`, the tangent cone `(4)` is a product of four distinct lines.
Consequently `delta_lambda` is reduced in `Rhat`, with four smooth formal
branches of multiplicity one.  The Delone--Faddeev algebra is finite free,
hence `S2`, and is generically etale.  Away from `(6)` it is etale in
codimension one; at each branch the discriminant valuation is one, so the
discriminant-index formula forces height-one index zero in the integral
closure.  It is therefore normal.

Modulo `m`, all four coefficients vanish and the special fibre is

```text
k[omega,theta]/(omega,theta)^2,                                 (16)
```

which is local.  Since `Rhat` is henselian, the finite algebra has no
nontrivial idempotent.  A connected normal finite algebra is a domain.
Finally, `(6)` is not a square in `Frac(Rhat)` because each of its four
height-one factors has odd valuation.  The connected generic cubic therefore
has nonsquare discriminant and Galois closure `S3`.

## 4. Exact boundary and next construction step

The formal packet now passes four gates simultaneously:

```text
one-place global discriminant target,
connected normal cubic order,
nonmonogenic index form,
generic S3 monodromy.                                           (17)
```

The missing implication is algebraization.  A finite polynomial row would
have to realize `(6)` exactly while retaining `(15)`; THM-3853 excludes only
maximum coefficient degree two.  Degree-three and interacting higher-degree
rows are the first live search.  Even after such a row is found, its maximal
etale open must still have only scalar units and admit the required plane
atlas.  Thus `(5)` is a positive local existence theorem and a sharp
algebraization gate, not a Jacobian counterexample.

Reproduction:

```bash
python3 04-computation/jc2_formal_inverse_discriminant_lift_thm3855.py
python3 -O 04-computation/jc2_formal_inverse_discriminant_lift_thm3855.py
```

Both modes byte-match the frozen 60-gate transcript.  **QED candidate,
pending independent hostile audit.**
