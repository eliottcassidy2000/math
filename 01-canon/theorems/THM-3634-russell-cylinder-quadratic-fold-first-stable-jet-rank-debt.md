---
id: THM-3634
title: "Russell-cylinder quadratic-fold first-stable-jet rank debt"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  retained-collision polynomial Q on
  the quadratic fold q=Q(x)+t^2, a hypothetical global regular target pair
  with nonzero constant source Jacobian has two independent zero-stable
  horizontal restrictions modulo constants and two independent first-stable
  sidecars.  This closes every pair with either rank at most one.  A sharp
  polynomial hostile exists only in the enlarged equal-triple-value ring;
  it is not claimed to lie in the actual target restriction ring.  The
  fully mixed rank-two/rank-two case remains OPEN.
source: root / audit_thm3629 THM-3630 fixed-polynomial boundary, 2026-08-21
audit: >
  PASS -- an independent hostile audit rederived the quadratic chain-rule
  determinant and both rank-two necessity arguments, independently expanded
  the enlarged-ring hostile, and verified byte-identical normal, optimized,
  and stored 28-gate transcripts plus all static gates.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
related:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3632-russell-cylinder-formal-pair-algebraization-triple-fibre-obstruction
script: 04-computation/jc2_russell_cylinder_quadratic_fold_first_stable_jet_rank_debt_thm3634.py
output: 05-knowledge/results/jc2_russell_cylinder_quadratic_fold_first_stable_jet_rank_debt_thm3634.out
script_sha256: bbb0930e99cdc66bc36d70de71c1f9fe5488e590677679f4faa1ad932c3a488f
output_sha256: 716e049f3aef6bd436585b36ce9b25b7967358be4fb1f8ff26ef37c07bbdbc8e
hash_basis: raw LF bytes
---

# THM-3634 -- Russell-cylinder quadratic-fold first-stable-jet rank debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
extracts an unbounded-support
necessary condition for a global pair on a quadratic stable fold.  It is a
global polynomial restriction argument, not a finite arbitrary-two-form
search.

All rings and closed points are over `C`.

## 0. Setup and statement

Use the exponent-two Danielewski surface and compiler

```text
Y_2=Spec R,       R=C[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3).       (1)
```

Let `Q in C[x]` be arbitrary, with no parity hypothesis, subject only to

```text
Q(-1)=Q(1)=-3,                    Q(0)=-3/4.            (2)
```

Take the quadratic fold before the Russell-cylinder isomorphism,

```text
E_Q(x,t)=(x,Q(x)+t^2,t),                              (3)
```

where the last target coordinate is `w=t`.  Because the Russell map is a
polynomial isomorphism, every pair of global regular functions on the
exponent-one target transports uniquely to

```text
F^#,G^# in R[w].                                       (4)
```

Write their stable expansions

```text
F^#=F_0+wF_1+w^2 F_(>=2),
G^#=G_0+wG_1+w^2 G_(>=2),             F_i,G_i in R.    (5)
```

Define the exact restriction homomorphism along the zero-stable source
curve

```text
gamma_Q:R -> C[x],
gamma_Q(H)=H(b(x,Q(x)),c(x,Q(x)),e(x,Q(x))).            (6)
```

Put

```text
U=gamma_Q(F_0),       V=gamma_Q(G_0),
A=gamma_Q(F_1),       B=gamma_Q(G_1).                   (7)
```

The theorem says that if the source pullbacks `f,g` of `F^#,G^#`
satisfy

```text
Jac_(x,t)(f,g)=kappa in C*,                             (8)
```

then necessarily

```text
dim_C span{[U],[V]}=2       in C[x]/C,                  (9)

dim_C span{A,B}=2           in C[x].                   (10)
```

Thus a survivor needs two independent zero-stable horizontal classes and
two independent first-stable sidecars.  The assertion holds for every `Q`
in `(2)`, even or non-even.

## 1. The exact first-stable determinant

The criticality of the quadratic fold is load-bearing:

```text
partial_t(Q(x)+t^2)|_(t=0)=0.                           (11)
```

Consequently the chain rule applied to `(5)` gives, on `t=0`,

```text
f_x=U',       f_t=A,       g_x=V',       g_t=B.         (12)
```

Equation `(8)` therefore has the exact zero-stable shadow

```text
                         U'B-AV'=kappa.                 (13)
```

No bounded target-degree assumption enters `(13)`.

At the three retained source points, direct substitution in `(1),(2)` gives

```text
(x,q)=(-1,-3),(0,-3/4),(1,-3)
             |-> (b,c,e)=(0,0,-3).                     (14)
```

Every polynomial in the image of `gamma_Q` therefore has one common value
at `x=-1,0,1`.  In particular,

```text
U(-1)=U(0)=U(1),             V(-1)=V(0)=V(1),
A(-1)=A(0)=A(1),             B(-1)=B(0)=B(1).           (15)
```

## 2. Zero-stable horizontal rank must be two

Suppose `(9)` failed.  Then there are `alpha,beta`, not both zero, and a
constant `delta` such that

```text
alpha U+beta V=delta.                                  (16)
```

Complete `(alpha,beta)` to a matrix `M in GL_2(C)` and apply that constant
linear change to the two target outputs.  Call the new source pullbacks
`h,k`, with `h=alpha f+beta g`.  Their Jacobian is the nonzero constant
`det(M)kappa`, while `(12),(16)` give

```text
h_x(x,0)=0.                                             (17)
```

Restricting the new Jacobian identity to `t=0` gives

```text
-h_t(x,0) k_x(x,0)=det(M)kappa.                         (18)
```

The two factors lie in `C[x]`.  A product equal to a nonzero constant forces
both to be polynomial units.  Hence

```text
k_x(x,0)=mu in C*,             k(x,0)=mu x+nu.          (19)
```

But `k(x,0)` is another constant linear combination of `U,V`, so `(15)`
forces equal values at `-1,0,1`.  The affine polynomial in `(19)` cannot do
that.  This contradiction proves `(9)`.

## 3. First-stable sidecar rank must be two

Suppose instead that `(10)` failed.  Then some nonzero constant row
`(alpha,beta)` satisfies

```text
alpha A+beta B=0.                                      (20)
```

Complete this row to `M in GL_2(C)` as before and let `h,k` be the resulting
output pair, with `h=alpha f+beta g`.  Equations `(12),(20)` give

```text
h_t(x,0)=0.                                             (21)
```

The zero-stable Jacobian identity is now

```text
h_x(x,0) k_t(x,0)=det(M)kappa.                          (22)
```

Again both factors are polynomial units.  Thus

```text
h_x(x,0)=mu in C*,             h(x,0)=mu x+nu,          (23)
```

contradicting the common triple value inherited from `(15)`.  This proves
`(10)`.

The conclusion is stronger than excluding a surface-only or stable-only
output.  It also excludes, for example,

- either output whose explicit stable dependence begins at order `w^2`;
- first-stable coefficients with one constant projective direction; and
- arbitrary surface support combined only with two constant affine
  `w`-coefficients.

## 4. Sharp hostile in the equal-triple-value enlargement

Put

```text
L=x(x^2-1).                                             (24)
```

The full ring of polynomials with one common value at `-1,0,1` is

```text
E=C+L C[x].                                             (25)
```

Indeed subtract the common value and use the three distinct roots.  For
every retained `Q`, equation `(15)` gives only the inclusion

```text
                         gamma_Q(R) subset E.           (26)
```

Inside this **enlarged test universe**, define

```text
U=(L/2)(3x^2-2),                  V=L,

A=L(225x^3/8-75x/4),             B=1+45xL/4.            (27)
```

They obey

```text
(U,V,A,B)|_(x=-1,0,1)=(0,0,0,1),

U'B-AV'=1,                                               (28)
```

and both ranks in `(9),(10)` are exactly two.  Therefore common triple
values plus the one-variable determinant identity do not yield a stronger
universal obstruction.

This is a hostile control, **not a target witness**.  No polynomial `Q` and
no elements of `R` mapping to the four polynomials in `(27)` are produced.
In particular, membership of `(27)` in any actual `gamma_Q(R)` is not proved.

## 5. Scope and OPEN boundary

What this theorem proves is a necessary rank debt for **global regular
decomposable pairs** on the quadratic fold.  It does not apply to a freely
chosen arbitrary target two-form, whose coefficients need not arise from
two functions or satisfy the stable-expansion relations `(5),(12)`.

The formal multi-germ pair of THM-3630 is also not contradicted.  Its three
source completions permit branchwise integration constants and are not the
restriction of one polynomial on the connected source plane.  Finite
Hermite approximation cannot detect that all-order globality debt.

The rank-two/rank-two case remains **OPEN**.  This theorem supplies
neither

- a no-pair theorem for arbitrary fully mixed target functions;
- one fixed non-even polynomial surviving all orders;
- a global target pair with constant source Jacobian; nor
- a Keller map or counterexample to `JC(2)`.

Conversely, any positive global pair in the remaining cell would already
give a polynomial Keller map `(f,g):A2->A2` identifying the three points in
`(14)`.  It would be a genuine `JC(2)` counterexample, not merely a formal or
arbitrary-two-form survivor.

## 6. Exact verification

The companion verifies the compiler relation and retained triple, the
quadratic critical derivative and determinant shadow, constant-output
Jacobian covariance, both unit-factor contradiction controls, a direct
equal-value enlargement control, every hostile value and factorization in
`(27),(28)`, both rank claims, and an active hostile mutation.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_quadratic_fold_first_stable_jet_rank_debt_thm3634.py
python3 -O 04-computation/jc2_russell_cylinder_quadratic_fold_first_stable_jet_rank_debt_thm3634.py
```

Both streams must be byte-identical to

```text
05-knowledge/results/jc2_russell_cylinder_quadratic_fold_first_stable_jet_rank_debt_thm3634.out
```

The frozen companion reports zero Python assertion statements.  Its exact
gates support the provisional proof but do not replace the pending
independent hostile audit.
