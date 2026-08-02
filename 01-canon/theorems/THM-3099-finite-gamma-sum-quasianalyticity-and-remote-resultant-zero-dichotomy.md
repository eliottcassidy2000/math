---
id: THM-3099
title: "Finite Gamma-sum quasianalyticity and remote-resultant zero dichotomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A finite
  real linear combination of eventually positive hypergeometric terms whose
  consecutive quotients are rational functions is either the exact zero
  sequence or has a fixed nonzero sign eventually.  Every one-parameter
  remote factorial physical resultant belongs to this class.  Consequently
  one exact good terminal extension of a fixed support rules out persistent
  identity and makes every sufficiently remote terminal extension good.
  The sign of the test value need not be the eventual sign.
source: multiscale-newton-flag-2026-08-01
audit: >
  An independent hostile audit rederived the ordered asymptotic scale,
  rational-shift evaluation determinant, polynomial inverse bound, and
  exclusion of all-orders flat nonidentities after choosing a basis from the
  original positive terms.  It checked that determinant-one elimination and
  full Newton expansion leave only finite products with eventually positive
  rational quotients, and independently accepted the finite-core escape
  corollary.  Normal, optimized, and stored output matched exactly; both LF
  hashes passed and no truth-bearing assert was present.  The dependency on
  THM-3097 and the requirement that an affine-clock test sample lie on that
  clock are explicit.
depends_on:
  - THM-3097-translated-support-monge-compactification-and-cofinite-bad-set-induction
related:
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3097-translated-support-monge-compactification-and-cofinite-bad-set-induction
script: 04-computation/gmc_finite_gamma_sum_quasianalyticity_thm3099.py
output: 05-knowledge/results/gmc_finite_gamma_sum_quasianalyticity_thm3099.out
script_sha256: 4f1417aaea8caeb1a99d6d9c3ba91588493e02a92ac3de6c941610041cb223b3
output_sha256: d93294df029997fe2b2119f1e197b5928e153e3e4d391bb00dc4cebc433cfdb7
hash_basis: LF-normalized bytes
---

# THM-3099 -- finite Gamma-sum quasianalyticity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The remote-terminal theorems identify a leading positive face when a lower
resultant is already nonzero.  A bad child has no such face, but its exact
remote resultant still has much more rigidity than an arbitrary sequence:
it is a finite sum of factorial products, and every summand has a positive
rational consecutive quotient.  This excludes infinitely recurring
sporadic cancellation.  What remains is an exact persistent identity.

The distinction is a useful holotopy boundary.  A single exact good fibre
does not determine the orientation at infinity, but it does prove that the
whole one-parameter section is not the zero section.

## 1. Positive rational-ratio terms

Fix an integer `n_0`.  Call a real sequence `T(n)`, `n>=n_0`, a **proper
positive hypergeometric term** if

```text
T(n)>0,
Q_T(n):=T(n+1)/T(n) belongs to R(n),
Q_T(n)>0                                                (1)
```

throughout the tail.  Removing finitely many initial indices always makes
the last two inequalities lawful, so the precise starting index is not
important.

Let

```text
F(n)=sum_(j=1)^J c_j T_j(n),             c_j in R.       (2)
```

Then exactly one of the following holds:

1. `F(n)=0` for every index in the common admissible tail `n>=n_0`; or
2. there are `N` and `epsilon in {+1,-1}` such that

   ```text
   epsilon F(n)>0                    for every n>=N.     (3)
   ```

In particular, a nonzero finite positive-ratio hypergeometric sum has only
finitely many integer zeros.

The first alternative is an **exact constant linear dependence** among the
sequences, not merely equality of their formal Stirling series.

## 2. The asymptotic scale

Write a rational quotient at infinity as

```text
Q_T(n)=Lambda n^delta
       (1+q_1/n+q_2/n^2+...),
Lambda>0,        delta in Z.                             (4)
```

Factoring its numerator and denominator, or applying Euler--Maclaurin to
`sum log Q_T(k)`, gives for every `M`

```text
log T(n)
 =delta log Gamma(n)+n log Lambda+beta log n+c
  +sum_(r=1)^M b_r/n^r+O(n^(-M-1)).                    (5)
```

All constants are real because `T` is eventually positive.  Formula `(5)`
can also be read directly from a Gamma factorization of the rational
function in `(1)`.

Order the finitely many summands first by `delta`, then by `Lambda`.
Different classes are separated superexponentially or exponentially.
Inside one class, after division by a fixed member, every ratio has a full
Poincare expansion

```text
C n^gamma (1+d_1/n+d_2/n^2+...),          C>0.          (6)
```

The exponents occurring in a finite sum lie in a finite union of sets
`gamma_j-Z_(>=0)`.  This union is reverse well ordered.  Thus a nonzero
formal combination has a first nonzero coefficient and already has a fixed
eventual sign.  The only possible obstruction would be a nonzero constant
combination flat to every order.  Section 3 rules that out exactly.

## 3. Rational-shift determinant excludes formal ghosts

Work inside one `(delta,Lambda)` class and discard exact constant linear
dependencies.  Choose `G_1,...,G_d` as a remaining subset of the original
positive terms that is a basis as sequences on the common tail.  Linear
independence of their evaluation vectors gives fixed
nonnegative shifts

```text
e_1,...,e_d                                               (7)
```

for which one evaluation matrix is nonsingular.  Normalize its columns at
a moving index:

```text
B_(i,j)(n)=G_j(n+e_i)/G_j(n)
          =product_(u=0)^(e_i-1) Q_(G_j)(n+u).           (8)
```

Every entry of `B(n)` is rational.  Its determinant is not the zero rational
function, because it is nonzero at the evaluation index used to choose
`(7)`.  Hence

```text
||B(n)^(-1)||=O(n^K)                                     (9)
```

for some finite `K`, outside finitely many indices.

Suppose `G=sum_j a_jG_j` were flat relative to one member of the class.
Every shifted value `G(n+e_i)` would still be flat, since a fixed shift costs
only a rational quotient.  But `(8)` gives

```text
B(n) diag(G_1(n),...,G_d(n)) a
 =(G(n+e_1),...,G(n+e_d))^T.                            (10)
```

Equations `(6),(9),(10)` would make every
`a_jG_j(n)/G_1(n)` flat.  A nonzero expression of the form `(6)` is never
flat, so every `a_j=0`, a contradiction.

Therefore a nonidentity combination in the maximal asymptotic class has a
first nonzero Poincare coefficient.  If the whole contribution of that
class is an exact identity, remove it and descend through the finitely many
classes.  This proves the dichotomy `(2)--(3)`.

This proof uses neither total positivity nor a zero theorem for arbitrary
P-recursive sequences.  The load-bearing facts are the positive quotient,
which removes root-of-unity oscillation, and the rational quotient, which
makes the inverse shift determinant only polynomially costly.

## 4. Factorial and Gamma monomials

Every sequence

```text
T(n)=A^n product_nu Gamma(a_nu n+b_nu)^e_nu,
A>0,
a_nu in Z_(>0),       e_nu in Z,                       (11)
```

with no tail pole is of the class `(1)`.  Indeed,

```text
T(n+1)/T(n)
 =A product_nu (a_nu n+b_nu)_(a_nu)^e_nu in R(n),     (12)
```

and it is positive eventually.  Gauss multiplication can alternatively
rewrite `(11)` as a geometric factor times products of unit-slope Gamma
functions.  Thus the theorem applies to every finite real sum of factorial
or fixed-shape product-Gamma monomials.

Exact Gamma multiplication identities merely place summands in alternative
`(1)` of the theorem.  They are not mistaken for asymptotically small
remainders.

## 5. Remote physical resultants

Fix a positive factorial support

```text
Lambda_0={lambda_1<...<lambda_(t-1)}                   (13)
```

and append one integer terminal exponent `C>lambda_(t-1)`.  Use

```text
f_a=s^a/a!,
F_r=L((sum_i x_i f_(lambda_i)+z f_C)^r),
1<=r<=t.                                                (14)
```

Because `F_1=sum_i x_i+z`, its determinant-one physical elimination is
independent of `C`.  Before elimination, and for each `k`-layer contribution
to a coefficient after elimination, the term has the form

```text
constant * (kC+b)!/(C!)^k,               0<=k<=r,      (15)
```

where `b` depends only on the fixed lower occurrences.  The homogeneous
resultant is a fixed polynomial in these coefficients.  Expanding its
finite Newton support expresses the physical resultant as

```text
R_(Lambda_0)(C)
 =sum_mu c_mu product_l (k_(mu,l)C+b_(mu,l))!
                    /(C!)^K_mu.                        (16)
```

Every nonzero summand in `(16)` has a positive rational consecutive
quotient.  Section 1 therefore gives the exact alternative

```text
R_(Lambda_0)(C)=0 for every admissible C,

or

R_(Lambda_0)(C)!=0 for every sufficiently large C.     (17)
```

Consequently:

> **One-sample remote certificate.**  If one exact admissible terminal
> exponent `C_0` makes `Lambda_0 union {C_0}` good, then every sufficiently
> remote terminal exponent makes it good.

The same conclusion holds on any fixed affine clock `C=an+b`, because its
consecutive clock quotient is again rational and eventually positive; the
one exact good sample used to rule out the zero sequence must lie on that
clock.

There is a finite-core consequence for the first unresolved factorial
width.  THM-3097 proves that the minimal bad four-prefix bank `C_4` and the
minimal bad five-prefix bank `C_5` are finite.  Suppose that every
`c in C_4` has at least one exact good one-slot extension.  For each fixed
`c`, `(17)` then leaves only finitely many bad terminal exponents.  Every bad
five-support has first bad prefix in `C_4` or `C_5`; hence

```text
[every c in C_4 has one good append]  =>  B_5 is finite. (18)
```

More generally, if SFC is known through width `b` and every member of
`C_(b+1)` has one good append, then `B_(b+2)` is finite.  This turns an
effective version of THM-3097's finite-prefix bank into a finite exact
certification program.  It does not control two or more free tail
coordinates over a bad child.

This result is complementary to THM-3069.  That theorem supplies an explicit
positive leading face whenever the child resultant is nonzero.  The present
theorem also applies when that face vanishes, but it replaces an effective
positive asymptotic by the dichotomy `exact persistent identity / eventual
nonvanishing`.  It does not identify which branch holds without one exact
good sample or a transverse jet.

## 6. Sharp boundaries

The sign at one sample need not be the sign at infinity:

```text
10*2^n-3^n                                                (19)
```

is positive at small positive `n` and negative eventually, although both
summands have positive constant quotient.  Thus the one-sample corollary is
about nonvanishing, not orientation or positivity.

Eventual positivity of the quotients is also necessary.  If a quotient of
`-1` is admitted, then

```text
1+(-1)^n                                                  (20)
```

has infinitely many zeros without being the zero sequence.  Arbitrary
P-recursive or oscillatory holonomic sequences are therefore outside the
theorem.

For several remote exponents moving independently, `(16)` becomes a
multivariate Gamma sum.  The one-dimensional proof applies along each fixed
affine ray, but gives no uniform orthant threshold and no conclusion on a
non-affine path.  It also does not prove SFC(4), identify the persistent
identity branch, preserve a preselected resultant sign, or improve
THM-3097's bad-prefix exponent without a child-detection sidecar.

## 7. Exact evidence

The exact companion checks positive rational shift quotients for factorial
monomials, exact proportionality identities, nonsingular rational-shift
evaluation matrices, the sign-changing positive-quotient hostile, the
root-of-unity boundary, every coefficient shape `(15)` in a finite bank,
and exact two-slot physical remote resultants.

Run

```text
python 04-computation/gmc_finite_gamma_sum_quasianalyticity_thm3099.py
python -O 04-computation/gmc_finite_gamma_sum_quasianalyticity_thm3099.py
```

Both modes must equal the stored transcript after LF normalization.

**QED.**
