---
id: THM-2842
title: "Ordered positive-cone Vandermonde multiplier observability"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Laguerre
  multipliers are exactly radial-variance jets.  On any
  k ordered nonzero positive adjacent-difference cones, the first k such
  jets form a strictly positive generalized Vandermonde decoder.  In the
  optimal D-dimensional Laguerre quotient, D-1 extra observations are
  sharply necessary in general: ell_(D-1) is a trace-zero unit invisible
  through the first D-1 scalar Krylov readouts.  A single selector works
  under the prior unit promise only in cardinal-idempotent form.  These
  jets remain external data and are not implied by ordinary Gaussian
  nullity.
source: root/ordered-cone-vandermonde-multiplier-access-2026-07-28
depends_on:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
related:
  - THM-2638-radial-height-graded-wick-decoder-and-laplace-forgetting-boundary
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - THM-2845-local-residue-versus-split-trace-unit-observability
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
script: 04-computation/gmc_laguerre_variance_jet_observability_thm2842.py
output: 05-knowledge/results/gmc_laguerre_variance_jet_observability_thm2842.out
script_sha256: a3e63877c3ad1162eeab5a30617d71022b2ff6ba497118fccd8bf40b730d610c
output_sha256: ac1729d28a1f19c098815c82269f4b267c3df2353bd1d1e2ebce7d7a5ecef290
hash_basis: LF-normalized bytes
---

# THM-2842 -- ordered positive-cone multiplier observability

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(s^n)=n!,                    f_n=s^n/n!,
d_n=f_(n+1)-f_n,                                      (1)
```

and let

```text
ell_r(s)=(-1)^r r! L_r(s)                             (2)
```

be the monic Laguerre polynomial.  Define the variance-deformed
factorial functional

```text
L_t(f)=t^(-1) integral_0^infinity f(s)e^(-s/t) ds,
t>0.                                                  (3)
```

The theorem has two complementary parts:

1. on ordered positive adjacent-difference cones, the first `k`
   variance jets give a strictly positive `k` by `k` decoder; and
2. on the unrestricted optimal Laguerre carrier, the full `D-1` extra
   observations are sharply necessary, even when the unknown class is
   nonzero at every positive Gauss--Laguerre node.

## 1. Variance jets are Laguerre multipliers

For every polynomial `f` and every `r>=0`,

```text
(partial_t^r L_t(f))|_(t=1)=L(ell_r f).                (4)
```

It is enough to check `f=s^n`.  The left side is

```text
partial_t^r(n!t^n)|_(t=1)=n! n^(falling r).            (5)
```

The standard Laguerre expansion gives the same value on the right.
Equivalently, `(4)` is the finite-difference identity

```text
sum_(j=0)^r
 (-1)^(r+j) r! binom(r,j)(n+j)!/j!
 =n! n^(falling r).                                   (6)
```

In particular,

```text
L_t(d_i)=t^(i+1)-t^i,
L(ell_r d_i)=r i^(falling (r-1)),          r>=1.       (7)
```

Thus the multiplier observations absent from scalar Gaussian nullity
have a precise physical meaning: they are radial-variance derivatives.

## 2. The ordered-cone Vandermonde theorem

Let `I_1,...,I_k` be finite nonempty subsets of nonnegative integers,
strictly ordered in the sense

```text
I_1<I_2<...<I_k

iff every i in I_a is smaller than every j in I_(a+1). (8)
```

For positive cone elements

```text
U_a=sum_(i in I_a) lambda_(a,i)d_i,
lambda_(a,i)>=0,             U_a!=0,                   (9)
```

form the jet matrix

```text
J_(r,a)=L(ell_r U_a),                 1<=r,a<=k.        (10)
```

Then

```text
det J
 =k! sum_(i_a in I_a)
      (product_a lambda_(a,i_a))
      product_(a<c)(i_c-i_a)
 >0.                                                    (11)
```

Indeed, for singleton labels `i_1<...<i_k`, `(7)` gives

```text
det_[r,a](r i_a^(falling (r-1)))
 =k! product_(a<c)(i_c-i_a),                           (12)
```

because the falling factorials are a monic polynomial basis.
Multilinearity in the columns gives `(11)`, and every surviving summand
is positive.

Consequently, for a known ordered family `U_1,...,U_k`, the observations

```text
L(ell_1 H),...,L(ell_k H)                              (13)
```

recover every coefficient in

```text
H=alpha_1U_1+...+alpha_kU_k.                           (14)
```

Every adjacent-difference combination has `L(H)=0`.  Therefore the monic
triangular change of basis between `1,ell_1,...,ell_k` and
`1,s,...,s^k` shows that the inserted bank

```text
L(sH),...,L(s^kH)                                     (15)
```

is equivalent.  Exactly `k` scalar linear observations are
dimension-minimal on a `k`-dimensional complex span.

For the two cones of THM-2830, write

```text
Lambda=sum_i lambda_i,       i_bar=sum_i i lambda_i/Lambda,
M=sum_j mu_j,                j_bar=sum_j j mu_j/M.     (16)
```

Support separation gives `i_bar<j_bar`.  If

```text
H=alpha U+beta V,
a=L(ell_1H),                 c=(1/2)L(ell_2H),          (17)
```

then the decoder is explicit:

```text
alpha=(j_bar a-c)/(Lambda(j_bar-i_bar)),
beta =(c-i_bar a)/(M(j_bar-i_bar)).                    (18)
```

This is rational and cutoff-free once the two cone elements are known.

## 3. The optimal quotient and its trace-zero units

Fix `D>=2` and let

```text
A_D=Q[s]/(ell_D).                                      (19)
```

THM-2815 gives its positive Gauss--Laguerre quotient readout

```text
Lambda_D([f])=sum_(i=1)^D w_i f(x_i)
```

and, for the degree-below-`D` representative, its weighted trace formula

```text
Lambda_D([f])=L(f)
             =Tr_(A_D/Q)(omega_D f),       w_i>0.       (20)
```

where `x_1,...,x_D` are the simple positive roots of `ell_D`.
The classes that are scalar-null but nonzero at every node are exactly

```text
Tr(omega_D f)=0,                 Norm_(A_D/Q)(f)!=0.    (21)
```

In other words, they are the trace-zero units for the weighted trace.
This locus is nonempty over `Q`: the explicit class

```text
f=ell_(D-1) mod ell_D                                  (22)
```

has zero readout by Laguerre orthogonality and is a unit because
successive Laguerre polynomials are coprime.

## 4. Sharp Krylov and Pade boundary

For the degree-below-`D` representative of `f`, put `u=[s] in A_D` and

```text
M_f^D(z)=sum_(r>=0)Lambda_D(u^r[f])z^r,
Q_D(z)=z^D ell_D(1/z)=product_i(1-x_i z).              (23)
```

Partial fractions in `(20)` give a unique numerator `N_f`, of degree
below `D`, such that

```text
M_f^D(z)=N_f(z)/Q_D(z).                                (24)
```

This is a quotient Krylov series, not the infinite original factorial
series.  Gauss--Laguerre exactness gives

```text
Lambda_D(u^r[f])=L(s^r f),                  0<=r<D,     (24a)
```

which is precisely the finite bank used below; beyond the exact horizon
the two sequences may differ.

The constant term and the simple pole residues give the exact
characterization

```text
L(f)=0                    iff z divides N_f,
f(x_i)!=0 for every i     iff gcd(N_f,Q_D)=1.          (25)
```

The factorial Hankel matrix `[L(s^(r+j))]_(0<=r,j<D)` is invertible.
Therefore every nonzero scalar-null class is exposed by at least one of

```text
L(sf),...,L(s^(D-1)f).                                 (26)
```

This length is sharp even under the stronger unit promise.  For `(22)`,

```text
L(s^r ell_(D-1))=0,                         0<=r<D-1,
L(s^(D-1)ell_(D-1))=((D-1)!)^2,                       (27)

N_(ell_(D-1))(z)=((D-1)!)^2 z^(D-1).                  (28)
```

Thus a class can be nonzero at every positive node while remaining
invisible for `D-1` consecutive scalar Krylov readouts.

## 5. The unique shape of one universal unit selector

After extending scalars to `C`, suppose the unknown class is promised to
be a trace-zero unit.  A single multiplier

```text
g in A_D tensor_Q C
```

detects every such class,

```text
Lambda_D(gf)!=0
 whenever Lambda_D(f)=0 and product_i f(x_i)!=0,       (29)
```

if and only if

```text
g=c+d lambda_i,                  d!=0,                 (30)
```

for some Gauss node `x_i` and some `c,d in C`, where `lambda_i` is its
Christoffel cardinal idempotent over the splitting field:

```text
lambda_i(x_j)=delta_(ij).                              (31)
```

The forward implication is linear algebra.  In node coordinates, let

```text
a=(w_1,...,w_D),             b=(w_1g(x_1),...,w_Dg(x_D)).
```

The common kernel of `a` and `b` contains no torus point precisely when
it lies in a coordinate hyperplane.  Over the infinite field `C`, this
is equivalent to some coordinate vector lying in `span(a,b)`, hence to
`g` being constant off one node, which is `(30)`.  Conversely, for
trace-zero `f`,

```text
Lambda_D((c+d lambda_i)f)=d w_i f(x_i)!=0.             (32)
```

Every nonconstant selector `(30)` has degree exactly `D-1` and depends
on the chosen cutoff/node.  Without the unit promise, the scalar-null
subspace has dimension `D-1`, so `D-1` additional linear observations
are necessary and `(26)` is optimal.

## 6. Sharp endogenous and homogeneous boundaries

The external bank `(13)` must not be mistaken for ordinary moment data.
THM-2830 obtains a different, endogenous cubic detector.  It is sharp:
with

```text
U=d_1,              V=d_2,
alpha=(-3+i sqrt(3))/2,
H=alpha U+V,                                           (33)
```

one has

```text
L(H)=L(H^2)=0,
L(H^3)=36+54i sqrt(3)!=0.                              (34)
```

Both rays are divisible by `s`.  Therefore, with

```text
h=H/s,                 P=W+Z h(ZW),        W=conj(Z),
```

for a standard complex Gaussian `Z`, charge balance shows that moments
one through five vanish and the sixth equals

```text
720+1080i sqrt(3).                                     (35)
```

The variance-jet route also respects the single-height boundary.  If a
radial cancellation is homogeneous of total height `A`, then

```text
L_t(G)=t^A L_1(G),                                     (36)
```

so a zero at `t=1` survives every variance derivative.  Variance access
is a genuine multiplier sidecar, not a universal bypass of homogeneous
cancellation, MISTAKE-211, arbitrary SFC, or HYP-8765.

## 7. Exact companion

The exact companion verifies:

1. the variance-jet identity and the adjacent-difference response;
2. point and positive-mixture Vandermonde determinants;
3. the monomial/Laguerre triangular basis change;
4. the sharp `ell_(D-1)` Krylov delay, quotient Pade numerator, and
   consecutive Laguerre coprimality through `D=12`;
5. root-free weighted traces, norms, and resultants through `D=8`;
6. finite-horizon agreement followed by the first original/quotient
   alias divergence;
7. `1,089` finite node-coordinate selector patterns, including all
   positive cardinal-selector cases; and
8. the sharp complex/Gaussian witness `(33)--(35)`.

All truth-bearing gates are explicit exceptions and all arithmetic is
integer or rational.  Reproduce with

```text
python 04-computation/gmc_laguerre_variance_jet_observability_thm2842.py
python -O 04-computation/gmc_laguerre_variance_jet_observability_thm2842.py
```

Both modes byte-match the stored transcript.

## 8. Independent hostile audit

An independent audit rederived every algebraic identity and caught one
load-bearing quotient distinction before promotion: the rational function
in `(24)` belongs to the quotient Krylov sequence, not to the infinite
original factorial sequence.  After the repair, the audit checked the
finite exact horizon and first alias, the selector iff over `C`, the
trace-zero unit and Pade boundaries, zero labels, the two-cone decoder,
and the complex moment-six witness.  It independently replayed normal and
optimized companions against the stored transcript and verified the
declared LF hashes.

**QED.**
