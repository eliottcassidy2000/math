---
id: THM-3104
title: "Multiplicatively infinitely divisible Stieltjes width-three common nullity"
status: >
  PROVED + VERIFIED-EXACT + CERTIFIED-ARB + INDEPENDENTLY HOSTILE-AUDITED.
  There is a unique parameter pair in an explicit positive
  rational box for which the full-support Stieltjes moment law
  exp(a n^(3/2)+b n^2) has a three-slot quadratic/cubic common null line.
  The law is multiplicatively infinitely divisible, every generalized
  Hankel minor is strict, and every higher logarithmic finite difference has
  the product-Gamma alternating sign.  Thus those soft positivity properties
  do not prove width-three detection.  This is not a product-Gamma,
  factorial, or Gaussian counterexample.
source: low-child-flag-extension-2026-08-02
audit: >
  An independent hostile audit rederived the Levy--Khintchine moment law,
  multiplicative infinite divisibility, full-support generalized-Hankel
  strictness, alternating curvature identity, exact width-three numerator
  typing, all four Arb/Taylor face signs, Poincare--Miranda existence,
  Jacobian uniqueness, Gram positivity, and the non-product-Gamma scope.
  Fresh normal and optimized executions byte-match the stored transcript;
  both LF hashes and the documentation checker pass.
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-3051-stieltjes-multiplier-gamma-flow-and-moving-lower-hankel-boundary
  - THM-3100-product-gamma-response-monge-compactification-and-bad-prefix-spectrum
script: 04-computation/gmc_stieltjes_multiplicative_id_width3_nullity_thm3104.py
output: 05-knowledge/results/gmc_stieltjes_multiplicative_id_width3_nullity_thm3104.out
script_sha256: 57a589b5f0e02f60a944796b779b027589b3164cca1e687cb6e8ecb70dfe3362
output_sha256: 35a6912d844c1acd0423a5b18439cb22bc7261be4c8dc9090bd7d64b08350e73
hash_basis: LF-normalized bytes
---

# THM-3104 -- a strict Stieltjes width-three nullity hostile

**PROVED + VERIFIED-EXACT + CERTIFIED-ARB + INDEPENDENTLY
HOSTILE-AUDITED.**

Strict Stieltjes positivity is not the missing width-three mechanism.  The
failure persists inside a two-parameter multiplicative convolution semigroup,
after imposing strict positivity of every generalized Hankel minor and the
entire alternating hierarchy of logarithmic finite differences.

The counterexample is nevertheless outside the factorial and finite
product-Gamma classes.  It therefore marks a proof boundary for THM-3100; it
does not refute that theorem or the open product-Gamma width-three problem.

## 1. The certified moment law

Put

```text
a_0=101128546502785376754010485709/10^31,
b_0=38728340591575303661743380503/(2*10^30),            (1)

h_a=10^(-12),                 h_b=10^(-13),             (2)

R=[a_0-h_a,a_0+h_a] x [b_0-h_b,b_0+h_b].               (3)
```

There is a unique pair `(a_*,b_*) in R` such that the moment sequence

```text
w_n=exp(a_* n^(3/2)+b_* n^2),             n>=0,         (4)
```

has the following property.  If `L(x^n)=w_n`, `f_n=x^n/w_n`, and

```text
U=f_1-f_0,                    V=f_2-f_1,                 (5)
```

then there are complex `(alpha,beta)!=(0,0)` for which

```text
L(alpha U+beta V)=0,
L((alpha U+beta V)^2)=0,
L((alpha U+beta V)^3)=0.                                 (6)
```

Thus the normalized support `{0,1,2}` is bad for this Stieltjes functional.
The two projective solutions in `(6)` are nonreal conjugates; no nonzero real
linear combination can have zero quadratic readout.

Numerically, only as orientation and not as proof,

```text
a_*=0.0101128546502785376754010485709...,
b_*=0.0193641702957876518308716902515....                (7)
```

## 2. A citation-free full-support Stieltjes representation

For `a>0`, take the Levy measure on the negative half-line

```text
nu_a(dz)=(3a/(4 sqrt(pi))) 1_(z<0)|z|^(-5/2) dz.         (8)
```

It is a Levy measure because `z^2 nu_a` is integrable at zero and `nu_a` is
integrable at infinity.  Its first absolute moment is also integrable away
from zero, so compensation by `z` on the whole half-line is lawful.  The
Levy--Khintchine formula gives an infinitely divisible real random variable
`Z_a` with, for every `t>=0`,

```text
log E exp(t Z_a)
 =integral_(-infinity)^0 (exp(tz)-1-tz) nu_a(dz)
 =(3a/(4 sqrt(pi))) Gamma(-3/2)t^(3/2)
 =a t^(3/2),                                               (9)
```

because `Gamma(-3/2)=4 sqrt(pi)/3`.  Convergence in `(9)` is direct: the
integrand is `O(z^2)` at zero and `O(|z|)` at negative infinity.

Let `G_b` be an independent centered Gaussian of variance `2b`, set

```text
Y_(a,b)=Z_a+G_b,                   X_(a,b)=exp(Y_(a,b)). (10)
```

Then

```text
E X_(a,b)^n=exp(a n^(3/2)+b n^2).                       (11)
```

Conditional on `Z_a`, the Gaussian density is positive everywhere.  Hence
`Y_(a,b)` has an everywhere-positive density on `R`, and `X_(a,b)` has an
everywhere-positive density on `(0,infinity)`.

The parameters add under independent multiplication:

```text
X_(a,b) X_(a',b')  has law X_(a+a',b+b').              (12)
```

In particular, every law in `(3)` is multiplicatively infinitely divisible.

## 3. Positivity properties that still do not detect the null line

Full support makes `(w_n)` a strict generalized-Hankel Stieltjes sequence.
Indeed, for strictly increasing nonnegative integer lists `p_i,q_j`, Andreief
gives

```text
det(w_(p_i+q_j))
 =1/k! integral det(x_l^p_i) det(x_l^q_j)
              product_l dmu(x_l)>0.                    (13)
```

On the ordered positive chamber both generalized Vandermonde determinants
have the same strict sign, and full support gives that chamber positive
measure.

There is also an all-order discrete-curvature law.  Write

```text
ell_n=log w_n=a n^(3/2)+b n^2.                         (14)
```

For every `n>=0` and `r>=2`,

```text
(-1)^r Delta^r ell_n>0.                                (15)
```

For the stable part this follows from the exact integral

```text
(-1)^r Delta^r n^(3/2)
 =1/Gamma(-3/2) integral_0^infinity
   exp(-nt)(1-exp(-t))^r t^(-5/2) dt>0.                (16)
```

The Gaussian part contributes `2b>0` at `r=2` and zero afterward.  Equivalently,
if `q_n=w_(n+1)/w_n`, then

```text
(-1)^(k+1) Delta^k log q_n>0,              k>=1.       (17)
```

This is the same alternating logarithmic-curvature direction as a finite
product-Gamma moment law.  It is still insufficient for `(6)`.

## 4. The two exact nullity equations

Use THM-2824's Gram and cubic tensors

```text
g11=L(U^2),       g12=L(UV),        g22=L(V^2),
t111=L(U^3),      t112=L(U^2V),
t122=L(UV^2),     t222=L(V^3).                          (18)
```

The quadratic and cubic forms have a common projective root exactly when

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22=0,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2=0.           (19)
```

For a general normalized moment vector `w_0=1,w_1,...,w_6`, exact symbolic
cancellation gives

```text
I1=N1/(w1^5 w2^3),               I2=N2/(w1^5 w2^4),   (20)
```

where `N1,N2` contain respectively `20,23` integer monomials.  Put

```text
F=N1,                         H=9N2-10N1.              (21)
```

After substituting `(4)`, every monomial has the exact form

```text
c product_(n=1)^6 w_n^e_n
 =c exp(a A_e+b B_e),
A_e=sum_n e_n n sqrt(n),          B_e=sum_n e_n n^2.   (22)
```

Thus every derivative needed below is a finite rigorously evaluable
exponential sum.  No ill-conditioned repeated substitution through `(18)` is
used in the interval proof.

## 5. The rigorous Poincare--Miranda box

At the rational centre `(1)`, the companion evaluates `(22)` at `256`-bit Arb
precision.  On each face it uses the centred second-order enclosure

```text
R_J <=1/2(M_aa h_a^2+2M_ab h_a h_b+M_bb h_b^2),       (23)
```

where each `M` is the Arb supremum of the corresponding second derivative on
the whole box.  The certified conservative rational bounds are

```text
a=a_0-h_a:   -7.13*10^(-16)<F<-4.00*10^(-16),
a=a_0+h_a:    4.00*10^(-16)<F< 7.13*10^(-16),         (24)

b=b_0-h_b:   -3.24*10^(-16)<H<-3.22*10^(-16),
b=b_0+h_b:    3.22*10^(-16)<H< 3.24*10^(-16).         (25)
```

The Taylor remainders in `(24)` and `(25)` are respectively less than
`1.8*10^(-26)` and `1.2*10^(-25)`.  Poincare--Miranda therefore gives a point
of `R` with `F=H=0`, hence `N1=N2=0`.

The same whole-box calculation gives

```text
partial_a F>1/2000,
det partial_(a,b)(F,H)>17/10^7.                        (26)
```

For each fixed `b`, `(24)` and `partial_aF>0` make `F=0` a unique graph
`a=phi(b)`.  Along that graph,

```text
d/db H(phi(b),b)
 =det partial_(a,b)(F,H)/(partial_aF)>0.               (27)
```

The two signs in `(25)` then make the common zero unique in `R`.

Finally, throughout `R`, direct interval evaluation gives

```text
0.0482<g11<0.0483,               0.0526<g22<0.0527,
0.000295<g11 g22-g12^2<0.000297.                       (28)
```

Thus the Gram form is positive definite, so THM-2824's divisibility criterion
is applicable and `(19)` really yields `(6)`.

## 6. Exact failure anatomy and scope

The first false implication is now explicit:

```text
full-support Stieltjes
 + strict generalized-Hankel total positivity
 + multiplicative infinite divisibility
 + every alternating logarithmic curvature inequality

 does not imply width-three quadratic/cubic detection.                 (29)
```

The missing coordinate is not another scalar Hankel minor.  THM-3100's
factorwise rational Gamma response and its marked-permutation coupling retain
structure absent from `(29)`.

This law cannot itself be a finite fixed-shape product-Gamma law.  The latter
has

```text
log w_n=A n log n+O(n),                                  (30)
```

up to a scale term, whereas `(4)` has `b_* n^2` with `b_*>0`.  Therefore this
theorem gives no product-Gamma or factorial bad support, no counterexample to
THM-3100, no Gaussian Moment Conjecture counterexample, and no SFC(3)
counterexample for the exponential/factorial functional.  It instead proves
that any successful width-three extension must use more than the soft moment
cone shared by these laws.

## 7. Reproduction

Run

```text
python 04-computation/gmc_stieltjes_multiplicative_id_width3_nullity_thm3104.py
python -O 04-computation/gmc_stieltjes_multiplicative_id_width3_nullity_thm3104.py
```

Both modes must byte-match the stored transcript after LF normalization.  The
script independently derives `(19)--(22)`, verifies the Levy normalization,
and certifies every face, derivative, Jacobian, and Gram bound quoted above.

**QED.**
