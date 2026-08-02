---
id: THM-3100
title: "Product-Gamma response-Monge compactification and bad-prefix spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every fixed finite product of positive-shape Gamma moments, translating
  an arbitrary fixed-width support sufficiently far makes its first-window
  resultant positive, uniformly over all internal gaps.  Above any fixed
  nonzero child resultant, every fixed-rank arbitrary-gap remote cluster is
  eventually nonzero with a threshold independent of its internal diameter.
  The response is strictly Monge, mixed coefficients obey a factorwise
  Jensen contraction, and every compactification face is a positive
  generalized alternant on the nodes r^A.  A new A-layer marked-permutation
  model proves every mixed adjacent-difference tensor of order at least two
  positive, hence widths one and two are unconditionally good.  Minimal bad
  prefixes are finite and fixed-width bad supports have cumulative count
  O_t(X^(t-3)); this improves to O_t(X^(t-4)) if width three is separately
  good.  Width three for arbitrary product-Gamma shapes remains open.
source: root/multiscale-newton-flag/product-gamma-compactification-2026-08-02
audit: >
  An independent immutable-file hostile audit rederived the product response
  and endpoint derivative, factorwise Jensen and covariance, r^A block limit
  and Schur floor, physical elimination orientation, q>=2 and q=1 remote
  asymptotics, A-layer deletion/reinsertion inclusion-exclusion and strict
  marked-cycle family, low widths, minimal-prefix counts, and the exact
  shape-two certificate boundary.  Normal, optimized, and stored companions
  matched byte-for-byte; both LF hashes and documentation checks passed.
depends_on:
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3093-arbitrary-gap-remote-cluster-monge-flag-compactification
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2853-gamma-adjacent-tensor-cycle-weighted-positivity
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3050-rational-product-gamma-radial-nullcone-and-critical-borel-order
  - THM-3097-translated-support-monge-compactification-and-cofinite-bad-set-induction
script: 04-computation/product_gamma_monge_compactification_thm3100.py
output: 05-knowledge/results/product_gamma_monge_compactification_thm3100.out
script_sha256: 80f3f39afe3415edde0b08d1b36fda7bb40633acbc89a9f56d944b36ae69ca7a
output_sha256: 5c454ed5315a42314abf1de38d00145a6b2fba0bff84f6cd4f7000aa9d03f612
hash_basis: LF-normalized bytes
---

# THM-3100 -- product-Gamma response-Monge compactification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The factorial response used by THM-3093 is not isolated.  Every finite
product-Gamma moment sequence has the same strict response ordering, with
one important change: a compactification block sees the nodes `r^A`, where
`A` is the number of Gamma factors.  This gives both a whole-support
translation face and the arbitrary-gap remote-cluster face above any
nonzero child.

The extension has a sharp logical boundary.  Widths one and two remain good,
and in fact all positive adjacent-difference tensors remain positive by a
new multilayer permutation model.  The factorial width-three atomic
orientation does not survive even the Gamma shape `2`.  Thus the resulting
bad-prefix spectrum starts at width three, not width four, unless a separate
product-Gamma width-three theorem is supplied.

## 1. Fixed product-Gamma setting

Fix once and for all

```text
A>=1,                    theta_1,...,theta_A>0,          (1)

w_n=product_(nu=1)^A (theta_nu)_n,
L_theta(s^n)=w_n,        f_n=s^n/w_n.                   (2)
```

All thresholds below may depend on `A`, the complete shape vector `theta`,
the width, and any fixed child support.  None is uniform while a shape or the
number of factors moves.

If `G_nu` are independent `Gamma(theta_nu,1)` variables and
`X=product_nu G_nu`, then

```text
L_theta(P)=E[P(X)].                                      (3)
```

The law of `X` is continuous, has support `(0,infinity)`, and charges every
nonempty open interval.  In particular

```text
L_theta(P^2)>0
```

for every nonzero real polynomial `P`.  Also `L_theta(f_n)=1`.

For a support

```text
0<=a_1<...<a_t,                                         (4)
```

write

```text
H=sum_(i<t)x_i(f_(a_i)-f_(a_t))
```

and let `R_t(a_1,...,a_t)` be the standard-orientation homogeneous
resultant of

```text
L_theta(H^2),...,L_theta(H^t).                          (5)
```

Positive row normalizations may be inserted without changing its zero set.
Set `R_1=1`.  Call a support **good** when its resultant is nonzero and
**bad** otherwise.

## 2. Exact response and factorwise Jensen contraction

Fix an integer base `N>=0`.  For a degree `r>=1` and an integer gap `h>=0`,
put

```text
V_r(h)
 =L_theta(f_(N+h)^r)/L_theta(f_N^r)
 =product_nu
   (theta_nu+rN)_(rh)/(theta_nu+N)_h^r.                (6)
```

For `u<v`, define the root response

```text
phi_(u,v)(r)=(1/r)log[V_r(v)/V_r(u)].                  (7)
```

Gamma interpolation gives the exact integral

```text
phi_(u,v)(r)
 =sum_nu integral_u^v
   [psi(theta_nu+r(N+x))-psi(theta_nu+N+x)]dx.          (8)
```

Its `r` derivative is

```text
sum_nu integral_u^v (N+x)psi'(theta_nu+r(N+x))dx>0.    (9)
```

The integrand can vanish only at the single endpoint `N=x=0`; every
positive-length integer edge still gives a strict integral.  Thus `(8)` is
strictly increasing for every `N>=0`, including degree one and base zero.
Equivalently,

```text
[V_a(v)/V_a(u)]^b < [V_b(v)/V_b(u)]^a,
                                      1<=a<b.          (10)
```

This is the product-Gamma response-Monge law.

For distinct gaps `0=h_1<...<h_q`, the normalized coefficient of
`v^alpha` in a degree-`r` all-remote form is

```text
Q_(r,alpha)
 =product_nu
   (theta_nu+rN)_(sum_t alpha_t h_t)
   /product_t(theta_nu+N)_(h_t)^alpha_t,               (11)

sum_t alpha_t=r.
```

For each factor, `s |-> log(theta_nu+rN)_s` is convex and vanishes at zero.
Jensen, after adding zero gaps if `sum alpha_t<r`, gives

```text
Q_(r,alpha)^r<=product_t V_r(h_t)^alpha_t.             (12)
```

Multiplying the factorwise inequalities proves `(12)` for the whole
product.  This unsigned inequality is the physical error control; it is not
a termwise resultant-sign comparison.

The carriers also satisfy exact common-offset covariance:

```text
V_(r,N)(C+h)=V_(r,N)(C)V_(r,N+C)(h).                  (13)
```

No rationality assumption on the shapes enters `(6)--(13)`.

## 3. The compactification blocks see `r^A`

Apply the staircase row and column potentials of THM-3093 to `(10)`.  They
expose exactly the same path of tied pure terms:

```text
{(i,i),(i,i+1):i<q} union {(q,q)}.                    (14)
```

All normalized pure weights are at most one, and `(12)` bounds every mixed
coefficient by its pure envelope.

Suppose an adjacent gap stays fixed while its block origin `M` tends to
infinity.  For a fixed local gap `d`, `(6)` gives

```text
lim_(M->infinity)
 [V_(r,M)(d)]^(1/r)=r^(A d).                           (15)
```

Consequently a fixed compactification block `I=[a,b]`, with local gaps
`d_j=h_j-h_a`, has positive within-block potential multipliers `Lambda_j`
and full line matrix

```text
L_(i,j)=(r_i^A)^d_j Lambda_j/
        ((r_i^A)^d_i Lambda_i).                       (16)
```

Put `x_i=r_i^A`.  Since `A>=1`, the positive nodes `x_i` are strictly
increasing.  The row and column `Lambda` factors cancel in the determinant.
The remaining generalized alternant is a Vandermonde times a Schur
polynomial, and the Schur polynomial contains the staircase monomial with
coefficient one.  Hence

```text
det L_I
 =det(x_i^d_j)_(i,j in I)/product_(i in I)x_i^d_i

 >=Vandermonde(x_a,...,x_b)
   /product_(i=a)^b x_i^(i-a)>0.                       (17)
```

If several adjacent gaps diverge, pass to a subsequence on which every edge
is fixed or tends to infinity.  Strict response contraction kills all
backward terms across a divergent edge, while `(12)` kills every mixed term
containing a killed column.  The limiting system is block upper triangular.
THM-3073 removes every arbitrary forward term and gives the normal resultant

```text
product_(I in P)(det L_I)^D>0,                         (18)
```

where `P` is the consecutive-block partition and `D` is the product of the
row degrees.  There are only finitely many such partitions, and `(17)` gives
a positive floor on each.  This is the composition-cube compactification.

## 4. Unconditional whole-support translation

Fix a width `t>=1` and arbitrary distinct integer gaps

```text
0=h_1<...<h_t.                                         (19)
```

Then there is

```text
N_0=N_0(A,theta,t)                                     (20)
```

independent of the entire gap vector, such that

```text
R_t(N+h_1,...,N+h_t)>0                 for N>=N_0.     (21)
```

For `t=1` this is the convention `R_1=1`.  For `t>=2`, use the uneliminated
degree system `1,...,t` in the `t` normalized monomials.  Its first form is

```text
F_1=v_1+...+v_t.                                       (22)
```

The response flag of Sections 2--3 applies with row degrees
`r_i=i`.  If uniform positivity failed, a countersequence with `N->infinity`
would have a fixed/divergent edge partition.  Every block origin tends to
infinity, including the first, so `(15)--(18)` make its normalized resultant
converge to one of finitely many strictly positive face values.  This is a
contradiction.  Restoring the carriers in `(6)` only multiplies by positive
numbers.

For completeness, eliminate `(22)` by

```text
v_i=y_i                    (i<t),
v_t=w-sum_(i<t)y_i.                                      (23)
```

This coordinate map has determinant one and `F_1=w`.  Moving `w` past the
`t-1` variables contributes `(-1)^(t-1)`, raised by resultant covariance to
the degree product `t!`.  It is therefore `+1` for every `t>=2`.  Restriction
to `w=0` is exactly the standard mean-zero system `(5)`.  Thus the positive
uneliminated face is the claimed physical resultant, not a formal corner.

The threshold `(20)` is existential.  Compactness over finitely many edge
partitions does not provide a useful computable box.

## 5. Arbitrary-gap remote clusters above a fixed child

Fix a physical child support of width `m>=1` with

```text
S_m!=0.                                                 (24)
```

Use `S_1=1`.  Fix `q>=1`, put

```text
k=m+q,       r_i=m+i                    (1<=i<=q),
D_z=product_i r_i=k!/m!,
rho_m=(m/(m+1))^m<1.                                  (25)
```

Append an arbitrary remote cluster

```text
C+h_1,...,C+h_q,        0=h_1<...<h_q,                (26)
```

to the fixed child.  There is a threshold depending on the fixed parameters,
child, `m`, and `q`, but independent of every internal gap, such that the
new resultant is nonzero for all larger `C`.

More precisely, the weight-independent fixed-pivot coordinate and atom
inventory of THM-3093 apply unchanged.  Replacing its one factorial ratio by
the product in `(6)` gives, after the staircase normalization,

```text
R_k(C;h)
 =S_m^D_z Etilde_(C,h)^m!
  product_(i=1)^q U_(r_i,C+h_i)^[k!/r_i]
  [1+O(C^K rho_m^(A C))].                             (27)
```

Here `K` is fixed and every carrier `U` is positive.  For `q>=2`, `Etilde`
has a positive uniform lower bound supplied by `(17)--(18)`.  For `q=1`, its
one-normal limit is `(-1)^k`, so its absolute value has a uniform positive
lower bound instead.  The estimate is uniform over gap vectors of arbitrary
diameter.  Each Gamma factor contributes the same Stirling exponential base
as the factorial factor, so their product raises the contraction from
`rho_m^C` to `rho_m^(A C)`; the shapes affect only the fixed polynomial
invoice.

For `q>=2`, `D_z` contains consecutive integers and is even.  Thus `(27)` is
eventually positive, even if the child orientation is negative.  For `q=1`,
the statement is eventual nonvanishing with inherited sign
`sign(S_m)^(m+1)`; width two over the empty child is positive directly.

This theorem does not repair a zero child.  When `S_m=0`, THM-3073 makes the
displayed leading face vanish and does not identify the next physical factor.

## 6. An `A`-layer marked-permutation positivity theorem

The low-width input admits an independent combinatorial proof.  Put

```text
d_n=f_(n+1)-f_n,
R_n=w_(n+1)/w_n=product_nu(theta_nu+n),
E_n=w_(n+1)d_n=s^(n+1)-R_n s^n.                       (28)
```

For `k>=2` and `n_1,...,n_k>=0`, set

```text
M=sum_i(n_i+1).
```

Expansion gives

```text
W_theta(n_1,...,n_k)
 :=L_theta(product_i E_(n_i))

 =sum_(S subset [k])(-1)^|S|
   product_(i in S)product_nu(theta_nu+n_i)
   product_nu(theta_nu)_(M-|S|).                      (29)
```

There is a literal positive model for `(29)`.  On the same `M` labelled
objects, distinguish marks `x_i` and pairwise-disjoint private sets `O_i` of
size `n_i`.  Take `A` permutation layers.  Give layer `nu` the cycle weight
`theta_nu^(number of cycles)`.  In layer `nu`, call `x_i` bad when it is a
singleton cycle or its predecessor lies in `O_i`.  Call `x_i`
**catastrophic** when it is bad in every layer.

For a fixed set `S` of catastrophic marks, delete those marks in every
layer.  In layer `nu`, each deleted `x_i` is restored either as a singleton
cycle, of weight `theta_nu`, or after one of the `n_i` private labels.  The
private insertion sites are disjoint, so the exact intersection weight is

```text
product_(i in S)product_nu(theta_nu+n_i)
product_nu(theta_nu)_(M-|S|).                         (30)
```

Inclusion-exclusion in the catastrophic events proves that `(29)` counts
the `A`-tuples in which every mark is good in at least one layer.  It is
therefore coefficientwise nonnegative in the shape variables.

Strictness for `k>=2` is visible.  In one fixed layer put all `k` marks in a
single directed cycle, permute the ordinary labels arbitrarily, and leave
every other layer arbitrary.  This admitted family has positive weight

```text
(k-1)! theta_1(theta_1)_(sum_i n_i)
product_(nu>1)(theta_nu)_M.                           (31)
```

Consequently

```text
L_theta(product_i d_(n_i))>0              for k>=2.  (32)
```

By multilinearity, `(32)` holds for every mixed moment of two or more
nonzero positive adjacent-difference cone elements.  This extends THM-2853
from one Gamma layer to every finite product.

Widths one and two follow.  A mean-zero polynomial on two normalized slots
is a scalar multiple of

```text
f_b-f_a=sum_(i=a)^(b-1)d_i,
```

whose second moment is strictly positive by `(32)` (equivalently by the
full-support integral `(3)`).  Thus

```text
C_1=C_2=empty                                           (33)
```

for the minimal-bad-prefix sets below.

## 7. Finite first-bad prefixes and the counting spectrum

For fixed `A,theta`, let `B_t` be the bad width-`t` supports and let `C_t`
consist of the bad supports whose every proper initial prefix is good.

Each `C_t` is finite.  Otherwise, `(21)` bounds the minimum exponent on an
infinite sequence.  Pass to a subsequence with a longest fixed proper prefix
`c` of width `j`; then the next exponent tends to infinity.  By definition of
`C_t`, the fixed child `c` is good.  Apply Section 5 to the remaining
`t-j` arbitrary gaps (including the one-normal case).  The resulting support
is eventually good, a contradiction.

Every bad support has a unique first bad prefix.  Therefore, for

```text
B_t(X)={a in B_t:0<=a_1<...<a_t<=X},                  (34)
```

one has the cylinder upper bound

```text
#B_t(X)
 <=sum_(m=3)^t sum_(c in C_m, c_m<=X)
      binom(X-c_m,t-m)
 =O_(A,theta,t)(X^(t-3)).                              (35)
```

The cylinders are an upper envelope: extending a bad prefix need not produce
a bad support.  If width three is separately proved good for the chosen
product-Gamma family, then `C_3` is empty and `(35)` improves to

```text
#B_t(X)=O_(A,theta,t)(X^(t-4)).                        (36)
```

More generally, goodness through width `b` gives exponent `t-b-1`.
Unconditionally, `(35)` says that bad supports have relative density
`O(X^-3)` among all width-`t` supports, and `B_3=C_3` is finite.

There is also a coefficientwise Hilbert-series form.  A prefix
`c in C_m` with last exponent `c_m` has full extension-cylinder series

```text
sum_(X>=0)binom(X-c_m,t-m)z^X
 =z^(c_m+t-m)/(1-z)^(t-m+1).                          (37)
```

Summing `(37)` over the finite first-bad-prefix bank recovers `(35)` and
identifies width three as the first possible boundary stratum.

## 8. The exact width-three boundary

The permutation theorem `(32)` does not orient a complex two-plane.  The
factorial proof in THM-2824 uses the stronger atomic inequality

```text
D=2t222 g12-3t122 g22>=0.                              (38)
```

It already fails for one Gamma factor of shape `theta=2` on the support
`{0,1,2}`.  With

```text
U=d_0,              V=d_1,
gij=L_theta(U_i U_j),
tijk=L_theta(U_i U_j U_k),                             (39)
```

exact arithmetic gives

```text
g11=1/2,       g12=1/2,       g22=5/6,
t111=1/2,      t112=1,        t122=13/6,      t222=16/3,

D=-1/12.                                                  (40)
```

This is a certificate hostile, not a moment hostile.  The two exact
quadratic/cubic divisibility remainders are

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22=-1/2,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2=-11/36.     (41)
```

Thus this triple is good, but the factorial atomic route cannot prove the
general product-Gamma statement.  No theorem here asserts arbitrary
product-Gamma width-three goodness.  Accordingly `(36)` is conditional and
must not be cited as the unconditional counting exponent.

## 9. Sharp scope and degeneration walls

- `A>=1` and every `theta_nu>0` are load-bearing.  At `A=0`, `w_n=1` is
  evaluation at one, every response node is `1`, and every alternant of
  order at least two vanishes.
- The shapes, factor count, child, and width are fixed before the translation
  or remote offset tends to infinity.  Already for one shape,
  `L_theta((f_1-f_0)^2)=1/theta`; hence no positive floor is uniform as
  `theta->infinity`, where the normalized law approaches evaluation.
- Repeated gaps repeat a column and kill the corresponding alternant.
- The remote theorem requires `S_m!=0`.  It is not a zero-child continuation
  theorem.
- The translation and remote thresholds are non-effective and need not agree.
- Rational shapes are unnecessary here.  THM-3050 needs rationality for its
  finite-place nullcone proof; the analytic inequalities here do not.
- Nothing here proves product-Gamma SFC in width three or higher, the Strong
  Factorial Conjecture, NC2, arbitrary-radial GMC(2), LRC(14), JC(2), or
  DC(2).

## 10. Exact companion

The dependency-free rational companion checks:

1. `1,620` strict response-Monge cells, including degree one and base zero;
2. `12,384` factorwise multivariate Jensen cells;
3. `972` staircase dual weights and every exposed tie;
4. `540` composition-face conditions with nodes `r^A`;
5. `1,458` exact common-offset covariance cells;
6. `63` strict adjacent tensors, six literal multilayer permutation
   enumerations, and the marked-cycle lower family;
7. `36` exact width-two variances;
8. the shape-two hostile `(40)--(41)`; and
9. all `A=0` and moving-shape boundary controls.

Run

```text
python 04-computation/product_gamma_monge_compactification_thm3100.py
python -O 04-computation/product_gamma_monge_compactification_thm3100.py
```

Both executions must byte-match the stored transcript after LF
normalization.  The finite checks audit the exact identities and sharp
boundaries; the compactness and permutation arguments supply the universal
quantifiers.

**QED.**
