---
id: THM-3065
title: "Reciprocal Beta Gregory-Newton wall and generalized Hankel sign regularity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a,b>0,
  every strict generalized-Hankel minor of
  m_k=(b)_k/(a)_k has the sign of product_(j<r)(a-b)_j.  The entire
  arbitrary-offset determinant is that universal product times an explicit
  strictly positive Selberg-Schur/Littlewood-Richardson factor.  Positive
  integer reciprocal gaps are exact Gregory-Newton degree walls; positive
  nonintegral gaps give nonzero sign-regular kernels with a sharp wall-crossing
  parity.  The same ordered-prefix calculus gives one rational-mesh carrier
  certificate for THM-3062.  Nonpositive dual prefixes force every order-two
  sign but not the higher checkerboard law: exact zero-cut hostiles fail at
  orders three and four.  Higher reciprocal-carrier signs and the strict-dual-
  prefix repair are finite-exact signals only.
source: root-and-gmc-tight-sequence-frontier-2026-08-01
audit: >
  An independent audit rederived the Gregory--Newton termination law, the
  one-sided alternant, the Selberg--Schur/Littlewood--Richardson positive
  factor after cancellation and analytic continuation, every integer-wall
  multiplicity, and the arbitrary-mesh prefix identity.  It separately
  checked the zero-cut H3/H4 hostiles and the summation-by-parts H2 survivor,
  then replayed normal and optimized modes byte-identically to the stored
  transcript and matched both declared LF hashes.  The strict-negative-prefix
  all-order repair remains explicitly OPEN.
depends_on:
  - THM-3053-beta-gamma-prefix-transport-and-multiplicative-holotopy-cone
  - THM-3062-four-slot-physical-transverse-augmentation-and-affine-tail-holotopy
script: 04-computation/gmc_reciprocal_beta_generalized_hankel_wall_thm3065.py
output: 05-knowledge/results/gmc_reciprocal_beta_generalized_hankel_wall_thm3065.out
script_sha256: 83faf0f4a46e7e9af6506cf8e4c25bf8d99f63aae0280e8c273c3a02ec51f1f6
output_sha256: f5e9ac33a111e9b58d626a04a4ecca5d706366d6b823df100861075ba2d9cb3d
hash_basis: LF-normalized bytes
---

# THM-3065 -- reciprocal Beta gaps are Gregory-Newton sign walls

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3053 separates the positive-prefix multiplicative cone from the full
Stieltjes cone.  Its smallest forward edge has moments `(a)_k/(b)_k` with
`a<b`.  Reversing that edge does not merely destroy positivity.  It produces
an exact sign-regular kernel whose signature changes only when the shape gap
crosses an integer.  The integers are not accidental Selberg poles: they are
the termination walls of the Gregory--Newton expansion.

The one-variable difference law, the one-sided alternant, and the full
two-sided generalized minor are three lifts of the same falling-factorial
object.

## 1. Universal generalized-minor factorization

Fix real `a,b>0` and put

```text
m_k=(b)_k/(a)_k,                     k>=0,             (1)
delta=a-b,                           d=b-a=-delta.     (2)
```

For strictly increasing nonnegative integer tuples

```text
u_1<...<u_r,                         v_1<...<v_r,      (3)
```

define the generalized-Hankel minor

```text
D_(u,v)(a,b)=det[m_(u_i+v_j)]_(i,j=1)^r.             (4)
```

Then there is an explicit factorization

```text
D_(u,v)(a,b)=Phi_(u,v)(a,b)
              product_(j=0)^(r-1)(delta)_j,          (5)

Phi_(u,v)(a,b)>0                    for every a,b>0.  (6)
```

Consequently every strict generalized minor of a fixed order has the same
sign, including arbitrary gaps on both the row and column sides:

```text
sgn D_(u,v)(a,b)
 =sgn product_(j=0)^(r-1)(a-b)_j.                    (7)
```

Formula `(5)`, rather than only `(7)`, is proved below.  In particular there
is no hidden genericity condition on the offsets.

## 2. Gregory--Newton is the elementary wall mechanism

Let `Delta f_k=f_(k+1)-f_k`.  Directly from `(1)`,

```text
Delta m_n=d m_n/(a+n).                                (8)
```

Induction, with the denominator shape shifting once at each difference,
gives the exact all-order law

```text
Delta^q m_n
 =d^(underline q) m_n/(a+n)_q,

d^(underline q)=d(d-1)...(d-q+1).                    (9)
```

Thus `d=m` a nonnegative integer is precisely Newton termination:

```text
Delta^(m+1)m_n=0,
m_n=(a+m)_n/(a)_n=(a+n)_m/(a)_m,                     (10)
```

so `(m_n)` is a polynomial sequence of exact degree `m`.  For nonintegral
`d>0`, no finite difference vanishes.  Moreover

```text
(a-b)_j=(-d)_j=(-1)^j d^(underline j),               (11)
```

so the universal determinant factor in `(5)` is exactly the product of the
Gregory--Newton difference coefficients, with the usual alternant sign
`(-1)^(r choose 2)`.  This is the structural reason integer and noninteger
reciprocal gaps behave differently.

## 3. Elementary one-sided alternant

Keep arbitrary rows `u_1<...<u_r`, but take consecutive columns beginning at
`t>=0`.  Then

```text
det[m_(u_i+t+j)]_(i=1..r,j=0..r-1)
 =V(u) product_(j=0)^(r-1)(a-b)_j
    product_(i=1)^r [m_(u_i+t)/(a+u_i+t)_(r-1)],     (12)

V(u)=product_(i<j)(u_j-u_i)>0.                       (13)
```

This is elementary.  Factor `m_(u_i+t)` from row `i` and write `x=u_i`.
After multiplying that row by `(a+x+t)_(r-1)`, column `j` becomes the degree
`r-1` polynomial

```text
P_j(x)=(b+x+t)_j(a+x+t+j)_(r-1-j).                   (14)
```

The determinant of the `P_j(x_i)` is `V(x)` times the determinant of their
coefficient matrix.  Evaluate instead at `x_k=-b-t-k`, `0<=k<r`.  The
evaluation matrix is triangular and its diagonal is

```text
(-k)_k(a-b)_(r-1-k)
 =(-1)^k k!(a-b)_(r-1-k).                            (15)
```

The Vandermonde of those evaluation nodes is
`(-1)^(r choose 2) product_(k=0)^(r-1)k!`; division leaves exactly
`product_(j<r)(a-b)_j`, proving `(12)`.

Taking `u_i=k+i-1` gives the contiguous closed form

```text
det[m_(k+i+j)]_(i,j=0)^(r-1)
 =product_(j=0)^(r-1)
   [j!(b)_(k+j)(a-b)_j/(a)_(k+r+j-1)].               (16)
```

Thus the familiar contiguous wall is only the consecutive-offset shadow of
the full generalized law.

## 4. The positive factor before analytic continuation

It remains to prove that arbitrary columns introduce no new sign.  This step
is made explicit before continuation.

First suppose `delta>0`.  The normalized Beta integral is

```text
m_k=1/B(b,delta)
    integral_0^1 t^(b+k-1)(1-t)^(delta-1)dt.          (17)
```

Put `p=u_1`, `q=v_1`, `alpha=b+p+q`, and define partitions

```text
lambda_i=u_(r+1-i)-p-(r-i),
mu_i    =v_(r+1-i)-q-(r-i),            1<=i<=r.      (18)
```

The generalized alternants are

```text
det[t_j^(u_i)]=product_j t_j^p Delta(t)s_lambda(t),
det[t_j^(v_i)]=product_j t_j^q Delta(t)s_mu(t).       (19)
```

Andreief applied to `(17)` therefore gives

```text
D_(u,v)=1/(r!B(b,delta)^r)
 integral_[0,1]^r Delta(t)^2 s_lambda(t)s_mu(t)
 product_j t_j^(alpha-1)(1-t_j)^(delta-1)dt_j.       (20)
```

Use the finite Littlewood--Richardson expansion

```text
s_lambda s_mu=sum_nu c_(lambda,mu)^nu s_nu,
c_(lambda,mu)^nu>=0.                                 (21)
```

For completeness, the needed Schur--Selberg evaluation is itself a
consequence of Andreief and the elementary formula `(12)`.  If
`k_i=nu_(r+1-i)+i-1`, apply Andreief to
`Delta(t)s_nu(t)=det[t_j^(k_i)]` and `Delta(t)`.  The resulting determinant
of Beta integrals is one-sided, so `(12)` gives

```text
I_nu(alpha,delta)
 :=integral_[0,1]^r s_nu(t)Delta(t)^2
       product_j t_j^(alpha-1)(1-t_j)^(delta-1)dt_j

 =S_r(alpha,delta)s_nu(1^r)
   product_(i=1)^r
   [(alpha+r-i)_(nu_i)
    /(alpha+delta+2r-i-1)_(nu_i)],                   (22)
```

where

```text
S_r(alpha,delta)
 =product_(j=0)^(r-1)
   [Gamma(alpha+j)Gamma(delta+j)Gamma(j+2)
    /Gamma(alpha+delta+r+j-1)].                       (23)
```

The Weyl alternant quotient gives the displayed `s_nu(1^r)`; no analytic
sign assertion is used in deriving `(22)`.

Substitute `(21)--(23)` into `(20)` and cancel the `Gamma(delta)` factors
*before* continuing.  The result is `(5)` with

```text
Phi_(u,v)(a,b)=C_r(a,b,alpha) P_(lambda,mu)(a,b),     (24)

C_r=1/r! [Gamma(a)/Gamma(b)]^r
    product_(j=0)^(r-1)
    [Gamma(alpha+j)Gamma(j+2)
     /Gamma(a+p+q+r+j-1)],                           (25)

P_(lambda,mu)
 =sum_nu c_(lambda,mu)^nu s_nu(1^r)
   product_(i=1)^r
   [(alpha+r-i)_(nu_i)
    /(a+p+q+2r-i-1)_(nu_i)].                         (26)
```

Every argument in `(25)` is positive for `a,b>0`.  Every summand in `(26)`
is nonnegative, at least one Littlewood--Richardson coefficient is positive,
and every displayed Pochhammer ratio is positive.  Hence `(6)` is literal.

Now hold `b,u,v` fixed and write `a=b+delta`.  After the cancellation just
performed, both sides of `(5),(24)--(26)` are analytic throughout the
connected interval `delta>-b`: all remaining Gamma arguments are positive
and every finite denominator is nonzero.  Equality on `delta>0` therefore
continues to every `a,b>0`.  Crucially, the continued factor is visibly the
positive expression `(24)--(26)`; positivity is not inferred from a
meromorphic continuation.

This proves `(5)--(7)`.

## 5. Exact integer ranks and nonintegral sign chambers

The factorization gives all boundaries without root finding.

If `a>b`, every `(a-b)_j` is positive, so `(m_k)` is strictly generalized-
Hankel totally positive.  This is also immediate from its nondegenerate Beta
law.

If `a=b`, then `m_k=1` and every minor of order at least two vanishes.

If `d=b-a=m` is a positive integer, `(10)` and `(5)` give

```text
rank[m_(u+v)]_(u,v>=0)=m+1,                           (27)

sgn D_(u,v)=(-1)^(r choose 2),       r<=m+1,
D_(u,v)=0,                           r>=m+2.          (28)
```

In fact every strict generalized minor through order `m+1` is nonzero, not
only one leading minor.

Let instead `m<d<m+1` with integer `m>=0`.  There are no zeros and

```text
sgn D_(u,v)=(-1)^E,
E=sum_(j=0)^(r-1) min(j,m+1).                         (29)
```

Equivalently,

```text
E=binom(r,2),                                  r<=m+2,
E=binom(m+2,2)+(r-m-2)(m+1),                  r>=m+2. (30)
```

At the integer wall `d=m`, every order `r>=m+2` minor has the universal zero
multiplicity

```text
r-m-1,                                               (31)
```

because precisely the factors `(a-b)_j`, `m+1<=j<=r-1`, contain `m-d`.
The sign across that wall flips exactly when `r-m-1` is odd.  These are the
sharp nonintegral-gap sign chambers.

## 6. Arbitrary shape meshes and the THM-3062 carrier

THM-3053's prefix proof does not require unit-spaced shapes.  For any

```text
0<alpha_0<...<alpha_N,
G_C=c^C product_(j=0)^N(alpha_j)_C^(n_j),            (32)
S_j=sum_(i=0)^j n_i,                                 (33)
```

the same cut argument gives a cancellation-free Gamma/forward-Beta product
if and only if every `S_j>=0`.  The canonical identity is

```text
product_j(alpha_j)_C^(n_j)
 =(alpha_N)_C^(S_N)
   product_(j=0)^(N-1)
   [(alpha_j)_C/(alpha_(j+1))_C]^(S_j).              (34)
```

The ratio in `(34)` is `Beta(alpha_j,alpha_(j+1)-alpha_j)`.  Conversely, a
cut after `alpha_j` counts Gamma residues to its left and forward Beta edges
crossing it, so its value is nonnegative.  This proves the arbitrary-mesh
extension, including necessity.

Apply it to THM-3062.  Remove only positive constants and a positive geometric
factor from the auxiliary base-width carrier `W_C B_C`.  Its ordered shape
inventory is

| shape | `n` | `n+1/4` | `n+1/3` | `n+1/2` | `n+2/3` | `n+3/4` | `n+1` |
|---|---:|---:|---:|---:|---:|---:|---:|
| multiplicity | `26` | `6` | `8` | `6` | `8` | `6` | `-14` |
| prefix | `26` | `32` | `40` | `46` | `54` | `60` | `46` |

This is one rational-mesh positive-prefix certificate, rather than a list of
unrelated Beta factors.

The literal carrier is `G_C=W_(C+h)B_C`.  Up to a positive constant
independent of `C`,

```text
W_(C+h)/W_C
 =product_(t=0)^(h-1)
   [(n+t+1)_C/(n+t)_C]^26
   [(n+t+2)_C/(n+t+1)_C]^20.                         (35)
```

Thus the fixed terminal gap adds `26` forward transfers from shape `n` to
`n+h` and `20` transfers from `n+1` to `n+h+1`.  For `h=1`, the nonzero
prefixes are

```text
6,14,20,28,34,26,46.                                 (36)
```

For `h>=2`, they are

```text
6,14,20,28,34,0,26,46.                               (37)
```

Hence the actual carrier has a single positive-prefix certificate with
terminal mass `46` for every fixed `h>=1`.  This recovers its strict
Stieltjes property and explains the net Gamma order `46` in THM-3062's
confluent alternant.

## 7. What reciprocity preserves and destroys

For a single forward Beta edge, reciprocity changes `(low)_k/(high)_k` into
`(high)_k/(low)_k`.  Sections 1--5 completely describe what replaces total
positivity: a sign-regular kernel, Newton rank termination at integral gap,
and nonterminating sign chambers at nonintegral gap.

There is also one universal statement for reciprocals of arbitrary strict
moment sequences.  If `(G_k)` is strictly generalized-Hankel TP of order two
and `H_k=1/G_k`, then for row gap `p>0`, column gap `q>0`, and base `x`,

```text
det[[H_x,H_(x+q)],[H_(x+p),H_(x+p+q)]]
 =[G_(x+p)G_(x+q)-G_x G_(x+p+q)]
   /[G_x G_(x+p)G_(x+q)G_(x+p+q)]<0.                (38)
```

Thus the reciprocal THM-3062 carrier is already excluded from the Stieltjes
cone at every strict generalized minor of order two.

The natural dual-prefix lift stops exactly there.  More generally, let

```text
H_k=c^k product_(j=0)^N(alpha_j)_k^(e_j),
0<alpha_0<...<alpha_N,                E_j=sum_(i<=j)e_i.    (39)
```

If every `E_j<=0` and the inventory is nonzero, then every generalized minor
of order two is negative.  Indeed, for the continuous ratio

```text
q(x)=H_(x+1)/H_x=c product_j(x+alpha_j)^(e_j),

d/dx log q(x)
 =E_N/(x+alpha_N)
  +sum_(j=0)^(N-1) E_j
   [1/(x+alpha_j)-1/(x+alpha_(j+1))]<0.               (40)
```

Strict decrease of `q` is strict log-concavity, and multiplying the intervening
ratios proves the arbitrary-gap order-two sign.  This is the strongest
survivor of the dual-prefix conjecture.

It does **not** lift to the checkerboard signature at higher order.  The
smallest three-shape zero-cut hostile found is

```text
(alpha_0,alpha_1,alpha_2)=(1/7,2/3,20),
(e_0,e_1,e_2)=(-1,1,-1),              E=(-1,0,-1),

H_k=(2/3)_k/[(1/7)_k(20)_k].                          (41)
```

Its leading order-three determinant has the wrong sign:

```text
det[H_(i+j)]_(i,j=0)^2
 =4914161/84138683904000
 =7^3*14327/[2^12*3^10*5^3*11^2*23] >0,             (42)
```

where the checkerboard signature at order three is negative.  A powered
version records the failure independently at even checkerboard parity:

```text
(e_0,e_1,e_2)=(-3,3,-1),              E=(-3,0,-1),
u=(0,1,2,7),                          v=(0,1,2,6).    (43)
```

Its order-four minor, which a checkerboard law would require to be positive,
is instead

```text
- 7^58 * 1953596693 * 4434246478119732782469492691
  /[2^39 * 3^89 * 5^14 * 11^4 * 13^2 * 19^3
    *23^3 * 29^6 * 43^6 * 71^3] < 0.                 (44)
```

Thus the first failed implication is

```text
nonpositive dual prefixes
  does not imply all-order checkerboard sign regularity.                (45)
```

The missing coordinate is quantitative: edge multiplicities, shape gaps,
and terminal reciprocal-Gamma mass remain visible beyond order two.  The two
nearest strict-prefix repairs `(-2,1,-1)` and `(-3,2,-1)` have prefixes
`(-2,-1,-2)` and `(-3,-1,-2)` and pass all `6,664` generalized-minor controls
through order five and offsets in `[0,6]`.  Independently, all `194` exponent
inventories in `[-4,4]^3` with every prefix strictly negative passed all
`194*1,225=237,650` arbitrary order-four controls on the same three shapes.
This supports, but does not prove, the sharper hypothesis that **strictly**
negative prefixes force the checkerboard law.  For THM-3062 this is still not
a complete explanation: the actual carrier has terminal mass `46`, while for
`h>=2` it also has the zero cut displayed in `(37)`.

The exact companion finds the checkerboard signature
`(-1)^(r choose 2)` for `26,656` reciprocal-carrier minors through order five,
two bases, four terminal gaps, and every row/column offset set in `[0,6]`.
This is a `FINITE-EXACT` signal only.  Hadamard products of individually
sign-regular kernels are not proved here to preserve that signature; no
all-order reciprocal-carrier claim is made.

## 8. Exact evidence and scope

The dependency-free exact companion checks:

- `1,152` Gregory--Newton identities and `162` integer terminations;
- `53,312` arbitrary-two-sided generalized minors through order five across
  forward, equal, positive-integer, and positive-noninteger gap chambers;
- `5,712` one-sided alternant equalities and `504` contiguous formulas;
- `72` hostile nonconsecutive minors immediately below, on, and above four
  integer walls, including the crossing parity `(31)`;
- the rational-mesh carrier prefixes, `120` actual-gap flow identities, and
  the `26,656` explicitly scoped reciprocal-carrier controls;
- the exact dual-prefix hostiles `(41)--(44)`, `441` surviving order-two
  controls, and `6,664` nearest strict-prefix repair controls.

Run

```text
python 04-computation/gmc_reciprocal_beta_generalized_hankel_wall_thm3065.py
python -O 04-computation/gmc_reciprocal_beta_generalized_hankel_wall_thm3065.py
```

Ordinary and optimized execution are byte-identical to the stored eleven-line
transcript.  The LF-normalized script has `16,831` bytes in `459` lines and
the output has `791` bytes in `11` lines; their hashes are pinned in the
frontmatter.

The theorem concerns one Gamma-ratio kernel at all orders, and the positive-
prefix factorization of the exact carrier.  It does not classify arbitrary
signed Gamma inventories outside the multiplicative cone, promote the finite
higher reciprocal-carrier signature, prove that reciprocal Hankel signs are
closed under Hadamard product, or add a physical GMC/SFC/NC2/LRC/Jacobian
consequence.

**QED.**
