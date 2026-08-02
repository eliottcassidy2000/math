---
id: THM-3063
title: "Terminal suspension, transverse resultant, and the five-slot physical tail"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE-FILE AUDIT
  REQUESTED.  For every k>=3, a fixed-low physical factorial cell with one
  moving terminal pair has an eventual positive resultant whenever its
  transverse (k-2)-slot resultant is nonzero.  The exact leading suspension
  is S^[k(k-1)] times the two top normal forms, and every fixed affine finite
  generalized-Hankel bank is uniformly positive along a carrier-to-physical
  holotopy.  For k=5 the transverse cell is an arbitrary three-slot cell, so
  THM-2824 proves the hypothesis unconditionally.  The first new obstruction
  is k=6, exactly the arbitrary four-slot first-window resultant.
source: kind-pasteur-2026-08-01-terminal-suspension
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-3040-formal-corner-resultant-width-quotient-and-all-order-bernoulli-law
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
  - THM-3062-four-slot-physical-transverse-augmentation-and-affine-tail-holotopy
related:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-3060-three-slot-physical-terminal-face-and-affine-tail-holotopy
script: 04-computation/gmc_terminal_suspension_five_slot_tail_thm3063.py
output: 05-knowledge/results/gmc_terminal_suspension_five_slot_tail_thm3063.out
script_sha256: 8bb9dda3069ced0122ea40dfeeb6138ea5259a2e1cd6115e19fd0cd5814645db
output_sha256: 3532437efae3e3d182b4b9a8078491064db84aa527a53f8f05c525c6ccb2290f
hash_basis: LF-normalized bytes
---

# THM-3063 -- terminal suspension and the five-slot physical tail

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3060 and THM-3062 expose the same phenomenon at consecutive transverse
depths.  The all-large physical face always has only two normal directions.
The information apparently lost in the other coordinates is not arbitrary:
it is exactly the lower-slot factorial resultant on the transverse plane.
This gives a suspension theorem.  It proves the full fixed-low terminal tail
through five slots and identifies, without termwise factorial comparison,
the first genuinely new obstruction at six slots.

## 1. Physical family and transverse sidecar

Let `L(s^j)=j!`, `f_j=s^j/j!`, and fix

```text
k>=3, n>=1, 0<=a_0<...<a_(k-3)<C, h>=1, M=C+h.       (1)
```

Use homogeneous coordinates

```text
(x_0,...,x_(k-3),z) in P^(k-2),
w=sum_(i=0)^(k-3)x_i,                                 (2)
```

and put

```text
Y_C=sum_i x_i(f_(n+a_i)-f_(n+M))
       +z(f_(n+C)-f_(n+M)),                           (3)

F_(r,C)=L(Y_C^r)/L(f_n^r),                 2<=r<=k,
R_C=Res(F_(2,C),...,F_(k,C)).                         (4)
```

This is the normalized intrinsic physical resultant, not a formal prime
corner or one selected Macaulay chart.

For `k>=4`, the transverse space is

```text
K={w=z=0} ~= P^(k-4).                                 (5)
```

It has vector dimension `k-3`.  Restrict the fixed layers:

```text
H_r=F_(r,0,C)|_K,                         2<=r<=k-2,
S=Res_K(H_2,...,H_(k-2)).                           (6)
```

The forms in `(6)` are independent of `C`.  There are `k-3` of them in
`k-3` homogeneous variables.  At the endpoints use

```text
k=3: K=P^(-1), S=1;
k=4: K=P^0,    S=H_2(K).                              (7)
```

The theorem assumes only

```text
S!=0.                                                  (8)
```

Thus the fixed-low scalar cell is the complete obstruction passed into the
moving tail.  No sign hypothesis on `S` is needed.

## 2. Normal filtration and the exact top pair

Split each physical form by the number of factors drawn from the moving pair
`{f_(n+C),f_(n+M)}`:

```text
F_(r,C)=sum_(q=0)^r F_(r,q,C).                        (9)
```

Every moving factor is linear in `(w,z)`, so exactly

```text
F_(r,q,C) belongs to (w,z)^q.                        (10)
```

For `r=k-1,k`, define

```text
U_(r,C)=L(f_(n+C)^r)/L(f_n^r)
       =(rn+1)_(rC)/(n+1)_C^r,
g_(r,C)=F_(r,C)/U_(r,C).                              (11)
```

The all-high binary normal blocks have full coefficientwise inverse-`C`
expansions and leading forms

```text
F_(r,r,C)/U_(r,C) -> L_r(w,z)^r,
L_r=-r^h w+(1-r^h)z.                                  (12)
```

Put

```text
p=k-1, q=k, Delta=k^h-(k-1)^h>0.                     (13)
```

The two limiting lines are independent.  Choose once a generic rational
complete flag in `K` for which the limiting lower Poisson boundary factors
are nonzero, and let `N~=P^1` be its final complementary normal line.  Define
the exact, flag-dependent boundary factor

```text
E_C=Res_N(g_(k-1,C)|_N,g_(k,C)|_N).                  (14)
```

This uses the full normalized boundary restrictions, not merely their
all-high pieces.  It is therefore the literal final Poisson factor.  Yet
`(12)` gives

```text
E_C=Delta^[k(k-1)](1+O(C^-1)),                       (15)
```

with a complete inverse-power expansion.  In particular `E_C!=0` for all
sufficiently large `C`.  The value in `(14)` depends on the auxiliary flag;
the final symbol below does not.

## 3. Uniform tube around the transverse space

For fixed `k,n,(a_i),h`, coefficient norms satisfy

```text
||F_(r,j,C)/U_(r,C)||
 <=poly(C)(j^j/r^r)^C,                 0<=j<r,        (16)
```

with `0^0=1` for the fixed layer.  The exact normal pair in `(14)` has a
uniformly nonzero binary resultant.  Hence its two powered directions are
jointly coercive.  At a common zero of `g_(k-1,C),g_(k,C)`, `(10),(16)` give
the two root inequalities

```text
c R^r<=poly(C)sum_(j<r)(j^j/r^r)^C R^j,
r=k-1,k,                                                (17)
```

where `R` is projective distance to `K`.  Since

```text
(j^j/r^r)^[1/(r-j)]<=1/r,                             (18)
```

the elementary polynomial root bound yields

```text
R<=poly(C)(k-1)^(-C).                                 (19)
```

On this tube every lower form satisfies, uniformly on both standard charts
of `K`,

```text
F_(r,C)=H_r+O(poly(C)rho_k^C),             2<=r<=k-2,
rho_k=((k-2)/(k-1))^(k-2)<1.                         (20)
```

Indeed a `j`-high layer costs at most
`(j/(k-1))^(jC)` after `(19)`; convexity puts the largest endpoint at
`j=k-2`, which also dominates `j=1`.  The number `rho_k` is a safe gap, not
a sharp one.  For `k=3` there is no lower layer and `(20)` is vacuous.

This is the load-bearing separation.  It controls whole contracted forms on
their common zero scheme and never compares equal-base Macaulay monomials
term by term.

## 4. Repeated Poisson contraction

Scale only the top two coefficient blocks in the resultant.  Their
multidegrees are

```text
mu_(k-1)=k!/(k-1),             mu_k=(k-1)!.           (21)
```

Repeated homogeneous Poisson formulas along the generic flag from `(14)`
are exact for all sufficiently large `C`.  At each stage, `(19)--(20)`
replace the next lower form on the finite intersection scheme by its
restriction to `K`; local multiplicities are retained.  The separated
source model is

```text
(H_2,...,H_(k-2)) in the K variables,
(A_(k-1),A_k)     in the two normal variables.        (22)
```

Its standard product formula is

```text
Res(22)=S^[k(k-1)] E^[(k-2)!].                        (23)
```

The exponents have a direct intersection meaning.  The normal complete
intersection has length `k(k-1)`, so it repeats the transverse resultant
that many times.  Conversely the lower degrees have product

```text
2*3*...*(k-2)=(k-2)!,                                (24)
```

so they repeat the binary normal resultant `(k-2)!` times.  Equivalently,
every exponent in `(23)` is forced by the coefficient multidegrees.  The
Poisson boundary factors telescope along the flag, leaving exactly `(23)`.

Consequently the actual physical resultant satisfies the refined formula

```text
R_C=S^[k(k-1)] E_C^[(k-2)!]
    U_(k-1,C)^[k!/(k-1)] U_(k,C)^[(k-1)!]
    (1+O(poly(C)rho_k^C)).                            (25)
```

This proof also covers multiple roots: first perturb the fixed generic flag
and lower coefficients so every Poisson fibre is reduced, apply the exact
identity, and specialize back by polynomial continuity.  No root is counted
without its local intersection multiplicity.

Combining `(15)` and `(25)` gives the flag-free leading face

```text
R_C~S^[k(k-1)] Delta^(k!)
     U_(k-1,C)^[k!/(k-1)] U_(k,C)^[(k-1)!].           (26)
```

Because `k(k-1)` and `k!` are even and every `U` is positive, `(8)` implies

```text
R_C>0 for every sufficiently large C.                 (27)
```

For `k=3`, `(14)` is the whole binary resultant and the lower Poisson product
is empty.  For `k=4`, `(25)` is THM-3062's exact `G_2^12 E_C^2` mechanism.

## 5. Exact exponential base and power

Gauss multiplication gives

```text
U_(r,C)=r^(rC) product_(j=1)^(r-1)
         [(n+j/r)_C/(n+1)_C],                         (28)

U_(r,C)~kappa_(r,n) r^(rC) C^[-(r-1)/2],
kappa_(r,n)=Gamma(n+1)^(r-1)/
             product_(j=1)^(r-1)Gamma(n+j/r)>0.       (29)
```

Put

```text
alpha_k=(k-2)! [2(k-1)^2-1]/2.                       (30)
```

Then `(26)` becomes

```text
R_C=A [k(k-1)]^(k! C) C^(-alpha_k)(1+O(C^-1)),       (31)

A=S^[k(k-1)] Delta^(k!)
  kappa_(k-1,n)^[k!/(k-1)]
  kappa_(k,n)^[(k-1)!]>0.                             (32)
```

All normalized coefficients in `(14)` are finite Gamma-ratio expressions,
and the remaining error in `(25)` is exponentially small.  Thus `(31)` has
a complete inverse-power expansion.

## 6. Unconditional five-slot theorem

Take `k=5`, so the fixed lows are arbitrary

```text
0<=a<b<c<C.                                           (33)
```

Now `K=P^1`, and

```text
S_3=Res_K(H_2,H_3).                                   (34)
```

THM-2824 proves that the positive-definite binary quadratic `H_2` and the
real binary cubic `H_3` have no common projective root for every three-slot
factorial support.  In the standard orientation

```text
S_3=leading(H_2)^3 |H_3(alpha)|^2>0.                 (35)
```

Therefore `(8)` is automatic and

```text
R_C=S_3^20 E_C^6 U_(4,C)^30 U_(5,C)^24
    (1+O(poly(C)(27/64)^C)),                          (36)

R_C~S_3^20(5^h-4^h)^120
     kappa_(4,n)^30 kappa_(5,n)^24
     20^(120C) C^(-93)>0.                             (37)
```

This proves the actual normalized five-slot physical terminal resultant is
eventually positive for every fixed triple of lower offsets and fixed gap.

The next case is exact rather than rhetorical.  At `k=6`, `(6)` is the
arbitrary four-slot resultant of degrees `2,3,4` on `P^2`.  THM-2843 gives
its nonnegative norm reduction but not universal nonvanishing.  Hence six
slots are blocked precisely by arbitrary four-slot first-window SFC; no new
terminal-tail cancellation is hiding before that sidecar.

## 7. Exact carrier and affine tail holotopy

Let

```text
A_k=k!(Harmonic(k)-1),
B_k=k!(k+1-2Harmonic(k)),
I_k=A_k+B_k=k!(k-Harmonic(k)).                        (38)
```

THM-3040's exact physical width flag at `t=1/n` is

```text
W_M=(n)_M^A_k (n+1)_M^B_k/n^(I_k M).                 (39)
```

Define the positive exact carrier

```text
B_C=S^[k(k-1)] Delta^(k!)
    U_(k-1,C)^[k!/(k-1)] U_(k,C)^[(k-1)!],
G_C=W_(C+h)B_C,                    Q_C=W_(C+h)R_C.    (40)
```

Formula `(28)` makes `B_C` a scaled product of Beta moments, while `(39)`
is a full-support product-Gamma moment.  Hence `(G_C)` is a strict Stieltjes
product.

Fix integers `s>=1,b_0>=0`, put `C_j=sj+b_0`, and retain the tail on which
all inequalities in `(1)` hold.  For

```text
T_j^(theta)=(1-theta)G_(C_j)+theta Q_(C_j),
0<=theta<=1,                                          (41)
```

every fixed finite generalized-Hankel bank is positive for all sufficiently
large `j`, uniformly in `theta`.  The actual/carrier ratio has a common full
symbol by `(25),(31)`.  The width flag supplies the first confluent
alternant

```text
(I_k s)^N V(u)V(v)/product_(ell=0)^(r-1)ell!>0,
N=r(r-1)/2,                                           (42)
```

and the fixed bank makes the remainder uniform.

For five slots,

```text
(A_5,B_5,I_5)=(154,172,326),                         (43)
```

and the carrier has thirty copies of each
`Beta(n+j/4,1-j/4)`, `j=1,2,3`, plus twenty-four copies of each
`Beta(n+j/5,1-j/5)`, `j=1,2,3,4`.  Its affine alternant order is `326s`.

## 8. Exact evidence and scope

The exact companion checks the normal-ideal filtration, safe layer gaps,
separated resultant exponent ledger, five-slot transverse binary resultants,
full normal-boundary resultants, Gauss/Beta products, `(154,172,326)` width
flag, and finite carrier Hankel minors.  It also builds the literal degree-11
physical Macaulay map over `F_1000003` on `36` mixed cells.  All maps have
shape `589x364` and full column rank.  This is a finite nonvanishing control,
not an intrinsic resultant value or a tail proof.  No unproved selected
`364x364` determinant chart is used as evidence.

Run

```text
python 04-computation/gmc_terminal_suspension_five_slot_tail_thm3063.py
python -O 04-computation/gmc_terminal_suspension_five_slot_tail_thm3063.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem is fixed-low, fixed-gap, and affine-clock.  It proves eventual
physical resultant positivity and every fixed finite Hankel bank, not one
all-order Stieltjes tail.  It does not handle moving lower offsets, varying
gap, nonlinear clocks, raw-chart persistence, or a wall-stripped core.  It
does not prove arbitrary four-slot SFC, the six-slot tail without `(8)`, or
any new NC2, GMC, LRC, or Jacobian statement.

**QED, pending independent immutable-file audit and status promotion.**
