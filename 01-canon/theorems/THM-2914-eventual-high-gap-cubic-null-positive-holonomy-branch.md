---
id: THM-2914
title: "Eventual high-gap cubic-null positive-holonomy branch"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING FINAL INDEPENDENT
  HOSTILE AUDIT.  For every sufficiently large integer r, the ordered
  support (0,1,2,r) has a locally unique positive high-shared cubic-null
  branch at scale ((3r)!/(r!)^3)^(-1/3).  Its normalized quartic endpoint
  determinant converges to an explicit positive factored limit, so the
  branch is not quartic-null.  The threshold is existential and no global
  branch uniqueness or arbitrary-support closure is claimed.
source: root/eventual-high-gap-cubic-null-branch-2026-07-29
depends_on:
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
related:
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
  - THM-2908-consecutive-four-slot-projective-resultant-closure
  - THM-2910-nonconsecutive-cubic-null-endpoint-holonomy-sign-reversal
script: 04-computation/gmc_eventual_high_gap_cubic_null_branch_thm2914.py
output: 05-knowledge/results/gmc_eventual_high_gap_cubic_null_branch_thm2914.out
script_sha256: c43a4f77b5afe3c29e6c49bfa63386453ecbb1f519c410df2ccd90c45f53ab72
output_sha256: 80d44c8930b325b02b22620b5f239af18c40a0560786bf2d25724eacada9e253
hash_basis: LF-normalized bytes
---

# THM-2914 -- eventual high-gap positive cubic-null branch

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING FINAL INDEPENDENT
HOSTILE AUDIT.**

Put

```text
L(s^k)=k!,                    f_j=s^j/j!,
d_j=f_(j+1)-f_j,              C_r=f_r-f_2.             (1)
```

For each integer `r>=3`, consider

```text
U_r=d_0+xC_r,                 V_r=d_1+yC_r.            (2)
```

Define

```text
T_k(r)=(kr)!/(r!)^k,
epsilon_r=T_3(r)^(-1/3),
H_r=epsilon_r^4 T_4(r)
   =T_4(r)/T_3(r)^(4/3)>0.                             (3)
```

Then there are a positive point `(xi_*,eta_*)` and an integer `R_0`
such that for every integer `r>=R_0` there is a unique point
`(xi_r,eta_r)` in one fixed neighborhood of `(xi_*,eta_*)` for which

```text
x_r=epsilon_r xi_r>0,          y_r=epsilon_r eta_r>0   (4)
```

and the binary quadratic moment of `U_r+zV_r` divides its binary cubic
moment.  Moreover,

```text
(xi_r,eta_r)=(xi_*,eta_*)+O(2^(-r)),                  (5)

J_r(x_r,y_r)/H_r
 -> (eta_*^2-2xi_*^2)
    (eta_*^2-4eta_*xi_*+2xi_*^2)>0,                  (6)
```

where `J_r` is THM-2872's quartic endpoint coefficient determinant.
Consequently `J_r(x_r,y_r)>0` for every sufficiently large `r`, and
this local positive cubic-null branch never extends to a quartic
multipole line.

The uniqueness in this statement is only inside the fixed scaled
neighborhood.  It says nothing about other positive branches or other
parameter scales.

## 1. Exact moment expansion and the finite error list

Write

```text
x=epsilon_r xi,               y=epsilon_r eta,
E(z)=d_0+z d_1.                                      (7)
```

For `m=2,3,4`, the binary moment is exactly

```text
M_(m,r)(z)
 =L((U_r+zV_r)^m)

 =sum_(k=0)^m binom(m,k) epsilon_r^k
    (xi+eta z)^k L(E(z)^(m-k) C_r^k).                (8)
```

There are no asymptotic interchanges in `(8)`.  If

```text
P(s)=sum_d p_d s^d
```

is any fixed polynomial, direct expansion of
`C_r=f_r-f_2` gives

```text
L(P C_r^k)
 =sum_(j=0)^k (-1)^(k-j) binom(k,j) 2^(-(k-j))
    sum_d p_d (jr+2(k-j)+d)!/(r!)^j.                 (9)
```

Every polynomial `P` needed in `(8)` has degree at most eight.
Therefore all errors belong to a finite list of factorial-ratio
sequences.

For

```text
M_(j,a)(r)=(jr+a)!/(r!)^j,                            (10)
```

the quadratic and cubic errors are finite sums of sequences of the
forms

```text
T_3^(-1/3) M_(1,a),
T_3^(-2/3) M_(2,a),
M_(j,a)/T_3,                     j<3.                 (11)
```

Their successive-term ratios tend to at most

```text
1/3,             4/9,             4/27,              (12)
```

respectively.  For the quartic moment divided by `H_r`, every lower
term has one of the forms

```text
T_3^((4-k)/3) M_(j,a)/T_4,       j<=k<=3,
M_(j,a)/T_4,                     j<4, k=4.            (13)
```

The largest successive-ratio limit in `(13)` is `81/256`.
Since the lists of `(j,k,a)` are finite, the ratio test yields, on
every compact set `K` in the `(xi,eta)` plane,

```text
||F_r-F_infinity||_(C^1(K))
 +||J_r(epsilon_r xi,epsilon_r eta)/H_r
      -J_infinity||_(C^0(K))
 <=C_K 2^(-r)                                         (14)
```

for all sufficiently large `r`.  Here `F_r=(I_(1,r),I_(2,r))` is the
pair of cubic-divisibility invariants.  Thus `(14)` is a uniform
coefficient estimate, not merely a pointwise Stirling heuristic.

For orientation only, the familiar leading scales are

```text
T_2/T_3^(2/3)
 =Theta(r^(1/6)(4/9)^r) ->0,

H_r=T_4/T_3^(4/3)
 =Theta(r^(-1/6)(256/81)^r) ->infinity.              (15)
```

Equation `(9)` supplies the rigorous proof used in `(14)`.

## 2. The limiting cubic-null point

Literal factorial evaluation gives the low quadratic and cubic
vectors

```text
(g_0,g_1,g_2)=(1,1,2),
(t_0,t_1,t_2,t_3)=(2,4,10,30).                       (16)
```

Equations `(8)--(14)` imply

```text
q_r(z) -> 1+2z+2z^2,

c_r(z) -> c_0(z)+(xi+eta z)^3,                       (17)
```

where `c_0` has coefficient vector `(2,4,10,30)` in the binary
binomial convention.  Hence

```text
F_(infinity,1)
 =-eta^3+6eta xi^2-4xi^3-14,

F_(infinity,2)
 =-2(eta^3-3eta^2xi+2xi^3+4).                       (18)
```

Eliminating `xi` gives

```text
Res_xi(F_(infinity,1),F_(infinity,2))

 =-256(2eta^9+18eta^6-729eta^3+54).                 (19)
```

Put `u=eta^3` and

```text
p(u)=2u^3+18u^2-729u+54.                             (20)
```

It has exactly one real root in each of

```text
(-25,-24),              (0,1),              (15,16). (21)
```

The two positive roots are the only possibilities with `eta>0`.
On `p(u)=0`, the other coordinate satisfies either exact form

```text
xi/eta=(2u^2+21u-603)/189
      =(u+9)/(2u-3).                                  (22)
```

For `0<u<1`, the first expression in `(22)` is negative, so that root
does not lie in the positive quadrant.  The root `u_*` in `(15,16)`
gives the unique positive-quadrant limiting point

```text
eta_*=u_*^(1/3),
xi_*=eta_*(2u_*^2+21u_*-603)/189,

6/7<xi_*/eta_*<35/27.                                (23)
```

Direct substitution verifies both equations `(18)`.  Their Jacobian is

```text
det D F_infinity
 =18(eta^2-2eta xi+2xi^2)^2>0                        (24)
```

away from the origin.  In particular it is nonzero at
`(xi_*,eta_*)`.

## 3. Uniform continuation by contraction

Let

```text
p_*=(xi_*,eta_*),             A=(D F_infinity(p_*))^(-1).
```

Choose `rho>0` so that the closed ball

```text
B=closed_ball(p_*,rho)
```

lies in the positive quadrant and

```text
sup_(w in B)||I-A D F_infinity(w)||<=1/4.             (25)
```

By `(14)`, for every sufficiently large integer `r`,

```text
sup_B||A(D F_r-D F_infinity)||<=1/4,
||A F_r(p_*)||<=rho/2.                               (26)
```

The map

```text
Phi_r(w)=w-AF_r(w)                                    (27)
```

is then `1/2`-Lipschitz on `B`, and `(26)` makes it map `B` into
itself.  Banach's fixed-point theorem gives one and only one
`p_r=(xi_r,eta_r)` in `B` satisfying `F_r(p_r)=0`.
The same estimate gives `(5)`.

The two polynomials in `(2)` are independent because the `f_0`
coefficient occurs only in `U_r`.  Since

```text
L(P^2)=integral_0^infinity P(s)^2 exp(-s) ds>0
```

for nonzero real `P`, their quadratic Gram determinant is positive.
Thus `F_r=0` is exactly quadratic divisibility of the cubic, not a
degenerate remainder artifact.

This proof establishes `exists R_0`.  Formula `(9)` makes `R_0`
effectively computable, but no numerical threshold is asserted here.

## 4. The positive quartic exit

After division by `H_r`, only the pure `C_r^4` term survives in the
quartic moment:

```text
A_(i,r)/H_r -> xi^(4-i) eta^i,              0<=i<=4. (28)
```

Combining `(28)` with `(16)` in THM-2872's endpoint determinant gives

```text
J_infinity
 =eta^4-4xi eta^3+8xi^3eta-4xi^4

 =(eta^2-2xi^2)(eta^2-4eta xi+2xi^2).                (29)
```

Let `tau=xi_*/eta_*`.  From `(23)`,

```text
1-2tau^2<=-23/49,

1-4tau+2tau^2<=-601/729.                             (30)
```

For the second inequality, the quadratic is convex and its values at
the endpoints in `(23)` are `-47/49` and `-601/729`.
Therefore

```text
J_infinity(p_*)
 >=(13823/35721) eta_*^4>0.                          (31)
```

Equations `(5)`, `(14)` and `(31)` prove `(6)` and eventual positivity
of `J_r`.  Quartic divisibility would force the two linear remainder
coefficients, hence `J_r`, to vanish by THM-2872.  This proves the
quartic exit.

## 5. Meaning and boundary

THM-2910 found a negative branch at top slot five and a positive branch
at top slot six.  The present theorem explains why positive endpoint
sign is stable in the high-gap limit: the high atom contributes a pure
rank-one cubic and an overwhelmingly larger rank-one quartic, while
every mixed low/high correction decays geometrically after scaling.

This is not a global branch classification.  In particular:

1. the fixed-point theorem gives uniqueness only in one scaled
   neighborhood of `p_*`;
2. other positive or signed branches and other parameter scales are not
   excluded;
3. no explicit value of `R_0` is proved;
4. middle, low and projective-boundary charts are untouched; and
5. there is no arbitrary-support, SFC(4), or GMC consequence.

The useful survivor is a new analytic tail module:

```text
finite exact support atlas + factorial-ratio tail stability
 -> global support-family coverage.                            (32)
```

For a complete SFC(4) result, the finite prefix and all other projective
charts would still need sign-blind full-remainder certificates.

## 6. Exact companion

The exact companion verifies:

1. the low moment vectors `(16)`;
2. the limiting invariants `(18)` and eliminant `(19)`;
3. all three real root boxes `(21)` and both selector forms `(22)`;
4. the positive ratio interval `(23)` and Jacobian `(24)`;
5. the endpoint factorization and quantitative floor `(29)--(31)`; and
6. every exponential base used to separate the finite error types.

It locks the four limiting polynomials by SHA-256 digests and uses
explicit `require` gates in ordinary and optimized mode.  Run

```text
python3 04-computation/gmc_eventual_high_gap_cubic_null_branch_thm2914.py
python3 -O 04-computation/gmc_eventual_high_gap_cubic_null_branch_thm2914.py
```

Both modes byte-match

```text
05-knowledge/results/gmc_eventual_high_gap_cubic_null_branch_thm2914.out.
```

The companion checks the finite algebraic core.  The uniform tail and
contraction argument are proved in Sections 1 and 3; the output does
not pretend to print an unproved numerical cutoff.
