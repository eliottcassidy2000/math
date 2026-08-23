---
id: THM-3896
title: "Cusp residual degree-four equality-seam emptiness"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  In the THM-3881 residual problem, the n=4
  THM-3884 equality seam is empty.  After the two proved gauge jets, the
  middle collision has the Kummer form q_0=x s^2 and r_4=rho x^3 s.  If
  s is not proportional to x, the degree-17 square-root equation fails at
  x-order three.  If s is proportional to x, the degree-18 and degree-17
  equations can lift, but the degree-16 equation has an unavoidable
  coefficient 2 sigma^6/81 at x^4 y^12.  This closes only the degree-four
  equality cell; JC(2) remains open.
source: tournament-jc-breakthrough / degree-four Kummer-shell scout, 2026-08-23
audit: >
  SELF-AUDITED PROVISIONAL CANDIDATE.  The focused exact companion rederives
  the n=4 Kummer collision, verifies the complete homogeneous
  parameterization, the A/B gauge-transport identities, the DQ high shell,
  both x-adic branches, and the two decisive coefficients.  A positive
  control reaches both leading square-root equations before failing at the
  claimed degree-16 coefficient.  Normal and optimized runs byte-match the
  frozen transcript.  Independent hostile audit must recheck the degree
  exhaustion, the UFD Kummer normalization, the h_2 divisibility gate, and
  the degree-16 valuation bound before status promotion.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
  - THM-3886-cusp-residual-equality-seam-second-layer-trichotomy
related:
  - THM-3894-cusp-residual-all-degree-gauge-kummer-parity-filtration
script: 04-computation/jc2_cusp_residual_degree_four_equality_seam_thm3896.py
output: 05-knowledge/results/jc2_cusp_residual_degree_four_equality_seam_thm3896.out
script_sha256: ec0e410392cab779713bf6e8dbb2cc7b1f8fdf7b4cc50849f61baeb14e8cb08f
output_sha256: ad402ebe3e3d2a6c588cdcd18d2942aa81a2bc3ca1956fbab74c68e012e8ab99
semantic_sha256: d18d8a8674ea353b7a387746b31af8414fe14b6fa6f927be9a7a23cfccba6b67
hash_basis: raw LF bytes
---

# THM-3896 -- the degree-four equality seam is empty

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Use the THM-3881 notation

```text
a=x+1,                     L=9x+4,
K=K_2+K_1+K_0,
K_2=y^2-15x^2,             K_1=-15x,             K_0=-4,
P=aL^2,                    Delta=a^3L^2-K^2,              (1)
r=aT+Kf,                   A=KT+aPf,
B=Pf^2-T^2.                                                (2)
```

Let `S(T,f)` be the residual in THM-3881 equation (21).  The provisional
claim is

```text
deg f=4,  deg T=5,  S(T,f) a square,
f_4=xq_0, T_5=-K_2q_0, q_0 nonzero
                         ==> contradiction.               (3)
```

Thus the candidate excludes the `n=4` equality seam.  It makes no assertion
about `deg T>=deg f+2` or the rest of the planar Jacobian problem.

## 1. Re-deriving the middle Kummer symbol

THM-3886 supplies the second gauge jet: there is a homogeneous quadratic
`q_1` such that

```text
f_4=xq_0,                  T_5=-K_2q_0,
f_3=q_0+xq_1,              T_4=-K_2q_1+15xq_0.           (4)
```

Let `r_4` be the next homogeneous part of `r`.  The leading forms retained
from THM-3884 are

```text
A_8=81x^5q_0,              B_11=81x^5q_0^2.              (5)
```

At total degree nineteen, all residual terms except `8AB+3r^2B` are lower,
and direct substitution gives

```text
S_19=243q_0^2x^5(r_4^2+216q_0x^5).                       (6)
```

Degree nineteen is odd.  Since `q_0` is nonzero, a square residual therefore
requires

```text
r_4^2=-216q_0x^5.                                         (7)
```

Unique factorization now gives, after absorbing a scalar square,

```text
q_0=xs^2,                  r_4=rho*x^3s,
deg s=1,                   rho^2=-216.                    (8)
```

Indeed, every prime other than `x` occurs to even order in `q_0`, while
`v_x(q_0)+5` is even.  Hence `q_0/x` is a square.  This is a self-contained
rederivation of the `n=4` specialization of the provisional THM-3894
passport; THM-3894 is related signal, not a proved dependency here.

## 2. Complete parameterization after the collision

The degree-four part of `r=aT+Kf`, using `(4)`, is

```text
r_4=xT_3+K_2(f_2-q_1)-15x^2q_1-4xq_0.                   (9)
```

Reduce `(9)` modulo `x`.  Since `K_2 mod x=y^2` and the right side of `(8)`
is divisible by `x`, one obtains `x | (f_2-q_1)`.  Thus there is a unique
homogeneous linear `q_2` for which

```text
f_2=q_1+xq_2,
T_3=-K_2q_2+15xq_1+4q_0+rho*x^2s.                       (10)
```

Write the unrestricted lower pieces as a linear `u`, a quadratic `tau_2`,
a linear `tau_1`, and a scalar `z`.  Every pair satisfying `(4),(8),(10)`
and the THM-3881 address is then, uniquely,

```text
f_4=xq_0,                    T_5=-K_2q_0,
f_3=q_0+xq_1,                T_4=-K_2q_1+15xq_0,
f_2=q_1+xq_2,                T_3=-K_2q_2+15xq_1+4q_0+rho*x^2s,
f_1=q_2+u,                   T_2=15xq_2+4q_1+tau_2,
f_0=z,                       T_1=4q_2+tau_1,
                              T_0=4z.                    (11)
```

Put

```text
W=q_0+q_1+q_2,
bar f=f-aW=u+z,
bar T=T+KW=rho*x^2s+tau_2+tau_1+4z.                     (12)
```

The full pair is `(T,f)=(bar T-KW,bar f+aW)`.  Therefore

```text
r=a bar T+K bar f,
A=bar A+Delta W,
B=bar B+2W bar A+Delta W^2,                              (13)
```

where

```text
bar A=K bar T+aP bar f,          bar B=P bar f^2-bar T^2. (14)
```

In particular,

```text
r_4=rho*x^3s,
r_3=rho*x^2s+x tau_2+K_2u,
r_2=tau_2+x tau_1+K_1u+K_2z.                             (15)
```

No lower coordinate has been discarded.

## 3. The `D Q` high shell

Define

```text
D=Delta W^2,                    Q=8Delta W+3r^2.          (16)
```

The relevant homogeneous pieces of `Delta` are

```text
Delta_5=81x^5,
Delta_4=315x^4-K_2^2,
Delta_3=475x^3+30xK_2.                                  (17)
```

Equation `(8)` gives the exact middle cancellation

```text
Q_8=8Delta_5q_0+3r_4^2=0.                               (18)
```

Since `deg bar A<=5`, `deg bar B<=6`, `deg W=3`, and `deg D<=11`, expansion
of `(13)` gives

```text
(8A+3r^2)B
=(Q+8 bar A)(D+2W bar A+bar B).                          (19)
```

After `(18)`, every term in `(19)` outside `DQ` has degree at most sixteen.
All other macros in `S` have degree at most fourteen.  Consequently

```text
S_d=(DQ)_d                         for d>=17,             (20)
S_18=D_11Q_7,
S_17=D_11Q_6+D_10Q_7.                                  (21)
```

The required components are

```text
D_11=Delta_5q_0^2,
D_10=2Delta_5q_0q_1+Delta_4q_0^2,
D_9 =Delta_5(q_1^2+2q_0q_2)
      +2Delta_4q_0q_1+Delta_3q_0^2,                     (22)
Q_7 =8(Delta_5q_1+Delta_4q_0)+6r_4r_3,
Q_6 =8(Delta_5q_2+Delta_4q_1+Delta_3q_0)
      +3(r_3^2+2r_4r_2).                                (23)
```

At degree sixteen, the same expansion gives the more precise formula

```text
S_16=D_11Q_5+D_10Q_6+D_9Q_7+8 bar A_5D_11.              (24)
```

The would-be `Q_8` terms vanish by `(18)`; every other correction has degree
at most fifteen.

## 4. A terminal line transverse to `x` dies at degree seventeen

Write

```text
s=alpha*x+beta*y.                                        (25)
```

First suppose `beta!=0`.  From `(17),(22),(23)`, the exact `x`-orders are

```text
v_x(D_11)=7,                 v_x(Q_7)=1,
v_x(D_10)=2.                                             (26)
```

The minimum terms are

```text
D_10=-beta^4*x^2*y^8+O(x^3),
Q_7 =-8beta^2*x*y^6+O(x^2).                              (27)
```

Hence `v_x(S_18)=8`.  If `S=G^2` and

```text
G=g_9+g_8+...                                             (28)
```

in homogeneous pieces, `S_18=g_9^2` forces `v_x(g_9)=4`.  But `D_11Q_6`
is divisible by `x^7`, while `(27)` gives the exact coefficient

```text
[x^3y^14] S_17=8beta^6!=0.                              (29)
```

Thus `v_x(S_17)=3`, contradicting the next square-root equation

```text
S_17=2g_9g_8,                                             (30)
```

whose right side is divisible by `x^4`.

## 5. The line `s=sigma*x` dies at degree sixteen

It remains to take

```text
s=sigma*x,                 sigma!=0,
q_0=sigma^2x^3.                                           (31)
```

Now

```text
v_x(D_11)=11,                v_x(Q_7)=3,
v_x(D_10)=6.                                             (32)
```

Moreover

```text
(Q_7/x^3) mod x=-8sigma^2y^4.                            (33)
```

If `S_18` is not a square, the branch is already empty.  Otherwise `(21),
(32),(33)` allow a homogeneous quadratic `h_2` such that

```text
Q_7=x^3h_2^2,
g_9=9sigma^2x^7h_2,
h_2^2 mod x=-8sigma^2y^4.                                (34)
```

Put

```text
M=sigma^2Delta_4+162x^2q_1.                              (35)
```

Then

```text
D_10=sigma^2x^6M,
S_17=sigma^2x^9Mh_2^2+81sigma^4x^11Q_6.                 (36)
```

Because `h_2 mod x` is a nonzero multiple of `y^2`, one has
`gcd(h_2,x)=1`.  The equation `S_17=2g_9g_8` is polynomial only if

```text
h_2 divides Q_6.                                         (37)
```

If `(37)` fails, the branch is empty at degree seventeen.  If it holds,
write `Q_6=h_2c_4`; then `(36)` determines

```text
g_8=(x^2Mh_2)/18+(9/2)sigma^2x^4c_4.                    (38)
```

Let `h_2 mod x=eta*y^2`, so `eta^2=-8sigma^2`.  Since

```text
M mod x=-sigma^2y^4,                                     (39)
```

equation `(38)` has the unavoidable leading term

```text
g_8=-(sigma^2eta/18)x^2y^6+O(x^3).                       (40)
```

On the other hand, `(22),(24),(31),(32),(37)` give

```text
v_x(D_9)>=3,
v_x(S_16)>=min(11,6,3+3,11)=6.                           (41)
```

Combining `(40),(41)` yields the exact coefficient

```text
[x^4y^12](S_16-g_8^2)=2sigma^6/81!=0.                   (42)
```

Thus `v_x(S_16-g_8^2)=4`, while the final square-root equation

```text
S_16-g_8^2=2g_9g_7                                      (43)
```

has right side divisible by `x^7`.  This contradiction closes the remaining
terminal line.

## 6. Positive controls

The last obstruction is not a consequence of an already-empty leading
shell.  For every linear `s`, take

```text
q_1=q_2=u=0,
tau_2=(17rho/18)xs=-204xs/rho.                            (44)
```

Then exact reduction by `rho^2=-216` gives

```text
Q_7/x=-8s^2K_2^2,                                        (45)
```

which is a square after adjoining `lambda` with `lambda^2=-8`.

On the exceptional line `s=x`, additionally take

```text
tau_1=(7rho/72)x=-21x/rho,                               (46)
```

while leaving `z` arbitrary.  Then

```text
Q_7=-8x^3K_2^2,
Q_6=6(40+rho*z)x^4K_2.                                  (47)
```

Thus `h_2=lambda K_2` satisfies `(34),(37)`: both the degree-eighteen and
degree-seventeen equations genuinely lift.  Reconstructing the full pair

```text
f=ax^3+z,
T=-Kx^3+rho*x^3+(17rho/18)x^2+(7rho/72)x+4z             (48)
```

and expanding the complete THM-3881 residual recovers

```text
[x^4y^12](S_16-g_8^2)=2/81.                              (49)
```

This is a positive control for every preceding layer and a hostile control
for the terminal obstruction.

## 7. Scope and reproduction

Sections 1--6 form a **provisional proof candidate**, not promoted canon.
An independent hostile audit is still required.  If promoted, the result
would close only `deg f=4`, `deg T=5` on the THM-3884 equality seam.  It does
not close other even Kummer terminals, the odd square passport, the
off-equality region, a Keller atlas, or JC(2).

Reproduce the exact packet with

```bash
python3 04-computation/jc2_cusp_residual_degree_four_equality_seam_thm3896.py
python3 -O 04-computation/jc2_cusp_residual_degree_four_equality_seam_thm3896.py
```

Both streams must byte-match
`05-knowledge/results/jc2_cusp_residual_degree_four_equality_seam_thm3896.out`.
