# A self-contained higher-order certificate for the endpoint phase

2026-10-19. This certificate is an independent control for the next-order audit. It strengthens the numerical enclosure of C_2; it does not strengthen the next-order error theorem.

## 1. Definitions and an exact residual identity

Let f(y)=y+y^2, a=-1/2, b>0, u_k=f^k(a), and v_j=g^j(b), where g is the nonnegative inverse of f. Let

h(y)=-1/y+log|y|,

D_0(y)=1/(1+y)+log(1+y)-1.

The endpoint phase is

C_R=h(b)-h(a)-sum_(k>=0) D_0(u_k)-sum_(j>=1) D_0(v_j),  R=b+1/2.

Set H(y)=h(y)+P(y), where

P(y)=-y/2+y^2/3-13y^3/36+113y^4/240-1187y^5/1800
     +877y^6/945-14569y^7/11760+176017y^8/120960.

Define the exact Abel residual

r(y)=H(f(y))-H(y)-1
    =D_0(y)+P(y+y^2)-P(y).

The logarithmic identity holds on both relevant real sides because 1+y>0. In particular, r extends smoothly through 0 and r(0)=0. Exact rational algebra gives

r'(y)=y^9 Q(y)/(1+y)^2,

where

Q(y)=(1745717+10608000y+30074982y^2+49440366y^3
      +50918352y^4+33678411y^5+13995135y^6
      +3344323y^7+352034y^8)/15120.

The accompanying script reconstructs P coefficient by coefficient by canceling the Abel residual through degree nine, then verifies this derivative identity exactly with SymPy. Thus the residual estimate is based on a rational identity, not on fitting a numerical convergence rate.

## 2. Analytic tail enclosure

Fix eta=1/100 and write Q(y)=sum_(i=0)^8 q_i y^i. For |y|<=eta,

|r'(y)| <= K |y|^9,

K=(sum_(i=0)^8 |q_i| eta^i)/(1-eta)^2.

Integrating from 0 gives

|r(y)| <= (K/10)|y|^10.

For any cutoff m with -eta<u_m<0<v_m<eta, both future tails stay inside this interval. The squared-increment identities are

sum_(k>=m) u_k^2=-u_m,
sum_(j>m) v_j^2=v_m.

Consequently,

sum_(k>=m) |r(u_k)| <= (K/10)|u_m|^9,
sum_(j>m) |r(v_j)| <= (K/10)v_m^9.

To connect the residual to the original normalization, telescope P(f(y))-P(y) on each half-orbit. Since P(y)->0 at y->0,

sum_(k>=m) r(u_k)=sum_(k>=m) D_0(u_k)-P(u_m),
sum_(j>m) r(v_j)=sum_(j>m) D_0(v_j)+P(v_m).

The finite h-identities then give the exact relation

C_R=H(v_m)-H(u_m)+2m
    -sum_(k>=m) r(u_k)-sum_(j>m) r(v_j).

Hence, with A_m=H(v_m)-H(u_m)+2m,

|C_R-A_m| <= (K/10)(|u_m|^9+v_m^9).

All bounds concern the residual itself, rather than rounded numerical differences between two large approximate phases.

## 3. Directed-interval evaluation

The script uses mpmath.iv with 80 decimal digits to propagate u and v from the exact rational decimal inputs. It verifies the cutoff condition by interval comparisons, evaluates A_m and K in directed interval arithmetic, and adds the full symmetric analytic tail interval. The positive inverse is evaluated as 2v/(1+sqrt(1+4v)), avoiding subtractive cancellation.

For b=3/2 and m=512, the interval arithmetic gives

C_2 in [-2.463836653392077114885601747333544260665627649968553426299030462009799029538349464708,
        -2.463836653392077114885581328452590992276863620693049760239525536347905589664417200762].

The analytic tail half-width is at most 1.020944047663419438202e-23 (rounded upward). A simple outward-rounded certificate, asserted by the script, is therefore

-2.463836653392077114885602 < C_2 < -2.463836653392077114885581.

The full endpoints above are computational interval output, not claims that all their displayed digits are uniquely determined digits of C_2. The short outward enclosure is the intended certificate. It lies strictly inside the initial note's valid coarse enclosure.

Certificates for b=0.05,0.25,10, a second cutoff for b=1.5, and the Abel-shift consistency check at b=3.75 are included in the captured output. All finite threshold-root checks are separately labeled NUMERICAL and are not used to certify this interval or the theorem's asymptotic remainder.

Reproduce with `python independent_phase_controls.py` using Python 3, SymPy 1.14.0, and mpmath 1.3.0. No network or packet imports are needed.
