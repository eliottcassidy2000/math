---
id: THM-3839
title: "Constant-tower bichromatic R2Z profiles do not enter the cubic pseudo-plane Darboux packet"
status: >
  PROOF CANDIDATE + VERIFIED-EXACT, PENDING INDEPENDENT HOSTILE AUDIT.  In
  the first r^2z extension of the cubic pseudo-plane canonical packet,
  suppose the preceding first-row rz^2 profile and the new first-row r^2z
  profile are nonzero.  If the two top Wronskians have distinct target colors
  lambda!=mu, their next bucket forces an 8/7 Kummer tower; the following
  bucket is either zero or its 8/5 rung.  When the common Kummer polynomial
  is constant, an exact factor-cofactor first integral and two target-color
  jets exclude every polynomial degree of the arm profile f.  The second
  tower coefficient is retained throughout: setting it to one is not a
  valid normalization.  Nonconstant Kummer towers, the aligned color,
  one-sided top profiles, higher slots, and the general planar Jacobian
  problem remain OPEN.
source: jc_zero_debt_lift / first-higher-canonical-slot lane, 2026-08-23
audit: >
  SELF-HOSTILE EXACT CANDIDATE.  An initial beta=1 proof was rejected after
  restoring the constant tower scalar exposed three nonlinear scaling
  seams.  The repaired companion keeps beta=b!=0 and checks the seams and
  the resulting all-b equations directly against the unreduced Poisson
  bracket.  Its 67 active gates verify the Poisson Casimir, monic canonical
  reduction, both top Wronskians, the 8/7 and optional 8/5 valuation laws,
  the complete factor-cofactor integral, the arm, every raw constant-tower
  bucket with b retained, the two tail-independent color jets, all terminal
  degree cases, and fixed degree 3--6 hostile replays.  Normal and optimized
  runs byte-match the frozen transcript.  Independent hostile audit remains.
depends_on:
  - THM-3821-cubic-pseudoplane-rz2-odd-ladder-terminal-riccati-gate
related:
  - THM-3836-cubic-factor-cofactor-darboux-packet
  - THM-3834-one-sided-second-row-r2z2-profile-nonentry
  - THM-3814-nodal-rz-kummer-profile-degree-gate
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
script: 04-computation/jc2_cubic_pseudoplane_constant_tower_bichromatic_r2z_thm3839.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_constant_tower_bichromatic_r2z_thm3839.out
script_sha256: dbfcc8082ef36a2686eba2fc2fdd896628348c6a21d13a783f20bc1678edcd0e
output_sha256: ddd9e81a83faa6b9bb181a100d86aa645aec2ed62ba8b26eacac3efc804490bc
semantic_sha256: 9150cf81c49c17d3d44dc7973b9c41183f1cc7fd7ef5a77b11ea287ce270ba09
hash_basis: raw LF bytes
---

# THM-3839 -- the constant bichromatic `r^2z` tower is empty

**PROOF CANDIDATE + VERIFIED-EXACT, PENDING INDEPENDENT HOSTILE AUDIT.**
Let `k` be an algebraically closed field of characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary profiles

```text
f,g,h,kappa,p,q,S,T,U,V in k[e]
```

consider the first `r^2z` extension of the THM-3821 packet:

```text
A=e^2-z/3+r g+z^2 f+rz p+rz^2 S+r^2z U,
C=e^3-e-ez/2+r h+z^2 kappa+rz q+rz^2 T+r^2z V.       (2)
```

Assume

```text
U!=0,        S!=0,
V=lambda U, T=mu S,        lambda!=mu.                (3)
```

The exact high buckets below force a polynomial `v` and constants
`alpha,beta in k*` such that

```text
U=alpha e^5v^8,             S=beta e^4v^7.             (4)
```

If `v` is constant, then

```text
{A,C}!=1.                                                (5)
```

This theorem closes the constant-`v` part of the genuinely bichromatic
branch (3).  It does not claim that `v` is always constant.

## 1. Two colors and the `8/7` tower

Reduce `{A,C}-1` by the monic relation

```text
z^3=r^2e+r.                                               (6)
```

The remainder has `z`-degree at most two and there is no quotient loss.  Its
two top Wronskians are

```text
[r^5z]=24e(UV'-VU'),
[r^5] =21e^2(ST'-TS').                                   (7)
```

Since `U,S` are nonzero, vanishing of (7) says

```text
(V/U)'=0,                    (T/S)'=0                    (8)
```

in `k(e)`.  The constant field is `k`, so (8) gives exactly the two colors
in (3), including the possibilities `lambda=0` or `mu=0`.

Put

```text
Delta=mu-lambda!=0,
M=kappa-mu f,       N=q-lambda p,       H=h-lambda g.   (9)
```

The next bucket is

```text
[r^4z^2]
 =3(lambda-mu)(7eS U'-8eU S'-3SU).                     (10)
```

Thus each irreducible prime `pi` gives a logarithmic valuation equation.  At
`pi=e` and `pi!=e`, respectively,

```text
7 ord_e(U)-8 ord_e(S)=3,
7 ord_pi(U)-8 ord_pi(S)=0.                              (11)
```

The nonnegative solutions are

```text
ord_e(U)=5+8m,      ord_e(S)=4+7m,
ord_pi(U)=8m,       ord_pi(S)=7m,       m>=0,           (12)
```

which proves (4).  Conversely (4) solves (10) identically, so the tower is
both necessary and exact.

After (4), `[r^4z]` is

```text
3(8eU N'-5eN U'+UN).                                   (13)
```

Hence either `N=0`, or the same prime-by-prime argument gives the optional
`8/5` rung

```text
N=gamma e^3v^5,                gamma in k*.             (14)
```

We retain both cases by allowing `gamma in k` below.

## 2. The factor-cofactor first integral

The arm bucket is

```text
D_mu f=2eM-1/12,             D_mu=3e^2-2mu e-1.        (15)
```

In particular,

```text
f(0)=1/12.                                               (16)
```

Substitute (4) and (14) into the complete pure `[r^4]` bucket and define

```text
K=7beta M+8alpha vH.                                    (17)
```

The bucket is exactly

```text
3e^5v^6(evK'-4ev'K-2vK).                               (18)
```

The domain permits cancellation of the displayed nonzero polynomial factor.
Moreover, (18) is precisely

```text
(K/(e^2v^4))'=0                                         (19)
```

up to a nonzero rational factor.  Therefore its complete integral in `k(e)`
is

```text
7beta M+8alpha vH=eta e^2v^4,          eta in k.        (20)
```

This factor-cofactor law is internal to the pseudo-plane profile calculation.
Its resemblance to THM-3836's intrinsic cubic factor-cofactor packet is a
structural connection, not a logical identification of the two objects.

## 3. The honest constant-tower equations

Now suppose `v in k*`.  Absorb its powers into the three constants and write

```text
U=a e^5,       S=b e^4,       N=n e^3,
a,b in k*,     n in k.                                     (21)
```

Set

```text
rho=a/b in k*,          delta=n/b in k.                    (22)
```

After renaming `eta`, equations (15) and (20) become

```text
D_mu f=2eM-1/12,
7M+8rho H=eta e^2.                                        (23)
```

It is essential that `b` remain in the calculation.  Setting `b=1` changes
three lower buckets and is not a valid normalization.  Put

```text
P_0=ep'-3p,          H_0=eH'-2H,          G_0=eg'-2g.     (24)
```

After cancelling only nonzero global polynomial factors, the next two
buckets are the exact equations

```text
E_2:=
 -3rho b e^4-8rho e^2f'-(8rho/Delta)e^2M'
 +16rho ef+(16rho/Delta)eM+7ep'-21p=0,                   (25)

G_0=
 (21eH_0+3rho delta b e^4-8rho lambda e^2-15delta P_0)
 /(21Delta e).                                            (26)
```

Equation (26) is an identity in `k(e)`, obtained from a polynomial bucket by
cancelling the nonzero polynomial `21bDelta e^4`.  If its right side is not
polynomial, then no polynomial `g` exists and the branch is already empty;
otherwise it is exactly the Euler derivative of `g`.  No evaluation at
`e=0` and no pointwise division is used.

Define

```text
R=24rho e^8-16rho lambda e^7-8rho e^6
  +7rho e^5H'-15rho e^4H+12e^5M'-22e^4M,                 (27)

Q=-4e^2Mf'+4e^2fM'-3eHp'+5epH'-Hp.                      (28)
```

The complete terminal `[r^3]` bucket, after (26), is

```text
E_3:=R+Q/b-5delta e^3G_0=0.                              (29)
```

The exact companion reconstructs (25), (26), and (29) directly from the
unreduced bracket with

```text
U=rho b e^5,       S=b e^4,       N=delta b e^3.         (30)
```

In particular, it detects the three terms that a simultaneous `b=1`
normalization would incorrectly suppress.

## 4. The two target-color jets

Equation (16) allows the expansion

```text
f=1/12+f_1e+f_2e^2+O(e^3).                              (31)
```

Solve (23) for `M,H`:

```text
M=(D_mu f+1/12)/(2e),
H=(eta e^2-7M)/(8rho).                                  (32)
```

The numerator in the first expression is divisible by `e` because of (16).
The first three coefficients of (25) determine the coefficients of `p`
through degree two; its degree-three coefficient is the Euler resonance and
is left arbitrary.  Substituting into (29) gives

```text
[e]E_3=
 -(6f_1+lambda)(6f_1+mu)/(36b(lambda-mu)),               (33)

[e^2]E_3=
 -(48f_1^2lambda+48f_1^2mu+48f_1f_2+16f_1lambda mu
   -12f_1+4f_2lambda+4f_2mu-lambda-mu)
  /(32b(lambda-mu)).                                     (34)
```

Both expressions are independent of `delta`, `eta`, the Euler-resonant
coefficient of `p`, and every coefficient of `f` of degree at least three.
Since `b(lambda-mu)!=0`, (33) chooses one of the two target colors:

```text
f_1=-c/6,                 c in {lambda,mu}.              (35)
```

Then (34) fixes the second jet:

```text
f_2=c^2/3+1/4.                                           (36)
```

These are exactly the first two nonconstant coefficients of

```text
-1/(12D_c)=1/(12(1+2ce-3e^2)).                           (37)
```

Thus the terminal does not merely say that a coefficient is nonzero: it
tries to make the arm profile choose one of the two incompatible target
colors.

## 5. Every polynomial degree is terminal

Let `d=deg f`.  We keep the small degrees separate because the fixed source
monomial `24rho e^8` can then be top.

### Degree zero

Here `f_1=f_2=0`.  Equations (33)--(34) force

```text
lambda mu=0,                  lambda+mu=0.               (38)
```

Thus `lambda=mu=0`, contradicting the bichromatic hypothesis.

### Degree one

Equation (33) chooses a color `c`.  Equation (34) would additionally require
`4c^2+3=0`.  More decisively, solve the complete polynomial equation (25),
retain its free Euler coefficient, and substitute into the complete (29).
For either color,

```text
[e^8]E_3=24rho!=0.                                      (39)
```

This identity is independent of `b,delta,eta` and the free coefficient.

### Degree two

Use both forced jets (35)--(36), solve the complete (25), and again substitute
into (29).  For either color the exact terminal coefficient is

```text
[e^8]E_3=24rho!=0.                                      (40)
```

### Degree at least three

Write `F` for the nonzero leading coefficient of `f`.  The arm, cofactor,
and (25), successively, force

```text
M =(3/2)F e^(d+1)+lower,
H =-(21/(16rho))F e^(d+1)+lower,
p =(12rho/(7Delta))F e^(d+2)+lower.                      (41)
```

The term involving `G_0` in (29) has degree at most `d+4`.  The `M,f` part
of `Q` has degree at most `2d+2`; every term of `R` has degree below
`2d+3` when `d>=3`.  Only the `H,p` part of `Q/b` reaches degree `2d+3`, and
its coefficient is

```text
[e^(2d+3)]E_3=-9(d-1)F^2/(2bDelta)!=0.                  (42)
```

This closes every `d>=3`.  The proof is uniform in `delta`; hence it covers
both `N=0` and the nonzero `8/5` rung.

## 6. Boundary and next design target

The exact result is

```text
U,S!=0, lambda!=mu, v constant  ==>  no pair in (2).    (43)
```

It does **not** close:

```text
v nonconstant;
lambda=mu;
U=0 or S=0;
profiles beyond the displayed r^2z slot;
the general planar Jacobian conjecture.                  (44)
```

The live positive design problem is therefore not another constant monomial
tower.  It must use a nonconstant Kummer divisor, merge the two target
colors, enter one-sidedly, or move to a later canonical slot.  In the first
option, equations (20) and (29) show exactly what has to be synchronized:
the zero divisor of `v`, the factor-cofactor payment, and the color selected
by the arm jet.

## Reproduction

```bash
python3 04-computation/jc2_cubic_pseudoplane_constant_tower_bichromatic_r2z_thm3839.py
python3 -O 04-computation/jc2_cubic_pseudoplane_constant_tower_bichromatic_r2z_thm3839.py
```

Both executions must byte-match
`05-knowledge/results/jc2_cubic_pseudoplane_constant_tower_bichromatic_r2z_thm3839.out`.
