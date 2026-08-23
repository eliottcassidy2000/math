---
id: THM-3885
title: "Cusp-residual f=0 arm dichotomy and quadratic closure"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.  In
  the THM-3881 rank-two residual equation, the full f=0 lane has an exact arm
  dichotomy: T(-1,y) is either zero or a constant c with c^3=-625/32.  Every
  total-degree-at-most-two pair in this lane is empty.  At the independent
  L=0 arm, all surviving square classes have an exact finite root-polarization
  grammar.  Nonlinear T=c+(x+1)U with deg U>=2, a Keller atlas, and JC(2)
  remain OPEN.
source: jc_zero_debt_lift / post-THM-3881 f=0 lane, 2026-08-23
audit: >
  PROVISIONAL EXACT PROOF CANDIDATE.  The companion verifies the f=0 residual,
  both arm specializations, the complete constant square Groebner basis and
  both positive controls, the quadratic address/odd-degree/one-color gates,
  and the L=0 gcd extraction and root polarization in 20 active gates.  Normal
  and optimized runs must byte-match the frozen output.  Independent audit
  must recheck the Mason degree bound, every constant boundary, the y-degree
  square argument, and both directions of the L=0 classification.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
related:
  - THM-3872-three-cusp-polarization-branches-and-minimal-affine-square-residual-gate
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
script: 04-computation/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.py
output: 05-knowledge/results/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.out
script_sha256: 004af773747738f3524935d076bdfebf5c07544b27a80b9021f52dd6ead74dab
output_sha256: 1701d8a3fd4f850036651953c0cf74628e2b21f70241ec07996698832441374e
semantic_sha256: 53ee84885a6a2415e44aca8523d1d2269eb1946016c58307ad151bcf2c98b1bf
hash_basis: raw LF bytes
---

# THM-3885 -- the f-zero lane starts above quadratic degree

**PROVED + VERIFIED-EXACT CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT.**
Work over an algebraically closed field `k` of characteristic zero and retain
the notation of THM-3881:

```text
D=k[x,y],                 a=x+1,
L=9x+4,                   K=y^2-15x^2-15x-4.             (1)
```

This theorem studies the exact subspace `f=0` of the rank-two pair `(T,f)`.
The address law becomes

```text
T(0,0)=0.                                                   (2)
```

It closes all total degrees at most two and gives two all-degree boundary
classifications.  It does not close the remaining nonlinear interpolation
problem.

## 1. Exact residual on `f=0`

Substitute `f=0` into THM-3881 equation (21).  Then

```text
r=aT,                       A=KT,
S(T,0)=L^4-6aL^2T^2-8KT^3-3a^2T^4.                      (3)
```

All statements below are necessary conditions for `(3)` to be a square in
`D`; the final `L=0` statement is an exact classification on that arm.

## 2. The `a=0` arm has only four constant addresses

Set `a=0`, equivalently `x=-1`, and write

```text
tau(y)=T(-1,y).
```

If `(3)` is a square, then for some `G in k[y]`,

```text
G^2+8(y^2-4)tau^3=625.                                    (4)
```

Suppose `tau` is nonconstant and put `n=deg tau>=1`.  The two nonconstant
summands in `(4)` are coprime because their sum is the nonzero constant `625`.
Their common degree is

```text
D_0=2+3n=2deg G;                                           (5)
```

if `D_0` is odd, this is already impossible.  Otherwise Mason--Stothers gives

```text
2+3n <= deg rad(G(y)(y^2-4)tau(y))-1
       <= (2+3n)/2+n+1.                                   (6)
```

Thus `n<=0`, a contradiction.  Therefore `tau=c` is constant.  Equation `(4)`
becomes

```text
G^2=625+32c^3-8c^3y^2.                                   (7)
```

For `c!=0`, a square root must be linear.  Its missing linear `y` coefficient
and nonzero quadratic coefficient force its constant term to vanish, hence
`32c^3+625=0`.  Conversely both possibilities are squares:

```text
c=0:                    G^2=25^2,
c^3=-625/32:            G^2=(25y/2)^2.                   (8)
```

Consequently every global square survivor satisfies exactly one of

```text
T(-1,y)=0,
T(-1,y)=c,                 c^3=-625/32.                   (9)
```

Equivalently,

```text
T=c+aU,                    U(0,0)=-c,                     (10)
```

where `c` is one of the four values in `(8)`.

## 3. Every quadratic cell is empty

Assume `deg T<=2`.  By `(10)`, the quotient has degree at most one, so write

```text
T=c+a(px+qy-c).                                             (11)
```

The last constant is forced by `(2)`.  At `x=0`, equation `(11)` gives
`T=qy`, and the specialization of `(3)` has degree at most five with

```text
[y^5]S(0,y)=-8q^3.                                        (12)
```

If `q!=0`, the specialization has odd degree five and is not a square.  Hence
`q=0`, so `T` depends only on `x`.  Now `(3)` has `y`-degree two, no linear
`y` coefficient, and

```text
[y^2]S=-8T(x)^3,             S(0,0)=256.                  (13)
```

If it were a square, its root would have the form `U_0(x)+yU_1(x)`.  The
missing linear term says `U_0U_1=0`; the nonzero constant part forces
`U_0!=0`, so `U_1=0`.  Equation `(13)` then forces `T=0`.

Thus

```text
f=0, deg T<=2, and S(T,0) square  ==>  T=0.               (14)
```

This includes the affine constant-span cells of THM-3872 and closes the first
new quadratic layer.

## 4. Exact root polarization at `L=0`

There is a second, independent boundary grammar.  Set

```text
x=-4/9,                 a=5/9,
K_0=y^2-8/27,           tau(y)=T(-4/9,y).                 (15)
```

Then

```text
S=-tau^3(8K_0+(25/27)tau).                                (16)
```

The zero polynomial `tau=0` is one solution.  Suppose `tau!=0`.  Normalize

```text
d=gcd(tau,K_0),          tau=d sigma,
K_0=de,                  gcd(sigma,e)=1.                  (17)
```

Because `K_0` is squarefree, `d,e` are coprime squarefree complementary
divisors.  Equation `(16)` is a square if and only if

```text
sigma=u^2,
8e+(25/27)u^2=v^2                                       (18)
```

for some `u,v in k[y]`.  Indeed, after removing the square `d^4`, the two
factors `sigma^3` and `8e+(25/27)sigma` are coprime; UFD parity proves the
forward direction, and the reverse direction is immediate because `-1` is a
square in `k`.

Put

```text
gamma=5/(3sqrt(3)).                                       (19)
```

The factors `v-gamma u` and `v+gamma u` are coprime: a common factor would
divide `u,v`, then `e`, contrary to `gcd(u,e)=1`.  Thus `(18)` is equivalent
to an ordered root partition

```text
e=e_-e_+,                 gcd(e_-,e_+)=1,
v-gamma u=lambda e_-,
v+gamma u=8lambda^(-1)e_+,                               (20)
```

with `lambda in k^*`.  Explicitly,

```text
u=(8lambda^(-1)e_+-lambda e_-)/(2gamma),
tau=d u^2.                                                (21)
```

Conversely every choice `(17),(20),(21)` makes `(16)` a square.  Since `K_0`
has only two simple roots, `(17),(20)` give finitely many root-assignment
shapes; only the nonzero scale `lambda` remains continuous.

## 5. Exact surviving problem

Any nonzero `f=0` square survivor must now have

```text
T=c+aU,                   c=0 or c^3=-625/32,
deg U>=2,                 U(0,0)=-c,                     (22)
```

and its second-arm restriction must obey the polarization `(17)-(21)`.
These conditions are necessary, not sufficient.  The nonlinear interpolation
between the two arms, the general `T,f!=0` equation, a polynomial-plane Keller
atlas, and JC(2) remain **OPEN**.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.py
python3 -O 04-computation/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.py
```

Both runs must byte-match
`05-knowledge/results/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.out`.
