---
id: THM-3885
title: "Cusp-residual f=0 arm dichotomy and cubic closure"
status: >
  PROVED + VERIFIED-EXACT; AUDITED BASE + CUBIC STRENGTHENING AWAITING
  FOCUSED HOSTILE AUDIT.  In the THM-3881 rank-two residual equation, the
  full f=0 lane has an exact arm dichotomy: T(-1,y) is either zero or a
  constant c with c^3=-625/32.  Every total-degree-at-most-three square pair
  in this lane is the base pair `T=0`.  At the independent L=0 arm, all
  surviving square classes have an exact finite root-polarization grammar.
  Nonlinear T=c+(x+1)U with deg U>=3, a Keller atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3881 f=0 lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).  The
  audit rederived the Mason common-degree/radical bound, exhausted every
  constant edge with both positive controls, and checked the odd-degree and
  missing-`y` square-root arguments in the quadratic cell.  It independently
  proved both UFD-parity directions and every unit/zero boundary in the
  `L=0` root-polarization grammar.  Here “empty” in the transcript means no
  nonzero `T`; the base pair remains.  CUBIC STRENGTHENING AWAITING FOCUSED
  HOSTILE AUDIT: the new coefficient recursion, the `p=0` seam, the reduced
  Groebner basis `<z,q^2>`, the affine `L=0` odd-degree gate, and the final
  one-color reduction must be independently rechecked.  The amended companion
  verifies these gates in addition to the audited base packet.  Normal and
  optimized runs must byte-match the frozen output.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
related:
  - THM-3872-three-cusp-polarization-branches-and-minimal-affine-square-residual-gate
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
script: 04-computation/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.py
output: 05-knowledge/results/jc2_cusp_residual_f_zero_arm_quadratic_thm3885.out
script_sha256: ec653744f276da842904249e93819f1f96ac870bcee9e578c143a6b089de056e
output_sha256: a6c1af7024c6155abb649e23c169826e419f58f3d990e8165f4fc794e2b7eca0
semantic_sha256: fa6375d39828c8eaf658b45830d8b4285472fab7043a03a1cbbab64ff35affd0
hash_basis: raw LF bytes
---

# THM-3885 -- the f-zero lane starts above cubic degree

**PROVED + VERIFIED-EXACT; AUDITED BASE + CUBIC STRENGTHENING AWAITING
FOCUSED HOSTILE AUDIT.**
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

It closes all total degrees at most three and gives two all-degree boundary
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

## 3. Every cubic cell has only the base point

Assume `deg T<=3`.  By `(10)`, the quotient has degree at most two.  Write its
complete coefficient universe as

```text
U=alpha x^2+beta xy+q y^2+gamma x+p y-c,
T=c+aU.                                                     (11)
```

The last constant is exactly the address `(2)`.  At `x=0`, put

```text
tau(y)=T(0,y)=p y+q y^2.                                   (12)
```

The specialized residual is

```text
S_0=256-96tau^2-8(y^2-4)tau^3-3tau^4.                     (13)
```

Suppose `S_0=G^2`.  Since `S_0(0)=256`, replace `G` by `-G` if necessary and
write

```text
G=16+g_2y^2+g_3y^3+g_4y^4.                                (14)
```

There is no `y` term because `[y]S_0=0`.  Successive comparison in degrees
two, three, and four uniquely gives

```text
g_2=-3p^2,
g_3=p(p^2-6q),
g_4=-(3/8)(p^4-8p^2q+8q^2).                               (15)
```

Put `z=p^2`.  The next three residual coefficients are

```text
[y^5](G^2-S_0)=-2p E_5,
[y^6](G^2-S_0)=E_6/4,
[y^7](G^2-S_0)=-(3/4)p E_7,                               (16)

E_5=3z^2-24zq-4z+48q^2,
E_6=13z^3-120z^2q+288zq^2+96zq-128q^3,
E_7=z^3-14z^2q+56zq^2-64q^3-32q^2.                       (17)
```

If `p=0`, equation `[y^6]=0` reads `-32q^3=0`, so `q=0`.  If `p!=0`, all
three `E_i` vanish.  Exact division over `Q[z,q]` gives the reduced basis

```text
Groebner(E_5,E_6,E_7; z,q, lex)=[z,q^2],                  (18)
```

contradicting `z=p^2!=0`.  Therefore `p=q=0` uniformly.

It remains to remove the mixed coefficient.  At `L=0`, equation `(11)` is

```text
tau_L(y)=A+By,
A=4c/9+80alpha/729-20gamma/81,
B=-20beta/81.                                              (19)
```

By the exact arm identity of Section 4,

```text
S_L=-tau_L^3(8K_0+(25/27)tau_L).                          (20)
```

If `B!=0`, this has degree five and leading coefficient `-8B^3`, so it is
not a square.  Thus `beta=0`.  For completeness, when `B=0` the arm is a
square exactly at

```text
A=0  or  A=64/25.                                         (21)
```

We have now forced `T=T(x)`.  Globally `(3)` has `y`-degree two, no linear
`y` coefficient, and

```text
[y^2]S=-8T(x)^3,             S(0,0)=256.                  (22)
```

If it were a square, its root would have the form `U_0(x)+yU_1(x)`.  The
missing linear term says `U_0U_1=0`; the nonzero origin value forces
`U_0!=0`, so `U_1=0`.  Equation `(22)` then forces `T=0`.  Hence

```text
f=0, deg T<=3, and S(T,0) square  ==>  T=0.               (23)
```

This includes the affine constant-span cells of THM-3872 and closes both the
quadratic layer and the first genuinely nonlinear cubic layer.

## 4. Exact root polarization at `L=0`

There is a second, independent boundary grammar.  Set

```text
x=-4/9,                 a=5/9,
K_0=y^2-8/27,           tau(y)=T(-4/9,y).                 (24)
```

Then

```text
S=-tau^3(8K_0+(25/27)tau).                                (25)
```

The zero polynomial `tau=0` is one solution.  Suppose `tau!=0`.  Normalize

```text
d=gcd(tau,K_0),          tau=d sigma,
K_0=de,                  gcd(sigma,e)=1.                  (26)
```

Because `K_0` is squarefree, `d,e` are coprime squarefree complementary
divisors.  Equation `(25)` is a square if and only if

```text
sigma=u^2,
8e+(25/27)u^2=v^2                                       (27)
```

for some `u,v in k[y]`.  Indeed, after removing the square `d^4`, the two
factors `sigma^3` and `8e+(25/27)sigma` are coprime; UFD parity proves the
forward direction, and the reverse direction is immediate because `-1` is a
square in `k`.

Put

```text
gamma=5/(3sqrt(3)).                                       (28)
```

The factors `v-gamma u` and `v+gamma u` are coprime: a common factor would
divide `u,v`, then `e`, contrary to `gcd(u,e)=1`.  Thus `(27)` is equivalent
to an ordered root partition

```text
e=e_-e_+,                 gcd(e_-,e_+)=1,
v-gamma u=lambda e_-,
v+gamma u=8lambda^(-1)e_+,                               (29)
```

with `lambda in k^*`.  Explicitly,

```text
u=(8lambda^(-1)e_+-lambda e_-)/(2gamma),
tau=d u^2.                                                (30)
```

Conversely every choice `(26),(29),(30)` makes `(25)` a square.  Since `K_0`
has only two simple roots, `(26),(29)` give finitely many root-assignment
shapes; only the nonzero scale `lambda` remains continuous.

## 5. Exact surviving problem

Any nonzero `f=0` square survivor must now have

```text
T=c+aU,                   c=0 or c^3=-625/32,
deg U>=3,                 U(0,0)=-c,                     (31)
```

and its second-arm restriction must obey the polarization `(26)-(30)`.
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
