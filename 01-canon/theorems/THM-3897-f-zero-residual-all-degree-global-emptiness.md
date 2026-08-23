---
id: THM-3897
title: "F-zero residual all-degree global emptiness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Using the
  now-proved THM-3895 high-y-degree
  input, the complete polynomial f=0 residual lane consists only of the
  base point T=0.  The quadratic-y
  channel collapses by a polynomial Pell/unit argument, the linear-y
  channel by odd degree, and the constant-y channel by its missing linear
  coefficient together with the origin address.  The f!=0 residual lane,
  a polynomial-plane Keller atlas, and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3895 three-channel descent, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (2026-08-23).  A structurally separate
  companion rederives the THM-3895 quartic covariant, 1,596-row degree
  ledger, four exact square ideals, and osculating a^(-3) pole before
  auditing the present three-channel descent.  It rechecks UFD parity,
  the unit-factor polynomial Pell argument, every zero-leading alternative,
  the sole use of T(0,0)=0, and the hostile rational k(x)[y] point in 1,719
  active gates.  Normal and optimized runs are byte-identical.  The primary
  companion checks every residual
  leading coefficient, the Pell factorization, the constant-channel origin
  value, and the rational-x hostile point in 15 active gates.  Normal and
  optimized runs byte-match the frozen output.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
  - THM-3895-f-zero-quartic-covariant-and-high-y-degree-emptiness
related:
  - THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality
script: 04-computation/jc2_f_zero_all_degree_global_emptiness_thm3897.py
output: 05-knowledge/results/jc2_f_zero_all_degree_global_emptiness_thm3897.out
script_sha256: c0a15aebe8542543c5cbad24f31d4bb3635a7a1ade31af1bd695619d554b7094
output_sha256: af62022edb2bf93bee90dfc5676f310b457eea0c29ed72f14e873b55ab47648b
semantic_sha256: 6d9be0bbf8c5242bf061d89a4bf2d18c42f81e373008006880946b07b4531240
independent_audit_script: 04-computation/jc2_f_zero_all_degree_global_emptiness_independent_audit_thm3897.py
independent_audit_output: 05-knowledge/results/jc2_f_zero_all_degree_global_emptiness_independent_audit_thm3897.out
independent_audit_script_sha256: 0762969b2d2df3720f9ac9f27d0d5d7e730e05910f3cb766e03c5976e462c71b
independent_audit_output_sha256: 69b9439abc5b6f4612bacf34c28f52d58b40ed1793c20158c33791770aae2035
independent_audit_semantic_sha256: 79d829a2562577865b20d68f4ab1151bbf1656d3e7ef5ac25e0a345ada94b693
hash_basis: raw LF bytes
---

# THM-3897 -- the polynomial f-zero residual has only its base point

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Put

```text
D=k[x,y],                 a=x+1,
L=9x+4,                   F=15x^2+15x+4,
K=y^2-F.                                                    (1)
```

Suppose `T,G in D` obey the `f=0` residual equation and its source address

```text
G^2=L^4-6aL^2T^2-8KT^3-3a^2T^4,
T(0,0)=0.                                                    (2)
```

Then

```text
T=0,                       G=+/-L^2.                        (3)
```

The proved THM-3895 reduces every polynomial solution of the square equation
in `(2)` to `deg_y(T)<=2`.  The present theorem is the final global descent
from those three remaining channels.

## 1. Inheritance and the exact remaining universe

THM-3885 supplies both equation `(2)` and its address.  THM-3888 identifies
the generic `k(x)[y]` curve with an equianharmonic elliptic surface and gives
the hostile rational-x point

```text
T_*=-2K/(3a^2),
G_*=4K^2/(3a^3)-L^2.                                      (4)
```

It satisfies the square equation exactly, but is not in `D` and also has
`T_*(0,0)=8/3`.  Thus both polynomial descent and the address in `(2)` are
essential.

By THM-3895, it remains only to write the complete three-channel universe

```text
T=u(x)y^2+v(x)y+w(x),             u,v,w in k[x].           (5)
```

The three channels die for different structural reasons.

## 2. Quadratic channel: a polynomial Pell equation has only the zero rung

The coefficient of `y^8` on the right side of `(2)` is

```text
-u^3(8+3a^2u).                                              (6)
```

If `u!=0`, this coefficient cannot vanish: its only formal alternative is
`u=-8/(3a^2)`, which is not a polynomial (equivalently, evaluation of
`3a^2u=-8` at `x=-1` is impossible).  Hence `(6)` is a nonzero square in
`k[x]`.

Now

```text
gcd(u,8+3a^2u)=1                                           (7)
```

because `8` is a unit.  Unique factorization and parity of valuations in
`(6)` imply, after absorbing nonzero constants into squares (the field is
algebraically closed), that

```text
u=r^2,                    s^2=8+3a^2r^2                   (8)
```

for some `r,s in k[x]`.  Choose `sqrt(3) in k`.  Equation `(8)` factors as

```text
(s-sqrt(3)ar)(s+sqrt(3)ar)=8.                              (9)
```

Both factors are units of `k[x]`; subtracting them shows that `ar` is a
constant.  Since `a=x+1` is nonconstant, this forces `r=0`, contrary to
`u!=0`.  Therefore

```text
u=0.                                                       (10)
```

This is the global denominator obstruction hidden by the generic point
`(4)`: over `k(x)`, `a` is a unit and the conclusion does not follow.

## 3. Linear channel: the residual has odd degree

With `u=0`, suppose `v!=0`.  The term `-8KT^3` has the unique highest
`y`-degree, namely five, and

```text
[y^5] RHS(2)=-8v^3!=0.                                    (11)
```

A nonzero polynomial square cannot have odd degree in `y`.  Hence

```text
v=0.                                                       (12)
```

## 4. Constant channel: the missing coefficient meets the address

It remains that `T=w(x)`.  If `w!=0`, the right side of `(2)` has degree two
in `y`, with

```text
[y^2] RHS(2)=-8w^3!=0,             [y] RHS(2)=0.           (13)
```

Write its square root as `p(x)y+q(x)`.  The quadratic coefficient makes
`p!=0`, while the missing linear coefficient gives `2pq=0`, hence `q=0`.
The constant coefficient must therefore vanish identically.  But the
address gives `w(0)=0`, and at `x=0` that constant coefficient is

```text
L(0)^4=4^4=256,                                            (14)
```

a contradiction.  Thus `w=0`.  Substitution into `(2)` gives
`G^2=L^4`, proving `(3)`.

## 5. Exact consequence and open boundary

The proof closes the **entire polynomial `f=0` lane**, not merely a bounded
degree cell.  In the THM-3881 residual coordinates, any counterexample-positive
square lift must therefore pay a genuinely nonzero sidecar `f`; another pass
through higher-degree `T` with `f=0` cannot help.

Nothing here classifies the `f!=0` residual equation, constructs a Keller
atlas, or proves `JC(2)`.  The next live problem is to understand whether
`f` can cancel the same quadratic Pell obstruction without creating a new
arm or denominator debt.

## 6. Reproduction

```bash
python3 04-computation/jc2_f_zero_all_degree_global_emptiness_thm3897.py
python3 -O 04-computation/jc2_f_zero_all_degree_global_emptiness_thm3897.py
python3 04-computation/jc2_f_zero_all_degree_global_emptiness_independent_audit_thm3897.py
python3 -O 04-computation/jc2_f_zero_all_degree_global_emptiness_independent_audit_thm3897.py
```

Each normal/optimized pair must byte-match its corresponding frozen output
in `05-knowledge/results/`.  The primary companion checks 15 active
identities and boundary controls; the independent companion checks 1,719
active gates and separately rebuilds the inherited high-y cutoff.  Their
exact algebra supports the displayed formulas; the UFD and unit arguments
above are the proof.
