---
id: THM-2045
title: "The smooth factorized family R=x(a-b x^r q^s) has no planar Jacobian mate"
status: >
  PROVED. For nonzero a,b and positive integers r,s, no polynomial Q in C[x,q]
  has Jacobian bracket {Q,R}=1 with R=x(a-b x^r q^s). In particular the first
  coordinate of THM-2044 cannot be de-stabilized to a planar Keller pair. The
  proof separates weighted Laurent sectors and shows the leading coefficient
  in the only constant-producing sector cannot vanish.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-1345
  - THM-2042
  - THM-2044
  - THM-2046
  - HYP-8802
  - HYP-8803
---

# THM-2045 -- no planar mate for the suspension coordinate

Let `a,b` be nonzero complex numbers, let `r,s` be positive integers, and put

```text
R=x(a-b x^r q^s).                                  (1)
```

There is no `Q in C[x,q]` satisfying

```text
R_x Q_q-R_q Q_x=1.                                  (2)
```

Write `g=gcd(r,s)`, `r=g r0`, `s=g s0`, and

```text
v=x^r0 q^s0,                 x^r q^s=v^g.           (3)
```

Grade a monomial `x^i q^j` by

```text
kappa=s i-r j.                                        (4)
```

Every term of the Hamiltonian derivation

```text
W=R_x partial_q-R_q partial_x                       (5)
```

raises `kappa` by exactly `r`. Consequently only the `kappa=-r` sector of `Q`
can contribute the constant in (2). The nonnegative integer solutions of
`s i-r j=-r` are

```text
(i,j)=(r0 t,1+s0 t),            t>=0,
```

so this sector has the form

```text
q f(v),                       f in C[v].             (6)
```

Direct differentiation gives

```text
W(v)=s0 (v/q)(a-bv^g),

W(q f(v))=[a-b(r+1)v^g]f(v)
           +s0 v(a-bv^g)f'(v).                     (7)
```

If `f` has degree `N` and leading coefficient `c_N!=0`, the coefficient of
`v^(N+g)` in (7) is

```text
-b (r+1+s0 N)c_N,                                  (8)
```

which is nonzero. Thus (7) cannot equal the constant one. No other weighted
sector can affect the constant term, proving the theorem.

For `(a,b,r,s)=(2,3,1,1)`, this applies to the `R` in THM-2044. The rank-two
Poisson counterexample therefore uses genuine symplectic stabilization:
deleting the extra canonical pair cannot leave a planar Jacobian counterexample
with this first coordinate. The wider family also supplies a two-parameter
Newton-edge exclusion stratum for planar JC.

THM-2046 upgrades this de-stabilization obstruction on the quantum side: the
same family cannot be the first multiplication-position image of any `A_2`
endomorphism whose second position is also a multiplication polynomial and
whose dual generators have differential order at most one.
