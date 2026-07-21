---
id: THM-2045
title: "The smooth factorized family R=x(a-bxq) has no planar Jacobian mate"
status: >
  PROVED. For nonzero a,b, no polynomial Q in C[x,q] has Jacobian bracket
  {Q,R}=1 with R=x(a-bxq). In particular the first coordinate of THM-2044
  cannot be de-stabilized to a planar Keller pair. The proof separates Laurent
  diagonal sectors and reduces the only possible constant-producing sector to
  a one-line differential equation.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-1345
  - THM-2042
  - THM-2044
  - HYP-8802
---

# THM-2045 -- no planar mate for the suspension coordinate

Let `a,b` be nonzero complex numbers and

```text
R=x(a-bxq).
```

There is no `Q in C[x,q]` satisfying

```text
R_x Q_q-R_q Q_x=1.                                  (1)
```

To prove this, localize at `x` and write `s=xq`. The Hamiltonian derivation in
(1) becomes

```text
W=R_x partial_q-R_q partial_x
 =x[b x partial_x+(a-bs)partial_s].                  (2)
```

Every monomial `x^i q^j` is `x^(i-j)s^j`. Thus the Laurent `x`-exponent
`k=i-j` grades `C[x,q]`, and (2) sends the `k`-sector to the `(k+1)`-sector.
Only the `k=-1` sector can contribute the constant on the right side of (1).
Its polynomial members have the form

```text
x^(-1)f(s),             with f(s) divisible by s.   (3)
```

Equation (1) on that sector is

```text
(a-bs)f'(s)-bf(s)=1,
```

or equivalently

```text
((a-bs)f(s))'=1.                                    (4)
```

Hence `(a-bs)f=s+C`. Polynomiality forces `C=-a/b` and `f=-1/b`, contradicting
the divisibility by `s` required in (3). No other sector can affect the
constant term, proving the theorem.

For `(a,b)=(2,3)`, this applies to the `R` in THM-2044. The rank-two Poisson
counterexample therefore uses genuine symplectic stabilization: deleting the
extra canonical pair cannot leave a planar Jacobian counterexample with this
first coordinate.
