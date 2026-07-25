---
id: THM-2180
title: "The first cubic pole stratum is empty for twice-odd Keller quartics"
status: >
  PROVED. For a planar Keller pair whose first component is
  V^2 z^4+beta z^3+gamma z^2+delta z+epsilon, every reduced mate degree
  n not divisible by four forces V|beta. At a hypothetical place with
  nu(beta)<nu(V), quadratic-deck depression has the normalized polar face
  Z^4-6Z^2+8Z-3=(Z-1)^3(Z+3). The top reduced Faber term is uniquely deepest
  on the original fibre, but its boundary value is
  4^n binom(n/4,n), which is nonzero. Thus the first pole stratum is empty.
  The linear coefficient of the canonical approximate root is polynomial;
  the remaining open condition is V|(4gamma-(beta/V)^2).
source: codex-2026-07-24-JC-quartic-pole-descent
depends_on:
  - THM-2129
  - THM-2158
related:
  - THM-2129
  - THM-2158
  - THM-2181
---

# THM-2180 -- the first quartic pole stratum is empty

Let

```text
P=V^2 z^4+beta z^3+gamma z^2+delta z+epsilon
```

belong to a planar polynomial Keller pair over `C`, with `V!=0`. Reduce the
mate by polynomial target shears, and let `n>0` be its remaining `z`-degree.
Assume

```text
4 does not divide n.                                  (1)
```

The twice-odd branch relevant to THM-2129 has `n=2 mod 4`, but the argument
uses only (1).

Then

```text
V divides beta.                                       (2)
```

## 1. The normalized polar face

Work over the quadratic extension of `C(x)` obtained by choosing `U^2=V`.
Put

```text
w=Uz,                 h=beta/(4U^3),
B_2=gamma/U^2,        B_1=delta/U,        B_0=epsilon,
Z=w+h.                                                   (3)
```

Direct depression gives

```text
P=Z^4+pZ^2+qZ+r,                                      (4)

p=B_2-6h^2,
q=B_1-2B_2h+8h^3,
r=B_0-B_1h+B_2h^2-3h^4.                              (5)
```

Fix an irreducible `pi|V`, extend its valuation to the quadratic field, and
write

```text
e=nu_pi(V)>0,                 b=nu_pi(beta).           (6)
```

Suppose for contradiction that `b<e`. Then

```text
nu_pi(h)=b-3e/2<0.                                    (7)
```

Since `beta,gamma,delta,epsilon` are polynomial,

```text
nu_pi(B_2/h^2)>=2e-2b>0,
nu_pi(B_1/h^3)>=4e-3b>0,
nu_pi(B_0/h^4)>=6e-4b>0.                              (8)
```

Consequently the normalized coefficients of (4) have the exact residue

```text
(p/h^2,q/h^3,r/h^4) -> (-6,8,-3).                    (9)
```

The hostile face is

```text
T^4-6T^2+8T-3=(T-1)^3(T+3).                          (10)
```

It is not the balanced square face `(T^2-1)^2`; this distinction is the
mechanism which kills the first pole stratum.

## 2. Reduced Faber normal form

For every `m>=0`, put

```text
E_m=Pol_Z(P^(m/4)).                                   (11)
```

THM-2129 proves that the Hamiltonian derivative of `E_m` has `Z`-degree at
most two. This also gives the usual triangular normal form for the mate:

```text
Q=H(P)+sum_(4 does not divide m, m<=n) c_m E_m,
H in C[T],                 c_m in C,                 (12)
```

with `c_n!=0`.

For completeness, the constants in (12) are forced, not assumed. If a
remaining summand has top term `a_d(x)Z^d`, then the coefficient of
`Z^(d+3)` in its Hamiltonian derivative is `-4a_d'(x)`. The Keller bracket
has `Z`-degree zero, so `a_d'=0`. Subtract the corresponding constant
multiple of the monic `E_d` and descend in `d`. Terms with `4|d` are powers
of `P` and collect into `H`.

## 3. The top boundary value cannot vanish

On the original fibre `z=0`, equation (3) gives `Z=h`. Both `Q(x,0)` and
`H(P(x,0))=H(epsilon(x))` are regular at `pi`. Weighted homogeneity gives

```text
E_m(h;p,q,r)
 =h^m E_m(1;p/h^2,q/h^3,r/h^4).                     (13)
```

Because `nu_pi(h)<0`, every term with `m<n` is strictly shallower than the
top `m=n` term. Regularity of (12) therefore forces the leading residue

```text
E_n(1;-6,8,-3)=0.                                    (14)
```

But (10), translated by `t=T-1`, becomes

```text
P_0(t+1)=t^4+4t^3.                                   (15)
```

THM-2129's exact boundary formula now gives

```text
E_n(1;-6,8,-3)
 =[u^n](1+4u)^(n/4)
 =4^n binom(n/4,n).                                  (16)
```

The generalized binomial coefficient can vanish only when `n/4` is an
integer in `{0,...,n-1}`. Condition (1) excludes this, so (16) is nonzero,
contradicting (14).

Thus `nu_pi(beta)>=nu_pi(V)` at every `pi|V`, which proves (2).

## 4. Exact remaining pole

THM-2158's canonical approximate root is

```text
H_0
 =Vz^2+[beta/(2V)]z
   +[4gamma V^2-beta^2]/(8V^3).                      (17)
```

Write `beta=V beta_1`, now justified by (2). Then

```text
H_0=Vz^2+(beta_1/2)z+[4gamma-beta_1^2]/(8V).         (18)
```

Hence the Keller equations have already made the linear coefficient
polynomial and reduced every remaining finite pole to depth at most
`nu_pi(V)`. Polynomiality of `H_0` is now exactly the single congruence

```text
V divides 4gamma-beta_1^2.                            (19)
```

This theorem does not prove (19). The balanced square-resonant stratum in
which (19) fails is the precise remaining nonmonic quartic pole problem.
QED.
