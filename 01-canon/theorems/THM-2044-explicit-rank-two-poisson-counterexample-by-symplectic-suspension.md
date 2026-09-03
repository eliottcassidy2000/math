---
id: THM-2044
title: "An explicit nonautomorphic rank-two Poisson endomorphism by symplectic suspension"
status: >
  PROVED. The four polynomials defined below satisfy the six canonical Poisson
  identities in C[x,q,p,z], and three explicitly displayed rational points form
  one exact fibre. The construction is a polynomial symplectic suspension of
  the THM-1300 three-dimensional Keller map. It disproves the Poisson
  Conjecture for two canonical pairs and gives a symplectic Keller
  counterexample in affine dimension four. It does not by itself disprove
  DC(2) or planar JC; see THM-2045 and HYP-8802.
source: codex-2026-07-21-DC2-JC2
related:
  - THM-1300
  - THM-2042
  - THM-2045
  - THM-2046
  - THM-4397
  - HYP-8802
  - HYP-8803
script: 04-computation/poisson_rank2_symplectic_suspension_codex_20260721.py
output: 05-knowledge/results/poisson_rank2_symplectic_suspension_codex_20260721.out
---

# THM-2044 -- rank-two Poisson counterexample

Use the convention

```text
{p,x}={z,q}=1
```

and all other generator brackets zero. Define `s=xq` and

```text
R  = x(2-3s),
D0 = ((1+3s)/2)p-3q^2z,
L  = 3x^2p+(2-6s)z.                                (1)
```

Put

```text
G(s)=252s^3+1008s^2+1379s+659,
g  = -q^2 G(s)/140,
ell=L+g.                                             (2)
```

In the following definitions `ell` is the polynomial (2), and

```text
y=q-x ell/3,             u=1+xy,
T=y+3xu^2 ell+3xy^2(4+3xy),
S=(u^3 ell+y^2u(4+3xy))/2.                          (3)
```

Finally set `D=D0+H`, where

```text
H=-ell B/1620,                                        (4)

B=2 ell^4 x^6(3s-2)
  +ell^3 x^4(-90s^2-30s+55)
  +ell^2 x^2(540s^3+720s^2-120s-270)
  +ell(-1620s^4-3780s^3-1890s^2+810s+540)
  +q^2(2430s^3+8100s^2+8640s+2430).
```

These are polynomials over `Q` (term counts after full expansion are
`2,35,246,78` for `R,T,D,S`). They obey

```text
{D,R}=1,                 {S,T}=1,
{R,T}={R,S}={D,T}={D,S}=0.                          (5)
```

Thus `(x,q,p,z) -> (R,T,D,S)` is a Poisson endomorphism.

## Why the construction works

Regard `ell` temporarily as an independent third variable and apply the shear

```text
(x,q,ell) -> (x, y=q-x ell/3, ell).                 (6)
```

Under (6), `(R,T,2S)` is exactly `(F3,F2,F1)` for the THM-1300 map

```text
F1=(1+xy)^3 ell+y^2(1+xy)(4+3xy),
F2=y+3x(1+xy)^2 ell+3xy^2(4+3xy),
F3=2x-3x^2y-x^3 ell.
```

Consequently

```text
det d(R,T,S)/d(x,q,ell)=1.                          (7)
```

For any polynomials `A,B` in `(x,q,ell)`, substitution of
`ell=L+g(x,q)` gives the Nambu identity

```text
{A,B}=-det d(R,A,B)/d(x,q,ell).                     (8)
```

Equation (7) therefore gives `{S,T}=1` and the two `R` cross-brackets.
The correction (4) is the polynomial connection making the vector field
`{D,-}` equal to differentiation in the `R` direction while holding `T,S`
fixed. Direct expansion gives the remaining three identities in (5).

The apparently arbitrary polynomial `g` kills the only base cohomology
remainder. Before shifting `L`, that remainder is

```text
q^3(54s^3+189s^2+222s+89)/2.
```

Writing `g=q^2 f(s)`, its cancellation is the one-variable identity

```text
(1-3s)f'(s)-21f(s)=54s^3+189s^2+222s+89,           (9)
```

whose unique cubic solution is `f=-G/140`. This is the promised reduction of
the symplectic identity to coefficient identities.

## Exact three-point fibre and nonautomorphy

The three distinct points

```text
P0=(0, 0, 1/24, -1/8),
P+=(1, 2/3, 224839/90720, -173417/60480),
P-=(-1, -2/3, 224839/90720, -173417/60480)          (10)
```

all map to

```text
(R,T,D,S)=(0,0,0,-1/8).                             (11)
```

This fibre contains exactly these three points. Indeed the source change

```text
(x,q,p,z) -> (x,q,d=D,ell=L+g)                      (12)
```

is a polynomial automorphism. At fixed `(x,q)`, the linear map
`(p,z)->(D0,L)` has determinant one, and the subsequent `g,H` changes are
triangular. In the coordinates (12), the map is simply

```text
(x,q,d,ell) -> (R(x,q),T(x,q,ell),d,S(x,q,ell)),    (13)
```

and (6) plus the target permutation/scaling identifies its three-dimensional
factor with THM-1300. The exactly-three fibre of THM-1300 therefore transports
to (10)-(11). In particular the Poisson endomorphism is not an automorphism.

Because (5) says its Jacobian matrix is symplectic, its ordinary Jacobian
determinant is one. Hence it is also a four-variable Jacobian counterexample.

## Scope

The result is rank two in the standard Poisson indexing: two canonical pairs,
four commutative generators. Direct replacement of brackets by Weyl
commutators is invalid; HYP-8802 records the nonzero ordering anomaly. The safe
Hamiltonian-dual/cotangent construction doubles to four Weyl pairs and lands in
conventional `A_4`, not `A_2`.

THM-2046 now also excludes the lower-dimensional first-order cotangent
pullback shortcut: with `R` as a multiplication position, a second
multiplication position and first-order dual momenta would force a planar
Keller mate, contradicting THM-2045. HYP-8803 identifies the surviving
nonfiltered problem as polynomial extension across the `x=0` divisor of an
exact localized Ore-Weyl chart.

THM-4397 proves that Long's later arXiv:2608.23777v1 presentation is exactly
this map up to an explicit polynomial Hamiltonian source translation and a
linear symplectic target quarter-turn.  The smaller presentation is useful for
computation, but it is the same right-left equivalence class and does not alter
the planar or `A_2` scope boundary.
