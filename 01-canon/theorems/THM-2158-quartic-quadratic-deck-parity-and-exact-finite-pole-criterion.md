---
id: THM-2158
title: "Quartic quadratic-deck parity and the exact finite-pole criterion"
status: >
  PROVED. Let a quartic over C[x] have leading coefficient V^2. After the
  customary quadratic base change U^2=V, monicization, and depression, the
  deck involution sends the depressed variable Z to -Z, fixes the even
  depressed coefficients, and negates the odd coefficient. Consequently the
  canonical square approximate root Z^2+p/2 and its linear remainder are
  deck invariant and already lie in C(x)[z]. Thus quadratic-deck descent is
  automatic; the remaining obstruction is finite-pole regularity. The
  approximate root is polynomial exactly when V divides the cubic
  coefficient beta and, writing beta=V beta_1, V divides
  4 gamma-beta_1^2, where gamma is the quadratic coefficient. Equivalently,
  V divides beta and V^3 divides 4 gamma V^2-beta^2. This is an exact
  two-jet criterion, not a proof that Keller equations force it, so the
  nonmonic quartic stratum and planar JC remain open.
source: codex-2026-07-22-JC2-quartic-deck-parity
depends_on: []
related:
  - THM-2129
  - THM-2136
  - THM-2141
  - THM-2147
---

# THM-2158 -- the quartic deck is not the pole obstruction

Let

```text
R=C[x],                    K=C(x),
```

and choose `V,beta,gamma,delta,epsilon in R` with `V!=0`. Consider

```text
P=V^2 z^4+beta z^3+gamma z^2+delta z+epsilon.        (1)
```

THM-2129 leaves (1) in the twice-odd quartic branch: the leading
coefficient is forced to be a square, but the usual monic coordinate is
introduced only after a quadratic base change. We prove that this quadratic
deck loses no information relevant to square completion. It is only a device
for displaying a canonical object which is already rational over the base.

## 1. Monic depression and deck parity

Work in the quadratic `K`-algebra

```text
L=K[U]/(U^2-V).                                      (2)
```

This formulation includes the split case. Since `V` is nonzero in `K`, the
algebra is etale, `U` is a unit, and it has the involution

```text
sigma(U)=-U,                  sigma|K=id.             (3)
```

Put

```text
w=Uz,
B_3=beta/U^3,        B_2=gamma/U^2,
B_1=delta/U,         B_0=epsilon.                    (4)
```

Then

```text
P=w^4+B_3 w^3+B_2 w^2+B_1 w+B_0.
```

Make the unique translation which kills the cubic term,

```text
Z=w+B_3/4=Uz+beta/(4U^3).                            (5)
```

Direct expansion gives

```text
P=Z^4+p Z^2+q Z+r,                                  (6)

p=gamma/U^2-3 beta^2/(8U^6),

q=delta/U-beta gamma/(2U^5)+beta^3/(8U^9),

r=epsilon-beta delta/(4U^4)
          +beta^2 gamma/(16U^8)-3 beta^4/(256U^12). (7)
```

Every parity is now visible:

```text
sigma(Z)=-Z,       sigma(p)=p,
sigma(q)=-q,       sigma(r)=r.                       (8)
```

Therefore both terms in the canonical square decomposition

```text
H_0=Z^2+p/2,
P=H_0^2+[qZ+(r-p^2/4)]                               (9)
```

are fixed by the deck involution. This is more than formal invariance: they
can be written explicitly over `K`.

Indeed, put

```text
zeta=z+beta/(4V^2),
Q_1=delta-beta gamma/(2V^2)+beta^3/(8V^4).           (10)
```

Then `Z=U zeta`, `q=Q_1/U`, and hence

```text
H_0
 =V zeta^2+(1/2)(gamma/V-3 beta^2/(8V^3))
 =V z^2+[beta/(2V)]z
      +gamma/(2V)-beta^2/(8V^3),                    (11)

qZ=Q_1 zeta in K[z].                                 (12)
```

Also `r-p^2/4` lies in `K`. Thus (9) is already a decomposition in
`K[z]`; no invariant-field theorem and no choice of a square root of `V`
is needed to descend it.

## 2. The square approximate root is canonical

The rational polynomial (11) is characterized without mentioning the deck.
Suppose

```text
H=V z^2+a z+c in K[z]                                (13)
```

and `deg_z(P-H^2)<=1`. Comparing the coefficients of `z^3` and `z^2`
forces

```text
2Va=beta,
a^2+2Vc=gamma.                                      (14)
```

Consequently

```text
a=beta/(2V),
c=gamma/(2V)-beta^2/(8V^3),                         (15)
```

so `H=H_0`. Conversely (14) verifies directly that `P-H_0^2` has
`z`-degree at most one. Hence `H_0` is the unique degree-two approximate
root with leading coefficient `V`.

This uniqueness is the exact descent statement. On the quadratic deck, the
monic depressed square root must be `Z^2+p/2`; over the base, it must be the
same uniquely characterized polynomial (11).

## 3. Exact finite-pole criterion

Because the coefficients of a polynomial in `R[z]` are independent, (11)
gives

```text
H_0 in R[z]
 iff beta/V in R
     and (4 gamma V^2-beta^2)/V^3 in R.              (16)
```

Since `R=C[x]` is a UFD, this is equivalently

```text
V divides beta,
V^3 divides 4 gamma V^2-beta^2.                     (17)
```

Writing `beta=V beta_1`, the second condition becomes

```text
V divides 4 gamma-beta_1^2.                         (18)
```

Thus the obstruction is exactly the first two coefficient jets of the
quartic at the finite divisor `V=0`. If (17) holds, both `H_0` and
`P-H_0^2` belong to `R[z]`, and the latter has `z`-degree at most one. If
(17) fails, deck invariance still holds but the canonical approximate root
has a genuine finite pole.

The two conditions are independent and sharp. With `V=x`:

```text
P_1=x^2 z^4+z^3
```

has

```text
H_0=xz^2+z/(2x)-1/(8x^3),                            (19)
```

so the first divisibility already fails. The polynomial

```text
P_2=x^2 z^4+xz^3
```

passes `V|beta` but has

```text
H_0=xz^2+z/2-1/(8x),                                 (20)
```

so the second jet still carries a pole. In the positive control

```text
P_3=x^2 z^4+xz^3+(1/4)z^2+z
   =(xz^2+z/2)^2+z,                                  (21)
```

both conditions hold and the canonical approximate root is polynomial.
These are algebraic controls, not asserted Keller components.

## 4. Consequence and exact scope

The remaining nonmonic-quartic route can now be factored as

```text
twice-odd leading equation
  -> leading coefficient V^2
  -> canonical H_0 in C(x)[z], with automatic deck descent
  -> prove the two finite-pole divisibilities (17)
  -> polynomial exact-square prefix
  -> apply THM-2147 under its weighted hypotheses.                  (22)
```

Accordingly, “quadratic-deck compatibility” is not a separate mathematical
obstruction. The actual missing lemma is that the Keller coefficient
equations force (17), or else provide a descent which strictly reduces the
finite pole divisor. This theorem does not prove either assertion. In
particular, deck invariance must not be confused with polynomial regularity,
as (19)--(20) demonstrate. The general quartic source-fiber stratum, JC(2),
and DC(2) remain open.

There is no faithful Tournament Analysis in this argument. The carrier is
the ordered coefficient jet `(V,beta,gamma)`, the quadratic involution, and
the pole divisor. Replacing these by pairwise orientations would erase the
valuation multiplicities in (17), which are exactly the surviving data.

QED.
