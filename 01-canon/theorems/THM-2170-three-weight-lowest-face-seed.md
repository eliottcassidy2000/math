---
id: THM-2170
title: "Explicit three-weight lowest-face seed for arbitrary radial coefficients"
status: >
  PROVED as an explicit specialization of THM-2022. For
  P=Z^p a(ZW)+b(ZW)+W^q c(ZW), the lowest balanced Wick face is determined
  by the radial orders of a,b,c. Its nonzero constant-term seed occurs at
  return level 1 or (p+q)/gcd(p,q), so the finite-place
  Kummer--Lucas--Frobenius amplifier closes every arbitrary-radial
  three-weight slice without a separate THM-2111 seed invocation and without
  leading-factorial dominance. This repairs the mechanism of THM-1515 and
  makes THM-2017's
  resonance bands effective refinements rather than NC2 obligations.
source: codex-2026-07-24-NC2-radial-face-synthesis
depends_on:
  - THM-2022
related:
  - THM-1515
  - THM-2014
  - THM-2017
  - THM-2018
  - HYP-8765
  - HYP-8766
---

# THM-2170 -- explicit three-weight lowest-face seed

Put `s=ZW` and let

```text
P=Z^p a(s)+b(s)+W^q c(s),             p,q>=1,          (1)
```

where `a,b,c` are nonzero complex polynomials. Write

```text
g=gcd(p,q),       p0=p/g,       q0=q/g,
r=p0+q0,                                              (2)
A=ord_s a,        B=ord_s b,     C=ord_s c,            (3)
alpha=[s^A]a,     beta=[s^B]b,   gamma=[s^C]c.         (4)
```

Thus `alpha beta gamma!=0`. The theorem computes the exact face seed needed
by THM-2022; it does not use radial degrees or an asymptotic channel order.

## 1. Balanced multiplicities and Wick height

Choose `x_+` terms from `Z^p a`, `x_0` terms from `b`, and `x_-` terms from
`W^q c`. Charge balance says

```text
p x_+=q x_-.                                          (5)
```

Hence every balanced charged block has

```text
x_+=q0 k,       x_-=p0 k                              (6)
```

for an integer `k>=0`, and has length `rk`. At the lowest radial monomial in
each polynomial, its common Wick exponent is

```text
D_ch=q0(p+A)+p0 C.                                    (7)
```

Each neutral letter contributes length one and Wick exponent `B`. Higher
radial monomials of the same charge strictly increase the Wick exponent.
Therefore every balanced channel with `k` charged blocks and `l` neutral
letters has

```text
length=rk+l,
Wick height>=k D_ch+l B,                              (8)
```

with equality exactly on the displayed lowest radial monomials.

For fixed length, the lowest face is determined by comparing the two
per-letter heights

```text
D_ch                 and                 rB.           (9)
```

This is the complete linear program because (5) leaves only the neutral ray
and the primitive charged ray.

## 2. The three explicit face seeds

Let `f_F(u)` be THM-2022's Laurent polynomial on the lowest face.

### Neutral side: `rB<D_ch`

The neutral ray has smaller height per letter. Thus

```text
f_F(u)=beta,
m0=1,
Q=CT_u(f_F)=beta!=0,
A0=B.                                                  (10)
```

### Charged side: `rB>D_ch`

The primitive charged ray has smaller height per letter. Thus

```text
f_F(u)=alpha u^p+gamma u^(-q),
m0=r,
Q=CT_u(f_F^r)
 =[r!/(q0!p0!)] alpha^q0 gamma^p0 !=0,
A0=D_ch.                                               (11)
```

Indeed the only balanced term in the `r`-th power uses `q0` positive-charge
letters and `p0` negative-charge letters.

### Tied side: `rB=D_ch`

All three lowest monomials lie on one face:

```text
f_F(u)=alpha u^p+beta+gamma u^(-q).                    (12)
```

There is no need to separate its colliding return channels. At level one,

```text
m0=1,
Q=CT_u(f_F)=beta!=0,
A0=B.                                                  (13)
```

Thus in every case a nonzero face seed exists with

```text
m0 in {1,r}.                                           (14)
```

No Duistermaat--van der Kallen or compound-root existence theorem is needed
to produce the seed for this three-charge support.

## 3. Whole-face amplification

Apply THM-2022 after the face has been exposed. Algebraically descend the
coefficient point while inverting every active coefficient and the explicit
`Q`. Choose a good residue characteristic `ell` only afterward. At moment
order

```text
M=ell*m0,                                              (15)
```

normalize by `(ell*A0)!`.

The three-case residue calculation of THM-2022 applies verbatim:

1. Kummer kills every multiplicity vector not divisible by `ell`;
2. strict face height kills every dilated off-face vector; and
3. Lucas and Frobenius preserve the complete tied face.

The normalized moment reduces in the chosen residue field to

```text
overline(E[P^(ell*m0)]/(ell*A0)!)=Qbar^ell!=0.         (16)
```

Consequently (1) cannot have all Gaussian moments zero.

The detecting characteristic `ell`, and hence the moment order (15), may
depend on the coefficients. Equation (14) is an explicit seed bound, not a
coefficient-uniform moment cutoff.

## 4. Degenerate omissions and nullcone statement

The omitted-polynomial cases require no artificial convention for
`ord_s(0)`.

- If `b!=0` and either charged side is absent, every balanced word is purely
  neutral; its lowest coefficient gives `m0=1,Q=beta`.
- If `b=0` but both charged sides are present, (11) is the only balanced
  face and gives `m0=r`.
- If `b=0` and at most one charged side is present, the polynomial is strict
  one-sided and all positive moments vanish.

Therefore, for arbitrary radial `a,b,c`, including zero polynomials,

```text
E[P^m]=0 for every m>=1
iff
b=0 and (a=0 or c=0).                                  (17)
```

This is the exact three-weight nullcone, for arbitrary charges `p,q`.

## 5. Repair and scope

For `p=q=1`, equation (17) is the conclusion claimed by THM-1515. Its
published leading-factorial proof is invalid by MISTAKE-202: lower channels
can have exponentially many allocations and cannot be discarded by
termwise factorial order. Equations (5)--(16) give the valid replacement:

```text
lowest s-adic face
 -> explicit nonzero whole-face seed
 -> prime chosen after the seed
 -> Kummer/Lucas/Frobenius preservation.               (18)
```

THM-2017's degree-gap and hyper-Bessel results remain sharper asymptotic
information, but no resonance offset is left open for NC2. This theorem is
a transparent specialization of the already general THM-2022, not a second
proof of full NC2 and not a coefficient-uniform detection theorem.
