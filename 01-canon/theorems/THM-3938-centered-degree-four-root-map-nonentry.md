---
id: THM-3938
title: "Centered degree-four repeated-root maps are scalar or fold a ramification arm"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the centered
  trace-zero linear-color binary-cubic
  grammar with deg A=3 and repeated-root degree four, local trace and the
  finite Riemann--Hurwitz budget leave five collision-free non-root-regular
  pole rows. Exact color division leaves one family in each row. Two
  exceptional parameter values give literal scalar endpoints and hence
  monogenic cubic orders. Every genuinely non-scalar survivor has an
  irreducible one-place implicit component occurring to exponent one in the
  order discriminant and a forced two- or three-address target fibre.
  Rank-three finite flatness coalesces those addresses at one source
  ramification point, so the source arm is non-unibranch and THM-3920
  forbids an affine-plane Keller open. This closes the stated centered
  degree-four grammar, not arbitrary root gauges, root degree at least five,
  other coefficient planes, or JC(2).
source: jc_zero_debt_lift / post-THM-3936 next root-degree stratum, 2026-08-23
audit: >
  FIVE INDEPENDENT HOSTILE AUDITS PASS (degree4_independent, root,
  jc3913_lattice_referee, jc_decic_lattice, and jc_tournament_response,
  2026-08-24). The first audit
  independently reconstructed the pole/Riemann--Hurwitz exhaustion, including
  the finite-infinity gauge and reversed row-F orientation; all color seams,
  roots at infinity, and pole cancellations; the scalar normalizations and
  coefficient ideals; every address factor and critical-value exclusion; the
  exponent-one maximal-order parity step; and the rank-three source
  non-unibranch bridge. The assertion-free 96-gate companion byte-matches in
  normal and optimized mode; frozen output, hashes, and documentation checks
  pass. The further hostile audits independently
  recovered the five color survivors, the exact two-address rather than
  three-address count in rows G--I, and the two scalar seams. They added the
  codimension-one maximality plus R1+S2 normality check in (25a)--(25b), so
  THM-3801 is applied to the actual maximal order rather than an auxiliary
  binary-cubic order.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3936-centered-degree-three-infinite-root-value-nonentry
related:
  - THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry
  - THM-3930-two-pole-linear-color-aligned-one-place-branch-packet
  - THM-3931-degree-two-pole-cubic-principal-ramification-no-atlas
  - THM-3941-all-degree-centered-cubic-pole-carrier-routing
script: 04-computation/jc2_centered_degree_four_root_map_thm3938.py
output: 05-knowledge/results/jc2_centered_degree_four_root_map_thm3938.out
script_sha256: d938b8af519e04d4c8d7957dce492fe28feec6ab63eb2f433d20bf2844bf65a3
output_sha256: deefbdc002977b2992d0138ae443856a7d0f16e4497fa4e191661eafbb0cc1dc
semantic_sha256: 8df0539a4608fb0a83929d851cae25cefdd0150bf0f3f6d828ba77e8a982ceac
hash_basis: raw LF bytes
---

# THM-3938 -- degree four has five doors, and every door returns to the same folded arm

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of
characteristic zero. Consider

```text
Phi=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3,                      (1)
```

and an irreducible one-place discriminant component with normalization
`A1_u`. Assume

```text
A=A(u), C=C(u) in k[u],       deg A=3,                    (2)
deg(t:P1_u -> P1_[U:V])=4,    t=U/V,                      (3)
```

and use the centered primitive incidence

```text
a(A)t^3-c(A)t-2d(A)=0.                                    (4)
```

Thus `Tr_(k(u)/k(A))(t)=0`. We also impose

```text
(a,C,c,d)=k[A,C],                                         (5)
```

and the normalization condition `k(A,C)=k(u)`. Separately, trace zero and
nonconstancy force `t notin k(A)`: an element of the base field would have
trace `3t`. Since `[k(u):k(A)]=3` is prime, `k(A,t)=k(u)`, and `(4)` is the
minimal polynomial after removal of coefficient content.

The theorem proves that no such component belongs to the finite completion
of a Keller `A2` open. The fixed linear-color and centered-root hypotheses
are essential scope; a root translation generally mixes the `C` slot and
does not preserve `(1)`.

## 1. Local trace leaves exactly five nonregular rows

The cubic polynomial map

```text
pi_A:P1_u -> P1_A
```

is totally ramified at infinity, so only two finite Riemann--Hurwitz units
remain. At a tame point with local equation `z=s^e`, trace keeps precisely
the Laurent exponents divisible by `e`. Hence an isolated finite pole obeys

```text
order 1: e=2 or 3;
order 2: e=3;
order 3: e=2;
order 4: e=3.                                             (6)
```

At polynomial infinity `e=3`; an order-three pole is impossible because
its leading exponent survives trace. Order four is allowed, while orders
one and two are allowed. If two distinct finite pole supports share an
`A`-value, both repeated roots limit to `[1:0]`. In the `U=1` chart

```text
Phi=a+Cw+cw^2+dw^3,
```

so both addresses have `a=C=0`. The exact-double/transversality and
rank-three fibre argument of THM-3936 then puts the two normalization
addresses on one source ramification point, making that source component
non-unibranch. We may discard every shared-address row before normalization.

For the remaining collision-free rows, the finite pole partition and its
minimum finite ramification cost are

```text
t(infinity) finite:
  4        costs 2;       3+1      costs 2;
  2+2, 2+1+1, 1+1+1+1 each cost at least 4.

m_infinity=1:
  3        costs 1;       2+1 and 1+1+1 cost at least 3.

m_infinity=2:
  2        costs 2;       1+1      costs 2.

m_infinity=3: forbidden by trace.
m_infinity=4: t is polynomial and THM-3929 applies.        (7)
```

Thus exactly five non-root-regular rows remain. Affine changes of `u,A`, a
diagonal root scaling, and the compensating nonzero color scaling give the
complete trace-normalized table

```text
row  pole grammar                 A(u)       t(u)

E    finite value; 4@e3          u^3        p/u+q/u^2+1/u^4

F    finite value; 3@e2+1@e2     u^3-3u     p/(u-1)-1/[2(u-1)^2]
                                             +1/(u-1)^3+q/(u+1)

G    m_inf=2; 2@e3               u^3        u^2+pu+q/u+k/u^2

H    m_inf=2; 1@e2+1@e2          u^3-3u     u^2+pu-2+q/(u-1)+k/(u+1)

I    m_inf=1; 3@e2               u^3+u^2    u+1/3+p/u-3k/(2u^2)+k/u^3.
                                                                  (8)
```

The nonzero exact-pole conditions are: the normalized top principal
coefficient in E/F is already one; `q!=0` in F; `k!=0` in G and I; and
`qk!=0` in H. The constants in `(8)` are forced by trace:

* in E, `Tr(u^-3)=3/A` kills the missing `u^-3` term, while trace of a
  finite value is three times that value; hence `t(infinity)=0`. The same
  constant-term argument forces the finite value in F to zero;
* in F, trace kills the constant and forces twice the order-two coefficient
  plus the order-three coefficient to vanish;
* in G, every displayed nonconstant term has trace zero;
* in H, `Tr(u^2)=6`, forcing the constant `-2`;
* in I, `Tr(u)=-1`, `Tr(u^-2)=2/A`, and `Tr(u^-3)=3/A`, forcing `1/3`
  and the ratio `-3k/2`.

This proves completeness of the five rational normal forms, including their
lower-order jets. In row F the apparently reversed orientation (triple pole
at `u=-1`, simple pole at `u=1`) is the same affine orbit under `u -> -u`,
followed by the allowed diagonal root/color rescaling; it is not a sixth row.

## 2. Exact color division

Eliminate `u` from `A-A(u)` and the denominator-cleared `T-t(u)`. The
primitive trace-zero resultant gives `a,c,d` in `(4)`, and the repeated-root
derivative gives

```text
C(u)=-(3a(A(u))t(u)^2+c(A(u)))/(2t(u)).                   (9)
```

The exact numerator/denominator resultants of `(9)` are

```text
E:  27 p^2 q^6,
F: -2^35 3^3 q^3(p+q)^2(2p+6q+1)^6,
G:  81 k^6(pq-k)^6,
H: -5184 k^3q^3[kp+k+2p^2+pq-q-2]^6,
I:  2^23 3^20 k^9(3p+2)^6(81k+18p+4)^3.                 (10)
```

Resultant zero is only a necessary first cut, so every seam must be divided
exactly.

### 2.1 Row E

On `p=0`, polynomial division leaves the nonzero remainder `3/q^4` unless
`q=0`. The seam `q=0` divides exactly. Hence

```text
t=(pu^3+1)/u^4,
C=-3u^8(pA+1)/2,
a=A^4, c=0, d=(pA+1)^3/2.                                (11)
```

At `p=0`, `d=1/2` is a literal scalar value of the binary index form. The
index criterion in THM-3801 therefore makes the cubic order monogenic, while
its Keller gate requires nonmonogenicity. From now on row E has `p!=0`.

### 2.2 Row F

On `p=-q`, the exact pseudo-remainder is

```text
-49152(4q+1)^7[(4q^2+5q+1)u-4q^2+q+1].                  (12)
```

Away from `q=-1/4`, its two coefficients cannot vanish together. At
`q=-1/4`, the seam overlaps the other factor. The second seam

```text
p=-3q-1/2                                                     (13)
```

divides for every `q`. After removing a harmless common scalar eight, put

```text
L_F=(4q+1)A+8q-2.
```

The survivor is

```text
a=(A-2)(A+2)^3,
c=3(2A-5)L_F^2/4,
d=-L_F^3/16,
C=3(u-1)^3(u+1)(u^4+6u^3-22u-21)L_F(A(u))/4.             (14)
```

The pole requires `q!=0`. At `q=-1/4`, `L_F=-4` and `d=4` is a literal
scalar index-form endpoint, hence monogenic and excluded by THM-3801's Keller
gate. The non-scalar row therefore has `q!=0,-1/4`.

### 2.3 Row G

The first seam in `(10)` is `k=pq`. Substitution leaves

```text
-3p^6(p^3-q)                                               (15)
```

as the final division remainder. Since `k!=0`,

```text
q=p^3, k=p^4, p!=0.                                       (16)
```

Put `L_G=A+p^3` and `Q_G=u^2-pu+p^2`. Then

```text
t=(u+p)^2Q_G/u^2,
C=-3u^4Q_G(u^2+3pu+p^2)/2,
a=A^2, c=3ApL_G^2, d=L_G^4/2.                             (17)
```

### 2.4 Row H

With `qk!=0`, the first seam is

```text
k=-(p-1)(2p+q+2)/(p+1).                                  (18)
```

The case `p=-1` would require `q=0` and is not live. Exact division after
`(18)` leaves

```text
-3(p-1)^3(p+1)^3[p^4+p^3-3p^2-5p-2q-2].                 (19)
```

The `p=1` factor makes `k=0`; hence the genuine survivor is

```text
q=(p-2)(p+1)^3/2,
k=-(p-1)^3(p+2)/2,
p notin {-2,-1,1,2}.                                      (20)
```

Put `L_H=A+p^3-3p` and `Q_H=u^2-pu+p^2-3`. Then

```text
t=(u+p)^2Q_H/[(u-1)(u+1)],
C=-3(u-1)(u+1)Q_H
  [p^2u^2-5p^2+3pu^3-11pu+u^4-4u^2-1]/2,
a=A^2-4,
c=-3L_H^2(-Ap+p^2+1),
d=L_H^4/2.                                                 (21)
```

### 2.5 Row I

On the first seam `p=-2/3`, division leaves

```text
-1024(81k-40)/729,
```

so `k=40/81`. On the other seam `k=-(18p+4)/81`, the exact remainder is

```text
4(9p+2)^4[81p^2-162pu-72p-12u-20]/243.                   (22)
```

The prefactor vanishes only at `p=-2/9`, which gives `k=0` and deletes the
triple pole. Killing the coefficient of `u` instead gives `p=-2/27`, where
the constant in brackets is nonzero. Thus the unique live row is

```text
p=-2/3, k=40/81.                                          (23)
```

Let

```text
L_I=27A-20,       Q_I=9u^2+15u+10.
```

Then

```text
t=(3u-2)^2Q_I/(81u^3),
C=-4u^3Q_I(27u^5+45u^4-30u^3-60u^2+16)/3,
a=216A^3,
c=-8(15A-4)L_I^2/243,
d=4L_I^4/19683.                                           (24)
```

Equations `(11)`, `(14)`, `(17)`, `(21)`, and `(24)` are the full
polynomial-color list. No census or bounded coefficient assumption enters.

## 3. Coefficient ideals and the two scalar seams

Every genuine survivor passes `(5)`:

* E: `gcd(A,pA+1)=1`;
* F: `L_F(2)=16q` and `L_F(-2)=-4`;
* G: `L_G(0)=p^3`;
* H:

  ```text
  L_H(2)=(p-1)^2(p+2),  L_H(-2)=(p+1)^2(p-2),
  ```

  nonzero precisely away from the four pole cancellations;
* I: `L_I(0)=-20`.

Thus coefficient primitivity does not kill the non-scalar rows. Conversely,
the two exceptional parameters

```text
E:p=0 gives d=1/2,       F:q=-1/4 gives d=4              (25)
```

are not merely failed division cases: they are actual scalar values of the
binary index form and are excluded by THM-3801's nonmonogenicity gate.

There is a necessary maximal-order check before invoking that gate. At the
row-E endpoint `p=0`, exact specialization gives

```text
J_E=(27A^8+8C^3)/8,                 Disc(O_E)=-2J_E.       (25a)
```

At the row-F endpoint `q=-1/4`, one has `L_F=-4` and

```text
Disc(O_F)=-16J_F.                                        (25b)
```

The pole-valuation generation argument of Section 4 remains valid on both
seams, so each `J` is irreducible. Its exponent is one; hence the displayed
finite-free cubic order is maximal at the sole height-one discriminant
prime and etale, therefore normal, at every other height-one prime. Finite
freeness gives `S2`, so `R1+S2` makes the whole order normal. Consequently
the scalar values in `(25)` belong to the actual maximal normalization, not
merely an auxiliary nonmaximal order. THM-3801 now applies exactly.

## 4. Every non-scalar implicit component is genuine ramification

For each surviving polynomial parametrization, define its implicit
resultant

```text
J_i(A,C)=Res_u(A-A_i(u),C-C_i(u)).                         (26)
```

Exact discriminant calculation gives

```text
Disc(Phi_E)=-2(pA+1)^3 J_E,
Disc(Phi_F)= L_F^3 J_F/4,
Disc(Phi_G)=-2L_G^4 J_G,
Disc(Phi_H)=-2L_H^4 J_H,
Disc(Phi_I)=-16L_I^4 J_I/19683.                           (27)
```

The `C`-degree of every `J_i` is three. The total degrees in the non-scalar
range are respectively

```text
11, 11, 8, 8, 10.                                        (28)
```

Their normalizations have one polynomial infinity. More importantly, every
`J_i` occurs to exponent **one** in the order discriminant. An
order-to-maximal-order correction is a square, so the maximal discriminant
still has valuation one at `J_i`. Hence `J_i` supports a genuine tame
residue-degree-one `(2,1)` ramification prime `E_i`.

The implicit equations are irreducible and the parametrizations are
birational. Indeed, over `k(A)` the cubic `(1)` is linear in `C`; a
factorization over `k(A)[C,T]` would have a `C`-independent factor dividing
both `T^2` and `aT^3+cT+d`, impossible because `d!=0`. Also `t notin k(A)`:
at a finite pole its order modulo the local ramification index is

```text
E:4 mod 3; F:3 mod 2; G:2 mod 3; H:1 mod 2; I:3 mod 2.   (29)
```

Thus `k(A,t)=k(u)`. At the generic exact double root, Euclid's algorithm
recovers `t` from the cubic and its derivative, so `k(A,C)=k(u)`. Therefore
`J_i` is the irreducible relation, not a parametrization power.

## 5. The address packets are on the source ramification arm

Each non-scalar row has a target point with multiple distinct normalization
addresses.

### E: three addresses

For `p!=0`, the three roots of

```text
u^3=-1/p
```

all map to `(A,C)=(-1/p,0)`. They are distinct because their discriminant is
`-27/p^2`.

### F: three addresses

For `q!=0,-1/4`, put

```text
A_0=(2-8q)/(4q+1).
```

The factor `L_F(A(u))` in `(14)` sends all three roots of `A(u)=A_0` to
`(A_0,0)`. The identities

```text
A_0-2=-16q/(4q+1),       A_0+2=4/(4q+1)
```

show that `A_0` is not a critical value, so the three addresses are distinct.

### G: two addresses

The identity

```text
(u+p)Q_G=u^3+p^3
```

shows that the two roots of `Q_G` map to `(-p^3,0)`. Their discriminant is
`-3p^2`, nonzero because `p!=0`.

### H: two addresses

Likewise

```text
(u+p)Q_H=u^3-3u+p^3-3p.
```

The two roots of `Q_H` map to `(-(p^3-3p),0)`. Their discriminant is
`-3(p-2)(p+2)`, and

```text
Q_H(1)=(p-2)(p+1),       Q_H(-1)=(p-1)(p+2),
```

so in the genuine range `(20)` they are distinct and avoid both pole
supports.

### I: two addresses

Finally,

```text
(3u-2)Q_I=27(u^3+u^2)-20.
```

The two roots of `Q_I` map to `(20/27,0)` and are distinct because
`Disc(Q_I)=-135`.

Now let `B` be the maximal normalization of the target plane in the cubic
field. It is a normal finite module over the regular two-dimensional target,
hence Cohen--Macaulay and finite flat of rank three. The normalization of
the residue-degree-one source ramification prime `E_i` is `J_i^nu=A1_u`.
Every point of `E_i` is non-etale and consumes local geometric fibre length
at least two. Two different `E_i` points cannot lie in the same rank-three
fibre. Consequently every address packet above coalesces at a single point
of `E_i`.

Thus `E_i` itself has two or three branches at one point. This is a source
ramification obstruction, not merely a collision of the branch image. A
Keller `A2` open must delete `E_i`, while THM-3920 says every irreducible
boundary curve of an affine-plane open in a normal affine surface is
unibranch. No non-scalar survivor admits a Keller plane atlas.

Combining this with the two scalar seams `(25)` and THM-3929's polynomial
root row proves the theorem.

## 6. Scope and next frontier

THM-3933 and THM-3936 close centered root degree three; this theorem closes
centered root degree four. The mechanism is now recursive: local trace
turns pole partitions into a finite grammar, color division leaves a few
exact families, and source-address multiplicity tests the survivors.

Still open are root degree at least five, noncentered gauges that do not
preserve the linear-color slot, nonlinear color directions, other binary
cubic coefficient planes, and JC(2).

## Reproduction

```bash
python3 04-computation/jc2_centered_degree_four_root_map_thm3938.py
python3 -O 04-computation/jc2_centered_degree_four_root_map_thm3938.py
```

The frozen transcript is
`05-knowledge/results/jc2_centered_degree_four_root_map_thm3938.out`.
