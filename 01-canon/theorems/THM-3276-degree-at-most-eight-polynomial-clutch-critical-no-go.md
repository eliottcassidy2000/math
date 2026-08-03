---
id: THM-3276
title: "Degree-at-most-eight polynomial clutch critical no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For
  either THM-3212 cubic accessory response pair and every polynomial clutch
  B of degree at most eight, the displayed first coordinate P_B has a
  critical point. If B is a unit on the owner divisor, at least 43 units of
  saturated critical-resultant multiplicity remain away from that divisor.
  Hence no member is a coordinate of a polynomial Keller pair. This excludes
  one explicit deformation family and proves no instance of JC(2).
source: root/2026-08-03
audit: >
  The assertion-independent exact companion pins THM-3225, THM-3237 and
  THM-3257, rebuilds both cubic accessory fields and response pairs, verifies
  the universal degree and boundary ledgers, and uses a fresh exact Laurent
  engine to eliminate the eight available simple-boundary jets successively.
  In each field it proves that the first two untunable Laurent coefficients
  have coprime numerators. An independent audit rederived the elimination,
  degree and boundary ledgers, proved the exact triangular slope
  16*j*c^4*v_1^3, and reconstructed the two coprime untunable pairs through
  separate rational-function and good-reduction engines. Normal, optimized
  and stored outputs agree byte-for-byte.
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
  - THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape
related:
  - THM-3237-degree-nine-jacobian-infinity-wall-and-square-root-escape
  - THM-3257-degree-eight-tuned-cubic-infinity-wall-and-three-root-critical-escape
  - THM-3265-degree-six-retuned-quintic-infinity-wall-and-five-root-critical-escape
script: 04-computation/jc_degree_at_most_eight_uniform_clutch_no_go_wildcard_20260803.py
output: 05-knowledge/results/jc_degree_at_most_eight_uniform_clutch_no_go_wildcard_20260803.out
script_sha256: acdd3b274f3e2692040219626c93302bd14e3bc729790c1a2d7aab5798f35790
output_sha256: 0472f414b5d55f04aaa0c165734daff6a3144fb456bd1ee651283b041e511c29
hash_basis: LF-normalized bytes
---

# THM-3276 -- degree-at-most-eight polynomial clutch critical no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Universe and statement

Let `K_i` be either cubic accessory field of
[THM-3212](THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch.md),
and embed it in an algebraically closed characteristic-zero field `K_0`.
Retain its response pair and owner divisor

```text
V=4SDT^2/Gamma^2,       A=2SET/Gamma,
g=S*T,                  2VA'-AV'=2V,
deg V=16,               deg A=8.                         (1)
```

For an arbitrary polynomial `B in K_0[x]` with `deg B<=8`, define

```text
P_B(x,z)=(V(x)z^2+B(x)z)^2+A(x)z+x.                     (2)
```

Then `P_B` has a critical point. More precisely, either:

1. `gcd(B,g)!=1`, and `(2)` has an explicit critical point over a common
   root; or
2. `gcd(B,g)=1`, and at least 43 units of the saturated critical-resultant
   multiplicity lie at roots `x_0` with `V(x_0)!=0`.

At a critical point of `P_B`, `Jac(P_B,Q)=0` for every polynomial `Q`.
Consequently no polynomial `B` in this degree range gives a constant-Jacobian
mate for `(2)`.

## 2. Gradient reduction and universal resultant

On `V!=0`, put

```text
y=Vz,              L=y^2+By,
J_B=2VB'-BV'.                                               (3)
```

The response identity in `(1)` reduces the two gradient equations exactly to

```text
R_1=2L(2y+B)+VA,
R_2=V^3+V^2y+J_B yL.                                      (4)
```

As proved in
[THM-3225](THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape.md),
these generate the localized gradient ideal and

```text
Res_y(R_1,R_2)=V^3 K_B,                                   (5)
```

where

```text
K_B=
 -A^3J^3+12A^2J^2V^2-4AB^3J^2V+4AB^2JV^2
 +24ABJV^3-48AJV^4-16AV^4-8B^4JV^2+8B^3JV^3
 +32B^2V^4-96BV^5+64V^6,                                 (6)
```

and `J=J_B`.

## 3. Uniform degree and owner-boundary ledger

If `d=deg B<8`, then `deg J_B<=d+15`. At `d=8`, the nominal top coefficient
of `J_B` is proportional to `2d-16` and cancels, so `deg J_B<=22`. Substitution
in all twelve terms of `(6)` gives

```text
deg K_B=96        for every d=0,...,8,                    (7)
```

and `64V^6` is the unique degree-96 term. Thus its leading coefficient cannot
cancel for any choice of `B`.

The monic passport boundaries are

```text
boundary_(4111)=S^3T^8x^9,
boundary_(3211)=S^3T^8x^6(x-1)^3,                        (8)
```

both of degree 44. The local response identity proves that `(8)` divides
`K_B` for arbitrary `B`; write

```text
K_B=boundary_i H_B.                                      (9)
```

Equations `(7)--(9)` give

```text
deg H_B=52.                                               (10)
```

At a root `alpha` of `T`, write

```text
V=v t^m+...,       A=2t/(2-m)+...,
m in {3,4,5,6},    b=B(alpha).                            (11)
```

If `b!=0`, direct substitution in `(6)` gives the exact first coefficient

```text
[t^(3m-1)]K_B=
 16m(m-1)/(m-2) b^5 v^3 !=0.                            (12)
```

The exponents `3m-1` are exactly those in `(8)`. Hence `H_B` has no `T` root
whenever `B` is a unit on `g`.

If instead `B(alpha)=0` at any root of `g`, then `V(alpha)=A(alpha)=0` and
`A'(alpha)!=0`. On that fibre

```text
P_z(alpha,z)=0,             P_x(alpha,z)=A'(alpha)z+1,
```

so

```text
(alpha,-1/A'(alpha))                                      (13)
```

is an explicit critical point. This proves the first lane of the theorem.

## 4. The complete degree-eight jet invoice at `S`

Assume now `gcd(B,g)=1`. At the simple root `s` of `S`, use `t=x-s` and
write

```text
V=v_1t+v_2t^2+...,
A=2t+a_2t^2+...,
B=c+b_1t+...+b_8t^8,              c!=0.                  (14)
```

The coefficient of `t^2` in `(6)` cancels identically, giving the inherited
`S^3` factor. After the earlier jets have been fixed, the coefficient of
`t^(j+2)` is affine in the new jet `b_j`, for `j=1,...,8`. In both exact
accessory fields its exact slope is

```text
16j c^4 v_1^3 !=0.                                      (15)
```

Therefore there is a unique successive tuning which kills the eight
coefficients

```text
[t^3]K_B,...,[t^10]K_B.                                  (16)
```

Working exactly in the Laurent field `K_i[c,c^(-1)]`, the tuned jet degree
profiles `(numerator degree, denominator degree)` are

```text
(1,0),(3,2),(4,3),(6,5),(7,6),(9,8),(10,9),(12,11).      (17)
```

No degree-eight jet remains after `(16)`. Put

```text
F_11(c)=[t^11]K_B,              F_12(c)=[t^12]K_B.       (18)
```

The exact two-field calculation gives

```text
                   numerator degree   denominator
F_11                       13              c^8
F_12                       15              c^10,

gcd(num(F_11),num(F_12))=1                               (19)
```

in each `K_i[c]`. Since `c!=0`, `(19)` says the two untunable coefficients
cannot vanish together over any algebraic closure of `K_i`. Thus every
degree-at-most-eight clutch satisfies

```text
ord_S(K_B)<=12,                  ord_S(H_B)<=9.           (20)
```

This is the load-bearing finite-dimensional obstruction: eight adjustable
jets can buy at most nine units of saturated contact at the simple owner
root.

## 5. Surviving critical points

In the coprime lane, `(12)` excludes every `T` root and `(20)` absorbs at
most nine of the 52 units in `(10)` at `S`. Hence at least

```text
52-9=43                                                     (21)
```

units remain at roots of `H_B` away from `g`, where `V!=0`.

The leading `y` coefficient of `R_1` is the constant four. Therefore each
such resultant root supports an affine common zero of `(4)`, and the exact
localized gradient reduction turns it into a critical point of `P_B`.
Multiplicity, not distinctness, is counted in `(21)`; only existence is used.

## 6. Failure anatomy and scope

The first failed implication is

```text
eight freely tuned B-jets -> absorption of all 52 residual units          FALSE.
```

This closes every `B`-only clutch through degree eight in the current
THM-3212 chart. Degree nine is the first degree which changes the infinity
face, as
[THM-3237](THM-3237-degree-nine-jacobian-infinity-wall-and-square-root-escape.md)
proves. Its presently known strict-transform tunings still leave 50 or more
critical points, but this theorem does not classify arbitrary degree-nine
clutches.

The result is a critical-point obstruction for the explicit family `(2)`.
It supplies no marked inverse cover, branchwise cofactor, polynomial second
coordinate or Jelonek component. It proves neither `JC(2)` nor `DC(2)` and
does not classify deformations in `C_0` or `E_0`.

The next lawful probe must change a different gradient coordinate or supply
the true inverse-cover/cofactor sidecar. Another `B` retuning below degree
nine is exhausted.

## 7. Exact reproduction

Run

```text
python3 04-computation/jc_degree_at_most_eight_uniform_clutch_no_go_wildcard_20260803.py
python3 -O 04-computation/jc_degree_at_most_eight_uniform_clutch_no_go_wildcard_20260803.py
```

and compare LF-normalized bytes with the declared output. The companion uses
exact rational and cubic-field arithmetic only, has no assertion node,
floating literal, randomness or optimization-sensitive branch, and pins six
direct inherited artifacts.

QED.
