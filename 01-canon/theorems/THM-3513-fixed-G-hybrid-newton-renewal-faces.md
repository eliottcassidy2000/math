---
id: THM-3513
title: "Fixed-G hybrid Newton renewal faces"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED.  For the fixed sporadic Keller polynomial
  G=L^43 N(J), the two renewal faces left open by THM-3506 are complete:
  max deg_z(G)=476 with face C_z x^410 z^476, and
  min(i-j-5k)=-1970 with face
  C_gamma z^271(27x^2z+y^3)^205, where both displayed scalars are explicit
  nonzero rationals.  Thus this fixed G has the full five-face packet
  A(271,99).  The next finite-sheet unit is left open here and subsequently
  closed by THM-3521.  THM-3522 subsequently proves fixed-chart renewal
  propagation and closes the complete packets of R_5 and R_6.  THM-3528
  subsequently proves the raw all-level polynomial-packet induction, and
  THM-3529 proves all complete packets are finite-sheet units.  A fifth image
  prime and every general Jacobian claim remain open at this theorem state.  The
  exact companion is independent of THM-3506's face script; an
  independent parent replay reproduced its semantic ledger exactly.  This
  is a replay audit, not a broader all-level or geometric proof audit.
source: codex/fixed-G-renewal/2026-08-16
depends_on:
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
related:
  - MISTAKE-413
  - MISTAKE-415
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
  - THM-3522-fixed-keller-five-face-renewal-propagation
script: 04-computation/keller_G_renewal_faces_independent_probe_20260816.py
output: 05-knowledge/results/keller_G_renewal_faces_independent_probe_20260816.out
script_sha256: f9e82f502026dfe499ebba9290295f98056d1b7dba7c893184d9871a032be01f
output_sha256: becaa80c075bd46e4193b216406c2152f3d5d8565f6116ba6db9b712409badaa
hash_basis: raw LF bytes; exact rational values use the ASCII numerator/denominator convention
---

# THM-3513 -- the two hybrid Newton limits renew the fixed `G` packet

**PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED.**

Let `F:C^3->C^3` be the fixed sporadic Keller map of THM-2473.  Use
THM-3495's normalization of the `66,146`-term polynomial `J`, and put

```text
G=L^43 N(J) in Q[x,y,z],
```

as in THM-3498.  For a monomial `x^i y^j z^k`, write

```text
gamma=i-j-5k.
```

Then

```text
max deg_z(G)=476,
in_max-z(G)=(3^1128/2^117) x^410 z^476,                (1)
```

and

```text
min gamma(G)=-1970,
in_min-gamma(G)
 =(3^513/2^117) z^271(27x^2z+y^3)^205.                (2)
```

Both are complete faces: there are no omitted equal-weight competitors.
Consequently the two renewal faces (6)--(7) in THM-3506 hold for this fixed
`G`, and `G` has the complete packet `A(271,99)`.

## 1. A singleton controls both limits

For a monomial of `J`, put

```text
delta_6=i-j-6k=gamma-k,
delta_8=i-j-8k=gamma-3k.                               (3)
```

Exact extraction from THM-3495's pinned coefficient ledger gives

```text
max k(J)=76,                 min gamma(J)=-314.         (4)
```

The complete `k=76` face is the singleton

```text
c_J x^66z^76,               c_J=2^15 3^171,            (5)
```

and the complete `gamma=-314` face has `34` terms.  Its intersection with
the face in (5) is exactly the monomial in (5).  Equations (3)--(4) therefore
give

```text
delta_6 >= -390,            delta_8 >= -542,           (6)
```

with equality in either inequality only at (5).  Independently, a direct
scan of all `66,146` terms gives the same two singleton initial forms.  The
next attained weights are `-389` and `-539`; hence no equal-weight term can
cancel (5) in either residual algebra.

## 2. The top `z`-face

Use target coordinates `(a,b,c)` temporarily, and scale

```text
(a,b,c)=(A,B,Cs^3),           w=q/s,        s -> infinity.
```

The leading inverse cubic and inverse coordinates are

```text
27A^2Cq^3-2=0,
(q_x,q_y,q_z)~(q/s,-s/q,-Cs^6/q^3).                   (7)
```

These identities follow by exact weighted extraction from the inverse
numerators of THM-3495.  The companion also reduces both coordinate
identities modulo the cubic in (7), so no choice of a formal root is used.

This scaling reads `delta_6`.  By (5)--(7), every inverse branch has

```text
J(q_x,q_y,q_z)
 ~c_J C^76q^-162s^390
 =c_J(27/2)^54A^108C^130s^390.                         (8)
```

The second equality uses `q^3=2/(27A^2C)`, so its coefficient is the same
on all three roots.  Since

```text
L^43~(27A^2C^2)^43s^258,                               (9)
```

the product of the three values in (8), multiplied by (9), is

```text
G(A,B,Cs^3)
 ~(3^1128/2^117)A^410C^476s^1428.                     (10)
```

THM-3498 proves that `G` is polynomial.  Therefore the powers of `s` in
this substitution are three times its `c`-degrees.  The complete generic
leading coefficient in (10) proves (1), after relabelling `(a,b,c)` as
`(x,y,z)`.

## 3. The gamma face and the nonmonic Vieta factor

Now scale

```text
(a,b,c)=(At,B/t,C/t^5),       w=qt,        t -> 0,
D=27A^2C+B^3.
```

Exact weighted extraction gives

```text
Dq^3-3Bq-2=0,
(q_x,q_y,q_z)~(qt,-1/(qt),-C/(q^3t^8)).                (11)
```

This scaling reads `delta_8`, so (5)--(6) imply that each branch contributes

```text
J(q_x,q_y,q_z)~c_JC^76q^-162t^-542.                   (12)
```

The cubic in (11) is nonmonic.  Retaining its leading coefficient, Vieta's
formula gives

```text
product(q)=2/D,
product(q^-162)=(D/2)^162.                             (13)
```

Also `L^43~(CD)^43t^-344`.  Taking the norm in (12) and using (13) gives

```text
G(At,B/t,C/t^5)
 ~(3^513/2^117)C^271D^205t^-1970.                     (14)
```

The power of `t` is precisely the target weight `i-j-5k`.  Since (14) is a
generic polynomial in `A,B,C`, it is the complete initial form, proving
(2).  In particular, the factor `D^205` is support data; it cannot be
discarded as a constant unit.

## 4. Consequence and boundary

For `(e,m)=(271,99)`, the renewal exponents in THM-3506 are

```text
2e-4m/3=410,       2e-2m/3=476,
e-2m/3=205,        -8e+2m=-1970.
```

Thus (1)--(2) are exactly the two faces missing from that theorem's fixed
`G` packet.  Together with THM-3506's three transported faces, they make the
five-face transform from `G` to

```text
R_5=L^271N(G)
```

lawful and give the first three output faces at `(1699,615)`.

This theorem itself does not determine `v_L(N(R_5))`: the two divergent
sheets are controlled by the exposed face, but the finite inverse sheet may
still vanish.  THM-3521 subsequently proves that finite sheet is a unit and
therefore closes the valuation.  THM-3522 then proves the fixed-chart
renewal implication and obtains the complete packets of both `R_5` and
`R_6`.  THM-3528 subsequently proves the raw all-level polynomial tower.
These results still prove no fifth image equation or prime, statement about
arbitrary Keller maps, or version of the general Jacobian conjecture.

Reproduce the exact ledger, both residual-cubic reductions, and both face
scalars with

```text
python 04-computation/keller_G_renewal_faces_independent_probe_20260816.py
python -O 04-computation/keller_G_renewal_faces_independent_probe_20260816.py
```

An independent parent replay passed with the same semantic ledger
`2887be137141414185ab305e7e0416754f73a009625924a5b6c6c1a268101dbd`.
This verifies reproducibility of the exact companion; it does not enlarge
the theorem's fixed-`G` scope.

**QED.**
