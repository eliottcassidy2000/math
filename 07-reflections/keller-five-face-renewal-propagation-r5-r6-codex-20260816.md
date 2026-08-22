# The two hybrid faces propagate: `R_5` and `R_6` inherit complete packets

**Status: PROVED in THM-3522 + VERIFIED-SYMBOLIC + TWICE INDEPENDENTLY
HOSTILE-AUDITED.**  The exact companion is
`04-computation/keller_five_face_renewal_propagation_probe_20260816.py`.
The `R_5` application uses proved THM-3506 and THM-3513.  The `R_6`
application additionally uses proved THM-3521.  No image-prime,
irreducibility, degree-`243`, all-level polynomiality, or general Jacobian
claim is made.

## 1. The apparent missing step

THM-3506 defines the five-face packet `A(e,m)` and proves that, for

```text
Q=L^e N(P),
e'=7e-2m,                 m'=3e-2m,                    (1)
```

the first three faces propagate whenever `Q` is polynomial.  It leaves the
top-`z` and minimum-`gamma` renewal faces as separate obligations.  THM-3513
then proves those two faces for the fixed

```text
G=L^43N(J),               (e,m)=(271,99).              (2)
```

The key observation is that THM-3513's hybrid calculation is not special to
the 66,146-term polynomial `J`.  Its only polynomial-specific input is one
monomial that is simultaneously the unique top-`z` term and the unique
lowest residual-weight term.  Every complete `A(e,m)` packet supplies exactly
such a monomial.

## 2. The unique common monomial

Put

```text
a=2e-4m/3,       b=2e-2m/3,       h=e-2m/3.            (3)
```

The two renewal faces of `P` are

```text
in_max-z(P)=c x^a z^b,
in_min-gamma(P)=d z^e(27x^2z+y^3)^h,                  (4)
```

with nonzero rational `c,d`.  The term of the second face indexed by
`r=h` is

```text
d 27^h x^(2h) z^(e+h)=d 27^h x^a z^b.                 (5)
```

It is the unique term of that face having `z`-degree `b`; completeness of
both faces therefore implies

```text
gamma >= -8e+2m,          k <= b,                      (6)
```

and equality in both inequalities occurs only at `x^a z^b`.  Consequently

```text
delta_6=gamma-k >=(-8e+2m)-b,
delta_8=gamma-3k>=(-8e+2m)-3b                         (7)
```

with a unique equality monomial in each case.  This is precisely the
singleton gate used in THM-3513.

## 3. Top-`z` propagation

Use THM-3513's exact inverse scaling

```text
(x,y,z)_target=(A,B,Cs^3),
(x,y,z)_inverse~(q/s,-s/q,-Cs^6/q^3),
27A^2Cq^3-2=0,
L~27A^2C^2s^6.                                         (8)
```

Write

```text
r=a-3b=-4e+2m/3.                                       (9)
```

The unique source monomial gives on each inverse sheet

```text
P~c(-1)^b C^b q^r s^(-a+6b).                          (10)
```

Since the product of the three roots in (8) is
`2/(27A^2C)`, multiplication over the three sheets and by `L^e` gives

```text
Q(A,B,Cs^3)
 ~c^3(-1)^(3b) 27^e(2/27)^r
   A^(2e-2r) C^(2e+3b-r) s^(6e-3a+18b).               (11)
```

Direct substitution of (1) and (3) gives

```text
2e-2r       =2e'-4m'/3,
2e+3b-r     =2e'-2m'/3,
6e-3a+18b   =3(2e'-2m'/3).                             (12)
```

The scalar in (11) is nonzero.  Because `Q` is polynomial, the complete
top-`z` face is therefore

```text
in_max-z(Q)
 =c_z x^(2e'-4m'/3) z^(2e'-2m'/3),       c_z!=0.       (13)
```

No coefficient ledger for `Q` is needed.

## 4. Minimum-`gamma` propagation and the nonmonic factor

Use the second exact scaling

```text
(x,y,z)_target=(At,B/t,C/t^5),
(x,y,z)_inverse~(qt,-1/(qt),-C/(q^3t^8)),
Dq^3-3Bq-2=0,       D=27A^2C+B^3,
L~CDt^-8.                                                (14)
```

The same singleton gives

```text
P~c(-1)^bC^bq^rt^(a-8b).                               (15)
```

Here the leading coefficient of the residual cubic is load-bearing:

```text
product(q)=2/D.                                         (16)
```

Thus

```text
Q(At,B/t,C/t^5)
 ~c^3(-1)^(3b)2^r
   C^(e+3b)D^(e-r)t^(-8e+3a-24b).                     (17)
```

Again the exponents simplify exactly:

```text
e+3b          =e',
e-r           =e'-2m'/3,
-8e+3a-24b    =-8e'+2m'.                              (18)
```

Since `27x^2z+y^3` evaluates to `Dt^-3`, equation (17) is exactly

```text
in_min-gamma(Q)
 =c_gamma z^e'(27x^2z+y^3)^(e'-2m'/3),   c_gamma!=0.   (19)
```

Equations (13) and (19) prove the renewal propagation lemma:

> If `P` has the two complete renewal faces of `A(e,m)` and
> `Q=L^eN(P)` is polynomial, then `Q` has the two renewal faces of
> `A(e',m')`.

Combined with THM-3506's first-three-face transform, a full packet propagates
whenever the next cleared norm is polynomial.

## 5. Exact calibration against `J -> G`

The source singleton used by THM-3513 is

```text
c_Jx^66z^76,          c_J=2^15 3^171,
e=43,                 r=66-3*76=-162.                  (20)
```

Equations (11) and (17) give respectively

```text
c_z=c_J^3 27^43(2/27)^(-162)=3^1128/2^117,
c_gamma=c_J^3 2^(-162)=3^513/2^117,                    (21)
```

exactly reproducing both independently proved THM-3513 scalars.  This is a
normalization check on the general calculation, not merely an exponent
comparison.

## 6. The complete `R_5` packet

THM-3513 proves that `G` has `A(271,99)`, and THM-3506 proves that

```text
R_5=L^271N(G)                                           (22)
```

is polynomial.  Therefore the lemma applies unconditionally to this fixed
`R_5`.  The new state and renewal exponents are

```text
(e_5,m_5)=(1699,615),
max deg_z(R_5)=2988,
min gamma(R_5)=-12362.                                  (23)
```

The common monomial of `G` has coefficient `3^1128/2^117`, so the exact
faces are

```text
in_max-z(R_5)
 =(3^7251/2^1369)x^2578z^2988,                         (24)

in_min-gamma(R_5)
 =(3^3384/2^1369)z^1699(27x^2z+y^3)^1289.             (25)
```

The coefficient of the common term in (25) is
`3^(3384+3*1289)/2^1369=3^7251/2^1369`, agreeing with
(24).  Thus `R_5` has the complete packet `A(1699,615)`; its renewal faces
are not an additional open computation once (22) and THM-3513 are combined.

## 7. The complete `R_6` packet

THM-3521 proves

```text
R_6=L^1699N(R_5)                                        (26)
```

is polynomial and `L`-coprime.  The same lemma therefore gives

```text
(e_6,m_6)=(10663,3867),
max deg_z(R_6)=18748,
min gamma(R_6)=-77570,                                  (27)
```

with exact faces

```text
in_max-z(R_6)
 =(3^46008/2^10493)x^16170z^18748,                     (28)

in_min-gamma(R_6)
 =(3^21753/2^10493)z^10663(27x^2z+y^3)^8085.           (29)
```

The two successive Cassini determinants are

```text
det((271,99),(1699,615))=-1536,
det((1699,615),(10663,3867))=12288=-8(-1536),          (30)
```

as required by the norm-face matrix determinant.

## 8. Scope and next gate

The propagation theorem changes the tower bottleneck.  Renewal is no longer
a fresh coefficient problem at every rung: once one full packet is known,
the only new algebraic gate needed to propagate another full packet is
polynomiality of the next cleared norm.  For the fixed tower, THM-3506 and
THM-3521 pay that gate through `R_6`.

What remains open is substantial:

- `R_5` and `R_6` are not proved irreducible or image equations;
- no fifth image component, image multiplicity, or degree-`243` separability
  follows;
- polynomiality and finite-sheet nonvanishing for the next norm are open;
- exact positive discriminant multiplicities and all-level induction remain
  open;
- none of this proves a counterexample to `JC(2)` or classifies Keller maps.

At the time of this derivation, the next sharp tower test was the finite-sheet
unit for `L^10663N(R_6)`, not another renewal-face extraction.  **Subsequent
status:** THM-3523 has now paid that gate, constructing polynomial,
`L`-coprime `R_7=L^10663N(R_6)` and transporting the complete packet
`A(66907,24255)`.  The corresponding next polynomiality test is the
finite-sheet unit needed for `L^66907N(R_7)`; image, irreducibility, and
all-level claims remain open.

## 9. Reproduction

Run

```text
python -B 04-computation/keller_five_face_renewal_propagation_probe_20260816.py
python -B -O 04-computation/keller_five_face_renewal_propagation_probe_20260816.py
```

The companion verifies every general linear-form identity, the unique-face
intersection, the THM-3513 scalar calibration, both new packets, both exact
coefficient pairs, and the Cassini determinants.  Two independent companions
then rederived the inverse ratios and nonmonic Vieta factors and supplied the
hostile audit used to promote THM-3522.
