---
id: THM-3239
title: "Universal sixth resonant PRS strict-transform wall and Pell atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  After dividing the universal sixth reciprocal pivot by the exceptional
  square rho_4^2 from THM-3231, the strict transform factors into inherited
  walls and one new primitive irreducible W28 of degree 28.  W28 is
  irreducible modulo 73 and equals 1 modulo 2.  Reconstructing rho_6 gives
  Pell content 70 and normalization gauge 239, and every fixed offset retains
  a finite explicit exceptional-prime atlas through row six.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  An independent theorem audit reconstructed both denominator exponent tables
  and the full numerator factorization, checked that the Rabin tests at 73
  are exactly the prime-divisor tests for degree 28, verified the fixed-offset
  exceptional-prime update and the p=43 strict-transform control, and replayed
  normal and optimized companions byte-for-byte against the stored output.
  No theorem, typing, or scope defect was found.
depends_on:
  - THM-3223-fourth-fifth-resonant-prs-primitive-walls-pell-content-and-pivot-resurrection
  - THM-3231-fraction-free-simple-pivot-wall-second-order-blowup-carry
related:
  - THM-3217-universal-resonant-degree-prs-wall-atlas-and-fixed-offset-exception-set
  - THM-3229-hasse-pluecker-simple-root-contact-gcd-flag-and-degree-termination
script: 04-computation/factorial_universal_sixth_strict_transform_wall_thm3239.py
output: 05-knowledge/results/factorial_universal_sixth_strict_transform_wall_thm3239.out
script_sha256: e0dcdaf148f7e9aeb2101f5b18ae326a071b59333d128936539ea02f33f8e2dd
output_sha256: 7934a8e881fd8f28147a879f364df1def7839737ab721b40b0b49b564ba50508
hash_basis: LF-normalized bytes
---

# THM-3239 -- universal sixth resonant PRS strict-transform wall and Pell atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3223 factors the fourth and fifth selected pivots.  THM-3231 then proves
that the sixth row has a universal exceptional square: if `rho_k` is the
constant coefficient of `R_k`, then `rho_4^2` divides `rho_6`.  The natural
next coordinate is therefore the strict transform `rho_6/rho_4^2`, not the
special-fibre zero.  It has a clean universal wall of its own.

## 1. Sixth row and strict-transform coordinate

Retain the degree-`d` reciprocal pair and fraction-free operator of
THM-3223:

```text
R_(-1)=b,       R_0=a,       R_(k+1)=E(R_k,R_(k-1)),

E(f,g)=z^(-2){f_0^2g-[f_0g_0+P_1(f,g)z]f}.              (1)
```

Write

```text
r=R_1, s=R_2, t=R_3, w=R_4, x=R_5, y=R_6,
rho_k=[z^0]R_k.                                          (2)
```

Thirteen input top coefficients suffice for `rho_6`.  Since `rho_4` is a
nonzero rational function of `d`, define

```text
Theta_6(d)=rho_6(d)/rho_4(d)^2.                          (3)
```

By THM-3231, `(3)` is the coefficient of the Rees strict transform across
the `rho_4=0` divisor, rather than an ad hoc cancellation.

## 2. The primitive degree-28 wall

Define

```text
W28(d)=
 169920813940688814080000d^28
 +331141682207614360879104000d^27
 -134162410011258685536379863040d^26
 +5300802398351890701896525021184d^25
 -684739297310367200983522643804160d^24
 +130331904744222160513651684933632000d^23
 -12943819225240695078283631960879267840d^22
 +713250673128714558996347298753011515392d^21
 -24834141237082462639059202312165067325440d^20
 +598783632792278487019889843452862948966400d^19
 -10626718017854951523648995520000467801210880d^18
 +144696459098683056904745037387918164398964736d^17
 -1556104532530719261165315051927889255805747200d^16
 +13492893612749552869825791823397116965185126400d^15
 -95734068097296356726622150832929787026834718720d^14
 +561609982081860252702393234504629255248835248128d^13
 -2743000014814555875946490646195207597676182896640d^12
 +11199447544319326718253356230653228142356647116800d^11
 -38279428000372095326340982875930990156643565568000d^10
 +109404450823229543746513469845082648703787025301504d^9
 -260447513492258947342576045172684040607636400701440d^8
 +512815846815398540657575958180230914285371601715200d^7
 -825986183003965800574110651732582518260652206080000d^6
 +1070580568879390678018966410506355769337873203200000d^5
 -1089659124692158529207880465753484653373043380800000d^4
 +838910469802753648130728862980511764057677888000000d^3
 -459389089498986441643744053260439590725236765000000d^2
 +159486965048208251611930715459002067591134125000000d
 -26393359344902929074282529688291558301498193359375.    (4)
```

For the ordered odd shifts

```text
k=(27,25,23,21,19,17,15,13,11,9,7,5),                  (5)
```

put

```text
D_Q=2^112 product_i (2d-k_i)^(q_i),
q=(1,2,3,4,5,6,11,16,27,38,65,92),                     (6)

D_6=2^168 product_i (2d-k_i)^(e_i),
e=(1,2,3,4,7,10,17,24,41,58,99,140).                   (7)
```

Using THM-3223's `H,J,U,W13`, exact reduction in `Q(d)` gives

```text
Theta_6
 =-(d-1)^46 H^16 J^6 U^4 W28/D_Q,                       (8)

rho_6
 =-(d-1)^70 H^24 J^10 U^4 W13^2 W28/D_6.               (9)
```

Equation `(9)` follows either by multiplying `(8)` by THM-3223's exact
`rho_4^2`, or by direct reduction of the sixth row.  Seven independent raw
moment expansions at every integer `14<=d<=20` reproduce both `(8)` and
`(9)` exactly.

## 3. W28 is a genuine fixed-offset wall

The polynomial `W28` is primitive.  Rabin's criterion modulo `73` gives

```text
X^(73^28)=X                                      mod W28,

gcd(W28,X^(73^14)-X)=gcd(W28,X^(73^4)-X)=1.       (10)
```

The exponents `14` and `4` correspond to the prime divisors `2,7` of `28`.
There is no degree drop modulo `73`, so `(10)` proves irreducibility over
`F_73`; Gauss's lemma proves irreducibility over `Q`.

Moreover,

```text
W28(d)=1                                           (mod 2). (11)
```

Every nonconstant coefficient in `(4)` is even and the constant coefficient
is odd.  Hence `W28` never vanishes at an integer.  It is a new primitive
chart divisor, but it creates no infinite fixed-offset degeneracy.

## 4. Pell content and gauge continue

Equation `(9)` gives the exact content orders through row six:

```text
ord_(d-1)(rho_1,...,rho_6)=1,2,5,12,29,70.               (12)
```

These are the first six positive Pell numbers.  No hidden `(d-1)` factor is
present in `W28`, because the factorization in `(8)` is primitive and
coprime.

Under common rescaling of the two input rows, the bidegree law
`E(cf,dg)=c^2dE(f,g)` gives

```text
1,1,3,7,17,41,99,239.                                   (13)
```

Thus the sixth row has gauge exponent `239`.  The Pell content and companion
gauge mechanisms observed through row five are exact through row six; no
all-depth recurrence for the primitive factors is asserted.

## 5. Every fixed offset has a finite row-six atlas

Let `Xi_s^[5]` be THM-3223's nonzero integer controlling all selected
coordinates through row five.  Define

```text
Xi_s^[6]
 =Xi_s^[5] W28(s)(2s-25)(2s-27).                         (14)
```

For every integer `s>=2`, equation `(11)` and parity of the two new linear
factors give

```text
Xi_s^[6] !=0.                                            (15)
```

Fix `s>=2`, put `d=p+s`, and suppose

```text
p>max(13,2s-4),                     p does not divide Xi_s^[6]. (16)
```

THM-3223 controls the earlier coordinates.  Equations `(7),(9)` show that
every numerator and denominator of `rho_6` reduces to a unit at `d=s`, and
the common normalization change has gauge exponent `239` from `(13)`.
Therefore the standard sixth pivot is also a `p`-adic unit.  Each fixed
offset has only the explicit finite exceptional-prime set through row six.

The condition is sufficient, not minimal: `Xi_s^[6]` deliberately retains
all earlier row factors so that one integer controls the whole atlas.

## 6. The p=43 wall is the strict-transform positive control

At `s=2,p=43`, `W13(2)=0 mod 43`, so the standard sixth pivot in `(9)` has
the inherited square wall `W13^2`.  THM-3231 proves the sharper local data

```text
v_43(rho_4)=1,                  v_43(rho_6)=2,

(rho_6/rho_4^2) mod 43=39,      rho_6/43^2 mod 43=19.    (17)
```

Thus `(8)` is exactly the coordinate that survives the apparent terminal
special fibre.  The new primitive `W28` does not vanish there.  This is the
sharp distinction between a standard scalar pivot atlas and its blowup:
`rho_6` dies to order two, while its strict transform is a unit.

## 7. Connection contract and scope

```text
source:       universal degree-d resonant factorial reciprocal pair;
operation:    sixth fraction-free row and exceptional-square division;
target:       primitive strict-transform divisor W28;
preserved:    exact rational factorization, p-adic valuation, fixed offset;
destroyed:    the remaining sixth row and any canonical root/face selector;
sidecar:      THM-3231 normal-cone coordinate on inherited W13 walls.
```

The theorem advances the selected-pivot arithmetic atlas by one row.  It
does not prove that `W28` controls the whole exterior state, that primitive
walls remain irreducible at every depth, or that a finite row always
separates arbitrary radial coefficients.  It supplies no new `NC(2)`,
`GMC(2)`, `JC(2)`, or `LRC(14)` conclusion.

## 8. Exact evidence

Run

```text
python 04-computation/factorial_universal_sixth_strict_transform_wall_thm3239.py
python -O 04-computation/factorial_universal_sixth_strict_transform_wall_thm3239.py
```

and compare LF-normalized bytes with the declared output.  The companion
hash-pins promoted THM-3223, promoted THM-3231, and the THM-3223 generator;
constructs thirteen top jets and all six rows in `Q(d)`; proves `(8),(9)`;
checks fourteen raw-moment coordinates by an independent integer expansion;
verifies the Rabin certificate, parity, Pell content, gauge, and both
denominator exponent tables; and uses no floating point, randomness,
optimizer, or assertion-sensitive check.

QED.
