---
id: THM-2282
title: "Thirteen-adic saturation and a unit shallow-blocker minor"
status: >
  PROVED + INDEPENDENTLY AUDITED. Every one of the 120 interior
  first-depth-one scalar profiles has two bounded scalar relations whose
  coefficient minor on c_1 and one coordinate outside the THM-2266 pair is
  nonzero modulo thirteen. A saturation descent in the torsion-free quotient
  relation lattice bounds the second relation by 4921, the scalar minor by
  3725197, and its fixed-section original-row lift by 7450394. The outside
  coordinate is not prescribed, no exact Fourier atom or owner transition is
  selected, no scalar profile is excluded, and LRC(14) remains open.
source: codex-2026-07-25-thirteen-adic-unit-minor
depends_on:
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
  - THM-2279-shallow-blocker-anchored-relation-minor
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
---

# THM-2282 -- thirteen-adic saturation and a unit anchored minor

**PROVED + INDEPENDENTLY AUDITED.**

THM-2279 forces a nonzero integer minor anchored at the shallow blocker
`c_1`, but its outside coefficient may be divisible by thirteen. The
repair is a finite saturation descent:

```text
all outside coefficients divisible by 13
  -> subtract a bounded multiple of the primitive shallow pair
  -> divide the complete relation by 13
  -> retain independence and crossing
  -> repeat until an outside coefficient is a 13-unit.          (1)
```

Infinite descent is impossible in the torsion-free quotient of the relation
lattice by the primitive pair. The height recurrence remains uniformly
bounded.

## 1. The primitive shallow pair

Use the scalar row

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3)                     (2)
```

in one of the `120` interior first-depth-one profiles

```text
c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

3<=b<=c-2,                    5<=c<=19.               (3)
```

Let

```text
Lambda_*={x in Z^9:x.w_*=0}.                          (4)
```

THM-2266 and the normalization recorded in THM-2279 give a primitive
relation

```text
p in Lambda_*                                         (5)
```

of one of the forms

```text
13d H-a c_1=0,                                        (6)

d c_1-13a q_i=0.                                      (7)
```

In both cases

```text
gcd(a,d)=1,              13 does not divide ad,
ad<=757.                                                 (8)
```

Write

```text
A=supp(p)={c_1,z},        B=[9]\A,                    (9)
```

where `z=H` in (6) and `z=q_i` in (7). The load-bearing coefficient facts
are

```text
p_k=0                       for k in B,

13 divides p_z,

13 does not divide p_(c_1),

0<|p_(c_1)|<=757,

P:=||p||_infinity<=9841.                              (10)
```

The scalar value `z` is a thirteen-unit, while the scalar value `c_1` is
divisible by thirteen.

THM-2275's adaptive cut supplies an independent relation

```text
r_0 in Lambda_*\Qp,                 ||r_0||_infinity<=2116,

L_A(r_0):=sum_(j in A)(r_0)_j(w_*)_j!=0.             (11)
```

The sharper height is `708` in the guard-owner case, but `2116` is the
uniform input used below.

## 2. One descent step

Let

```text
r in Lambda_*\Qp                                     (12)
```

be any current relation. If some

```text
k in B
```

satisfies

```text
13 does not divide r_k,                               (13)
```

the descent stops.

Suppose instead that

```text
13 divides r_k                 for every k in B.      (14)
```

Reduce `r.w_*=0` modulo thirteen. Every term indexed by `B` vanishes
because its coefficient is divisible by thirteen. The `c_1` term vanishes
because `13|c_1`. The only remaining term is

```text
r_z z=0 mod 13.                                       (15)
```

Since `z` is a thirteen-unit,

```text
13 divides r_z.                                       (16)
```

The coefficient `p_(c_1)` is invertible modulo thirteen. Choose the unique
residue class `n mod 13` satisfying

```text
r_(c_1)-n p_(c_1)=0 mod 13                           (17)
```

and take its centered representative

```text
-6<=n<=6.                                             (18)
```

Equations (10), (14), (16), and (17) show coordinate by coordinate that

```text
r-np in 13 Z^9.                                       (19)
```

Define

```text
r'=(r-np)/13.                                         (20)
```

Then `r'` is an integral relation:

```text
r' in Lambda_*.                                       (21)
```

It remains independent of `p`. Indeed,

```text
r' in Qp
  -> r=13r'+np in Qp,                                 (22)
```

contradicting (12).

The crossing partial sum is also retained. Because `p` is supported on
`A` and is itself a relation,

```text
L_A(p)=p.w_*=0.                                       (23)
```

Therefore

```text
L_A(r')=L_A(r)/13.                                    (24)
```

Starting from (11), every relation in the descent continues to cross
`A|B`.

## 3. Termination in the quotient lattice

The vector `p` is primitive in `Z^9`. Hence

```text
Z^9/Zp
```

is a free abelian group. The inclusion `Lambda_* subset Z^9` induces an
injective map

```text
Q:=Lambda_*/Zp  ->  Z^9/Zp,                           (25)
```

because its kernel is

```text
(Lambda_* intersection Zp)/Zp=0.                     (26)
```

Thus `Q` is a finitely generated torsion-free abelian group, hence free.

At each descent step,

```text
[r]=13[r']                 in Q.                      (27)
```

If the process never stopped, the initial nonzero class would satisfy

```text
[r_0] in intersection_(N>=1) 13^N Q.                 (28)
```

Choose a basis of the free group `Q`. An integer coordinate divisible by
every power of thirteen is zero, so

```text
intersection_(N>=1) 13^N Q={0}.                      (29)
```

Equations (28)--(29) would give `[r_0]=0`, contrary to (11). Therefore
after finitely many steps there is a relation

```text
r_* in Lambda_*\Qp
```

and a coordinate `k in B` such that

```text
13 does not divide (r_*)_k.                          (30)
```

## 4. Uniform height invariant

Put

```text
R=||r||_infinity.
```

Equations (18), (20), and (10) give the recurrence

```text
||r'||_infinity<=(R+6P)/13.                          (31)
```

Use

```text
M=4921.                                               (32)
```

The initial height in (11) satisfies `2116<=M`, and

```text
P<=9841,
M>=ceil(P/2).                                         (33)
```

Consequently,

```text
(M+6P)/13<=M.                                         (34)
```

Induction through (31) proves

```text
||r_*||_infinity<=4921.                               (35)
```

The integer recurrence can sharpen `4921` by one, but no such refinement
is needed here.

## 5. The thirteen-unit anchored minor

Use the relation rows `(p,r_*)` and the columns `(c_1,k)` from (30).
Since `k in B`, equation (10) gives `p_k=0`, so

```text
Delta_(c_1,k)
 =det [[p_(c_1),0],[(r_*)_(c_1),(r_*)_k]]
 =p_(c_1)(r_*)_k.                                    (36)
```

Both factors on the right are thirteen-units by (10) and (30). Therefore

```text
Delta_(c_1,k)!=0 mod 13.                              (37)
```

Equations (10) and (35) give the explicit scalar bound

```text
0<|Delta_(c_1,k)|
 <=757*4921
 =3725197.                                            (38)
```

Thus the reductions of the two relation rows have rank at least two over
`F_13`, witnessed by a minor with the prescribed shallow-blocker column.

## 6. Fixed-section original-row lift

THM-2203 lifts scalar relation coefficients by

```text
(x_H,x_rest) |->(2x_H,x_rest)                         (39)
```

on the fixed nine-coordinate original-row section. This preserves
independence and multiplies the minor in (36) by two only when `k=H`.
Hence the original row has a minor on `c_1` and the corresponding outside
coordinate satisfying

```text
Delta_original(c_1,k)!=0 mod 13,

0<|Delta_original(c_1,k)|
 <=2*3725197
 =7450394.                                            (40)
```

The factor two is itself a thirteen-unit, so (37) survives. The lifted
relation rows have uniform height at most

```text
max(2*9841,2*4921)=19682.                             (41)
```

## 7. Frontier effect and exact stopping boundary

THM-2279's missing residue coordinate is now repaired:

```text
120 interior scalar profiles
  -> a fixed-section rank-two relation packet
  -> a nonzero c_1-anchored Plucker coordinate over F_13.       (42)
```

This is stronger than a nonzero integer determinant and stronger than an
unanchored rational rank statement. It supplies an algebraic transversality
sidecar for THM-2269's marked residue energy, THM-2273's gap ancestry, and
THM-2278's labelled root-fibre characters.

It does not yet identify an exact global Fourier atom. Restriction to a
successor gap mixes the diagonal integer Fourier modes inside one residue
class, as THM-2278 explicitly records. Nor does (37) prescribe the outside
label `k`, prove that its relation coefficient corresponds to a positive
ancestry return, or turn a root-fibre vector energy floor into a scalar
cover contradiction.

The connection and loss ledger is

```text
source:
  a primitive shallow pair with a unit c_1 coefficient and an
  independent adaptive crossing relation;

target:
  a bounded c_1-anchored relation minor nonzero over F_13;

map:
  whenever every outside coefficient vanishes mod 13, cancel the c_1
  residue with a centered multiple of the pair and divide by 13;

preserved:
  integrality, the scalar relation equation, independence, crossing,
  the shallow blocker label, and a uniform height invariant;

destroyed:
  the original crossing row, the number of descent steps, the terminal
  outside label, integer Fourier frequency, and root-owner ancestry;

cheapest hostile probes:
  all outside coefficients initially divisible by 13, maximal pair
  height P=9841, k=H in the fixed-section lift, and infinite descent;

needed sidecar:
  identify the terminal unit column with a labelled root character and
  prove that its exact frequency survives the gap/owner restriction.      (43)
```

No scalar profile is excluded, and LRC(14) remains open. QED.
