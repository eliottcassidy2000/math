---
id: THM-2822
title: "Sextic response centered-lift mod-three Faber obstruction"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Each of THM-2817's two sextic response carriers has a canonical
  centered polynomial exact-square lift.  In every constant-coefficient
  full even Faber bank, its first two flux conditions leave only rows
  j=2 mod 3 and force R_Q/q to lie in q^2 C[q^4], whereas the carrier
  requires R_Q/q=1.  Independently, the zero section makes a constant
  Jacobian impossible.  This excludes only the centered lift, not another
  multiplier or arbitrary Keller-chart entry of the abstract response.
source: root/sextic-centered-lift-faber-obstruction-2026-07-28
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2817-sextic-e3-maximal-pole-power-chebyshev-accessory-classification
related:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2811-special-image-linear-intertwiner-rigidity-and-jc-degree-wall
script: 04-computation/jc_sextic_centered_lift_mod3_faber_obstruction_thm2822.py
output: 05-knowledge/results/jc_sextic_centered_lift_mod3_faber_obstruction_thm2822.out
script_sha256: 818c1fd790483237cd88239d28dc118fd66898307129f92d3b2c4c68f4f581dd
output_sha256: 92fb93af04be369111eff473b23a6bcf77330405d3e1c4a873b05f13150444f4
hash_basis: LF-normalized bytes
---

# THM-2822 -- the canonical centered lifts miss every Faber degree

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2817 leaves two abstract sextic response carriers.  Both pass all
universal response tests, so there is no contradiction at that quotient.
There is nevertheless a natural polynomial source lift, and that lift is
impossible in two independent ways.  The more informative proof is an
all-degree mod-three support split in the full Faber gauge.

## 1. The two exact response packets

Write `T_pole` for the product of the four distinct pole factors, to avoid
confusing it with the quartic source invariant.  In both carriers

```text
E^2=D+A_0,                    C=-6A_0,

F=E^2/D,                     G=C E/(2D T_pole),
V=4D T_pole^2/C^2,           M=E T_pole.              (1)
```

For the power carrier, in coordinate `x`,

```text
D=x^3(x^3-1),                T_pole=x(x^3-1),
E=x^3-1/2,                   A_0=1/4,       C=-3/2,

V=(16/9)x^5(x^3-1)^3,
G=-3(2x^3-1)/[8x^4(x^3-1)^2],
M=x(x^3-1)(x^3-1/2).                                (2)
```

For the Chebyshev carrier, in the centered coordinate `y`, put
`T_3(y)=4y^3-3y`.  Then

```text
D=(y^2-1/4)^2(y^2-1),
T_pole=(y^2-1/4)(y^2-1),
E=T_3(y)/4,                    A_0=1/16,      C=-3/8,

V=(256/9)(y^2-1/4)^4(y^2-1)^3,
G=-3T_3(y)/(64D T_pole),
M=(T_3(y)/4)T_pole.                                  (3)
```

Direct substitution gives, in either row,

```text
F'=2G,                         F=VG^2,
2VG'+V'G=2,

P_resp:=VG=(2/C)M,
V-(4/C^2)M^2=-(4A_0/C^2)T_pole^2,
1/(1-F)=-D/A_0.                                      (4)
```

Thus the power and Chebyshev square defects have degree eight and are,
respectively,

```text
-(4/9)T_pole^2,               -(16/9)T_pole^2.       (5)
```

This records why neither carrier can be rejected by the response-side
tests of THM-2796.

On the quadratic deck put

```text
U=(2/C)T_pole sqrt(D),         R=E/sqrt(D).           (6)
```

The same identities give

```text
U^2=V,                         R^2=F,
R'=1/U,                        UG=R.                  (7)
```

## 2. The canonical centered polynomial lift

Choose the source multiplier to be the response carrier itself:

```text
A_src=P_resp=VG=(2/C)M,
H=Vz^2,                        L=A_src z,
P=H^2+L=V^2z^4+A_src z.                              (8)
```

With `w=Uz`, the normalized source polynomial is

```text
P=w^4+qw,                      q=A_src/U=UG=R.        (9)
```

Hence the quartic invariants are

```text
d=0,                           s=0,
T_F=q^2=F.                                           (10)
```

Normalize a hypothetical nonzero constant Jacobian to one.  THM-2784's
nonsplit transfer says that its full Faber representative has only even
rows and satisfies

```text
Phi_Q=0,                       Psi_Q in C,
R_Q'=1/U.                                            (11)
```

Both `R_Q` and `R` are deck-anti-invariant, so `(7)` and `(11)` leave no
additive constant:

```text
R_Q=R=q,                       K:=R_Q/q=1.            (12)
```

## 3. Exact mod-three support of every even Faber row

Index the even rows by

```text
m_j=4j-2,                      alpha_j=j-1/2,
j>=1.                                                (13)
```

At `(10)`, the coefficient series in THM-2760 becomes simply

```text
(1+q t^3)^alpha_j=sum_(n>=0)c_n t^n.                 (14)
```

Therefore

```text
c_n =
  binom(j-1/2,n/3) q^(n/3),    if 3 divides n,
  0,                           otherwise.             (15)
```

Using

```text
Phi_j=4c_(m_j+1),              Psi_j=4c_(m_j+2),
R_j=4c_(m_j+3)                                      (16)
```

gives the disjoint support table

```text
Phi_j !=0  iff j=1 mod 3,
  Phi_j=4 binom(j-1/2,(4j-1)/3) q^((4j-1)/3);

Psi_j !=0  iff j=0 mod 3,
  Psi_j=4 binom(j-1/2,4j/3) q^(4j/3);

R_j !=0    iff j=2 mod 3,
  R_j=4 binom(j-1/2,(4j+1)/3) q^((4j+1)/3).          (17)
```

Every displayed half-integral binomial coefficient is nonzero.  This is an
identity for every `j`, not a finite interpolation.

## 4. Full-bank contradiction

Let

```text
Q=sum_(j=1)^J a_j E_(4j-2),          a_J!=0,          (18)
```

be any constant-coefficient full even Faber bank.  A nonconstant function
on a curve is transcendental over `C`, so distinct powers of `q` are
linearly independent over `C`.

By `(17)`, `Phi_Q=0` kills every coefficient with `j=1 mod 3`.
Likewise `Psi_Q in C` kills every coefficient with `j=0 mod 3`; all its
displayed powers are positive, so the constant is actually zero.  Thus a
nonzero survivor uses only `j=3t+2`, and

```text
R_Q/q
 =sum_t a_(3t+2) b_t q^(4t+2)
 =q^2 H(q^4),                                        (19)
```

where every `b_t` is nonzero and `H` is a nonzero polynomial.  In
particular,

```text
R_Q/q !=1.                                           (20)
```

This contradicts `(12)`.  It also shows that any nonzero top row on this
root-zero boundary has `j=2 mod 3`, hence reduced degree

```text
4j-2=6 mod 12.                                       (21)
```

Equation `(21)` is compatible with THM-2760's exceptional projective lane;
the new information is that the complete constant-coefficient bank still
cannot realize the centered carrier.

## 5. Independent zero-section proof

There is a shorter proof that does not use the Faber gauge.  From `(8)`,

```text
P_x(x,0)=0,                       P_z(x,0)=A_src(x),  (22)
```

and `A_src` has degree seven for both carriers.  If a polynomial `Q`
satisfied `Jac(P,Q)=kappa!=0`, then at `z=0`,

```text
-A_src(x) Q_x(x,0)=kappa,                            (23)
```

which is impossible in `C[x]`.  Thus the canonical centered lift has no
polynomial constant-Jacobian mate, independently confirming `(20)`.

## 6. Exact boundary: what remains open

The conclusion is deliberately not an arbitrary chart exclusion.  In the
general polynomial exact-square source

```text
H=Vz^2+Bz+C_0,                  L=A_src z+E_0,
P=H^2+L,                                            (24)
```

the zero section is

```text
P_x(x,0)=2C_0 C_0'+E_0',
P_z(x,0)=2C_0B+A_src.                                (25)
```

A necessary Keller condition is only

```text
gcd(2C_0 C_0'+E_0', 2C_0B+A_src)=1.                 (26)
```

This boundary is real: already `C_0=0,E_0=x` makes the first entry in
`(26)` equal to one.  That is a control on the argument, not a Keller
example.

Changing the multiplier is equally load-bearing.  The response quotient
fixes only

```text
A_src K=(2/C)M,                 K=R_Q/q,              (27)
```

whereas

```text
d=C_0-B^2/(4V),
s=A_src B/(2V)-E_0,
T_F=A_src^2/V.                                       (28)
```

Thus another `A_src`, or nonzero `B,C_0,E_0`, changes the Faber trajectory
and need not obey `(14)`.  The first live inherited test begins with the
reduced degree-twenty-six bank.  A legal constant translation
`P_c=P+c` changes the next coefficient by
`c_22 -> c_22-(13/2)c`, so choose `c=2c_22/13`.  In that normalized
coordinate the bank is

```text
Q=E_26+c_18E_18+c_14E_14+c_10E_10+c_6E_6+c_2E_2,            (29)
```

with `(28)`, the two flux equations, `(27)`, polynomial reconstruction,
and `(26)` imposed simultaneously.  No symmetry of the abstract
power/Chebyshev quotient may be assumed to lift through these choices.

Therefore this theorem proves neither arbitrary entry of the two abstract
carriers nor `JC(2)` or `DC(2)`.

## 7. Exact companion

The companion verifies with exact rational and symbolic arithmetic:

1. all identities `(1)--(7)` for both carriers;
2. the centered lift and zero-section formulas;
3. direct series extraction against `(15)` for `j=1,...,30`;
4. the three disjoint support classes and the exponent pattern in `(19)`;
5. sample rows
   `K_2=-q^2/4`, `K_5=9q^6/512`, and
   `K_8=-195q^10/131072`; and
6. the lower-term gcd boundary after `(26)`.

It contains no Python `assert` node.  Run

```text
python 04-computation/jc_sextic_centered_lift_mod3_faber_obstruction_thm2822.py
python -O 04-computation/jc_sextic_centered_lift_mod3_faber_obstruction_thm2822.py
```

Normal, optimized, and stored transcripts agree exactly.
