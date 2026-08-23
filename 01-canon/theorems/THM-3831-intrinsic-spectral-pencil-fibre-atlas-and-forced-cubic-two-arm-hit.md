---
id: THM-3831
title: "Intrinsic spectral-pencil fibre atlas and forced cubic two-arm hit"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The five THM-3827
  spectral fibres are classified inside the actual
  nonlinear cubic surface.  Each of the two roots of 7a^2+3 has one Laurent
  G_m component and only the minus square-root sign.  Each of the three roots
  of 3a^3+7a^2+1 has exactly two comaximal Laurent G_m components, C=0 and
  C=-1/a^2, carrying the minus and plus signs.  Consequently, in the
  irreducible-h branch of any plane atlas, the mandatory disconnected fibre
  occurs at a cubic slope and the atlas meets both intrinsic components.  No
  plane atlas or Keller pair is constructed.
source: jc_quartic_c3_construct / THM-3827 disconnected spectral-fibre lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary, 2026-08-23).
  The audit rederived the quotient presentations, checked that k is a unit,
  and found that lift-law reconstruction alone did not certify saturation at
  D=0.  The repaired proof uses THM-3811's actual q-root chart, proves that
  every nonzero spectral fibre avoids its missing companion divisor, and
  classifies the full localized chart including D=0 points.  It also checked
  comaximality, both W signs, and the global sign transfer under plane base
  change.  The strengthened exact companion checks the five
  slopes and every required unit denominator; derives the cubic factor
  z(k+a^2 z); reconstructs all three original Delone--Faddeev multiplication
  laws, both lift laws, the SL2 law, and the different in the universal cubic
  product and the quadratic Laurent chart; verifies comaximality; and labels
  both cubic components and the quadratic component by the normalized
  discriminant root.  Normal and optimized replay agree with the frozen
  output.  Normal and optimized runs byte-match the frozen transcript and the
  raw hashes agree.
depends_on:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
related:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3830-coordinate-cross-bichromatic-split-nonentry
script: 04-computation/jc2_intrinsic_spectral_pencil_fibre_atlas_thm3831.py
output: 05-knowledge/results/jc2_intrinsic_spectral_pencil_fibre_atlas_thm3831.out
script_sha256: 9d307f8e9089962aca7cb3af1d1db35cf671c2b64afc7ded67785f74bf67631d
output_sha256: 92fd33bb830f81293206762b4965e33a66f3084e36ec60f244bfe0bb1f9b7722
semantic_sha256: 990a421380ed0135caf5d8b99630c2fedc2610dc4f0f34bb073eafe9c800f9ae
hash_basis: raw LF bytes
---

# THM-3831 -- the mandatory sign split is one of three cubic two-arm fibres

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.
Let `U=Spec B` and use the intrinsic functions of THM-3822, writing `k` for
its `k_0`.  Thus

```text
Ck-mh=1,                                                        (1)
D(7h^2+3k^2)=1+2Ck,                                            (2)
hD(9h+14k)=km+3hC^2,                                           (3)
```

and

```text
A=hD,             omega=kD,             theta=(m-14A)/3.       (4)
```

The exact converse in THM-3822 says that `(1)--(4)` reconstruct the three
original cubic multiplication laws and the different
`D=Comega-3Atheta-14A^2`.

Put

```text
Q(a)=7a^2+3,                    P(a)=3a^3+7a^2+1,
b(a)=(a+1)(2a+1)(3a-1).                                      (5)
```

The two roots of `Q` and the three roots of `P` are distinct and nonzero.
They are exactly the five spectral slopes from THM-3827.  For any one of
them define

```text
B_a=B/(h-ak),                         z=Ck.                     (6)
```

## 1. Every spectral fibre has an invertible Laurent parameter

In `B_a`, equation `(1)` reads

```text
k(C-am)=1.                                                       (7)
```

Thus `k` is a unit and

```text
C=z/k,                       m=(z-1)/(ak).                     (8)
```

This unit is what prevents hidden affine-line or vertex components in the
classification below.

There is also an intrinsic sign label on every spectral fibre.  Let

```text
q=7h^2+3k^2,
B_0=k^5-7h^2k^3-3h^2k^2-6h^3k^2-7h^4,
W=(h^2q^2D+B_0)/k.                                             (9)
```

Division by `k` is legal in `B_a`.  The THM-3822 quadratic in `D` gives

```text
W^2=H(h,k).                                                     (10)
```

At a spectral slope, the completed square of THM-3827 has `A_5=0`, hence

```text
W^2=(kB_3)^2,
B_3=(h+k)(2h+k)(3h-k).                                        (11)
```

The coefficient `b(a)` is nonzero at all five slopes, so `kB_3=b(a)k^4` is
a unit on `B_a`.  The signs `W=+kB_3` and `W=-kB_3` are therefore intrinsic
component labels.

The lift laws by themselves must not be treated as a presentation at `D=0`.
To classify the actual quotient, use the first THM-3811 root chart

```text
T=K[A,C,u]/(A G(u)-C(C+u^2)),       G(u)=u^3+7u+3,
J=A(3u^2+7)-2Cu.                                             (11a)
```

Its open `Spec T[J^-1]` lies in `U` and has

```text
D=AJ,                 h=1/J,                 k=u/J.           (11b)
```

The only companion divisor missed by this chart is `P_1`.  On `P_1` one has
`h=0` and `k=C^-1`, so `h-ak=-a/C` is a unit because every spectral slope is
nonzero.  Thus the **entire** fibre `B_a` lies in `(11a)--(11b)`.  There
`h=ak` is exactly `u=1/a`, and specialization gives

```text
A P(a)/a^3=C(C+1/a^2),
J=A Q(a)/a^2-2C/a.                                          (11c)
```

These equations, with `J` inverted, are the saturated scheme-theoretic fibre
presentation used below.  In particular, they retain rather than guess the
points where `D=AJ` vanishes.

## 2. The two quadratic slopes have one minus component

Assume

```text
Q(a)=0.                                                         (12)
```

Equation `(2)` immediately forces `z=-1/2`.  Equations `(1)` and `(3)` then
give

```text
C=-1/(2k),                    m=-3/(2ak),
D=(14k+3)/(4(9a+14)k^3).                                     (13)
```

Here `a(9a+14)` is nonzero at both roots of `Q`.  More decisively, `(11c)`
has `J=-2C/a`, so localizing at `J` makes `C` a unit, while `P(a)!=0`
determines `A` uniquely.  Hence the actual root chart gives

```text
B_a=K[C,C^-1]=K[k,k^-1]            if Q(a)=0.                    (15)
```

The displayed formulas `(13)` and `(4)` recover every generator and include
the point `C=-1/a^2`, where `A=D=0` and `k=a^2/2`.  Thus `(15)` has not lost
an exceptional graph point.

Direct substitution in `(9)` gives

```text
W=-kB_3.                                                        (16)
```

Thus each quadratic spectral fibre is one `G_m`, is connected even scheme-
theoretically, and carries only the minus sign.

## 3. The three cubic slopes have two opposite-sign components

Now assume

```text
P(a)=0.                                                         (17)
```

Then `Q(a)!=0` and

```text
b(a)=-Q(a).                                                     (18)
```

Equation `(2)` and `(8)` give

```text
D=(1+2z)/(Q(a)k^2).                                            (19)
```

Substitution into `(3)`, followed only by multiplication by the displayed
units, gives

```text
k[N_0(a)+N_1(a)z]-3a^2Q(a)z^2=0,                              (20)
N_0=3P,                     N_1=3b.
```

Using `(17),(18)`, equation `(20)` becomes

```text
z(k+a^2z)=0.                                                    (21)
```

This compatibility calculation is a useful check, but exactness comes from
the saturated chart `(11c)`: since `P(a)=0`, it gives

```text
C(C+1/a^2)=0.                                                 (21a)
```

The ideals are comaximal because `(1+a^2C)-a^2C=1`.  On `C=0`, inverting
`J=AQ(a)/a^2` leaves `K[J,J^-1]`; on `C=-1/a^2`, inverting
`J=AQ(a)/a^2+2/a^3` leaves the same Laurent ring.  Since `k=1/(aJ)`, these
are both `K[k,k^-1]`.  Equivalently, with `z=Ck`, `(21a)` is `(21)`, and
there is an exact isomorphism

```text
B_a=K[k,k^-1,z]/(z(k+a^2z))
   =K[k,k^-1] x K[k,k^-1].                                  (23)
```

The universal assignments

```text
h=ak,       C=z/k,       m=(z-1)/(ak),
D=(1+2z)/(Q(a)k^2)                                           (22)
```

satisfy all original cubic laws and reconstruct the different modulo
`P(a)` and `(21)`.  The root-chart proof above is what makes this converse
exhaustive at `D=0`; the lift-law calculation alone would not.

The two intrinsic components are

```text
U_a^-: z=0,                  C=0,
        m=-1/(ak),           D=1/(Q(a)k^2);                    (25)

U_a^+: z=-k/a^2,             C=-1/a^2,
        m=(-k/a^2-1)/(ak),   D=(1-2k/a^2)/(Q(a)k^2).           (26)
```

Each is exactly `G_m`; formulas `(22)` and `(4)` give the remaining
coordinates.  The second chart component contains `A=D=0`, `k=a^2/2`, so the
root-chart argument explicitly shows that `(25),(26)` exhaust the full fibre
including its exceptional point.  Finally `(9)` gives the exact labels

```text
W=-kB_3 on U_a^-,                  W=+kB_3 on U_a^+.            (27)
```

Thus the three cubic spectral fibres, and only those three among the five,
are intrinsically disconnected into the two square-root signs.

## 4. Consequence for a plane atlas

Let

```text
psi:A2 -> U                                                        (28)
```

be a dominant etale plane atlas and suppose the pulled-back `h` is
irreducible.  THM-3827 then supplies a spectral slope `a` for which the plane
fibre `h-ak=0` has nonempty plus and minus sign factors.

In the source function field, the polynomial square root used by THM-3827
equals either `psi^*W` or `-psi^*W`; the sign is global because the source
field is a domain.  A global swap does not change the assertion that both
signs occur.  Equations `(15),(16)` show that neither quadratic slope can
support both signs.  Hence

```text
P(a)=0,                                                         (29)
```

and the base-changed fibre meets both target components `(25),(26)`.  In
particular, a surviving atlas in the irreducible-`h` branch must cover points
on both the `C=0` and `C=-1/a^2` Laurent arms of one cubic spectral fibre.

This is a necessary image/incidence passport, not a construction.  The
existence of `(28)`, its Keller equation, and the reducible-`h` branch remain
open.  **QED.**
