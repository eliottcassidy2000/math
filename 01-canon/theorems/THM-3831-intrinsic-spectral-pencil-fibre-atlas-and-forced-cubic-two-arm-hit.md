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
  verified the quadratic and cubic reconstruction formulas including D=0
  points, checked comaximality and both W signs, and audited the global sign
  transfer under plane base change.  The 34-gate exact companion checks the five
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
script_sha256: 4dc7892f2e36d467c68e8268bbce48125e41fa968f913781a7ddc61c0dfd7ae9
output_sha256: 79217dcd7e35d491266ff3a0d2d22b6502e068dbf4b98d4f1bd37854961633d9
semantic_sha256: 292840fe6e1ee62b404f211b32ad51cb69a1f74078981d3f53346370ab0e0f9f
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

Here `a(9a+14)` is nonzero at both roots of `Q`.  Formulas `(4)` express
`A,omega,theta` in `K[k,k^-1]`, so there is a surjection

```text
K[k,k^-1] -> B_a.                                               (14)
```

Conversely, substitute `(13)` and `(4)` into all three original cubic laws
and the definition of `D`.  They vanish modulo `Q(a)`, so the assignments
define an inverse to `(14)`.  Therefore

```text
B_a = K[k,k^-1]                    if Q(a)=0.                    (15)
```

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

As in the quadratic case, `(8),(19),(4)` express every generator of `B_a`
in the ring on the right below.  Conversely, the universal assignments

```text
h=ak,       C=z/k,       m=(z-1)/(ak),
D=(1+2z)/(Q(a)k^2)                                           (22)
```

satisfy all original cubic laws and reconstruct the different modulo
`P(a)` and `(21)`.  Hence there is an exact isomorphism

```text
B_a = K[k,k^-1,z]/(z(k+a^2z)).                                 (23)
```

The ideals `(z)` and `(k+a^2z)` are comaximal because their difference after
multiplying `(z)` by `a^2` is the unit `k`.  The Chinese remainder theorem
therefore gives

```text
B_a = K[k,k^-1] x K[k,k^-1].                                  (24)
```

The two intrinsic components are

```text
U_a^-: z=0,                  C=0,
        m=-1/(ak),           D=1/(Q(a)k^2);                    (25)

U_a^+: z=-k/a^2,             C=-1/a^2,
        m=(-k/a^2-1)/(ak),   D=(1-2k/a^2)/(Q(a)k^2).           (26)
```

Each is exactly `G_m`; formulas `(22)` and `(4)` give the remaining
coordinates and show that `(25),(26)` exhaust the full fibre, including any
points with `D=0`.  Finally `(9)` gives the exact labels

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
