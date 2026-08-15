---
id: THM-3441
title: "The weighted quartic has one identity and one transposition Jelonek component"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For THM-3438's quartic Keller map G, the exact Jelonek set is V(C) union
  V(L), with L an explicit irreducible surface.  Generically C=0 loses one
  sheet with identity inertia and even primitive clearing exponent 12, while
  L=0 loses two sheets with transposition inertia and odd exponent 5.  This
  refutes HYP-9027's Keller successor claiming odd inertia on every Jelonek
  component, while verifying parity=inertia-sign for both components.  It
  does not classify quartic Keller boundaries.  THM-3448 subsequently
  realizes genuine Keller C3 inertia in degree five and excludes it only
  inside the weighted quartic family.
source: Socrates independent quartic boundary audit, integrated by codex2, 2026-08-15
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
related:
  - THM-3448-weighted-keller-cyclic-jelonek-inertia-family
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
script: 04-computation/jc_weighted_quartic_jelonek_inertia_thm3441.py
output: 05-knowledge/results/jc_weighted_quartic_jelonek_inertia_thm3441.out
script_sha256: adace19f41b7816d273a47556013253e5d05014f517c0212743e9d74dcc5093f
output_sha256: 8e9629df938785d0c25b0a2ba35c29698c992bb8c785351055c01f58bad73a04
semantic_sha256: 84d29019ac0b6cc51caf6f87f75de74c5b655f1d8a0254c259e44d9b224290cf
transcript_payload_sha256: 6e88ebd254430d90958b954d17308527b7860fe6c5fa04602bb5d83040e68cb9
---

# THM-3441 -- the weighted quartic has one identity and one transposition Jelonek component

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. The discriminant and the chart-index trap

For THM-3438's quartic map `G`, write target coordinates `(A,B,C)` and

```text
P=BC,                 Q=AC^2,
f(w)=w^4-w^2+2Pw-Q.                                  (1)
```

Direct calculation gives

```text
Disc_w(f)=-16D(P,Q),
D(P,Q)=27P^4-36P^2Q-P^2+16Q^3+8Q^2+Q.               (2)
```

After substitution,

```text
D(BC,AC^2)=C^2 L(A,B,C),                              (3)

L=16A^3C^4+8A^2C^2-36AB^2C^2
  +A+27B^4C^2-B^2.                                   (4)
```

The target change `(A,B,C)->(P,Q,C)` has Jacobian `-C^3`; it is birational
only over `C!=0`.  Thus `C^2` in `(3)` cannot be read automatically as finite
ramification.  At `C=0`,

```text
f(w)=w^2(w-1)(w+1).                                  (5)
```

The double `w=0` is a collision of primitive-coordinate values in a
nonmaximal monogenic order.  Since `det JG=-6`, no factor of `(3)` is the image
of a finite critical locus.

The surface `L=0` is irreducible.  Over `C!=0` it is the pullback of the
irreducible discriminant curve parametrized by

```text
P=r-2r^3,             Q=r^2-3r^4.                    (6)
```

That parametrization is birational; generically

```text
r=(12PQ-P)/(18P^2-4Q-1).                              (7)
```

Moreover `partial_A L=1` at `C=0`, so the closure introduces no hidden power
or additional component.

## 2. Exact finiteness locus

For a root of `(1)`, put

```text
gamma=P-w+2w^3=(1/2)f'(w).                            (8)
```

The source reconstructs as

```text
x=C/gamma,
y=(w-gamma)/(3C),
z=gamma(7gamma-3gamma^2-4w)/(3C^2).                   (9)
```

On the complement of `V(CL)`, both `C` and `Disc_w(f)` are units.  Hence
`gamma` is a unit in the finite quartic algebra defined by `(1)`, and `(9)`
is regular there.  This proves

```text
S_G subseteq V(C L).                                  (10)
```

Both components are actually asymptotic.

### The component `C=0`

Fix any target `(A_0,B_0,0)` and a parameter `s`.  Set

```text
w_s=-1-B_0s+(3B_0^2-A_0)s^2/2,
P_s=B_0s,
Q_s=w_s^4-w_s^2+2P_sw_s,
(A_s,B_s,C_s)=(Q_s/s^2,B_0,s).                        (11)
```

Then the target tends to `(A_0,B_0,0)`, `f(w_s)=0`, and `(9)` gives

```text
x_s~-s,              y_s->B_0,             z_s~2/s^2.
```

Thus every point of `C=0` is nonproper.  At a generic point `A!=B^2`, the
four roots have integral-power expansions

```text
w_-~-1-BC,           w_+~1-BC,
w_(1,2)~C(B +- sqrt(B^2-A)).                           (12)
```

A meridian therefore permutes no roots: generic infinity inertia is the
identity.  Exactly the `w_-` sheet escapes; the other three remain finite.
For the control target `(A,B,C)=(0,1,0)`, the finite preimages are

```text
(0,1,-29/2),          (1,-1/3,7/3),          (-1,1,5). (13)
```

### The component `L=0`

On a dense open, take `C=c!=0` and `(P,Q)` from `(6)`, so

```text
B=(r-2r^3)/c,         A=(r^2-3r^4)/c^2,
f(w)=(w-r)^2(w^2+2rw+3r^2-1).                         (14)
```

An escaping deformation is

```text
w_s=r+s,              gamma_s=s,
P_s=gamma_s+w_s-2w_s^3,
Q_s=w_s^4-w_s^2+2P_sw_s.                              (15)
```

With target `(Q_s/c^2,P_s/c,c)`, reconstruction yields

```text
x_s=c/s,
y_s=r/(3c),
z_s=-s(3s^2-3s+4r)/(3c^2).                            (16)
```

The two colliding branches escape.  Transversely to `L=0`,

```text
w-r~+-sqrt(t/(6r^2-1)),                               (17)
```

for a local equation `t` of `L`.  Generic inertia is therefore one
transposition, of cycle type `2 1^2`.

Equations `(10)--(17)` prove the exact set equality

```text
S_G=V(C) union V(L)=V(C L).                            (18)
```

## 3. Raw and generic-primitive discriminant ledgers

Along `C=0`, the special primitive `w` remains bounded on the escaping sheet
and separately fails to distinguish two finite residual sheets at `w=0`.
Its raw discriminant orders are

```text
v_C(Disc_w f)=2,                 v_L(Disc_w f)=1.       (19)
```

The first satisfies `2=inertia_defect+2i=0+2*1`: it is the monogenic order's
index square from `(5)`, not a critical or inertia factor.  The second has
`i=0` and is the true transposition factor.

For the intrinsic infinity ledger take a Zariski-generic affine source
primitive

```text
T=alpha x+beta y+delta z,        alpha delta!=0.        (20)
```

Along `C=0`, one `T`-root has valuation `-2`, while three are finite and
residually distinct.  Hence

```text
v_C(Disc q_T)=2*3*(-2)=-12.                            (21)
```

The primitive clearing coefficient has order `rho_C=2`; for the integral
quartic `N_T`,

```text
v_C(Disc N_T)=6rho_C+v_C(Disc q_T)=0,
E_C=6rho_C-v_C(Disc N_T)=12.                           (22)
```

Thus `E_C=12` belongs to the **generic infinity-separating primitive** `(20)`,
not to the bounded raw primitive `w`.

Along `L=0`, two roots have valuation `-1/2` and two are finite.  Five root
pairs involve an escaping root, so

```text
v_L(Disc q_T)=2*5*(-1/2)=-5.                          (23)
```

Here `rho_L=1`, and

```text
v_L(Disc N_T)=6-5=1,
E_L=6rho_L-v_L(Disc N_T)=5.                            (24)
```

Thus the exact component ledger is

| component | generic lost sheets | infinity inertia | raw `w` order | generic-primitive `E` |
|---|---:|---|---:|---:|
| `C=0` | 1 | identity `1^4` | 2 (even index square) | 12 (even) |
| `L=0` | 2 | transposition `2 1^2` | 1 (odd) | 5 (odd) |

In both rows, primitive clearing parity equals the sign of infinity inertia.

## 4. Verdict and loss boundary

HYP-9027's Keller-restricted successor “every Jelonek component has odd
infinity inertia” is **REFUTED** by `C=0`.  The corrected theorem survives on
this first Keller test:

```text
cleared-discriminant parity = infinity-inertia sign.   (25)
```

The target-coordinate factor `C^2`, finite critical branching, and infinity
inertia are three different objects.  Saturating away `C` because the
`(P,Q,C)` chart degenerates would wrongly delete a genuine Jelonek component;
treating `C^2` as ramification would wrongly assign it a transposition.

This theorem does not classify other quartic Keller maps, identify all special
strata inside `V(CL)`, or settle `JC(2)`/`DC(2)`.  THM-3448 subsequently gives
a degree-five Keller `C_3` component and proves that no weighted quartic has
one; the arbitrary quartic/two-jet lane remains open.  The present theorem
describes exactly the boundary of the explicit weighted quartic `G`.

## 5. Exact companion and audit

The deterministic companion verifies `(2)--(5)`, irreducibility controls for
`L`, the reconstruction identities, both asymptotic deformations, the
finite-fibre control `(13)`, the valuation/clearing table, and hostiles
separating chart index from finite and infinity ramification.  Normal and
optimized runs reproduce the stored transcript byte-for-byte after LF
normalization:

```text
python -B 04-computation/jc_weighted_quartic_jelonek_inertia_thm3441.py
python -B -O 04-computation/jc_weighted_quartic_jelonek_inertia_thm3441.py
```

The script, transcript, semantic, and transcript-payload hashes are recorded
in the frontmatter.
