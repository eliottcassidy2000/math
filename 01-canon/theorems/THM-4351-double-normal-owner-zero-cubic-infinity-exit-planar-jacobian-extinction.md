---
id: THM-4351
title: "Double-normal-owner-zero cubic infinity-exit planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4344/4350 + VERIFIED-EXACT + TWICE
  INDEPENDENTLY HOSTILE-AUDITED. In the inherited reduced (2,3), exact-
  weight-twelve seam, Z=beta_11=zeta_3=W=xi_10=alpha_11=Theta=0 and
  U*K!=0 are extinct, with no condition on eta. All 24 exact supports have
  fan M,E00. Off eta^2=4KU the carrier is smooth genus two of form order six.
  On equality it has an exact A5 contact whose fixed nonzero cubic coefficient
  forces critical depths 1,2,3; their invariant tails have genera 2,1,1,
  two proved attachments, and positive primitive orders 34,16,10. U=0, K=0,
  seam entry, JC(2), and DC(2) are not claimed.
source: root + F00 collision/referee agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4350-first-normal-owner-cubic-infinity-exit-planar-jacobian-extinction
related:
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
  - THM-4352-even-self-similar-cusp-reciprocal-parity-and-attachment-law
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_thm4351.py
primary_output: 05-knowledge/results/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_thm4351.out
primary_script_sha256: b7058b43b9dc2dea1c3ad0807aa15a167f8927dab5f981ae3c2876a4f7aef847
primary_output_sha256: eb65b3915da1534fc05d8bf2486ceb16dd042b4d342e756ef61a9bcfa706bc59
referee_script: 04-computation/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_referee_thm4351.py
referee_output: 05-knowledge/results/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_referee_thm4351.out
referee_script_sha256: a7e25fa5bc7884a8292eb5787a0fa7c6466d2b75bcd0d42fe886efed40a1f0df
referee_output_sha256: 0b5e4757bbb16bdfe5490620071407f8d58baeaf08f589e304e9200b9ebf821c
hash_basis: raw LF bytes
audit: >
  PASS AFTER PRIMITIVE-BASE REPAIR. The import-free primary reconstructs the
  complete source, 24-support fan, exact faces and edges, generic carrier,
  doubled-edge factorization, A5 intersection length, reciprocal Morse
  chart, all realized critical depths, primitive and common base tails,
  both attachment charts, graph conservation, and form orders. THM-4350's
  independently written primary and global referees reproduce the fan,
  critical series, common-base tails, and incidence. The theorem records
  68/64/60 only as common-cover orders; the invariant primitive orders are
  34/16/10. A second import-free hostile referee independently repeats the
  support/fan, primitive chart, exact critical series, tail, attachment,
  surface, order, and graph ledgers in 293 exact checks. Both normal/optimized
  pairs byte-match their frozen outputs.
---

# THM-4351 -- Double-normal-owner-zero cubic infinity-exit planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4344/4350 + VERIFIED-EXACT + TWICE
INDEPENDENTLY HOSTILE-AUDITED. THE ENTIRE DISPLAYED `E00` CORNER IS EXTINCT.
`U=0`, `K=0`, SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OUTSIDE THIS THEOREM.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam, with the
complete sixteen-term source and

```text
e=-1376/135,                    K=2848/45-(7/6)Delta.       (1)
```

Write `u=upsilon_5` and `alpha=alpha_11`. Impose

```text
Z=beta_11=zeta_3=W=xi_10=alpha=Theta=0,
U*K!=0.                                                    (2)
```

Then no nonautomorphic planar Keller pair lies on `(1)--(2)`, for arbitrary
`eta,u,Phi,Delta` subject only to `(1)`. The conclusion is relative to the
inherited toroidal, good-target, and proper-flat interfaces. It neither
proves entry into this seam nor proves `JC(2)` or `DC(2)`.

The inheritance pass was:

- closest proved mechanism: THM-4350's exact four-owner atlas already
  identifies the omitted `M,E00` fan;
- canonical hostile: THM-4344's even `A5` contact, where the affine tail
  loses one genus unit unless both boundary attachments are retained;
- corrected near miss: common-cover orders are not primitive orders. Both
  remain positive here, but their numerical values must not be conflated;
- least-used sidecar: the fixed coefficient `e`, which terminates what first
  appears to be an unbounded critical-value hierarchy.

The Anchor / Niche / Wildcard board was

```text
E00 carrier | doubled outer edge | critical series | primitive ray
two marked ends | graph cycle | reciprocal depth | natural block.          (3)
```

## 2. Exact two-face fan and the generic carrier

On `(2)`, the seam has three exact support-status classes, represented by
`Delta=0,3968/63,1`; the presence bits of `u,Phi,eta` are independent. Thus
there are exactly

```text
3*2^3=24                                                  (4)
```

realizable active supports. Literal lower-hull enumeration gives

```text
24/24: M,E00.                                             (5)
```

Up to torus monomials, the exact faces are

```text
M   =(P-S^2)(UP^6-1),

E00 =S^2(1-UP^6-eta*S*P^4-K*S^2P^2).                     (6)
```

Their polygon ledger is

| object | vertices | `(2Area,boundary,interior)` |
|---|---|---:|
| `M` | `(0,1),(2,0),(2,6),(0,7)` | `(24,14,6)` |
| `E00` | `(2,0),(4,2),(2,6)` | `(12,10,2)` |
| global | `(0,1),(2,0),(4,2),(2,6),(0,7)` | `(36,12,13)` |

The outer schemes and the internal scheme are

```text
X-1; 1-KX^2; K+eta X+UX^2; X-1; U-X^6;
M--E00: 1-UX^6.                                          (7)
```

Under `U*K!=0`, all fixed schemes are reduced. The only variable edge
discriminant is

```text
L00=eta^2-4KU.                                           (8)
```

Put `Y00=2KSP+eta P^3`. Modulo the `E00` equation,

```text
Y00^2=(eta^2-4KU)P^6+4K.                                (9)
```

When `L00!=0`, the branch polynomial in `(9)` has degree six and is
squarefree. Hence the smooth projective `E00` normalization has genus two.
The `M` face has seven rational components and twelve internal nodes; its six
attachments to `E00` give

```text
(V,E,b1)=(8,18,11),                 11+2=13.             (10)
```

The global source-infinity packet is `(11,7,7,2,2,1)` and satisfies
`sum(e_i-1)=24=2*13-2`. This checks, but does not prove, completeness.

The `E00` supporting plane is `(1/3,1/6,-2/3)`. The inherited density formula
therefore gives

```text
base=6,        order=6*(5/6-1/3-1/6+2/3)=6>0.           (11)
```

Thus the generic genus-two carrier is Keller-constant; every other component
is rational.

## 3. Collision factorization and exact `A5` incidence

Suppose `L00=0`. Set

```text
a=-eta/(2K),              U=Ka^2,              a!=0,
lambda^2=K.                                                (12)
```

Then the three incidence identities are

```text
E00/S^2
 =[1-lambda(SP-aP^3)][1+lambda(SP-aP^3)],

K+eta X+UX^2=K(1-aX)^2,

1-UX^6=(1-lambda*a*X^3)(1+lambda*a*X^3).                (13)
```

The two cubics in the last line are reduced and disjoint. Hence the six
`M--E00` nodes split `3+3` between the two rational signs.

At reciprocal infinity put `x=P^-1` and `v=S/P^2`. The signs become

```text
x^3-lambda(v-a)=0,            x^3+lambda(v-a)=0.         (14)
```

Their intersection quotient is

```text
C[[x,v-a]]/(v-a,x^3),                                   (15)
```

of length three. Thus `(14)` is exactly a two-branch `A5` contact with
`delta=3`, not merely a doubled edge. Before resolving its normal series,

```text
(V,E,b1)=(9,18,10),                 10+delta(A5)=13.      (16)
```

## 4. Complete reciprocal chart and forced depth

Use the primitive `E00` face scale

```text
Q=sigma^6,       s=sigma^-2 S,       p=sigma^-1 P,
G=sigma^4 F_Q.                                             (17)
```

Define

```text
H6=UP^6+eta SP^4+KS^2P^2,       H5=uP^5+Phi SP^3.       (18)
```

Direct expansion of every surviving source row gives, without ellipsis,

```text
G=(S^2-sigma^3P)
 [1-H6-sigma H5-sigma^2Delta P^4-e sigma^3P^3
    -(8/3)sigma^4P^2+3sigma^5P]-sigma^6S^2/2.            (19)
```

Now put

```text
x=P^-1,                v=S/P^2,                rho=sigma*x.
```

Multiplication by `x^10` transforms `(19)` exactly into

```text
(v^2-rho^3)
 [x^6-A(v)-rho B(v)-rho^2Delta-e rho^3
      -(8/3)rho^4+3rho^5]-rho^6v^2/2=0,                 (20)

A(v)=U+eta v+Kv^2,                B(v)=u+Phi v.
```

At `(12)`, `A(v)=K(v-a)^2`. The prefactor in `(20)` is the unit `a^2` at
`(v,rho)=(a,0)`, and the Hessian is `2K`. Division and the parameter-
preserving formal Morse lemma therefore give one completed equation

```text
y^2=x^6-psi(rho),                  rho=sigma*x.          (21)
```

The critical-point approximation `v0=a-rho*Phi/(2K)` has derivative error
starting only in degree `rho^9`. Consequently

```text
psi(rho)
 =(u+Phi*a)rho
 +(Delta-Phi^2/(4K))rho^2
 +e rho^3+O(rho^4).                                      (22)
```

Because `e=-1376/135!=0`, the critical depth is forced to be

```text
r=ord_rho psi in {1,2,3}.                               (23)
```

All three depths occur. With `a=1`, hence `U=K,eta=-2K`, take respectively

```text
(Delta,Phi,u)=(0,0,1),       (1,0,0),       (0,0,0),    (24)
```

and determine `K` from `(1)`.

## 5. Tails and the primitive-base firewall

Write `psi=rho^r C(rho)`, with `C(0)=c0!=0`. A convenient common dominating
base is

```text
sigma=tau^[2(6-r)],        x=tau^[2r]X,
y=tau^[6r]Y.                                                (25)
```

After removing the horizontal square `X^[2 floor(r/2)]`, the reduced
exceptional function field is

```text
Y0^2=X^[r mod 2](X^(6-r)-c0).                            (26)
```

Every branch polynomial in `(26)` is squarefree. However, `(25)` is not the
minimal base. Put `g_r=gcd(r,6-r)`. The primitive weights are

```text
(ord_t sigma,ord_t x,ord_t y)
 =((6-r)/g_r, r/g_r, 3r/g_r).                           (27)
```

Their gcd is one, and direct substitution produces exactly the same
function field `(26)`, not a ramified root-coordinate cover. The complete
table is

| `r` | tail equation | genus | persistent delta | primitive weights | primitive surface | primitive form order | common-cover order |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `Y0^2=X(X^5-c0)` | 2 | 0 | `(5,1,3)` | `A5` | 34 | 68 |
| 2 | `Y0^2=X^4-c0` | 1 | 1 | `(2,1,3)` | `A2` | 16 | 64 |
| 3 | `Y0^2=X(X^3-c0)` | 1 | 1 | `(1,1,3)` | `A1` | 10 | 60 |

The common base is a further ramification of indices `2,4,6`, respectively.
Thus `68,64,60` are valid common-cover orders but are not the invariant
primitive numbers. This distinction is the MISTAKE-540 firewall.

## 6. Two attachments and exact genus conservation

On the common base set

```text
z=1/X=tau^12/rho,                 h=y/x^3.              (28)
```

The exact complementary chart is

```text
z*rho=tau^12,
h^2=1-z^(6-r)C(rho).                                     (29)
```

At `z=rho=tau=0`, the points `h=+1,-1` are distinct and etale. Normalizing
the Morse coordinate by `y=lambda(v-a)+O(rho)` shows that they land one on
each sign in `(13)`. Each sign already has three nodes into the connected
`M` subcurve. Adding the tail vertex and its two distinct attachment edges
therefore closes exactly one graph cycle:

```text
(9,18,10) -> (10,20,11).                                (30)
```

The two attachment surfaces are `A11` on the common base. On the primitive
bases they are `A5,A2,A1` for `r=1,2,3`. Resolving `z*rho=t^N` replaces one
contracted edge by `N` edges and `N-1` rational vertices, so resolving both
surfaces leaves `b1=11`.

For all three rows,

```text
floor(r/2)+g_tail+1=3=delta(A5),                         (31)

b1+g_tail+persistent_delta=13.                           (32)
```

Here `(10,20,11)` is the **contracted normalized component graph**. The
persistent horizontal delta is separate and must not also be counted as a
dual-graph edge.

## 7. Differential orders

Let `Phi_root=x^10G`. Since `P=x^-1` and `S=v/x^2`, exact differentiation
gives

```text
G_S=x^-8(Phi_root)_v,          -dP/G_S=x^6 dx/(Phi_root)_v.   (33)
```

Morse preparation makes `(Phi_root)_v` a unit times `y` on `(21)`. Combining
this with the face order six from `(11)` gives

```text
phi^*eta_0=unit*sigma^6*x^6 dx/y.                        (34)
```

On the common base `(25)`, its order is

```text
6*2(6-r)+6*2r+2r-6r=72-4r,                              (35)
```

namely `68,64,60`. On the primitive base `(27)`, the order is

```text
(36-2r)/gcd(r,6-r),                                     (36)
```

namely `34,16,10`. Every nonrational tail is therefore Keller-constant on
the invariant minimal chart.

## 8. Proper-flat extinction

If `L00!=0`, the only nonrational component is the smooth genus-two `E00`
carrier of positive order six. If `L00=0`, both primary signs are rational;
the only nonrational component is the genus-two or genus-one tail in `(26)`,
of positive primitive order. The `M` components and every surface-resolution
chain are rational. Equations `(7)` and `(29)` show there is no omitted
collision owner.

After finite base change, normalization, and regularization, every special
component map to the inherited good elliptic target is constant. With actual
positive component multiplicities on a common dominating base, the inherited
proper-flat identity gives, for a positive-degree target line bundle `L`,

```text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0.       (37)
```

This contradicts the positive generic response degree of a nonautomorphic
Keller pair and proves `(1)--(2)`. **QED.**

## 9. Reciprocal fundamental domain and natural-number compression

An abstract even `A5` split admits `r=1,...,5`, with reciprocity

```text
r <-> 6-r.                                               (38)
```

The fixed nonzero `e` in `(22)` restricts the actual corner to `r<=3`.
Therefore `(23)` contains exactly one representative from each reciprocal
orbit

```text
{1,5},             {2,4},             {3}.              (39)
```

This is not a chosen convention: the source coefficient selects a canonical
fundamental domain.

The general oriented even-cusp index at `g=3` is

```text
n=(g-1)^2+r=4+r=2T_2+1+(r-3).                           (40)
```

Using `T(-3)=T(2)`, its center is also `2T(-3)+1=7`.
Reciprocity negates the offset `r-3` and reflects the five-number block
`5,6,7,8,9`. The forced depths occupy the half-block `5,6,7`. On the quotient
`h=min(r,6-r)=r in {1,2,3}`; replacing the `h`th odd square by the natural
rank `h` is lossless only inside this fixed `A5` block.

The rank alone still forgets the coefficient `c0`, tail equation, primitive
base, and the two sign addresses in `(29)`. The lawful packet is

```text
(h; reciprocal fixed/oriented flag, c0, primitive ray,
 two marked attachment addresses).                                  (41)
```

This also explains why a tournament edge is too small. Reciprocity is an
involution with the fixed type `r=3`, while the tail has two distinct ends.
Forcing the fixed point into an orientation or forgetting either endpoint
deletes the graph unit in `(31)`.

## 10. Reproduction and remaining boundaries

Run from the repository root:

```bash
python3 -B 04-computation/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_thm4351.py
python3 -B -O 04-computation/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_thm4351.py
python3 -B 04-computation/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_referee_thm4351.py
python3 -B -O 04-computation/jc2_m12_double_normal_owner_zero_cubic_infinity_exit_extinction_referee_thm4351.py
```

Both normal/optimized pairs byte-match their frozen outputs. Each import-free
script reconstructs the source rather than importing THM-4350's fan data. The
THM-4350 primary and global-referee certificates independently contain the
same `E00` fan, critical series, common-base tails, and `3+3` incidence.

What is closed is exactly `(1)--(2)`. The next genuinely new fan problems are
the intersections `U=0` and `K=0`; seam entry and `JC(2)`/`DC(2)` remain
upstream global questions.
