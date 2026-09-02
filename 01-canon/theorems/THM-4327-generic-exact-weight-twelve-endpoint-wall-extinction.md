---
id: THM-4327
title: "Generic exact-weight-twelve endpoint-wall extinction"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4103/4230/4297 +
  VERIFIED-EXACT.  In the inherited reduced (2,3) seam, good-target
  differential order excludes (i) U=0, WZ nonzero when the two displayed
  outer-face discriminant implications hold, with Lambda arbitrary, and (ii) Z=0,
  U beta_11 zeta_3 nonzero, with W and all remaining lower coefficients
  arbitrary.  Every positive-genus face and every audited collision tail is
  Keller-constant; all other components are rational, so proper-flat degree
  conservation gives a contradiction.  Remaining endpoint subwalls include
  U=0 with WZ=0, the repeated (2,5) and deepest (2,3) cusps, and on Z=0
  with U nonzero the beta_11=0 and zeta_3=0 owner strata.
  Seam entry, JC(2), and DC(2) remain open.
source: root + u_zero_wall + z_zero_wall / planar-Jacobian endpoint-wall session, 2026-09-01
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
related:
  - THM-4232-weight-eleven-u-zero-primitive-cm-planar-jacobian-exclusion
  - THM-4248-weight-eleven-z-zero-owner-descent-planar-jacobian-exclusion
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
  - THM-4298-weighted-face-source-normal-unimodular-visibility-transform
  - THM-4299-d-zero-square-face-elliptic-splitting-and-off-corner-extinction
  - THM-4328-seam-covariant-student-stein-face-visibility
  - THM-4334-beta-zero-exact-weight-twelve-endpoint-wall-extinction
mistake_firewall:
  - MISTAKE-487
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through the audited
  general face/edge-model use in THM-4045.  The endpoint-specific primitive
  heights, face and edge regularity, local contacts, and component ledgers
  are checked below and in the exact certificates.
scripts:
  - 04-computation/jc2_m12_u0_endpoint_extinction_thm4327.py
  - 04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py
outputs:
  - 05-knowledge/results/jc2_m12_u0_endpoint_extinction_thm4327.out
  - 05-knowledge/results/jc2_m12_z0_endpoint_extinction_thm4327.out
script_sha256:
  - 389a6a7aae958a61aecb7b9340f6875973b2443409f9bb872555f5dd2378bb60
  - 6aa9087afd413833d3168b27efe3e65e779ab796a9b537ef1e624a6380fac551
output_sha256:
  - 94203d89e320ca237d733086e8d4755797aba98872eedb97be8d0ea32d2d431c
  - dfc9b586362d04ff6ef4a1e1cbb9c7d78e33ed27807106736a65e4bd90bd2c80
hash_basis: raw LF bytes
audit: >
  PASS.  Two standard-library Fraction certificates independently own the
  two endpoint branches.  The U-wall path checks 7,500 simultaneous
  row/collision hostiles, five owner complexes, an additional 1,024-mask
  owner audit, polygons, genera, edges, packets, differential orders, the
  Lambda-zero A23 specialization/order table, and the sharp hostiles.  The
  Z-wall path checks 24,576 theorem hostiles and a 131,072-mask over-atlas,
  both Lambda-zero identities, symbolic normal forms, generic-point
  derivatives, edge gates, component graphs, polygons, genera, packets,
  and exact differential orders.  Normal and optimized streams byte-match the frozen
  transcripts; the U-wall audit mode is also optimization-stable.
---

# THM-4327 -- Generic exact-weight-twelve endpoint-wall extinction

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4103/4230/4297 +
VERIFIED-EXACT.  THIS CLOSES THE TWO DISPLAYED OPEN SUBSTRATA, NOT THE WHOLE
`U=0` OR `Z=0` WALL.  SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement, notation, and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam.  Retain the
complete residual source through weight twelve from THM-4230 and put

```text
U=[p^6]H,          W=[p^3y^2]H,       Z=[y^4]H,
u=[p^5]H=upsilon_5,                  a=[p^4y]H=alpha_11,
beta=[py^3]H=beta_11,                zeta=[y^3]H=zeta_3,
d=[p^4]H=Delta,     eta=[p^3y]H,      e=[p^3]H=-1376/135,
D=W^2-4UZ,          Lambda=U+W+Z,
E5=a^2-4Wu,         E3=eta^2-4eW.                       (1)
```

No lower coefficient omitted from THM-4230's source is silently set to
zero.

> **Theorem U (generic `U=0` endpoint wall).**  No nonautomorphic planar
> Keller pair exists on
>
> ```text
> U=0,                 WZ!=0,
> u!=0  implies E5!=0,
> u=a=d=0 implies E3!=0.                                  (2)
> ```
>
> Subject only to the displayed conditions, `Lambda=W+Z` and all remaining
> lower coefficients are arbitrary.

> **Theorem Z (generic `Z=0` endpoint wall).**  No nonautomorphic planar
> Keller pair exists on
>
> ```text
> Z=0,                 U*beta*zeta!=0.                           (3)
> ```
>
> Subject only to the displayed conditions, `W` and all remaining lower
> coefficients are arbitrary.

The inheritance pass was:

- closest proved mechanism: THM-4297/4299's good-differential extinction
  followed by proper-flat degree conservation;
- canonical hostiles: THM-4218's hidden elliptic face, MISTAKE-531's double
  top-edge root, and MISTAKE-540's genus inflation on a ramified toric cover;
- corrected near miss: the weight-eleven primitive-CM exclusions
  THM-4232/4248 do not transfer--the `U=0` main component has an explicit
  `j=0` quotient;
- least-used sidecar: the exact supporting plane of every face, which gives
  the divisorial order of the good target differential before any Hom
  decomposition.

The live concept board was

```text
{endpoint owner, lower hull, invariant face field, good differential,
 repeated contact, graph genus, proper-flat degree}.                   (4)
```

## 2. The common good-differential formula

Write the exact source equation as

```text
F_Q=(s^2-p)(1-QH)-Q s^2/2.                              (5)
```

Let a lower face have supporting plane

```text
h(r,l)=a_s r+b l+c,                                     (6)
```

and choose `Q=sigma^L`, with `L` clearing the source face and the good
elliptic target.  Put

```text
s=sigma^(-L a_s)S,             p=sigma^(-L b)P,
F_Q=sigma^(Lc)G.                                        (7)
```

THM-4103's Keller residue identity and the chain rule give

```text
Q ds/(F_Q)_p
 =sigma^[L(1-a_s-b-c)] dS/G_P.                          (8)
```

The original target differential is `sigma^(L/6)` times the invariant
differential `eta_0` on the good elliptic model.  Therefore

```text
ord_face(phi^*eta_0)=L(5/6-a_s-b-c).                    (9)
```

At the generic point of every positive-genus face below, the certificates
and displayed equations show that `G_P` is nonzero in the face function
field.  Thus `dS/G_P` is a generic relative-differential generator.  If
`(9)` is positive, the specialized component pulls `eta_0` back to zero and
its map to the characteristic-zero elliptic curve is constant.  Every
rational face, toric subdivision, and point-blowup component is also
constant.  The proper-flat interface of THM-4045/4103/4230 then conserves
the positive generic map degree, contradicting a special fibre on which
every component has degree zero.

This mechanism is used below only after a complete face, contact, and genus
inventory.  A positive order on a list of primary faces alone is not treated
as a theorem about an unaudited degenerate contact.

### 2.1 Endpoint regular-model applicability

The face/edge theorem imported and audited in THM-4045 applies to the new
complexes, for the following endpoint-specific reasons.

First, after the bases used below, every supporting plane becomes an
integral affine height graph.  The coefficient of the height coordinate in
its three-dimensional normal is `-1`, so the normal is primitive and every
face has multiplicity one.  The complete scaled slope triples
`L(a_s,b,c)` are

```text
U wall A--C, L=60: M=(5,10,-10), H5=(0,12,-12),
                    V4a=(-15,15,-15), V3a=(-40,20,-20);
U wall D, L=24:    M=(2,4,-4), V4W=(-3,6,-6);
U wall E, L=12:    M=(1,2,-2), V3W=(-4,4,-4);
Z wall W!=0,L=84: M=(7,14,-14), V7=(12,12,-24),
                    T3=(28,0,-56);
Z wall W=0,L=36:  M=(3,6,-6), V9=(4,6,-8),
                    T3=(12,0,-24).                      (9a)
```

All these bases are divisible by six, so the good target is integral too.

Second, every irreducible face branch is torus-smooth under the theorem
gates.  On the `U` wall, the two nonconstant exponent vectors of

```text
C=1-WS^2P^5-ZS^4P^4
```

have determinant `-12`; the two logarithmic derivatives cannot vanish
simultaneously.  The `H5` face has `PF_P=5` on the face, the `V4a,V3a`
faces are linear in `S`, and `V4W` has nonzero `S`-derivative on the torus.
For `V3W`, the monomial change

```text
(S,P)->(V=SP,X=P^(-1))                                  (9b)
```

has determinant `-1`, and `(15)` is a smooth cubic exactly when `E3!=0`.
On the `Z` wall, the relevant pairs of nonconstant exponent vectors have
determinants

```text
-12 for the W!=0 main branch, -7 for V7,
-18 for the W=0 genus-six branch.                       (9c)
```

Hence their logarithmic gradients have no common torus component; the `T3`
face has nonzero `S`-derivative on its equation.  These same calculations
show that `G_P` is a nonzero residue-field element at every positive-genus
generic point.

Third, the complete symbolic edge schemes are listed in Sections 4, 7, and
8.  Their derivatives and resultants are nonzero under precisely the active
coefficient and discriminant gates.  All internal edges are primitive or
have the displayed squarefree polynomial; the general model therefore adds
only rational toric chains.  The sole nonsimple schemes retained by a
theorem are the `Lambda=0` double roots handled in Section 6.

The reducible main faces are also explicit.  When `Lambda!=0`, each
`R:S^2=P` intersection with the other main branch is transverse.  In local
branch parameters the exact scaled equation is

```text
u_branch v_branch=sigma^L*(unit),                       (9d)
```

so its resolution is an `A_(L-1)` rational chain.  On `Z=W=0`, the six
factors of `1-UP^6` are distinct and the same calculation applies to all
twelve `R--line` nodes.  The dual-graph and Pick ledgers below account for
the entire generic genus, ruling out an unlisted positive-genus component.

Finally, the torus model is the inherited generic source, not a quotient:
`t=p-s^2` and `x=s/t` invert `(x,t)->(s=xt,p=t+s^2)`.  On `t=0` the torus
equation `(5)` becomes `-Q s^2/2`, so `t` cannot vanish there.  Thus the
torus curve is a dense open of the smooth generic completion.  Dokchitser's
face/edge construction, the explicit node resolutions `(9d)`, and the
Section 6 contact audit supply a proper flat model whose complete component
inventory is exactly the one used below.  This is the missing applicability
lemma; the order table alone would not suffice.

## 3. The `U=0` lower hull is an exact five-case atlas

Under `(2)`, every coefficient specialization lies in exactly one row of
the following table:

| case | coefficient gate | lower faces besides `M` |
|---|---|---|
| `A` | `u!=0`, `E5!=0` (`a` arbitrary) | `H5` |
| `B` | `u=0`, `a*d!=0` | `H5,V4a` |
| `C` | `u=d=0`, `a!=0` | `H5,V3a` |
| `D` | `u=a=0`, `d!=0` | `V4W` |
| `E` | `u=a=d=0`, `E3!=0` | `V3W` |

The exact planes `(a_s,b,c)` are

```text
M   =( 1/12, 1/6,-1/6),       H5 =(   0,1/5,-1/5),
V4a =(-1/4, 1/4,-1/4),        V3a=(-2/3,1/3,-1/3),
V4W =(-1/8, 1/4,-1/4),        V3W=(-1/3,1/3,-1/3).     (10)
```

The finite universe consists of all fifteen allowed labelled rows through
weight twelve after deleting `U p^6`.  A row `p^i y^j` contributes

```text
(j+2,i+j,1),                 (j,i+j+1,1),               (11)
```

together with the fixed points `(2,0,0),(0,1,0),(2,0,1)`.  The certificate
independently toggles every optional lower row and independently deletes
each active aggregate point among

```text
(2,3,1),(2,4,1),(2,5,1),(3,4,1),(3,5,1),(4,5,1).      (12)
```

These deletions deliberately overcount coefficient cancellations.  The
exact hostile populations in cases `A--E` are

```text
4500, 1080, 720, 720, 480,                 total 7500,  (13)
```

and every trial returns precisely `(10)`.  A second owner-mask audit crosses
all `1,024` coefficient states and recovers the same five complexes.

## 4. `U=0` faces, genera, and sharp edge boundary

Up to torus monomials, the face cores are

```text
g_M=(S^2-P)C,
C=1-W S^2P^5-Z S^4P^4,                         genus 4,

H5=-1+P^5(u+aS+W S^2),                         genus 2,

V4a=-1+dP^4+aSP^5,                             genus 0,
V3a=-1+eP^3+aSP^5,                             genus 0,

V4W=-1+dP^4+W S^2P^5,                         genus 2,
V3W=-1+eP^3+eta SP^4+W S^2P^5.                genus 1.  (14)
```

The positive-genus normal forms include

```text
H5:   Y=(2WS+a)P^3,          Y^2=E5 P^6+4WP,
V4W:  Y=SP^3,                W Y^2=P-dP^5,
V3W:  X=1/P, V=SP,            (2WV+eta)^2=4WX^3+E3.     (15)
```

Thus `E5!=0` is exactly the extra reduced-face gate in case `A`; in cases
`B,C`, it equals `a^2!=0`.  The main face has polygon
`conv((0,0),(2,5),(4,4))`, with Pick ledger `(12,6,4)`.  The remaining Pick
ledgers reproduce the genera in `(14)`.

The outer edge schemes have the common prefix

```text
X-1,             1-ZX^4,             (X-1)(WX+Z),       (16)
```

followed respectively by

```text
A: W+aX+uX^2, u-X^5;
B: W+aX, a+dX, d-X^4;
C: W+aX, a+eX, e-X^3;
D: W+dX, d-X^4;
E: e+eta X+WX^2, e-X^3.                                (17)
```

The internal schemes are `1-WX`, and additionally `1-aX` in `B,C`.
Exact discriminants and resultants show that the case gates, `E3!=0` in
case `E`, and `WZ!=0`
make every active root simple, except the intentionally retained
`Lambda=0` collision in `(16)`.  The excluded gate `E5=0,u!=0` has local
model

```text
W(S-S0)^2=z^5,                                         (18)
```

a repeated `(2,5)` cusp with a possible genus-two tail.  It is not covered
by this theorem.

The other excluded discriminant `E3=0` in case `E` makes `(15)` the
cuspidal cubic `Y^2=4WX^3`.  Its normalization is rational, but a smoothing
can create a genus-one tail; that repeated `(2,3)` face also requires its
own differential audit.

For `Lambda!=0`, `R:S^2=P` and `C` meet in twelve transverse points.
Together with all adjacent primitive intersections, the dual graph has
first Betti number eleven.  The complete global ledgers are

| case | `(2Area,boundary,genus)` | positive face genera | graph `b1` |
|---|---:|---:|---:|
| `A` | `(46,14,17)` | `4+2` | `11` |
| `B` | `(45,13,17)` | `4+2` | `11` |
| `C` | `(44,12,17)` | `4+2` | `11` |
| `D` | `(44,12,17)` | `4+2` | `11` |
| `E` | `(42,12,16)` | `4+1` | `11` |

Hence no positive-genus component is missing from `(14)`.

## 5. All `U=0` face orders are positive

The minimal clearing bases in cases `A--E` are

```text
L=60,60,60,24,12.                                       (19)
```

Every plane in `(10)` has `b+c=0`.  Formula `(9)` gives

```text
M:       45,45,45,18,9;
H5:      50;
V4a:     65;                 V3a: 90;
V4W:     23;                 V3W: 14.                   (20)
```

All are strictly positive.  At a generic face point, `G_P` is nonzero:

```text
H5:       P F_P=5;
V_k:      P F_P=5+(k-5)c_kP^k, not identically zero;
V3W:      P F_P=5-2eP^3-eta SP^4, equal to 3 at S=0 on the face;
C:        G_P=(S^2-P)C_P, with both factors nonzero.    (21)
```

Consequently all curves in `(14)` are Keller-constant, and the rational
components are constant as well.

The old weight-eleven primitive-CM mechanism cannot replace `(20)`.  The
main genus-four component always has the explicit `j=0` quotient

```text
T=SP,             X=P^2,             Y=WP^3+2ZT^2,
Y^2=W^2X^3+4Z.                                             (22)
```

In case `E`, `V3W` itself is a smooth `j=0` elliptic curve under `E3!=0`.
Positive differential order, rather
than absence of elliptic Hom, is the decisive obstruction.

## 6. Both endpoint `Lambda=0` contacts are extinct

MISTAKE-531 forbids importing a simple-root response packet at this wall,
so no such import is used.  Put `Z=-U-W`, set

```text
A=2U+W,
```

and use the top-infinity variables

```text
b=1/S,                         r=P/S^2.                  (23)
```

The general main closure and its derivative along `R:r=1` are

```text
Ctilde=b^12-Ur^6-Wr^5-Zr^4,
partial_r Ctilde(0,1)=-(6U+5W+4Z)=-A.                  (24)
```

The exact Lambda-zero identity is

```text
Ur^6+Wr^5+Zr^4
 =(A/2)(r^6-r^4)-(W/2)r^4(r-1)^2.                     (25)
```

On Theorem U, `U=0,Z=-W`, so `A=W!=0` and the left side of `(25)` is
`Wr^4(r-1)`.  On Theorem Z, `Z=0,W=-U`, so `A=U!=0` and it is
`Ur^5(r-1)`.  In both cases `(24)` gives one length-twelve `A23` contact,
and `(25)` puts the germ in THM-4297's repeated-contact chart with
`U_eff=A/2!=0`.  THM-4297's proved four-step critical ladder and
terminal `b^12q` splitter therefore give the five possible tail orders

```text
6s+2b, (5s+9b)/2, 2s+4b, (3s+7b)/2, s+3b,             (27)
```

for positive divisorial valuations `s,b`.  Each is positive.  THM-4297 also
proves that the correction term in `(25)` first enters after the complete
critical ladder; the endpoint certificates check both specialized contacts
and all five resulting inequalities.  The A23 resolutions and all
possible positive-genus tails are therefore constant.  Replacing twelve
nodes by the length-twelve contact preserves the arithmetic-genus ledgers in
Sections 4 and 7, completing both theorems for arbitrary `Lambda`.

There is no base mismatch in this import: every endpoint base `L` is a
multiple of twelve.  If THM-4297 uses `Q=sigma^12` and the present chart uses
`Q=Sigma^L`, put `sigma=Sigma^(L/12)`; all local orders are multiplied by
the positive ramification index `L/12`.

For reference only, when `Lambda!=0` the audited source-infinity packets
and full/finite response totals are

```text
A--C: (11,11,5,5,2,2,2,2,1),             41/33;
D:    (11,11,9,2,2,2,2,1),               40/32;
E:    (11,11,4,4,2,2,2,2,1),             39/31.        (28)
```

They are not transported to `Lambda=0` and are not needed for extinction.

## 7. The `Z=0`, `W!=0` face complex

Assume `(3)` and first take `W!=0`.  Every one of the `16,384` simultaneous
row/collision hostiles has exactly the three planes

```text
M =(1/12,1/6,-1/6),
V7=(1/7,1/7,-2/7),
T3=(1/3,0,-2/3).                                        (29)
```

Their face cores are

```text
(S^2-P)(1-UP^6-WS^2P^5),
S^2(1-WS^2P^5-beta S^3P^4),
S^2(1-zeta S^3P^3-beta S^3P^4).                        (30)
```

The first two positive-genus normalizations have genus three; explicit
function-field forms are

```text
y^2=W x(x^6-U),
P^7=(beta^2/W^3) x^3/(1-x)^2.                          (31)
```

The `T3` face is rational.  Exact torus Jacobians and `G_P` resultants show
that the two genus-three components are reduced and that `G_P` is nonzero at
their generic points.

With `L=84`, formula `(9)` gives the exact orders

```text
63, 70, 98.                                             (32)
```

The global polygon and its Pick ledger are

```text
(0,1),(2,0),(5,3),(5,4),(4,5),(0,7),
(2Area,boundary,genus)=(46,14,17).                      (33)
```

When `Lambda!=0`, the face graph has vertices `R,C,V7,T3`; its twelve `R--C` nodes and the
primitive `C--V7,V7--T3` intersections give `E=14,V=4`, hence `b1=11`.
For `Lambda=0`, Section 6 replaces those twelve nodes by the length-twelve
contact and its audited tails without changing arithmetic genus.  The two
face genera sum to six, so no positive-genus component is missing.

The nontrivial edge schemes reduce to

```text
X-1, 1-zeta X^3, -beta X-zeta, -W X-beta,
(X-1)(UX+W), U-X^6, 1-WX, 1-beta X.                    (34)
```

All schemes except `(X-1)(UX+W)` are squarefree under
`U*W*beta*zeta!=0`.  The remaining scheme is squarefree when `Lambda!=0`
and is the Section 6 double root when `Lambda=0`.  In the simple-root case,
the corresponding source-infinity packet and response controls are

```text
(11,11,7,4,2,2,2,1),                    finite/full 34/40.    (35)
```

They saturate the genus ledger but are not needed for the constant-component
contradiction.

## 8. The `Z=W=0` replacement complex

When `W=0`, all `8,192` hostile masks instead have planes

```text
M =(1/12,1/6,-1/6),
V9=(1/9,1/6,-2/9),
T3=(1/3,0,-2/3).                                        (36)
```

The main core `(S^2-P)(1-UP^6)` splits into rational factors.  The sole
positive-genus normalization is

```text
1-UP^6-beta S^3P^4=0,                         genus 6,  (37)
```

while the `T3` face remains rational.  With `L=36`, their good-differential
orders are

```text
27, 28, 42.                                             (38)
```

The outer and internal edge schemes are

```text
X-1, 1-zeta X^3, -beta X-zeta, -U X-beta,
U(X-1), U-X^6, 1-UP^6, 1-beta X.                       (38a)
```

They are squarefree under `(3)`.

The global polygon and complete boundary controls are

```text
(0,1),(2,0),(5,3),(5,4),(2,6),(0,7),
(2Area,boundary,genus)=(45,13,17),
packet=(17,11,4,2,2,2,1),              finite/full responses 33/39.   (39)
```

Here the graph has `R`, the six rational factors of `1-UP^6`, the
genus-six face, and `T3`: twelve `R--line`, six `line--V9`, and one
`V9--T3` intersections give `E=19,V=9`, hence `b1=11`.  The genus-six face
plus that graph accounts for all genus seventeen.
Equations `(37)--(38)` make the only positive-genus component constant; all
others are rational.  Sections 7--8 prove Theorem Z.

## 9. What the larger `Z=0` atlas says--and does not say

As a signal for the remaining wall, the Z-certificate also drops the
`beta` and `zeta` gates while retaining `Z=0,U!=0`.  Its hostile
over-atlas has

```text
131072 masks,             51 lower-face complexes,
32 supporting planes,     min(5/6-a_s-b-c)=3/4>0.       (40)
```

Thus every primary face in that larger atlas still has positive formal
good-differential order.  This is **FINITE-EXACT evidence, not an exclusion
theorem**: at `beta=0` or `zeta=0`, the face complex changes and the
resulting contacts have not all been inventoried.  By contrast,

```text
U+W=0  =>  (X-1)(UX+W)=U(X-1)^2,                       (41)
```

is covered by Section 6.  MISTAKE-531 still forbids transporting the
simple-root response packet across it.

The precise remaining endpoint inventory after this theorem is:

```text
U=0:  WZ=0; or u!=0 with E5=0; or u=a=d=0 with E3=0;
Z=0:  U=0; or beta=0; or zeta=0.                        (42)
```

Intersections among these walls may refine the list further.  Neither
Theorem U nor Theorem Z asserts entry into the exact-weight-twelve seam.
THM-4334 subsequently closes the part `beta=0,U*W*zeta!=0`, including its
`Lambda=0` contact.  The live residual inventory is maintained in
`00-navigation/CURRENT-FRONTIER.md`.

## 10. The source-normal/Stein gate is strictly weaker

THM-4328 gives exact points on each extinct endpoint stratum which pass the
THM-4308 projected row-eight gate and THM-4315 row-nine bracket equation.
For example, its `Z=0` control has

```text
beta=1, zeta=-3/2,
U=-5200771070534/1042743375,
W=10155617023591/463441500,
U+W=70597468930183/4170973500 !=0.                     (43)
```

It lies inside Theorem Z and is therefore geometrically impossible despite
passing those finite algebraic gates.  This is the sharp separation:

```text
finite source-normal bracket data
        do not determine endpoint Newton components or their map degrees. (44)
```

The new connection has source a supporting plane, target the divisorial
order `(9)`, map the face scaling `(7)`, preserved predicate positivity of
the pulled-back good form, destroyed information lower-row labels and
contact gluing, required sidecar the full face/contact/genus inventory, and
cheapest decisive test the repeated-root hostiles `(18)` and `(41)`.

## 11. Reproduction and honest scope

Run from the repository root:

```bash
python3 -B 04-computation/jc2_m12_u0_endpoint_extinction_thm4327.py
python3 -B -O 04-computation/jc2_m12_u0_endpoint_extinction_thm4327.py
python3 -B 04-computation/jc2_m12_u0_endpoint_extinction_thm4327.py --audit
python3 -B -O 04-computation/jc2_m12_u0_endpoint_extinction_thm4327.py --audit
python3 -B 04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py
python3 -B -O 04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py
```

The primary normal and optimized streams must byte-match the frozen outputs;
both U audit modes must also terminate with `RESULT=PASS`.

What is proved, relative to the inherited toroidal/proper-flat interfaces,
is exactly `(2)` and `(3)`.  The response packets are consistency controls,
not transported through repeated roots.  The larger atlas `(40)`, seam
entry, the residual subwalls `(42)`, `JC(2)`, and `DC(2)` remain open.
