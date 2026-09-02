---
id: THM-4337
title: "Zeta-zero exact-weight-twelve endpoint-wall extinction"
status: >
  PROVED RELATIVE TO THM-4292/4297/4327 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. In the inherited reduced (2,3), exact-weight-twelve seam, the
  endpoint stratum Z=zeta_3=0 and U*beta_11!=0 is extinct with W,
  K=[y^2]H, every other lower coefficient, and Lambda=U+W arbitrary. For
  W!=0 the positive-genus faces have genera 3 and 3; for W=0 the remaining
  positive-genus face has genus 6. In both cases graph b1=11 completes genus
  17, all good-form orders are positive, and the Lambda=0 A23 contact splits
  before every t^6 correction. Exact-seam entry, JC(2), and DC(2) remain open.
source: root + jc zeta-zero agent / planar-Jacobian continuation session, 2026-09-02
depends_on:
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
related:
  - THM-4292-lambda-zero-repeated-face-keller-extinction
  - THM-4297-general-lambda-zero-central-and-tail-keller-extinction
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4334-beta-zero-exact-weight-twelve-endpoint-wall-extinction
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
mistake_firewall:
  - MISTAKE-487
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_z0_zeta0_endpoint_extinction_thm4337.py
primary_output: 05-knowledge/results/jc2_m12_z0_zeta0_endpoint_extinction_thm4337.out
primary_script_sha256: 8d43b5cbd789ff497b7fdd15d98c47e20b8ff9747cd748edd15245ee3427cec2
primary_output_sha256: 04bfaf696913cfa166e07e19d4c8d45a22a58b8a844892fcda96895d5451934e
independent_audit_script: 04-computation/jc2_m12_z0_zeta0_endpoint_extinction_independent_audit_thm4337.py
independent_audit_output: 05-knowledge/results/jc2_m12_z0_zeta0_endpoint_extinction_independent_audit_thm4337.out
independent_audit_script_sha256: 0a0a31fb8b69dcee245f56d8411d76db466a59eb51e285cff0d4adacc50a7a5d
independent_audit_output_sha256: 96a8502ea82ba69798306670839a7fbed8a519bdc2bfdcddc2e922243ebe080f
hash_basis: raw LF bytes
audit: >
  PASS. The primary path exhausts conservative 16,384-state and 4,096-state
  support over-atlases and checks all four W/K complexes, invariant face
  fields, Pick and graph genera, edges, positive orders, and the Lambda-zero
  contact. An import-free clean-room path generates all 293 affine candidates
  through the 28-point master support, tests independently cancellable
  multiply-owned points, and reconstructs the polygons, Kummer genera, edge
  coefficient sequences, Euler remainders, graphs, orders, and local splitter.
  A separate hostile referee audited the owner/deletion distinction, toric
  orbits, K-row timing, local/global base conversion, and proper-flat bridge.
  Normal and optimized streams byte-match both frozen outputs.
---

# THM-4337 -- Zeta-zero exact-weight-twelve endpoint-wall extinction

**PROVED RELATIVE TO THM-4292/4297/4327 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE `Z=zeta_3=0`, `U*beta_11!=0` ENDPOINT STRATUM IS EXTINCT FOR
ARBITRARY `W`, `K`, AND `Lambda`; EXACT-SEAM ENTRY, `JC(2)`, AND `DC(2)`
REMAIN OPEN.**

## 1. Statement and inheritance

Work over an algebraically closed field of characteristic zero in the reduced
`(2,3)`, exact-weight-twelve source family inherited from THM-4230 and
THM-4327.  Put

```text
Z=0,                 zeta_3=0,                 U*beta_11!=0,          (1)
```

and retain every other allowed lower coefficient.  In particular

```text
W and Lambda=U+W are arbitrary,              K=[y^2]H is arbitrary. (2)
```

Below `beta` abbreviates the nonzero coefficient `beta_11`.

The theorem is:

> No nonautomorphic planar Keller pair lies on `(1)`.  Every component in
> each of the four `W`/`K` special-fibre strata below maps constantly to the
> good elliptic target, contradicting proper-flat conservation of positive
> generic degree.

This is relative to the inherited exact seam.  It proves neither exact-seam
entry nor `JC(2)` nor `DC(2)`.

The inheritance pass is:

- closest proved mechanism: THM-4327's exact `Z=0` face/edge model and
  good-differential extinction;
- canonical hostile: MISTAKE-531's `Lambda=0` double top root, which requires
  a fresh local tail calculation rather than a transported response packet;
- corrected near miss: MISTAKE-540--all cyclic genera below are computed in
  their actual invariant function fields;
- least-used sidecar: the optional weight-six owner `K`, which inserts a
  rational face and an attachment without changing graph genus.

The live board was

```text
{zeta owner, K face, W split, invariant face field, edge orbit,
 A23 splitter, good differential, graph genus}.                           (3)
```

## 2. An exact four-stratum lower atlas

Use the complete sixteen-row universe and lift `(i,j)` to

```text
(j+2,i+j,1),                       (j,i+j+1,1),              (4)
```

as in THM-4327, together with the fixed hull points
`(2,0,0),(0,1,0),(2,0,1)`.  In particular, the residual point `(2,0,1)`
is nonoptional.  Rows `Z=(0,4)` and `zeta_3=(0,3)` are absent; rows
`p,p^2,p^3,U=(6,0),beta_11=(1,3)` are required.  First suppose `W!=0` and
require its row `(3,2)` as well.  Eight rows remain optional.  Independently
delete the six inherited collision/deletion points

```text
(2,3,1),(2,4,1),(2,5,1),(2,6,1),(3,4,1),(3,5,1).           (5)
```

This is a conservative over-atlas: on `zeta_3=0`, for example, `(3,4,1)` is
owned only by `eta`, and other listed points can also have only one active
owner after optional-row deletion.  Thus its artificial deletion bits need
not describe realizable coefficient cells; they enlarge the support universe
being tested.  All `2^8 2^6=16,384` hostile states nevertheless have exactly
one of

```text
K=0:  (M,V7),                     K!=0: (M,V7,Kf).           (6)
```

When `W=0`, the point `(2,6,1)` is uniquely owned by `U`.  We use the smaller
over-atlas that does not artificially delete either this point or the
eta-only point `(3,4,1)`, and delete only

```text
(2,3,1),(2,4,1),(2,5,1),(3,5,1)                            (7)
```

gives `2^8 2^4=4,096` states, every one with exactly

```text
K=0:  (M,V9),                     K!=0: (M,V9,Kf).           (8)
```

The supporting planes are

```text
M =(1/12,1/6,-1/6),             V7=(1/7,1/7,-2/7),
V9=(1/9,1/6,-2/9),              Kf=(1,-1/2,-2).             (9)
```

The equality triangle of `Kf` is

```text
(2,0,0),(4,2,1),(5,4,1).                                  (10)
```

Its middle point is uniquely owned by `K`, proving that the extra face is
present exactly when `K!=0`.  The independent certificate reconstructs the
row universe without importing THM-4327's hull engine and checks the four
complexes separately.

For `W!=0`, use global base `Q=Sigma^84`; for `W=0`, use `Q=Sigma^36`.
The integral slope triples are

```text
L=84:  M=(7,14,-14), V7=(12,12,-24), Kf=(84,-42,-168),
L=36:  M=(3, 6, -6), V9=( 4, 6, -8), Kf=(36,-18, -72).      (11)
```

Every height normal has coefficient `-1`, so all face multiplicities are
one.  Both bases are divisible by six and make the good elliptic model
integral.

## 3. Face fields, genera, and generic derivatives

### 3.1. The `W!=0` strata

Up to torus monomials, the positive-genus face cores are the two already
visible in THM-4327:

```text
G_M=(S^2-P)C,                  C=1-UP^6-WS^2P^5,
G_7=S^2B7,                     B7=1-WS^2P^5-beta S^3P^4.   (12)
```

Write `R:S^2=P`.  For `C`, the birational coordinates

```text
x=P^(-1),                 y=WS/P,
P=x^(-1),                 S=y/(Wx)                          (13)
```

give

```text
y^2=W x(x^6-U),                                           (14)
```

a squarefree degree-seven hyperelliptic model of genus three.  For `B7`, put

```text
x=WS^2P^5,                   1-x=beta S^3P^4.
```

Then

```text
P^7=(beta^2/W^3)x^3/(1-x)^2,
S=W(1-x)P/(beta x).                                        (15)
```

Thus `(15)` is the actual face field.  Its Kummer valuations `(3,-2,-1)`
are prime to seven, so the cover is connected and

```text
2g(B7)-2=-14+3(7-1)=4,                    g(B7)=3.          (16)
```

When `K!=0`, the additional face is

```text
G_K=S^2 A,                     A=1-KS^2P^2-beta S^3P^4.    (17)
```

The unimodular coordinate `V=SP` gives

```text
A=1-KV^2-beta V^3P,
P=(1-KV^2)/(beta V^3),                 S=V/P.              (18)
```

Hence `Kf` is irreducible and rational.  It contributes an attachment but no
positive genus.

The torus exponent determinants for `C,B7,Kf` are respectively

```text
-12,                         -7,                         2, (19)
```

so their logarithmic gradients have no common torus zero.  Their exact
Euler remainders are

```text
P C_P-5C=-(UP^6+5),
P (B7)_P-4B7=-(WS^2P^5+4),
P A_P-4A=2KS^2P^2-4.                                    (20)
```

The first two are not zero in their respective function fields by the
birational models `(14)--(15)`; the third is not zero because `V` is the
nonconstant parameter in `(18)`.  Since

```text
(G_M)_P|C=(S^2-P)C_P,        (G_7)_P|B7=S^2(B7)_P,
(G_K)_P|Kf=S^2A_P,                                           (21)
```

with every displayed prefactor nonzero at the relevant generic point, the
full `G_P` is nonzero on every face component used below.

### 3.2. The `W=0` strata

Now

```text
G_M=(S^2-P)(1-UP^6),
G_9=S^2B9,                    B9=1-UP^6-beta S^3P^4.       (22)
```

The main face is `R` plus the six rational factors of `1-UP^6`.  The actual
`B9` field is

```text
S^3=(1-UP^6)/(beta P^4).                                  (23)
```

It is a connected cyclic-three cover of the `P`-line.  The six simple finite
zeros, the order-four pole at zero, and order two at infinity are all
nontrivial modulo three, so Riemann--Hurwitz gives genus six.  The optional
`Kf` remains the rational curve `(18)`.  Its edge coefficients never require
`W`.

Here the torus determinants for `B9,Kf` are `-18,2`, and

```text
P (B9)_P-4B9=-(2UP^6+4).                                  (24)
```

At a boundary point `S=0,B9=0`, `(24)` equals `-6`, proving generic
nonvanishing.  Thus the only positive-genus component has genus six and
nonzero `G_P`.

## 4. Pick and graph completeness

The exact polygon ledgers `(2Area,boundary,interior)` are

```text
W!=0:
 M=(36,10,14), C=(12,8,3), B7=(7,3,3), Kf=(2,4,0),
 global K=0: (43,11,17),             global K!=0: (45,13,17);

W=0:
 M=(24,14,6), B9=(18,8,6), Kf=(2,4,0),
 global K=0: (42,10,17),             global K!=0: (44,12,17). (25)
```

When `W!=0` and `Lambda!=0`, the `K=0` components `R,C,B7` have twelve
`R--C` nodes and one `C--B7` node.  Thus `(V,E,b1)=(3,13,11)`.  When
`K!=0`, one primitive `B7--Kf` attachment changes this to `(4,14,11)`.  Hence

```text
3+3+11=17.                                                (26)
```

For `W=0,K=0`, the components are `R`, the six rational main-face lines,
and `B9`; their twelve `R--line` and six `line--B9` nodes give
`(V,E,b1)=(8,18,11)`.  The optional primitive `B9--Kf` attachment gives
`(9,19,11)`.  Therefore

```text
6+11=17.                                                  (27)
```

Every value matches the corresponding global Pick genus in `(25)`.  No
positive-genus component is omitted, and toric subdivisions add only
rational chains.

## 5. Complete edge schemes and regular-model boundary

For `W!=0`, the `K=0` schemes are

```text
X-1, 1-beta X, -WX-beta, (X-1)(UX+W), U-X^6, 1-WX.        (28)
```

When `K!=0`, replace the exposed edge `1-beta X` by

```text
1-KX^2,                     -beta X-K,                    (29)
```

and retain `1-beta X` as the primitive internal `B7--Kf` edge.  Thus there
are eight active schemes.  Under

```text
U*W*beta*K*Lambda!=0                                        (30)
```

they are active, avoid toric corners, and are squarefree.  When `K=0`, omit
the `K` schemes and its factor from `(30)`.  The top discriminant is
`Lambda^2`; every other derivative loses squarefreeness or a relative-
interior root only on an explicitly excluded owner wall.

For `W=0`, the `K=0` schemes are

```text
X-1, 1-beta X, -UX-beta, U(X-1), U-X^6, 1-UX^6.           (31)
```

For `K!=0`, again insert `(29)` and the internal `1-beta X`.  All are active,
corner-free, and squarefree under `U*beta!=0`, with `K!=0` only in the
`Kf` stratum.  Here `Lambda=U!=0`, so no repeated top root occurs.

The primitive normals `(11)`, torus smoothness `(19),(24)`, edge lists
`(28)--(31)`, and the exact node charts inherited from THM-4327 verify the
face/edge regular-model hypotheses.  Numerical equality of root coordinates
on different edges does not identify their strata: attachments are labelled
by relative interiors of distinct toric divisors.

## 6. `Lambda=0`: `beta_11` forces an early splitter

The only possible `Lambda=0` case has `W=-U!=0`, hence belongs to the
`L=84` regime.  Use the global base `Q=Sigma^84` and the local contact base

```text
Q=sigma^12,                   sigma=Sigma^7,
t=sigma b,                    s=v(sigma)>0,
nu=v(b)>0,                    q=r-1,
b=1/S,                        r=P/S^2.                     (32)
```

Exactly as an identity of the specialized top polynomial,

```text
Ctilde=b^12-U r^5(r-1),       partial_r Ctilde(0,1)=-U.    (33)
```

Thus `R` and `C` have one length-twelve `A_23` contact.  No simple-root
packet is transported through it.

The following local calculation is restated because THM-4297's global
statement has a different coefficient gate.  Put `A=2U+W=U!=0`.  The
prepared equation is

```text
F=q ell+q^2V(q,t)-t^12/2,                 V(0,0)=A.        (34)
```

Before repetition, every generically separable face has Keller-form order at
least

```text
3s+5nu>0,                                                   (34a)
```

the local THM-4292/4297 calculation which depends only on `V(0,0)=A!=0`.

For `nu<s`, any coefficient/`b^12` cancellation ends in discriminant
`x^2+2A`, hence only simple roots.  On `nu=s`, after lower rows vanish, write
`c` for the `t^6` coefficient of `ell`; the discriminant is
`(c-X^6)^2+2A`.  A nonzero multiple root forces `A=0`.  Its possible root at
`X=0`, characterized by `c^2+2A=0` and therefore `c!=0`, is exactly passage
to `nu>s`.

On that repeated range the only difference from the `W=0` quadratic model
is the exact term

```text
Delta F=-(W/2)r^4q^3.                                      (35)
```

After `q=t^6y` and division by `t^12`, it begins at `t^6`.  Deep repetition
requires

```text
c_1=alpha_11+beta_11=0.                                   (36)
```

The `K=[y^2]H` row occurs as `K t^6r^2` in `Hhat`.  After the substitution
and division it becomes

```text
Kyr^2=Ky+2Kt^6y^2+Kt^12y^3.
```

Its leading `Ky` is absorbed into the balanced coefficient `c`, and the
balanced repeated case has `c!=0`, so this absorption cannot annihilate the
quotient below.  Its remaining terms, and all inherited higher-`q` terms,
begin at `t^6`; `(35)` is the only new `W`-dependent correction.

But `beta_11!=0`, so `alpha_11=-beta_11!=0`.  The first moving-critical
coefficient is therefore

```text
C_1=alpha_11/c^2!=0.                                      (37)
```

There are two nonzero candidate splitters, with gaps

```text
d_1=s+nu,                       d_b=6(nu-s).               (38)
```

Their exact comparison is

```text
1<nu/s<7/5:  d_b<d_1,        first order 6s+2nu>0;
nu/s=7/5:    d_b=d_1,         common order
                                  6s+2nu=(5s+9nu)/2>0;
nu/s>7/5:    d_1<d_b,        first order (5s+9nu)/2>0.     (39)
```

At equality the critical-value face is
`X(C_1-(1/c)X^5)`.  Every nonzero root is simple; at such a root the
normalized equations

```text
X-X_0=unit*T^2,                  F_y=unit*T
```

cancel the apparent derivative zero.  The root `X=0` passes to the next
ratio case rather than creating an unlisted component.  Away from equality,
the earlier of `b^12q` and `C_1t` is already nonzero, so no refinement can
pass it.  Both gaps precede every `t^6` correction, including the `K`
remainder and `(35)`:

```text
6(s+nu)-d_b=12s>0,               6(s+nu)-d_1=5(s+nu)>0.   (40)
```

This exhausts every local valuation range.  The base identity
`sigma=Sigma^7` gives `s=7v(Sigma)>0`; substituting this in either positive
form in `(39)` preserves strict positivity.  All exceptional divisorial
components over `(33)`, including any
positive-genus tail, are Keller-constant.  The optional `Kf` attachment lies
in a different toric-edge orbit and remains simple, so it creates no second
contact.

At contact, the arithmetic-genus formulas are

```text
K=0:  (3+3)+12+1-3+1=17,
K!=0: (3+3)+12+1+1-4+1=17,                               (41)
```

where the extra `1` in the second line is the `B7--Kf` node.  These agree
with `(25)` and leave no unaccounted tail.

## 7. Good-form extinction and degree conservation

For a lower plane `h=a_s r+b l+c`, THM-4327 gives

```text
ord_face(phi^*eta_0)=L(5/6-a_s-b-c).                      (42)
```

The exact orders in the two `W` regimes are

```text
W!=0, L=84:  M=63, V7=70, Kf=196;
W=0,  L=36:  M=27, V9=28, Kf=84.                          (43)
```

Every value is positive.  Equations `(20)--(24)` show that the relative
differential generator is nonzero at every positive-genus generic point.
Thus `C,B7`, or `B9` maps constantly to the good elliptic curve.  Every
remaining primary component and ordinary toric chain is rational, and
Section 6 handles the sole repeated contact and all exceptional tails.

After finite base change and regularization, let `M_phi` be the actual
pullback of a positive-degree bundle on the good target.  Proper-flat
intersection, with all fibre multiplicities retained, gives

```text
deg(M_phi|generic fibre)=sum_i m_i deg(M_phi|X_i)=0.       (44)
```

This contradicts the positive response degree of a nonautomorphic Keller
pair and proves the theorem.

## 8. Sharp boundary and the next contact object

The face and local mechanisms first lose force at `beta_11=0`.  With
`W!=0`, the exact beta-zero atlas has only `M` and the optional plane

```text
T=(1/2,0,-1),                         order at L=12: 16.    (45)
```

Its face core contains the cubic edge polynomial

```text
A(P)=K+Theta P+xi_10 P^2+W P^3.                            (46)
```

Primary differential order remains positive, but a repeated nonzero root is
an unaudited contact.  The exact hostile

```text
(K,Theta,xi_10,W)=(1,-1,-1,1),
A(P)=(P-1)^2(P+1)                                         (47)
```

shows the issue is real.  At `W=beta_11=0`, the hostile atlas expands to
twelve complexes and nine planes; its first repeated edge can already be
`U+alpha_11X+xi_10X^2=(X+1)^2`.

THM-4335's LRC owner-address work suggests a precise organization for `(46)`,
not a theorem dependency.  Its symmetric labelled-root separation is

```text
Disc(A)=xi_10^2 Theta^2-4W Theta^3-4xi_10^3K
        -27W^2K^2+18Wxi_10Theta K.                         (48)
```

The regimes `KW!=0,Disc(A)!=0` (three etale torus attachments),
`KW!=0,Disc(A)=0` (interior contact), and `K=0` or `W=0` (a root exits through
a toric endpoint) must be kept distinct.  The required sidecar is the
Laurent-saturated labelled root algebra together with local intersection
multiplicities; support alone and the scalar discriminant alone each lose
necessary information.

After this theorem, the exact `Z=0,U!=0` residual is

```text
beta_11=0 and W*zeta_3=0.                                 (49)
```

The `U=0` repeated cusps and owner intersections, exact-seam entry, `JC(2)`,
and `DC(2)` also remain open.

## 9. Reproduction and honest scope

The primary certificate reuses the pinned THM-4327 `Z=0` hull primitives;
the clean-room audit must not.  Run both in normal and optimized modes and
byte-match their frozen outputs.

```bash
python3 -B 04-computation/jc2_m12_z0_zeta0_endpoint_extinction_thm4337.py
python3 -B -O 04-computation/jc2_m12_z0_zeta0_endpoint_extinction_thm4337.py
python3 -B 04-computation/jc2_m12_z0_zeta0_endpoint_extinction_independent_audit_thm4337.py
python3 -B -O 04-computation/jc2_m12_z0_zeta0_endpoint_extinction_independent_audit_thm4337.py
```

What is proved here is exactly `(1)--(2)`, relative to the inherited
face/edge and proper-flat interfaces.  Packets are consistency controls only
and are never transported through `(33)`.  None of the open claims in
Section 8 is promoted by this theorem. **QED.**
