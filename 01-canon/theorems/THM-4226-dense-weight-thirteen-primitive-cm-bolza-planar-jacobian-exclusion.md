---
id: THM-4226
title: "Dense weight-thirteen primitive-CM Bolza planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4222 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED. In the inherited b=d=0 reduced (2,3) seam, the
  complete exact-M=13 locus with nonzero p^5y, p^2y^3, p^6, y^4
  coefficients and separated main roots contains no nonautomorphic planar
  Keller pair. The positive-genus inventory is a primitive-CM genus-six
  component plus a Bolza genus-two tail, both Hom-orthogonal to the good
  j=0 target. The five coefficient walls, other cells, seam entry, JC(2),
  and DC(2) remain OPEN.
source: codex-jc-m13-dense-audit-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4222-dense-weight-eleven-primitive-cm-planar-jacobian-exclusion
related:
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
  - THM-4220-weight-ten-zeta-zero-genus-two-planar-jacobian-exclusion
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through the audited use
  in THM-4045, supplies the face/edge model and rational toric chains.
  J. S. Milne, "Complex Multiplication" course notes (2020), Sections
  1.9--1.10 and Proposition 3.13, supplies primitive CM-pair simplicity.
  All face arithmetic, CM characters, Bolza identification, regular-model
  ledgers, and the planar-Jacobian consequence are checked below.
script: 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py
output: 05-knowledge/results/jc23_dense_weight13_primitive_cm_bolza_candidate.out
script_sha256: 280d2c2515fe4d1b3f5f229947319d0922b3982d83ad88fd39cd498b78e8404b
output_sha256: fd9b297486e923fe6946eb7a98fb5bd09f7d22080fc354c32202f993c4d1dc96
independent_audit_script: 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.out
independent_audit_script_sha256: 0eac79869bfbef796774531a1cc3607e1badbef8b165db473b18bf6b2a90280c
independent_audit_output_sha256: 774e328d0fda469e51adc956daef4df4e72799ae8700e6890f63e2601ed16d58
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy certificate executes 16,913 exact checks; a
  standard-library-only implementation independently executes 113 checks.
  Both cover the monomial universe, hulls, collisions, Pick and packet
  ledgers, edge schemes, primitive slopes, regular charts, CM type, Bolza
  tail, and genus inventory. Normal, optimized, and fixed-hash-seed streams
  byte-match their frozen outputs. The proof/citation layer is explicit.
---

# THM-4226 -- dense weight-thirteen primitive-CM Bolza exclusion

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4222 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Put

```text
A=[p^5y]H,       B=[p^2y^3]H,
U=[p^6]H,        Z=[y^4]H.                              (1)
```

> **Theorem.** The complete exact-weight-thirteen locus
>
> ```text
> A*B*U*Z*(A+B) != 0                                   (2)
> ```
>
> contains no nonautomorphic planar Keller pair.

Every other lower coefficient is arbitrary subject only to the inherited
normal form. There is no hidden `Z+B`, `A+U`, or `A-B` condition. The five
walls in `(2)`, other cells, seam entry, `JC(2)`, and `DC(2)` remain open.

The closest mechanism is THM-4222's primitive-CM total-`Hom` obstruction.
The hostile example is THM-4218's `j=0` positive-genus tail. The corrected
near miss is that THM-4222's internal slope probes cannot be reused: the new
`BD` edge has primitive transverse step `(3,2)` and two, not five,
intermediate components. The least-used sidecar is THM-4012's exact Bolza
splitting and CM mismatch.

## 2. Complete source and lower hull

Use `s=XT`, `p=T+s^2`, `y=sp`, `t=p-s^2=T`. The complete polynomial is

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+zeta_3 y^3
  +upsilon_5 p^5+xi_10 p^2y^2+alpha_11 p^4y+beta_11 py^3
  +U p^6+omega_12 p^3y^2+Z y^4+A p^5y+B p^2y^3,

K=2848/45-(7/6)Delta.                                  (3)
```

Enumerating `0<2i+3j<=13` and deleting the forbidden `y,py` rows gives
exactly these eighteen monomials. For `Q=q^-1`, set

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                              (4)
```

On its torus `p=s^2` is impossible, so `t=p-s^2` and `X=s/t` recover the
actual source function field.

A monomial `p^iy^j` contributes valued endpoints

```text
(j+2,i+j,1),                  (j,i+j+1,1).             (5)
```

The three lower planes are

```text
M: nu=(r+2l-2)/13,
V: nu=(l-1)/6,
T: nu=(r+l-2)/8.                                      (6)
```

The two endpoint gaps above them are

```text
M: ((13-2i-3j)/13, (13-2i-3j)/13),
V: ((7-i-j)/6,       (6-i-j)/6),
T: ((8-i-2j)/8,      (9-i-2j)/8).                     (7)
```

All gaps are nonnegative. Equality occurs exactly at

```text
M: p^5y,p^2y^3;       V: p^6,p^5y;       T: y^4,p^2y^3. (8)
```

Thus the required rows own every face independently of all lower
coefficients. Of the eight coincident support points, only `(3,6,1)`, with
coefficient `B-A`, lies on a lower face, and it is the midpoint of a retained
edge. Hence `A=B` is harmless. The additional unit `A+B` instead separates
the main intersections below.

## 3. Faces, edges, and genus ledgers

The face polygons and Pick data are

| face | polygon | `(2Area,boundary,interior)` |
|---|---|---:|
| `M` | `(0,1),(2,0),(5,5),(1,7)` | `(39,5,18)` |
| `V` | `(0,1),(1,7),(0,7)` | `(6,8,0)` |
| `T` | `(2,0),(6,4),(5,5)` | `(8,6,2)` |

Their equations, up to units and monomials, are

```text
g_M=(S^2-P)(1-A S P^6-B S^3P^5)=R*C,
g_V/P=-1+P^6(U+A S),
g_T/S^2=1-Z S^4P^4-B S^3P^5.                          (9)
```

The vertical face is rational. The tail is torus-smooth because, after
`x=SP`, differentiation of `1-Zx^4-Bx^3P^2` with respect to `P` is
`-2Bx^3P!=0`.

The global polygon and packet are

```text
(0,1),(2,0),(6,4),(5,5),(1,7),(0,7),
(2Area,boundary,g)=(53,15,20),                         (10)

packet=(12,12,8,6,2,2,2,2,1),
sum=47,                         defect=38=2*20-2.       (11)
```

The six outer and two internal edge schemes are

```text
A_0B_0: X-1,                    B_0C_0: 1-ZX^4,
C_0D_0: -Z-BX,                 D_0E_0: (X-1)(AX+B),
E_0F_0: A+UX,                  F_0A_0: U-X^6,
A_0E_0: AX-1,                  B_0D_0: 1-BX.           (12)
```

Their nonconstant discriminants are `-256Z^3,(A+B)^2,46656U^5`.
Therefore `(2)` is exactly the reduced corner-avoiding gate. In particular,
`A=B` gives `A(X-1)(X+1)` and introduces no collision; `Z+B` and `A+U` are
not units in this model.

The only main-face torus intersections satisfy

```text
P=S^2,                         1-(A+B)S^13=0.           (13)
```

Their gradient determinant is `-13(A+B)S^12`, so there are thirteen distinct
transverse points. The exponent determinant on `C` is `-13`, and `A*B!=0`
makes `C` smooth. The tail derivative, the eight edge schemes, and derivative
`AP^6` on `V` cover every other torus and compactified boundary point.

## 4. Multiplicity-one regular model

Take

```text
Q=sigma^312.                                             (14)
```

The integral lower heights and primitive normals are

```text
M:24r+48l-48,           (24,48,-1),
V:52l-52,               (0,52,-1),
T:39r+39l-78,           (39,39,-1).                    (15)
```

Every face has multiplicity one. The six outer three-dimensional edge gcds
equal the planar lengths `1,4,1,2,1,6`; both internal graph edges are
primitive. The outer slopes `24,-39,-39,-24,-52,0` add no chain.

On `A_0E_0`, the sequence `-48>-49>-50>-51>-52` adds three intermediate
rational components. On `B_0D_0`, the primitive functional

```text
L*=-5r+3l+10,       L*(B_0)=0,       L*(3,2)=1          (16)
```

gives adjacent slopes `120,117`; `120>119>118>117` adds exactly two. All
have multiplicity one.

In the main chart set

```text
s=sigma^-24S,       p=sigma^-48P,
H_M=sigma^312H(sigma^-24S,sigma^-48P),
U_0=S^2-P,          V_0=(1-H_M)/S^2.                   (17)
```

The exact scaled equation is

```text
U_0V_0=sigma^312/2.                                    (18)
```

Thus the thirteen points in `(13)` are `A_311` smoothings whose resolutions
add only rational multiplicity-one chains.

The side charts

```text
V: s=S,             p=sigma^-52P;
T: s=sigma^-39S,    p=sigma^-39P                       (19)
```

reduce exactly to the other faces in `(9)` after multiplying `(4)` by
`sigma^52,sigma^78`. Face smoothness, the reduced edge schemes, primitive
denominators, and these charts meet Dokchitser's face/edge regularity gate
through THM-4045.

The moving outer restrictions are

```text
q-1/2=K W^2+zeta_3 W^3+Z W^4,
q=H(p,0),                         deg_p H(p,0)=6.       (20)
```

They are squarefree over `C(q)`: a common derivative root is constant over
the coefficient field and cannot have transcendental image `q`. The generic
completion is smooth and is the actual source completion by the inverse
after `(4)`.

## 5. The two positive-genus components

The nonrational main normalization is

```text
C:1-A S P^6-B S^3P^5=0.                               (21)
```

Put `x=ASP^6`. Then

```text
P^13=(B/A^3)x^3/(1-x).                                 (22)
```

This cyclic degree-thirteen cover has branch residues `(3,12,11)` and genus
six. Its Newton triangle has interiors

```text
(1,2),(1,3),(1,4),(1,5),(2,4),(2,5).                  (23)
```

The deck action

```text
tau:(S,P) |-> (zeta_13^7S,zeta_13P)                   (24)
```

has holomorphic character set

```text
Phi={5,6,9,10,11,12}.                                  (25)
```

Its negatives complete all twelve nontrivial characters, so
`Q[tau]=Q(zeta_13)` has degree `12=2g`. Exact multiplication by the twelve
units gives stabilizer `{1}`. The CM type is primitive; Milne Proposition
3.13 makes `J(C)` simple. Since it has dimension six,

```text
Hom(J(C),E)=0                                           (26)
```

for every elliptic curve `E`.

For the tail, put

```text
x=SP,                         Y=x^2P.                  (27)
```

Its normalization becomes

```text
B Y^2=x-Zx^5=x(1-Zx^4).                                (28)
```

The quintic is squarefree under `(2)` and has genus two. Choose `u,v` with
`Z=-u^4`, `B=v^2/u`; scaling identifies it with the Bolza curve

```text
y^2=x(x^4+1).                                          (29)
```

THM-4012 proves that its Jacobian is isogenous to `E_8000^2` and that
`Hom(E_8000,E_0)=0` for `E_0:Y^2=X^3+1`. Hence

```text
Hom(J(Bolza),E_0)=0.                                   (30)
```

## 6. Genus completeness and degree zero

After contracting rational paths, the core has vertices `R,C,Bolza,V`,
thirteen `R--C` paths, one `C--Bolza` path, and one `C--V` path. Thus

```text
b_1=15-4+1=12,
g_special=6+2+12=20.                                  (31)
```

This equals the global Pick genus and proves component completeness. Every
other component is rational.

Since `q=sigma^-312`, scale the inherited target by

```text
A_target=sigma^-104X,                C_target=sigma^-156Y. (32)
```

The exact target equation becomes

```text
Y^2=X^3+1-(3a^2/4)sigma^208X-(a^3/4)sigma^312,         (33)
```

with smooth special fibre `E_0`. A hypothetical Keller pair induces a finite
nonconstant generic map to this target. Resolve its rational extension over
the regular source model. Riemann--Hurwitz makes every rational component
constant; `(26)` and `(30)` make both positive-genus components constant.

For a relative degree-one target line bundle, proper-flat degree conservation
and the multiplicity-one ledger give

```text
deg(phi_generic)=sum_i m_i deg(phi_i^*L)=0,             (34)
```

contradicting generic finiteness. This proves the theorem.

## 7. Verification and firewall

The primary certificate checks all eighteen monomials, `2^14` lower-row
subsets, `2^8` collision patterns, every polygon and edge scheme, exact unit
gate, primitive slopes, charts, thirteen `A_311` models, CM characters,
Bolza identity, target scaling, and the complete genus inventory. The
independent standard-library audit reconstructs those ledgers without SymPy.

Replay with

```bash
python3 -B 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py
python3 -B -O 04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py
PYTHONHASHSEED=4213 python3 -B \
  04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate.py

python3 -B \
  04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.py
python3 -B -O \
  04-computation/jc23_dense_weight13_primitive_cm_bolza_candidate_independent_audit.py
```

All streams byte-match the frozen outputs. The scripts do not prove Milne's
theorem, THM-4012's Hom calculation, the inherited Keller-map interface, or
proper-flat degree conservation; those are the explicit cited/proof layer.
Nothing here closes a wall in `(2)`, another cell, seam entry, `JC(2)`, or
`DC(2)`. **QED.**
