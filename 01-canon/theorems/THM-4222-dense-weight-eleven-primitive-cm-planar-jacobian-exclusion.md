---
id: THM-4222
title: "Dense weight-eleven primitive-CM planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3992/3997/4012/4045/4218 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED.
  In the inherited b=d=0 reduced (2,3) seam, the complete exact-M=11
  locus with nonzero p^4y, py^3, p^5, y^3 coefficients and separated
  weight-eleven top roots contains no nonautomorphic planar Keller pair.
  Its only positive-genus special component is a genus-five curve with
  primitive Q(zeta_11)-CM and hence simple Jacobian; every component map
  to the good j=0 target is constant, contradicting degree conservation.
  The five named coefficient walls, other cells, seam entry, JC(2), and
  DC(2) remain OPEN.
source: codex-planar-jacobian-weight-eleven-session-20260826
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
related:
  - THM-4217-complete-mixed-off-antidiagonal-delta-zero-planar-jacobian-exclusion
  - THM-4220-weight-ten-zeta-zero-genus-two-planar-jacobian-exclusion
external: >
  Tim Dokchitser, "Models of curves over DVRs," arXiv:1807.00025v2,
  Definitions 3.7, 3.9, and 3.12 and Theorem 3.14, through the audited use
  in THM-4045, supplies the general face/edge model and rational toric
  chains. J. S. Milne, "Complex Multiplication" course notes (2020),
  Sections 1.9--1.10 and Proposition 3.13, supplies the standard theorem
  that primitive CM-pairs classify simple CM abelian varieties up to
  isogeny; https://www.jmilne.org/math/CourseNotes/CM.pdf. The exact face
  arithmetic, CM characters and stabilizer,
  regular-model ledger, and planar-Jacobian consequence are proved here.
script: 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222.py
output: 05-knowledge/results/jc23_dense_weight11_primitive_cm_exclusion_thm4222.out
script_sha256: f8d92c59e754eee75c2a8e101ad1c2f597372406a3fefadfafd717638b384609
output_sha256: bbf4161da4c0f290c384a8a662fff1231588526d53ca57d5c538a3899438ee47
independent_audit_script: 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_dense_weight11_primitive_cm_exclusion_thm4222_independent_audit.out
independent_audit_script_sha256: fe9ed3bcb4651168be37201794b8b5bd8261b28b317a4350a5cbca94099def57
independent_audit_output_sha256: de5d236c496c9e80d42d56864e46cfa72abd45086b177569df6ddfe7a3f5fcdb
---

# THM-4222 -- dense weight-eleven primitive-CM planar Jacobian exclusion

**PROVED RELATIVE TO THM-3992/3997/4012/4045/4218 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.**

## 1. Statement and scope

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Put

```text
A=[p^4y]H,       B=[py^3]H,       U=[p^5]H,       Z=[y^3]H.       (1)
```

> **Theorem.** The complete exact-weight-eleven locus
>
> ```text
> A*B*U*Z*(A+B) != 0                                             (2)
> ```
>
> contains no nonautomorphic planar Keller pair.

Every other lower coefficient is arbitrary. In particular, `K,Phi,Delta,
Theta,eta,Xi` may vanish or cancel. The theorem starts inside the inherited
seam; it proves neither entry into that seam nor `JC(2)`.

The closest mechanism is THM-4220's degree-zero specialization. The
canonical hostile is THM-4218's hidden elliptic side face. The corrected near
miss is to treat the order-eleven attachment difference as a genus-one degree
invoice: killing one torsion class in a genus-five Jacobian does not force a
curve-map degree to be divisible by eleven. The decisive sidecar is instead
the complete cyclotomic CM type and its stabilizer.

## 2. Complete source and universal lower hull

Use `s=XT`, `p=T+s^2`, `y=sp`, and `t=p-s^2=T`. Enumerating
`0<2i+3j<=11` and deleting exactly the forbidden `y,py` rows gives

```text
p,p^2,p^3,y^2,p^2y,p^4,py^2,p^3y,y^3,
p^5,p^2y^2,p^4y,py^3.                                      (3)
```

Thus there is no ellipsis in

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+Z y^3
  +U p^5+Xi p^2y^2+A p^4y+B py^3,

K=2848/45-(7/6)Delta.                                      (4)
```

The weight-eleven part in `(s,p)` is `A s p^5+B s^3p^4`. For a generic
pencil value set `Q=q^-1` and

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                                  (5)
```

On its torus, `p=s^2` is impossible because it would give
`F_Q=-Qs^2/2!=0`; hence `t=p-s^2`, `X=s/t` recovers the actual source
function field.

A term `p^iy^j` contributes valued endpoints

```text
(j+2,i+j,1),                       (j,i+j+1,1).           (6)
```

For `(r,l)` the only lower planes under `(2)` are

```text
M: nu=(r+2l-2)/11,       V: nu=(l-1)/5,       T: nu=(r-2)/3. (7)
```

The two endpoint gaps above `M` are both `(11-2i-3j)/11`; above `V` they
are `(6-i-j)/5,(5-i-j)/5`; above `T` they are `(3-j)/3,(5-j)/3`.
Consequently equality off `M` occurs only for the second endpoints of
`p^5,p^4y` on `V` and the first endpoints of `y^3,py^3` on `T`.

The retained terms already have exactly the three faces `(7)`, and every
optional endpoint lies above their lower envelope. The certificate checks all
`2^9=512` lower-support subsets and all 32 deletion patterns among the five
coincident support points

```text
(2,3,1),(2,4,1),(2,5,1),(3,4,1),(3,5,1).               (8)
```

The last point has coefficient `B-A` and is only the midpoint of a retained
edge. Thus `A=B` is harmless when `(2)` holds.

The face polygons and Pick ledgers are

| face | polygon | `(2Area,boundary,interior)` |
|---|---|---:|
| `M` | `(0,1),(2,0),(5,4),(1,6)` | `(33,5,15)` |
| `T` | `(2,0),(5,3),(5,4)` | `(3,5,0)` |
| `V` | `(0,1),(1,6),(0,6)` | `(5,7,0)` |

Their equations, up to nonzero monomials, are

```text
g_M=(S^2-P)(1-A S P^5-B S^3P^4)=R*C,
g_T/S^2=1-(SP)^3(Z+B P),
g_V/P=-1+P^5(U+A S).                                  (9)
```

Both side curves are rational: after `T_0=SP`, `g_T` is linear in `P`, and
`g_V` is linear in `S`. Only `C` can have positive genus.

## 3. Polygon, packet, and all eight edge schemes

The global polygon and Pick ledger are

```text
(0,1),(2,0),(5,3),(5,4),(1,6),(0,6),
(2Area,boundary,g)=(41,13,15).                          (10)
```

Using the intrinsic index `e=u+v-c` on an edge with primitive inward normal
`(u,v)` and support constant `c`, and omitting the affine divisor `s=0`, the
labelled outer inventory is

```text
AB: e1 on R;        BC: e2,e2,e2 on T;       CD: e4 on T;
DE: e10 on R and e10 on C;                   EF: e5 on V. (11)
```

Therefore

```text
packet=(10,10,5,4,2,2,2,1),
n_full=36,                    defect=28=2*15-2.         (12)
```

The six outer and two internal special edge schemes, up to units, monomials,
and reversal, are

```text
AB X-1,          BC 1-ZX^3,       CD Z+BX,
DE (X-1)(AX+B),  EF A+UX,         FA U-X^5,
AE AX-1,         BD 1-BX.                                  (13)
```

The nonconstant discriminants are `-27Z^2,(A+B)^2,5^5U^4`; all schemes are
reduced and avoid corners under exactly `(2)`. There is no extra wall
`U+A=0` or `Z+B=0`.

The `BC` restriction is the prime cubic carrier

```text
q-1/2=K W^2+Z W^3.                                    (14)
```

It has three geometric index-two branches, so removing it from the origin
fibre changes the packet sum from `36` to `30`. Both are divisible by three;
THM-4218's congruence is genuinely unavailable. The proof below does not use
either response degree.

## 4. Exact regular model

Take the common base change

```text
Q=sigma^330.                                           (15)
```

The integral lower heights and primitive normals are

```text
M:30r+60l-60       (30,60,-1),
V:66l-66           (0,66,-1),
T:110r-220         (110,0,-1).                         (16)
```

Every face therefore has multiplicity one. The six outer three-dimensional
edge gcds are `1,3,1,2,1,5`, equal to their planar lengths; both internal
gcds are one. The outer slope ledger

```text
30,-110,-110,-30,-66,0                                (17)
```

adds no chain. The determinant-one internal sequences are

```text
AE: -60>-61>...>-66,       five intermediate components;
BD: -90>-91>...>-110,      nineteen intermediate components. (18)
```

Every inserted component is rational of multiplicity one.

The branches `R,C` meet exactly where

```text
P=S^2,                         1-(A+B)S^11=0.           (19)
```

Their gradient determinant is `-11(A+B)S^10`; hence these are eleven
distinct transverse torus nodes. In the main chart set

```text
s=sigma^-30S,       p=sigma^-60P,
H_sigma=sigma^330H(sigma^-60P,sigma^-90SP),
U_0=S^2-P,          V_0=(1-H_sigma)/S^2.               (20)
```

The exact scaled equation is

```text
U_0V_0=sigma^330/2.                                    (21)
```

Thus each point of `(19)` is an `A_329` smoothing; its resolution inserts
329 rational multiplicity-one curves along that path.

There are no other singularities. On `C`, the two derivative equations have
coefficient determinant `-11AB`; on `V`, differentiation by `S` gives
`AP^5`; on `T`, coordinates `(T_0=SP,P)` give derivative `-BT_0^3`.
The edge schemes `(13)` cover every compactified boundary point. The exact
side charts

```text
s=S,p=sigma^-66P;               s=sigma^-110S,p=P       (22)
```

reduce to `g_V,g_T` after multiplying `(5)` by `sigma^66,sigma^220`.
Dokchitser's face/edge theorem, through its audited use in THM-4045, now
supplies the proper regular semistable model after the explicit resolutions
`(21)`.

The generic outer restrictions are linear apart from the separated `DE`
roots, the prime cubic `(14)`, and the affine polynomial `H(p,0)-q` of degree
five. The last is coprime to its derivative over `C(q)`, since a derivative
root is constant and cannot have transcendental image `q`. Thus the generic
completion is smooth and, by the inverse after `(5)`, is the actual source
completion.

After contracting rational paths, the core has vertices `R,C,T,V`, eleven
`R--C` paths, one `C--T` path, and one `C--V` path. Therefore

```text
b_1=13-4+1=10,                    g_special=5+10=15.    (23)
```

This matches `(10)` and proves component and packet completeness. Every
component except `C` is rational.

## 5. Primitive CM and the Hom obstruction

The only nonrational normalization is

```text
C:1-A S P^5-B S^3P^4=0.                               (24)
```

Put `x=ASP^5`. The exact identity

```text
P^11=(B/A^3)x^3/(1-x)                                 (25)
```

presents `C` as a cyclic degree-eleven cover of the `x`-line, with primitive
branch residues `3,10,9` at `0,1,infinity`. Riemann--Hurwitz gives `g(C)=5`.
Its Newton triangle `(0,0),(1,5),(3,4)` has interior points

```text
(1,2),(1,3),(1,4),(2,3),(2,4).                         (26)
```

For the covariant Jacobian convention used here, the deck action

```text
tau:(S,P) |-> (zeta_11^6S,zeta_11P)                    (27)
```

acts on the tangent/regular-differential CM type with characters

```text
Phi={4,5,8,9,10}.                                      (28)
```

The contravariant convention replaces `Phi` by `-Phi`; its stabilizer and
the primitivity conclusion are unchanged.

Their negatives are `{1,2,3,6,7}`, so `H^1(C)` contains every nonzero
character once and `Q[tau]=Q(zeta_11)` has degree `10=2g`. Exact
multiplication by the ten units modulo eleven gives

```text
{u in (Z/11)^*:uPhi=Phi}={1}.                          (29)
```

Thus the CM-pair is primitive. Milne, *Complex Multiplication*, Sections
1.9--1.10 and Proposition 3.13, is the explicit cited input that primitive
CM-pairs classify simple CM abelian varieties up to isogeny. The in-repo
calculation `(25)--(29)`, not the citation, supplies this curve's CM type and
primitivity. Hence `J(C)` is simple.

Since `dim J(C)=5`, any nonzero homomorphism `J(C)->E` to an elliptic curve
would have a positive-dimensional connected kernel, contradicting simplicity.
Thus

```text
Hom(J(C),E)=0                                           (30)
```

for every elliptic curve `E`; every morphism `C->E` is constant.

The side paths attach above `x=1` and `x=0`. The identity

```text
div(x/(1-x))=11(P_0-P_1)                               (31)
```

shows their difference has exact order eleven (order one would give a
degree-one function on a genus-five curve). This label is not used as a
false genus-one degree invoice; `(30)` is stronger.

## 6. Good target reduction and degree zero

Since `q=sigma^-330`, put

```text
A_target=sigma^-110X,                  C_target=sigma^-165Y. (32)
```

The inherited target becomes

```text
Y^2=X^3+1-(3a^2/4)sigma^220X-(a^3/4)sigma^330,         (33)
```

with smooth special fibre `E_0:Y^2=X^3+1`.

A hypothetical Keller pair induces a finite nonconstant generic morphism to
this target. Resolve its rational extension over the regular source model.
Every exceptional component is rational. The components `R,T,V`, every
toric and `A_329` chain, and every exceptional component map constantly to
`E_0` by Riemann--Hurwitz; `(30)` makes `C` constant as well.

For a relative target line bundle `L` of degree one, proper flat degree
conservation and the multiplicity-one ledger give

```text
deg(phi_generic)=sum_i m_i deg(phi_i^*L)=0.             (34)
```

This contradicts generic finiteness and proves the theorem.

## 7. Sharp boundary and verification

Each factor in `(2)` has a distinct role:

```text
A=0:    loses the p^4y top vertex and two-term CM face;
B=0:    loses the py^3 top vertex and two-term CM face;
U=0:    loses the V owner and degree-five affine endpoint;
Z=0:    loses the T owner and cubic carrier;
A+B=0:  collides the DE roots and moves R--C contact to the boundary. (35)
```

These are sharp walls for this certificate and remain **OPEN**, not asserted
survivor loci. Conversely `A=B`, `K=0`, and every wall involving only
`Phi,Delta,Theta,eta,Xi` are included. The next contraction atlas is finite:
with highest pure `p` owner `c_kp^k`, the replacement side is the rational
curve `-1+c_kP^k+ASP^5`; when `Z=0,K!=0`, it is the rational curve
`1-K(SP)^2-B(SP)^3P`. No such wall is claimed here.

The primary certificate checks the complete monomial universe, 512 support
subsets, 32 collision hostiles, all face and Pick ledgers, packet, eight edge
schemes, primitive normals and edge denominators, slope sequences, exact
main/side charts, eleven `A_329` models, genus completeness, cyclic-cover
identity, CM characters and stabilizer, generic edges, target scaling, and
degree-zero inventory. It marks Milne Proposition 3.13 as a cited input.
The clean-room audit independently reconstructs the facets, side charts,
chain sequences, Chevalley--Weil character spectrum, CM-type stabilizer,
Milne gates, and degree-conservation ledger.

Replay with

```bash
python3 -B 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222.py
python3 -B -O 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222.py
PYTHONHASHSEED=4222 python3 -B 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222.py

python3 -B 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222_independent_audit.py
python3 -B -O 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222_independent_audit.py
PYTHONHASHSEED=4222 python3 -B 04-computation/jc23_dense_weight11_primitive_cm_exclusion_thm4222_independent_audit.py
```

All three streams byte-match the frozen output. **QED.**
