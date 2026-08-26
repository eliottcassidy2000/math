---
id: THM-4176
title: "Complete repeated-top wall planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/
  4157/4171/4173 + VERIFIED-EXACT + INDEPENDENTLY RESULTANT-AUDITED. On the live
  exact-weight-nine reduced (2,3) seam, the whole repeated-top wall
  zeta=-eta with eta*Delta!=0 contains no nonautomorphic planar Keller pair.
  Rows A and D are inherited from THM-4173 and THM-4157. Exact B/C
  resultant-divisor towers, the row-independent Hessian bridge, and the
  THM-4171 carrier-orbit lemma remove every remaining discriminant,
  projected-fibre-separation, and leading-degree wall in rows B and C.
source: codex-frontier-synthesis-creative-20260826ax
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4157-repeated-top-edge-wall-planar-jacobian-exclusion
  - THM-4171-row-a-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4173-repeated-top-row-a-complete-planar-jacobian-exclusion
related:
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
script: 04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.py
output: 05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.out
independent_audit_script: 04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.gp
independent_audit_output: 05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.out
script_sha256: 928d77ae3a9815abaa38cf7ca6e29716f5d3f1575fc46e68ec93481baf59cc98
output_sha256: 7b32ef45ab66500a181215a1c14d6beb4d03d1af8fed2569fd63c30cdd199d79
independent_audit_script_sha256: 37d5f1da666e2ed023987d9e0ab8cbd376bb5e53dfd9fec40422ad0c65a0b734
independent_audit_output_sha256: a6edc0762e692e054dea56553cb60edbc2f587fdfd0852655d2e4b4b4d39eb59
hash_basis: raw LF bytes
primary_audit: >
  PASS. The exact SymPy certificate reconstructs the complete B/C source,
  proves the row-independent Hessian bridge, checks both finite-T X-infinity
  sidecars, and recomputes both resultants and their coefficient towers
  through the terminal coprimality gate. Its degree-to-length and monodromy
  checks use the inherited universal T=0 identity proved and frozen in
  THM-4157/4173; it does not directly recompute that pair. Normal, optimized,
  and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT. PARI/GP independently reconstructs only the B/C resultant towers.
  It verifies the common leading wall, the row-B double-root and terminal
  tower, the row-C terminal specialization, and the terminal
  gcd(H,Delta*U*V*N)=1 (script variables gcd(H,D*A*b*N)=1). It does not
  independently audit the Hessian bridge, T=0 pair, X-infinity sidecars,
  packets, carrier transport, or monodromy.
---

# THM-4176 -- complete repeated-top wall planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/
4147/4157/4171/4173 + VERIFIED-EXACT + INDEPENDENTLY RESULTANT-AUDITED;
JC(2) AND DC(2) REMAIN OPEN.**

## 1. Exact theorem and inheritance pass

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. For every coefficient tuple

```text
(Delta,Theta,Phi,eta) in C^4,             eta*Delta!=0,                (1)
```

put

```text
zeta=-eta,
K=2848/45-(7/6)Delta,
Ctop=Delta+Theta,
Btop=-1376/135+K=7(2048-45Delta)/270.                  (2)
```

Use the complete normalized source

```text
P=T+X^2T^2,                            Y=XTP,

G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4+Theta P Y^2
  +eta P^3Y-eta Y^3.                                      (3)
```

> **Theorem.** No coefficient tuple satisfying `(1)--(3)` is the normalized
> exact-weight-nine source of a nonautomorphic planar Keller pair in the
> inherited reduced `(2,3)` seam.
>
> Equivalently, the complete repeated-top wall
>
> ```text
> zeta=-eta,                         eta*Delta!=0                       (4)
> ```
>
> is empty of nonautomorphic planar Keller realizations.

The quantifier is over **every** complex coefficient point in `(4)`. No
hypothesis is imposed on a coordinate eliminant discriminant, on the value of
an eliminant at `T=-1/6`, on either independent source projection, or on the
leading coefficient of a nominal residual degree.

The four exhaustive, pairwise disjoint structural rows are

```text
A: Ctop!=0;
B: Ctop=0, Phi!=0;
C: Ctop=Phi=0, Btop!=0;
D: Ctop=Phi=Btop=0.                                      (5)
```

In row `D`, `(2)` forces

```text
Delta=2048/45,       Theta=-2048/45,       K=1376/135. (6)
```

The closest proved source/packet/carrier mechanism is THM-4157. THM-4173,
using THM-4171 on its inner wall, already proves all of row `A`. THM-4157
already proves all of row `D`, including its endpoint wall `J_D=0` and its
common nonreduced wall `H_30=0`. The new content is the exhaustive exact
degree-drop analysis in rows `B,C` and its use with the row-independent local
bridge below.

The canonical hostile is MISTAKE-421: several reduced critical points may
share one projected `T`-coordinate, so full etaleness does not make a chosen
coordinate eliminant squarefree. MISTAKE-450 and MISTAKE-445 require separate
degree-preserving and degree-drop strata and forbid calling a fixed top-slot
coefficient the leading coefficient after it vanishes. MISTAKE-509 requires
the three conjugates of the prime cubic carrier to respond together.

## 2. Source-target contract and retained data

The statement is relative to the inherited Keller/elliptic-pencil contract:

```text
Keller source map:
  (A,C_K): A^2_C -> A^2_C,              det D(A,C_K)=1;

target pencil:
  E=q, with precisely the two inherited Morse nodes at values 0 and 1/2;

source pencil:
  the normalization of the generic fibre G=q of (3), together with all
  affine inverse points of the two nodal target values;

map of generic fibres:
  the map induced by (A,C_K), after compactification, resolution, relative
  normalization, and restriction to a finite morphism;

retained sidecars:
  complete infinity packet, residue-field degree of the nonrational place,
  total finite-carrier index beta, affine Morse inverse count L, finite
  origin index, full packet defect, and handle-support overlap;

discarded data:
  ordering and coordinates of the three carrier conjugates, ordering of
  affine critical points, and the value of a chosen projection coordinate.
                                                                    (7)
```

The proof order is load-bearing. First resolve and shrink to proper smooth
families with a finite relative map. Only then delete the target origin, the
complete cubic carrier and its full inverse image, and the finite Hurwitz
exceptional set. The resulting restriction is finite etale. Parallel Milnor
cores transport every affine Morse inverse point as a distinct fixed sheet.
No assertion that a punctured family itself is proper is used.

At an affine critical point of a hypothetical realization `G=E(A,C_K)`,

```text
grad G=D(A,C_K)^t grad E=0.
```

Since `det D(A,C_K)=1`, the second-derivative chain rule at a target node is

```text
Hess(G)=D(A,C_K)^t Hess(E)D(A,C_K),
det Hess(G)=det Hess(E)!=0.                              (8)
```

Thus every affine critical point of a hypothetical Keller realization is
Morse, and every such point lies over one of the two target nodes.

## 3. Row-independent Morse-resultant divisor lemma

> **Lemma.** Let `k` be algebraically closed of characteristic zero, let
> `G in k[X,T]` satisfy `T | G_X`, and put
>
> ```text
> f=G_X/T,                         h=G_T.               (9)
> ```
>
> Fix a coefficient specialization for which `Res_X(f,h)` is nonzero.
> Suppose that, for every `t in k^*`, the projective `X`-closures of `f=0`
> and `h=0` have no common point at `X=infinity`. If every common zero of
> `f,h` on `T!=0` is a Morse critical point of `G`, then:
>
> 1. `V(f,h) intersect {T!=0}` is finite and reduced;
> 2. for every `t in k^*`,
>
>    ```text
>    ord_(T=t) Res_X(f,h)
>      =sum_(z:T(z)=t) i_z(f,h)
>      =#{z:T(z)=t};                                      (10)
>    ```
>
> 3. if
>
>    ```text
>    Res_X(f,h)=T^v U(T),                  U(0)!=0,      (11)
>    ```
>
>    then the number of distinct common zeros on `T!=0` is the **actual**
>    polynomial degree `deg U`, counted after specialization.

**Proof.** Direct differentiation gives the polynomial identity

```text
T det D(f,h)=det Hess(G)+f G_XT.                        (12)
```

At a common zero of `f,h` with `T!=0`, `(12)` and the Morse hypothesis give
`det D(f,h)!=0`. Hence the local complete intersection has length one. The
no-common-`X`-infinity hypothesis identifies the `t`-adic resultant valuation
with the sum of the finite local intersection multiplicities, proving
`(10)`. Since `U(0)!=0`, every finite root of `U` lies on `T!=0`, and the sum
of its root multiplicities is its actual degree. QED.

This lemma distinguishes two reduced points with one `T`-coordinate from one
doubled source point. A repeated projected root is harmless in the first
case; the second case would make `(12)` force a zero Keller Hessian. Therefore
neither `disc_T U!=0` nor separation from another known fibre is a hypothesis.

For `(3)`, independently of the structural row,

```text
f(X,0)=-X,                         h(X,0)=-(X^2+6)/2.   (13)
```

The actual critical ideal is `(Tf,h)`, so it restores exactly

```text
T=0,      X^2=-6,      G=0,        det Hess(G)=+6.    (14)
```

The factor `(6T+1)^2` occurring below is exactly

```text
T=-1/6,   X^2=6,       G=1/2,      det Hess(G)=-6.    (15)
```

If a residual factor also vanishes at `T=-1/6`, `(10)--(12)` show that its
extra valuation comes from additional distinct points in that fibre; it is
not lost affine length.

## 4. Inherited complete rows A and D

### 4.1 Row A

Put

```text
D_A=4Theta K^2-27eta^2.                                (16)
```

THM-4173 proves `D_A!=0` using

```text
Res_X(f,h)=T^42(6T+1)^2Q_19(T),
Q_19(0)=-12288 Ctop^6,
[T^19]Q_19=-1458 Ctop eta^4D_A^2,                     (17)
```

the lemma above, and the complete packet/carrier response. It gives
`L=19+2+2=23`. THM-4171 proves the entire `D_A=0` locus. Its four exhaustive
source-resultant strata have residual degrees and lengths

```text
18,17,16,15       ->       L=22,21,20,19.              (18)
```

All retain genus ten and packet `(7,7,4,2,2,2,1)`. Thus row `A` is complete,
with no residual discriminant or `Q_19(-1/6)` hypothesis.

### 4.2 Row D

THM-4157 gives genus seven, packet

```text
(7,4,2,2,2,1),                        defect=12,       (19)
```

and generic residual degree twelve, hence `L=16`. On

```text
J_D(eta)=22143375eta^2+15510536192=0,                 (20)
```

the actual residual degree is eleven and `L=15`. On the remaining common
discriminant wall `H_30(eta)=0`, the exact linear `X`-subresultant isolates a
unique finite point above a repeated `T`-root; its local intersection length
is at least two, so `(12)` forces `det Hess(G)=0`, contradicting `(8)`.
Projection switching removes the projection-only factors `H_24,H_36`.
Consequently THM-4157 proves all of row `D` for `eta!=0`.

## 5. Complete row-B degree tower

Assume

```text
Ctop=0,                         Phi*eta*Delta!=0.       (21)
```

Set

```text
U=105Delta-5696,
V=825Delta-22784,
J_BC=Delta U^2+54675eta^2.                             (22)
```

Exact elimination in the complete source `(3)` gives

```text
deg_X(f,h)=(6,7),
[X^6]f=7T^6(Phi+eta T),
[X^7]h=T^6(7Phi+8eta T).                              (23)
```

The two leading coefficients in `(23)` never vanish simultaneously for
`T!=0` when `Phi*eta!=0`; hence the lemma has no finite-`T` point at
`X=infinity` anywhere in row `B`.

The exact resultant identity is

```text
Res_X(f,h)=T^30(6T+1)^2Q_B(T),
Q_B(0)=-(3*7^5/2^5)Phi^5,
[T^17]Q_B=eta^3 J_BC^2/22500.                         (24)
```

In the quotient by `J_BC`, put

```text
F=180Delta eta V^2-4Delta U^2V Phi-1215eta U^2Phi^2.
                                                                    (25)
```

Then

```text
[T^16]Q_B = Delta U^4F/373669453125        mod J_BC,
disc_Phi(F)=16Delta U^2V^2J_BC.                        (26)
```

On `J_BC=0`, hypotheses `(21)` imply `U!=0`, and the quadratic leading
coefficient `-1215eta U^2` is a unit. Therefore `(25)` is exactly

```text
F=-1215eta U^2(Phi-Phi_*)^2               mod J_BC,
Phi_*=-2Delta V/(1215eta).                             (27)
```

Thus the actual degree is sixteen on `J_BC=0,Phi!=Phi_*`.

On `J_BC=0,Phi=Phi_*`, define

```text
H(Delta)=187425Delta^3-25939920Delta^2
         +1215936512Delta-18687983616,                 (28)

S=397166535676406250000.
```

Exact reduction gives

```text
[T^15]Q_B=-3Delta^3U^6H(Delta)^2/(S eta^3).           (29)
```

Hence the actual degree is fifteen when `H(Delta)!=0`. On `H(Delta)=0`, the
next coefficient is

```text
[T^14]Q_B=-Delta^3U^4N(Delta)/(S eta^3),              (30)
```

where

```text
N(Delta)=
 1161862920421875Delta^8
-436262053412193750Delta^7
+72428665164732450000Delta^6
-6925998948093233280000Delta^5
+416163890453413588992000Delta^4
-16047146699676377898024960Delta^3
+386666237774835302050824192Delta^2
-5306381690338021964332924928Delta
+31649810610496845164669042688.                        (31)
```

The exact terminal coprimality gate is

```text
gcd(H, Delta*U*V*N)=1.                                (32)
```

Thus every root of `H` is an admissible unit for the prefactors in `(30)`,
`Phi_*!=0` there, and `(30)` is nonzero. There is no fifth stratum. The
complete row-B degree/length tower is therefore

| exact coefficient stratum | `deg Q_B` | `L=deg Q_B+2+2` |
|---|---:|---:|
| `J_BC!=0` | `17` | `21` |
| `J_BC=0, Phi!=Phi_*` | `16` | `20` |
| `J_BC=0, Phi=Phi_*, H!=0` | `15` | `19` |
| `J_BC=0, Phi=Phi_*, H=0` | `14` | `18` |

The first `+2` is the divisor `(6T+1)^2`, counted with multiplicity and made
reduced by the lemma; the second `+2` is the restored pair `(14)`. Since
`Q_B(0)!=0`, no residual point is hidden in the discarded `T=0` Sylvester
artifact.

## 6. Complete row-C degree tower

Assume

```text
Ctop=Phi=0,                 eta*Delta*Btop!=0.          (33)
```

Use `U,V,J_BC` from `(22)`. Exact differentiation gives

```text
deg_X(f,h)=(6,7),
[X^6]f=7eta T^7,
[X^7]h=8eta T^7.                                      (34)
```

Thus no common point at `X=infinity` occurs for any finite `T!=0`. The exact
resultant is

```text
Res_X(f,h)=T^32(6T+1)^2Q_C(T),
Q_C(0)=2*3^6 eta Btop^5,
[T^15]Q_C=eta^3J_BC^2/22500.                          (35)
```

Modulo `J_BC`,

```text
[T^14]Q_C=
 4Delta^2eta U^4V^2/8303765625.                       (36)
```

On `J_BC=0`, hypotheses `(33)` imply `U!=0`; otherwise
`J_BC=54675eta^2!=0`. Thus `(36)` vanishes only when

```text
V=0,                      Delta=22784/825.             (37)
```

This value is not the row-D value `2048/45`, so `Btop!=0` remains valid.
Equation `J_BC=0` then forces

```text
eta^2=-739213574144/187171875!=0.                     (38)
```

At the two complex points `(37)--(38)`, exact reduction gives

```text
[T^13]Q_C=
 -39082296781894638834860007697799970816
  --------------------------------------------------- eta !=0.
              8994856609344482421875                  (39)
```

There is no further degree drop. The complete row-C tower is

| exact coefficient stratum | `deg Q_C` | `L=deg Q_C+2+2` |
|---|---:|---:|
| `J_BC!=0` | `15` | `19` |
| `J_BC=0, V!=0` | `14` | `18` |
| `J_BC=0, V=0` | `13` | `17` |

Again `Q_C(0)!=0`, and neither a residual discriminant nor
`Q_C(-1/6)!=0` is used.

## 7. Complete packets, carrier, and uniform contradictions

THM-4157's strict transforms depend only on the structural units in `(5)`.
The new `J_BC`, `Phi_*`, and `H` walls do not change them. Finiteness of the
critical scheme follows from `(24)` or `(35)` and their nonzero constants,
together with the finite `T=0` row `(13)--(14)`. THM-3827 then supplies
geometric connectedness, and THM-4103 plus Riemann--Hurwitz makes the
displayed packet complete.

The row data are

| row | complete packet | genus | full defect | finite `(n,beta)` | finite origin index | full `n` |
|---|---|---:|---:|---:|---:|---:|
| `A` | `(7,7,4,2,2,2,1)` | `10` | `18` | `(19,3)` | `15` | `25` |
| `B` | `(6,6,4,2,2,2,1)` | `9` | `16` | `(17,3)` | `13` | `23` |
| `C` | `(5,5,4,2,2,2,1)` | `8` | `14` | `(15,3)` | `11` | `21` |
| `D` | `(7,4,2,2,2,1)` | `7` | `12` | `(12,3)` | `9` | `18` |

The nonrational place in every row is

```text
q-1/2=K W^2-eta W^3.                                  (40)
```

Because `eta!=0`, `(40)` is one prime separable cubic closed point over
`C(q)`, including when `K=0`. It either responds at the origin as a whole or
splits after base change into three conjugate carrier meridians, each a
transposition. Thus `m=beta=3`; MISTAKE-509 forbids selecting an arbitrary
geometric subset.

Let `X,Y` be the handle permutations and `L=r_0+r_1` the affine Morse inverse
count. Fixed-sheet transport gives

```text
|supp X|+|supp Y|<=2n-L.                              (41)
```

The carrier-orbit lemma of THM-4171 applies in every finite row because
`n>m+1=4`. With

```text
kappa=|supp X intersect supp Y|,
```

it gives

```text
kappa<=n+m-L,
ind(mu_O)<=2kappa+m<=2(n+m-L)+m.                      (42)
```

For the full response, `X,Y` themselves generate transitively, so

```text
kappa<=n-L,
ind(mu_O)=ind([X,Y])<=2(n-L).                         (43)
```

Use the smallest length in each complete degree tower. Larger lengths only
strengthen the bounds:

| row | `L_min` on length-controlled strata | finite ceiling from `(42)` | required finite origin index | full ceiling from `(43)` | required full defect |
|---|---:|---:|---:|---:|---:|
| `A` | `19` | `2(19+3-19)+3=9` | `15` | `2(25-19)=12` | `18` |
| `B` | `18` | `2(17+3-18)+3=7` | `13` | `2(23-18)=10` | `16` |
| `C` | `17` | `2(15+3-17)+3=5` | `11` | `2(21-17)=8` | `14` |
| `D` | `15` | `2(12+3-15)+3=3` | `9` | `2(18-15)=6` | `12` |

Every ceiling is strictly smaller than the required origin index. Thus no
finite or full response can be transitive on any length-controlled stratum.
The row-D `H_30=0` locus was already excluded locally by `(8),(12)` before
this monodromy comparison. Rows `(5)` are exhaustive, so the theorem follows.
**QED.**

## 8. Exact scope and firewalls

This theorem closes precisely the complete coefficient locus `(4)` inside
the inherited exact-weight-nine live reduced `(2,3)` seam. It closes every
row-`A,B,C,D` coefficient, leading-degree, projected-root-collision, and
universal-fibre-collision wall compatible with `eta*Delta!=0`.

It does **not** assert any of the following:

1. `eta=0`, where the repeated-top degree-three carrier and exact top row
   change;
2. `Delta=0`, where the common Newton polygon contracts;
3. entry into the reduced seam or a different reduced cell;
4. residual weight at least ten;
5. validity of the finite-separable-carrier transport outside the inherited
   THM-4147/4157 contract;
6. an arithmetic period-three orbit associated to the cubic closed point;
7. `JC(2)` or `DC(2)`.

The source projections `Res_s(Acrit,Ccrit)` used in THM-4157 are not
hypotheses in the new B/C argument. Their endpoints and discriminants may
vanish. The B/C proof uses only the actual specialized `T`-resultant divisor,
the finite-`T` `X`-infinity sidecar, the Keller--Morse bridge, and the complete
packet/carrier contract.

## 9. Exact artifacts and replay

The primary exact SymPy certificate directly reconstructs `(3)`, proves the
row-independent bridge `(12)`, checks both B/C finite-`T` `X`-infinity
sidecars, computes both symbolic resultants, verifies the encoded coefficient
tower identities `(22)--(39)`, and checks the terminal coprimality gate. Its
degree-to-length and four worst monodromy checks use the universal `T=0` pair
`(13)`, whose exact symbolic proof and frozen replay are inherited from
THM-4157/4173; the primary certificate does not directly recompute `(13)`.

```text
04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.py
raw-LF sha256:
928d77ae3a9815abaa38cf7ca6e29716f5d3f1575fc46e68ec93481baf59cc98

05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.out
raw-LF sha256:
7b32ef45ab66500a181215a1c14d6beb4d03d1af8fed2569fd63c30cdd199d79
```

The independent PARI/GP certificate is deliberately narrower. It independently
reconstructs only the B/C resultant towers and verifies the common leading
wall, row-B double-root and terminal tower, row-C terminal specialization, and
`gcd(H,Delta*U*V*N)=1` (script variables `gcd(H,D*A*b*N)=1`). It does not
independently audit the Hessian bridge, `T=0` pair, `X`-infinity sidecars,
packets, carrier transport, or monodromy.

```text
04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.gp
raw-LF sha256:
37d5f1da666e2ed023987d9e0ab8cbd376bb5e53dfd9fec40422ad0c65a0b734

05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.out
raw-LF sha256:
a6edc0762e692e054dea56553cb60edbc2f587fdfd0852655d2e4b4b4d39eb59
```

Replay the canonical artifacts against their frozen outputs with

```bash
python3 -B \
  04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.py \
  | diff -u \
      05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.out -
python3 -B -O \
  04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.py \
  | diff -u \
      05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.out -
PYTHONHASHSEED=4176 python3 -B \
  04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.py \
  | diff -u \
      05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.out -
gp -fq \
  04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.gp \
  | diff -u \
      05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.out -
```

The first three commands byte-match the same primary output, including
`checks=26` and `verdict=B_AND_C_CLOSE`. The GP command byte-matches the
independent output `PARI_GP_INDEPENDENT_BC_BRIDGE_ACCEPT`. Verify the four
raw-LF bindings with

```bash
sha256sum \
  04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.py \
  05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176.out \
  04-computation/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.gp \
  05-knowledge/results/jc23_repeated_top_rows_bc_complete_exclusion_thm4176_independent.out
```

The inherited row-A and row-D proof and replay bundles remain bound by
THM-4173 and THM-4157 respectively; the universal `T=0` identity used by the
primary length audit is bound there rather than duplicated here.
