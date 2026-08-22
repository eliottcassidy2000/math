---
id: THM-3289
title: "Affine transverse C0-E0 coupled clutch critical no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For either THM-3212 cubic accessory response pair, every simultaneous
  affine C_0,E_0 clutch with B=1 leaves a critical point.  The genuinely
  coupled lane has a degree-53 saturated resultant; the finite clutch
  controls all T contact, and exact cubic-field PRS calculations cap the
  remaining S contact at three, leaving at least 50 units of off-owner
  critical-resultant root multiplicity.  This is multiplicity, not a claim
  of 50 distinct critical points.
source: root/creative-synthesis-recover/2026-08-03
audit: >
  The primary exact scout uses a literal 6-by-6 Sylvester determinant and
  direct number-field Euclidean algorithms.  The independent audit does not
  import or execute it: it rebuilds both response pairs, uses SymPy's builtin
  resultant, derives the generic S jets with coefficient-list arithmetic,
  and recomputes the final wall ideals with its own rational-triple cubic
  arithmetic and dense PRS.  Both routes verify the 40-term quotient, unique
  degree-97 infinity term, four finite-gate T rows, q3/q4 hostile walls,
  terminal-unit q5/q6 ideals, rad(V)=ST, disjoint degree-44 boundary,
  degree-53 residual, and squarefree owner-disjoint controls.  Normal,
  optimized, and stored transcripts agree; both sources have zero assertion
  nodes and zero floating literals.
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
  - THM-3279-affine-transverse-clutch-critical-no-go
related:
  - THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape
script: 04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.py
output: 05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.out
script_sha256: f63ff06e3f5ed30f3f6bc5be99756c347d6af5f8e9b220ce8336abff2cd2ca31
output_sha256: 1aef4341650cdfaf1043a8699e3a1725a0100af6d9848d99dfa924b6f054dba1
audit_script: 04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.py
audit_output: 05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.out
audit_script_sha256: b2fa8c96854549ccb9e515485214c119b685b31456fb7c53c5e2bd83f7933831
audit_output_sha256: 48d50289c98c9dd17e099497d21cd9648cac27a097b339a1f75a1e13d8fd8837
hash_basis: LF-normalized bytes
---

# THM-3289 -- affine transverse C0-E0 coupled clutch critical no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let `K_i` be either cubic accessory field of
[THM-3212](THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch.md),
embedded in an algebraically closed characteristic-zero field `K_0`.  Retain
the response pair and owner divisor

```text
V=4SDT^2/Gamma^2,       A=2SET/Gamma,
g=ST,                   2VA'-AV'=2V,
rad(V)=g,
deg V=16,               deg A=8.                         (1)
```

For arbitrary `c_0,d,e_0,k in K_0`, put

```text
C=c_0+d x,              E=e_0+k x                       (2)
```

put

```text
P_(C,E)(x,z)=(V(x)z^2+z+C(x))^2+A(x)z+E(x).             (3)
```

Then every polynomial `(3)` has a critical point.  The constant `e_0` is
gradient-inert, and the proof uses the following disjoint parameter lanes:

1. if `d=0`, this is the complete constant-clutch lane of THM-3212;
2. if `d!=0` and `k=0`, this is THM-3279; and
3. if `dk!=0`, either the finite clutch fails and gives an explicit critical
   point over `g`, or at least `50` units of saturated critical-resultant
   root multiplicity remain away from `g`.

At any critical point of `(3)`, `Jac(P_(C,E),Q)=0` for every polynomial `Q`.
Thus `(3)` cannot be one coordinate of a polynomial Keller pair.  The theorem
constructs no mate or inverse cover and proves no case of the planar Jacobian
conjecture.

## 2. Exact finite clutch

At a root `alpha` of `g`, THM-3212 gives

```text
V(alpha)=A(alpha)=0,       A'(alpha)!=0.                 (4)
```

There

```text
P_z=2(z+C(alpha)),
```

whose unique zero is `z=-C(alpha)`.  At this point the square contributes
nothing to `P_x`, and

```text
P_x(alpha,-C(alpha))=k-A'(alpha)C(alpha).                (5)
```

Consequently the finite clutch passes exactly when

```text
Delta=k-A'C                                                (6)
```

is a unit modulo `g`.  If it is not, `(alpha,-C(alpha))` is an explicit
critical point over any common root.  The rest of the proof assumes that `(6)`
is a unit.

## 3. Localized gradient pair and the 40-term quotient

On `V!=0`, set

```text
y=Vz,                  L=y^2+y+CV.                       (7)
```

After multiplying by powers of `V`, the first gradient equation is

```text
R_1=2L(2y+1)+VA.                                          (8)
```

The second reduces, using `(1)`, to

```text
R_2=V^3k+V^2y+L(-V'y+2V^2d).                             (9)
```

The exact sign ledger is

```text
raw-R_2=(V'y/2)R_1.                                      (10)
```

A literal `6 x 6` Sylvester determinant gives

```text
Res_y(R_1,R_2)=V^3 K_(C,k).                              (11)
```

The quotient has `40` monomials.  If `K_(C,0)` denotes the exact 20-term
factor in THM-3279, then the coupled correction is

```text
K_(C,k)=K_(C,0)
 +k(12A^2V^2(V')^2+64ACV^4V'd+16ACV^2(V')^2
     +96AV^5d-16AV^3V'd-24AV^3V'-4AV(V')^2
     -128CV^5d+64CV^5-32CV^3V'
     +32V^4d+32V^4+8V^2V')
 +k^2(48AV^4V'+128CV^6d+32CV^4V'
      -32V^5d-96V^5-8V^3V')
 +64V^6k^3.                                              (12)
```

For affine `C`, the unique degree-`97` monomial in `(12)` is

```text
128 C V^6 d k^2.                                        (13)
```

Since the leading coefficient of `C` is `d`, the leading coefficient of
`K_(C,k)` in the live lane `dk!=0` is

```text
128 d^2 k^2 V_lead^6 !=0.                               (14)
```

Thus

```text
deg K_(C,k)=97.                                          (15)
```

## 4. The four T rows are exactly the finite gate

At a root of `T`, write

```text
V=v t^m+...,        A=a t+...,       C=c+dt,
a=2/(2-m),          m in {3,4,5,6}.                     (16)
```

Direct substitution in `(12)` gives

```text
ord_t K_(C,k)>=3m-1,
[t^(3m-1)]K_(C,k)
 =16m(m-1)/(m-2) v^3 (k-ac).                            (17)
```

Here `a=A'(alpha)`, so the last factor is exactly `(6)` at this root.  In the
finite-gate lane it is nonzero.  The four orders are `8,11,14,17`, precisely
those encoded by

```text
boundary_(4111)=S^3T^8x^9,
boundary_(3211)=S^3T^8x^6(x-1)^3.                       (18)
```

The four owner factors are pairwise disjoint.  Equation `(17)` gives the
stated local order at every `T` root, while direct expansion at the simple
`S` root makes the first three coefficients vanish universally.  These local
orders, rather than the later rational control, prove divisibility by the
degree-`44` boundary in `(18)`.  Hence

```text
K_(C,k)=boundary_i H_(C,k),        deg H_(C,k)=53.       (19)
```

Equation `(17)` also proves that `H_(C,k)` has no root on `T` after the
finite clutch passes.

## 5. The two live S walls

At the simple root of `S`, use `t=S` and write

```text
V=v_1t+v_2t^2+v_3t^3+v_4t^4+v_5t^5+...,
C=c+dt.                                                   (20)
```

The response equation in `(1)` recursively determines the jet of `A`.
Truncated exact series arithmetic first gives

```text
q_3=[t^3]K_(C,k)
 =(8/3)v_1^2(2c-k)(6cv_1^2+3kv_1^2+4v_2).               (21)
```

Since `A'(S=0)=2`, the branch `k=2c` is exactly finite-clutch failure and is
excluded.  The only live `q_3` wall is therefore

```text
k=-2c-4v_2/(3v_1^2).                                    (22)
```

On `(22)`, the next coefficient is

```text
q_4=64(3cv_1^2+v_2)/(45v_1)
 (30c^2v_1^4+10cv_1^2v_2+30dv_1^3
  -15v_1^3+18v_1v_3-16v_2^2).                          (23)
```

The apparent exceptional slope is not live.  If

```text
3cv_1^2+v_2=0,                                          (24)
```

then `(22)` gives `k=2c`, again violating the finite clutch.  Otherwise
`q_4=0` has the unique solution

```text
d=(15v_1^3-18v_1v_3+16v_2^2
   -30c^2v_1^4-10cv_1^2v_2)/(30v_1^3).                 (25)
```

This closes every exceptional-slope case before any PRS calculation.

## 6. The q5/q6 cubic-field obstruction

After imposing `(22)` and `(25)`, define

```text
U=3cv_1^2+v_2,

Q_5=
 -210c^2v_1^4v_2+630cv_1^5+504cv_1^3v_3
 -518cv_1^2v_2^2+105v_1^3v_2-360v_1^2v_4
 +750v_1v_2v_3-400v_2^3.                               (26)
```

Then

```text
q_5=-32 U Q_5/(315v_1^2),                               (27)
q_6=-32 Q_6/(4725v_1^3),                                (28)
```

where

```text
Q_6=
 12600c^4v_1^8v_2
 +37800c^3v_1^9-37800c^3v_1^7v_3-13650c^3v_1^6v_2^2
 +189000c^2v_1^7v_2+48600c^2v_1^6v_4
 +32310c^2v_1^5v_2v_3-100420c^2v_1^4v_2^3
 +18900cv_1^8-9450cv_1^6v_3+97125cv_1^5v_2^2
 -21000cv_1^5v_5-28800cv_1^4v_2v_4-504cv_1^4v_3^2
 +172056cv_1^3v_2^2v_3-124946cv_1^2v_2^4
 +1575v_1^6v_2-16200v_1^5v_4+25875v_1^4v_2v_3
 -1025v_1^3v_2^3-7000v_1^3v_2v_5+9720v_1^3v_3v_4
 -23640v_1^2v_2^2v_4-15882v_1^2v_2v_3^2
 +73098v_1v_2^3v_3-37168v_2^5.                         (29)
```

As polynomials in `c`, `UQ_5` and `Q_6` have degrees `3` and `4`.  Exact
Euclidean division in both cubic accessory fields gives

```text
(4,3) -> (3,2) -> (2,1) -> (1,0),
gcd(UQ_5,Q_6)=1.                                        (30)
```

After removing the already excluded slope factor `U`, the sharper live
profile is

```text
(4,2) -> (2,1) -> (1,0),
gcd(Q_5,Q_6)=1.                                         (31)
```

The characteristic-zero monic-pair digests for `(UQ_5,Q_6)` are

```text
4111  92c368a8101e6ccebc69c4961bda10b41a85dd3f367bdbfeb7ee04747421c319
3211  f7268035c556bfdc85afd670405986ddbe6839d2ddfb3a13900cba2baabc901a
```

Thus `q_5` and `q_6` cannot vanish simultaneously over any extension of
either accessory field.  It follows that

```text
ord_S K_(C,k)<=6,             ord_S H_(C,k)<=3.          (32)
```

Only the generated ideal and its unit gcd are load-bearing here.  Raw PRS
representatives and intermediate normalizations may differ by field units;
the displayed degree profiles record these two exact routes and are not
treated as canonical coefficient lists (MISTAKE-360).

## 7. Residual count and hostile controls

Equations `(17)`, `(19)`, and `(32)` leave at least

```text
53-3=50                                                   (33)
```

units of resultant-root multiplicity away from `g`.  This is multiplicity,
not a count of `50` distinct roots or critical points.  Homogenizing in `y`,
the leading `y`-coefficient `4` of `R_1` excludes a common projective root at
`y=infinity`.  Consequently every distinct off-owner root supporting this
multiplicity yields an affine common zero of `(8),(9)`, hence a critical point
of `(3)`.

The rational coupled control

```text
C=1+x,                   E'=1                            (34)
```

passes the finite gate in both accessory fields.  Exact global division gives

```text
(deg K,deg boundary,deg H)=(97,44,53),
gcd(H,g)=1,                 gcd(H,H')=1.                 (35)
```

Its monic residual digests are

```text
4111  ae98818930edc795b9a7a73a957ddf2f2e70dee40695bb3ea489c9a8c16780d4
3211  621f57adfa31900c88fa44d0443d9baf81d199a201a4b9963d99914d134b1d21
```

These controls show that the degree-53 invoice and boundary separation are
attained; the universal boundary divisibility comes from the disjoint local
orders above, not from the controls.

## 8. Loss ledger and audit record

The inheritance map is

```text
(constant C lane) union (constant E lane)
       --> simultaneous affine (C,E) lane.              (36)
```

It preserves the response identity `(1)`, the owner divisor `g`, the exact
finite clutch `(6)`, and the localized gradient interpretation.  It changes
the infinity degree from `96` to `97` and adds one possible unit of residual
`S` contact.  The needed sidecar is the pair of univariate S-wall
polynomials `(UQ_5,Q_6)`; their cheapest decisive test is the exact cubic-field
PRS `(30)`.

The primary companion uses exact rational and number-field arithmetic, has no
floating literals or Python `assert` nodes, and agrees in normal and optimized
modes with its stored transcript.  The independent audit rebuilds the response
pairs and uses a builtin resultant rather than the primary literal Sylvester
determinant; coefficient-list jets rather than the primary symbolic series;
and independent rational-triple cubic arithmetic and dense PRS rather than the
primary field-polynomial Euclidean route.  It separately verifies the finite
hostile, exceptional-slope hostile, `rad(V)=ST`, disjoint local boundary
orders, two wall ideals, and squarefree global controls.
