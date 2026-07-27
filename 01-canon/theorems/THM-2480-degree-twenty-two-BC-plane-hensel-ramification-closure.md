---
id: THM-2480
title: "Degree-twenty-two B-C plane Hensel-ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane D=E=W=0 is empty. For B,C nonzero, the
  root-free invariants lambda=C^2/B^3 and r=C/(By) give a 32-term plane
  curve of bidegree (8,5). Its Newton quadrilateral has only two line
  and two quadratic small Minkowski summands. The common line ansatz has
  unit coefficient ideal. For a quadratic factor, the squarefree fixed
  L_5 section kills the first Hensel deformation, and the six scaled
  equations through order three have unit ideal. Thus every physical
  fibre is absolutely irreducible. The v-discriminant is
  constant*lambda^7*W_5^2*K_30, where W_5 is exactly the moving
  first-flux wall. The complete exceptional divisor includes one
  exponent-four collision factor, but it is wall-disjoint; the only
  collision/wall overlap has exponent three. A Sylvester-valuation bound
  and ramification parity force normalization genus at least six in
  every fibre. No rational trajectory survives. Together with the nine
  earlier mixed-plane closures and the axes, this closes all support-at-
  most-two coefficient strata in the inherited chart. It does not close
  higher mixed strata, JC(2), or DC(2).
source: codex-2026-07-27-degree-twenty-two-BC-hensel-ramification
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
  - THM-2428-degree-twenty-two-B-axis-trigonal-ramification-closure
related:
  - THM-2429-degree-twenty-two-CW-plane-hyperelliptic-family-closure
  - THM-2437-degree-twenty-two-DW-plane-quartic-ramification-closure
  - THM-2463-degree-twenty-two-BD-plane-square-lift-closure
  - THM-2468-degree-twenty-two-BW-plane-square-lift-closure
  - THM-2469-degree-twenty-two-CD-plane-coprime-weight-ramification-closure
  - THM-2470-degree-twenty-two-EW-plane-coprime-weight-ramification-closure
  - THM-2472-degree-twenty-two-DE-plane-coprime-weight-ramification-closure
  - THM-2475-degree-twenty-two-CE-plane-hensel-ramification-closure
  - THM-2476-degree-twenty-two-BE-plane-hensel-ramification-closure
script: 04-computation/jc2_degree22_bc_plane_hensel_ramification_thm2480.py
output: 05-knowledge/results/jc2_degree22_bc_plane_hensel_ramification_thm2480.out
script_sha256: 5281789eb1ae9c42301ee12d7565261b14aa5be898074924ac86a73d1ad0a1c3
output_sha256: d856fde91e049fc88e058edf134af45258cd4483d99c60adfedc21055415c529
hash_basis: working-tree bytes (LF)
---

# THM-2480 -- the degree-twenty-two B-C plane is empty

**PROVED + VERIFIED-EXACT.**

This is the last of the ten support-two coefficient planes.  Both active
coefficients move the first-flux denominator, the eliminant has a new Newton
quadrilateral, and the ramification collision divisor has one exponent-four
component.  The fixed-section Hensel operation still proves irreducibility,
and exact wall separation keeps the exponent-four collision harmless.  The
exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
D=E=W=0
    => contradiction.                                           (1)
```

Thus the complete `B,C` support-two coefficient plane is empty.

## 1. Coprime-weight invariant coordinates

Use the target-translated coordinates of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `B=0` and `C=0` are the closed `C` and `B` axes of THM-2425
and THM-2428.  Suppose

```text
B,C!=0.                                                (3)
```

First take `y!=0` in `C(x)` and define

```text
lambda=C^2/B^3 in C*,             r=C/(By).           (4)
```

Then

```text
B/y^2=r^2/lambda,                  C/y^3=r^3/lambda.  (5)
```

Put

```text
v=u/y^2,                           zeta=Z/y^3.         (6)
```

Multiplying the first two normalized fluxes by `lambda` gives

```text
F_1
 =1331[616r^2+lambda(63-1089v)]zeta
   +4[(-745360v+6160)r^2+(2342560v-58080)r^3
      +lambda(922383v^2-25410v+63)]
 =0,                                                        (7)

F_2
 =r^2[65591680zeta+1443016960v^2-71554560v+98560]
  +r^3[-206145280zeta+449771520v-1239040]
  +lambda[15944049zeta^2+(-162339408v+2236080)zeta
           -1190488992v^3+147581280v^2-1219680v+672]
 =0.                                                        (8)
```

The open first-flux chart is

```text
616r^2+lambda(63-1089v)!=0.                         (9)
```

Equation (7) reconstructs `zeta` uniquely.

## 2. The 32-term eliminant

Exact elimination gives

```text
Res_zeta(F_1,F_2)=28344976 P_lambda(r,v),              (10)
```

where

```text
P_lambda
 =(55873616691200v-1385296281600)r^8
  +(-24889156526080v+558316380160)r^7

  +[-49388286182400lambda v^2+5714347161600lambda v
     -111318451200lambda+34222590223360v^2
     +3959638538240v-44411530240]r^6

  +lambda(59714927838720v^2-3653718942720v
           +64184440320)r^5

  +lambda(-149234938081152v^3-9500102156160v^2
           +695599766400v-6017413248)r^4

  +lambda^2(-44449457564160v^3+4836786704640v^2
             -113342423040v+1317254400)r^3

  +lambda^2(206782580709936v^4+6246495741024v^3
             -1509756494400v^2+34466937120v
             -193496688)r^2

  -567lambda^3 L_5(v),                                  (11)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (12)
```

The fixed section and initial gap are

```text
P_lambda(0,v)=-567lambda^3L_5(v),

[r]P_lambda=0.                                        (13)
```

The quintic `L_5` is squarefree.

## 3. The Newton quadrilateral and all factor shapes

The Newton polygon is

```text
N=conv{(0,0),(8,0),(8,1),(0,5)}.                     (14)
```

Its primitive directed edges have lengths `(8,1,4,5)`.  A lattice
Minkowski summand has edge lengths

```text
(a,b,c,d)=(2c,b,c,b+c),

0<=b<=1,                    0<=c<=4.                 (15)
```

The nonzero constant and `v^5` coefficients normalize every factor to be
monic in `v` and remove translations.  Since factor `v`-degrees add to five,
one factor has degree one or two.  The complete small inventory is

```text
degree one: (b,c)=(1,0),(0,1),

degree two: (b,c)=(1,1),(0,2).                        (16)
```

The two line supports are both contained in

```text
v+ar^2+br+c.                                           (17)
```

Exact division gives

```text
I_1=(1) in Q[a,b,c,lambda].                            (18)
```

For degree two, `(0,2)` is the triangle `2 Delta`, while `(1,1)` is its
quadrilateral sub-support.  It therefore suffices to exclude the larger
triangle.

## 4. The scaled Hensel obstruction

Let a proposed monic quadratic factor have fixed section

```text
q_0=v^2+cv+h.                                         (19)
```

It divides the squarefree `L_5`, so it is coprime to its cubic cofactor.
The missing order one in (13) and Hensel uniqueness force

```text
q_1=0.                                                (20)
```

The maximal triangular support permits

```text
q_2=av+f,                       q_3=e.                (21)
```

The quadrilateral support is the specialization `e=0`.  Scale away the
nonzero physical ratio by putting

```text
A=lambda a,             F=lambda f,             E=lambda e. (22)
```

To state the exact finite obstruction compactly, set `lambda=1` in (11),
truncate modulo `r^4`, and call the result `H(r,v)`.  Put

```text
Q(r,v)=v^2+cv+h+(Av+F)r^2+Er^3.                      (23)
```

For `i=0,1` and `j=0,2,3`, define

```text
G_ij=[v^i r^j] rem_v(H,Q).                            (24)
```

Direct coefficient comparison before specialization verifies that the
original equations are `lambda^3G_i0` at order zero and
`lambda^2G_ij` at orders two and three.  Thus no physical parameter was
lost.  Exact Groebner reduction gives

```text
(G_10,G_00,G_12,G_02,G_13,G_03)
 =(1) in Q[A,F,E,c,h].                                (25)
```

No quadratic divisor of the fixed section lifts even through order three.
Equations (18),(25) exhaust (16), so `P_lambda` is absolutely irreducible
for every `lambda!=0`.  Let `C_lambda` be its smooth projective
normalization.  Projection to the `r`-line has degree five.

## 5. Exact wall and ramification divisor

Equality in (9) gives

```text
v=(616r^2+63lambda)/(1089lambda).                     (26)
```

Define

```text
W_5(r,lambda)
 =745360r^5-71148r^4+43560lambda r^3
   +5082lambda r^2+945lambda^2.                       (27)
```

Direct substitution yields the exact rational wall identity

```text
P_lambda(r,(616r^2+63lambda)/(1089lambda))
 =256W_5(r,lambda)^2/(9lambda).                       (28)
```

The squared factor is therefore excluded wall mass.  The quintic
discriminant is

```text
Disc_v(P_lambda)
 =c lambda^7 W_5(r,lambda)^2 K_30(r,lambda),          (29)
```

for `c in Q*`.  Normalize `K_30` to be primitive with

```text
K_30(0,lambda)=38724883599623115234375lambda^13.      (30)
```

It has degree thirty in `r`, degree thirteen in `lambda`, and 87 nonzero
terms.  Its leading coefficient is

```text
[r^30]K_30
 =927499340872432064135168000000 L_0(lambda),         (31)

L_0
 =49121386875lambda^3+35843727150lambda^2
   +34155999339lambda+1686616064.                     (32)
```

Put

```text
L_1=1782lambda+245,                 L_2=22275lambda+2744,

A_5
 =134711704729495154443359375lambda^5
  +1800011771181150399288281250lambda^4
  +6418659057813101346623578125lambda^3
  +1620235334410707585796706250lambda^2
  +103287017473928512574894775lambda
  -384735691349720065722248.                          (33)
```

Exact factorization gives

```text
Disc_r(K_30)
 =c' lambda^325 L_1^2 L_2^4 A_5^3
      R_18 P_20^2 Q_20^3,                            (34)
```

where `R_18,P_20,Q_20` are the unique primitive positive-leading
irreducible factors of degrees `18,20,20` occurring with exponents `1,2,3`.

The wall collision divisor is

```text
Res_r(K_30,W_5)
 =c'' lambda^52 A_5^3 S_8(lambda),                    (35)
```

where

```text
S_8
 =12337829835953864392149281489742232500000000000000lambda^8
  -44369823212208566246105819275522866775312500000000lambda^7
  +16885864307881167661079424840506375193448798828125lambda^6
  -21179519183140752458837472696065162224372050000000lambda^5
  -9294048284435692495345979933901524446649735250000lambda^4
  +131589285074186063043261783346562971089053066250lambda^3
  +205953213712224149574057621432756223317422854200lambda^2
  +942027416127503933058469688031656737095805128lambda
  -1311283001340114972135581249229837964617938539.   (36)
```

The eight polynomials

```text
L_0,L_1,L_2,A_5,R_18,P_20,Q_20,S_8                  (37)
```

are squarefree and pairwise coprime.  In particular, the exponent-four
collision factor `L_2` is wall-disjoint.  The sole collision/wall overlap is
the displayed exponent-three factor `A_5`.

At each of the three roots of `L_0`, the `r^29` coefficient, the specialized
degree-30 discriminant, and the wall resultant are nonzero.  The
leading-degree specialization identity therefore makes all three degree-29
fibres squarefree and wall-disjoint.

## 6. Uniform ramification floor

The Sylvester valuation lemma says that discriminant valuation `e` bounds

```text
deg gcd(k_0,partial_rk_0)<=e.                         (38)
```

At an `L_1` root, at least `30-2*2=26` roots of `K_30` are simple; at an
`L_2` root at least `30-2*4=22` are simple.  Neither parameter meets the
wall.  The factors `R_18,P_20,Q_20` are also wall-disjoint and have exponent
at most three.

At an `A_5` root, at least `30-2*3=24` roots are simple and at most the five
roots of `W_5` must be discarded, leaving at least nineteen.  At an `S_8`
root, `K_30` is squarefree and at most five roots are lost, leaving 25.  The
degree-drop and generic fibres are squarefree and wall-disjoint.  These cases
exhaust every `lambda in C*`, so uniformly

```text
at least 19 simple finite roots of K_30 lie off W_5.   (39)
```

Each contributes one unit of tame ramification: (29) has valuation one, and
normalization changes order discriminant valuation by twice the index.

For the degree-five projection,

```text
2g(C_lambda)-2=-10+R.                                 (40)
```

Total ramification `R` is even.  The nineteen visible units force

```text
R>=20,                    g(C_lambda)>=6.             (41)
```

## 7. Genus and trajectory closure

The rational pair `(r,v)` from (4),(6) would give a nonconstant rational map
from `P^1` to the positive-genus normalization, impossible.  Hence `r,v`
are constant.  Since `r=C/(By)` is nonzero, `y` and then `u` are constant.
Equation (7) reconstructs constant `zeta`, so `Z,T,q` are constants,
contradicting the genuine deck.

It remains to treat `y=0`.  Put

```text
D_0=1331(616B-1089u),              n=-9370240Cu.      (42)
```

The open chart says `D_0!=0`, and the first flux gives

```text
D_0Z-n=0,                         Z=n/D_0.            (43)
```

The second flux is

```text
15944049Z^2-206145280CZ
 +1443016960Bu^2-1190488992u^3=0.                    (44)
```

After (43) and multiplication by `D_0^2`, (44) is a nonzero
constant-field polynomial of degree five in `u`; its leading coefficient is
nonzero.  Thus `u`, and then `Z,T,q`, are constants, giving the same deck
contradiction.  If `u=0`, (43) instead gives `Z=0`, directly contradicting
`Z=T^2` with `T=q^2!=0`.  Together with the two axes, this proves (1).

## 8. Complete support-two closure and scope

The ten mixed planes are now closed by

```text
(C,W)  THM-2429,       (D,W)  THM-2437,
(B,D)  THM-2463,       (B,W)  THM-2468,
(C,D)  THM-2469,       (E,W)  THM-2470,
(D,E)  THM-2472,       (C,E)  THM-2475,
(B,E)  THM-2476,       (B,C)  THM-2480.              (45)
```

Together with THM-2423/2425/2428 on the axes, this proves:

```text
in the inherited genuine degree-22 branch with mathcal A!=0,
support(B,C,D,E,W)<=2
    => no trajectory.                                 (46)
```

The next frontier is therefore genuinely higher-dimensional: coefficient
support at least three.  The reusable operation is now clear.  Form a
root-free Bezout coordinate, classify the Newton summands, lift fixed-section
factors only to the first allowed sparse order, identify the entire wall
square exactly, and use the exceptional collision divisor rather than a
generic-only genus count.

Branches outside the inherited reduction, split/even short edges, support
at least three, and integral order raising remain open.  Nothing here proves
`JC(2)` or `DC(2)`.

## 9. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_bc_plane_hensel_ramification_thm2480.py
python3 -O 04-computation/jc2_degree22_bc_plane_hensel_ramification_thm2480.py
```

The companion reconstructs (7)--(13), the Newton quadrilateral and all small
summands, the line and scaled low-order Hensel ideals, the rational wall
identity, the complete discriminant and wall divisors, every coprimality and
degree-drop control, the stratum-specific branch floors, and the `y=0`
quintic.  Normal and optimized transcripts byte-match the stored output.  All
truth-bearing checks use explicit exceptions and remain active under
optimized Python.  **QED.**
