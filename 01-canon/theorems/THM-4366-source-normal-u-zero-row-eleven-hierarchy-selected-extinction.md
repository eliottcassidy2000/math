---
id: THM-4366
title: "Source-normal U-zero row-eleven hierarchy-selected extinction"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4364 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On the full U=0 source pullback, the Phi=0 boundary
  has genuine row-nine points but a fixed nonzero row-ten bracket residual.
  On Phi!=0, row-ten bracket plus joint P_2/P_3 depth leaves exactly six
  strict source points, three Phi^2 roots times a sign sheet, with affine-
  eight terminal fibres; a row-eleven bracket residual coprime to their
  cubic kills all six. The new opposite-order L_(10,11,3)(A) hierarchy
  consumer is distinct but equals -4/3 times the old C consumer on this
  graph, so it is an alternative certificate, not an extra cut. No all-row
  lift, seam-entry, Keller-pair, JC(2), or DC(2) conclusion is asserted.
source: root + hierarchy-consumer scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
  - THM-4364-sharp-binomial-diagonal-annihilator-hierarchy
related:
  - THM-4357-source-normal-row-eight-endpoint-pullback-stratification
  - THM-4361-source-normal-beta-zero-row-ten-joint-depth-extinction
  - THM-4362-tetrahedral-diagonal-annihilator-and-sharp-depth-threshold
mistake_firewall:
  - MISTAKE-354
  - MISTAKE-522
  - MISTAKE-540
primary_script: 04-computation/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_thm4366.py
primary_output: 05-knowledge/results/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_thm4366.out
primary_script_sha256: 3d942317f566f2937ee876c2384a611d631d72bd5b7d880ce1b3af57fc774027
primary_output_sha256: e9d10af6bac654520cd626c7cda344d53a2f3d0601519d6ca1727c03fb32d118
independent_referee_script: 04-computation/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_independent_referee_thm4366.py
independent_referee_output: 05-knowledge/results/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_independent_referee_thm4366.out
independent_referee_script_sha256: ce32acf13d4744cbdde8ff00fd0a923904e55f69a549ca7805116162ed99a6c3
independent_referee_output_sha256: 26bcaa8877c72e8fa85fda78abee3e3356d35ec3dc75d7bcdfbbbe58b82194f4
hash_basis: raw LF bytes
audit: >
  PASS WITH RESIDUAL-SUPPORT AND RANK WORDING REPAIRS. The primary performs
  117 exact symbolic checks. A 461-assertion import-free referee independently rebuilds
  the source branches, four projected-depth systems, hierarchy consumers,
  strict-gate certificates, and row-eleven obstruction. Normal/-O/hash-
  seeded/frozen LF streams match.
---

# THM-4366 -- Source-normal U-zero row-eleven hierarchy-selected extinction

**PROVED FINITE-ROW RELATIVE TO THM-4308/4315/4364 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THIS CLOSES THE DECLARED `U=0` FINITE SOURCE PULLBACK
THROUGH ROW ELEVEN. IT DOES NOT ASSERT AN ALL-ROW LIFT, POLYNOMIAL
TERMINATION, SEAM ENTRY, A KELLER PAIR, `JC(2)`, OR `DC(2)`.**

## 1. Source pullback and statement

Work over THM-4308's algebraically closed characteristic-zero, source-normal,
residual-weight-at-most-twelve universe. Its fixed coordinates include

```text
Delta=896/15,       Theta=512/75,       K=-32/5,
upsilon_5=-731648/2025,                  zeta_3=-3Phi/2. (1)
```

Impose

```text
U=0  iff  xi_10=237757952/54675.                         (2)
```

Then

```text
W=-13(820125Phi^2+13056802816)/9841500,

Z=-13(32805000Phi^2+36905625Phi eta-1600270057472)
   /265720500.                                           (3)
```

The theorem is:

```text
no point of (2)--(3) passes the row-eleven bracket while
retaining the required row-ten projections
pi_10(A) in pi_10(P_2),       pi_10(C) in pi_10(P_3).     (4)
```

The two `Phi` strata make `(4)` sharp as a staged statement. The `Phi=0`
boundary survives row nine and dies at row ten. The `Phi!=0` branch has six
genuine row-ten bracket-and-joint-depth source points and dies one row later.

The inheritance pass was:

- closest proved mechanism: THM-4315's row-nine source equation and exact
  tangent selection, followed by THM-4361's joint-depth consumer;
- canonical hostile: each depth block can be compatible while their shared
  coordinates are not, and a new functional can be proportional on the
  selected source graph despite being independent in the ambient module;
- corrected near miss: the `Phi=0` row-ten residual support has two entries,
  although one fixed nonzero entry already proves extinction;
- least-used sidecar: THM-4364's opposite-order `A` functional on diagonal
  intercept eleven.

The live board was

```text
U-zero source graph | Phi boundary | row-ten joint cokernel | sign sheet
opposite-order hierarchy consumer | affine terminal fibre | row-eleven gcd. (5)
```

## 2. Literal source rows

The next source rows, retained without specialization before solving, are

```text
G_9 =(20U+10W+4Z)x^6 +(10alpha+6beta)x^7
     +(5upsilon_5+4xi_10)x^8 +(eta+zeta_3)x^9,

G_10=(15U+10W+6Z)x^8 +(5alpha+4beta)x^9
     +(upsilon_5+xi_10)x^10,

G_11=(6U+5W+4Z)x^10 +(alpha+beta)x^11.                  (6)
```

All row and depth conclusions below come from the literal Hasse-truncated
systems attached to `(6)`; no specialization is promoted coefficientwise.

## 3. The Phi-zero boundary is nonempty at row nine

At `Phi=0`, the row-nine Student equation is

```text
E_9=(-11/243)F_0(eta),
F_0(eta)=5646560625eta^2+379697122115584.                 (7)
```

The quadratic is squarefree with nonzero constant, so it has two geometric
roots. The row-nine joint-depth terminal system has rank three, with pivots
`theta_(9,7),theta_(9,8),theta_(9,9)`, and leaves an affine-seven fibre;
`alpha,beta` remain source-free. Thus this boundary is genuinely nonempty at
row nine.

After selecting the seven remaining directions with `G_10`, the possible
nonzero residuals have support `{6,8}`. In particular,

```text
r_8=-11388581281792/332150625 !=0,                       (8)
```

independently of `alpha,beta`; the other residual is

```text
r_6=(321853955625alpha eta-34327712203669504)/215233605000.
```

Hence every point on `(7)` dies before row ten. Taking `alpha=beta=0` at
either root is a positive hostile to premature row-nine extinction: the
source equations and row-nine depth hold, and `W,Z,upsilon_5` and
`alpha^2-4W upsilon_5` are all nonzero.

## 4. The Phi-nonzero row-ten bracket graph

Assume `Phi!=0`. Row nine solves

```text
alpha=(13553385750Phi^2-69677820000Phi eta-5646560625eta^2
       -379697122115584)/(11293121250Phi).               (9)
```

Its literal bracket selector has rank seven on the old row-eight terminal
space. The new row-nine joint-depth system has rank three with the same three
pivots and leaves an affine-seven terminal fibre, with no extra source
equation.

The row-ten bracket selector has rank seven on that fibre. Its two residuals
solve successively for `eta` and `beta`; after substituting `(9)`, the source
graph is

```text
eta=(1752089427968-32805000Phi^2)/(36905625Phi),

alpha=(145629074958046875Phi^4-6583704417821122560000Phi^2
       -26093447590576484415176704)/(23154427662890625Phi^3),

beta=(-319030749548372281640625000Phi^6
      +14659146525631027935375360000000Phi^4
      +101452811911656563438652405841920000Phi^2
      +434321509795518334240224474125496745984)
     /(7690757619746410400390625Phi^5),                  (10)

W=-13(820125Phi^2+13056802816)/9841500,
Z=-225611776/30375.
```

This is a one-parameter row-ten bracket locus before joint depth.

## 5. Joint depth leaves exactly six source points

The exact row-ten projection matrices are

```text
pi_10(P_2): 88 rows, 193 columns, rank 68, left nullity 20;
pi_10(P_3): 99 rows, 304 columns, rank 83, left nullity 16. (11)
```

Each standalone terminal system has coefficient/augmented rank `3/3`, with
pivots `theta_(10,8),theta_(10,9),theta_(10,10)`. The combined coefficient
rank is also three. The two selected solutions agree in coordinates 9 and 10,
while their coordinate-8 difference is

```text
4Q(Phi^2)/(23072272859239231201171875Phi^5),             (12)
```

where

```text
Q(Y)=373891487235896675830078125Y^3
 -15097287707154073014589440000000Y^2
 -101452811911656563438652405841920000Y
 -434321509795518334240224474125496745984.               (13)
```

Thus the combined augmented rank is three on `Q(Phi^2)=0` and four away
from it. Joint depth is compatible exactly on that cubic locus.

The polynomial `Q` is squarefree with nonzero constant. Therefore its three
distinct nonzero `Y` roots each have two square roots `Phi`, giving exactly
six geometric source points. The terminal coefficient rank is three on
eleven fresh row-ten coordinates, so each has an affine-eight terminal fibre.

## 6. The new hierarchy consumer is distinct but redundant here

THM-4364 proves that

```text
L_(10,11,3)(A)
 =35a_(6,1)-20a_(7,3)+10a_(8,5)-4a_(9,7)+a_(10,9)       (14)
```

annihilates `pi_10(P_2)`. After selecting the `P_3` terminal solution, direct
evaluation gives

```text
L_(10,11,3)(A)
 =-2Q(Phi^2)/(23072272859239231201171875Phi^5).          (15)
```

In the opposite order, select `P_2` first. The known functional

```text
L_(10,10,3)(C)
 =56c_(5,0)-35c_(6,2)+20c_(7,4)-10c_(8,6)
  +4c_(9,8)-c_(10,10)                                   (16)
```

annihilates `pi_10(P_3)` and reads

```text
L_(10,10,3)(C)
 =Q(Phi^2)/(15381515239492820800781250Phi^5).            (17)
```

Consequently `(15)=(-4/3)(17)`. The two functionals live on distinct modules
and are opposite-order certificates, but their source ideals on `(10)` are
both `(Q)`. The new functional is not an additional cut. A finite exact scan
of every THM-4364-admissible row-ten hierarchy functional finds these as the
only nonzero opposite-order evaluations.

The rational point `Phi=1` on `(10)` is a hostile to confusing either
standalone system with the joint one: it satisfies the bracket and each depth
block separately, with all strict gates, but `Q(1)!=0`.

## 7. Strict gates and the sign sheet

Put

```text
W_7=1-2Y,
D_7=-Y^4+2Y^3-Y+1,
Q_7=-Y^3+3Y^2+3Y+2.                                    (18)
```

Here `W_7` is the reduction of the cleared `W` factor, and `D_7` is the
reduction of the primitive cleared numerator of
`alpha^2-4W upsilon_5`. Exact identities in `F_7[Y]` are

```text
(Y-1)Q_7+(2Y^2+3Y+1)Q_7'=1,
3Q_7+(2Y^2+2Y+2)W_7=1,
(-Y^3-Y^2+Y-3)Q_7+Y^2D_7=1.                            (19)
```

The reductions retain their displayed degrees. Hence `Q` is squarefree and
coprime to both strict-gate factors in characteristic zero. Together with
`Q(0)!=0`, the fixed nonzero `Z` in `(10)`, and nonzero `upsilon_5`, this
proves every one of the six points satisfies all declared strict gates.

The involution

```text
Phi -> -Phi,          (eta,alpha,beta,zeta_3) ->
                      -(eta,alpha,beta,zeta_3)           (20)
```

fixes `xi_10,W,Z` and exchanges the two sign sheets over each `Y` root.

## 8. Row-eleven bracket extinction

Choose the row-ten `P_2` terminal solution; modulo `Q(Phi^2)` it also
satisfies `P_3`. The row-eleven bracket selector has rank eight on the
surviving affine-eight row-ten terminal fibre. Solving those eight directions
leaves

```text
r_9=8Q(Phi^2)/(84598333817210514404296875Phi^5),

r_8=-R(Phi^2)/(71525718603423911567894897460937500Phi^6), (21)
```

where

```text
R(Y)=6846329377771290182382546697998046875Y^4
 -713835723041306505264998768716800000000000Y^3
 -2754991513504883058403611855707575418880000000Y^2
 -31916203206707002973657986739896008646412206080000Y
 -156854967149983010817497418504735580308619473018945536. (22)
```

Modulo seven,

```text
R_7=-2Y^4+3Y^3+3Y^2+3Y-2,

(3Y^3-2Y^2-3Y)Q_7+(2Y^2-2Y+3)R_7=1.                  (23)
```

Thus `gcd(Q,R)=1`. On the six-point row-ten carrier, `r_9=0` and `r_8!=0`.
No survivor passes the row-eleven bracket, and no row-eleven terminal fibre
exists. Row-eleven projected-depth calculation is unnecessary.

## 9. Audit and scope

The 117-check primary reconstructs both source branches, row-nine equations,
all bracket selectors, four depth matrices, both hierarchy evaluations, the
finite hierarchy scan, strict-gate certificates, sign involution, and
row-eleven gcd obstruction. The 461-assertion import-free referee rebuilds
these objects without importing the primary. Normal, optimized, isolated,
hash-seeded, and frozen LF streams agree.

The theorem is confined to exact finite Hasse rows in THM-4308's declared
source-normal residual-weight universe. It does not prove full `P_2/P_3`
membership, extension beyond finite rows, polynomial termination, source-
chart entry for an arbitrary hypothetical counterexample, seam entry, a
Keller pair, `JC(2)`, or `DC(2)`.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_thm4366.py
python3 -B -O 04-computation/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_thm4366.py
python3 -B 04-computation/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_independent_referee_thm4366.py
python3 -B -O 04-computation/jc2_source_normal_u_zero_row11_hierarchy_selected_extinction_independent_referee_thm4366.py
```
