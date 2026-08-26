---
id: THM-4183
title: "P-only Delta-zero planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4176/4180
  + VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED. On the live exact-weight-
  nine b=d=0 reduced (2,3) seam, the whole P-only Delta-zero wall zeta=0,
  Delta=0, eta!=0 contains no nonautomorphic planar Keller pair. The two
  exhaustive rows Theta!=0 and Theta=0 have critical lengths 24 and 23 and
  complete packets (8,5,4,3,2,2,1) and (8,7,4,2,2,1). At Theta=0 one
  critical root escapes in each projection and the index-three/index-five
  faces blow down to one primitive index-seven face. Exact full and prime-
  quadratic-carrier response bounds exclude both rows without discriminant,
  Q_d(-1/6), or Phi hypotheses. Separately, exact-M=9 source completeness
  reclassifies eta=zeta=0 as a filtration exit, so THM-4176/4180 already
  cover the complete repeated-top exact-M=9 locus. Entry, other P/B walls,
  other cells, M>=10, JC(2), and DC(2) remain OPEN.
source: jc-eta-zero-normal-20260826
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
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
  - THM-4180-repeated-top-delta-zero-planar-jacobian-exclusion
related:
  - THM-4045-live-two-three-max-seven-hidden-elliptic-tail-no-go
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
script: 04-computation/jc23_p_only_delta_zero_complete_exclusion_thm4183.py
output: 05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_thm4183.out
independent_audit_script: 04-computation/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.py
independent_audit_output: 05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.out
script_sha256: da9405e43824f91611768e6d9ac6d554fb8a3a45f7f5190ee7183e78260aadd6
output_sha256: b8fb656c76e108b5ec205f9fdcd5b7a24535760f3bede17d77e9e21a58ac4a34
independent_audit_script_sha256: a0a214a023bd270214de60a4d8818431474116b45501de73357184814390e186
independent_audit_output_sha256: 41b3458b1cabd5d3b43d92b5ff1a2f6631a8d772f51cf7d78eb9113f91668955
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone normalized (X,T) certificate reconstructs the complete
  source, proves the Morse-resultant identity, computes the symbolic
  three-parameter resultant and its direct Theta-zero specialization,
  verifies both reciprocal normal coefficients, rebuilds every Newton face,
  and checks both packet and response ledgers. Normal, optimized, and fixed-
  hash executions byte-match the frozen output.
independent_audit: >
  ACCEPT. A separate rational (s,p) implementation imports no primary code.
  It derives a different polynomial critical pair, checks its Hessian bridge
  by ideal reduction, recomputes p-resultants and the linear reciprocal
  escape gauge, reconstructs the polygons by supporting inequalities, and
  obtains the same lengths, packets, and strict response bounds. Normal,
  optimized, and fixed-hash executions byte-match the frozen output.
---

# THM-4183 -- P-only Delta-zero planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4176/4180
+ VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED; JC(2) AND DC(2) REMAIN
OPEN.**

## 1. Theorem and inheritance

Work over `C` in the live `b=d=0` reduced `(2,3)` seam at exact residual
weight nine. Put

```text
P=T+X^2T^2,                  Y=XTP,
K_0=2848/45,                 epsilon=-1376/135.          (1)
```

Use the complete normalized source

```text
G=-X^2T/2-3P+(8/3)P^2+epsilon P^3+K_0Y^2+Phi P^2Y
  +Theta P Y^2+eta P^3Y.                                (2)
```

> **Theorem.** For every
>
> ```text
> (Phi,Theta,eta) in C^3,                    eta!=0,     (3)
> ```
>
> polynomial `(2)` is not the normalized exact-weight-nine source of a
> nonautomorphic planar Keller pair in the inherited reduced seam.
> Equivalently, the complete coefficient wall
>
> ```text
> zeta=0,                    Delta=0,          eta!=0    (4)
> ```
>
> is empty of nonautomorphic planar Keller realizations.

There are exactly two exhaustive, disjoint rows:

```text
P_A: Theta!=0;                         P_0: Theta=0.    (5)
```

The coefficient `Phi` is unrestricted in both rows. No residual
discriminant, selected projection-fibre separation, or value condition at
`T=-1/6` is imposed.

The closest proved mechanism is THM-4147's P-only quadratic-carrier response.
The hostile is `Theta=0`: its nominal degree-twenty critical endpoint dies
and two Newton edges disappear. The corrected near miss is MISTAKE-509: the
two geometric carrier branches are one quadratic closed point and must
respond together. The least-used sidecar is the first reciprocal normal
coefficient, computed independently in `u=1/T` and `v=1/p` below.

## 2. Source completeness and the Keller contract

In normalized rational coordinates

```text
s=XT,                 p=T+s^2,              y=sp,
t=p-s^2=T,            K=2848/45-(7/6)Delta,             (6)
```

THM-3992/3997 and THM-4007 give

```text
G=-s^2/(2t)+H(p,y),
H=-3p+(8/3)p^2+epsilon p^3+K y^2+Phi p^2y
  +Delta p^4+Theta p y^2+eta p^3y+zeta y^3.            (7)
```

For `wt(p)=2,wt(y)=3`, the only weight-nine monomials are `p^3y` and
`y^3`. Thus exact weight nine is exactly

```text
(eta,zeta)!=(0,0).                                      (8)
```

On `(4)`, equation `(7)` specializes to `(2)`, `K=K_0`, and `eta p^3y`
keeps exact weight nine. There is no omitted term and no extra structural
row.

The target/source contract is the one proved and retained by THM-4147:

```text
source: normalization of the generic fibre G=q, together with all affine
        inverse points of the two nodal target values;
target: the inherited elliptic generic fibre E=q over C(q), with precisely
        its two Morse nodes and E_q(C(q))={O};
map:    the finite morphism induced by a hypothetical Keller pair after
        compactification, resolution, relative normalization, and shrinking;
kept:   complete infinity packet, residue-field degree, affine fixed-sheet
        supply L, finite-carrier index beta, and origin defect;
lost:   coordinates/order of conjugate carrier points and the value of a
        chosen source projection coordinate.                            (9)
```

At an affine critical point of a hypothetical realization `G=E(A,C)`, the
unit Jacobian gives

```text
Hess(G)=D(A,C)^t Hess(E)D(A,C),
det Hess(G)=det Hess(E)!=0.                              (10)
```

Every affine critical point is therefore Morse and lies above one of the
two target nodes. This is the mechanism that replaces every residual
squarefreeness assumption.

## 3. Complete critical divisors and both strict-transform gauges

Set

```text
f=G_X/T,                         h=G_T.                 (11)
```

The exact polynomial identity

```text
T det D(f,h)=det Hess(G)+f G_XT                         (12)
```

implies that every common zero on `T!=0` is reduced under `(10)`. In both
rows `(5)`,

```text
deg_X(f,h)=(8,9),
[X^8]f=[X^9]h=9eta T^8.                                (13)
```

Hence no finite `T!=0` critical point is lost at `X=infinity`.

Exact symbolic elimination gives

```text
Res_X(f,h)=T^56(6T+1)^2 Q_20(T),

Q_20(0)=-(3^15/2^7)eta^7,
[T^20]Q_20=72900 K_0^4 eta^5 Theta^4.                  (14)
```

Thus `P_A` has actual residual degree twenty. Direct specialization and a
second resultant computation agree on `P_0`:

```text
Res_X(f|_(Theta=0),h|_(Theta=0))
  =T^56(6T+1)^2Q_19(T),

Q_19(0)=-(3^15/2^7)eta^7,
[T^19]Q_19=944784 K_0^6 eta^7.                         (15)
```

All four endpoints in `(14)--(15)` are units in their stated rows.

The degree drop is one genuine strict-transform escape, not an unexamined
wall. Put `u=1/T` and

```text
Qhat(u,Theta)=u^20Q_20(1/u).
```

At `(u,Theta)=(0,0)`, its first two live rows are

```text
Qhat=72900K_0^4eta^5Theta^4
     +944784K_0^6eta^7u+terms in (Theta*u,u^2).          (16)
```

Since the `u` coefficient is a unit, exactly one projected root escapes and

```text
u=1/T= -25Theta^4/(324K_0^2eta^2)+higher terms.         (17)
```

The independent source chart gives a deliberately different gauge. Define

```text
A=(-sp+t^2H_s)/p,
C_0=s^2+2t^2H_p,
B=(C_0+sA)/t^2.                                        (18)
```

These are polynomials and

```text
t^2G_s=pA,                    2t^2G_p=t^2B-sA.          (19)
```

Their `s`-degrees are `(5,2)`, with leading rows

```text
LC_s(A)=2p(K_0+Theta p),
LC_s(B)=8p(Theta p+3K_0/4).                            (20)
```

The two rows in `(20)` cannot vanish together for `p!=0`, so this projection
also loses no finite open-chart point at `s=infinity`. Differentiating `(19)`
at an `(A,B)` zero gives

```text
p det D(A,B)=2t^2 det Hess_(s,p)(G).                   (21)
```

The independent certificate verifies the exact ideal reduction behind
`(21)`. Its resultants are

```text
Res_s(A,B)=p^2R_20(p),
[p^20]R_20=-65610eta^6Theta,       R_20(0)=-31104K_0^2;

Theta=0:
Res_s(A,B)=p^2R_19(p),
[p^19]R_19=-78732K_0eta^6,         R_19(0)=-31104K_0^2. (22)
```

For `v=1/p`, the reciprocal initial form is

```text
v^20R_20(1/v)=-65610eta^6Theta-78732K_0eta^6v
              +terms in (Theta*v,v^2),

v= -5Theta/(6K_0)+higher terms.                        (23)
```

Equations `(17)` and `(23)` are projection-specific gauges; no rootwise
identification between them is claimed. They independently certify the same
one-root loss.

The factor `T^56` in `(14)--(15)` and `p^2` in `(22)` are coordinate
degree-drop artifacts. The normalized chart restores exactly

```text
T=0,       X^2=-6,       G=0,       det Hess(G)=+6;
T=-1/6,    X^2= 6,       G=1/2,     det Hess(G)=-6.     (24)
```

The second pair accounts for `(6T+1)^2`; the first pair is absent from the
`(f,h)` divisor. By `(12)`, additional roots of `Q_d` at `T=-1/6` are
additional distinct Morse points rather than hidden multiplicity. Therefore

```text
P_A: L=20+2+2=24;                 P_0: L=19+2+2=23.    (25)
```

No discriminant or `Q_d(-1/6)` condition occurs.

## 4. Delta-zero polygons and the Theta blowdown

In `(s,p)` coordinates, equation `(2)` is `G=-s^2/(2t)+H`, where

```text
H=-3p+(8/3)p^2+epsilon p^3+K_0s^2p^2+Phi sp^3
  +Theta s^2p^3+eta sp^4.                              (26)
```

For `Q=q^-1`, put

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                              (27)
```

Exact support collection gives

| row | Newton polygon, counterclockwise | `(2Area,B,g_Pick)` |
|---|---|---:|
| `P_A` | `(0,1),(2,0),(4,2),(4,3),(3,4),(1,5),(0,4)` | `(28,10,10)` |
| `P_0` | `(0,1),(2,0),(4,2),(3,4),(1,5),(0,4)` | `(27,9,10)` |

The seven exact faces in `P_A`, in order, are

```text
s^2(1-Q/2)-p,
s^2((1-Q/2)-K_0Q(sp)^2),
-Qp^2s^4(K_0+Theta p),
-Qp^3s^3(eta p+Theta s),
Qeta p^4s(p-s^2),
Qp^4(epsilon+eta sp),
p(-1+Q(-3p+(8/3)p^2+epsilon p^3)).                    (28)
```

The last face has `s=0,p!=0`, hence `X=0,T=p`: it consists of affine source
points and is not part of the infinity packet. In `P_0`, the third and
fourth faces of `(28)` disappear and are replaced by the primitive face

```text
-Qp^2s^3(K_0s+eta p^2).                               (29)
```

Its edge vector is `(-1,2)`, and its torus equation is linear in the
primitive edge coordinate. Thus it is one smooth branch. This is the exact
geometric mechanism: the index-three and index-five faces blow down to one
index-seven face. The boundary loses one degree but no genus or defect.

The edge ledgers are

| row | primitive inward data `(u,v;c)` | indices |
|---|---|---|
| `P_A` | `(1,2;2),(-1,1;-2),(-1,0;-4),(-1,-1;-7),(-1,-2;-11),(1,-1;-4)` | `1,(2,2),3,5,8,4` |
| `P_0` | `(1,2;2),(-1,1;-2),(-2,-1;-10),(-1,-2;-11),(1,-1;-4)` | `1,(2,2),7,8,4` |

The length-two face is one prime quadratic closed point:

```text
q-1/2=K_0W^2,                       W=sp.              (30)
```

It is separable and irreducible over `C(q)` because `q-1/2` has odd
valuation at `q=1/2`. All other displayed infinity faces have a single
rational simple torus root. In particular, the diagonal replacement has

```text
W=-epsilon/eta=1376/(135eta)!=0.                       (31)
```

The coefficient `Phi` occurs on no boundary face. After deleting the finite
set of vertical critical values in the pencil base, all face schemes are
smooth. THM-4103's residue identity therefore gives the candidate packets

```text
P_A: (8,5,4,3,2,2,1),       defect=18,   sum=25;
P_0: (8,7,4,2,2,1),         defect=18,   sum=24.       (32)
```

Equations `(14)--(15)` make the affine critical scheme finite. THM-3827's
closed-polynomial factor theorem then makes the geometric generic source
connected: a nontrivial composition would put a complete curve in that
finite critical scheme. Pick gives normalization genus at most ten, while
Riemann--Hurwitz over the elliptic target and defect eighteen give genus at
least ten. Hence both rows have genus ten and `(32)` is complete. There is
no hidden affine ramification or further boundary branch.

## 5. Quadratic-carrier responses

THM-4120 gives `E_q(C(q))={O}`. Every rational place in `(32)` therefore
maps to the target origin. Let `L_2/C(q)` be the quadratic residue field of
`(30)`. If that point has a finite horizontal target image with residue
field `M`, then

```text
C(q) subset M subset L_2.                              (33)
```

Prime degree and the absence of a finite `C(q)`-point force `M=L_2`.
Consequently the two geometric branches respond together as two conjugate
simple carrier points, with total transposition index `beta=2`. The only two
responses are

| row | finite response `(n,beta)` | full response `n` |
|---|---:|---:|
| `P_A` | `(21,2)` | `25` |
| `P_0` | `(20,2)` | `24` |

The proper-smooth/finite-etale construction and parallel Milnor-core
transport are exactly THM-4147's finite-separable-carrier lemma. Its
hypotheses have now been verified in both rows: connected source, complete
packet, one separable prime carrier, and finite reduced affine critical
scheme. If `X,Y` are the two handle permutations and `r_0+r_1=L`, then

```text
#Fix(X)>=r_0,              #Fix(Y)>=r_1,
|supp X|+|supp Y|<=2n-L.                               (34)
```

## 6. Uniform monodromy contradictions

For a full response, `X,Y` generate transitively. Hence their supports cover
all `n` sheets, their overlap is at most `n-L`, and THM-4147's commutator
bound gives

```text
ind([X,Y])<=2(n-L).                                    (35)
```

The origin meridian has index equal to the complete packet defect eighteen.
But in both rows

```text
P_A: 18<=ind([X,Y])<=2(25-24)=2;
P_0: 18<=ind([X,Y])<=2(24-23)=2,                       (36)
```

a contradiction.

For a finite response, the handles and the two carrier transpositions must
generate transitively. If at least one handle is nonidentity, their total
generator index is at most

```text
2n-L-1+beta.                                           (37)
```

The exact bounds are

```text
P_A: 2*21-24-1+2=19<20=n-1;
P_0: 2*20-23-1+2=18<19=n-1.                           (38)
```

If both handles are identities, total index is only `beta=2<n-1`. No finite
response is transitive. Equations `(36),(38)` exclude both exhaustive rows,
proving the theorem. **QED.**

## 7. Repeated-top exact-weight-nine scope repair

Equation `(8)` repairs the apparent `eta=0` residual in the repeated-top
route. On

```text
zeta=-eta,                                             (39)
```

exact weight nine is equivalent to `eta!=0`; at `eta=0`, both weight-nine
coefficients vanish. Therefore the exact-M=9 repeated-top locus has the
exhaustive disjoint partition

```text
eta!=0,Delta!=0         or         eta!=0,Delta=0.     (40)
```

THM-4176 excludes the first row and THM-4180 excludes the second. Thus:

> **Scope corollary.** Inside the live exact-weight-nine reduced `(2,3)`
> seam, THM-4176 and THM-4180 already exclude the **complete** repeated-top
> locus `zeta=-eta`. There is no additional `eta=0` exact-M=9 coefficient
> wall.

This is a filtration statement, not a transverse realization theorem. The
formal quotient

```text
(eta(P^3Y-Y^3))/eta=P^3Y-Y^3                         (41)
```

is a normal vector to the lower-weight fibre, not the source polynomial of
an `eta=0` Keller pair. There is no base realization on which to linearize
the Keller equation. Lower-weight coefficient points are governed by the
exact `M<=8` canon and are not promoted here as new consequences.

## 8. Exact scope and firewalls

THM-4183 closes exactly `(4)` in the inherited live seam. It crosses all
`Phi`, critical-discriminant, `Q_d(-1/6)`, and `Theta` specializations on
that wall. It does **not** close:

1. P-only `Delta!=0` coefficient or critical walls not already in THM-4147;
2. the mixed B chamber or its coefficient/critical walls;
3. another reduced cell or entry into this seam;
4. exact residual weight at least ten;
5. a transport setup outside the inherited proper-smooth finite map;
6. an arbitrary subset of the quadratic carrier orbit;
7. `JC(2)` or `DC(2)`.

The reciprocal gauges `(17),(23)` record projection-specific loss only. They
do not identify roots between projections, preserve coordinates of escaping
points, or provide a Keller deformation through the lower-weight fibre.

## 9. Exact artifacts and replay

```text
04-computation/jc23_p_only_delta_zero_complete_exclusion_thm4183.py
sha256 da9405e43824f91611768e6d9ac6d554fb8a3a45f7f5190ee7183e78260aadd6

05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_thm4183.out
sha256 b8fb656c76e108b5ec205f9fdcd5b7a24535760f3bede17d77e9e21a58ac4a34

04-computation/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.py
sha256 a0a214a023bd270214de60a4d8818431474116b45501de73357184814390e186

05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.out
sha256 41b3458b1cabd5d3b43d92b5ff1a2f6631a8d772f51cf7d78eb9113f91668955
```

Replay both certificates in normal, optimized, and fixed-hash modes:

```bash
python3 -B \
  04-computation/jc23_p_only_delta_zero_complete_exclusion_thm4183.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_thm4183.out -
python3 -B -O \
  04-computation/jc23_p_only_delta_zero_complete_exclusion_thm4183.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_thm4183.out -
PYTHONHASHSEED=4183 python3 -B \
  04-computation/jc23_p_only_delta_zero_complete_exclusion_thm4183.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_thm4183.out -

python3 -B \
  04-computation/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.out -
python3 -B -O \
  04-computation/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.out -
PYTHONHASHSEED=4183 python3 -B \
  04-computation/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.py \
  | diff -u \
      05-knowledge/results/jc23_p_only_delta_zero_complete_exclusion_independent_audit_thm4183.out -
```

All three primary streams print `checks=49` and
`verdict=P_ONLY_DELTA_ZERO_CLOSES`. All three independent streams print
`checks=91` and `verdict=INDEPENDENT_P_ONLY_DELTA_ZERO_ACCEPT`.
