---
id: THM-2476
title: "Degree-twenty-two B-E plane Hensel-ramification closure"
status: >
  PROVED + VERIFIED-EXACT. In the open first-flux chart of the genuine
  nonsplit polynomial exact-square-prefix degree-twenty-two branch, the
  complete coefficient plane C=D=W=0 is empty. For B,E nonzero, the
  root-free invariants lambda=E^2/B^5 and r=E/(B^2y) give a 25-term
  plane curve of bidegree (10,5). Its Newton triangle reduces every
  factorization to one linear or quadratic v-factor. The linear ideal
  is unit. For a quadratic factor, the squarefree fixed L_5 section
  and the missing odd Hensel orders force the complete six-variable
  factor ansatz, whose exact coefficient ideal is unit. Thus every
  physical fibre is absolutely irreducible. The v-discriminant is
  constant*lambda^14*W_5^2*K_30, where W_5 is exactly the moving
  first-flux wall. A complete exceptional-parameter factorization,
  a Sylvester-valuation bound, and ramification parity force
  normalization genus at least six in every fibre. No rational
  trajectory survives. Together with the axes, this closes the ninth
  of ten support-two planes. It does not close the B-C plane, higher
  mixed strata, JC(2), or DC(2).
source: codex-2026-07-27-degree-twenty-two-BE-hensel-ramification
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
script: 04-computation/jc2_degree22_be_plane_hensel_ramification_thm2476.py
output: 05-knowledge/results/jc2_degree22_be_plane_hensel_ramification_thm2476.out
script_sha256: 58d28830676f2f97b105c59b53c29cf9c11c9f4f56d27cb196efdf0201c03e59
output_sha256: b2edc8a888074f7fb653d58e76269de7b69f424f6686f3e5f07fe648183bd6ac
hash_basis: working-tree bytes (LF)
---

# THM-2476 -- the degree-twenty-two B-E plane is empty

**PROVED + VERIFIED-EXACT.**

The fixed-section Hensel operation introduced in THM-2475 persists when the
weight-two coefficient also moves the first-flux denominator.  Here the
formal quadratic lift is allowed at order two, but sparsity forces its odd
orders to vanish and the complete order-four ansatz is inconsistent.  The
exact conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
C=D=W=0
    => contradiction.                                           (1)
```

Thus the complete `B,E` support-two coefficient plane is empty.

## 1. Coprime-weight invariant coordinates

Use the target-translated coordinates of THM-2411:

```text
y=11s,                   u=dT,                   Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The cases `B=0` and `E=0` are the closed `E` and `B` axes of THM-2425
and THM-2428.  Suppose

```text
B,E!=0.                                                (3)
```

First take `y!=0` in `C(x)` and define

```text
lambda=E^2/B^5 in C*,             r=E/(B^2y).        (4)
```

The lossless identities are

```text
B/y^2=r^2/lambda,                  E/y^5=r^5/lambda^2. (5)
```

Put

```text
v=u/y^2,                           zeta=Z/y^3.         (6)
```

Multiplying the first two normalized fluxes by `lambda^2` gives

```text
F_1
 =1331lambda[616r^2+lambda(63-1089v)]zeta
   +4[lambda(-745360v+6160)r^2-3748096r^5
      +lambda^2(922383v^2-25410v+63)]
 =0,                                                        (7)

F_2
 =lambda r^2[65591680zeta+1443016960v^2
              -71554560v+98560]
   -239878144r^5
   +lambda^2[15944049zeta^2+(-162339408v+2236080)zeta
              -1190488992v^3+147581280v^2-1219680v+672]
 =0.                                                        (8)
```

The open first-flux chart is now

```text
616r^2+lambda(63-1089v)!=0.                         (9)
```

Equation (7) reconstructs `zeta` uniquely there.

## 2. The 25-term eliminant

Exact elimination gives

```text
Res_zeta(F_1,F_2)=28344976lambda^2 P_lambda(r,v),      (10)
```

where

```text
P_lambda
 =126434012626944r^10+22755800252416r^9

  +lambda(-50286255022080v+2299591827456)r^7

  +lambda(34222590223360v^2+3959638538240v
           -44411530240)r^6

  +lambda^2(44449457564160v^2-5142912445440v
             +35419507200)r^5

  +lambda^2(-149234938081152v^3-9500102156160v^2
             +695599766400v-6017413248)r^4

  +lambda^3(206782580709936v^4+6246495741024v^3
             -1509756494400v^2+34466937120v
             -193496688)r^2

  -567lambda^4 L_5(v),                                  (11)
```

and

```text
L_5(v)
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                         (12)
```

The fixed squarefree section and odd gaps are

```text
P_lambda(0,v)=-567lambda^4L_5(v),

[r]P_lambda=[r^3]P_lambda=0.                          (13)
```

## 3. Uniform absolute irreducibility

For every `lambda!=0`, the polynomial `P_lambda` is absolutely irreducible
in `C[r,v]`.

Its Newton polygon is again

```text
N=conv{(0,0),(10,0),(0,5)}=5 Delta,

Delta=conv{(0,0),(2,0),(0,1)}.                        (14)
```

Every nontrivial lattice Minkowski summand is `k Delta` up to translation.
The nonzero constant and `v^5` coefficients normalize factors to be monic
in `v`, so a factorization has a factor of `v`-degree one or two.

The complete line ansatz

```text
v+ar^2+br+c                                            (15)
```

has exact coefficient ideal

```text
I_1=(1) in Q[a,b,c,lambda].                            (16)
```

For a quadratic factor, let `q_0=v^2+cv+h`.  It divides the squarefree
`L_5`, so it is coprime to its cubic cofactor.  The missing order one in
(13) and Hensel uniqueness force

```text
q_1=0.                                                (17)
```

Newton support permits the general order-two term

```text
q_2=av+f.                                             (18)
```

At order three, both `P_3` and the cross-terms involving `q_1` vanish.
The same coprime argument gives

```text
q_3=0.                                                (19)
```

Finally `2 Delta` permits only a scalar term `q_4=d`.  Hence every possible
quadratic factor has the forced complete form

```text
q=v^2+cv+h+(av+f)r^2+dr^4.                            (20)
```

Exact division of (11) gives

```text
I_2=(1) in Q[a,c,d,f,h,lambda].                       (21)
```

Thus neither factor degree occurs.  Let `C_lambda` be the smooth projective
normalization of (11); projection to the `r`-line has degree five.

## 4. Exact ramification divisor

Solving equality in (9) gives

```text
v=(616r^2+63lambda)/(1089lambda).                     (22)
```

Define

```text
W_5(r,lambda)
 =702768r^5+23716r^4-1694lambda r^2-315lambda^2.      (23)
```

Direct substitution yields

```text
P_lambda(r,(616r^2+63lambda)/(1089lambda))
 =256W_5(r,lambda)^2.                                 (24)
```

Thus the moving square is exactly the excluded first-flux wall.

The quintic discriminant factors as

```text
Disc_v(P_lambda)
 =c lambda^14 W_5(r,lambda)^2 K_30(r,lambda),         (25)
```

for `c in Q*`.  Normalize `K_30` to be primitive with

```text
K_30(0,lambda)=-5532126228517587890625lambda^14.      (26)
```

It has degree thirty in `r`, degree fourteen in `lambda`, and 56 nonzero
terms.  Its leading coefficient is

```text
[r^30]K_30
 =10284112691593526727130742784 L_0(lambda),          (27)

L_0
 =980708488959375lambda^2+32820181258788lambda
   +120472576000.                                     (28)
```

Put

```text
A_3
 =19731972545646643418400000000lambda^3
  +346409095257857532829648104lambda^2
  +1415693788892415142663014lambda
  -981468600381938943169.                             (29)
```

Exact factorization gives

```text
Disc_r(K_30)
 =c' lambda^379 A_3^3 P_12^2 Q_12^3 R_16,           (30)
```

where `P_12,Q_12,R_16` are the unique primitive positive-leading
irreducible factors of degrees `12,12,16` with exponents `2,3,1`.

The wall collision divisor is

```text
Res_r(K_30,W_5)
 =c'' lambda^56 A_3^3 S_5(lambda),                    (31)
```

where

```text
S_5
 =131955965143315044357267878722928640000000000000000lambda^5
  +3095105736402260153329374221431569404424300000000lambda^4
  +53687721220771954218792719815073181002030599351104lambda^3
  -6341516493144875869324270617788470732795917450816lambda^2
  +235228509291405819426976849459832940179373568841lambda
  +178772563473178835054245823352852460725645500.    (32)
```

The six polynomials

```text
L_0,A_3,P_12,Q_12,R_16,S_5                           (33)
```

are squarefree and pairwise coprime.  Root and wall collision overlap only
through `A_3`.

At each of the two roots of `L_0`, the `r^29` coefficient, the specialized
degree-30 discriminant, and the wall resultant are nonzero.  The
leading-degree specialization identity therefore makes both degree-29
fibres squarefree and wall-disjoint.

## 5. Uniform ramification floor

The Sylvester valuation lemma used in THM-2469/2470/2472/2475 says that a
discriminant valuation `e` bounds

```text
deg gcd(k_0,partial_rk_0)<=e.                         (34)
```

Every nonzero exceptional factor in (30) is squarefree and has exponent at
most three.  Coprimality with `L_0` makes the leading coefficient a unit.
Thus a degree-30 exceptional fibre has gcd degree at most three, and its
nonsimple roots consume total multiplicity at most six.  At least 24 roots
of `K_30` remain simple.

At an `A_3` root, discarding all five possible wall roots leaves at least
nineteen simple roots off the wall.  The `P_12,Q_12,R_16` roots have no wall
collision.  At an `S_5` root, `K_30` is squarefree and at most five roots are
lost.  The two degree-drop fibres and the generic fibre are squarefree and
wall-disjoint.  Uniformly,

```text
at least 19 simple finite roots of K_30 lie off W_5.   (35)
```

Each contributes one unit of tame ramification: (25) has valuation one, and
normalization changes order discriminant valuation by twice the index.

For the degree-five projection,

```text
2g(C_lambda)-2=-10+R.                                 (36)
```

Total ramification `R` is even.  The nineteen visible units force

```text
R>=20,                    g(C_lambda)>=6.             (37)
```

## 6. Genus and trajectory closure

The rational pair `(r,v)` from (4),(6) would give a nonconstant rational map
from `P^1` to the positive-genus normalization, impossible.  Hence `r,v`
are constant.  Since `r=E/(B^2y)` is nonzero, `y` and then `u` are constant.
Equation (7) reconstructs constant `zeta`, so `Z,T,q` are constants,
contradicting the genuine deck.

It remains to treat `y=0`.  Put

```text
D_0=1331(616B-1089u),              k=14992384E.       (38)
```

The open first-flux chart says `D_0!=0`, and the original first flux is

```text
D_0Z-k=0,                         Z=k/D_0.            (39)
```

The second flux is

```text
15944049Z^2+1443016960Bu^2-1190488992u^3=0.          (40)
```

Substitution and multiplication by `D_0^2` give the constant-field quintic

```text
15944049k^2
 +(1443016960Bu^2-1190488992u^3)D_0^2=0.             (41)
```

Its leading coefficient is nonzero, so `u`, and then `Z,T,q`, are constants,
again contradicting the deck.  Together with the axes, this proves (1).

## 7. Scope and structural lesson

This theorem closes the ninth of ten support-two planes.  Only

```text
(B,C)                                                   (42)
```

remains open.

The moving denominator does not destroy the ramification method: its equality
section can be substituted exactly, and the whole apparent discriminant
square is again identified as excluded wall mass.  On the irreducibility
side, the right object is not an unrestricted factor coefficient ideal but
the unique sparse Hensel lift of a fixed-section divisor.  The first missing
allowed coefficient gives the cheapest obstruction.

Higher mixed strata, branches outside the inherited reduction, split/even
short edges, and integral order raising remain open.  Nothing here proves
`JC(2)` or `DC(2)`.

## 8. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_be_plane_hensel_ramification_thm2476.py
python3 -O 04-computation/jc2_degree22_be_plane_hensel_ramification_thm2476.py
```

The companion reconstructs (7)--(13), the Newton triangle, both factor
obstructions, the moving-wall identity (22)--(24), the complete
discriminant and wall divisors, all coprimality and degree-drop controls, the
branch/genus floors, and the `y=0` quintic.  Normal and optimized transcripts
byte-match the stored output.  All truth-bearing checks use explicit
exceptions and remain active under optimized Python.  **QED.**
