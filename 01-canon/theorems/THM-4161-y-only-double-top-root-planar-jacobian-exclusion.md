---
id: THM-4161
title: "Complete Y-only nontriple double-top-root planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
  + VERIFIED-EXACT + INDEPENDENT SOURCE-CHART AUDIT. On the Y-only
  eta=Delta=0 exact-weight-nine wall, every top cubic with one double and one
  simple root, zeta!=0, and I_C!=0 is excluded. The J_top!=0 stratum has
  (g,L,packet)=(9,21,(8,5,3,2,2,2,1)); the J_top=0 ordinary-node stratum has
  (8,20,(8,3,2,2,2,2,2,1)); its unique cusp has
  (8,18,(8,3,3,2,2,2,1)). THM-4164 subsequently closes the triple-root
  locus off `I_C=0`; the common `I_C=0,Disc(C)=0` intersections, other cells,
  entry, M>=10, JC(2), and DC(2) remain OPEN.
source: codex-lrc14-jc-sharp-fronts-20260825b
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
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
related:
  - THM-4157-repeated-top-edge-wall-planar-jacobian-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
script: 04-computation/jc23_y_only_double_top_root_exclusion_thm4161.py
output: 05-knowledge/results/jc23_y_only_double_top_root_exclusion_thm4161.out
secondary_wall_script: 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161.py
secondary_wall_output: 05-knowledge/results/jc23_y_only_double_top_root_secondary_wall_thm4161.out
secondary_wall_independent_audit_script: 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161_independent_audit.py
secondary_wall_independent_audit_output: 05-knowledge/results/jc23_y_only_double_top_root_secondary_wall_thm4161_independent_audit.out
independent_audit_script: 04-computation/jc23_y_only_double_top_root_exclusion_thm4161_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_y_only_double_top_root_exclusion_thm4161_independent_audit.out
supplementary_direct_root_script: 04-computation/jc23_y_only_double_top_root_wall_thm4161.py
supplementary_direct_root_output: 05-knowledge/results/jc23_y_only_double_top_root_wall_thm4161.out
supplementary_direct_root_independent_audit_script: 04-computation/jc23_y_only_double_top_root_wall_thm4161_independent_audit.py
supplementary_direct_root_independent_audit_output: 05-knowledge/results/jc23_y_only_double_top_root_wall_thm4161_independent_audit.out
script_sha256: a56f09c2caa143023158920149c934720f19d0cb63b7d15f34473e0c785f9034
output_sha256: 17536d51e1cb6028a891a858786402971c2d8936f86063136d7a5b55395c6206
secondary_wall_script_sha256: 4ee43862a1dcab4f3ab78f82575c09d4fb180f60aab38efd288844b89d2fae06
secondary_wall_output_sha256: 2076d218fc5675643203d921fe45ff236391a347e4792ff01056b878ace55de2
secondary_wall_independent_audit_script_sha256: 41ae504d2f3b876ec1886373d3288b6210eca77de55b4f5b6fcab9f17319513b
secondary_wall_independent_audit_output_sha256: 7812dde796acae90ed3cc4736597c7cfdf523ee0e011983447a35ff6e63c53ba
independent_audit_script_sha256: 20b1848044b6cdec82c58aa8281fa14086476f8fe686e693ea91537202ce86c7
independent_audit_output_sha256: e39b5a0bd1924cb60a6d901c505482780cdc2da96e9e8fea77fc8f676fc9e801
supplementary_direct_root_script_sha256: 6f2181240b8fd5897f922f3c26067f4c91031f869b77e8bb14dfc89df438cbfd
supplementary_direct_root_output_sha256: 6ace665c7d3d9449fd04b0cad9913ef3c1bce955b6ae04513a5fb9102195269b
supplementary_direct_root_independent_audit_script_sha256: 9846b3d1e659af51ffe5f2f3bc5dd841c8184c0274e01f617c334413bfc289c3
supplementary_direct_root_independent_audit_output_sha256: 8ecdced043a929d942cc952496abac6b2d9a902d8fa9c3a3df5948462b767bc9
semantic_sha256: 69d972a9d14186b328695dc08804eebffb089449648ff3cd20ddec53c1fd0b25
secondary_wall_semantic_sha256: 8943b86c9bc0ea5ec9643939d57e5ae6f714482b55591a550b9abb3c29bb29ef
secondary_wall_independent_semantic_sha256: 1e9ac64e08b9d18ec7560f8e9704c57b8e4d25e10968050ffd13ebe10c75b079
independent_semantic_sha256: 632743f7b69862cb9cc39c0abc80060087ed08c5c45e23141c773e83ec0169be
supplementary_direct_root_semantic_sha256: 9e00cdbf68489a26e824458a177af856cda10244c613f062662896b96b4556ff
supplementary_direct_root_independent_semantic_sha256: ca74eeb5e783553cd866a708b8f08d119c7fd18ea561ac10762940e5cb38755d
hash_basis: raw LF bytes
primary_audit: >
  PASS. The generic-wall certificate reconstructs both critical projections,
  resolves the smooth quadratic tangency, and verifies both carrier
  responses. Normal, optimized, and hash-seeded outputs byte-match.
secondary_wall_audit: >
  PASS. Exact quadratic-field reduction exhausts J_top=0, proves the unique
  node-to-cusp transition, reconstructs two independent critical projections
  on both strata, and closes the cusp equality with the carrier-orbit lemma.
  Normal, optimized, and hash-seeded outputs byte-match.
secondary_wall_independent_audit: >
  ACCEPT. A separate source-chart referee uses only the alternative
  critical pair (A,C_0), a disjoint ordinary-node control, and independent
  node/cusp, packet, universal-point, and carrier-orbit reconstructions.
independent_audit: >
  ACCEPT. A clean-room source-chart referee recovers the double root from its
  linear subresultant, uses the alternative critical pair (A,C_0), a disjoint
  control, and independently reconstructs the generic packet and responses.
supplementary_direct_root_audit: >
  ACCEPT. Two further exact certificates use the reciprocal-free direct-root
  chart C=zeta(W-r)^2(W-u). They independently reconstruct I_C, both degree-17
  endpoint resultants, the quadratic top contact, ramification index five,
  packet, genus, and response contradictions on the J_C*I_C!=0 subchart.
---

# THM-4161 -- complete Y-only nontriple double-top-root exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
+ VERIFIED-EXACT + INDEPENDENT SOURCE-CHART AUDIT; JC(2) REMAINS OPEN.**

## 1. Statement and intrinsic chart

Put

```text
kappa=1376/135,                 K0=2848/45,
C(W)=zeta W^3+Theta W^2+Phi W-kappa,
I_C=4 Theta K0^2-27 zeta^2.
```

Define the linear first nonzero subresultant row

```text
A_C=6 Phi zeta-2 Theta^2,
B_C=-9 kappa zeta-Phi Theta,
S_1(W)=zeta(A_C W+B_C),
J_C=15 A_C^2+356 B_C^2.
```

> **Theorem.** The coefficient locus
>
> ```text
> eta=Delta=0, zeta!=0, Disc(C)=0,
> C has no triple root, I_C!=0
> ```
>
> contains no nonautomorphic planar Keller pair.

The inheritance pass is exact.  The closest proved mechanism is THM-4155's
critical-length and labelled-carrier obstruction.  Its canonical hostile is
the top-face collision `D_C=0`, where the degree-seventeen endpoint used there
vanishes.  The corrected near miss is THM-4157: its repeated slanted edge has
the same numerical length sidecar on one row, but a different polygon,
normalization loss, and packet, so it cannot be transported to this wall.
The least-used relevant sidecars are the strict transform at the collided
horizontal root and the orbit count of the handle subgroup; together they
separate the node from its unique cusp and repair the finite cusp equality.

The chart covers the stated locus. A cubic with nonzero constant and leading
coefficient, zero discriminant, and no triple root has a unique double root
`r` and a unique remaining root `s`; both are nonzero. Put `a=r^-1` and
`b=s^-1`. Factoring the constant term gives uniquely

```text
C(W)=kappa(aW-1)^2(bW-1),
zeta=kappa a^2 b,
Theta=-kappa(a^2+2ab),
Phi=kappa(2a+b),
a*b*(a-b)!=0.
```

Conversely these formulas give exactly such a cubic when `a!=b`. On the
discriminant wall `S_1` is proportional to the gcd of `C,C'`, and
`r=-B_C/A_C=1/a`. Direct substitution gives

```text
J_C=(14339490709504/332150625) a^2(a-b)^4(15a^2+356).
```

Thus the intrinsic gate `J_C!=0` is exactly

```text
J=15a^2+356 !=0
```

away from the already excluded triple-root row. Sections 2--4 first treat
`J_C!=0`; Section 5 closes `J_C=0` completely under the same nontriple and
`I_C!=0` hypotheses. The displayed `I_C` is exactly the inner-endpoint factor
of THM-4155, not a new symbol.

## 2. The `J_C!=0` source and exact critical length

Use

```text
s=XT, p=T+s^2, t=p-s^2, P=T+X^2T^2, Y=XTP,

H=-3p+(8/3)p^2-(1376/135)p^3+K0 s^2p^2
  +Phi sp^3+Theta s^2p^3+zeta s^3p^3,
G=-s^2/(2t)+H.
```

For the THM-4155 source-critical pair

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2,
```

exact elimination after the `(a,b)` substitution gives

```text
Res_s(A,B)=p^6 R_17(p),
R_17(0)=-46656 zeta I_C,

[p^17]R_17=
 -(3289935900927224469054816256/252226880859375)
  a^11 b^5(a-b)^3(15a^2+356).
```

The leading `s`-rows remain `3zeta p^2` and `9zeta p^2`, so no critical
intersection is lost at `s=infinity` for `p!=0`.

Independently reconstruct the normalized polynomial and put

```text
f=G_X/T, h=G_T.
```

The second exact elimination is

```text
Res_X(f,h)=T^56(6T+1)^2 Q_17(T),
Q_17(0)=-(3^15/2^7)zeta^7,

[T^17]Q_17=
 -(210555897659342366019508240384/9308590679915771484375)
  a^9 b^3(a-b)^3(15a^2+356)
  (17415a^3b^2+1013888a+2027776b)^2.
```

The last squared factor is a nonzero scalar multiple of `I_C/a`. All
endpoints are therefore live on the stated locus. The `T^56` factor is
a Sylvester artifact. As in THM-4155, the actual critical ideal restores

```text
T=0, X^2=-6, G=0, det Hess(G)=+6,
T=-1/6, X^2=6, G=1/2, det Hess(G)=-6.
```

Consequently a hypothetical Keller realization has reduced affine critical
scheme of exact length

```text
L=17+2+2=21.
```

No discriminant of `R_17` or `Q_17` is required. Projected-coordinate
collisions preserve scheme length, while Keller Hessian congruence makes the
actual critical points Morse. The exact control `(a,b)=(1,2)` has

```text
(Theta,Phi,zeta)=(-1376/27,5504/135,2752/135),
J=371, I_C=-9051394048/10935,
```

and both residual polynomials are squarefree of degree seventeen.

## 3. Horizontal-edge strict transform

Clear the generic-fibre equation as

```text
F=q(s^2-p)-(s^2-p)H-s^2/2.
```

Set `u=1/p`, multiply by `u^4`, and center `v=s-1/a`. The exact boundary
polynomial `K=u^4F(s,u^-1)` begins

```text
K=-(1376/135)a(a-b)v^2+[8(15a^2+356)/(45a^2)]u
  +terms of total (v,u)-order at least two not already displayed.
```

Both displayed coefficients are nonzero on the stated locus. Hence the
closure is smooth at the collided top root and

```text
u=[172a^3(a-b)/(3(15a^2+356))]v^2+O(v^3).
```

On the fibre, differentiation of `F=u^-4K` at fixed `s` gives

```text
F_p=-u^-2 K_u,                 omega=ds/F_p=-u^2 ds/K_u.
```

Thus `ord_v(omega)=4`, so the two former rational index-three punctures merge
into one rational index-five puncture. The remaining simple root gives one
rational index-three puncture.

The Newton polygon itself stays

```text
(0,1),(2,0),(5,3),(3,4),(0,4),
(2Area,B,I)=(27,11,9).
```

Because the collision point is smooth, there is no delta loss: the
normalization genus remains nine. All other faces are unchanged from
THM-4155. The complete packet is therefore

```text
(8,5,3,2,2,2,1),
rational (8,5,3,1) + one cubic orbit (2,2,2).
```

Its defect is

```text
7+4+2+3=16=2*9-2,
```

so Riemann--Hurwitz also proves completeness and rules out hidden affine
ramification.

The carrier face remains

```text
q-1/2=K0 W^2+zeta W^3.
```

Since `zeta!=0`, the rational function on the right has degree three.
Therefore `C(W)/C(q)` has prime separable degree three. The finite-separable
carrier lemma from THM-4147/4155 applies without alteration.

## 4. The two responses and strict contradictions

Every rational puncture maps to the target origin by THM-4120. Prime-degree
carrier transport leaves exactly

```text
full: n=23;
finite: n=17, beta=3.
```

For the full response, the two handle supports cover all sheets and

```text
|supp X_0 intersect supp X_1|<=n-L=2.
```

The commutator-index lemma gives index at most four, whereas the origin
meridian has packet defect sixteen. Thus `4<16` is a contradiction.

For the finite response the exact generator-index capacities are

```text
both handles nonidentity: (2n-L-2)+beta=14,
exactly one identity:     (2n-L-1)+beta=15,
both identities:                         beta=3.
```

All are below the transitive minimum `n-1=16`. This is the second contradiction and excludes the `J_C!=0` stratum.

## 5. Exact closure of `J_C=0`

Since

```text
J_C=(14339490709504/332150625)a^2(a-b)^4(15a^2+356),
```

every nontriple point on `J_C=0` has

```text
a^2=-356/15,                 b=ca.
```

Conversely these formulas parameterize the whole nontriple `J_C=0` wall.
Here `c=0` is exactly `zeta=0`, while `c=1` is exactly the triple-root row.
Exact reduction in `Q(c)[a]/(15a^2+356)` gives

```text
zeta=-489856ac/2025,
Theta=489856(2c+1)/2025,
Phi=1376a(c+2)/135,

I_C=(1986636480512/20503125)(387c^2+80c+40).
```

Thus the theorem-local wall universe is

```text
c*(c-1)*(387c^2+80c+40)!=0.                         (J1)
```

### 5.1 Generic ordinary-node row

Reducing the two ambient resultants coefficientwise modulo `15a^2+356`
gives, away from `43c-33=0`,

```text
Res_s(A,B)=p^6 R_16(p),
[p^16]R_16=
 -(156350393282470363396924806083371672901760279642112
   /1077383180545806884765625)
  a c^5(c-1)^2(43c-33),
R_16(0)=-46656 zeta I_C,

Res_X(f,h)=T^56(6T+1)^2 Q_16(T),
[T^16]Q_16=
 (451470050926386584017169136106653470203416750403666143346688
  /66269166070884740352630615234375)
 a c^3(c-1)^2(43c-33)(387c^2+80c+40)^2,
Q_16(0)=-(3^15/2^7)zeta^7.
```

The leading source-coordinate rows remain nonzero, and the normalized chart
has the same universal factors as Section 2. Hence

```text
L=16+2+2=20.
```

The exact algebraic control

```text
a=2sqrt(-1335)/15,             c=2
```

makes both degree-sixteen residuals squarefree over `Q(sqrt(-1335))`.

On the boundary, with `v=s-1/a`, the vanished `u`-linear row exposes

```text
K=[-489856(c-1)/2025]v^2-(16a/3)uv-3u^2+O_3(v,u),
disc_quad(K)=-45568(43c-33)/675.
```

Therefore `c=33/43` is the only quadratic-form degeneracy. For every other
`c` in `(J1)` the point is an ordinary node. Each normalization branch has
`u` of order one, `K_u` of order one, and
`omega=-u^2ds/K_u` of order one; the two branches have indices `(2,2)`.
The node has delta invariant one, so the genus drops from arithmetic genus
nine to eight. The complete packet and responses are

```text
packet=(8,3,2,2,2,2,2,1),       defect=14=2*8-2,
rational=(8,3,2,2,1) + cubic carrier=(2,2,2),
full n=22,                       finite (n,beta)=(16,3).
```

For the full response, `n-L=2` gives commutator index at most four, below
defect fourteen. For the finite response, the exact capacities are
`(13,14,3)`, all below `n-1=15`.

### 5.2 The unique cusp `c=33/43`

At `c=33/43`, put

```text
y=v-(135a/2848)u.
```

The exact strict transform begins

```text
K=(22784/405)y^2-(q+1161/80)u^3+higher Newton terms.
```

This is one rational `(2,3)` cusp. Indeed, with `t=y/u`, its normalization
has `u=(22784/405)t^2/(q+1161/80)+...`; no quadratic residue extension is
introduced. In the original coordinates `ds` has order one, `K_u` order
three, and `u^2` order four, so `ord(omega)=2` and the cusp branch has index
three. Its delta invariant is one, hence the genus remains eight.

At this point the degree-sixteen rows vanish by two further degrees, not one:

```text
Res_s(A,B)=p^6 R_14(p),
[p^14]R_14=
 (830271724986997884692231133070215530174115479552
  /886735128021240234375)a,
R_14(0)=5247065779999061573632a/18984375,

Res_X(f,h)=T^56(6T+1)^2 Q_14(T),
[T^14]Q_14=
 -(3973312205663515142507987786054543929167803303620264264254619648
   /490882711636183261871337890625)a,
Q_14(0)=
 -171036442361017141608902091573694117707776a/15016937255859375,

Q_14(-1/6)=
 -34951230658071442505568170678029386006915977428361031253491712a
  /2347876792391803819849491119384765625.
```

Both residuals are squarefree over `Q(sqrt(-1335))`; the two embeddings give
the two conjugate coefficient points. Thus

```text
L=14+2+2=18,
packet=(8,3,3,2,2,2,1),        defect=14=2*8-2,
rational=(8,3,3,1) + cubic carrier=(2,2,2),
full n=21,                      finite (n,beta)=(15,3).
```

The full response has overlap at most `n-L=3`, so commutator index at most
six is still below fourteen. The crude finite capacity has a one-identity
equality and must not be cited as strict.

Use the following carrier-orbit refinement. Let `tau_1,...,tau_m` be
transpositions and suppose `n>m+1` and
`<X,Y,tau_1,...,tau_m>` is transitive. If `U=supp X union supp Y`, then
`H=<X,Y>` has at most `m+1` orbits, hence rank at least `n-m-1>0`. Since
`H` fixes the complement of `U` and has at least one orbit in `U`, its rank
is at most `|U|-1`; therefore `|U|>=n-m`. Fixed-sheet transport gives
`|supp X|+|supp Y|<=2n-L`, whence

```text
|supp X intersect supp Y|<=n+m-L.                    (J2)
```

For `(n,m,L)=(15,3,18)`, the hypothesis `15>4` holds and `(J2)` forces zero
handle overlap. The punctured-torus relation and permutation-index
subadditivity then give

```text
ind(mu_O)<=ind([X,Y])+sum ind(tau_i)<=0+3=3.
```

But the finite rational packet `(8,3,3,1)` gives `ind(mu_O)=11`. This final
contradiction closes the cusp and therefore all of `J_C=0` in theorem scope.

## 6. Distinction from THM-4157 and stopping walls

THM-4157 treats `zeta=-eta`, `Delta!=0`, and a repeated slanted length-two
edge in an arithmetic-genus-eleven polygon. Its row `B` loses one
normalization genus and has two index-six branches. Here `eta=Delta=0`; two
roots collide on a horizontal length-three edge, the strict transform is
smooth, genus remains nine, and there is one index-five branch. The two
routes happen to share the numerical sidecars

```text
g=9, defect=16, L=21, full n=23, finite (n,beta)=(17,3),
```

but their packets and collision maps are different and neither theorem
implies the other.

The combined calculation stops at exactly

```text
I_C=0: inner critical endpoint wall;
C triple: a=b, requiring a cubic strict transform.
```

There is no remaining `J_C` exception: its ordinary-node and cusp strata are
both included above.

THM-4159 subsequently closes `I_C=0` away from the top collision, and
THM-4164 closes the triple-root locus away from `I_C=0`. Their common
`I_C=0,Disc(C)=0` intersections remain outside all three theorems.

## 7. Supplementary direct-root chart audit

An independent reciprocal-free chart writes

```text
C(W)=zeta(W-r)^2(W-u),       zeta=kappa/(r^2u),
P_I=2027776r^3u+1013888r^2u^2+17415,
I_C=-44032P_I/(273375r^4u^2).
```

Under

```text
r*u*(r-u)*(356r^2+15)*P_I != 0,
```

this is exactly the `J_C*I_C!=0` part of Sections 1--4 after
`(a,b)=(1/r,1/u)`. Direct elimination, without reciprocal substitution,
independently gives

```text
Res_s(A,B)=p^6R_17,
R_17(0)=3877634048P_I/(50625r^6u^3),

Res_X(f,h)=T^56(6T+1)^2Q_17,
Q_17(0)=-72965752821794209792/(56953125r^14u^7).
```

Both leading coefficients contain `(r-u)^3(356r^2+15)`; the normalized one
also contains `P_I^2`. The direct top chart has unit `L_z` and
`z~(W-r)^2`, hence `ord(omega)=4` and ramification index five. It reproduces

```text
(g,L,packet)=(9,21,(8,5,3,2,2,2,1))
```

and all full and finite response contradictions. A second clean-room script
reconstructs the same endpoint factors, direct-root chart, and the control
`(r,u)=(1,2)`. This is supplementary corroboration of the generic subchart;
Section 5, not this audit, closes `J_C=0`.

## 8. Exact artifacts and replay

```text
python3 -B 04-computation/jc23_y_only_double_top_root_exclusion_thm4161.py
python3 -B -O 04-computation/jc23_y_only_double_top_root_exclusion_thm4161.py
PYTHONHASHSEED=211 python3 -B 04-computation/jc23_y_only_double_top_root_exclusion_thm4161.py

python3 -B 04-computation/jc23_y_only_double_top_root_exclusion_thm4161_independent_audit.py
python3 -B -O 04-computation/jc23_y_only_double_top_root_exclusion_thm4161_independent_audit.py
PYTHONHASHSEED=271 python3 -B 04-computation/jc23_y_only_double_top_root_exclusion_thm4161_independent_audit.py

python3 -B 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161.py
python3 -B -O 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161.py
PYTHONHASHSEED=313 python3 -B 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161.py

python3 -B 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161_independent_audit.py
python3 -B -O 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161_independent_audit.py
PYTHONHASHSEED=347 python3 -B 04-computation/jc23_y_only_double_top_root_secondary_wall_thm4161_independent_audit.py

python3 -B 04-computation/jc23_y_only_double_top_root_wall_thm4161.py
python3 -B -O 04-computation/jc23_y_only_double_top_root_wall_thm4161.py
PYTHONHASHSEED=127 python3 -B 04-computation/jc23_y_only_double_top_root_wall_thm4161.py

python3 -B 04-computation/jc23_y_only_double_top_root_wall_thm4161_independent_audit.py
python3 -B -O 04-computation/jc23_y_only_double_top_root_wall_thm4161_independent_audit.py
PYTHONHASHSEED=127 python3 -B 04-computation/jc23_y_only_double_top_root_wall_thm4161_independent_audit.py
```

This theorem does not prove `JC(2)` or `DC(2)`. **QED.**
