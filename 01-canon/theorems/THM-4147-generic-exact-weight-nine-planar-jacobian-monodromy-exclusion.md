---
id: THM-4147
title: "Generic exact-weight-nine and anti-diagonal planar Jacobian monodromy exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4141/4143
  + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The critical-open
  exact-M=9 chambers P,Y,B and the generic anti-diagonal subchamber A on the
  live b=d=0 reduced (2,3) seam contain no nonautomorphic planar Keller pair.
  Their genera are 10,11,11,10; critical lengths 24,24,25,23; and packets
  (8,5,4,3,2,2,1), (11,8,2,2,2,1), (8,8,4,2,2,2,1), and
  (7,7,4,2,2,2,1). Prime quadratic/cubic carrier transport and sharp merger
  and commutator inequalities exclude every finite and full response.
  Remaining coefficient/critical walls, entry, other cells, M>=10, JC(2),
  and DC(2) remain OPEN.
source: codex-frontier-synthesis-creative-20260825ag
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4141-delta-d-collision-wall-boundary-monodromy-exclusion
  - THM-4143-two-term-collision-wall-critical-boundary-monodromy-exclusion
related:
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
script: 04-computation/jc23_weight9_generic_monodromy_thm4147.py
output: 05-knowledge/results/jc23_weight9_generic_monodromy_thm4147.out
independent_audit_script: 04-computation/jc23_weight9_generic_monodromy_thm4147_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_weight9_generic_monodromy_thm4147_independent_audit.out
script_sha256: e4666241bcee9df55632b7b484e8e06d7a855ed2cabbaedda233a4d6ec056691
output_sha256: 2c9772af2e66d9e6f0b312904ca72c94349822cf12884033a767601c36751096
independent_audit_script_sha256: f0be79df97daa9dc589c96a6c410f2a80dd3d99147e53cc4cbe1ecabc64b8ce2
independent_audit_output_sha256: 32b34a68bf25df98d5de0e41eca857c72ed305cb6335d8702e24187b15077e6a
semantic_sha256: e88dc5ce513e71cf18701434e3c5c14d9e5e22a06655cbc05e89580da4df1a17
extension_script: 04-computation/jc23_exact_weight_nine_generic_antidiagonal_exclusion_thm4147.py
extension_output: 05-knowledge/results/jc23_exact_weight_nine_generic_antidiagonal_exclusion_thm4147.out
extension_script_sha256: 2cbf9a8efb3e048770e8e0dc4815fe0b7c452dbf745c4b4adb56cd7142d8c9f1
extension_output_sha256: 08c53522d7adb4bf0ee44547f1ca517a3a4b906441e202e61a775bbc23f30faf
extension_semantic_sha256: f6dbab29b287dca1829d0a1470030d45b2cd71113b39e3fbf50f3969f8d49d12
extension_independent_audit_script: 04-computation/jc23_weight9_antidiagonal_generic_independent_audit_thm4147.py
extension_independent_audit_output: 05-knowledge/results/jc23_weight9_antidiagonal_generic_independent_audit_thm4147.out
extension_independent_audit_script_sha256: c57d3fb5652ac118ffb923da302f5cc42f8699b4a3175ac2aa7a9f7511fed108
extension_independent_audit_output_sha256: 279bf80b41ef152c717043dff4502878839d19d172a5a7c08fed2325422e5954
extension_independent_semantic_sha256: 6b219119e493f6391e52af1f517d2141e374b4ce7848b19292ee3713342758c6
supplementary_eta_script: 04-computation/jc23_weight9_eta_only_supplementary_crosscheck_thm4147.py
supplementary_eta_output: 05-knowledge/results/jc23_weight9_eta_only_supplementary_crosscheck_thm4147.out
supplementary_eta_script_sha256: 96960d6120176820a11387fac81efe02a5f807b3f5612cd0f1ddb8e25f1e97a6
supplementary_eta_output_sha256: ae712c0cdb751fbbc2b5e285c66bf3d18e08ed96dba8d3d6a3abb71f3e29a0d4
supplementary_eta_independent_script: 04-computation/jc23_weight9_eta_only_supplementary_crosscheck_thm4147_independent_audit.py
supplementary_eta_independent_output: 05-knowledge/results/jc23_weight9_eta_only_supplementary_crosscheck_thm4147_independent_audit.out
supplementary_eta_independent_script_sha256: 57d0b714b240e8c904fe3e2d1a7451033c66fff13f14adfb997c3142329061dc
supplementary_eta_independent_output_sha256: 0af28bd29d2a9921205b8423db2f1f68d35abd42c8602f7ceee9a65057623e8a
supplementary_eta_semantic_sha256: 8d00d0a84d02b99a33b3c71ea39c97fcbf19a14b22c5a6e558c17b5c80aa6cfe
hash_basis: raw LF bytes
primary_audit: >
  PASS. The exact certificate enumerates the complete weight-nine monomial
  universe, reconstructs all three valued lower supports, Newton polygons,
  boundary packets, and two critical projections, restores both collapsed
  coordinate strata, freezes three rational controls, and checks the finite
  and full permutation budgets plus the genuine eta+zeta=0 hostile. Normal,
  optimized, and hash-seeded executions byte-match the frozen output.
independent_audit: >
  ACCEPT after three incorporated scope repairs. A clean-room SymPy referee
  expands the source equation instead of importing the primary support
  routine, uses a disjoint rational control, recovers every polygon, packet,
  critical length, and response budget, and exhaustively verifies the
  commutator-overlap lemma on 533,417 ordered pairs through S_6. It also
  identified the source-completeness dependency, the cubic-carrier Hurwitz
  deletion, and the dominated valued-support point repaired below. Normal,
  optimized, and hash-seeded executions byte-match.
extension_audit: >
  PASS after scope repair. The anti-diagonal certificate binds Delta!=0,
  residual discriminant and Q_19(-1/6), derives the ordinary boundary node,
  the two index-seven branches, critical length 23, and both strict
  permutation deficits. Normal, optimized, and hash-seeded replays agree.
extension_independent_audit: >
  ACCEPT. A clean-room source-chart implementation uses the disjoint control
  Delta=-64/105, K=64, Phi=1, Theta=2, eta=3, reconstructs both resultants,
  derives the differential orders rather than importing the packet, and
  reproduces the cubic carrier and degree 25/19 contradictions. A separate
  symbolic face audit verifies that every non-top Newton face remains
  generically squarefree under (1a), including K=0 and Theta=0.
supplementary_eta_audit: >
  ACCEPT AS REDUNDANCY ONLY. Two disjoint eta-row implementations recover
  critical length 24, genus 10, the complete packet and the quadratic-carrier
  contradictions. They add no scope beyond the P/Y/B theorem.
---

# THM-4147 -- generic exact-weight-nine and anti-diagonal exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4141/4143
+ VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, INCLUDING THE GENERIC
ANTI-DIAGONAL EXTENSION; JC(2) REMAINS OPEN.**
Work over `C` in the live `b=d=0` reduced `(2,3)` seam.

## 1. The theorem and inheritance

> **Theorem.** None of the following three critical-open exact maximum
> residual-weight-nine chambers contains a nonautomorphic planar Keller pair:
>
> ```text
> P: zeta=0, eta*Delta*K*Theta!=0;
> Y: eta=0, zeta*Delta!=0;
> B: eta*zeta*Delta*(eta+zeta)!=0.                     (1)
> ```
>
> In every row the critical-open condition of Section 3 is also imposed.
> The coefficient `Phi` is unrestricted.

Section 9 proves a fourth dense-open exclusion on the repeated top-edge wall:

```text
A: zeta=-eta,
   eta*Delta*(Delta+Theta)*(4Theta K^2-27eta^2)!=0.     (1a)
```

It has its own exact critical-open condition. The exceptional divisors and
critical discriminant inside that wall remain outside the theorem.

The inheritance pass is deliberately label-sensitive:

- the closest proved mechanism is THM-4141/4143's exact critical-length
  versus labelled boundary-response obstruction;
- the canonical hostile is `eta+zeta=0`, where the top edge repeats and
  the critical strict transform changes;
- the corrected near miss is MISTAKE-509: an arbitrary geometric subset of
  boundary points is not a response because its residue field and horizontal
  image have been forgotten;
- the least-used sidecar is the pair `(L,beta)`: affine critical fixed-sheet
  supply together with total finite-carrier permutation index.

THM-4053 is related but is not used outside its exact-`M=8` scope. Source
completeness at weight nine comes instead from THM-3992/3997, THM-4007, and
the explicit weight enumeration below.

## 2. Exact source completeness

Put

```text
s=XT,                 p=T+s^2,                 y=sp,
P=T+X^2T^2,           Y=XTP,                   t=p-s^2,
K=2848/45-(7/6)Delta.                             (2)
```

THM-3992/3997 gives the normalized source form

```text
G=-s^2/(2t)+H(p,y),       H=lambda*p+R(p,y),
R in (p^2,y),             [y]R=[py]R=0.             (3)
```

For a monomial `p^i y^j`, the residual weight is `2i+3j`. Enumerating
`0<2i+3j<=9` and deleting exactly the two forbidden rows `y,py` leaves

```text
p, p^2, p^3, y^2, p^2*y, p^4, p*y^2, p^3*y, y^3.    (4)
```

THM-4007 supplies the forced lower row, including the displayed relation for
`K`. In the fixed gauge the complete exact-`M=9` polynomial is therefore

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4+Theta P Y^2
  +eta P^3Y+zeta Y^3.                                  (5)
```

There is no ellipsis in `(5)`. Exact weight nine says
`(eta,zeta)!=(0,0)`; the only new terms beyond weight eight are
`eta p^3y+zeta y^3`.

For a generic pencil value put `Q=q^-1`. Clearing `t` gives

```text
F_Q(s,p)=(s^2-p)(1-QH)-Q s^2/2.                        (6)
```

The primary certificate uses the collapsed **valued lower support** of
`(6)`, not its literal three-dimensional expansion. The term `-Qs^2/2`
has the same `(s,p)` exponent as the valuation-zero term `s^2`; its
height-one lift `(2,0,1)` is dominated and disappears when the common
coefficient's `Q`-adic valuation is taken. The independent audit expands
`(6)` first, combines every coincident coefficient, and obtains the same
lower hull.

## 3. The exact critical-open condition and affine length

Set

```text
R(T)=Res_X(G_X/T,G_T).                                  (7)
```

The **critical-open condition** in a chamber means that the appropriate row
below holds, both endpoints and the value at `T=-1/6` of the displayed
residual are nonzero, and its `T`-discriminant is nonzero:

```text
P: R=T^56(6T+1)^2 Q_20(T);
Y: R=T^56(6T+1)^2 Q_20(T);
B: R=T^56(6T+1)^2 Q_21(T).                             (8)
```

This is a Zariski-open condition. It is nonempty: the primary exact controls
use

```text
Delta=1, K=5591/90, Phi=11/7, Theta=19/11,
eta=23/13, zeta=29/17                                  (9)
```

with the absent top coefficient set to zero in `P,Y`. A disjoint
independent control uses

```text
Delta=-64/105, K=64, Phi=1, Theta=2,
(eta,zeta)=(1,0),(0,2),(1,2).                          (10)
```

Symbolically on the corresponding rows, the residual endpoint formulas are

```text
P: Q_20(0)=-(3^15/2^7)eta^7,
   [T^20]Q_20=72900 eta^5 K^4 Theta^4;
Y: Q_20(0)=-(3^15/2^7)zeta^7;
B: Q_21(0)=-(3^15/2^7)(eta+zeta)^7.                   (11)
```

For `T!=0`, both critical polynomials have common leading `X`-row
`9(eta+zeta)T^8`, nonzero in all three chambers, so no critical point is
lost at `X=infinity`. The exponent `56` in `(8)` is a Sylvester
degree-drop artifact:

```text
(G_X/T)(X,0)=-X,              G_T(X,0)=-(X^2+6)/2.     (12)
```

The actual critical ideal restores

```text
T=0,      X^2=-6,      G=0,                            (13)
```

while `(6T+1)^2` is the other universal pair

```text
T=-1/6,   X^2=6,       G=1/2.                          (14)
```

At the controls, an independently formed rational `(s,p)` critical pair gives

```text
P: Res_s(A,C)=p^4 R_20(p);
Y: Res_s(A,C)=p^8 R_20(p);
B: Res_s(A,C)=p^8 R_21(p),                             (15)
```

with nonzero endpoints, nonzero value at `T=-1/6`, and squarefree residuals
at the controls. Its
collapsed `p=0` and `t=0` strata restore the same two universal pairs.
Thus the critical-open affine lengths are exactly

```text
P: L=24,                 Y: L=24,                 B: L=25. (16)
```

For a hypothetical Keller realization `G=E(A,C)`, the Hessian congruence
at the two target nodes makes every affine critical point Morse. If `r_i`
counts affine inverse images of the node with value `i/2`, then

```text
r_0+r_1=L.                                               (17)
```

## 4. Newton polygons and complete labelled packets

A term `p^i y^j` in `H` contributes valued endpoints

```text
(j+2,i+j,1),                 (j,i+j+1,1).               (18)
```

Exact lower-hull enumeration gives:

| chamber | Newton polygon, counterclockwise | `(2Area,B,g_Pick)` | infinity packet |
|---|---|---:|---|
| `P` | `(0,1),(2,0),(4,2),(4,3),(3,4),(1,5),(0,5)` | `(29,11,10)` | `(8,5,4,3,2,2,1)` |
| `Y` | `(0,1),(2,0),(5,3),(3,4),(0,5)` | `(30,10,11)` | `(11,8,2,2,2,1)` |
| `B` | `(0,1),(2,0),(5,3),(1,5),(0,5)` | `(31,11,11)` | `(8,8,4,2,2,2,1)` |

For an oriented edge with primitive inward normal `(u,v)` and supporting
constant `c`, THM-4103's Keller residue identity gives

```text
ord(omega)=u+v-c-1,              e=u+v-c.               (19)
```

The length-four vertical edge `s=0,p!=0` consists of affine source points,
not source punctures, and is excluded from the packets. The only
nonrational edge is

```text
P:   q-1/2=K W^2;
Y,B: q-1/2=K W^2+zeta W^3.                             (20)
```

The first equation is irreducible because `q-1/2` has odd valuation at
`q=1/2`, so it is one separable quadratic closed point over `C(q)`. In the
second equation, `C(W)/C(q)` has degree three because the polynomial map
`W |-> 1/2+K W^2+zeta W^3` has degree three. Thus it is one separable cubic
closed point. In chamber `B`, the
length-two top edge is

```text
-(S^2-P)(eta P+zeta S^2),                               (21)
```

whose constant roots `S^2/P=1,-eta/zeta` are distinct exactly when
`eta+zeta!=0`.

The packet defects are `18,20,20`, equal to `2g_Pick-2`. The finite
critical scheme and THM-3827's closed-polynomial factor theorem imply that
the geometric generic source is connected: a factorization
`G=R(G_0)`, `deg R>1`, would make a complete curve lie in the finite
critical scheme. Its normalization genus is at most the Pick genus, whereas
Riemann--Hurwitz over the elliptic target forces it to be at least that genus
from the displayed boundary defect. Equality follows. Hence the genera and
packets in the table are complete: there is no hidden affine ramification or
unrecorded genus loss.

## 5. Prime carrier response

THM-4120 proves the target-only statement

```text
E_q(C(q))={O}.                                           (22)
```

THM-4122 supplies the affine-line normalization of a horizontal
nonproperness component. Let `K_0=C(q)`, and let `L/K_0` be the residue extension of the
nonrational boundary point; its prime degree is `ell=2` in `P` and
`ell=3` in `Y,B`. If its image is finite, the residue field `M` of
that horizontal target image is intermediate:

```text
K_0 subset M subset L.                                  (23)
```

Prime degree leaves `M=K_0` or `M=L`. The first alternative would be a
finite `K_0`-point of `E_q`, forbidden by `(22)`. Thus `M=L`: the
boundary maps birationally to one degree-`ell` horizontal carrier.
Separability gives exactly `ell` distinct conjugate target points, each
with a transposition meridian. Consequently

| chamber | finite response `(n,beta)` | full response `n` |
|---|---:|---:|
| `P` | `(21,2)` | `25` |
| `Y` | `(20,3)` | `26` |
| `B` | `(21,3)` | `27` |

Here `beta=ell` is the total permutation index of the finite carrier
meridians. Polynomiality prevents an affine source point from mapping to the
projective target origin. Hence the full degree is the sum of the complete
infinity packet, while the finite degree removes the `ell` index-two
punctures from the origin fibre.

## 6. Finite-separable-carrier transport

The cubic case cannot import THM-4138's quadratic conclusion by analogy.
The needed statement is the following source-independent strengthening of
its local/proper mechanism.

> **Finite-separable-carrier transport lemma.** Let
> `varphi_q:C_q->E_q` be the induced map between normalized generic fibres
> of the inherited elliptic pencil. Assume that `C_q` is geometrically
> connected, the displayed boundary packet is the complete ramification
> divisor, and the finite branch locus away from the origin is one separable
> horizontal carrier of degree `ell`. If `r_i` is the number of affine
> Keller-Morse inverse images of target node `o_i`, then there is a smooth
> reference fibre with handle loops `delta_0,delta_1` and carrier meridians
> `mu_1,...,mu_ell` such that, for `X=rho(delta_0)` and
> `Y=rho(delta_1)`,
>
> ```text
> #Fix(X)>=r_0,              #Fix(Y)>=r_1.              (24)
> ```
>
> Every `rho(mu_j)` is a transposition, and `X,Y` together with these
> carrier meridians generate a transitive action on the sheets.

**Proof.** First compactify the pencil graph, resolve it, and take the
relative normalization. The generic map is nonconstant between smooth
projective curves and hence finite. Shrink the `q)-line so that both
resolved families are proper and smooth, the resolved morphism is fibrewise
quasifinite, and therefore the proper morphism is finite. Delete the finite
set where a family is singular, the map changes degree or has a vertical
exceptional component, the carrier is singular, ramified over the base,
self-intersecting or incident with the origin, or the complete ramification
packet specializes or acquires an extra Hurwitz collision.

For the cubic in `(20)`,

```text
Disc_W(KW^2+zeta W^3-(q-1/2))
 =(q-1/2)(4K^3-27zeta^2(q-1/2)).                       (25)
```

Thus the deletion explicitly contains `q=1/2`--a double collision, triple
when `K=0`--and, for `K!=0`,
`q=1/2+4K^3/(27zeta^2)`. This changes no generic degree.

On the remaining open base `U`, the reduced carrier `S_U->U` is finite
etale of degree `ell` and is disjoint from the origin section `O_U`.
Completeness says that every relative ramification point lies over
`D=O_U union S_U`. Only now delete the complete marked divisor and its
full inverse image:

```text
varphi^circ:
 Cbar_U-varphi^(-1)(D) -> Ebar_U-D.                     (25a)
```

This restriction is finite etale. The proper-smooth construction precedes
divisor deletion; the punctured families themselves are not proper.

At an affine inverse `z` of target node `o_i`, target Morse coordinates
give `E-q_i=uv`. Since the Keller map is a local biholomorphism at `z`,
pullback gives the identical source equation. A Milnor core in a nearby
smooth fibre therefore has one closed degree-one lift in the inverse
neighborhood of `z`. The carrier meets the Milnor annulus in finitely many
points. A parallel core avoids all of them even when carrier branches collide
in the central fibre, because the core lies in a nearby smooth fibre.
Pairwise disjoint inverse neighborhoods give distinct closed lifts.

Choose a contractible tree in `U` joining nearby punctured neighborhoods of
the two nodal values to one reference value. Over it the proper marked family
and the finite-etale cover `(25a)` trivialize. Transport preserves all
closed lifts; sheet labels change only by simultaneous conjugation. Choose
the inherited distinguished paths first and transport both parallel Milnor
cores in this one marked trivialization. Small simultaneous pushes avoid the
carrier throughout, preserve their intersection number one, and therefore
produce the same handle generators used below rather than unrelated loop
classes. This proves `(24)`.

The inherited target vanishing cycles from THM-4130 have intersection number
one. Simultaneous small parallel pushes represent them by loops avoiding
`D` and intersecting once. A thin regular neighborhood is a once-punctured
torus whose complementary disk contains `O` and all carrier points. Hence

```text
pi_1(E_q-{O,Q_1,...,Q_ell})
 = <delta_0,delta_1,mu_O,mu_1,...,mu_ell |
    [delta_0,delta_1] mu_O mu_1 ... mu_ell=1>.          (25b)
```

Eliminating `mu_O` proves the generation claim. Each carrier point has one
simple index-two boundary preimage, so its meridian is a transposition.
Finally, deleting finitely many points from the connected projective source
does not disconnect its Riemann surface. Thus the punctured cover is
connected and its monodromy action is transitive. QED.

For a full response there is no finite carrier; the same proof applies with
`ell=0` and `D=O_U`.

Combining `(17)` and `(24)` gives

```text
|supp X|+|supp Y|<=2n-L.                                (26)
```

## 7. The two permutation gates

### 7.1 Finite response

For every nonidentity permutation `sigma`,

```text
ind(sigma)<=|supp(sigma)|-1.                            (27)
```

A generator can merge at most its index many current orbit blocks; hence a
transitive generating system on `n` letters has total index at least
`n-1`. If at least one handle is nonidentity, `(26),(27)` bound the
total index of the handles and carrier meridians by

```text
2n-L-1+beta.                                            (28)
```

Every row has `L>n+beta`, so `(28)<n-1`. If both handles are identities,
the total index is only `beta<n-1`. Explicitly,

```text
P: 19<20,                 Y: 18<19,                 B: 19<20. (29)
```

Thus every finite response is impossible.

### 7.2 Full response

For arbitrary permutations `X,Y`, put

```text
k=|supp(X) intersect supp(Y)|.
```

Then

```text
ind([X,Y])<=2k.                                         (30)
```

Indeed, let `D=supp(Y)` and restrict `X` to `D`. It is the identity
on `D-supp(X)` and has only `k` nonidentity arrows. The directed
path/cycle decomposition of this partial injection can be closed to a
permutation `Xhat` with `ind(Xhat)<=k` and

```text
Xhat Y Xhat^-1=X Y X^-1.
```

Triangle inequality for Cayley transposition distance and conjugacy
invariance give

```text
ind([X,Y])
 <=ind(Xhat)+ind(Y Xhat^-1 Y^-1)
 =2ind(Xhat)<=2k.                                       (31)
```

In a full response, `X,Y` alone generate transitively, so their supports
cover all `n` sheets. Equation `(26)` yields

```text
k=|supp X|+|supp Y|-n<=n-L.                             (32)
```

The origin meridian is `[X,Y]^-1`; completeness of the boundary packet
makes its permutation index equal to the packet defect. The contradictions
are

```text
P: 18<=ind([X,Y])<=2(25-24)=2;
Y: 20<=ind([X,Y])<=2(26-24)=4;
B: 20<=ind([X,Y])<=2(27-25)=4.                         (33)
```

Thus every full response is impossible, proving the theorem.

## 8. Hostiles, scope, and exact audit

The support cancellations `K=-1376/135`, `Theta=Delta`, and
`eta=zeta` were tested exactly and do not change their admitted polygons;
`eta=zeta` retains `T^56(6T+1)^2Q_21`. In contrast,
`eta+zeta=0` gives at the hostile control

```text
R=T^42(6T+1)^2Q_19(T),                                  (34)
```

so the repeated top edge is a genuine strict-transform wall. Sections 1--7
do not cross it; Section 9 treats its dense open subchamber `A` separately.

The theorem does not cross the remaining coefficient contractions excluded
in `(1),(1a)`, a zero leading endpoint, `Q_d(-1/6)=0`, a residual
discriminant, another reduced cell, entry into this seam, exact weight at
least ten, `JC(2)`, or `DC(2)`.

The primary certificate performs `15,104` exact checks. Normal, optimized,
and fixed-hash-seed runs byte-match. The clean-room audit reconstructs the
source polynomial before valuation, uses the disjoint control `(10)`, and
checks all `533,417` ordered pairs in `S_1,...,S_6`; `22,818` nonzero
pairs attain equality in `(30)`.

The primary semantic digest is

```text
e88dc5ce513e71cf18701434e3c5c14d9e5e22a06655cbc05e89580da4df1a17.
```

## 9. Generic anti-diagonal extension

On `A`, the exact top row contracts as

```text
eta(P^3Y-Y^3)=eta X T^2 P^3.                          (35)
```

Let `Q_19` be the residual in the `(X,T)` critical projection. Add to
`(1a)` the critical-open conditions

```text
Disc_T(Q_19)!=0,                    Q_19(-1/6)!=0.     (36)
```

The endpoint identities are

```text
Res_source=p^6 R_19(p),
R_19(0)=46656 eta(4Theta K^2-27eta^2),
[p^19]R_19=1327104 eta^5(Delta+Theta)^4;

Res_X=T^42(6T+1)^2 Q_19(T),
Q_19(0)=-12288(Delta+Theta)^6,
[T^19]Q_19=-1458(Delta+Theta)eta^4
                  (4Theta K^2-27eta^2)^2.             (37)
```

Thus `(1a),(36)` give nineteen reduced residual critical points, while the
two Morse pairs at `T=0` and the two at `T=-1/6` remain separated. Hence

```text
L_A=19+2+2=23.                                        (38)
```

The exact rational control

```text
(Delta,Theta,Phi,eta,zeta)=(1,19/11,11/7,23/13,-23/13)
```

has squarefree `Q_19` with nonzero value at `-1/6`, so this open chamber is
nonempty.

Before normalization the Newton polygon is the `B` polygon and has Pick
genus eleven. At its repeated top point, put

```text
s=z^-1,                         p=(1-a)z^-2.
```

If `L(a,z)=z^11 F_Q(z^-1,(1-a)z^-2)`, its tangent cone is

```text
Q a(eta a-(Delta+Theta)z).                             (39)
```

The two factors are distinct on `(1a)`, so this is one ordinary node with
branches `a~0` and `a~((Delta+Theta)/eta)z`. Moreover

```text
omega=Q ds/(F_Q)_p = Q z^7 dz/L_a.                    (40)
```

On each branch `L_a` has order one in `z`; therefore `omega` has order six
and both normalized branches have ramification index seven. Replacing the
two raw index-eight branches gives the candidate packet

```text
(7,7,4,2,2,2,1),                 defect=18.           (41)
```

The finite critical scheme and THM-3827 first give geometric connectedness,
exactly as in Section 4. The ordinary node gives normalization genus at most
ten. Conversely the displayed defect and THM-4103's Keller residue identity
give, by Riemann--Hurwitz over the elliptic target, `2g-2>=18`, so `g>=10`.
Consequently `g=10`, `(41)` is complete, and there is no hidden ramification.

The nonrational boundary equation is the separable cubic

```text
-eta W^3+K W^2=q-1/2.                                 (42)
```

For completeness, a direct face audit gives the other non-top equations

```text
2q+(1-2q)W,                 eta+Delta W,
H_0(p)-q,
H_0(p)=Delta p^4-(1376/135)p^3+(8/3)p^2-3p.           (42a)
```

The first two are squarefree torus labels when `eta*Delta!=0`; the vertical
quartic is generically squarefree over `C(q)` because
`gcd(H_0-q,H_0')=1`.  The cubic in `(42)` has discriminant

```text
(q-1/2)(4K^3-27eta^2(q-1/2)),                         (42b)
```

so it remains generically irreducible and separable even when `K=0`.
Likewise `Theta=0` changes no non-top face.  The exceptional values displayed
in `(42b)` and the finitely many critical values of `H_0` are deletions in the
pencil base, not new coefficient gates.

THM-4120/4122 and Section 6's finite-separable-carrier lemma apply without
change. If the cubic point maps to the origin, `n=25`; otherwise its three
conjugate index-two branches give `n=19,beta=3`. In the full response,

```text
ind([X,Y])<=2(25-23)=4<18.                            (43)
```

In the finite response, the maximum merger capacity is

```text
2*19-23-1+3=17<18=n-1.                               (44)
```

Both responses are impossible. This proves the generic anti-diagonal
extension. The divisors omitted in `(1a),(36)`, including the exact drop
strata recorded by the certificate, remain open. **QED.**

## 10. Supplementary eta-row cross-check

The two `eta_only_supplementary_crosscheck` artifacts give a disjoint exact
audit of row `P`; they add no theorem scope. On the narrower open locus where
the saturated critical scheme has twenty points and the `T=-1/6` fibre is
exactly `X^2=6`, the witness

```text
(Delta,Phi,Theta,eta)=(2,5,7,11),       K=2743/45
```

reconstructs `(g,L,packet)=(10,24,(8,5,4,3,2,2,1))`. The quadratic carrier
can then use THM-4138 directly. A one-point support intersection makes the
degree-25 commutator a three-cycle, while degree 21 has merger capacity at
most nineteen. These independent mechanisms corroborate row `P` but are not
dependencies for rows `Y,B,A`.

Replay with:

```bash
python3 -B 04-computation/jc23_weight9_generic_monodromy_thm4147.py
python3 -B -O 04-computation/jc23_weight9_generic_monodromy_thm4147.py
PYTHONHASHSEED=271828 python3 -B 04-computation/jc23_weight9_generic_monodromy_thm4147.py
python3 -B 04-computation/jc23_weight9_generic_monodromy_thm4147_independent_audit.py
python3 -B -O 04-computation/jc23_weight9_generic_monodromy_thm4147_independent_audit.py
PYTHONHASHSEED=271828 python3 -B 04-computation/jc23_weight9_generic_monodromy_thm4147_independent_audit.py
python3 -B 04-computation/jc23_exact_weight_nine_generic_antidiagonal_exclusion_thm4147.py
python3 -B 04-computation/jc23_weight9_antidiagonal_generic_independent_audit_thm4147.py
python3 -B 04-computation/jc23_weight9_eta_only_supplementary_crosscheck_thm4147.py
python3 -B 04-computation/jc23_weight9_eta_only_supplementary_crosscheck_thm4147_independent_audit.py
```
