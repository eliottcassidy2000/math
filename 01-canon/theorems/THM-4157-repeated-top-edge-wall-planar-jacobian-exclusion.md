---
id: THM-4157
title: "Repeated top-edge wall planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147
  + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the exact-weight-nine
  repeated top-edge wall zeta=-eta with eta*Delta!=0, the four exhaustive
  structural rows have packets (7,7,4,2,2,2,1), (6,6,4,2,2,2,1),
  (5,5,4,2,2,2,1), and (7,4,2,2,2,1). Nonempty critical opens in the first
  three rows contain no nonautomorphic planar Keller pair, and the deepest
  row is excluded for every eta!=0. The latter closure combines projection
  switching, an exact quadratic endpoint wall, and a degree-thirty
  nonreduced-critical-point/Hessian obstruction. The first-three-row walls
  left outside this theorem are subsequently closed by THM-4173/4176.
  Eta=0, Delta=0, entry, other cells, exact weight at least ten, JC(2), and
  DC(2) remain OPEN.
source: codex-frontier-synthesis-creative-20260825aq
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
related:
  - THM-4141-delta-d-collision-wall-boundary-monodromy-exclusion
  - THM-4143-two-term-collision-wall-critical-boundary-monodromy-exclusion
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4171-row-a-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4173-repeated-top-row-a-complete-planar-jacobian-exclusion
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
script: 04-computation/jc23_weight9_repeated_top_wall_thm4157.py
output: 05-knowledge/results/jc23_weight9_repeated_top_wall_thm4157.out
symbolic_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_symbolic_endpoints_thm4157.out
factor_script: 04-computation/jc23_weight9_repeated_top_wall_factors_thm4157.py
factor_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_factors_thm4157.out
j_wall_script: 04-computation/jc23_weight9_repeated_top_wall_J_wall_thm4157.py
j_wall_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_J_wall_thm4157.out
h30_script: 04-computation/jc23_weight9_repeated_top_wall_H30_thm4157.py
h30_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_H30_thm4157.out
independent_newton_audit_script: 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_newton_audit.py
independent_newton_audit_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_thm4157_independent_newton_audit.out
independent_quotient_audit_script: 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_quotient_audit.py
independent_quotient_audit_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_thm4157_independent_quotient_audit.out
independent_projection_audit_script: 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_projection_audit.py
independent_projection_audit_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_thm4157_independent_projection_audit.out
independent_j_wall_audit_script: 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_J_wall_audit.py
independent_j_wall_audit_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_thm4157_independent_J_wall_audit.out
independent_h30_point_audit_script: 04-computation/jc23_weight9_repeated_top_wall_H30_thm4157_independent_point_audit.py
independent_h30_point_audit_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_H30_thm4157_independent_point_audit.out
independent_h30_tower_audit_script: 04-computation/jc23_weight9_repeated_top_wall_H30_thm4157_independent_tower_audit.py
independent_h30_tower_stage1_output: 05-knowledge/results/jc23_weight9_repeated_top_wall_H30_thm4157_independent_tower_stage1_audit.out
script_sha256: a69244d03e52aa33f843696d58c505fda0f7d0618d4fa0d6515eead728221927
output_sha256: 39a645fcaf3756d9eb0838557c81d020994163b91d5632f852b12385f4eab485
symbolic_output_sha256: be28a7806f09c656f6ad7e7f79ff2d7566858cda5b5a5ed90891ffa8a131836e
factor_script_sha256: 079852fd424591858c3560150ba5250233d299cfe2f2285050641642dbf2f513
factor_output_sha256: 91b91809c99ebeecf26af27fdb7168311fa5124a90f71a8f942c1495133cf483
j_wall_script_sha256: dd4a6ead99d3013b260913114561ef1027ef0273d642502492b15def2f87bb90
j_wall_output_sha256: 1a577bdb1f44cc2e2f27927ca01b19527e9e1603767aba54509baf8495d99d10
h30_script_sha256: 90d7715649fd7d7f55cfc11c87429311a843a929e3f53820a785bfbb70f70073
h30_output_sha256: ea0e093bfeb9abeb9b6878dd23751e51922f3bfa06231815944f445875823f1f
independent_newton_audit_script_sha256: 191919f9d0fd470e6e450ab024618a6aa9d22782d063d656c6174041717d877d
independent_newton_audit_output_sha256: 80286b92a33fa71fb2026500c07d5f66dc47a6aad6f8ca726e35947ead47335d
independent_quotient_audit_script_sha256: d9d58785ac4953893fe79f86117e690a9862dc528304f2aaa655c0e6c563fb71
independent_quotient_audit_output_sha256: 71f647c17088e620d9c61a92dc57032e434cea75d74b155763c0e223bf7482ca
independent_projection_audit_script_sha256: 86e4515bbdd78612787ad972187e7a9aa677df7e9af609cce9d79a2910ad4b76
independent_projection_audit_output_sha256: 63367d33e0c740536acae3941de5f5052fce805df0d9972750e86fd275f89966
independent_j_wall_audit_script_sha256: a0d093e95fabf6f5f5fa9948de2539232dadf2e20d0cf66ad7fae053f80a1e49
independent_j_wall_audit_output_sha256: af66df80cd2bd5bcfcbd2ba8285151afbb7d17f14f4de0126ba1a608d237bc56
independent_h30_point_audit_script_sha256: 754555871dbd082bd21a912b4216bf55ce8e389037bff2887ac2445da13be0c7
independent_h30_point_audit_output_sha256: be82d02a296cf6dd40bcb7daf8278fbb49212b58f8de0771fdbc237d92779843
independent_h30_tower_audit_script_sha256: 248dc4f3b727a602270d67a2457c701a1442ced5c567b60ac7c15b96b3c400ab
independent_h30_tower_stage1_output_sha256: 760f9556a740874fcf9708a609972ef4f251eca6c9406c49b3c2c841b416c9e1
semantic_sha256: 6267e49f6c6e5edecb2ae99bcb88f6111010f2b02663c84652ad30654334af2d
factor_semantic_sha256: 08f3d689c392ea14d3d87df61068cb35267056d2a2f8d9318617cc568657d2f8
j_wall_semantic_sha256: b373e61184f208ce7315f655d74834fd0182c8614e381766cca1d41b186ef908
h30_semantic_sha256: f400db7367d8ff42e6bddebe74a5952ed2af8f21999bc57202c68dc2aa77fe94
independent_h30_point_semantic_sha256: 3534b84eac3679630eeaed120f68478cc2f6a16248620323957a80832b1995e0
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact symbolic and rational-control certificates reconstruct the
  complete wall source, all four strict transforms, Newton packets, two
  critical projections, universal critical points, carrier sidecars, and
  finite/full monodromy deficits. Companion certificates factor the deepest
  discriminants, close J=0 over Q(sqrt(-30)), and close H_30=0 by an exact
  linear subresultant, local intersection multiplicity, and Hessian identity.
  Normal, optimized, and hash-seeded executions byte-match their frozen
  outputs.
independent_audit: >
  ACCEPT. Four clean-room implementations independently recover the Newton
  hierarchy, quotient critical lengths, a primitive alternative projection,
  the common H_30 factor, and the entire J=0 response. A separate hostile
  audit proves H_30 irreducible modulo 347 and verifies the field tower,
  exact double root, unit/infinity gates, linear subresultant, local
  intersection length, Hessian identity, and point-level zero-Hessian
  contradiction. No hidden saturation or projection loss remains.
---

# THM-4157 -- repeated top-edge wall planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147
+ VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; JC(2) REMAINS OPEN.**
Work over `C` in the live `b=d=0` reduced `(2,3)` seam.

## 1. Decision

On the exact-weight-nine repeated top-edge wall

```text
zeta=-eta,                    eta!=0,                    (1)
```

put

```text
C=Delta+Theta,
B=epsilon+K,       epsilon=-1376/135,
K=2848/45-(7/6)Delta.                                 (2)
```

Assume `Delta!=0`, so the common Newton polygon below does not contract.  The
wall has four exhaustive structural strata:

```text
A: C!=0;
B: C=0, Phi!=0;
C: C=Phi=0, B!=0;
D: C=Phi=B=0.                                         (3)
```

In row `D`, equation `(2)` forces

```text
Delta=2048/45,       Theta=-2048/45,       K=1376/135. (4)
```

On a nonempty Zariski-open critical locus in **each** row, no nonautomorphic
planar Keller pair exists.  The exact normalized packets, affine critical
lengths, and responses are:

| row | normalized packet | `g`, defect | `L` | finite `(n,beta)` | full `n` |
|---|---|---:|---:|---:|---:|
| `A` | `(7,7,4,2,2,2,1)` | `10,18` | `23` | `(19,3)` | `25` |
| `B` | `(6,6,4,2,2,2,1)` | `9,16` | `21` | `(17,3)` | `23` |
| `C` | `(5,5,4,2,2,2,1)` | `8,14` | `19` | `(15,3)` | `21` |
| `D`, `J*H_30!=0` | `(7,4,2,2,2,1)` | `7,12` | `16` | `(12,3)` | `18` |
| `D`, `J=0` | `(7,4,2,2,2,1)` | `7,12` | `15` | `(12,3)` | `18` |

Rows `A,B,C` and the `J*H_30!=0` part of row `D` have the same two strict
margins:

```text
finite: L=n+beta+1;
full:   n-L=2, hence ind([X,Y])<=4<packet defect.       (5)
```

Thus the commutator-overlap gate, together with the prime cubic carrier
sidecar, excludes all four nonempty critical-open subchambers.  On `J=0`, the
critical length drops to `15`; the sharpened finite/full inequalities in
section 5.1 still exclude the slice.  On `H_30=0`, section 5.2 constructs a
nonreduced affine critical point and hence a zero Hessian, contradicting the
Keller--Morse congruence.  Thus row `D` is completely excluded in scope.  The
row-`A,B,C` coefficient/resultant/discriminant walls outside this theorem are
listed in section 9 and are subsequently closed by THM-4173/4176.

## 2. Complete source and type ledger

Source completeness is inherited from the reduced `(2,3)` construction of
THM-3992 (`THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual`),
the Hasse coefficient repair of THM-3997
(`THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go`), and the
third source-normal row of THM-4007
(`THM-4007-live-two-three-third-normal-row-five-weight-floor`).  This route does
**not** use THM-4053 as a source-completeness assertion.  THM-4147
(`THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion`)
canonically performs the resulting explicit weight-nine monomial enumeration;
the repeated-wall contraction and every strict-transform/critical calculation
below are independent of THM-4147's generic-chamber conclusion.

The canonical independent computation
`04-computation/jc23_exact_weight_nine_generic_antidiagonal_exclusion_thm4147.py`
(commit `19ce3486d3`) is the primary row-`A` cross-check.  Its independent
`(s,p)` and `(X,T)` eliminations give the same packet
`(7,7,4,2,2,2,1)`, `L=23`, finite `(n,beta)=(19,3)`, full `n=25`, and exact
endpoints.  The four-row primary calculation below retains row `A` only as an
independent replay; rows `B,C,D` and the two row-`D` wall closures are not
claimed from that canonical script.

After the scalar normalization used in THM-4130, set

```text
P=T+X^2T^2,                         Y=XTP,
K=2848/45-(7/6)Delta.                                 (6)
```

At exact residual weight nine the complete normalized source polynomial is

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta P^4+Theta P Y^2
  +eta P^3Y+zeta Y^3.                                 (7)
```

The only weight-nine residual monomials are the last two.  No higher-weight
tail is being silently truncated into this statement: after the inherited
lower coefficients are inserted, direct enumeration of residual monomials
`p^i y^j` with `2i+3j<=9` leaves precisely `p^3y,y^3` on the new top row.
In the rational chart

```text
s=XT,       p=T+s^2,       t=p-s^2,       y=sp,         (8)
```

wall `(1)` has the exact contraction

```text
eta p^3y-eta y^3=eta s p^3(p-s^2)=eta s p^3t.          (9)
```

The typed objects used later are:

```text
source:  normalization of the generic fibre G=q, plus its affine critical
         fibres over the two nodal target values;
target:  the generic elliptic target fibre E=q over C(q);
map:     restriction of the hypothetical Keller map (A,C), followed by
         compactification/normalization and finite-separable transport;
data kept: residue-field degree, ramification index, affine fixed-sheet
         supply L, finite-carrier index beta, and origin commutator defect;
data lost: carrier coordinates, ordering of conjugates, and coefficients
         within a connected critical-open chamber.                    (10)
```

## 3. Common Newton polygon and exact strict transform

Write

```text
gamma=-1/2, lambda=-3, alpha=8/3, epsilon=-1376/135,
H=lambda p+alpha p^2+epsilon p^3+K s^2p^2+Phi sp^3
  +Delta p^4+Theta s^2p^3+eta sp^3(p-s^2).             (11)
```

For `Q=q^-1`, the source fibre is

```text
F_Q=(s^2-p)(1-QH)+gamma Qs^2.                          (12)
```

Exact support collection gives in every row of `(3)`

```text
Newt(F_Q)=conv{(0,1),(2,0),(5,3),(1,5),(0,5)},
2Area=31,       boundary=11,       g_arith=11.          (13)
```

The length-two edge `(5,3)--(1,5)` has repeated face root.  Use its actual
toric chart

```text
s=1/z,       p=(1-a)/z^2,       r=1-a,
L(a,z)=z^11 F_Q(1/z,r/z^2).                             (14)
```

The exact strict transform is

```text
L=a z^9+gamma Qz^9+Qeta a^2r^3
  -Qa r^3(Delta r+Theta)z-QPhi a r^3z^2
  -Qa(epsilon r^3+Kr^2)z^3-Qalpha ar^2z^5
  -Qlambda arz^7.                                      (15)
```

Since `Delta r+Theta=C-Delta a`, successive vanishing of `C,Phi,B`
produces the complete local hierarchy:

| row | first Newton balance | normalization branches | local `delta` | top indices |
|---|---|---|---:|---|
| `A` | `Qa(eta a-Cz)` | `a=(C/eta)z+...`, `a=-z^8/(2C)+...` | `1` | `(7,7)` |
| `B` | `Qa(eta a-Phi z^2)` then `-QPhi az^2-Qz^9/2` | `a=(Phi/eta)z^2+...`, `a=-z^7/(2Phi)+...` | `2` | `(6,6)` |
| `C` | `Qa(eta a-Bz^3)` then `-QBaz^3-Qz^9/2` | `a=(B/eta)z^3+...`, `a=-z^6/(2B)+...` | `3` | `(5,5)` |
| `D` | `Q(eta a^2-z^9/2)`, `gcd(2,9)=1` | `z=tau^2`, `a=c tau^9+...`, `c^2=1/(2eta)` | `4` | `(7)` |

For the first three rows the two smooth branches meet to order `1,2,3`.
For row `D`, the `(2,9)` plane branch has delta invariant four.  The Keller
residue identity is, up to a nonzero sign,

```text
omega=Q ds/(F_Q)_p=Q z^7 dz/L_a.                       (16)
```

In rows `A,B,C`, `ord_z L_a=1,2,3`, giving differential orders `6,5,4`
and indices `7,6,5` on both branches.  In row `D`, pulling `(16)` back by
`z=tau^2,a=c tau^9+...` gives order `14+1-9=6`, hence index seven.

The other infinity edges are squarefree and contribute

```text
rational indices 4 and 1,
one cubic closed point with three geometric indices 2,2,2.              (17)
```

Subtracting the local delta invariant from `g_arith=11` gives the upper
genus bounds `10,9,8,7`.  The packets in section 1 have defects
`18,16,14,12`, so Riemann--Hurwitz for the nonconstant map to the elliptic
target forces the matching lower bounds.  Equality follows.  Thus every
displayed packet is complete: there is no hidden genus loss or affine
ramification.

Geometric connectedness is supplied by the closed-polynomial factor argument
of THM-3827 after the finite critical scheme below is known.  A disconnected
generic fibre would make `G=R(G_0)`, and a root of `R'` would put a whole curve
inside `Crit(G)`, contradicting finiteness.

## 4. Exact critical-open loci and endpoint cascade

Put

```text
f=G_X/T,                         h=G_T,
Acrit=t^2 G_s/p,                 Ccrit=2t^2 G_p.        (18)
```

For row `i`, the **two-projection critical-open condition** means

```text
Res_X(f,h)=T^v(6T+1)^2 Q_d(T),
Res_s(Acrit,Ccrit)=p^8 R_d(p),                         (19)
```

with the displayed degrees, nonzero leading and constant coefficients,
`Q_d(-1/6)!=0`, and

```text
disc_T(Q_d)*disc_p(R_d)!=0.                            (20)
```

The exact generic rows are

| row | `(v,d)` in `(19)` | exact `Q_d(0)` |
|---|---:|---|
| `A` | `(42,19)` | `-12288 C^6` |
| `B` | `(30,17)` | `-(3*7^5/2^5) Phi^5` |
| `C` | `(32,15)` | `2*3^6 eta B^5` |
| `D` | `(35,12)` | `-(3^11/2^5) eta^5` |

These endpoint formulae explain, rather than merely detect, the four
structural strata.  Each endpoint vanishes exactly when the next local Newton
coefficient must be used; row `D` terminates because `eta!=0`.

Four exact rational controls prove that all four critical-open loci are
nonempty:

| row | `(Delta,Theta,Phi,eta,K,B)` |
|---|---|
| `A` | `(1,19/11,11/7,23/13,5591/90,14021/270)` |
| `B` | `(1,-1,11/7,23/13,5591/90,14021/270)` |
| `C` | `(1,-1,0,23/13,5591/90,14021/270)` |
| `D` | `(2048/45,-2048/45,0,23/13,1376/135,0)` |

At those controls, exact endpoint and squarefreeness data are:

| row | `Q_d(0)` | `[T^d]Q_d` | `R_d(0)` | `[p^d]R_d` |
|---|---|---|---|---|
| `A` | `-8957952000000/1771561` | `-11206045421019265531550064784/407151596119125` | `1325533905629568/604175` | `6918774197944320000/5436100813` |
| `B` | `-483153/32` | `85821754813159273972508/352960408125` | `-70369952239872/54925` | `30531977307612/371293` |
| `C` | `12463005381719088344323/12793950000` | `85821754813159273972508/352960408125` | `-70369952239872/54925` | `994981781335325982843932/18796708125` |
| `D` | `-1140178853421/11881376` | `84349672224497501284348663943/231577323770812500` | `-3875767847803456/2471625` | `-1809463840379127/62748517` |

In every row `gcd(Q_d,Q_d')=gcd(R_d,R_d')=1`, and `Q_d(-1/6)!=0`.
The frozen output records exact normalized coefficient digests.

The `T`-power in `(19)` is a Sylvester degree-drop artifact.  Direct
specialization always gives

```text
f(X,0)=-X,                    h(X,0)=-(X^2+6)/2.        (21)
```

The factor `(6T+1)^2` is the pair

```text
T=-1/6,       X^2=6,       G=1/2,       det Hess(G)=-6,             (22)
```

and the actual ideal `(Tf,h)` restores the omitted pair

```text
T=0,          X^2=-6,      G=0,         det Hess(G)=6.              (23)
```

In the independent chart, `(22)` lies at `p=0`, while `(23)` is omitted at
`t=0`.  Therefore the full affine critical lengths are exactly

```text
A: 19+2+2=23,   B: 17+2+2=21,
C: 15+2+2=19,   D: 12+2+2=16.                         (24)
```

For a Keller realization, Hessian congruence through the determinant-one map
makes every source critical point Morse.  If `r_0,r_1` are the inverse counts
over the two target nodes, then `r_0+r_1=L`.

## 5. Deepest-row projection portfolio

Row `D` has only the free coefficient `eta`.  Here the two projections can be
combined more sharply than condition `(20)`.  Put

```text
J(eta)=22143375 eta^2+15510536192.                     (25)
```

Exact symbolic elimination gives

```text
Q_12(0)=-(3^11/2^5)eta^5,
[T^12]Q_12=eta^3 J(eta)^2/3690562500,
R_12(0)=-(64/1125)eta J(eta),
[p^12]R_12=-531441 eta^7.                              (26)
```

The two discriminants factor over `Q[eta]` as

```text
disc_T Q_12 = unit * eta^32 * H_30(eta) * H_36(eta)^2,
disc_p R_12 = unit * eta^40 * H_30(eta) * H_24(eta)^2. (27)
```

The subscripts are degrees in `eta`; all three factors are even.  Their monic
coefficient digests are

```text
H_30: 7f65856a5fa69e3a35511e1c7eaf84f0c517397d9d08e4275965cbc5760f2e2b
H_36: 0311e867c9d077d409761a0f47c84967d304479650a2d6c6b9b532b198fcec15
H_24: 21267b85a00cd6ea8f70e78df570ee669d258c194bd0a36dda3c16b83f862341
```

For direct wall calculations, the companion extractor

```text
04-computation/jc23_weight9_repeated_top_wall_factors_thm4157.py
```

reconstructs both discriminants from the source, verifies the common
degree-thirty factor, substitutes `u=eta^2`, and prints the complete primitive
integer and monic rational coefficient lists in descending powers of `u`.
The corresponding monic-`u` digests (distinct from the monic-`eta` convention
above) are

```text
H_30(u): c26b5e4c5c390de1fe0d9400c9be54dd8592828f027423f5541c62d24dea511a
H_36(u): 482d125c032bafbbbfd75aa8338974bb8093382d366d9592dc0839b89d98dff6
H_24(u): 6d9a9379c7a3ab158360d3d09aae7acdd2ab8adfb498f2d68e87bf7f28b056f1
```

Moreover,

```text
gcd(H_24,H_36)=gcd(H_24,numer(Q_12(-1/6)))=gcd(J,H_30)=1.            (28)
```

Consequently the union of the two valid projection charts proves the row-`D`
critical length on the larger exact open set

```text
eta*J(eta)*H_30(eta)!=0.                               (29)
```

The factors `H_36` and `H_24` are projection-coordinate collisions, not
common obstructions: if one projection fails, the other remains squarefree.
This leaves `J=0` and `H_30=0` after the two-projection open calculation.  The
two exact closures are sections 5.1 and 5.2.

### 5.1 Exact closure of the endpoint wall `J=0`

Work over the exact quadratic field

```text
K_J=Q[eta]/(22143375 eta^2+15510536192)=Q(sqrt(-30)),
eta=(88064/18225)sqrt(-30).                            (29a)
```

Both complex points of `J=0` are represented by its two embeddings.  Rebuild
the row-`D` source `(7)`, put `f=G_X/T` and `h=G_T`, and project its residual
critical ideal by `(X,T)->T`.  Coefficientwise reduction of the exact identity

```text
Res_X(f,h)=T^35(6T+1)^2 Q_12(T)
```

modulo `J` kills the degree-twelve coefficient and no further coefficient.
The specialized polynomial `Q_11 in K_J[T]` has

```text
[T^11]Q_11=
 1887494744535166760726994878464 eta/12258226409765625,
Q_11(0)=-7518022905104433152 eta/2767921875,
Q_11(-1/6)=
 -2983872572026179618984706441216 eta/402131117372361328125. (29b)
```

Exact algebraic-field factorization gives one factor of degree eleven with
multiplicity one; equivalently `Q_11` is irreducible and squarefree over
`K_J`.  The frozen certificate prints all twelve coefficients in the basis
`1,eta`, so this conclusion is not based on numerical root separation.

The projection loses the `X`-coordinate.  Its required infinity sidecar is

```text
deg_X(f)=6,  [X^6]f=7T^7 eta;
deg_X(h)=7,  [X^7]h=8T^7 eta.                          (29c)
```

Since `eta!=0`, `Q_11(0)!=0`, and both leading coefficients in `(29c)` are
nonzero at every residual root, no finite-`T` intersection is lost at
`X=infinity`.  Squarefreeness therefore gives residual length eleven.  The
universal pairs are independently restored as

```text
T=-1/6, X^2=6:   two points, det Hess(G)=-6;
T=0,    X^2=-6:  two points, det Hess(G)= 6.
```

Thus the full affine critical length on `J=0` is exactly

```text
L=11+2+2=15.                                           (29d)
```

The boundary packet and response degrees do not change on this coefficient
slice.  For the full response, `n=18` and transported fixed sheets give
overlap at most `n-L=3`; hence the commutator-overlap lemma gives
`ind([X,Y])<=6<12`, while the origin meridian has packet defect `12`.

For the finite response, `(n,beta)=(12,3)`.  If both handles are nonidentity,
the sum of generator indices is at most

```text
2n-L-2+beta=10<n-1=11,
```

contradicting transitivity.  If either handle is the identity, then
`[X,Y]=1`; the surface relation would express the origin permutation, of
cycle type `(7,4,1)` and index `9`, as a product of the three carrier
transpositions, whose index is at most `3`.  This is also impossible.  Hence
the entire `J=0` slice is excluded.  Combining it with `(29)` and
`gcd(J,H_30)=1` proves row `D` is excluded whenever

```text
eta*H_30(eta)!=0.
```

### 5.2 Exact nonreducedness closure of `H_30=0`

Write the even irreducible factor as

```text
H_30(eta)=h_15(eta^2),              u=eta^2,            (29e)
```

and work first over `K_H=Q[u]/(h_15)`, then in the quadratic tower
`L_H=K_H[eta]/(eta^2-u)`.  Irreducibility of `H_30` in the exact factorization
implies both that `h_15` is irreducible and that the tower has degree two.
The typed source is the affine critical scheme
`Spec L_H[X,T]/(f,h)`, the target is the `T`-line, and the map is projection
`(X,T)->T`.  This projection forgets `X`; the linear subresultant below
recovers it, while `(29c)` is the infinity sidecar.

The exact parity reduction of the primary residual and the primitive linear
`X`-subresultant is

```text
Q_12(T,eta)=eta Q_u(T,u),
S_1=A(T,u)X+eta B(T,u),                               (29f)
deg_T(Q_u,A,B)=(12,12,11).
```

Here `S_1` is obtained from the penultimate symbolic subresultant after
removing its exact content

```text
T^35 eta(6T+1)/6.                                     (29g)
```

In `K_H[T]`, exact gcd computation gives

```text
deg gcd(Q_u,dQ_u/dT)=1.                               (29h)
```

Let its root be `t_0`.  Exact coprimality checks give

```text
gcd(h_15,Q_u(0)Q_u(-1/6))=1,
gcd(h_15,Res_T(Q_u,A))=1,
gcd(H_30,eta J)=1.                                    (29i)
```

Thus `t_0(6t_0+1)` and `A(t_0)` are units.  The specialized linear
subresultant is nonzero, so the fibre over `t_0` has exactly one finite common
zero

```text
x_0=-eta B(t_0)/A(t_0)
```

of `f=G_X/T` and `h=G_T`.  The leading rows `(29c)` are units at `t_0`, so no
second common zero is hidden at `X=infinity`.

The load-bearing local elimination fact is the following.  When the two
`X`-leading coefficients are units at `t_0`, the `t_0`-adic valuation of
`Res_X(f,h)` equals the sum of the local intersection multiplicities of the
finite common zeros over `t_0`.  Equation `(29h)` makes `t_0` a repeated root
of the residual resultant.  Since there is exactly one point above it, that
point has local intersection multiplicity at least two.  Hence

```text
det D(f,h)(x_0,t_0)=0.
```

Direct differentiation gives the polynomial identity

```text
T det D(f,h)=det Hess(G)+f G_XT.                       (29j)
```

At the critical point, `f=0` and `t_0` is a unit, so
`det Hess(G)(x_0,t_0)=0`.

This cannot occur for a Keller realization.  Indeed, for `G=E(A,C)` with
`det D(A,C)=1`, the identity `grad G=D(A,C)^t grad E` sends every affine
critical point to one of the two target nodes.  There the second-derivative
chain rule reduces to

```text
Hess(G)=D(A,C)^t Hess(E)D(A,C),
```

whose determinant is nonzero by THM-4130 and THM-4147.  The exact point above
therefore contradicts the Keller--Morse congruence under every complex
embedding of `L_H`.  Thus `H_30=0` is excluded.  Together with `(29)` and
section 5.1, this proves the complete row-`D` statement

```text
eta!=0  ==>  no row-D Keller realization.              (29k)
```

## 6. Prime cubic carrier and exact response degrees

The nonrational boundary edge is always

```text
q-1/2=K W^2-eta W^3.                                  (30)
```

Because `eta!=0`, the rational function on the right has degree three, so,
with `K_0=C(q)` and `L=C(W)`, one has `[L:K_0]=3`.  This is one separable
cubic closed point and remains true when `K=0`.  The other punctures are
rational.  THM-4120 gives
`E_q(C(q))={O}`, hence all rational punctures map to the target origin.

If the cubic packet is finite and `M` is the residue field of its target
image, functoriality gives

```text
K_0 subset M subset L.                                      (31)
```

Primality gives `M=K_0` or `M=L`.  The first alternative would be a forbidden
finite `K_0`-rational target point.  Hence `M=L`.  Equivalently, after using
THM-4122 to normalize the horizontal image, the boundary maps birationally to
one degree-three carrier and yields three distinct conjugate transpositions.
Its total carrier index is

```text
beta=3.                                                   (32)
```

If it maps to `O`, all three conjugates do so.  Thus the rational-index sums
and exhaustive responses are

| row | rational index sum = finite `n` | full `n` |
|---|---:|---:|
| `A` | `7+7+4+1=19` | `19+3*2=25` |
| `B` | `6+6+4+1=17` | `23` |
| `C` | `5+5+4+1=15` | `21` |
| `D` | `7+4+1=12` | `18` |

The carrier polynomial has exact discriminant

```text
Disc_W(-eta W^3+K W^2-(q-1/2))
 =(q-1/2)(4K^3-27eta^2(q-1/2)).                        (33)
```

Thus `q=1/2` is a carrier collision value (a triple collision if `K=0`), and
when `K!=0` the second value is `q=1/2+4K^3/(27eta^2)`.  Both belong to the
finite deleted Hurwitz set; the proof never treats the cubic as three global
sections across them.

The transport input is the finite-separable-carrier lemma proved in THM-4147,
applied to the new complete packets above.  Explicitly, first resolve the compactified pencil graph
and shrink the base until the two families are proper and smooth and the
resolved map is fibrewise quasifinite; proper plus quasifinite makes it finite.
Only then delete `O`, the complete cubic carrier when finite, its full inverse
image, and the finite extra Hurwitz discriminant.  The restriction is finite
etale and separable in characteristic zero.  Parallel Milnor cores avoid carrier/node
collisions.  Ehresmann transport then conjugates sheet labels but preserves
the distinct fixed sheets contributed by every affine node inverse.  This is
the source-wall-independent transport gate used in THM-4138/4141.  The proof
does not assume a coordinate formula for the cubic target carrier.

## 7. The two monodromy exclusions

Let `X,Y` be handle permutations around the two nodal target fibres.  Fixed
sheet transport gives

```text
#Fix(X)>=r_0,       #Fix(Y)>=r_1,
|supp X|+|supp Y|<=2n-L.                               (34)
```

### Finite response on the critical-open rows

The three carrier meridians have total index `beta=3`.  If at least one of
`X,Y` is nonidentity, the maximum total generator index is

```text
ind(X)+ind(Y)+beta <=2n-L-1+beta=n-2<n-1,              (35)
```

because `L=n+beta+1` in each of the four critical-open rows.  If both handles
are identities, total index is only three.  A transitive system on `n` letters
needs at least `n-1` orbit mergers, so every critical-open finite response is
impossible.  The `J=0` finite response uses the sharpened two-case argument in
section 5.1; `H_30=0` is already excluded locally in section 5.2.

### Full response on the critical-open rows

Now `X,Y` alone generate transitively.  Their supports cover all `n` letters,
and `(33)` gives

```text
k=|supp(X) intersect supp(Y)|<=n-L=2.                  (36)
```

For arbitrary permutations, THM-4147's commutator-overlap lemma gives

```text
ind([X,Y])<=2k.                                        (37)
```

Indeed, the partial injection `X|supp(Y)` differs from the identity on only
the `k` overlap letters.  Extend its directed path/cycle decomposition to a
permutation `Xhat` with `ind(Xhat)<=k` and the same conjugation of `Y`.
Cayley transposition distance then gives

```text
ind([X,Y])
 <=ind(Xhat)+ind(Y Xhat^-1 Y^-1)<=2k.                  (38)
```

The origin meridian is `[X,Y]^-1`, whose index is the complete packet defect.
Thus rows `A,B,C,D` would require respectively

```text
18,16,14,12 <= ind([X,Y]) <=4,                         (39)
```

all impossible.  The primary certificate independently checks `(37)` through `S_5`
in addition to the canonical general proof.  Again, the two row-`D` wall
closures are handled separately in sections 5.1 and 5.2.

## 8. Firewalls and lost data

The cubic object in `(30)` is a **residue-degree-three boundary place**.  It
is not an arithmetic period-three orbit, and no map from the rational
three-cycle construction of THM-4139/4146 to `(30)`, the critical scheme, or
the target carrier is asserted.  Arithmetic three-cycles remain fully
firewalled.

The monodromy compression intentionally forgets carrier coordinates and the
ordering of the three conjugates.  The sidecar that must survive is exactly

```text
(closed-point degree 3, beta=3, L, packet defect).      (40)
```

For the full response, the additional invariant is support overlap `k`, not
the raw source degree.

## 9. Honest residual walls

This theorem itself does not cross:

1. `eta=0`, which leaves exact weight nine/repeated-top scope;
2. `Delta=0`, which contracts the common Newton polygon;
3. in rows `A,B,C`, any zero of a leading resultant coefficient,
   `Q_d(-1/6)`, `disc_T Q_d`, `R_d(0)`, `[p^d]R_d`, or `disc_p R_d`;
4. failure of the inherited finite-separable transport setup, including a
   nondeleted carrier/node or extra Hurwitz collision;
5. exact weight at least ten, another reduced cell, or entry into this
   reduced `(2,3)` seam.

The structural equalities `C=0`, `Phi=0`, and `B=0` are **not** silently
discarded: they are rows `B,C,D` and have been treated separately.  Row `D`
has no remaining coefficient/resultant/discriminant wall in the present
scope: `H_36,H_24` are removed by projection switching, `J=0` by section 5.1,
and `H_30=0` by section 5.2.  Conversely, this complete deepest-row result
does not itself close any row-`A,B,C` wall in item 3. THM-4173 subsequently
closes row `A`, and THM-4176 closes rows `B,C`; those successors leave only
items 1, 2, 4, and 5 open on the present route.

## 10. Discharged continuation

THM-4171/4173 discharged the row-`A` continuation below. THM-4176 then exports
the same resultant-multiplicity/Hessian mechanism to every row-`B,C` degree
wall. The following preserves the then-decisive computation that led there.

After complete row `D`, the sharpest explicit remaining coefficient wall was
the row-`A` critical-infinity wall exposed by the canonical independent
anti-diagonal computation.  Put

```text
D_A=4Theta K^2-27eta^2.                                (41)
```

Its exact endpoints are

```text
Q_19(0)=-12288 C^6,
[T^19]Q_19=-1458 C eta^4 D_A^2,
R_19(0)=46656 eta D_A,
[p^19]R_19=1327104 eta^5 C^4.                         (42)
```

Thus `D_A=0` is common to both endpoint ledgers rather than a one-projection
coordinate collision.  It is nonempty already at the exact rational control

```text
Delta=2, eta=7, Phi=5, K=2743/45,
Theta=2679075/30096196,                                (43)
```

where `C!=0`.  The decisive computation is to specialize the affine critical
ideal, its projective `X`-infinity chart, and the independent chart at `p=0`
to `D_A=0`, saturating by

```text
eta*C*T*(6T+1),                                        (44)
```

The exact test is binary: either the lost resultant multiplicity is supported
at a unique finite point, so the local lemma of section 5.2 gives an immediate
Hessian contradiction, or it lies at critical infinity and its strict
transform supplies the next required sidecar.  This is the first row-`A`
wall not settled by the canonical open-chamber computation at that stage.

## 11. Exact artifacts and replay

The primary certificate checks all four rational controls, both critical
projections, every strict transform, the deepest-row discriminant portfolio,
and `15,017` commutator hostiles. Its symbolic mode proves the endpoint
formulas rather than sampling them. Three companion certificates reconstruct
the deepest factorization, the `J=0` field calculation, and the `H_30=0`
local nonreducedness certificate.

```text
04-computation/jc23_weight9_repeated_top_wall_thm4157.py
05-knowledge/results/jc23_weight9_repeated_top_wall_thm4157.out
05-knowledge/results/jc23_weight9_repeated_top_wall_symbolic_endpoints_thm4157.out

04-computation/jc23_weight9_repeated_top_wall_factors_thm4157.py
05-knowledge/results/jc23_weight9_repeated_top_wall_factors_thm4157.out

04-computation/jc23_weight9_repeated_top_wall_J_wall_thm4157.py
05-knowledge/results/jc23_weight9_repeated_top_wall_J_wall_thm4157.out

04-computation/jc23_weight9_repeated_top_wall_H30_thm4157.py
05-knowledge/results/jc23_weight9_repeated_top_wall_H30_thm4157.out
```

The four clean-room implementations independently cover the strict-transform
hierarchy, quotient algebra, primitive second projection, common deepest
factor, and `J=0` response:

```text
04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_newton_audit.py
04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_quotient_audit.py
04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_projection_audit.py
04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_J_wall_audit.py
```

A fifth clean-room audit proves `H_30` irreducible by Rabin's test modulo
`347`, proves the repeated `T`-root is exactly double, and separately verifies
at the resulting algebraic point that

```text
f=h=det D(f,h)=det Hess(G)=0.                           (45)
```

This independently closes the only fact imported by the primary `H_30`
certificate from its locked factorization and confirms the local-length
argument pointwise.

```text
04-computation/jc23_weight9_repeated_top_wall_H30_thm4157_independent_point_audit.py
05-knowledge/results/jc23_weight9_repeated_top_wall_H30_thm4157_independent_point_audit.out
04-computation/jc23_weight9_repeated_top_wall_H30_thm4157_independent_tower_audit.py
05-knowledge/results/jc23_weight9_repeated_top_wall_H30_thm4157_independent_tower_stage1_audit.out
```

Replay the canonical certificates with

```bash
python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157.py
python3 -B -O 04-computation/jc23_weight9_repeated_top_wall_thm4157.py
PYTHONHASHSEED=314159 python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157.py --symbolic-endpoints

python3 -B 04-computation/jc23_weight9_repeated_top_wall_factors_thm4157.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_J_wall_thm4157.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_H30_thm4157.py

python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_newton_audit.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_quotient_audit.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_projection_audit.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_thm4157_independent_J_wall_audit.py
python3 -B 04-computation/jc23_weight9_repeated_top_wall_H30_thm4157_independent_point_audit.py
H30_AUDIT_STAGE1_ONLY=1 python3 -B 04-computation/jc23_weight9_repeated_top_wall_H30_thm4157_independent_tower_audit.py
```

Normal, optimized, and fixed-hash-seed runs byte-match each frozen output.
The primary default and symbolic runs pass `15,194` and `15,200` exact
checks. The raw LF-byte and semantic hashes in the front matter bind the
canonical artifacts.

This theorem proves only the repeated-wall strata stated in Section 1. It
does not prove JC(2) or DC(2). **QED.**
