---
id: THM-4155
title: "Generic Y-only Delta-zero weight-nine planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147
  + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the Y-only exact-M=9
  wall eta=Delta=0, the open
  locus zeta*(4Theta*K0^2-27zeta^2)*Disc(C)!=0, where K0=2848/45 and
  C(w)=zeta*w^3+Theta*w^2+Phi*w-1376/135, contains no nonautomorphic
  planar Keller pair. Its normalized genus is 9, affine critical length is
  22, and labelled infinity packet is (8,3,3,3,2,2,2,1).
  THM-4159/4161/4164/4165 subsequently close every zeta!=0 inner/top
  boundary. The zeta=0 wall, other coefficient cells, entry, M>=10, JC(2),
  and DC(2) remain OPEN.
source: codex-lrc14-planar-jc-breakthrough-20260825
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
script: 04-computation/jc23_generic_y_only_delta_zero_weight_nine_exclusion_thm4155.py
output: 05-knowledge/results/jc23_generic_y_only_delta_zero_weight_nine_exclusion_thm4155.out
script_sha256: 4e9e8ca0c1a3f864003ee5869882f15439669858cbaf22096b02724a376d8e32
output_sha256: 8d78bc243e4027c76f6116d35c5dcd16191e45221ef4a80bfa98240bae976513
semantic_sha256: e24081ff7259ce7c61abba1c5476bf2cf68588b934ef0d19ee461100d3e97d17
independent_audit_script: 04-computation/jc23_weight9_y_only_delta_zero_wall_thm4155_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_weight9_y_only_delta_zero_wall_thm4155_independent_audit.out
independent_audit_script_sha256: 6814f463adcd61d75c626135eab140e53965e9e41312ad8ea689d5a7836e2eb5
independent_audit_output_sha256: 65a06be05dc43c74fdded1690327d58863760e0d924427077df84007979b812b
independent_semantic_sha256: e1a7880dff5092afe8b2de6ca30d412d00656cec6e73744d09655c7e9441358c
hash_basis: raw LF bytes
primary_audit: >
  PASS. The standalone certificate computes the general symbolic critical
  resultants independently in source and normalized coordinates, checks both
  endpoint factorizations, restores the four universal critical points,
  reconstructs the valued Newton polygon and every edge label, verifies a
  nonempty squarefree exact control, and checks both strict monodromy
  contradictions. Normal, optimized, and hash-seeded executions byte-match
  the frozen output.
independent_audit: >
  ACCEPT. A clean-room SymPy referee imports no primary code. It uses the
  alternative source pair (A,C_0), rather than the primary pair (A,B), and
  obtains p^8*R18. It expands the generic fibre directly to recover every
  face polynomial, separates the rational top cubic from the prime cubic
  carrier, checks an independent squarefree control, and exhausts the
  commutator-support containment through S_6. Its 1,600,317 live checks pass
  in normal, optimized, and hash-seeded runs.
---

# THM-4155 -- generic Y-only Delta-zero weight-nine exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147
+ VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; JC(2) REMAINS OPEN.** Work
over `C` in the live `b=d=0` reduced `(2,3)` seam.

## 1. Statement and inheritance

Put

```text
K0=2848/45,
C(W)=zeta W^3+Theta W^2+Phi W-1376/135,
D_C=Disc_W(C),
I_C=4Theta K0^2-27zeta^2.                              (1)
```

> **Theorem.** The exact-weight-nine coefficient locus
>
> ```text
> eta=0,                 Delta=0,
> zeta*I_C*D_C != 0                                      (2)
> ```
>
> contains no nonautomorphic planar Keller pair.

This is a genuine codimension-one continuation of chamber `Y` in THM-4147,
not an assertion about all of its boundary. The inheritance pass is:

- THM-3992/3997 and THM-4007 give the complete normalized source and force
  `K=2848/45-(7/6)Delta`, hence `K=K0` on this wall;
- THM-4147 supplies the closest proved mechanism: exact critical length,
  complete labelled boundary response, finite-separable-carrier transport,
  and the two permutation inequalities;
- the hostile examples are the two factors `I_C=0` and `D_C=0`; the first
  destroys a critical endpoint and the second collides the three new
  index-three places;
- the least-used sidecar is the total finite-carrier permutation index
  `beta`, retained below rather than suppressing the carrier after base
  change.

The coefficient `Phi` is unrestricted except through the displayed cubic
discriminant. No discriminant of the residual polynomial in `T` is assumed.
The exact control

```text
Theta=2,               Phi=3,               zeta=5       (3)
```

has

```text
I_C=63521957/2025,      D_C=-10233928/135,                (4)
```

and its residual critical polynomial is squarefree, so the open locus `(2)`
is nonempty.

## 2. Complete wall source and two critical projections

Use the source and normalized coordinates

```text
s=XT,          p=T+s^2,          t=p-s^2=T,
P=T+X^2T^2,                     Y=XTP.                    (5)
```

On `(2)`, the complete exact-`M=9` source from THM-4147 specializes to

```text
G=-s^2/(2t)+H(p,s),

H=-3p+(8/3)p^2-(1376/135)p^3+K0 s^2p^2
  +Phi sp^3+Theta s^2p^3+zeta s^3p^3.                   (6)
```

There is no omitted weight-nine term: `eta=0` deletes `p^3y` and
`Delta=0` deletes `p^4`, while `zeta y^3=zeta s^3p^3` keeps the maximum
residual weight equal to nine.

Define the polynomial source-critical pair

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2.                                         (7)
```

Direct differentiation gives

```text
t^2 G_s=pA,        2t^2 G_p=t^2B-sA.                    (8)
```

The leading `s`-rows are `3zeta p^2` and `9zeta p^2`; hence `(2)` excludes
a source-critical point at `s=infinity` for `p!=0`. Exact symbolic
elimination gives

```text
Res_s(A,B)=p^6 R_18(p),

R_18(0)=-46656 zeta I_C,
[p^18]R_18=-236196 zeta^5 D_C.                          (9)
```

Thus every endpoint in this projection is live under `(2)`.

Independently, in normalized coordinates put

```text
f=G_X/T,                  h=G_T.                         (10)
```

These have `X`-degrees `(8,9)` and the same leading row
`9zeta T^8`. A second general-parameter symbolic elimination gives

```text
Res_X(f,h)=T^56(6T+1)^2 Q_18(T),

Q_18(0)=-(3^15/2^7)zeta^7,
[T^18]Q_18=-(6561/4)zeta^3 I_C^2 D_C.                  (11)
```

The factor `T^56` is a Sylvester degree-drop artifact. The actual critical
ideal restores the two points

```text
T=0,                    X^2=-6,             G=0,         (12)
```

whose Hessian determinant is `+6`. The factor `(6T+1)^2` gives the other
two universal points

```text
T=-1/6,                 X^2=6,              G=1/2,       (13)
```

whose Hessian determinant is `-6`. Consequently the affine critical-scheme
length is exactly

```text
L=18+2+2=22.                                             (14)
```

Under a hypothetical Keller realization `G=E(A_0,C_0)`, the inherited
Hessian congruence makes every affine critical point Morse. If `r_i` is the
number above the target node of value `i/2`, then

```text
r_0+r_1=22,                    r_0,r_1>=2.               (15)
```

A repeated root of `Q_18` can only merge projected `T`-values; it does not
reduce the critical-scheme length. This is why `(2)` contains no hidden
residual-discriminant condition.

## 3. Newton polygon, genus, and labelled packet

For `F_Q=(s^2-p)(1-QH)-Qs^2/2`, a monomial `p^i y^j` of `H` contributes
the valued endpoints

```text
(j+2,i+j,1),                   (j,i+j+1,1).              (16)
```

After coincident coefficients are combined, the lower Newton polygon is

```text
(0,1), (2,0), (5,3), (3,4), (0,4),
(2Area,B,I)=(27,11,9).                                  (17)
```

For an edge with primitive inward normal `(u,v)` and level `c`, THM-4103
gives ramification index `e=u+v-c`. The complete edge ledger is

| edge | length | `(u,v;c)` | `e` | label |
|---|---:|---:|---:|---|
| `(0,1)--(2,0)` | 1 | `(1,2;2)` | 1 | rational |
| `(2,0)--(5,3)` | 3 | `(-1,1;-2)` | 2 | one cubic closed point |
| `(5,3)--(3,4)` | 1 | `(-1,-2;-11)` | 8 | rational |
| `(3,4)--(0,4)` | 3 | `(0,-1;-4)` | 3 | three rational points |

The vertical edge `s=0` has length three and consists of affine source
points, not punctures. It is therefore omitted from the infinity packet.

The length-three top face is exactly `C(W)`. Since the constant and leading
coefficients are nonzero and `D_C!=0`, it has three distinct nonzero roots
in the algebraically closed constant field `C`; these give three rational
index-three places. The index-two face is

```text
q-1/2=K0 W^2+zeta W^3.                                  (18)
```

The rational function on the right has degree three. Therefore
`C(W)/C(q)` has prime degree three and is separable. The complete labelled
packet is

```text
(8,3,3,3,2,2,2,1),

rational (8,3,3,3,1) + one cubic orbit (2,2,2).          (19)
```

Its boundary defect is

```text
7+3(3-1)+3(2-1)=16=2*9-2.                               (20)
```

Pick gives normalization genus at most nine. As in THM-4147, the finite
critical scheme and THM-3827's closed-polynomial factor theorem make the
geometric generic source connected: a nontrivial polynomial factorization
would put a complete curve in that finite scheme. Riemann--Hurwitz over the
elliptic target and the displayed boundary defect give genus at least nine.
Hence the normalized genus is exactly nine, `(19)` is complete, and there
is no hidden affine ramification or further genus loss.

## 4. The cubic carrier and its two responses

THM-4120 gives `E_q(C(q))={O}`. Thus every rational place in `(19)` maps to
the target origin. For the cubic place, let `L_3/C(q)` be its residue
extension and `M` the residue field of a finite horizontal target image.
Then

```text
C(q) subset M subset L_3.                               (21)
```

Prime degree gives `M=C(q)` or `M=L_3`. The first alternative would produce
a forbidden finite `C(q)`-point of the target. Thus a finite image has the
full cubic residue field and splits geometrically into three conjugate
carrier points, each with a transposition meridian. There are exactly two
responses:

```text
full:       n=8+3+3+3+2+2+2+1=24;
finite:     n=8+3+3+3+1=18,        beta=3.              (22)
```

Here `beta` is the total permutation index of the three carrier meridians.

For completeness, the carrier transport used here is the
finite-separable-carrier lemma proved in THM-4147, not a cosmetic deletion of
the cubic orbit. Make the finite base change splitting `(18)`, normalize and
resolve the two pencil families, and delete the finite set of singular,
degree-drop, section-collision, origin-collision, and Hurwitz-collision
parameters. The three carrier sections are then distinct. They cannot agree
generically: agreement would descend their image through a proper
intermediate residue field, contrary to `(21)` and THM-4120.

After the proper smooth construction, delete the origin and carrier sections
and their complete inverse image. The remaining morphism is finite etale. A
small parallel Milnor core at either target node avoids the finitely many
marked sections, even if a carrier specializes to the node, and its
degree-one lifts transport to a common smooth reference fibre. Hence, for
the two handle permutations `X_0,X_1`, one has

```text
#Fix(X_i)>=r_i.                                          (23)
```

The two handles together with the three carrier transpositions generate a
transitive action because the connected source minus finitely many marked
points remains connected. In the full response the same argument has no
carrier sections and the two handles generate transitively. This preserves
the needed sidecar instead of treating finite and full responses as the same
cover.

## 5. Two strict monodromy contradictions

Equation `(15)` and `(23)` imply, in either response,

```text
|supp X_0|+|supp X_1| <= 2n-L.                          (24)
```

For the full response, transitivity makes the two supports cover all `24`
sheets. Their intersection therefore has size at most

```text
n-L=24-22=2.                                             (25)
```

THM-4147's commutator-index lemma gives

```text
ind([X_0,X_1]) <= 2|supp(X_0) intersect supp(X_1)| <= 4. (26)
```

The origin meridian is the inverse commutator, while packet `(19)` has index
equal to its defect `16`. Thus `4<16`, excluding the full response. The
standalone certificates additionally check the underlying disagreement
containment and its support consequence `6<23`.

For the finite response, transitivity requires total generator index at
least `n-1=17`. A nonidentity permutation has index at most support minus
one. Adding the three carrier transpositions gives the exact capacities

```text
both handles nonidentity:      (2n-L-2)+beta=15;
exactly one identity:          (2n-L-1)+beta=16;
both identities:                              beta=3.    (27)
```

Every capacity is strictly below `17`, so the finite response is impossible.
This proves the theorem.

## 6. Exact stopping walls and replay

The proof deliberately stops at exactly the factors detected by both the
critical and boundary computations:

```text
zeta=0:       no longer the Y-only exact-weight-nine locus;
I_C=0:        R_18(0)=0 and [T^18]Q_18=0;
D_C=0:        [p^18]R_18=[T^18]Q_18=0 and top-face roots collide. (28)
```

The last two walls need new critical and boundary normalizations; neither is
excluded here. No claim is made about another reduced cell, entry into this
seam, exact weight at least ten, `JC(2)`, or `DC(2)`.

The primary script performs `58` live checks and the clean-room audit performs
`1,600,317`. Normal, optimized, and fixed-hash-seed runs byte-match their
frozen outputs. Replay with

```bash
python3 -B 04-computation/jc23_generic_y_only_delta_zero_weight_nine_exclusion_thm4155.py
python3 -B -O 04-computation/jc23_generic_y_only_delta_zero_weight_nine_exclusion_thm4155.py
PYTHONHASHSEED=73 python3 -B 04-computation/jc23_generic_y_only_delta_zero_weight_nine_exclusion_thm4155.py
python3 -B 04-computation/jc23_weight9_y_only_delta_zero_wall_thm4155_independent_audit.py
python3 -B -O 04-computation/jc23_weight9_y_only_delta_zero_wall_thm4155_independent_audit.py
PYTHONHASHSEED=127 python3 -B 04-computation/jc23_weight9_y_only_delta_zero_wall_thm4155_independent_audit.py
```

**QED.**
