---
id: THM-4159
title: "Complete Y-only inner-resultant wall planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
  + VERIFIED-EXACT + INDEPENDENT SOURCE-PAIR AUDIT. The complete Y-only
  eta=Delta=0 exact-weight-nine inner wall I_C=0 with zeta*Disc(C)!=0
  contains no nonautomorphic planar Keller pair. Its three exhaustive
  critical strata have lengths 21,20,19; all retain genus 9 and packet
  (8,3,3,3,2,2,2,1). The carrier-orbit lemma excludes every finite response.
  THM-4161/4164 subsequently close the nontriple/triple top collisions off
  `I_C=0`; the common `I_C=0,Disc(C)=0` intersections, zeta=0, other cells,
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
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
script: 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159.py
output: 05-knowledge/results/jc23_y_only_inner_resultant_wall_exclusion_thm4159.out
independent_audit_script: 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_y_only_inner_resultant_wall_exclusion_thm4159_independent_audit.out
script_sha256: d97f30ad1321916f6ee04afd6db35d03409ed9a6d952700bb5da7bfb34cc8c8e
output_sha256: f28dac46e3199ae53fc0fca20794d64bd85b2af35d438d78c7660e749b996a0d
independent_audit_script_sha256: 4dacce587f017c3c88162f96e7a107f49b2360ad0979c8e483fa15baf136db6a
independent_audit_output_sha256: 590623feadebaaaf69abf4a2cc6cd0432236e3346a45f386051cc9c77579d034
semantic_sha256: 27bc93095e0f13030049e97584d769e7b78d47a0a81a0945458c250c5e617f9a
independent_semantic_sha256: f2feb98cb9cd02dd5dd2c4346b5c0a0dbc992e6db96fe7282bf3150d24032a69
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact symbolic elimination reconstructs the complete parameter map,
  the three critical-length strata, both endpoint ladders, all terminal gcd
  firewalls, the unchanged boundary packet, and the full and finite
  monodromy contradictions. Normal, optimized, and hash-seeded outputs
  byte-match.
independent_audit: >
  ACCEPT. A clean-room source-coordinate referee uses the alternative
  critical pair (A,C_0), resultant rather than gcd terminal firewalls,
  disjoint squarefree controls, and independent face/response reconstruction.
---

# THM-4159 -- complete Y-only inner-resultant wall exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
+ VERIFIED-EXACT + INDEPENDENT SOURCE-PAIR AUDIT; JC(2) REMAINS OPEN.**

Work over `C` in the live `b=d=0` reduced `(2,3)` seam. Put

```text
K0=2848/45,
C(W)=zeta W^3+Theta W^2+Phi W-1376/135,
D_C=Disc_W(C),
I_C=4Theta K0^2-27zeta^2.
```

## Theorem

The exact-weight-nine coefficient locus

```text
eta=Delta=0,             I_C=0,             zeta*D_C != 0
```

contains no nonautomorphic planar Keller pair.

More precisely, there is an exhaustive, disjoint critical-length
stratification. Define the unique parameter `u` and the two strict-transform
coordinates

```text
zeta=(5696/135)u,        Theta=3u^2,
u=3zeta/(2K0),

J=8544Phi-22784u-1215u^3,
S=2460375u^4-204543360u^2+5580439552.
```

Then

| coefficient stratum | affine critical length |
|---|---:|
| `J != 0` | `L=21` |
| `J=0, S!=0` | `L=20` |
| `J=S=0` | `L=19` |

These three rows cover every point of the stated locus. In every row the
normalization genus is `9`, the complete infinity packet is

```text
(8,3,3,3,2,2,2,1),
```

and both its full and finite cubic-carrier responses are impossible.

The inheritance pass is exact. The closest proved mechanism is THM-4155's
critical-length and labelled-carrier argument. Its canonical hostile is
`I_C=0`, where both critical endpoints vanish. The corrected near miss is
the old finite-response capacity: at `L=19` it is no longer strict. The
least-used sidecar is the number of orbits of the handle subgroup
`<X,Y>` before adjoining the three carrier transpositions. Retaining that
orbit count produces the strict replacement in Sections 4--5.

## 1. Parameter-map coverage

If `I_C=0` and `zeta!=0`, set `u=3zeta/(2K0)`. Then `u!=0`,

```text
zeta=(2K0/3)u=(5696/135)u,
Theta=27zeta^2/(4K0^2)=3u^2.
```

Conversely these formulas give `I_C=0`. Hence the map loses no point and has
no sign quotient. For every `(u,Phi)`, exactly one of `J!=0`,
`J=0,S!=0`, and `J=S=0` holds. On `J=0`, necessarily

```text
Phi=u(1215u^2+22784)/8544.
```

Thus the displayed three rows are exhaustive and disjoint.

## 2. Two critical eliminations and the three strict transforms

Use

```text
s=XT,       p=T+s^2,       t=p-s^2,
P=T+X^2T^2,                Y=XTP,

H=-3p+(8/3)p^2-(1376/135)p^3+K0s^2p^2
  +Phi sp^3+Theta s^2p^3+zeta s^3p^3,

A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2.
```

The two source pairs `(A,B)` and `(A,C_0)` are independent exact
eliminations. The former gives

```text
J != 0:       Res_s(A,B)=p^7 R_17(p),
J=0,S!=0:     Res_s(A,B)=p^8 R_16(p),
J=S=0:        Res_s(A,B)=p^9 R_15(p),
```

while the latter gives the same residual degrees with exceptional powers
`p^9,p^10,p^11`. Off `J=0`, the first residual has endpoints

```text
R_17(0)=(8305770496/1125)u^2 J,
[p^17]R_17=-(23983352712374779904/759375)u^5 D_C.
```

On `J=0`, the next strict transform has

```text
R_16(0)=(2916352/16875)u^3 S,
[p^16]R_16=-(23983352712374779904/759375)u^5 D_C.
```

Independently, with `f=G_X/T` and `h=G_T`, normalized elimination gives

```text
Res_X(f,h)=T^56(6T+1)^2 Q_d(T),
d=17,16,15 in the three rows.
```

For `J!=0`,

```text
Q_17(0)=-(1519777094677765052956672/56953125)u^7,
[T^17]Q_17
 =-(23983352712374779904/13839609375)u^3 J^2 D_C.
```

For `J=0,S!=0`,

```text
Q_16(0)=-(1519777094677765052956672/56953125)u^7,
[T^16]Q_16
 =-(2956854296576/3113912109375)u^3 S^2 D_C.
```

It remains to ensure that `J=S=0` has no hidden intersection where the
degree falls again. Work in the reduced algebra

```text
Q[u]/(S).
```

The polynomial `S` is squarefree and coprime to `u`. The terminal normalized
endpoints are nonzero rational multiples of

```text
u(2005507674605933782764615u^2
  +316908385228357703794041472),

u(401085u^2-16287712),
```

and the terminal source endpoints are nonzero rational multiples of

```text
u(369170566011315u^2-23248683486112768),

u(3904455285u^2-155035505152).
```

Exact Euclidean algorithms give gcd `1` between `S` and each of these four
factors. On `J=0`, write

```text
D_C=-u^2 B(u)/99781787520000,

B(u)=30267225703125u^8+2043284356800000u^6
    +264381824212992000u^4+6498574373014732800u^2
    +498260889496415371264.
```

One also has `gcd(S,B)=1`. Thus every root of `S`, not merely a sample root,
has `D_C!=0` and nonzero source and normalized endpoints. This proves the
last residual degree is exactly `15`.

The common leading `X`-rows remain nonzero because `u!=0`; hence no omitted
affine critical point lies at `X=infinity`. The two universal pairs remain

```text
T=0,     X^2=-6,       G=0,       Hess(G)=+6,
T=-1/6,  X^2=6,        G=1/2,     Hess(G)=-6.
```

They restore four points. Under a hypothetical Keller realization the
inherited Hessian congruence makes the entire affine critical scheme Morse.
Consequently the exact lengths are `17+4=21`, `16+4=20`, and `15+4=19`.
No residual-`T` discriminant hypothesis is needed: projected root collisions
do not reduce scheme length under the Keller-Morse condition.

## 3. Boundary packet and responses

The direct generic-fibre expansion has Newton polygon

```text
(0,1),(2,0),(5,3),(3,4),(0,4),
(2Area,B,I)=(27,11,9).
```

Because `u!=0`, both `zeta` and `Theta` remain live. Because `D_C!=0`, the
constant top cubic has three distinct roots. The index-two face is

```text
q-1/2=K0W^2+zeta W^3.
```

It defines a separable prime cubic residue extension over `C(q)`. The
complete packet and labels are therefore unchanged from THM-4155:

```text
(8,3,3,3,2,2,2,1)
= rational (8,3,3,3,1) + cubic orbit (2,2,2).
```

Its defect is `16=2*9-2`; the inherited connectedness and
Riemann--Hurwitz/Pick argument gives exact genus `9` and packet completeness.

The full response has degree `24`. In the finite response, prime residue
transport forces the full cubic carrier rather than a quotient carrier.
After the finite splitting base change it gives three distinct target
punctures on the same `18`-sheet permutation action, each with a
transposition meridian. Hence

```text
full:    n=24, origin index=16;
finite:  n=18, origin packet=(8,3,3,3,1), origin index=13,
         three carrier transpositions, total carrier index beta=3.
```

The finite origin index is exactly

```text
(8-1)+3(3-1)+(1-1)=13.
```

No residue quotient is taken in the `+3`: each of the three separated
carrier meridians acts on this same `18`-sheet cover and has Cayley index
one. If products cancel, their product index only decreases, which
strengthens the upper bound below.

## 4. Carrier-orbit lemma

> **Lemma.** Let `n>m+1`, and let `X,Y,tau_1,...,tau_m` generate a transitive
> permutation action on `n` letters, with every `tau_j` a transposition. Put
> `a=|supp X|`, `b=|supp Y|`, and
> `k=|supp X intersect supp Y|`. Suppose transported fixed sheets give
> `a+b<=2n-L`. Then
>
> ```text
> |supp X union supp Y| >= n-m,
> k <= n+m-L.
> ```

**Proof.** Put `H=<X,Y>`. Adding one transposition can merge at most two
current `H`-orbits, so it reduces the orbit count by at most one. Since the
action after adjoining all `m` transpositions is transitive, `H` has at most
`m+1` orbits. Thus

```text
n-#Orb(H) >= n-m-1.                         (a)
```

Let `U=supp X union supp Y`. The inequality `n>m+1` and transitivity after
only `m` transpositions force `U` to be nonempty. Every point outside `U` is
an `H`-fixed singleton, while `U` contributes at least one further orbit.
Therefore

```text
n-#Orb(H) <= |U|-1.                         (b)
```

Equations `(a),(b)` give `|U|>=n-m`. Finally

```text
k=a+b-|U| <= (2n-L)-(n-m)=n+m-L.
```

Every inequality is in the needed direction. QED.

For the punctured elliptic target, the exact fundamental-group relation is

```text
[X,Y] mu_O tau_1 ... tau_m=1.
```

The inherited commutator lemma and the triangle inequality for permutation
index give

```text
ind([X,Y]) <= 2k,
ind(mu_O) <= 2k+m <= 2(n+m-L)+m.           (c)
```

## 5. Monodromy contradictions

In the full response `X,Y` alone generate transitively, so their support
union is all `24` sheets. The fixed-sheet inequality gives

```text
k<=24-L,
ind(mu_O)=ind([X,Y])<=2(24-L).
```

For `L=21,20,19`, the upper bounds are respectively `6,8,10`, all strictly
below the exact origin index `16`.

In the finite response use the lemma with `(n,m)=(18,3)`. It gives

```text
|supp X union supp Y|>=15,
k<=21-L.
```

Equation `(c)` gives the respective origin-index ceilings

```text
L=21: 3,          L=20: 5,          L=19: 7.
```

All are strictly below the exact finite origin index `13`. Thus neither
response exists in any of the three exhaustive coefficient strata. This
proves the theorem.

## 6. Exact scope

This closes the complete `I_C=0` stopping wall of THM-4155 under
`zeta*D_C!=0`. It does not cross the top-face collision `D_C=0`, the
`zeta=0` coefficient wall, another reduced cell, seam entry, exact residual
weight at least ten, `JC(2)`, or `DC(2)`.

## 7. Exact artifacts and replay

```bash
python3 -B 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159.py
python3 -B -O 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159.py
PYTHONHASHSEED=181 python3 -B 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159.py

python3 -B 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159_independent_audit.py
python3 -B -O 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159_independent_audit.py
PYTHONHASHSEED=223 python3 -B 04-computation/jc23_y_only_inner_resultant_wall_exclusion_thm4159_independent_audit.py
```

Both normal/optimized/hash-seeded output comparisons pass exactly.

**QED.**
