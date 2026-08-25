---
id: THM-4130
title: "Theta-only extremal-seam critical-count monodromy obstruction"
status: >
  PROVED RELATIVE TO THM-3996/4053/4103/4120/4126 + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT on the smooth nonresonant theta-only exact
  maximum-weight-eight seam. That seam contains no nonautomorphic planar
  Keller pair. The result does not treat any collision wall, maximum residual
  weight at least nine, another reduced cell, or JC(2) in general.
source: codex-frontier-synthesis-creative-20260825k
depends_on:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4126-jc23-extremal-target-vertical-nonproperness
related:
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4124-planar-keller-integral-degree-ratio-all-vertex-shear
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
script: 04-computation/jc23_theta_only_extremal_seam_critical_monodromy_thm4130.py
output: 05-knowledge/results/jc23_theta_only_extremal_seam_critical_monodromy_thm4130.out
independent_audit_script: 04-computation/jc23_theta_only_extremal_seam_critical_monodromy_thm4130_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_theta_only_extremal_seam_critical_monodromy_thm4130_independent_audit.out
script_sha256: 9cc246bb07894f548ef5ff4814ca92e9384a77de54617e3311df341458420fb3
output_sha256: 020f79bcc07384e02623053adcaee20edb240ee5bef363f1f3bbe86056cafefd
semantic_sha256: 62f23741001fa94b0dbaddd0381d80e8a5c6f11f7be5df59fedc14ebb218bb58
independent_audit_script_sha256: 322756789ad0c47d421130b60db7dabe78066a8dd972ee707dcd7a6d57cbf8fa
independent_audit_output_sha256: 3231ec6d8429a9ffc25cac0fb94cf6679e1348434cddcf1aaf019d9d93affac5
independent_semantic_sha256: 62f23741001fa94b0dbaddd0381d80e8a5c6f11f7be5df59fedc14ebb218bb58
hash_basis: raw LF bytes
primary_audit: >
  PASS. The normalized (X,T) calculation factors the critical resultant,
  freezes every coefficient of its degree-sixteen residual, separates its
  T^42 degree-drop artifact from the two actual T=0 points, verifies the two
  universal Morse pairs and the Phi=0 value obstruction, and exhausts all
  seventeen shared-pivot support ledgers. Normal, optimized, and two
  hash-seeded executions byte-match the frozen output.
independent_audit: >
  ACCEPT. A clean-room rational (s,p) calculation imports no primary code,
  derives different critical equations, eliminates s, identifies its p^2
  degree-drop artifact, counts the p=0 and omitted t=0 strata separately,
  reconstructs the two Picard--Lefschetz generators, and rechecks the
  permutation obstruction with a dictionary implementation. Normal,
  optimized, and two hash-seeded executions byte-match; semantic hashes agree.
---

# THM-4130 -- critical monodromy empties the theta-only extremal seam

**PROVED RELATIVE TO THM-3996/4053/4103/4120/4126 + VERIFIED-EXACT +
INDEPENDENTLY VERIFIED-EXACT on one reduced seam; JC(2) OPEN.** Work over
`C`. Retain the exact reduced `(2,3)` cell, target pencil, and every
hypothesis of THM-4120 and THM-4126.

## 1. Theorem and inheritance

> **Theorem.** There is no nonautomorphic planar Keller pair in THM-4053's
> smooth nonresonant theta-only exact maximum-weight-eight survivor
>
> ```text
> delta=0,              theta!=0,              Delta_V!=0.        (1)
> ```

The closest proved mechanism is the complete generic-fibre cover packet:
THM-4103 gives the connected genus-eight source and ramification indices,
and THM-4120 forces degree `21` and puts the whole branch divisor over the
target origin. THM-4126 supplies the vertical nodal nonproperness
classification needed to transport sheets without loss away from the two
singular pencil values.

The canonical corrected near miss is THM-4124: `Q^2-cP^3` is not a target
automorphism, so the following coefficient reduction cannot be manufactured
by a nonlinear target shear. The least-used relevant sidecar is the affine
critical scheme of THM-4053's exact source polynomial. Its length, combined
with the two target vanishing cycles, is the missing invariant.

The proof first forces the exact ledgers

```text
#Crit_aff(E(A,C))=20,
r_0+r_1=20,                 2<=r_0,r_1<=18,             (2)
S_F=N_0 union N_1,
m_0+m_1=22,                 3<=m_0,m_1<=19,             (3)
```

where `r_i` is the number of affine preimages of the node of `N_i` and
`m_i=21-r_i`. The two vanishing-cycle permutations would then have supports
meeting in one point and would have a three-cycle as commutator. This
contradicts the inherited boundary cycle type

```text
(7,3,3,3,2,2,1).                                      (4)
```

## 2. Exact normalized polynomial

THM-3992/4053 give, on `(1)`,

```text
gamma=-a^3/2,                 lambda=-3/a^2,
alpha=8/(3a^7),               epsilon=-1376/(135a^12),
kappa=2848/(45a^12),          a!=0.                    (5)
```

Choose a square root of `a^5` and put

```text
X=a^(5/2)x,       T=t/a^5,       Phi=a^(29/2)phi,
Theta=a^17 theta, P=T+X^2T^2,    Y=XTP.                 (6)
```

After division by `a^3`, the entire source polynomial is

```text
calG=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
     +(2848/45)Y^2+Phi P^2Y+Theta P Y^2.               (7)
```

There is no omitted lower term. The inherited gates become

```text
Theta!=0,
W:=135Phi^2+5504Theta=135a^29 Delta_V!=0.              (8)
```

Changing the chosen square root changes the signs of `X,Phi` together and
does not change any assertion below.

## 3. The degree-sixteen eliminant and twenty affine critical points

Since `calG_X` is divisible by `T`, define

```text
f=calG_X/T,                       h=calG_T.              (9)
```

As polynomials in `X`, their generic degrees are `7,8`, and both leading
coefficients are

```text
8Theta T^7.                                             (10)
```

Direct exact elimination gives

```text
Res_X(f,h)=-T^42 Theta^3(6T+1)^2 Q_16(T),              (11)
deg_T Q_16=16,
[T^16]Q_16=
  (16842242073296896/20503125)Theta^2 W,
Q_16(0)=12288Theta^3.                                  (12)
```

All seventeen coefficient polynomials of `Q_16` are frozen in the primary
audit output. Equations `(8),(12)` make both endpoint coefficients nonzero.
The `T^42` in `(11)` is not a critical multiplicity: specialization gives

```text
f(X,0)=-X,                  h(X,0)=-(X^2+6)/2,          (13)
```

which have no common root. It is exactly the Sylvester degree-drop artifact
from the fall of the generic `X`-degrees. For `T!=0`, `(10)` and `Theta!=0`
exclude a common root at `X=infinity`.

The actual critical ideal is `(Tf,h)`. At `T=0`, equation `(13)` gives the
two points

```text
T=0,                    X^2=-6,                         (14)
```

both with value `calG=0`. At `T=-1/6`, exact reduction modulo `X^2-6`
gives

```text
f=h=0,                  calG=1/2,                       (15)
```

so this is a second universal pair.

It remains to turn resultant multiplicity into point count. Write

```text
E(U,V)=V^2-U^3+(3a^2/4)U+a^3/4.                        (16)
```

Because `calG` is a nonzero scalar-coordinate version of `E(A,C)` and
`det D(A,C)=1`,

```text
grad(E o F)=DF^t grad E.                               (17)
```

Thus every affine critical point lies over one of the two target nodes

```text
o_0=(-a/2,0),                   o_1=(a/2,0).            (18)
```

At such a point the gradient term in the second-derivative chain rule
vanishes, hence

```text
Hess(E o F)=DF^t Hess(E)DF,
det Hess(E o F)=det Hess(E)=+6a at o_0, -6a at o_1.    (19)
```

Both are nonzero. Therefore the affine critical scheme is reduced. At a
solution with `T!=0`, the Jacobian determinant of `(f,h)` is
`T^-1 det Hess(calG)`, so its local intersection multiplicity is one.
Equations `(11),(12)` consequently count

```text
2+16=18 distinct critical points with T!=0.             (20)
```

The two points in `(14)` are also simple by `(19)`. This proves

```text
#Crit_aff(calG)=18+2=20.                                (21)
```

The universal pairs `(14),(15)` lie over `o_0,o_1`, respectively. If

```text
r_i=#F^-1(o_i),                                         (22)
```

then `(21)` and `(17)` give `(2)`.

### 3.1 A sharp coefficient sub-obstruction

If `Phi=0`, the line `X=0` supplies critical `T`-coordinates

```text
1376T^2-240T+135=0.                                    (23)
```

Eliminating `T` against `z=calG(0,T)` gives the exact critical-value
polynomial

```text
29584z^2+14680z+10935.                                 (24)
```

Its values at the only permitted critical values are

```text
z=0: 10935,                  z=1/2: 25671.             (25)
```

Both are nonzero, contradicting `(17),(18)`. Hence already `Phi!=0` on any
hypothetical realization. This sub-obstruction is not needed for the final
monodromy contradiction.

The noncollision gate in `(8)` is load-bearing. For the exact hostile

```text
(Phi,Theta)=(5504,-743040)                              (26)
```

one has `Theta!=0,W=0`, and the specialized `Q_16` has degree `15`, not
`16`. No statement here crosses that collision wall.

## 4. Fibre defect and finite transport

THM-4120/4126 identify the global mapping degree with `21`. Since each
`r_i>=2` and `r_0+r_1=20`, each `r_i<=18<21`. THM-3996's exact fibre-defect
criterion therefore puts both node values in `S_F`. The three alternatives
in THM-4126 then force

```text
S_F=N_0 union N_1.                                     (27)
```

This proves `(3)`.

Put

```text
B=C-{0,a^3/2}.                                         (28)
```

THM-4126 also shows that the restriction of `F` over the target open
`E^-1(B)` has no nonproper value. It is quasifinite and etale because `F` is
Keller, and it is proper there by the definition of `S_F`. Hence it is a
finite etale cover of degree `21`. This is the load-bearing transport fact;
the generic degree alone would not prevent a sheet from escaping to infinity
while a vanishing cycle is moved through the pencil.

For a smooth `q in B`, THM-4103/4120 give the connected cover

```text
varphi_q:C_q -> E_q,          deg varphi_q=21,          (29)
```

branched only over the target origin `O`, with

```text
varphi_q^*(O) having indices (7,3,3,3,2,2,1).          (30)
```

After deleting the finite inverse image of `O`, the source remains connected,
so the monodromy action on the `21` sheets is transitive.

## 5. Affine node preimages inject into fixed sheets

> **Fixed-sheet lemma.** Let `delta_i` be the vanishing cycle of the target
> node `o_i`, transported to one smooth reference fibre. If `Xi_i` is the
> monodromy permutation of `(29)` along `delta_i`, then
>
> ```text
> #Fix(Xi_i)>=r_i.                                      (31)
> ```
>
> Distinct affine preimages give distinct fixed sheets. The assertion is an
> injection; it makes no claim that every fixed sheet comes from an affine
> node preimage.

### Proof

Let `z in F^-1(o_i)`. The target holomorphic Morse lemma gives coordinates
`u,v` near `o_i` with

```text
E-q_i=uv.                                               (32)
```

The Keller identity makes `F` a local biholomorphism at `z`. Pulling back
`u,v` therefore gives source coordinates with the **same** base equation
`uv=q-q_i`; there is no ramified replacement `q-q_i=tau^ell`. For nearby
nonzero `q-q_i`, the local Milnor annulus maps biholomorphically to the target
annulus. Its vanishing circle consequently has a closed degree-one lift, so
the corresponding sheet is fixed.

Choose disjoint inverse neighborhoods around the distinct points of
`F^-1(o_i)`. Their local annuli give distinct sheets. Finite-etale transport
over `(28)` is a bijection on sheets, so distinct sheets cannot merge while
they are carried to the reference fibre. Moreover, the vanishing circles
avoid `O`, and THM-4120's full boundary exhaustion says every source-boundary
point lies over `O`; no sheet over these circles is an uncounted point at
source infinity. Possible additional fixed sheets only strengthen `(31)`.
QED.

Local etaleness and preservation of the base parameter are both essential.
The sharp hostile

```text
(u,v) |-> (u^ell,v^ell),       uv=tau,
xy=tau^ell                                               (33)
```

has an `ell`-cycle after one target-base turn. Its Jacobian
`ell^2(uv)^(ell-1)` is not a unit, and it changes the base parameter. Thus it
does not satisfy the lemma's hypotheses.

## 6. The two vanishing cycles generate the punctured torus

After a complex scaling set `a=1` and choose the real reference value
`q_*=1/4`. The target cubic is

```text
V^2=U^3-(3/4)U+q-1/4.                                  (34)
```

At the two nodal values its root packets are

```text
q=0:   (-1/2,-1/2,1),
q=1/2: (-1,1/2,1/2).                                   (35)
```

At `q_*`, lift the two adjacent real root intervals. Their lifts
`delta_0,delta_1` are simple vanishing cycles, avoid `O`, and intersect once.
Their regular neighborhood is a once-punctured torus whose complementary
disk contains `O`. Hence

```text
pi_1(E_(q_*)-{O})=<delta_0,delta_1>                    (36)
```

is free of rank two. With a choice of orientation, the boundary loop about
`O` is

```text
[delta_0,delta_1] or [delta_0,delta_1]^-1,              (37)
```

up to conjugacy. A matrix cross-check uses the two positive twists

```text
[[1,1],[0,1]],             [[1,0],[-1,1]].             (38)
```

Their product has trace `1`, the `II*` trace, while the general trace is
`2-(delta_0 dot delta_1)^2`; this independently confirms intersection
absolute value one.

Let

```text
X=rho(delta_0),                 Y=rho(delta_1).         (39)
```

By `(30),(37)`, either `[X,Y]` or its inverse has cycle type `(4)`. Inverse
and conjugate conventions do not change cycle type. This audits the
commutator convention without choosing an orientation gauge.

## 7. Support-intersection obstruction

> **Permutation lemma.** Let `X,Y in S_21` generate a transitive action.
> Suppose
>
> ```text
> #Fix(X)>=r_0,       #Fix(Y)>=r_1,       r_0+r_1=20.   (40)
> ```
>
> If their commutator is nontrivial with cycle type `(4)`, then `(40)` is
> impossible.

### Proof

Put `A=supp(X),D=supp(Y)`. Equation `(40)` gives

```text
|A|+|D|=42-#Fix(X)-#Fix(Y)<=22.                        (41)
```

A common fixed point is precisely a point outside `A union D`; it would be a
singleton orbit of both generators. Transitivity therefore gives

```text
A union D={1,...,21}.                                  (42)
```

Thus the set of shared **nonfixed** points satisfies

```text
|A intersect D|=|A|+|D|-21<=1.                         (43)
```

This treats shared fixed and shared nonfixed points separately. If one
support were empty, `[X,Y]` would be the identity, contrary to `(4)`. If both
were nonempty and the intersection in `(43)` were empty, the two disjoint
supports would be separate invariant subsets, contrary to transitivity.
Therefore

```text
A intersect D={p}.                                     (44)
```

Any nontrivial cycle of `X` not containing `p` lies in `A-{p}`, where `Y`
is the identity, and would itself be a proper invariant orbit. Hence `X` is
one cycle through `p`; the same argument applies to `Y`. Write these cycles
as `alpha,beta`. Their other entries are disjoint, and direct multiplication
in the convention `[X,Y]=XYX^-1Y^-1` gives

```text
[alpha,beta]=(p,alpha(p),beta(p)),                     (45)
```

up to cyclic notation. Thus the commutator is one three-cycle and fixes the
other eighteen sheets. Its ramification is `2`, whereas `(4)` has
ramification

```text
(7-1)+3(3-1)+2(2-1)=14.                                (46)
```

This is the contradiction. QED.

Equations `(41)--(44)` also force equality throughout: there are no extra
fixed sheets after all,

```text
#Fix(X)=r_0,        #Fix(Y)=r_1,
|A|=m_0,           |D|=m_1,          m_0+m_1=22.       (47)
```

Notice that this surjectivity conclusion comes **after** transitivity; the
local lemma asserted only the safe injections `(31)`.

The numerical bounds alone are sharp. On `21` letters take

```text
X=(1 2 3 4 5 6 7 8 9 10 11),
Y=(11 12 13 14 15 16 17 18 19 20 21).                 (48)
```

Their supports have sizes `11,11`, meet only at `11`, their fixed counts sum
to `20`, and they generate a transitive action. But

```text
[X,Y]=(1 12 11)                                        (49)
```

in the convention above. This hostile realizes every support premise and
shows that the inherited branch packet `(4)`, not the scalar count by itself,
is the decisive obstruction.

Combining the fixed-sheet lemma with the permutation lemma proves the theorem.

## 8. Independent critical-count route

The independent audit does not call the primary code. Away from `t=0`, use

```text
s=XT,                  p=T+s^2,                  t=p-s^2. (50)
```

With `kappa=2848/45`, and initially `tp!=0`, the two critical equations reduce
to

```text
A=-s+(p-s^2)^2p(2kappa s+Phi p+2Theta sp),
B=-6+(32/3)p-(2752/45)p^2+6kappa s^2p
    +7Phi sp^2+8Theta s^2p^2.                          (51)
```

Independent elimination gives

```text
Res_s(A,B)=-(2/373669453125)p^2 R_16(p),               (52)
deg_p R_16=16,
[p^16]R_16=34012224000000 Theta^5 W,
R_16(0)=23277095392665600000.                          (53)
```

All seventeen coefficients of this different residual are frozen in the
independent output. At `p=0`, the specialized equations in `(51)` are
`A=-s,B=-6`, so the `p^2` in `(52)` is again a degree-drop artifact, not a
finite critical point. For `p!=0`, the leading `s`-coefficients are

```text
2p(45Theta p+2848)/45,
8p(15Theta p+712)/15.                                  (54)
```

Their nonzero roots are distinct, so the homogenized equations have no
common root at `s=infinity`. Also, `t=0` in `(51)` would force `s=p=0`, which
is outside this stratum.

Before dividing by `p`, the direct `p=0` equation is

```text
2calG_p=1/s^2-6,                                       (55)
```

giving exactly two points `s^2=1/6`, `t=-1/6`, of value `1/2`. Thus `(52)`
gives `16+2=18` points with `t!=0`; restoring the two `t=0` points in `(14)`
again gives `20`. This is genuinely a different projection and explicitly
guards both omitted-coordinate strata.

## 9. Audit and boundary

Reproduce the primary and independent audits with

```bash
python3 04-computation/jc23_theta_only_extremal_seam_critical_monodromy_thm4130.py
python3 04-computation/jc23_theta_only_extremal_seam_critical_monodromy_thm4130_independent_audit.py
```

Each script uses explicit exceptions rather than `assert`. Normal,
`python3 -O`, `PYTHONHASHSEED=0`, and `PYTHONHASHSEED=8675309` executions
byte-match their frozen outputs. The two coordinate routes emit the same
semantic hash for the critical, support, branch, and verdict ledger.

This theorem empties only the smooth nonresonant theta-only exact `M=8` seam
inside the inherited reduced `(2,3)` cell. It does not treat the theta-only
collision wall `Delta_V=0`, the delta-only collision wall, the two-term
collision wall, maximum residual weight at least nine, any other pole-depth
cell, arbitrary planar Keller maps, DC(2), or JC(2) as a whole. THM-4134 later
reduces the `Delta_V=0` wall and excludes its full-boundary degrees `20,19`,
but leaves finite horizontal-carrier alternatives at degrees `16,15`.
