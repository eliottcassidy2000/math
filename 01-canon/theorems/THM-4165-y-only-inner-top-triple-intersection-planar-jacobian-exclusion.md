---
id: THM-4165
title: "Complete Y-only inner/top intersection planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
  + VERIFIED-EXACT + INDEPENDENT SOURCE-PAIR/LOCAL AUDIT. The entire Y-only
  eta=Delta=0 exact-weight-nine intersection zeta!=0, I_C=D_C=0 contains no
  nonautomorphic planar Keller pair. Its four exhaustive strata have
  (L,g)=(20,9),(19,9),(19,9),(19,8), with respectively smooth-double,
  smooth-double, triple, and ordinary-node top boundary. Eta=zeta=0 exits
  exact M=9 (MISTAKE-514); other cells, entry, M>=10, JC(2), and DC(2)
  remain OPEN.
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
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4164-y-only-triple-top-root-planar-jacobian-exclusion
script: 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165.py
output: 05-knowledge/results/jc23_y_only_inner_top_intersection_exclusion_thm4165.out
independent_audit_script: 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_y_only_inner_top_intersection_exclusion_thm4165_independent_audit.out
script_sha256: 5b12569e9cf12753026566a7886fe47419f2ba267f6c807bfe05167d11521eee
output_sha256: dbc10308c5ee706934630f125d137d579ce1912cb26382197cfac3f4061f69a2
independent_audit_script_sha256: 5deac85bcd8ccae317fc26afd35d34653ba57ce916229f0fc5c5d6605b97ddd6
independent_audit_output_sha256: a4eaf246e5f582ff3356cf85fb09f2fc6f31e0c4430d4183aee3808779ad2bb2
semantic_sha256: 65a99510c4e6cf2349ed4870a7bd3ddaf5c795d6f26ce6883236dd25567bbaa9
independent_semantic_sha256: 9a66f5791221cad1c16e59fbe890fc094e1711a35016f718d1a8c0049393c072
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact symbolic elimination constructs the irreducible coefficient
  chart and its surjective repeated-root chart, the AB and normalized
  projections, all actual first/last-row numerator and denominator gcd
  firewalls, coordinate-infinity gates, universal values, four boundary
  strata, packets, genera, and full/finite response inequalities. Normal,
  optimized, and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT within the stated independence scope. The alternative source pair
  (A,C_0), resultant rather than gcd endpoint firewalls, a separately derived
  repeated-root/local chart, node tangent-cone audit, and response arithmetic
  reproduce the theorem. It also replays the common normalized projection as
  a cross-check, but that replay is not claimed as a second independent
  normalized derivation.
weakest_link_audit: >
  The weakest inherited link is THM-4147/4155 finite-separable cubic-carrier
  and fixed-sheet transport after the new boundary normalization. The carrier
  face is unchanged, every new local branch is normalized explicitly, and the
  no-quotient +3 carrier-meridian accounting is audited below. This theorem
  remains relative to that inherited transport statement.
---

# THM-4165 -- complete Y-only inner/top intersection exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
+ VERIFIED-EXACT + INDEPENDENT SOURCE-PAIR/LOCAL AUDIT; JC(2) REMAINS OPEN.**
Work over `C` in the live `b=d=0` reduced `(2,3)` seam.

## 1. Inheritance pass and theorem

Put

```text
K0=2848/45,
C(W)=zeta W^3+Theta W^2+Phi W-1376/135,
D_C=Disc_W(C),
I_C=4 Theta K0^2-27 zeta^2.
```

> **Theorem.** The exact-weight-nine coefficient locus
>
> ```text
> eta=Delta=0,       zeta!=0,       I_C=D_C=0
> ```
>
> contains no nonautomorphic planar Keller pair.

The closest proved mechanisms are THM-4159 on `I_C=0,D_C!=0`, THM-4161 on
the nontriple part of `D_C=0,I_C!=0`, and THM-4164 on its triple-root part.
The canonical hostile is their missing codimension-two collision: a critical
endpoint and a top boundary root collide simultaneously. The corrected near
miss is therefore to transport either adjacent theorem unchanged. On this
intersection both the critical residual degree and the local boundary packet
can change, and they must be recomputed together. The least-used sidecar is
the pair

```text
(actual first/last critical rows, repeated-root residue differential),
```

followed in the finite response by the handle-subgroup orbit count. Keeping
all three pieces produces strict contradictions on every subwall.

## 2. Complete coefficient and repeated-root charts

The equation `I_C=0` with `zeta!=0` has the unique chart

```text
zeta=(5696/135)u,       Theta=3u^2,       u=135zeta/5696!=0.       (1)
```

There is no sign quotient. Direct substitution gives

```text
D_C=u F(u,Phi)/12301875,                                      (2)

F=-2076192000Phi^3+110716875Phi^2u^3
  -285684019200Phi u^2+13541904000u^5
  -61429478588416u.
```

The polynomial `F` is irreducible, hence squarefree, in `Q(u)[Phi]`, and

```text
Disc_Phi(F)
 =-4936024008000000u^2(2460375u^4+44643516416)^3.              (3)
```

Thus `Q(u)[Phi]/(F)` is a field; the quotient identities below are not being
read in a silently reducible algebra.

Let `r` be the repeated root of `C`. The equations `C(r)=C'(r)=0` give

```text
R(u,r)=11392u r^3+405u^2r^2+1376=0,
Phi=Phi_r=-2ru(2848r+135u)/45.                                (4)
```

Conversely every point of `F=0,u!=0` has a unique repeated root (also at a
triple root), so `(4)` lifts every point. The exact elimination identity

```text
Res_r(R,45Phi+2ru(2848r+135u))=-5696u^2F                      (5)
```

certifies that coverage algebraically. The polynomial `R` is irreducible in
`Q(u)[r]`. Its two partial derivatives cannot vanish together on `R=0` when
`u!=0`: they would force both `11392r+810u=0` and
`34176r+810u=0`, hence `r=0`, contrary to `R=0`. Thus `(4)` is a smooth
repeated-root chart of the whole intersection.

## 3. Four exhaustive and disjoint strata

Define the source bottom endpoint and the critical top endpoint

```text
J=8544Phi-22784u-1215u^3,

E=-3244050Phi^2u-51904800Phi u^2+4185329664Phi
  +2460375u^5-595175040u^3-22321758208u.                      (6)
```

On `(4)`, put

```text
alpha=u(5696r+135u)/45=C''(r)/2,
beta=8(356r^2+15)/45.
```

Exact reduction modulo `R` gives

```text
E(u,Phi_r)=u^2(5696r+135u)^3(356r^2+15)/15.                  (7)
```

Hence `E=0` means exactly a triple top root (`alpha=0`) or the nodal
boundary condition (`beta=0`). The three finite exceptional parameter
algebras are

```text
A_T=2460375u^4+44643516416,
Phi_T=405u^3/5696,                  r_T=-135u/5696;

B_N=36905625u^4-4721414400u^2+239958900736,
Phi_N=(2025u^2+489856)(6075u^2-489856)/(230688000u),
r_N=-(6075u^2-489856)/(170880u);

B_J=30267225703125u^8+2043284356800000u^6
    +264381824212992000u^4+6498574373014732800u^2
    +498260889496415371264,
Phi_J=u(1215u^2+22784)/8544.                                 (8)
```

These are biconditional parameter maps, not positive controls. Namely,

```text
Res_r(R,5696r+135u)       =-5696A_T,
Res_r(R,356r^2+15)        =356B_N,
Res_r(R,J(u,Phi_r))       =-(16222208/3375)u^3B_J,

R(u,r_T)=A_T/32444416,
R(u,r_N)=43B_N/(38448000u^2),
356r_N^2+15=B_N/(82022400u^2),
F(u,Phi_J)=-uB_J/8111104.                                  (9)
```

Each of `A_T,B_N,B_J` is squarefree and coprime to `u`; they are pairwise
coprime. Therefore the following four rows are exhaustive and disjoint,
including every possible intersection of `E` and `J`:

| stratum | exact condition | top boundary | `L` | `g` |
|---|---|---|---:|---:|
| smooth generic | `E*J!=0` | smooth double-root tangency | 20 | 9 |
| smooth `J` | `J=0` | smooth double-root tangency | 19 | 9 |
| triple | `5696r+135u=0` | smooth cubic tangency | 19 | 9 |
| node | `356r^2+15=0` | ordinary node, not a cusp | 19 | 8 |

In particular, `J=0` meets neither the triple nor node row, and the triple
and node rows do not meet each other.

## 4. Two source projections, one normalized projection, and exact lengths

Use

```text
t=p-s^2,
H=-3p+(8/3)p^2-(1376/135)p^3+K0s^2p^2
  +Phi sp^3+Theta s^2p^3+zeta s^3p^3,

A=(-sp+t^2H_s)/p,
C_0=s^2+2t^2H_p,
B=(C_0+sA)/t^2.
```

The source degrees and leading `s`-rows are

```text
(deg_s A,deg_s B)=(6,3),       (LC_s A,LC_s B)=(3zeta p^2,9zeta p^2),
(deg_s A,deg_s C_0)=(6,7),     (LC_s A,LC_s C_0)=(3zeta p^2,6zeta p^2).
                                                                    (10)
```

Since `zeta!=0`, no finite nonzero-`p` intersection is lost at
`s=infinity`. The primary `(A,B)` resultant and the independent `(A,C_0)`
resultant give

| stratum | `(A,B)` projection | `(A,C_0)` projection |
|---|---|---|
| smooth generic | `p^7 R_16(p)` | `p^9 S_16(p)` |
| smooth `J` | `p^8 R_15(p)` | `p^10 S_15(p)` |
| triple or node | `p^7 R_15(p)` | `p^9 S_15(p)` |

Independently form

```text
P=T+X^2T^2,       Y=XTP,
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K0Y^2+Phi P^2Y+Theta PY^2+zeta Y^3,
f=G_X/T,          h=G_T.
```

Then

```text
(deg_X f,deg_X h)=(8,9),       LC_X(f)=LC_X(h)=9zeta T^8,    (11)

Res_X(f,h)=T^56(6T+1)^2 Q_d(T),
d=16 on the smooth-generic row and d=15 on all three exceptional rows.
                                                                    (12)
```

Thus there is likewise no unrecorded `X=infinity` loss. Before exceptional
specialization, the actual first and last source rows on `F=0` are

```text
[p^16]R=(1534934573591985913856/622782421875)u^6E,
R(0)   =(8305770496/1125)u^2J.                                (13)
```

For the normalized projection they are

```text
[T^16]Q=(47966705424749559808/38306957530517578125)u^5N,
Q(0)=-(1519777094677765052956672/56953125)u^7,

N=(108000/u)EJ^2                  in Q(u)[Phi]/(F).            (14)
```

The irreducibility gate in Section 2 makes `(14)` a field identity. It also
shows why both `E=0` and `J=0` must be strict-transformed.

The terminal assertion is all-roots exact. Extracting the *actual* first and
last coefficients after reduction in each exceptional wall algebra gives,
up to a nonzero rational scalar,

```text
wall       source AB and AC0: (first,last)
           normalized:        (first,last)

triple    (u(405u^2+88064), u(2025u^2-979712))
          (u(36585u^2+427154432), u^3)

node      ((20614295925117675u^2-789558891060057728)/u,
            u(247725u^2+21063808))
          ((1199320411266696917025u^2+56388834945845741080192)/u,
            u(153378225u^2-12933178112))

J         (u(18030796678125u^6+5214383198496000u^4
              +118847166205132800u^2+14477208077889175552),
            u^3(2460375u^4-204543360u^2+5580439552))
          (u(112249545005345596875u^6+4039772647342545648000u^4
              +349920353377066350182400u^2
              +9283134759621226372530176),
            u^7).                                                (15)
```

For each of the twelve displayed actual coefficients, its numerator and
denominator are separately nonzero and coprime to the corresponding one of
`A_T,B_N,B_J`. The primary certificate proves this by Euclidean gcd; the
source-pair/local referee uses nonzero resultants. Hence the residual degrees
in `(12)` and the source table are exact at every exceptional root.

The actual critical ideal restores the four universal points

```text
T=0,     X^2=-6,       G=0,       det Hess(G)=+6,
T=-1/6,  X^2=6,        G=1/2,     det Hess(G)=-6.              (16)
```

They contribute four points. Under a hypothetical Keller realization, the
inherited Hessian congruence makes the full affine critical scheme Morse;
projected-coordinate collisions therefore do not reduce its length. Thus
`L=d+4`, giving exactly the four lengths in Section 3. No residual
discriminant hypothesis is used.

## 5. Boundary normalization, packets, and genera

For a generic fibre parameter `q`, clear its equation as

```text
F_q=(s^2-p)(q-H)-s^2/2,
z=1/p,       W=s,       K=z^4F_q(W,z^-1).
```

At the repeated root put `v=W-r`. The exact rows needed for the local
normalization are

```text
K=alpha v^2+zeta v^3+beta z+gamma vz+delta z^2
  +terms of higher relevant order,

gamma=5696r/45,
delta=-(2848r^4+120r^2+135)/45.                              (17)
```

On a smooth-double row, `alpha*beta!=0`, so `z` has order two in `v`.
Since `K_z` is a unit, the Keller residue differential

```text
omega=-z^2 dW/K_z
```

has order four and produces index `5`; the remaining simple top root has
index `3`. On the triple row, `alpha=0,beta!=0`, so `z` has order three,
`ord(omega)=6`, and the single top place has index `7`.

On the node row set `y=u/r`. The equations `R=beta=0` reduce exactly to

```text
R_y=91125y^2+2563200y+174388736=0,
alpha=-y(135y+5696)/1068,       delta=-3,
gamma^2-4alpha delta=-(135y+2848)^2/12015.                   (18)
```

The polynomial `R_y` is coprime to `y`, `135y+5696`, and `135y+2848`.
Thus `alpha` and the tangent discriminant are nonzero. The tangent cone in
`(v,z)` has two distinct lines, so this is an ordinary node and no cusp
survives on the inner wall. On both branches `z` has order one and `K_z`
has order one; each therefore gives index `2`. The remaining simple root
again gives index `3`.

The non-top packet inherited from THM-4155 is `(8,2,2,2,1)`. The Newton
polygon remains

```text
(0,1),(2,0),(5,3),(3,4),(0,4),       (2Area,B,I)=(27,11,9).
```

The smooth double and triple tangencies have delta zero. The ordinary node
has delta one. Consequently the complete packets are

| stratum | complete packet | defect | genus |
|---|---|---:|---:|
| smooth generic or `J` | `(8,5,3,2,2,2,1)` | 16 | 9 |
| triple | `(8,7,2,2,2,1)` | 16 | 9 |
| node | `(8,3,2,2,2,2,2,1)` | 14 | 8 |

In every row the defect is `2g-2`. This both saturates Riemann--Hurwitz and
rules out an omitted affine branch contribution.

## 6. Carrier-orbit lemma and both monodromy responses

The prime cubic carrier remains the unchanged face

```text
q-1/2=K0W^2+zeta W^3.                                      (19)
```

Because `zeta!=0`, it is a separable degree-three extension of `C(q)`, hence
prime. THM-4120 sends every rational infinity place to the target origin.
The THM-4147/4155 finite-separable transport therefore separates the carrier
into three target punctures on the same finite permutation action.

We use the following lemma with all inequalities and directions explicit.

> **Carrier-orbit lemma.** Let `n>m+1`, and suppose
> `X,Y,tau_1,...,tau_m` generate a transitive permutation action on `n`
> letters, with each `tau_j` a transposition. Put
>
> ```text
> a=|supp X|,       b=|supp Y|,
> U=supp X union supp Y,       k=|supp X intersect supp Y|.
> ```
>
> If transported fixed sheets give `a+b<=2n-L`, then
>
> ```text
> |U|>=n-m,       k<=n+m-L.                                (20)
> ```

**Proof.** Let `H=<X,Y>`. Adjoining one transposition merges at most two
current `H`-orbits, hence lowers the orbit count by at most one. Transitivity
after all `m` additions gives

```text
#Orb(H)<=m+1,       so n-#Orb(H)>=n-m-1.                     (21)
```

The hypotheses `n>m+1` and transitivity imply `U` is nonempty; otherwise
`m` transpositions could merge at most `m` of the `n` singleton gaps. Every
point outside `U` is an `H`-fixed singleton, while nonempty `U` contributes
at least one more orbit. Therefore

```text
#Orb(H)>=(n-|U|)+1,       so n-#Orb(H)<=|U|-1.               (22)
```

Combining `(21)` and `(22)` in the displayed directions gives `|U|>=n-m`.
Finally

```text
k=a+b-|U|<=(2n-L)-(n-m)=n+m-L.
```

This proves `(20)`. QED.

For the punctured elliptic target, the exact relation is

```text
[X,Y] mu_O tau_1 ... tau_m=1.
```

The inherited commutator bound and permutation-index triangle inequality
give

```text
ind([X,Y])<=2k,       ind(mu_O)<=2k+m.                       (23)
```

In the full response, `X,Y` themselves generate transitively, so `U` is all
`n` sheets and `k<=n-L`. In the finite response, the carrier packet
`(2,2,2)` is removed from the origin packet and becomes `m=3` distinct
transposition punctures. There is no quotient loss in the `+3`: each
separated meridian acts on the same finite cover and has Cayley index one.
Any cancellation in their product only decreases the upper bound `(23)`.

The exact response ledger is

| stratum | `L` | full `(n,ind(mu_O),2(n-L))` | finite origin packet | finite `(n,ind(mu_O),k_max,2k_max+3)` |
|---|---:|---|---|---|
| smooth generic | 20 | `(23,16,6)` | `(8,5,3,1)` | `(17,13,0,3)` |
| smooth `J` | 19 | `(23,16,8)` | `(8,5,3,1)` | `(17,13,1,5)` |
| triple | 19 | `(22,16,6)` | `(8,7,1)` | `(16,13,0,3)` |
| node | 19 | `(22,14,6)` | `(8,3,2,2,1)` | `(16,11,0,3)` |

Here `k_max=n+3-L`. In every finite row `n` is `16` or `17`, so the lemma's
load-bearing hypothesis `n>m+1=4` holds. Every full ceiling is strictly below
the listed full origin index, and every finite ceiling is strictly below its
finite origin index. Thus neither response exists on any of the four
exhaustive strata, proving the theorem.

## 7. Weakest link, scope, and replay

The weakest inherited link is the finite-separable carrier and fixed-sheet
transport from THM-4147/4155 after a changed boundary normalization. This
proof does not hide that dependency: `(17)--(18)` normalize every new local
branch, `(19)` shows the carrier face itself is unchanged, and Section 6
audits its three meridians without quotienting them. The source-pair/local
referee independently checks those new local inputs. Its normalized
elimination is a replay of the common coordinate model and is not advertised
as an independent normalized-projection proof.

This theorem closes exactly

```text
eta=Delta=0,       zeta!=0,       I_C=D_C=0.
```

It makes no claim on the lower-weight filtration exit `eta=zeta=0`, another
reduced cell, seam entry, exact residual weight at least ten, `JC(2)`, or
`DC(2)`.

Primary replay:

```text
python3 -B 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165.py
python3 -B -O 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165.py
PYTHONHASHSEED=241 python3 -B 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165.py
```

Independent source-pair/local replay:

```text
python3 -B 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165_independent_audit.py
python3 -B -O 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165_independent_audit.py
PYTHONHASHSEED=251 python3 -B 04-computation/jc23_y_only_inner_top_intersection_exclusion_thm4165_independent_audit.py
```

Both normal/optimized/hash-seeded outputs byte-match the canonical artifacts.

**QED.**
