---
id: THM-4289
title: "A23 blowdown observers and the Kahler-dualizing quotient"
status: >
  PROVED FORMAL-LOCAL RELATIVE TO THM-4272/4279/4280/4284/4288 +
  FINITE-EXACT AUDIT PASS. For every resolution triple
  D_s: b*z*(z-b^s)=0, equal-value function descent has an exact length-s
  obstruction, and ambient Kahler extension of normalized regular one-forms
  has the differentiated length-s obstruction. Regular dualizing forms impose
  no such condition: their exact loss is conductor/J_f ~= k[[b]]/(b^s), with
  hostile b*z*eta=(0,0,db). At the A_m endpoint the analogous quotient has
  length m-1. Blowdown predicates factor through a linear observer exactly
  when its kernel lies below the requested contact depth; this gives sharp
  failures of the complexified 124 and 247 observers and a finite order-state
  quotient. No ambient-Kahler membership for the Keller differential, raw
  descent, or JC(2) is proved.
source: root/planar-jc-higher-order-20260830
depends_on:
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4279-four-channel-formal-log-hasse-observer-for-e0-hom-at-fat-contact
  - THM-4280-integral-three-channel-fat-contact-observer-and-sharp-five-jet-bound
  - THM-4284-a23-conductor-defect-and-degree-shell-first-character-nondescent
  - THM-4288-a23-partial-normalization-relative-differential-and-etale-base-change-obstruction
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel
  - THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
external_method_source: >
  Lingyuan Ye and Yiqi Xu, "Failure of Higher-Order Truth within
  Intuitionistic Propositional Logic," arXiv:2608.26874v1, is used only for
  the observer-congruence proof-design lens. Its theorem does not transfer to
  the Jacobian problem. The paper and two repairable v1 typos are audited in
  05-knowledge/reference/CORE-PAPERS-HEYTING-HIGHER-ORDER-TRUTH-2026-08-29.md.
primary_script: 04-computation/jc23_a23_blowdown_kahler_dualizing_observer_thm4289.py
primary_output: 05-knowledge/results/jc23_a23_blowdown_kahler_dualizing_observer_thm4289.out
primary_script_sha256: 47ae505c1569ce96008855c5c6025b9b6d2a826d00372f054ba7e939ce888c8b
primary_output_sha256: 157c94029c1990947d02a9738863de991dcd225af94422443c5baf417124e283
hash_basis: raw LF bytes
audit: >
  PASS. A dependency-free finite-jet audit over F_1000000007 checks the two
  obstruction maps, their kernels/ranks and differentiation compatibility for
  every 1<=s<=11 at N=2s+5; finite cyclic conductor/Jacobian quotients; the
  A_12 endpoint; nonvacuous dualizing lifts; minimal hostiles; and inherited
  observer-matrix kernels. State sets and partial-normalization lengths remain
  dependency results, printed as such. Ordinary and optimized output streams
  are byte-identical.
---

# THM-4289 -- A23 blowdown observers and the Kahler-dualizing quotient

**PROVED FORMAL-LOCAL RELATIVE TO THM-4272/4279/4280/4284/4288 +
FINITE-EXACT AUDIT PASS. THE KELLER AMBIENT-FORM BRIDGE AND `JC(2)` REMAIN
OPEN.**

## 1. Exact function descent on every resolution triple

Let `k` be a characteristic-zero field, `R=k[[b]]`, and, for `s>=1`, put

```text
D_s=Spec k[[b,z]]/(b*z*(z-b^s)).                       (1)
```

Its normalized branches are

```text
E: b=0,                    R_0: z=0,                  C: z=b^s.  (2)
```

Write

```text
V_s={(e(z),r(b),c(b)):e(0)=r(0)=c(0)=a}.              (3)
```

Define

```text
Phi_s(e,r,c)=[c(b)-r(b)-e(b^s)+a]
             in bR/b^(s+1)R.                          (4)
```

Then restriction to the normalization has the exact sequence

```text
0 -> O_(D_s) -> V_s --Phi_s--> bR/b^(s+1)R -> 0.      (5)
```

In particular, equal-value branch functions descend to `D_s` if and only if

```text
c-r-e(b^s)+a in b^(s+1)R.                             (6)
```

To prove this, start with

```text
F_0(b,z)=e(z)+r(b)-a.                                  (7)
```

It has the required restrictions on `E` and `R_0`. Its error on `C` is the
representative in `(4)`. A function vanishing on both `E` and `R_0` is a
multiple of `bz`, and on `C` this restricts to `b^(s+1)`. Conversely every
multiple of `b^(s+1)` is supplied by `bzH(b)`. This proves exactness at
`V_s`; the triples `(0,0,b^j)`, `1<=j<=s`, prove surjectivity.

The case `s=1` is exactly THM-4288's final circuit
`c_C'-c_R'-c_E'`. Formula `(4)` recovers every higher blowdown character,
not just its first-order shadow.

## 2. Ambient Kahler forms are the differentiated obstruction

Consider arbitrary regular normalized branch forms

```text
e(z) dz,                    r(b) db,                  c(b) db.  (8)
```

They are restrictions of one ambient form

```text
A(b,z) db+B(b,z) dz                                       (9)
```

if and only if

```text
Psi_s(e,r,c)
 =[c-r-s*b^(s-1)e(b^s)] in R/b^sR
 =0.                                                     (10)
```

Thus there is an exact cokernel

```text
coker(Omega^1_(k[[b,z]]) ->
      Omega^1_E direct-sum Omega^1_R0 direct-sum Omega^1_C)
 ~=R/b^sR.                                              (11)
```

Indeed, `r(b)db+e(z)dz` supplies the `E` and `R_0` restrictions. Its
restriction to `C` is

```text
[r(b)+s*b^(s-1)e(b^s)]db.                              (12)
```

A correction which vanishes on `E` and `R_0` changes the `C` restriction by
a multiple of `b^s`; conversely `zH(b)db` supplies every such multiple. This
proves `(10)--(11)`.

Differentiation gives an isomorphism

```text
d:bR/b^(s+1)R -> (R/b^sR)db,
  sum_(j=1)^s a_j b^j |-> sum_(j=1)^s j*a_j b^(j-1)db. (13)
```

Characteristic zero is load-bearing. If the branch functions in `(3)` are
formal coordinates of maps to a smooth curve, then `(13)` takes `Phi_s` to
`Psi_s` applied to their differentials. Consequently:

> **Differential descent criterion.** After translating an `E_0`-valued
> branch map to the origin and applying the formal logarithm, equal-value map
> descent to `D_s` is equivalent to ambient-Kahler extension of its pulled-back
> invariant differential.

This is an iff criterion on the displayed formal-local object. In
characteristic `p`, `b^p` is a sharp warning: ordinary differentiation loses
its descent obstruction.

## 3. The exact dualizing firewall

Put

```text
f=b*z*(z-b^s),
eta=Res(db wedge dz/f).                                (14)
```

With the sign convention in `(14)`, its normalized restrictions are

```text
eta_E=-dz/z^2,
eta_R0=-db/b^(s+1),
eta_C= db/b^(s+1).                                    (15)
```

The conductor of `O_(D_s)` in its normalization and the plane Jacobian ideal
are

```text
c_s=(z(z-b^s), b(z-b^s), bz),
J_f=(f_b,f_z) subset c_s.                              (16)
```

Inside rational differentials,

```text
nu_*Omega^1_(D_s normalization)=c_s*eta,
im(Omega^1_(D_s)/torsion)=J_f*eta.                    (17)
```

The first equality follows directly from `(15)`: a multiplier gives regular
forms on the three branches exactly when its restrictions are divisible by
`z^2,b^(s+1),b^(s+1)`, which is precisely the conductor. The second follows
from

```text
db=f_z*eta,                    dz=-f_b*eta.             (18)
```

There is an exact quotient

```text
c_s/J_f ~= R/(b^s),                                    (19)
```

generated by the class of `bz`. To see it without a length formula, set

```text
a_0=z(z-b^s),       r_0=b(z-b^s),       c_0=bz.
```

Then

```text
f_b=a_0-s*b^(s-1)c_0,             f_z=r_0+c_0.         (20)
```

Modulo `J_f`, therefore, `a_0=s*b^(s-1)c_0` and `r_0=-c_0`. The conductor
relations give `b^s c_0=0`. Sharpness is certified by the normalized-component
functional

```text
h |-> [(h_C+h_R+s*h_E(b^s))/b^(s+1)] mod b^s,          (21)
```

which kills `J_f` and sends `bz` to `1`.

The smallest hostile is consequently

```text
bz*eta=(0,0,db).                                       (22)
```

It has a regular numerator and regular pullback on every normalized branch,
so it is an honest dualizing/Poincare-residue section. But `(10)` sends it to
`1`; it is not ambient Kahler. Therefore

```text
regular dualizing residue  does not imply  ambient-Kahler descent. (23)
```

The exact missing sidecar is not the word "regular" but the ideal membership

```text
h in J_f=(f_b,f_z),                                    (24)
```

or, equivalently, vanishing of `(21)`. For this weighted-homogeneous plane
curve, `delta=s+2` and `mu=tau=2s+2`, so the numerical identity
`length(c_s/J_f)=tau-delta=s` agrees with the explicit module proof; the
numerical identity is not needed for it.

## 4. The raw two-branch endpoint and the complete A23 ladder

For comparison, return to THM-4284's raw contact

```text
A_m=k[[b,q]]/(q(q-b^m)),                               (25)
```

with branches `R_0:q=0` and `C:q=b^m`. Normalized forms
`r(b)db,c(b)db` come from an ambient form in `k[[b,q]]` if and only if

```text
c-r in b^(m-1)R.                                       (26)
```

For equal-value branch functions, raw descent is `L_C-L_R in b^mR`.
Differentiation identifies these two conditions in characteristic zero.

If

```text
eta_m=Res(db wedge dq/[q(q-b^m)]),                     (27)
```

then

```text
eta_m|R_0=-db/b^m,               eta_m|C=db/b^m,
c_m=(q-b^m,q),                   J_m=(f_b,f_q),
c_m/J_m ~=R/(b^(m-1)).                                 (28)
```

The quotient is generated by `q`, and

```text
q*eta_m=(0,db)                                         (29)
```

is the endpoint dualizing hostile. Thus for `m=12`, ambient Kahler forms see
exactly the eleven post-value characters which dualizing regularity erases.

Blowing up `(25)` once in the chart `q=bz` gives the reduced total transform
`D_(m-1)`. Successive blowups give

```text
D_(m-1),D_(m-2),...,D_1.                               (30)
```

If the maps on `E` and `R_0` are constant and the formal log on `C` is `L`,
then `(6)` becomes

```text
descent to D_s iff L in b^(s+1)R.                     (31)
```

Read in increasing `s`, these are the successive conductor characters. The
hostile `L=b^s` passes every earlier `D_j`, `j<s`, and first fails at `D_s`.
At `A_12`, the actual orders `1,2,4` from THM-4280 therefore first fail at
`D_1,D_2,D_4`, respectively, while `D_11` recovers full raw order-twelve
descent. There is no additional unexplained terminal condition.

## 5. Observer congruence and the Ye--Xu transfer boundary

Let `V` be any linear subspace of based length-twelve formal logs and put

```text
F^jV=V intersect b^jR/b^12R,
P_j(v) iff v in F^jV.                                  (32)
```

For a linear observer `q:V->W`, the exact fibre-purity test is

```text
P_j is constant on q-fibres iff ker(q) subset F^jV.    (33)
```

Necessity compares each kernel vector with zero. For sufficiency, two points
in one fibre differ by a kernel vector, and membership in the subspace `F^jV`
is then unchanged. The containment in `(33)`, not equality with `F^jV`, is
the correct condition.

Apply `(33)` to THM-4280's two arithmetically minimal triples on the
chosen-embedding complexification `M_sigma`:

```text
q_124=(c_1,c_2,c_4),       ker q_124=C(kappa f+g),     ord=7;
q_247=(c_2,c_4,c_7),       ker q_247=C(kappa f-g),     ord=1.  (34)
```

Therefore

```text
q_124 preserves P_j exactly for j<=7 and first fails at P_8;
q_247 preserves only P_1 and first fails at P_2.       (35)
```

The zero vector paired with the corresponding kernel vector is sharp. In
the curve notation of Section 4, `P_(s+1)` is exactly descent to `D_s`.
Thus the `q_124` hostile first fails on curve `D_7`, while the `q_247`
hostile first fails on curve `D_1`. These vectors exponentiate to formal
`E_0` germs but are not actual global `C_0->E_0` maps. On the actual
Eisenstein Hom lattice, THM-4280 proves both triples injective, so every
`P_j` factors through either one.

For the whole predicate family the coarsest order quotient is

```text
tau(v)=min(12,ord_b(v)),             P_j(v) iff tau(v)>=j. (36)
```

Its exact state sets on the four relevant consumers are

| consumer | order states |
|---|---|
| arbitrary based `12`-jets | `1,2,3,4,5,6,7,8,9,10,11,12` |
| chosen-embedding complex Hom span | `1,2,4,7,12` |
| actual `O`-Hom lattice | `1,2,4,12` |
| degree-`34/42` shells | `1` |

Here state `12` includes zero/full order-twelve descent. On the final row,
THM-4280 gives `c_1!=0` for every candidate, so common-tangent descent is
already uniformly false. Refining its finite observer cannot repair the
Keller gap; the missing assertion is that the specialized Keller response is
governed by the raw/ambient consumer at all.

This is the typed application of the finite-observer mechanism isolated in
the Ye--Xu paper audit. The fixed `A_23` filtration is a finite chain, with no
infinite antichain and no completion escape. The paper's higher-order topos
obstruction does **not** activate. What transfers is only the discipline that
a quotient must preserve both its successor/blowdown operation and its target
predicate.

## 6. A dimension firewall for partial descent

For `1<=ell<m`, put

```text
J_(ell,m)=b^ellR/b^mR.                                 (37)
```

On this arbitrary branch-difference space, full descent to `A_m` is the zero
predicate. A linear observer determines that predicate if and only if it is
injective, so at least `m-ell` scalar channels are necessary and sufficient.
THM-4288 gives instead

```text
length Omega_(A_ell/A_m)=min(m-ell,2ell).              (38)
```

Thus cotangent-length-many scalar channels are dimensionally insufficient
whenever `m>3ell`. At `m=12`:

| `ell` | channels needed for arbitrary full descent | relative-`Omega` length |
|---:|---:|---:|
| 1 | 11 | 2 |
| 2 | 10 | 4 |
| 4 | 8 | 8 |

At `ell=4` the dimension no-go disappears, but no positive factorization
follows. This explains why relative differentials prove ramification while
usually failing to reconstruct the full conductor gluing.

## 7. Keller consequence and strict scope

THM-4288 left two possible local bridges: extend the resolved map over a
surface neighbourhood of the exceptional curve, or prove that its invariant
differential is ambient Kahler on the raw terminal chart. Equations
`(17)--(24)` sharpen the second to one explicit ideal-membership test.

They also refute a tempting shortcut. A Poincare residue identity naturally
lands in the dualizing module. Even a residue with a regular conductor
numerator need not land in `J_f*eta`; `(22)` is the minimal witness. Therefore
dualizing, logarithmic, reflexive, and normalized-branch regularity cannot be
silently promoted to raw ambient-Kahler regularity.

The cheapest unresolved local test is now:

```text
restrict the Keller pullback differential to the terminal triple;
write it as h*eta;
test the single class [h] in c_s/J_f.                  (39)
```

Vanishing closes the corresponding blowdown layer; nonvanishing identifies
the exact hostile character. THM-4290 subsequently supplies an independent
cyclic-equivariance bypass on the exact-weight-twelve lane. It does not turn
the dualizing form in this theorem into an ambient one.

No claim here constructs `(39)`, descends the resolved Keller map to the raw
`A_23` surface, proves a Cartesian etale pullback, or proves `JC(2)`.

## 8. Exact audit

The optimization-safe script

```text
04-computation/jc23_a23_blowdown_kahler_dualizing_observer_thm4289.py
```

constructs the truncated maps `Phi_s` and `Psi_s`. For every `1<=s<=11` at
`N=2s+5` it verifies their surjectivity, that the ambient restriction images
are their full kernels, and that differentiation intertwines them. It also
recovers

```text
rank(function restrictions)=3N-2-s,
rank(ambient-form restrictions)=3N-s,
length(c_s/J_f)=s.                                    (40)
```

It separately verifies cyclic quotient bases generated by `bz` and `q`, their
dual quotient functionals and annihilating powers, the `A_12` endpoint length
eleven, every nonzero dualizing numerator visible in the truncation, the two
minimal hostiles, and the exact observer-matrix kernels and failure depths in
`(34)--(35)`. The order-state sets and Section
6 lengths are inherited from THM-4280/4288 and are only labelled dependency
readouts in the output. Its frozen output is

```text
05-knowledge/results/jc23_a23_blowdown_kahler_dualizing_observer_thm4289.out
```

**QED.**
