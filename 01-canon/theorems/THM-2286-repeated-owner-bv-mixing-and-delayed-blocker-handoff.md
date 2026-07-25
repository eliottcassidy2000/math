---
id: THM-2286
title: "Repeated-owner BV mixing and delayed blocker handoff"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT HOSTILE AUDIT ACCEPTED. In every
  one of the fifteen repeated-first
  scalar profiles (1,1,c), 5<=c<=19, some labelled exclusive blocker
  owner has mass at least 5229541/593783190. The selected owner has
  explicit depth lambda_j in {1,c}; the guard-safe, unit-mask-free
  residual outside that same owner has mass at least 2593/90090.
  If S is the sum of the nine scalar coefficients, then at every time
  k no earlier than that owner's expiration and satisfying
  13^k >= (2375132760/5229541)S, a set of mass at least
  13560199813/106987855174200 travels from the exclusive owner to a
  state serviced by another blocker and by neither the marked owner,
  guard absorber, nor any unit mask. Thus, together with THM-2288, every
  one of the 165 remaining first-depth-one profiles forces a positive
  ancestry-compatible delayed blocker-only handoff. This excludes no
  profile and does not place the handoff at prescribed expiration.
source: codex-2026-07-25-repeated-owner-bv-handoff
depends_on:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
related:
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2268-two-shell-private-owner-trident-and-raw-carry-cocycle-no-go
  - THM-2271-expiration-support-forces-a-weighted-owner-absorber-cut
  - THM-2288-shallow-owner-bv-mixing-and-delayed-blocker-handoff
script: 04-computation/lrc14_repeated_owner_bv_delayed_handoff_thm2286.py
output: 05-knowledge/results/lrc14_repeated_owner_bv_delayed_handoff_thm2286.out
script_sha256: 129939aa70a5e838f314cc3738cb2d656c3c4818ffc453b970b8e4f51da721af
output_sha256: 64788559e0c93d142e59979a8e1f446fa206b164cdbdcb3b20463155adaf2171
hash_basis: working-tree bytes (LF)
---

# THM-2286 -- delayed blocker-only handoff on every repeated row

Use the scalar five-unit/three-blocker notation

```text
T(x)=13x mod 1,

D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},

A_0=C_H minus union_(i=1)^5 D_(q_i).                 (1)
```

Assume the scalar cover

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(r=1)^3 D_(c_r)             (2)
```

almost everywhere and one of the fifteen repeated-first profiles

```text
c_1=13u_1,        c_2=13u_2,        c_3=13^c u_3,

5<=c<=19.                                             (3)
```

Here `H,q_1,...,q_5,u_1,u_2,u_3` are thirteen-units, `H` is odd, and
the usual distinctness assumptions hold. Put

```text
lambda_1=lambda_2=1,             lambda_3=c,

S=H+sum_(i=1)^5 q_i+sum_(r=1)^3 c_r.                 (4)
```

For a blocker label `j`, define its exact exclusive-owner set and its
blocker-only target by

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(r!=j)D_(c_r),

R_j=A_0 minus D_(c_j).                               (5)
```

Then there is a label `j in {1,2,3}` such that, for every integer `k`
satisfying

```text
k>=lambda_j+1,

13^k >= (2375132760/5229541)S,                       (6)
```

the ancestry-aware delayed return set

```text
H_(j,k)=E_j intersection T^(-k)(R_j)                 (7)
```

obeys

```text
measure(H_(j,k))
 >=13560199813/106987855174200
 >0.                                                 (8)
```

Almost every source point in (7) is serviced by `c_j` alone among all
nine labels. Almost every target point is serviced by one or both of the
other blockers and by none of the guard absorber, five unit masks, or
`c_j`. Thus (8) is a positive-mass, ancestry-compatible
blocker-to-other-blocker handoff.

The selected label is not asserted to be shallow. Its depth is retained
exactly by (4):

```text
j=1 or 2  -> lambda_j=1,
j=3       -> lambda_j=c.                             (9)
```

This is why the first condition in (6), rather than a universal time two,
is part of the theorem.

## 1. A labelled source exists on every repeated profile

THM-2263 proves the compatible repeated-row owner ledger

```text
sum_(j=1)^3 measure(E_j)
 >=5229541/197927730.                                (10)
```

The bound is uniform over all `c=5,...,19`; its unique worst profile for
that ledger is `(1,1,5)`, with the compatible equality carrier given by
common dilates of

```text
(1,2,13^4)       or       (1,2,2*13^4)               (11)
```

after the common factor thirteen is removed.

Choose a label `j` of maximal exclusive mass. Then

```text
measure(E_j)
 >=e_0:=5229541/593783190.                           (12)
```

This pigeonhole step is intentionally label-faithful. No union image is
substituted for the chosen `E_j`, and (9) records the actual expiration
time of that label.

## 2. The same label has a positive blocker-only target

For a blocker `q` of exact thirteen-adic depth `d>=1`, reduce the
guard/blocker pair as

```text
A=H/gcd(H,q),          B=q/gcd(H,q)=13^d v.          (13)
```

Then `A` is odd, `gcd(A,v)=1`, and neither `A` nor `v` is divisible by
thirteen. The exact fold in THM-2080 gives the parity-sharp cap

```text
measure(C_H intersection D_q)<=G_d,

G_d=
  5/49+5/(49*13^d),             d odd,
  5/49+5/(294*13^d),            d even.              (14)
```

For completeness, the finite part behind (14) has no profile
assumption. THM-2080's elementary square bound leaves only `Av<=2` in
the odd case and `Av<=14` in the even case. Exact evaluation gives the
unique endpoints `(A,v)=(1,1)` and `(1,6)`, respectively. In particular,
over every possible selected depth in (9),

```text
G_(lambda_j)<=G_1=10/91.                             (15)
```

Since `A_0 subset C_H`, equations (5), (15), and THM-2263's residual
floor give

```text
measure(R_j)
 >=measure(A_0)-measure(C_H intersection D_(c_j))

 >=961/6930-10/91

 =eta_0:=2593/90090.                                 (16)
```

This works for the same label selected in (12), including when `j=3`.
The depth-one cap is merely the uniform worst debit; a selected deep
owner has a strictly better target floor.

The target `R_j` retains the information needed later. By definition it
lies inside `C_H`, outside all five unit masks, and outside `D_(c_j)`.
The global cover (2) therefore makes another blocker available almost
everywhere on `R_j`.

## 3. BV mixing preserves ancestry

Let `P` be the normalized Perron operator of `T`:

```text
P f(y)=(1/13)sum_(r=0)^12 f((y+r)/13).               (17)
```

For every circle-BV function,

```text
Var(Pf)<=Var(f)/13.                                  (18)
```

Indeed, apply the triangle inequality along a cyclic partition. The
thirteen lifted inverse-branch partitions occupy the thirteen consecutive
subintervals of the source circle, so their total variation is at most
`Var(f)`; the normalization supplies the factor `1/13`. Iterating,

```text
Var(P^k f)<=13^(-k)Var(f).                           (19)
```

The operator preserves integrals. If a circle-BV function `g` has mean
zero, then

```text
||g||_infinity<=Var(g).                              (20)
```

Apply (19)--(20) to `f=1_(E_j)`:

```text
P^k 1_(E_j)
 >=measure(E_j)-13^(-k)Var(1_(E_j))                  (21)
```

almost everywhere.

Every boundary point of `E_j` belongs to the combined endpoint bank of
the guard, five unit combs, and three blockers. A comb with coefficient
`a` has at most `2a` boundary points, and the guard has at most `2H`.
Consequently

```text
Var(1_(E_j))<=2S.                                    (22)
```

The numerical factor in (6) is exactly

```text
4/e_0=2375132760/5229541.                            (23)
```

Thus (6), (12), and (21)--(23) imply

```text
P^k 1_(E_j)>=e_0/2                                  (24)
```

almost everywhere.

Finally, the exact transfer identity is

```text
measure(E intersection T^(-k)R)
 =integral_R P^k 1_E dmeasure.                       (25)
```

Use (16) and (24):

```text
measure(H_(j,k))
 >=(e_0/2)eta_0

 =13560199813/106987855174200,                       (26)
```

which proves (8). The inequality holds at every sufficiently large
integer time, not merely along a subsequence.

To obtain literal labelled transition pieces, remove the null exceptional
set of (2) and partition each target point by the least available label
among the other two blockers. Their two masses sum to (26), so one target
label receives at least half that amount. This refinement is optional:
the intrinsic cut

```text
{c_j} | {the other two blockers}                     (27)
```

already has the full ancestry weight in (26).

## 4. Exact stopping boundary at prescribed expiration

THM-2263 also gives the chosen owner's support floor at its prescribed
expiration:

```text
measure(T^(lambda_j+1)(E_j))
 >=(169/20)e_0

 =5229541/70270200

 =1/14+210241/70270200

 <1/7.                                               (28)
```

Thus the repeated branch crosses one oriented half-comb threshold but not
one full danger-comb capacity. Equation (28) does **not** imply

```text
T^(lambda_j+1)(E_j) intersection R_j
has positive measure.                               (29)
```

The BV proof establishes (29) only after the coefficient-dependent mixing
horizon (6). THM-2261's expiration-surjectivity example explains why the
fixed-time implication cannot be recovered from raw owner support, and
the finite root-word controls in THM-2288 explain why terminal phase or
carry cannot be dropped when trying to transport the delayed handoff back
to expiration.

Accordingly the theorem proves neither a scalar-profile exclusion nor
LRC(14). Its exact gain is the removal of a branch gap:

```text
150 strict profiles (THM-2288)
 + 15 repeated-first profiles (this theorem)
 = 165 first-depth-one profiles

all force a positive delayed blocker-only handoff.                 (30)
```

What remains is to charge that handoff against the cover, or to select it
at a bounded ancestry/gap horizon. An exact ordered selector state,
rather than a static owner-overlap graph, is the natural retained object.

## 5. Exact verification

The companion checks the fifteen-row census, every possible selected-owner
depth, the parity-sharp target caps, the source and target floors, the BV
horizon normalization, the positive handoff constant, and the
half-comb/full-comb stopping inequalities in (28). Reproduce with

```bash
python3 04-computation/lrc14_repeated_owner_bv_delayed_handoff_thm2286.py
python3 -O 04-computation/lrc14_repeated_owner_bv_delayed_handoff_thm2286.py
```

Normal and optimized transcripts are byte-identical to the stored output.
QED.
