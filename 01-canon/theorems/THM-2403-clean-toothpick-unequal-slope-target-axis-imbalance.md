---
id: THM-2403
title: "Clean toothpicks force a global fully-all-safe unequal-slope target current"
status: >
  REPAIRED PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT REAUDIT
  PENDING. In the standard owner-pivot packet, omit q_* from the pivot
  and choose two distinct ordinary helpers. The resulting lawful target
  dipoles move one ordinary unit factor with each target blocker. The
  global packet retaining all six unit-safe and all three blocker-safe
  present factors is zero at the untwisted slice by the scalar cover.
  On the positive clean-toothpick set, one unequal-slope orbit has
  shifted mass at least nine per parent sheet; the deep successor probe
  doubles this before the 1/13 root Jacobian. Consequently every one of
  the twelve nonzero first-target colours survives, and for each colour
  some exact fixed triangle has nonzero target, a nonzero deep colour,
  and gcd(m,91)=1. This is a lawful fully-all-safe present packet, not
  the preselected positive shallow-owner current. It does not give an
  all-coordinate 91-unit address, preserve a previously marked triangle
  or owner phase, exclude a row, or prove LRC(14).
source: codex-2026-07-26-clean-target-axis-imbalance
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2384-owner-pivot-primal-dual-obstruction-and-two-probe-repair-corner
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
  - THM-2397-clean-root-same-parent-charged-role-partition
related:
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
  - THM-2401-common-filter-endpoint-or-first-death-certificate
script: 04-computation/lrc14_clean_target_axis_imbalance_thm2403.py
output: 05-knowledge/results/lrc14_clean_target_axis_imbalance_thm2403.out
script_sha256: 083b51185d5b04c9c5561b37158248c1efe575c2d96045e47402842da407d376
output_sha256: 39970491339e89ac95f931e14ab8e05d6c5ba6b16e842d2b9ffd102872cfb150
hash_basis: working-tree bytes (LF)
---

# THM-2403 -- a clean toothpick fires a global fully-all-safe target current

**REPAIRED PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT REAUDIT
PENDING.**

THM-2397 identifies the exclusive `q_*` root word as literal
single-factor deletion support at one common bare endpoint. THM-2400
then isolates the lawful obstruction: a true owner-pivot target
covector moves different factors with unequal slopes.

That unequal-slope branch is positive. Restore the `q_*`-safe factor,
move one other ordinary unit-safe factor together with its target
blocker, and retain every other present safe factor. The scalar cover
makes the untwisted packet identically zero, while a clean toothpick
forces positive shifted mass:

```text
global fully-all-safe present packet
  + one lawful owner-pivot target dipole
  + positive clean-toothpick mass
  -> anchored nonflat rational C_13 mass bank
  -> all twelve nonzero first-target colours
  -> for each colour, some exact nonzero-target fixed triangle
     with a 91-unit deep multiplier.                           (1)
```

The clean set is used only as a positive lower-bound region inside the
global packet. It is never inserted as a frozen projector.

## 1. Two exact banks on one clean toothpick

Put `G=F_13`. Let

```text
Q_*,Q_i,Q_2,Q_3,Q_4 subset G,          |Q_j|=2,

H subset G,                             |H|=4,       (2)
```

be the five ordinary and one guard danger masks on a clean inverse-root
fibre. They cover all thirteen roots. Their total incidence is

```text
4+5*2=14,                                             (3)
```

so exactly one root is double-covered and the other twelve are
single-covered.

Choose `Q_i` distinct from the distinguished ordinary label `Q_*`.
Remove these two labels temporarily and put

```text
B=G minus (H union Q_2 union Q_3 union Q_4).          (4)
```

Thus `B` is the set missed by the other four masks, and

```text
B subset Q_* union Q_i.
```

The literal `q_*`-deletion support at the base shift is

```text
A_0=B minus Q_i.
```

The unique double root gives

```text
(|A_0|,|B|) in {(1,3),(2,3),(2,4)}.                  (5)
```

This is the three-type deletion lemma used in the first candidate
version of this theorem.

Now restore the `q_*`-safe factor and put

```text
C=B minus Q_*.
```

Then

```text
empty != C subset Q_i,             |C| in {1,2}.    (6)
```

Indeed, at most one `Q_i` root can meet one of the other four masks. If
the second met `Q_*`, that would create a second double root. If neither
meets the other four, at most one can meet `Q_*`, again because the
double root is unique. Hence at least one `Q_i` root lies in `C`.

Let `alpha in F_13^*` and translate the retained ordinary factor:

```text
F_s=C minus (Q_i+alpha s),             s in F_13.    (7)
```

Because `C subset Q_i`,

```text
F_0=empty.                                            (8)
```

For each fixed root in `C`, exactly two shifts put it in
`Q_i+alpha s`. Therefore

```text
sum_s |F_s|=11|C|.                                   (9)
```

## 2. One ordinary blocker gate leaves a sharp gap nine

Let

```text
g_s in {0,1},       g_0=1,       #{s:g_s=0}<=2.     (10)
```

Since `|F_s|<=|C|`, equation (9) gives

```text
sum_s g_s|F_s|
 >=11|C|-2|C|
 =9|C|
 >=9.                                               (11)
```

The universal constant nine is sharp in the labelled physical
root/gate universe. Take `alpha=1` and the clean cover

```text
Q_*={0,1},             Q_i={2,3},

Q_2={2,4},             Q_3={5,6},       Q_4={7,8},

H={9,10,11,12}.                                      (12)
```

Then `B={0,1,3}` and `C={3}`. If `g_2=g_3=0` and all
other `g_s=1`, the fully masked mass vector is

```text
(g_s|F_s|)_(s=0)^12
 =(0,0,0,0,1,1,1,1,1,1,1,1,1),                    (13)
```

whose total is nine. The two killed shifts are supplied by one strict
ordinary danger comb: at phase `z=5/26`,

```text
||z-s/13||<1/14

iff s in {2,3}.                                      (14)
```

For comparison, before restoring `Q_*`, the bank

```text
A_s=B minus (Q_i+alpha s)
```

satisfies

```text
sum_s g_s|A_s|-13|A_0|>=1,                          (15)
```

with sharp typewise floors `14,1,10`. Equation (15) remains a valid
finite lemma, but it is not used to localize the analytic packet:
freezing the clean cell would also freeze a moving blocker gate.

## 3. The global lawful owner-pivot packet

Return to one of the `165` scalar rows. Choose the positive shallow
source label `j` from THM-2349, let `a` be the other low blocker, and
put

```text
c=c_3.
```

In the THM-2309 owner pivot choose

```text
u_0=q_*,
```

and choose distinct ordinary graft helpers

```text
k_a=q_i,                 k_c!=q_i,q_*.
```

The true target covectors are the THM-2384 dipoles

```text
eta=e_a-e_(q_i),          ell=e_c-e_(k_c)            (16)
```

modulo thirteen. This is a dual statement. THM-2309's displayed primal
quotient alone would not justify (16).

Let `E_(s,t)` be the actual present packet which retains

```text
all six unit/guard safe factors
and all three blocker safe factors,                  (17)
```

with every factor shifted according to the lawful target action
`s eta+t ell`. Let `Q(Rx)` be any fixed nonnegative rational transported
terminal word, written in its nine-factor THM-2334 form, with

```text
R=13^k,                 k>=1.
```

Because `13|R`, the word is target-neutral. At the untwisted slice the
scalar cover gives the global pointwise identity

```text
E_(0,0)(x)=0                    almost everywhere.  (18)
```

Put

```text
Delta_r(x)=d(c x-r/13),

H_+(r,s,t)
 =integral_T E_(s,t)(x)Q(Rx)Delta_r(x) dx.          (19)
```

This is exactly the THM-2365 lawful table, with a fully-all-safe
present packet rather than the positive shallow-owner rectangle. It is
nonnegative. Since the `c`-safe factor moves with `t`,

```text
H_+(t,s,t)=0                         for all s,t     (20)
```

away from the null strict endpoints.

Define

```text
M_s=sum_r H_+(r,s,0).                               (21)
```

Equation (18) gives the anchored zero

```text
M_0=0.                                              (22)
```

## 4. Clean mass forces the shifted global packet to fire

Let `S` be THM-2392's clean set and

```text
delta=mu(S)>0.
```

For `y in S`, write

```text
x_r=(y+r)/13.
```

At `t=0`, the dipole `eta` moves only `q_i` and blocker `a`.
All other present factors are unchanged. The unit-root support after
restoring `q_*` is precisely the set `F_s` from (7).

The moving blocker-safe factor is independent of the root:

```text
g_s(y)
 =1-d(c_a x_r-s/13)
 =1-d(C_a y-s/13).                                  (23)
```

Since `y` lies outside all three quotient blockers,

```text
g_0(y)=1.
```

The exact root count, almost everywhere away from strict endpoints, is

```text
sum_s d(C_a y-s/13)=2-d(c_a y)<=2.                  (24)
```

Thus (10)--(11) apply pointwise on `S`.

The delayed word is root-constant:

```text
Q(Rx_r)=Q(13^(k-1)y).                               (25)
```

The deep successor count is also exact on `S`:

```text
sum_r Delta_r(x)=2-d(13c x),

d(13c x_r)=d(cy)=0.                                 (26)
```

Put

```text
rho_R
 =integral_S Q(13^(k-1)y)dy.                        (27)
```

Disintegrating `dx` over the thirteen inverse sheets and using
(11), (25), and (26) gives

```text
sum_s M_s
 >=(1/13) integral_S
      2 Q(13^(k-1)y) sum_s g_s(y)|F_s(y)| dy

 >=18 rho_R/13.                                     (28)
```

Every contribution from outside `S` is nonnegative, so it can only
strengthen (28). No indicator of `S` occurs in (19).

For the unfiltered word `Q=1`,

```text
rho_R=delta.
```

For any fixed positive rational THM-2305 terminal word, THM-2397's
two-BV estimate gives

```text
rho_R>0
```

at every sufficiently large clock. Thus (28) applies to an actual
fixed delayed terminal word, not only to the unfiltered control.

The corresponding endpoint packet without the deep successor probe
already obeys

```text
sum_s mu(E_(s,0)Q(Rx))>=9rho_R/13.                  (29)
```

## 5. All twelve first-target colours and exact floors

Every `M_s` is rational. With `zeta=exp(2*pi*i/13)`, put

```text
Mhat(b)=(1/13)sum_s M_s zeta^(b s).                 (30)
```

If `Mhat(b)=0` for some `b!=0`, the rational polynomial

```text
P(X)=sum_s M_s X^s
```

is divisible by `Phi_13(X)=1+X+...+X^12`. Since
`deg P<=12`, every `M_s` would be equal. Equations (22) and (28)
exclude this whenever `rho_R>0`. Hence

```text
Mhat(b)!=0                 for every b!=0.          (31)
```

Let

```text
Gamma=sum_s M_s-13M_0=sum_s M_s.
```

The sharp anchored variance and maximum-mode estimates are

```text
sum_(b!=0)|Mhat(b)|^2
 >=Gamma^2/2028
 >=27rho_R^2/28561,

max_(b!=0)|Mhat(b)|
 >=Gamma/156
 >=3rho_R/338.                                      (32)
```

For the endpoint bank in (29), the corresponding floors are

```text
energy>=27rho_R^2/114244,

max mode>=3rho_R/676.                               (33)
```

THM-2396 gives

```text
delta>=1/26754
```

uniformly and `delta>=66/4459` in the common-core chain. For `Q=1`,
the deep-probe floors in (32) become

```text
universal:
  energy>=3/2271477008164,
  max mode>=1/3014284;

common core:
  energy>=117612/567869252041,
  max mode>=99/753571.                              (34)
```

The endpoint floors in (33) become

```text
universal:
  energy>=3/9085908032656,
  max mode>=1/6028568;

common core:
  energy>=29403/567869252041,
  max mode>=99/1507142.                             (35)
```

No positive lower bound is claimed for each individual exact triangle.

## 6. Every first-target colour yields a 91-unit deep triangle

Use the normalized full transform

```text
B_+(a,b,h)
 =1/13^3 sum_(r,s,t)
    H_+(r,s,t)zeta^(a r+b s+h t).                  (36)
```

At the slice `t=0`, put

```text
J(a,b)
 =1/13^2 sum_(r,s)
    H_+(r,s,0)zeta^(a r+b s)
 =sum_h B_+(a,b,h).                                 (37)
```

Equations (21) and (30) give

```text
J(0,b)=Mhat(b)/13.                                  (38)
```

Meanwhile the diagonal zero (20) at `r=t=0` gives

```text
sum_a J(a,b)
 =1/13 sum_s H_+(0,s,0)zeta^(b s)
 =0.                                                (39)
```

Fix `b!=0`. By (31) and (38), `J(0,b)!=0`; by (39), some
`a!=0` has

```text
J(a,b)!=0.
```

Equation (37) then gives some `h` with

```text
B_+(a,b,h)!=0.                                      (40)
```

THM-2365's exact target typing is

```text
q=(b,a+h).                                          (41)
```

Thus `q!=0` because its first coordinate is the prescribed
`b!=0`.

The absolutely convergent `m`-then-`X` expansion from THM-2365 gives

```text
B_+(a,b,h)
 =sum_(m=a mod 13) sum_X A_(X,m)(b,a+h).            (42)
```

Hence (40) supplies some exact `m,X` with

```text
A_(X,m)(q)!=0,

m=a mod 13.                                         (43)
```

Since `a!=0`, thirteen does not divide `m`. The centred deep danger
coefficient vanishes at every nonzero multiple of seven, so every live
term in (43) also has `7` not dividing `m`. Therefore

```text
gcd(m,91)=1.                                        (44)
```

Equations (31), (40), and (43)--(44) hold separately for every
`b in F_13^*`.

## 7. The repaired localization and exact scope

The first THM-2403 candidate integrated the finite deletion gap over a
fixed clean set and then called the result a lawfully shifted packet.
That was not correctly typed: the clean-set indicator already contains
the unshifted quotient-blocker safe conditions. Multiplying it by the
shifted blocker gate freezes or duplicates the factor that the target
action is meant to move.

The repair is structural:

```text
wrong:
  insert 1_S as a fixed packet factor, then move c_a;

right:
  use the global target orbit E_(s,t),
  use S only to lower-bound its nonnegative shifted mass.        (45)
```

This repair also reveals the stronger fully-all-safe bank (7)--(11).
There is no owner/status/root-translate partition cost.

The conclusion has exact scope:

- `E_(s,t)` is a lawful physical present-factor packet retaining all
  six unit/guard safes and all three blocker safes;
- `Q(Rx)` is a target-neutral fixed transported terminal word;
- `H_+` supplies a genuine THM-2334 fixed triangle with nonzero target
  and its own `91`-unit deep multiplier;
- the untwisted all-safe packet is identically zero, so this is not the
  preselected positive shallow-owner rectangle;
- the nonzero target residues may cancel after the target address is
  forgotten, exactly as THM-2365's target-line identity requires;
- no all-coordinate `91`-unit relation address, preselected marked
  triangle, owner phase, or canonical positive-owner endpoint is
  identified.

Thus the unequal-slope branch is no longer an analytic hostile, and the
fully-all-safe present packet has an exact target-current realization.
The remaining bridge is alignment with the positive owner/terminal
marked current and its all-coordinate address mask.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 8. Exact companion

The dependency-free exact companion:

- verifies clean-cover examples for all three deletion types;
- exhausts `13,728` reduced deletion configurations and `1,084,512`
  ordinary gates;
- exhausts `23,166` restored fully-masked configurations and
  `1,830,114` gates;
- obtains the sharp full-mask table `11,10,9` and `22,20,18`;
- verifies the physical gap-nine control (12)--(14);
- checks the prime-cyclotomic reduction, anchored variance arithmetic,
  and all fractions in (34)--(35);
- verifies all `2,197` first-target phase-projector cases; and
- checks the finite transform normalizations in (37)--(41).

Run

```bash
python3 04-computation/lrc14_clean_target_axis_imbalance_thm2403.py
python3 -O 04-computation/lrc14_clean_target_axis_imbalance_thm2403.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_clean_target_axis_imbalance_thm2403.out
```

Every truth-bearing finite check raises explicitly, so optimized mode
executes the same audit.
