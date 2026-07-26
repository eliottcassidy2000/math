---
id: THM-2397
title: "Clean-root same-parent charged role partition"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. On a
  fixed THM-2392 clean-root owner/status/translate cell, the complement
  of the distinguished q_* singleton/adjacent mask has a measurable
  disjoint partition among the guard and four lower q roles. For every
  nonzero target-root colour, some role has strictly negative real
  charged product with q_*; one fixed role works on at least four
  colours. Summed over roles and colours, one fibrewise mixed
  coefficient is at most -rho/845 and one aggregate same-parent
  charged product is at most -rho^2/845. THM-2396 makes these uniform
  on the last septimal lane and sharpens the common-core floors to
  33/97964230 and 1089/11357385040820. On the owner-resolved cell the
  same product carries every septimal colour. This is a derived
  least-role Boolean current at the clean parent, not a translated
  single-factor probe, terminal current, row exclusion, ledger
  decrement, or proof of LRC(14).
source: codex-2026-07-26-clean-root-charged-role
depends_on:
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
related:
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2383-polarized-complete-subcube-gram-tomography
script: 04-computation/lrc14_clean_root_charged_role_thm2397.py
output: 05-knowledge/results/lrc14_clean_root_charged_role_thm2397.out
script_sha256: b7e0f12ca9848105d7eaa3ca4c3ddff48b3b85f30df8f4d67d6d5610201bc25b
output_sha256: 55f3af504e55af3594b5e55c8ff3774fafb9c789ddbf058ee32995b77a89f529
hash_basis: working-tree bytes (LF)
---

# THM-2397 -- a clean root has a charged lower-role neighbour

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2392 gives a positive parent set carrying a fixed singleton or
adjacent two-root `q_*` word. THM-2396 now makes that parent mass
uniform. The missing target-overlap question at this stage has a
pointwise answer: every root outside the exclusive `q_*` mask is
covered by at least one of the other five scalar roles.

A deterministic least-role assignment turns those five covers into a
partition. Its Fourier sum is the complement word, hence the exact
negative of the `q_*` word at every nonzero colour. No quadratic
tomography or external pair twist is needed to prove that this
same-parent charged product is nonzero.

## 1. Fixed clean-root cell

Retain THM-2392's clean-parent set `S`. For `y in S`, write

```text
x_r=(y+r)/13,                         r in F_13.       (1)
```

Let `q_*` be the distinguished top septimal ordinary label. Partition
`S` by:

```text
one deterministic low-blocker owner:      2 choices;

q_* singleton/adjacent status:             2 choices;

root translate:                           13 choices.  (2)
```

On one fixed cell `Y`,

```text
rho=mu(Y)>=delta/52.                                   (3)
```

The exclusive `q_*` root mask `A` is constant on `Y`. It is either one
singleton or two adjacent roots. Put

```text
W=1-A.                                                 (4)
```

With the normalized root-target transform

```text
a_k=(1/13)sum_r A(r) zeta^(-kr),

w_k=(1/13)sum_r W(r) zeta^(-kr),          zeta^13=1,  (5)
```

THM-2392 proves, for every `k!=0`,

```text
w_k=-a_k!=0.                                          (6)
```

The exact nonzero-colour energies are

```text
singleton:  sum_(k!=0)|a_k|^2=12/169;

adjacent:   sum_(k!=0)|a_k|^2=22/169.                 (7)
```

## 2. Least-role partition without fixing the full root word

Order the five roles

```text
L={H} union {q_i:q_i!=q_*}.                           (8)
```

For each `y in Y` and root `r`, let `chi_i(y,r)` be the activity bit of
role `i in L`. Define

```text
B_i(y,r)
 =W(r) chi_i(y,r) product_(j<i)(1-chi_j(y,r)).        (9)
```

These are measurable Boolean masks. More importantly,

```text
sum_(i in L) B_i(y,r)=W(r)                            (10)
```

pointwise. Indeed, if `r` is not a `q_*` root, the cover gives at least
one active role in `L`. If `r` is the one `q_*` double root excluded
from `A`, the second member of the double lies in `L`. Thus every
`W`-root has a least active role. This argument does not partition by
the `648648000` possible complete labelled root words.

Write

```text
b_(i,k)(y)
 =(1/13)sum_r B_i(y,r)zeta^(-kr).                    (11)
```

Equations (6) and (10) give, for every `y in Y` and `k!=0`,

```text
sum_i b_(i,k)(y)=-a_k.                               (12)
```

## 3. Fibrewise signed grouping

Define the coefficient-derived same-parent mixed grouping

```text
G_i(k)
 =integral_Y a_k conjugate(b_(i,k)(y))dy.            (13)
```

Because `a_k` is fixed on `Y`, equation (12) gives

```text
sum_i G_i(k)=-rho |a_k|^2<0                          (14)
```

for every `k!=0`. Therefore every nonzero root-target colour has a
strictly negatively correlated lower role:

```text
for every k!=0, some i has

Re G_i(k)<=-rho |a_k|^2/5<0.                         (15)
```

All masks are real, so the products at `k` and `-k` are complex
conjugates and have the same real part. Choose one witness role for
each of the six conjugate colour pairs. By pigeonhole among five roles,
one fixed role works for at least two pairs, hence at least four
distinct nonzero colours.

Summing (14) over all five roles and twelve colours, then using (7),
gives one pair `(i,k)` with

```text
-Re G_i(k)
 >=rho/845                                          (16)
```

uniformly in the singleton/adjacent alternatives. In the adjacent
case alone, the stronger summed floor is

```text
-Re G_i(k)>=11rho/5070.                              (17)
```

The per-colour adjacent floor is also explicit:

```text
-Re G_i(k)
 >=4rho sin^2(pi/26)/(5*169)                         (18)
```

for some role `i`.

## 4. Aggregate currents have genuine common target support

The grouping in (13) is linear in parent mass. There is also an exact
aggregate charged product. Put

```text
D(k)=integral_Y a_k dy=rho a_k,

R_i(k)=integral_Y b_(i,k)(y)dy.                      (19)
```

Then

```text
sum_i R_i(k)=-rho a_k,

 D(k)conjugate(R_i(k))=rho G_i(k),

sum_i D(k)conjugate(R_i(k))
 =-rho^2 |a_k|^2.                                   (20)
```

Thus the same conclusions as (15) hold for the aggregate currents:
for every `k!=0`, `D(k)` and at least one `R_i(k)` are simultaneously
nonzero, with negative real charged product. One fixed role overlaps
`D` on at least four colours. Moreover the same pair `(i,k)` which
realizes either fibrewise bound realizes its aggregate version after
multiplication by `rho`. Summing over all pairs gives

```text
some (i,k):

-Re[D(k)conjugate(R_i(k))]>=rho^2/845.               (21)
```

This is the exact same-target support overlap tested abstractly by
THM-2380, on the clean-root `F_13` line and in one common parent gauge.
It proves the overlap algebraically; it does not manufacture
THM-2380's separate physical translated-correlation measurement.

## 5. Uniform last-lane and common-core constants

The quantitative localization in THM-2396 and the complementary
branches of THM-2393 give, throughout the last septimal lane,

```text
delta>=1/26754,

rho>=1/1391208.                                      (22)
```

Consequently (16) and (21) give

```text
-Re G_i(k)>=1/1175570760,

-Re[D(k)conjugate(R_i(k))]
 >=1/1635463445878080                               (23)
```

for some lower role and nonzero colour.

On the literal common-core chain, THM-2396 gives

```text
delta>=66/4459,

rho>=33/115934.                                     (24)
```

The corresponding sharper floors are

```text
-Re G_i(k)>=33/97964230,

-Re[D(k)conjugate(R_i(k))]
 >=1089/11357385040820.                             (25)
```

All fractions are exact and reduced.

## 6. Owner-resolved septimal tensor

Instead use THM-2392's owner-resolved cell `Y_(7x13)`. It fixes the
septimal excess address `d`, the root status, and the root translate,
and in the common-core branch has mass

```text
rho_o>=delta/338>=33/753571.                         (26)
```

For every `ell in F_7`, its fixed septimal coefficient is

```text
j_ell=(1/7)zeta_7^(-ell d).                          (27)
```

Repeating Sections 2--4 on this cell gives one `(i,k)` with

```text
-Re G_i(k)>=33/636767495,

-Re[D(k)conjugate(R_i(k))]
 >=1089/479849517974645.                            (28)
```

Multiplication by `j_ell` yields a nonzero derived aggregate
`F_7 x F_13` coefficient for every `ell`, of magnitude at least

```text
1089/3358946625822515.                              (29)
```

The phase in (29) is fixed by the owner-resolved septimal address; it
does not depend on a choice of full labelled thirteen-root word.

## 7. Scope and next target

The theorem proves

```text
uniform positive clean-root mass
  -> same-parent q_*/lower-role charged overlap
  -> one role live on at least four target colours.              (30)
```

The masks `B_i` are lawful positive Boolean subsets, labelled by actual
guard/lower-`q` roles. They are nevertheless **derived least-role
selectors**, not translations of one single danger factor. The theorem
does not identify `R_i` with a terminal blocker-word current, preserve
the endpoint phase through later factors, or make (29) a separately
measurable physical pair-twist observable.

Thus THM-2380's disjoint-support hostile is impossible at this
clean-root stage, but terminal filtration can still destroy the
alignment. The next sharp target is:

```text
transport one B_i-labelled charged edge through the terminal word
without changing its parent gauge or root-target colour.                 (31)
```

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 8. Exact companion

The dependency-free exact companion:

- exhausts all `63` nonempty local activity subsets and verifies that
  every nonexclusive-`q_*` root has a least role among the other five;
- checks the singleton/adjacent nonzero-colour energy ledger;
- verifies the `5*12=60` pair pigeonhole and the fixed-role
  four-colour conjugate-pair pigeonhole;
- checks every rational floor in (22)--(29); and
- retains separate fibrewise, aggregate, and septimal normalizations.

Run

```bash
python3 04-computation/lrc14_clean_root_charged_role_thm2397.py
python3 -O 04-computation/lrc14_clean_root_charged_role_thm2397.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_clean_root_charged_role_thm2397.out
```

after LF normalization. Every assertion remains active under optimized
Python.
