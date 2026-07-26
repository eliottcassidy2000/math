---
id: THM-2395
title: "Common-core successor shell and forced hole-escape tax"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. In
  THM-2394's no-clean common-core residual, the mixed-septimal successor
  shells of q_*/7 and c_3/7 intersect in mass at least 5/637. Hence
  type-II base mass at least 1216/8281 remains high-safe after
  multiplication by thirteen. A covariant type-II hole must land in
  type III, whose total capacity is at most 24/169=1176/8281.
  Therefore a base set of mass at least 40/8281 undergoes a genuine
  hole-status escape. Its selected physical roots have mass at least
  40/57967 and satisfy R(x)=1, R(13x)=0,
  c_1(x)=C_1(13x)=1, and c_2(x)=C_2(13x)=0. The actual c_1-factor
  translation has successor overlap zero, so THM-2362 gives nonzero
  mode sum at least 440/753571, one real mode at least 110/2260713,
  and energy at least 48400/1703607756123. A fixed current septimal
  pair gives an F_7^* x F_13^* role-jet tensor. A strict physical
  type-II/type-III two-cycle shows that local carry-labelled covariance
  remains possible. This is a positive successor-role jet, not a row
  exclusion, ledger decrement, or proof of LRC(14).
source: codex-2026-07-26-common-core-successor-shell
depends_on:
  - THM-2362-thirteen-shift-successor-statistic-and-role-jet-floor
  - THM-2394-common-core-six-address-transversal-and-labelled-hole-trichotomy
related:
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
script: 04-computation/lrc14_common_core_successor_escape_thm2395.py
output: 05-knowledge/results/lrc14_common_core_successor_escape_thm2395.out
script_sha256: 219d8aed4f27b7cf4bd1b55a4125f793faa73ec0eea9dcfe3ce734fb9aa7c87e
output_sha256: d568693c4e96e381fcba501edb0aa8c5ee639db8457d6e0cbc26a581ba7eedf6
hash_basis: working-tree bytes (LF)
---

# THM-2395 -- a forced successor-role escape

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2394 turns the last no-clean common-core boundary into three labelled
seven-root states. This theorem follows those states through the actual
map

```text
T(y)=13y mod 1.                                      (1)
```

The central fact is a capacity mismatch:

```text
retained type II:       at least 1216/8281,

all type III:           at most 1176/8281.           (2)
```

Hole covariance sends II only to III. The difference in (2) is therefore
a positive, canonically typed failure of covariance rather than an
unlocated mass surplus.

## 1. The corrected high-safe base

Use the divided base variables from MISTAKE-264's repair of THM-2393.
Put

```text
u=q_*/7,

r=C_3/7,

v=c_3/7=13r.                                        (3)
```

The last-lane valuations give

```text
nu_7(u)=0,                  nu_7(v)>0.               (4)
```

The high-safe base is

```text
Y_0=D_u^c intersection D_r^c intersection D_v^c.    (5)
```

Let `E_II,E_III` be the type-II and type-III base sets in the THM-2394
trichotomy. That theorem proves

```text
mu(E_II)>=3335/8281,

mu(E_III)<=24/169=1176/8281.                         (6)
```

## 2. A mixed-septimal shell lemma

For any positive speed `w`, define its successor shell

```text
S_w=D_w^c intersection D_(13w).                     (7)
```

The common-line overlap is `1/91`, so

```text
mu(S_w)=1/7-1/91=12/91.                             (8)
```

> **Mixed-shell lemma.** If `7` does not divide `u` and `7` divides
> `v`, then
>
> ```text
> mu(S_u intersection S_v)>=5/637.                  (9)
> ```

**Proof.** Disintegrate over the seven roots

```text
y_s=(z+s)/7.
```

Because `7|v`, the `S_v` bit is constant on the roots and equals the
`S_(v/7)` bit at `z`. Its base mass is `12/91`.

Both `D_u` and `D_(13u)` occupy one generic root. Thus `S_u` occupies
one root unless those two addresses coincide, in which case it occupies
none. Root disintegration and

```text
mu_z{the two addresses coincide}
 =7mu(D_u intersection D_(13u))
 =1/13                                                (10)
```

give

```text
7mu(S_u intersection S_v)
 >=12/91-1/13
 =5/91.
```

Division by seven proves (9). QED.

Equations (8)--(9) imply

```text
mu(S_u union S_v)
 <=24/91-5/637
 =163/637
 =2119/8281.                                        (11)
```

## 3. At least `1216/8281` of type II survives the shift

For `y in Y_0`, the middle high condition is automatic after shifting:

```text
1_(D_r)(Ty)=1_(D_(13r))(y)=1_(D_v)(y)=0.            (12)
```

The other two conditions can fail only through

```text
D_u(Ty)=D_(13u)(y),

D_v(Ty)=D_(13v)(y).
```

Since `D_u(y)=D_v(y)=0` on `Y_0`,

```text
Y_0 minus T^(-1)(Y_0)
 subset S_u union S_v.                              (13)
```

Consequently

```text
mu(E_II intersection T^(-1)(Y_0))
 >=3335/8281-2119/8281
 =1216/8281.                                        (14)
```

No independence of the successor shells from `E_II` is used.

## 4. The covariant-hole automaton

At a base `y`, write the core addresses as

```text
(a,b,c)
 =(r_h(y),r_(13h)(y),r_(169h)(y)),
```

and let `d` be the unique `K`-hole. Put

```text
kappa=floor(13y) mod 7.
```

If the physical hole root `(y+d)/7` is sent by `T` to address `d'`, then

```text
d'=kappa-d.                                         (15)
```

THM-2394's carry law gives the first two successor core addresses

```text
a'=kappa-b,

b'=kappa-c.                                         (16)
```

Assume for the moment that the image address `d'` is again the successor
`K`-hole. Direct substitution in the three cases of THM-2394 gives the
exact automaton

```text
I   -> I or II,

II  -> III,

III -> I or II.                                     (17)
```

For example, in type II one has `d=b!=c`, so

```text
d'=a'!=b'.
```

The only successor type whose hole is `a'` rather than `b'` is type III.
The other two lines of (17) are identical one-line checks.

Let `E_cov` be the subset of (14) on which the hole is covariant. By
(17),

```text
T(E_cov) subset E_III.
```

For every measurable set `A`,

```text
mu(A)<=mu(T(A)),
```

because `A subset T^(-1)(T(A))` and `T` preserves Haar measure. Hence

```text
mu(E_cov)<=mu(E_III)<=1176/8281.                    (18)
```

Subtracting (18) from (14) forces a noncovariant type-II carrier

```text
mu(E_escape)>=40/8281.                              (19)
```

## 5. Physical residual descent and owner transport

For `y in E_escape`, select the current type-II hole root

```text
x=(y+b(y))/7.                                       (20)
```

Root disintegration gives a physical set `F_escape` with

```text
mu(F_escape)>=40/(7*8281)=40/57967.                 (21)
```

The image root has successor address `a'`. Since the transition is not
covariant, the successor type is I or II, whose unique hole is `b'`.
Equation (16) and current `b!=c` give `a'!=b'`. The successor
six-address transversal therefore has value one at the image root.

Writing the full scalar residual as

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),
```

both fibres have `q_*` safe, so `R=1-K`. Thus on `F_escape`,

```text
R(x)=1,                       R(Tx)=0.               (22)
```

The common-core labels are equally explicit:

```text
1_(D_(c_1))(x)=1,             1_(D_(c_2))(x)=0,

1_(D_(C_1))(Tx)=1,            1_(D_(C_2))(Tx)=0.    (23)
```

In words, the positive residual loses its hole status while the middle
owner transports from actual `c_1` to quotient `C_1`.

At `Tx`, exactly one of the guard and four lower ordinary roles is
active. Partitioning by that label gives one fixed successor role on a
physical set of mass at least

```text
(40/57967)/5=8/57967.                                (24)
```

## 6. The actual `c_1` factor has a zero-drift role jet

The carrier in (21) lies in the pure word

```text
{R=1} intersection D_(c_1)
 minus (D_(c_2) union D_(c_3)).                     (25)
```

On `{R=1}` with the other two blockers safe, the scalar cover makes the
displayed `c_1` factor logically redundant. Therefore THM-2362 applies
to translation of the **actual named factor**, not an auxiliary
duplicate probe.

Push the remaining carrier forward in the `c_1` coordinate. Since

```text
c_2=13c_1
```

and `c_2` is safe on (21), THM-2362's successor overlap is exactly

```text
rho_+=0.                                             (26)
```

With `rho=mu(F_escape)>=40/57967`, its formulas give

```text
sum_(k!=0) Mtilde(k)
 >=11rho/13
 >=440/753571,                                      (27)

some k!=0 has
 Re Mtilde(k)>=11rho/156
             >=110/2260713,                         (28)

sum_(k!=0)|Mtilde(k)|^2
 >=121rho^2/2028
 >=48400/1703607756123.                             (29)
```

Thus (21) is a canonical zero-drift positive role jet.

There is also a fixed-phase transverse refinement. Partition
`E_escape` by its current ordered type-II pair `(b,c)`. One of the
`42` cells has base mass at least

```text
20/173901
```

and selected physical mass at least

```text
20/1217307.                                         (30)
```

On that cell THM-2394's septimal hole/double phases are constant.
Applying (28) to (30) and multiplying by the normalized septimal
coefficient or charged product gives, for every `ell!=0`, a common
nonzero factor-translation mode `k!=0` with

```text
|F_7^* x F_13^* role-jet tensor|
 >=55/332324811,                                    (31)

|charged role-jet tensor|
 >=55/2326273677.                                   (32)
```

These are factor-coloured coefficients, not terminal target currents.

## 7. A strict exclusive-owner two-cycle

The covariant branch of (17) is physically real. Take

```text
y_0=1/24,                    y_1=13/24,

T(y_0)=y_1,                  T(y_1)=y_0,

h=1,

H=29,

(q_1,q_2,q_3,q_4)=(27,55,83,71),

q_*=14,

C_3=2*7^2*13^4=2798978,

c_3=13C_3=36386714.                                (33)
```

At both bases, the guard roots are `{5,6}`, the four lower ordinary
roots are `{1},{2},{3},{4}`, and all three high words are safe. Thus

```text
K=(0,1,1,1,1,1,1).                                  (34)
```

The core address triples are

```text
y_0: (a,b,c)=(0,1,0),       type III, hole 0,

y_1: (a,b,c)=(6,0,6),       type II,  hole 0.       (35)
```

The selected physical holes `1/168` and `13/168` form a two-cycle under
`T`. Every inequality is strict, so this is an open local chamber, not
an endpoint artifact. It is a lawful physical comb hostile to any
pointwise claim that exclusive owner type must persist. It is not a
global scalar cover or an LRC counterexample.

Equation (17) also shows that every one-state covariant loop is type I,
where the two actual low blockers share the hole. Hence (35) is the
minimal exclusive-owner covariant cycle.

## 8. Scope

THM-2395 proves

```text
no-clean common core
 -> positive owner-labelled residual descent
 -> nonzero actual-factor successor jet.            (36)
```

It does not identify the `F_13` factor-translation colour with the
canonical relation target or preserve the terminal endpoint phase.
Those are exactly the losses recorded in THM-2362. The two-cycle (33)
also prevents a pointwise closure from the carry automaton alone.

The next target is either a target/terminal intertwiner for (31), or a
whole-`49`-orbit incompatibility using both the top-danger bin and all
six top-safe bins. No row is excluded, the scalar ledger remains `165`,
and LRC(14) remains open.

## 9. Exact companion

The dependency-free exact companion:

- verifies `mu(S_w)=12/91`, the mixed shell floor `5/637`, and two
  direct physical shell controls;
- checks every rational constant in (11), (14), and (18)--(32);
- exhausts the three-state carry-labelled automaton (17);
- verifies all physical root words, valuations, owner types, and
  transitions in the strict two-cycle (33)--(35);
- checks that every displayed hostile predicate is away from its strict
  endpoint; and
- runs every assertion under ordinary and optimized Python.

Run

```bash
python3 04-computation/lrc14_common_core_successor_escape_thm2395.py
python3 -O 04-computation/lrc14_common_core_successor_escape_thm2395.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_common_core_successor_escape_thm2395.out
```

after LF normalization.
