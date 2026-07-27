---
id: THM-2559
title: "Target-informed chord and universal old repair packet"
status: >
  PROVED + VERIFIED-EXACT.  On every one of the 165 live first-depth-one
  owner rows, the already known target-active guard/unit role k_a can be
  inserted into the root selector itself.  Its nonempty failure mask is
  disjoint from the occupied A_0 mask.  Lexicographically selecting one
  root from each and stratifying by their nonzero difference produces a
  positive fixed-slope occupied-to-k_a-failure chord.  The Cayley endpoint
  sum is pointwise minus the owner weight.  At the selected head, every
  blocker status is inherited from the occupied source; in particular the
  deepest blocker and the blocker paired with k_a remain safe.  Applying
  THM-2379 therefore gives a positive primitive deep-by-repair coefficient
  on every row, with exact relative floor 11/24336 for an ordinary k_a and
  1/2704 for a guard k_a.  The selected head is at the old predecessor
  horizon, the chord slope need not be the physical guard slope, and the
  paired blocker is retained but not co-shifted.  No genuinely later root,
  lawful paired target dipole, Hall arrival, scalar-row exclusion, or
  LRC(14) proof follows.
source: lrc-semantic-frontier-2026-07-27-target-informed-chord
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2531-prime-necklace-guard-boundary-selector
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2558-sparse-owner-fibre-all-slope-target-role-forcing-boundary
related:
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
script: 04-computation/lrc14_target_informed_chord_old_repair_thm2559.py
output: 05-knowledge/results/lrc14_target_informed_chord_old_repair_thm2559.out
script_sha256: 39a9b6121b619433a1c000ef7131dcc57284bbc009c1ba72be47ec77a1c75dcd
output_sha256: 23124d9a7f34f48f8ce8202a41134c3d5ae55e73f6aa926d7779db39df927b71
hash_basis: working-tree bytes (LF)
---

# THM-2559 -- use the known target role to select its own old head

**PROVED + VERIFIED-EXACT.**

THM-2558 shows that a selector which sees only the owner mask can miss a
target-active failure on 202 rotation necklaces.  That is the wrong
information order once THM-2350/2461 has already named the unique
target-active member `k_a` of the five-role bank.  Retain its literal failure
mask as a sidecar and select a root from that mask directly.

This simple reversal has two exact consequences:

```text
occupied owner mask + known disjoint k_a-failure mask
  -> positive fixed-slope occupied-to-failure chord
  -> positive deepest-comb by k_a-complement repair coefficient.       (1)
```

The construction is uniform over all `165` valuation rows and bypasses the
all-slope blind layer completely.  It is deliberately an **old-predecessor**
theorem.  The final paired target action and future-root map remain absent.

## 1. The active direct fibre and its exact masks

Fix one live row and use THM-2349/2537 notation

```text
F(x)=1_(E_j)(x)1_(Q_(j,sigma))(13^k x),       k>=2,

W(x)=F(x)g(x),                                g in {0,1},       (2)
```

where `g` is the sufficiently late rational Boolean owner block.  On the
direct predecessor fibre

```text
iota_r(z)=(z+r)/13,                            r in F_13,       (3)
```

put

```text
e_r(z)=F(iota_r(z)).                                           (4)
```

Let `B` be the rational base set on which

```text
g(z)=1,                         e(z) is nonconstant.          (5)
```

THM-2537 equation (49) gives

```text
integral_B g(z) dz>0.                                         (6)
```

All statements below ignore the finite set of strict-comb endpoints.

The source blocker is `j`, the other nondeep blocker is `a`, and the deepest
blocker is `c=c_3`.  Every blocker coefficient has positive 13-adic
valuation.  Hence for `c_h=13^(lambda_h)u_h`, `lambda_h>=1`,

```text
c_h iota_r(z)
 =13^(lambda_h-1)u_h z                 mod 1,                 (7)
```

independent of `r`.  The word in (2) is also root-constant because

```text
13^k iota_r(z)=13^(k-1)z                 mod 1.               (8)
```

On a fibre in `B`, some root belongs to `F`.  Thus all root-constant blocker
and word factors in the definition of `F` equal one there.  The remaining
varying mask is exactly

```text
e_r(z)=1_(A_0)(iota_r(z)),

A_0=C_H minus union_(i=1)^5 D_(q_i).                         (9)
```

In particular `e(z)` is nonempty and is safe for all six guard/unit roles.
This is the direct-mask typing repaired in THM-2558; no Perron-rebased sparse
owner root is used.

Choose the THM-2350 pivot as in THM-2461, with `q_*=k_b`.  The role `k_a` is
the unique target-active member of the five roles other than `k_b`.  Let

```text
L_a=1,       if k_a is an ordinary unit role,
L_a=2,       if k_a is the guard role,                         (10)
```

and define its literal failure mask on (3) by

```text
T_a(z)={r:d_(L_a)(k_a iota_r(z))=1}.                         (11)
```

Here `d_L(y)=1_(||y||<L/14)`.  THM-2379's translated-tooth identity gives

```text
|T_a(z)| in {1,2},             L_a=1,
|T_a(z)| in {3,4},             L_a=2.                         (12)
```

Thus `T_a(z)` is always nonempty and proper.  Since (9) is safe for `k_a`,

```text
T_a(z) subset e(z)^c.                                      (13)
```

This disjointness, not a capacity comparison, is the key input.

## 2. A rational target-informed chord

Fix once and for all a nonzero orientation gauge `tau_0`; the physical guard
gauge from THM-2526 is a canonical choice.  For any nonempty proper mask `m`,
let `alpha_(tau_0)(m)` be THM-2531's unique lexicographic marker.  Define on
`B`

```text
s(z)=alpha_(tau_0)(e(z)),
t(z)=alpha_(tau_0)(1_(T_a(z))),
delta(z)=t(z)-s(z) in F_13^*.                              (14)
```

The first marker is occupied, the second lies in `T_a`, and (13) makes
`delta(z)` nonzero.  Every input factor is a rational step function.  The
two masks are therefore piecewise constant off finitely many rational walls;
the finite marker lookup makes `s,t,delta` rational measurable step maps.

For each fixed `delta!=0`, put

```text
B_delta={z in B:delta(z)=delta}                             (15)
```

and define the source and head packets by

```text
S_delta(iota_r(z))
 =g(z)1_(B_delta)(z)1_(r=s(z)),

T_delta(iota_r(z))
 =g(z)1_(B_delta)(z)1_(r=t(z)).                            (16)
```

They are rational Boolean step functions.  Since `dx=dz/13` on every
inverse branch,

```text
mu(S_delta)=mu(T_delta)
 =1/13 integral_(B_delta)g(z)dz,                            (17)

sum_(delta!=0)mu(S_delta)
 =1/13 integral_B g(z)dz>0.                                (18)
```

Consequently at least one **fixed** nonzero chord slope `delta` has positive
source and head mass.  Moreover

```text
S_delta<=W,                    T_delta W=0,                 (19)
```

and every point of `T_delta` has literal `k_a` failure.  The slope in (18)
is selected after a finite rational stratification.  It need not equal the
physical guard slope `tau_0`; no guard-oriented current is inferred merely
from that equality of notation.

The construction is affine-covariant when the orientation gauge and both
masks are transported together.  It uses the named target-role mask as
input.  It is therefore not a contradiction to THM-2558's target-blind
all-slope hostile.

## 3. The Cayley endpoint sum has a fixed sign

Let

```text
p^W_r(z)=g(z)[e_r(z)-1/13 sum_u e_u(z)]                    (20)
```

be the centred root profile.  Use THM-2537's convention

```text
(P_delta p)_v=p_(v+delta),

C_delta=P_delta-P_delta^2+...+P_delta^11-P_delta^12.       (21)
```

Its signless-incidence identity is

```text
(I+P_delta)C_delta=P_delta-I.                              (22)
```

On `B_delta`, `t=s+delta`, `e_s=1`, and `e_t=0`.  Evaluating (22) at `s`
therefore gives the pointwise identity

```text
(C_delta p^W)_s+(C_delta p^W)_t
 =p^W_t-p^W_s=-g(z).                                      (23)
```

Equivalently, for THM-2533's circle lift `mathcal C_delta`,

```text
-integral_(support(S_delta) union support(T_delta))
       mathcal C_delta W(x)dx

 =1/13 integral_(B_delta)g(z)dz
 =mu(S_delta)=mu(T_delta)>0                                (24)
```

for every positive stratum.  There is no sign cancellation and no
singleton/pair or 202-necklace exception.  The finite referee checks (23)
on all `1,577,940` ordered pairs of disjoint nonempty source and failure
masks.

## 4. The head inherits the exact THM-2379 anchors

Fix a positive slope stratum and write

```text
w=T_delta,                         rho=integral w>0.         (25)
```

At the occupied source root, membership in `E_j` says

```text
j dangerous,                 a safe,                 c_3 safe. (26)
```

All three statuses are root-constant by (7), so (26) remains true at the
selected head even though its guard/unit failure makes `A_0`, and hence
`E_j`, false there.  Combining (11), (16), and (26) gives

```text
support(w) subset {d_1(c_3 x)=0}
                         intersection {d_(L_a)(k_a x)=1},   (27)

support(w) subset {d_1(a x)=0}.                             (28)
```

Equation (27) is exactly THM-2379 hypothesis (3), with

```text
c=c_3,             v=k_a,             alpha=beta=0.        (29)
```

Define its complement-repair table

```text
K^+_(r,q)
 =integral_T w(x)
    d_1(c_3x-r/13)u_(L_a)(k_ax-q/13)dx,                    (30)

Khat^+(A,B)
 =1/13^2 sum_(r,q)K^+_(r,q)zeta^(Ar+Bq).                   (31)
```

THM-2379 gives some `A,B in F_13^*` with

```text
Re Khat^+(A,B)
 >=(13-2L_a)rho/[13^2*12^2]                                (32)

 =11rho/24336,                         L_a=1,

 =rho/2704,                            L_a=2.               (33)
```

The deep and repair probe multipliers selected by this coefficient are
units modulo `91`.  Hence every live row has a positive old-predecessor
deep-by-`k_a` complement-repair coefficient, not merely a generic
first-failure label.

Equation (28) retains the blocker paired with `k_a` as a literal safe
sidecar.  This is stronger typing than bare THM-2379, but it does **not**
yet realize the lawful target dipole.  In (30) the `k_a` complement is
translated while the paired blocker inside `w` is held fixed.  THM-2350's
target action must co-translate those two factors in opposite directions.
That polarized two-factor table is new mathematics and is not supplied by
(30)--(33).

## 5. Exact gain and exact stopping boundary

The theorem removes two former local uncertainties on all `165` rows:

```text
which old empty head carries k_a?       resolved by (14)--(18);

does that head support a positive
deep-by-repair primitive coefficient?   resolved by (27)--(33).           (34)
```

It does not remove the chronological and target-action debts:

```text
same-predecessor chord
  -/-> old-to-new owner transition;

static paired-blocker sidecar
  -/-> oppositely co-shifted target dipole;

old root t(z)
  -/-> THM-2545 genuinely later root b;

positive signed repair coefficient
  -/-> nonnegative Hall atom.                              (35)
```

THM-2555 makes the third separation sharp: an old ancestry root and the
future immediate digit are independent coordinates until a semantic
intertwiner is proved.  Thus (24) does not prove THM-2537 equation (56) in
its later-field sense, and (32) is not a Hall diagonal.  No scalar row is
excluded and LRC(14) remains open.

The next faithful object is the three-factor polarized head table which
translates the deepest probe, the `k_a` complement, and its paired blocker
with the exact THM-2350 dipole signs while retaining (25).  A nonzero target
character there would close the paired-action debt; a circulant table would
exhibit its sharp survivor.

## 6. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_target_informed_chord_old_repair_thm2559.py
python3 -O 04-computation/lrc14_target_informed_chord_old_repair_thm2559.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_target_informed_chord_old_repair_thm2559.out.
```

The referee uses only integers and `Fraction` arithmetic.  It checks:

- all `1,577,940` disjoint nonempty source/failure-mask pairs and their
  nonzero target-informed chord;
- the sign in (23) for every pair and the complete twelve-slope stratum
  histogram;
- blocker root constancy on all `165` valuation profiles and positive-clock
  word constancy;
- THM-2379's exact ordinary/guard translated-tooth counts on a common exact
  chamber refinement;
- the THM-2558 blind mask `e={0,1,4}`, which misses failure root `3` at every
  old all-slope head but is selected by the new chord `0 -> 3`; and
- both exact relative floors in (33).

The rational measurability argument, inheritance of (26), and application of
THM-2379 are symbolic proofs above, not finite extrapolations. **QED.**
