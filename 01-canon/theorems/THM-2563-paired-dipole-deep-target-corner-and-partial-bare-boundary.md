---
id: THM-2563
title: "Paired-dipole deep corner and the missing-left-residue boundary"
status: >
  PROVED + VERIFIED-EXACT.  Let w be any positive old-head packet which
  keeps the deepest blocker c_3 safe and makes the target-active role k_a
  dangerous, as every THM-2559 packet does.  Insert a translated deepest
  danger probe and, on the moving partial-bare endpoint, both complete
  THM-2350 safe dipoles a--k_a and c_3--k_b.  The resulting nonnegative
  13^3 table has total mass at least 63 rho, because the two paired safe
  shift sets have at least 7 and 9 common shifts and at most one role is
  the guard.  Its r=0 and s=0 coordinate planes and its r=t diagonal are
  zero.  Double finite-character summation therefore forces a coefficient
  with nonzero deepest colour A and nonzero a--k_a moving-endpoint colour
  B.  The normalized full-table real-part floor is 7 rho/35152, improving
  to 9 rho/35152 when both graft roles are ordinary.  This closes the
  previously missing paired-blocker co-shift on an old-head partial-bare
  cospan, uniformly on all 165 rows.  Every live displayed deep-probe
  multiplier is coprime to 91.  It does not produce a completed THM-2334
  target coefficient: the
  fixed selector/head leg is not target-co-shifted, so its left residue is
  absent, and the danger-to-safe transition is off-diagonal under THM-2452
  complete-mask restoration.  No terminal word at the head, future root,
  scalar-row exclusion, or LRC(14) conclusion follows.
source: lrc-semantic-frontier-2026-07-27-paired-dipole-corner
depends_on:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2559-target-informed-chord-and-universal-old-repair-packet
related:
  - THM-2561-positive-interval-physical-six-comb-blind-necklace-hostile
script: 04-computation/lrc14_paired_dipole_deep_corner_thm2563.py
output: 05-knowledge/results/lrc14_paired_dipole_deep_corner_thm2563.out
script_sha256: 6f3bde89c134f7e3db430f27e015dfd795d26cf9d52f932768e132850f3b294c
output_sha256: 631d6fa087cf97ae0f7180ed2a2f46bcd3090dc663480e283064b401c1702338
hash_basis: working-tree bytes (LF)
---

# THM-2563 -- the paired repair corner is positive, but still partial-bare

**PROVED + VERIFIED-EXACT.**

THM-2559 puts the unique target-active failure on a positive old head and
retains both paired blockers as safe sidecars.  THM-2461 identifies the
remaining local defect correctly: translating the `k_a` safe complement
alone is not the lawful target dipole because its blocker `a` must move in
the opposite direction.  The same issue applies to the designated graft
`k_b=q_*` and its paired deepest blocker `c=c_3`.

Completing both pairs gives a stronger positive corner:

```text
fixed old head:          c safe, k_a dangerous;

moving partial endpoint:
  a safe  <----> k_a safe,
  c safe  <----> k_b safe;

translated deep probe + two anchored zero planes
  -> one primitive-deep / nonzero paired-endpoint colour.       (1)
```

This is the exact object requested at the end of THM-2559, but only at the
partial-bare cospan level of THM-2445.  Section 5 proves why calling it a
completed target coefficient would still lose one residue.

## 1. The fully paired three-parameter table

Put `p=13` and, for `L in {1,2}`, write

```text
d_L(y)=1_(||y||<L/14),             u_L=1-d_L.          (2)
```

Let `a,c,k_a,k_b` be the two THM-2350 blocker--graft pairs

```text
eta=e_a-e_(k_a),             ell=e_c-e_(k_b),

c=c_3,                       k_b=q_*.                 (3)
```

Both blockers use `L=1`.  Put `L_a,L_b in {1,2}` for the role types of
`k_a,k_b`; `L=2` means that role is the unique guard.  Since the roles are
distinct, at most one of `L_a,L_b` equals two.

Let `w>=0` be a rational step weight of mass

```text
rho=integral_T w(x)dx>0                                  (4)
```

whose support obeys

```text
d_1(c x)=0,                    d_(L_a)(k_a x)=1.          (5)
```

The THM-2559 application has the additional old-head anchors

```text
d_1(a x)=0,                    d_1(j x)=1,                (6)
```

but they are not needed for the quantitative argument.

Use THM-2365's sign convention and define, for `r,s,t in F_13`,

```text
R_(r,s,t)
 =integral_T w(x)
    d_1(c x-r/13)

    u_1(a x-s/13) u_(L_a)(k_a x+s/13)

    u_1(c x-t/13) u_(L_b)(k_b x+t/13) dx.           (7)
```

The `s` pair is exactly the `eta` dipole and the `t` pair exactly the `ell`
dipole.  Thus every positive cell contains both blockers and both graft
roles safe on the moving endpoint.  No member of either pair on that moving
endpoint is frozen; the old-leg anchors inside `w` remain fixed.

## 2. Exact capacity and three zero loci

THM-2379's translated-tooth identity is

```text
sum_(q in F_13)d_L(y-q/13)=2L-d_L(13y).             (8)
```

Consequently the numbers of safe translates of the two factors in one
paired axis are at least `11` and `13-2L`.  Inclusion--exclusion on the
thirteen shifts gives

```text
sum_s u_1(a x-s/13)u_L(kx+s/13)
 >=(13-2)+(13-2L)-13
 =11-2L.                                             (9)
```

Hence the `s`-pair capacity is at least `11-2L_a`, the `t`-pair capacity
is at least `11-2L_b`, and

```text
sum_r d_1(c x-r/13)=2-d_1(13c x)>=1.                (10)
```

The three shift variables are independent.  Tonelli and (9)--(10) give

```text
M:=sum_(r,s,t)R_(r,s,t)
 >=(11-2L_a)(11-2L_b)rho.                           (11)
```

There is only one guard, so

```text
(11-2L_a)(11-2L_b)
 =81,                     both roles ordinary,

 =63,                     exactly one is the guard. (12)
```

In particular `M>=63rho>0` uniformly.

The anchors (5) give two whole zero planes:

```text
R_(0,s,t)=0                 for every s,t,

R_(r,0,t)=0                 for every r,t.          (13)
```

The repeated `c` factor in (7) also gives the THM-2365 diagonal zero

```text
R_(t,s,t)=0                 for every s,t.          (14)
```

Equations (11), (13), and (14) show directly that the table is nonzero and
cannot be a circulant function of `r-t`: a circulant table with the whole
`s=0` plane zero would vanish identically.

## 3. A simultaneous deep and paired-endpoint colour

Let `zeta=exp(2 pi i/13)` and use the normalized full transform

```text
Rhat(A,B,C)
 =1/13^3 sum_(r,s,t)R_(r,s,t)zeta^(Ar+Bs+Ct).       (15)
```

Character orthogonality and the two zero planes in (13) give

```text
sum_(A!=0,B!=0)Rhat(A,B,0)
 =M/13^3.                                           (16)
```

Indeed `sum_(A!=0)zeta^(Ar)=-1` for `r!=0`; both minus signs multiply to
one, and the zero-coordinate terms vanish.  There are `12^2` summands, so
some `A,B in F_13^*` satisfy

```text
Re Rhat(A,B,0)
 >=M/[13^3*12^2]                                    (17)

 >=9rho/35152,             both roles ordinary,

 >=7rho/35152,             exactly one guard.       (18)
```

Equivalently, marginalize the second dipole without normalizing,

```text
K_(r,s)=sum_t R_(r,s,t),

Khat(A,B)=1/13^2 sum_(r,s)K_(r,s)zeta^(Ar+Bs).      (19)
```

Then `Khat(A,B)=13 Rhat(A,B,0)` and the corresponding floors are

```text
9rho/2704,                    both ordinary,

7rho/2704,                    exactly one guard.     (20)
```

The colour `A` is a nonzero deepest-probe residue.  At the Poisson--Abel
boundary, every contributing exact deep multiplier obeys

```text
m=A mod 13,

d_hat_1(m)=sin(pi m/7)/(pi m)!=0.                   (21)
```

Thus `13` and `7` divide no live `m`, and every displayed deep multiplier
is coprime to `91`.  As in THM-2379, this statement controls the displayed
probe multiplier only, not all Fourier indices hidden in `w`.

## 4. Uniform application to the 165 old heads

Take a positive slope stratum from THM-2559 and put

```text
w=T_delta.
```

At its old target-informed head,

```text
c_3 safe,                 a safe,                 j dangerous,

k_a dangerous.                                           (22)
```

The three blocker statuses are root-constant, so (5)--(6) hold.  Choose
the THM-2350 pivot with `k_b=q_*`; equations (3) and (7) are then the two
actual paired dipoles.  Sections 1--3 apply on every one of the `165` rows.

If desired, partition `w` by the five remaining guard/unit truth bits.  One
of at most `32` pieces has mass at least `rho/32` and is a fully
truth-resolved old-head atom: all blocker bits and all six role bits,
including the separately retained `q_*=k_b` core bit, are then fixed.  This
is a refinement of, not a synonym for, THM-2452's seven-bit complete-mask
bank.  The relative bounds (18)--(20) remain unchanged on that piece.  The
refinement does not solve endpoint matching; it makes the mismatch in
Section 5 literal.

## 5. Why this is not yet a completed target current

Equation (7) repairs the exact defect named in THM-2461: the moving `k_a`
complement no longer travels without `a`, and the moving `k_b` complement
no longer travels without `c`.  It is therefore a lawful **moving-endpoint
dipole table** and a positive partial-bare cospan in THM-2445's sense.

The fixed factor `w`, however, is the old selector/head leg.  It is not
co-shifted in (7).  In an atomic endpoint refinement, let `u` be a left
index coming from that fixed leg and `v` a moving endpoint index.  The
`B`-character in (15) sees only the moving residue (up to the global choice
of endpoint/Fourier sign)

```text
plus_or_minus eta.v.                                  (23)
```

A completed THM-2334 target character must see the difference

```text
eta.(u-v).                                            (24)
```

Thus (23) implies (24) only after a left-residue sidecar such as
`eta.u=0` is proved.  THM-2559's target-informed singleton selector is not
known to be target-neutral or to possess a covariant target action; in its
native root chart, a coordinate index `u` is not even part of the retained
data.  The missing information is therefore the left endpoint residue, not
another right paired factor.

THM-2452 does not add it for free.  On the fixed head, `k_a` is dangerous;
on every positive moving repair cell in (7), `k_a` is safe.  After complete-
mask refinement these are different Boolean atoms.  They form an
off-diagonal endpoint transition, precisely the kind whose fixed-frequency
terms can be nonzero while its full-`X` recombination cancels.  The matched-
mask idempotent slide preserves diagonal atoms, not this danger-to-safe
transition.

The loss is sharp even in the one-axis finite group algebra.  Put

```text
K(s)=1_(s!=0),

c_0=12/13,               c_v=-1/13 for v!=0.        (25)
```

Then

```text
K(s)=sum_v c_v zeta^(-vs)                           (26)
```

is nonnegative, has an anchored zero at `s=0`, and has every nonzero
one-sided character.  Regard each summand as an endpoint pair with equal
left and right residues `u=v`.  Its canonical difference residue is then
zero term by term, and

```text
sum_v c_v zeta^((u-v)s)=sum_v c_v=0                 (27)
```

for every `s`.  A positive one-sided boundary table can therefore coexist
with complete zero-target cancellation.

## 6. Positive physical control on the blind interval

THM-2561 supplies a fully physical check.  On its nested owner interval

```text
J=(4117/399854,4129/399854),

x=(z+3)/13,                  z in J,                (28)
```

take

```text
c=13^5,              a=169,

k_a=95,              k_b=93,                       (29)
```

and let `w` be the branch interval in (28).  Its exact mass is

```text
rho=6/2599051.                                       (30)
```

Both graft roles are ordinary.  The exact referee splits (28) at all `22`
relevant strict-comb walls and integrates all `13^3` cells of (7).  It finds

```text
positive cells:                    1330,

minimum (deep,s-pair,t-pair)
capacities:                        (1,10,10),

sum_(r,s,t)R_(r,s,t)
 =2110/4826809
 =(7385/39)rho,                                      (31)

max_(r,s,t)R_(r,s,t)=rho/6.                          (32)
```

Both zero planes and the diagonal in (13)--(14) vanish exactly.  This
control is stronger than the universal capacity floor, but it remains the
same local non-cover example as THM-2561.

## 7. Exact frontier after the paired completion

The proved implication is

```text
target-informed old head
  + both paired moving-endpoint dipoles
  + primitive deepest probe
  -> positive deep / paired-endpoint partial-bare colour.    (33)
```

The invalid promotion is

```text
paired moving-endpoint colour
  -/-> left-minus-right target residue
  -/-> completed matched endpoint
  -/-> terminal-word target fibre
  -/-> genuinely later semantic root.                       (34)
```

The next exact object is a **charged endpoint transition cospan** retaining
the old selector's left residue together with the moving repair residue.
Either a nonzero difference character survives, or its full-`X`
recombination gives a concrete off-diagonal cancellation law.  Only after
that step could a THM-2365/2452 completed-endpoint landing be attempted;
common-endpoint and semantic typing would still require audit.  The
sheet/carry/future-root and Hall obligations remain subsequent.  No scalar
row is removed and LRC(14) remains open.

## 8. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_paired_dipole_deep_corner_thm2563.py
python3 -O 04-computation/lrc14_paired_dipole_deep_corner_thm2563.py
```

Both executions byte-match the stored output.  The referee checks the exact
ordinary/guard tooth counts, paired capacity products and rational floors;
all basis cells of the double-zero character identity; a nontrivial exact
integer table; the missing-left-residue hostile (25)--(27); and the complete
rational physical table (28)--(32).  The symbolic Tonelli argument and the
partial-bare versus completed-target typing are proved above, not inferred
from the controls. **QED.**
