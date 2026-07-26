---
id: THM-2400
title: "Clean-parent root-gauge quotient and target-slope boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the
  complete THM-2392 clean-parent set, common root relabelling cancels
  pointwise in every literal-role charged product. Thus no
  singleton/adjacent status or root-translate cell is needed: each of
  the five actual guard/lower-q roles survives all twelve target
  colours and has a signed Hermitian-pair floor delta/1014, while some
  actual role/colour has floor delta/845. Fixing only one of thirteen
  owner/digit categories retains a nonzero F_7 x F_13 coefficient with
  mass at least delta/13. A lawful target covector descends to this common-root
  gauge at the labelled factor level exactly when all unit-factor
  slopes eta_i/w_i agree; unequal slopes deform collision status and
  admit a sharp two-factor hostile. No lawful target current,
  canonical-endpoint transport, row exclusion, or LRC(14) conclusion
  is asserted.
source: codex-2026-07-26-clean-parent-root-gauge
depends_on:
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
related:
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2384-owner-pivot-primal-dual-obstruction-and-two-probe-repair-corner
  - THM-2397-clean-root-same-parent-charged-role-partition
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
script: 04-computation/lrc14_clean_parent_root_gauge_thm2400.py
output: 05-knowledge/results/lrc14_clean_parent_root_gauge_thm2400.out
script_sha256: 62663a4ecb47a4966729c13fabbdcc9298bf016ee69c1a80f196e553e50aa3e1
output_sha256: 4b9aaa9265764fc0a1e3a52298b44bcd3097f68f2c28ef3089a7037e09161ebe
hash_basis: working-tree bytes (LF)
---

# THM-2400 -- quotient the harmless root gauge before paying for a cell

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2397 fixed the singleton/adjacent status and root translate of the
exclusive `q_*` mask. That is necessary for a separately integrated
linear current, but not for the joint charged product. A common relabelling
of the thirteen roots rotates both Fourier factors by the same phase, which
cancels in their product.

Quotienting this harmless gauge before pigeonholing removes the entire
`2*13` support-category tax:

```text
complete clean parent S
  -> every literal role survives every nonzero target colour;

fix only the transverse owner/digit category
  -> one F_7 x F_13 packet of mass at least delta/13.            (1)
```

This does not make the physical target action a common relabelling. Its
factor slopes generally differ. Section 6 isolates that exact boundary.

## 1. The complete clean parent

Let `S` be THM-2392's positive clean-parent set and put

```text
delta=mu(S)>0,

m_A=integral_S |A_y|dy.                            (2a)
```

Since `|A_y|` is one or two, `m_A>=delta`.

For generic `y in S`, use the inverse roots

```text
x_r=(y+r)/13,                         r in F_13.       (2)
```

Let `A_y` be the exclusive `q_*` mask. It is a singleton when `q_*`
belongs to the unique double pair and a two-root mask otherwise:

```text
|A_y| in {1,2}.                                      (3)
```

For the five other actual roles

```text
L={H} union {q_i:q_i!=q_*},
```

write `C_i(y,r)` for the danger indicator. Their root counts are fixed:

```text
n_H=4,                    n_(q_i)=2.                 (4)
```

Exclusivity gives, pointwise,

```text
A_y(r)C_i(y,r)=0.                                   (5)
```

Use normalized transforms

```text
a_y(k)=(1/13)sum_r A_y(r)zeta^(-kr),

c_(i,y)(k)=(1/13)sum_r C_i(y,r)zeta^(-kr),          (6)
```

where `zeta^13=1`.

## 2. Common root gauge cancels

Let `t=t(y) in F_13` be any measurable choice of root origin and relabel
every mask by

```text
A_y^t(r)=A_y(r+t),

C_(i,y)^t(r)=C_i(y,r+t).                            (7)
```

Then

```text
a_y^t(k)=zeta^(kt)a_y(k),

c_(i,y)^t(k)=zeta^(kt)c_(i,y)(k),                   (8)
```

and therefore

```text
a_y^t(k)conjugate(c_(i,y)^t(k))
 =a_y(k)conjugate(c_(i,y)(k)).                      (9)
```

The joint product is root-gauge invariant even when the gauge varies with
`y`. Define it on the whole clean parent:

```text
G_i(k)
 =integral_S a_y(k)conjugate(c_(i,y)(k))dy.         (10)
```

No status or root-translate refinement appears in (10).

## 3. Rationality forces all twelve colours

Define the normalized circular cross-correlation

```text
K_i(s)
 =(1/13)integral_S sum_r
    A_y(r+s)C_i(y,r)dy.                             (11)
```

The clean set and all masks are finite Boolean combinations of rational
comb intervals, so every `K_i(s)` is a nonnegative rational number.
Equation (5) gives

```text
K_i(0)=0.                                           (12)
```

On the other hand,

```text
sum_s K_i(s)
 =n_i m_A/13
 >=delta n_i/13>0.                                 (13)
```

Let `P_i(X)=sum_s K_i(s)X^s`. A direct substitution gives

```text
P_i(zeta^(-k))=13G_i(k).                            (14)
```

If `G_i(k)=0` for one `k!=0`, then the minimal polynomial

```text
Phi_13(X)=1+X+...+X^12
```

divides `P_i`. Since both degrees are at most twelve, `P_i` is a rational
multiple of `Phi_13`; all `K_i(s)` are equal. Equations (12)--(13)
contradict this. Hence

```text
G_i(k)!=0
 for every actual role i and every k!=0.            (15)
```

This is the direct two-sided mechanism of THM-2398, now applied before
fixing a root gauge.

### 3a. Relative translation, not the common-root diagonal

The distinction from a common-endpoint current is already exact at this
level.  Let `f(y)>=0` be any root-constant rational step function and put

```text
K_i^f(s)
 =(1/13)integral_S f(y)sum_r
    A_y(r+s)C_i(y,r)dy.                            (15a)
```

If `f` has positive mass on `S`, the proof of (12)--(15) is unchanged:
`K_i^f(0)=0`, its total mass is positive, and every nonzero Fourier
colour survives.  In contrast, simultaneous common-root relabelling
gives the diagonal

```text
J_i^f(s)
 =(1/13)integral_S f(y)sum_r
    A_y(r+s)C_i(y,r+s)dy
 =0                                                (15b)
```

for every `s`, directly from (5).  Thus the complete translate bank in
(11) is a **relative endpoint cross-correlation**, not a relabelled
THM-2370 lawful co-shift current. Equation (15b) concerns only the
equal-slope common-root diagonal; it does not identify every THM-2370
target action with that diagonal.

The exact sufficient service is equally transparent.  A lawful polarized
two-endpoint bank

```text
P_i^f(s,t)
 =(1/13)integral_S f(y)sum_r
    A_y(r+s)C_i(y,r+t)dy
 =K_i^f(s-t)                                       (15c)
```

recovers the entire charged spectrum. For this route, either a common
rational nonflat circulant as in THM-2398 or a lawful
relative-shift/polarization bank would suffice. The diagonal `s=t`
alone is the sharp zero hostile (15b).

## 4. Signed floors on the whole parent

Normalized Parseval and (5) give, pointwise in `y`,

```text
sum_(k!=0)a_y(k)conjugate(c_(i,y)(k))
 =-|A_y|n_i/169.                                   (16)
```

After integration,

```text
sum_(k!=0)G_i(k)
 =-n_i m_A/169
 <=-delta n_i/169.                                 (17)
```

Thus every actual role has a negative Hermitian pair and some `k!=0`
with

```text
-Re G_i(k)>=delta n_i/(169*12).
```

In particular,

```text
every lower-q role:       -Re G_i(k)>=delta/1014;

the guard role:           -Re G_H(k)>=delta/507.    (18)
```

The whole role bank is sharper. Pointwise, the sum over the five actual
role masks is `1-A_y` in the singleton status and `1-A_y+D_y` in the
two-root status, where `D_y` is the unique double root and
`A_yD_y=0`. The two normalized Parseval ledgers are

```text
-12/169,                    singleton status;

-24/169,                    two-root status.        (19)
```

Consequently

```text
sum_i sum_(k!=0)Re G_i(k)
 =-12m_A/169
 <=-12delta/169.
```

Among the `5*12=60` literal role/colour pairs, some one satisfies

```text
-Re G_i(k)>=delta/845.                              (20)
```

No aggregate product of separately integrated currents is inferred:
the root phase cancels inside the joint integral (10).

## 5. Exact uniform and owner/digit floors

THM-2396 and the complementary branches give

```text
delta>=1/26754                                      (21)
```

throughout the last septimal lane. Equations (18) and (20) therefore
give

```text
every lower-q role:       -Re G_i(k)>=1/27128556;

some literal role/colour: -Re G_i(k)>=1/22607130.   (22)
```

On the common-core chain,

```text
delta>=66/4459,
```

so the corresponding floors are

```text
every lower-q role:       -Re G_i(k)>=11/753571;

some literal role/colour: -Re G_i(k)>=66/3767855.   (23)
```

THM-2392 partitions `S` into exactly thirteen transverse categories:
one simultaneous-owner category with `d=0`, and the twelve exclusive
`(d,owner)` categories with `d!=0`. Choose a largest category `S_o`;
its mass satisfies

```text
rho_o>=delta/13.                                    (24)
```

On it, the septimal coefficient

```text
j_ell=(1/7)zeta_7^(-ell d)
```

is fixed and nonzero for every `ell`. Repeating Sections 2--4 on `S_o`
shows

```text
j_ell G_i^o(k)!=0
 for every actual role i, every k!=0, and every ell. (25)
```

For every actual role `i`, some Hermitian pair `k_i!=0` obeys

```text
universally:
  rho_o>=1/347802,
  -Re G_i^o(k_i)>=1/352671228,
  |j_ell G_i^o(k_i)|>=1/2468698596;

common core:
  rho_o>=66/57967,
  -Re G_i^o(k_i)>=11/9796423,
  |j_ell G_i^o(k_i)|>=11/68574961.                  (26)
```

Equation (25), but not the numerical floor (26), holds for every
`k!=0`. For some literal role/colour pair, the sharper bank floors are

```text
universally:
  |G_i^o(k)|>=1/293892690,
  |j_ell G_i^o(k)|>=1/2057248830;

common core:
  |G_i^o(k)|>=66/48982115,
  |j_ell G_i^o(k)|>=66/342874805.                   (27)
```

These are joint same-parent coefficients. They are not physical
`91`-unit terminal atoms.

## 6. The exact target-slope boundary

The root-gauge quotient must not be confused with a lawful target
covector. Let the six labelled unit-factor speeds have nonzero residues

```text
w_j in F_13^*
```

and let a target parameter `s` shift factor `j` by coefficient `eta_j`.
On the inverse-root coordinate, its word has the form

```text
C_(j,s)(r)
 =C_(j,0)(r-s alpha_j),

alpha_j=eta_j/w_j in F_13.                          (28)
```

At the **labelled factor-action level**, this is a common root
relabeling exactly when

```text
alpha_1=...=alpha_6.                                (29)
```

If the common value is `alpha`, subtracting the source-translation gauge
`alpha w` from `eta` kills every unit-factor coordinate. Equivalently,
the target class has a representative supported only on the blocker
coordinates. This is the pure-blocker dual boundary.

Unequal slopes genuinely can deform the exclusive word rather than
translate it. On `F_13`, take a two-root failed-factor mask and one
two-root competing mask

```text
D={0,1},                    C={1,2},

alpha_D=0,                  alpha_C=1.              (30)
```

The exclusive deletion supports

```text
A_s=D\ (C+s)
```

have sizes

```text
|A_0|=1,                    |A_12|=0,

|A_1|=2.                                            (31)
```

A common translate preserves cardinality, so this lawful two-slope
action cannot descend to the root-gauge orbit. The example is sharp:
if both slopes are equal, every `A_s` is exactly a translate of `A_0`.

There is also no operator on the collapsed support `A_0` alone which
recovers unequal-slope transport.  Replacing `C` by

```text
C'={1,3}
```

leaves the same collapsed support `D\C=D\C'={0}` at `s=0`, but the two
families differ at `s=12`: the first support is empty and the second is
`{1}`.  The missing datum is the labelled factorization, not another
invariant of the bare exclusive word.

Equation (29) is an iff for descent of the **labelled factor action**.
Special evaluated words can have accidental redundancies, so unequal
slopes are not claimed to deform every individual live cell.

### 6a. The standard owner-pivot target has no common-slope direction

The pure-blocker alternative is already impossible in the standard
THM-2384 owner-pivot packet.  Its two true target duals are

```text
eta_a=e_a-e_(k_a),             eta_b=e_b-e_(k_b),   (32a)
```

where `k_a,k_b` are distinct members of the six unit labels.  A
combination

```text
eta=u eta_a+v eta_b
```

has unit coordinates `-u` at `k_a`, `-v` at `k_b`, and zero at the
other four unit labels.  If all six slopes in (28) agreed, any one of
those four zero coordinates would force the common slope to be zero.
The `k_a,k_b` coordinates would then force `u=v=0`.  Therefore

```text
no nonzero owner-pivot target covector
has a pure-blocker representative.                 (32b)
```

This is a target-dual statement, not the false primal/dual
identification excluded by THM-2384.  It collapses the standard
owner-pivot continuation to the unequal-slope branch.

## 7. Consequence and remaining bridge

The theorem separates two operations which had been conflated:

```text
common root relabelling:
  harmless gauge, quotient before pigeonholing;

lawful target covector:
  generally unequal factor slopes, changes collision data.      (32)
```

It removes the status/translate tax from the clean-parent joint carrier
and reduces the transverse tax from `338` to `13`. It does not turn the
joint coefficient into a lawfully co-shifted THM-2370 packet, preserve
the target quotient/diagonal-zero/Poisson--Abel typing, or move from the
common bare endpoint to the canonical fully masked owner endpoint.

The next exact target is no longer root-gauge alignment.  In the
standard owner-pivot packet, (32b) removes the second alternative
altogether:

```text
control the unequal-slope deformation in (28) while retaining a
lawful relative endpoint bank and the canonical terminal filter. (33)
```

Only outside, or after changing, that standard owner-pivot packet could
a pure-blocker representative remain an option.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 8. Exact companion

The dependency-free exact companion:

- exhausts every disjoint size-`1/2` exclusive mask and size-`2/4`
  actual-role mask on `F_13`;
- verifies common-translate invariance of the full correlation bank;
- checks that the common-root diagonal is identically zero and the polarized
  two-endpoint bank depends only on the relative shift;
- checks the zero base shift, positive total, and all signed ledgers;
- verifies every rational floor in (22)--(27);
- checks the thirteen owner/digit categories; and
- replays the equal-slope positive control and unequal-slope
  cardinality hostile (30)--(31).

Run

```bash
python3 04-computation/lrc14_clean_parent_root_gauge_thm2400.py
python3 -O 04-computation/lrc14_clean_parent_root_gauge_thm2400.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_clean_parent_root_gauge_thm2400.out
```

All truth-bearing executable checks remain active under optimized Python.

Two independent hostile audits rederived the root-gauge cancellation,
the rational cross-correlation argument, all signed floors, the exact
thirteen-category pigeonhole, the relative/diagonal distinction, the
factor-slope iff, and the owner-pivot no-common-slope corollary. Both
audits replayed normal and optimized execution against the stored
transcript and independently reproduced the recorded LF hashes. QED.
