---
id: THM-3682
title: "LRC pure-role to lawful-target leak tariff"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER TYPING
  REPAIR.  The nonzero-character energy of any covariant profile is at least
  the squared positive part of sqrt(frozen-role energy)-sqrt(leak energy).
  THM-2362 therefore gives a sharp 121 rho^2/2028 gate and THM-3674 gives
  exact drift tariffs for a current-scale present-Q auxiliary table.  The
  original THM-2365 handoff instead contains Q(Rx): because 13 divides R,
  its redundant word factor is target-fixed, so the present-Q table is not
  that delayed row.  A translated danger-weight hostile attains complete
  cancellation.  The abstract tariff and auxiliary realization are proved;
  the target/word intertwiner, covering-row leak bound, and LRC(14) remain open.
source: kps-s197 / pure-role versus diagonal-coshift comparison, 2026-08-21
audit: >
  PASS AFTER ONE LOAD-BEARING SCALE REPAIR -- Lovelace independently checked
  every Fourier and THM-3674 normalization, the sharp interval hostile, and
  the root-energy constants.  The audit caught that the first draft moved a
  redundant factor inside the delayed word Q(Rx); periodicity actually fixes
  it under the lawful target co-shift.  The repaired theorem applies directly
  to an explicitly constructed current-scale present-Q table and records the
  missing target/word-action intertwiner for the original delayed current.
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2362-thirteen-shift-successor-statistic-and-role-jet-floor
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-3674-sharp-successor-variance-drift-and-target-energy-tariff
related:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-3670-lrc-successor-mass-transfer-and-thirteen-number-gate
script: 04-computation/lrc_pure_role_lawful_target_leak_tariff_thm3682.py
output: 05-knowledge/results/lrc_pure_role_lawful_target_leak_tariff_thm3682.out
script_sha256: 7c50823d4479ef02099f5b8fd7cf1d5ea526caa8dea068e03125301b01bdbf97
output_sha256: dfa7c760cc357bab621e25274d69e162fd6b922364b713a22858faf3a9eb8dfa
semantic_sha256: 8b596be28ef491b71c1a8bb4255ffee1e60a3086c4d82a8feebccfcfc7accb7d
hash_basis: raw LF bytes for files; canonical exact hostile ledger for semantic hash
---

# THM-3682 -- the role mode survives exactly when covariance cannot pay for it

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED AFTER TYPING
REPAIR.**  This theorem identifies the missing coordinate between
THM-2362's factor-coloured role and a covariant target profile.  It is an
exact Hilbert-space comparison.  Its physical realization below is a new
current-scale present-word table, not the original delayed handoff current.

## 1. The universal finite-group tariff

Let `p>=2`, let `zeta=exp(2*pi*i/p)`, and for a complex sequence
`X=(X_s)_(s in F_p)` use the normalized transform

```text
X_hat(k)=1/p sum_s X_s zeta^(-ks),

E(X)=sum_(k!=0)|X_hat(k)|^2
    =1/p sum_s |X_s-X_bar|^2.                       (1)
```

Suppose `R` is a frozen one-factor role profile, `C` is the corresponding
lawful covariant profile, and define the **leak**

```text
L_s=C_s-R_s.                                       (2)
```

Projection away from the constant sequence gives `C^o=R^o+L^o`.  Therefore

```text
E(C)=E(R)+E(L)+2 Re <R^o,L^o>,                     (3)

sqrt(E(C)) >= max(0,sqrt(E(R))-sqrt(E(L))).        (4)
```

Equation (4) is just the reverse triangle inequality in the mean-zero
subspace, so it is valid for arbitrary complex profiles and is sharp.  In
particular,

```text
E(L)<E(R)  implies  C is nonconstant.              (5)
```

If `E(L)<=theta E(R)` with `0<=theta<1`, then

```text
E(C)>=(1-sqrt(theta))^2 E(R).                       (6)
```

No sign or phase hypothesis is hidden in this estimate.

## 2. Exact current-scale realization and the delayed-word boundary

Put `p=13`.  Fix a THM-2305 pure word `Q=Q_(j,{a})`, with `a` the nondeep
target blocker and `c` the deepest target blocker.  At the **current scale**
write its indicator as

```text
1_Q(x)=d(c_a x)G_0(x).                               (7)
```

The scalar cover makes the displayed `a`-danger factor redundant, so

```text
support(G_0) subset {x:d(c_a x)=1}.                 (8)
```

Fix the second target coordinate at `t=0` and include the exact successor
weight

```text
chi(x)=2-d(13c x),                 1<=chi<=2.        (9)
```

Push the remaining base density forward through `y=c_a x`.  The finite-root
average gives an integrable nonnegative density `w_0(y)` with

```text
support(w_0) subset D={d=1},
rho=int_T w_0(y)dy>0.                              (10)
```

With the sign inherited from THM-2365, put

```text
A_s(y)=d(y-s/13),
R_s=int_T w_0(y)A_s(y)dy.                          (11)
```

Changing `s` to `-s` does not change any energy.  Hence THM-2362 gives

```text
E(R)>=121 rho^2/2028.                              (12)
```

Construct a **present-Q auxiliary table** by setting the delayed factor to
`W=1` and co-shifting every current-scale factor of `Q` in the chosen lawful
two-target chart.  Write

```text
E^Q_(s,0)(x)=d(c_a x-s/13)G_s(x).                  (13)
```

Delete only the moving `a`-factor and push `G_s(x)chi(x)dx` forward through
`y=c_a x`; explicitly,

```text
w_s(y)=1/|c_a| sum_(c_a x=y)G_s(x)chi(x).           (14)
```

The auxiliary successor row is now exactly

```text
C_s=S_Q(s,0)=int_T w_s(y)A_s(y)dy.                 (15)
```

Thus the covariant leak has the concrete physical formula

```text
L_s=int_T (w_s(y)-w_0(y))A_s(y)dy.                 (16)
```

For this present-Q table, the source and target use the same current-scale
word, named danger factor, successor weight, and quotient-dual coordinate.
Freezing the other factors destroys only their covariance; `(16)` records
exactly that lost information.

The original THM-2365 handoff table is different.  It has

```text
F_(s,0)(x)=E_(s,0)(x)Q(Rx),               13|R.     (17)
```

The redundant danger factor lies inside `Q(Rx)`.  Under the lawful target
co-shift it becomes

```text
d(Rc_a x-Rs/13)=d(Rc_a x),                         (18)
```

so it is fixed by periodicity.  Pushing through `y=c_a x` produces `d(Ry)`,
not `A_s(y)`.  Therefore `(15)` does **not** type the original delayed row.
Relating `S_Q` to that row requires a target/word-action intertwiner or
Bockstein sidecar not proved here.

Combining `(4)` and `(12)` yields the usable gate

```text
E(L)<121 rho^2/2028

  => E(C)>=(11rho/sqrt(2028)-sqrt(E(L)))^2>0.       (19)
```

Consequently the present-Q successor marginal is nonconstant and its lawful
auxiliary drift is positive.  No conclusion about the delayed handoff row
follows without the intertwiner after `(18)`.

## 3. Sharp target-energy invoice

Let `S_Q(s,t)` be the complete present-Q auxiliary successor table and keep
the row `C_s=S_Q(s,0)`.  The global mean minimizes squared distance on this
row, so

```text
Var(S_Q)
 =1/p^2 sum_(s,t)|S_Q(s,t)-S_Q_bar|^2
 >=E(C)/p.                                         (20)
```

Retaining the deep target complement in `G_s` and inserting the shifted deep
danger probe gives the same exact diagonal zero `H_Q(t,s,t)=0` as in
THM-2365.  The constant `1/p` is sharp: take every other row to be constant
at the mean of `C`.  THM-3674 now gives, at `p=13`,

```text
D_Q   >=E(C)/2028,
E_dt,Q >=E(C)/26364.                               (21)
```

Here `E_dt,Q` is the nonzero-deep, nonzero-target energy of this auxiliary
tensor.  Under the strict hypothesis in `(19)`, set

```text
Gamma=(11rho/sqrt(2028)-sqrt(E(L)))^2.
```

Then the entirely explicit consequence is

```text
Var(S_Q)>=Gamma/13,
D_Q>=Gamma/2028,
E_dt,Q>=Gamma/26364.                               (22)
```

For the auxiliary table, the remaining task is to keep the other factors'
covariant leak from becoming exactly antiparallel.  For the original handoff,
an earlier task comes first: construct the target/word intertwiner in `(18)`.

## 4. The cancellation threshold is attained physically

Let

```text
D=(-1/14,1/14) mod 1,             w_0=d,
A_s(y)=d(y-s/13).                                  (23)
```

The frozen overlap profile is exactly

```text
R_0=1/7,
R_1=R_12=6/91,
R_s=0 otherwise.                                   (24)
```

Its nonzero-character energy is

```text
E(R)=2508/1399489
    >121/99372=121 rho^2/2028.                     (25)
```

Now covary the remainder by the same translation:

```text
w_s=A_s.
```

Then

```text
C_s=int A_s(y)^2dy=1/7             for every s,

E(C)=0,             L^o=-R^o,      E(L)=E(R).       (26)
```

This is an exact interval-set hostile, not a rounding or endpoint artifact.
It attains equality in `(4)` and proves there is no universal positive
constant `c` for which `E(C)>=cE(R)` follows from nonnegativity, danger
support, or a large role floor alone.  The hostile is not asserted to be a
scalar-cover owner packet; its job is to certify the sharp abstract boundary.

## 5. Quadratic root energy: gain and remaining sidecar

THM-2305 supplies, on the same pure word and for every nonzero root character
`kappa`, a positive energy density

```text
W_kappa=1_Q |M_kappa|^2,

int W_kappa>p_a^2.                                  (27)
```

Multiplying by the successor weight `(9)` only increases its mass.  Deleting
the redundant named danger factor and freezing the root-energy weight gives
the rigorously typed quadratic **role** floor

```text
E(R_kappa)>121 p_a^4/2028.                          (28)
```

The root-character action is not automatically the lawful target action.  If
one supplies a covariant quadratic family `w_(kappa,s)` whose base is this
frozen weight and whose `s` coordinate is the THM-2365 quotient-dual action,
then Sections 1--3 apply to its leak.  In that conditional lane, the concrete
sufficient test

```text
E(L_kappa)<=121 p_a^4/2028                          (29)
```

forces a nonconstant lawful **quadratic-energy** profile.  If
`p_a>rho_j/3`, the smaller source-only threshold gives the stronger sufficient
leak hypothesis

```text
121 rho_j^4/164268.                                 (30)
```

This does not by itself land the original linear marked current.  The
integral of `|M_kappa|^2` over a moving word is an auxiliary tensor-square
period, not the modulus square of the integrated current.  A polarization or
action-intertwining sidecar must still show that nonzero quadratic target
energy cannot be made entirely from cross terms whose linear target fibres
cancel.  THM-2337's vertical-jet hostile shows that this distinction is
load-bearing.

The theorem therefore makes the quadratic proposal decidable without
overclaiming it:

```text
ordinary pure word:
  present-Q table + bound E(L) below the role floor -> auxiliary drift;

original delayed word:
  first construct a target/word-action intertwiner;

root-energy word:
  first construct the target/root covariant family,
  then bound E(L_kappa) below (29) -> quadratic drift,
  then prove a polarization/intertwiner sidecar -> linear target current.
                                                               (31)
```

## 6. Exact companion and scope

The dependency-free companion verifies the exact danger-overlap profile,
the role floor, the covariant equality hostile, polarization and Cauchy on
`64` deterministic rational controls, the sharp row-to-table invoice, and
both THM-3674 constants.  Reproduce with

```bash
python -B 04-computation/lrc_pure_role_lawful_target_leak_tariff_thm3682.py
python -O -B 04-computation/lrc_pure_role_lawful_target_leak_tariff_thm3682.py
```

Both executions agree with the pinned LF transcript and pass `169` active
gates.  The theorem does not identify the present-Q auxiliary table with the
original delayed handoff, bound either leak on a hypothetical scalar cover,
exclude THM-3670's thirteen-number recirculation locus, or identify quadratic
and linear targets.  The independent audit found and repaired exactly this
current-scale/delayed-scale seam.  LRC(14) remains open.  **QED.**
