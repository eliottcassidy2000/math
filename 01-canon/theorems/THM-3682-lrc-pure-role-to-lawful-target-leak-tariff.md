---
id: THM-3682
title: "LRC pure-role to lawful-target leak tariff"
status: >
  PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.  On a pure
  THM-2305 target word, freeze every factor except the redundant named danger
  factor and compare that forced THM-2362 role profile with the lawful
  THM-2365 co-shift row.  Their difference is an exact covariant-leak
  sequence.  The nonzero-character energy of the lawful row is at least the
  squared positive part of sqrt(role energy)-sqrt(leak energy).  Thus any
  leak energy below 121 rho^2/2028 forces target drift, with explicit sharp
  THM-3674 tariffs.  A translated danger-weight hostile attains equality at
  complete cancellation, so positivity and even quadratic root energy do
  not remove the leak sidecar.  This is a conditional transfer theorem, not
  a covering-row leak bound or LRC(14) proof.
source: kps-s197 / pure-role versus diagonal-coshift comparison, 2026-08-21
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

**PROVED + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.**  This theorem
identifies the missing coordinate between THM-2362's factor-coloured role and
THM-2365's lawful target action.  It is an exact Hilbert-space comparison,
not an assertion that the missing term is small on a hypothetical cover.

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

## 2. Exact typing on a pure LRC word

Put `p=13`.  Fix a THM-2305 pure word `Q_(j,{a})`, with `a` the nondeep
THM-2365 target blocker and `c` the deepest target blocker.  The scalar cover
makes the displayed `a`-danger factor redundant: after deleting it, the
remaining nonnegative density is still supported where that factor equals
one.

Fix the second target coordinate at `t=0` and include the exact successor
weight

```text
chi(x)=2-d(13c x),                 1<=chi<=2.        (7)
```

Push the remaining base density forward through `y=c_a x`.  The finite-root
average gives an integrable nonnegative density `w_0(y)` with

```text
support(w_0) subset D={d=1},
rho=int_T w_0(y)dy>0.                               (8)
```

With the sign inherited from THM-2365, put

```text
A_s(y)=d(y-s/13),
R_s=int_T w_0(y)A_s(y)dy.                           (9)
```

Changing `s` to `-s` does not change any energy.  Hence THM-2362 gives

```text
E(R)>=121 rho^2/2028.                              (10)
```

Now let every other present factor undergo the actual lawful quotient-dual
co-shift.  Delete only the moving `a`-factor and push the resulting remainder
forward in the same way; call its density `w_s`.  The actual successor row is

```text
C_s=S(s,0)=int_T w_s(y)A_s(y)dy.                   (11)
```

Thus the covariant leak has the concrete physical formula

```text
L_s=int_T (w_s(y)-w_0(y))A_s(y)dy.                 (12)
```

This is the promised map.  The source and target use the same word, named
danger factor, successor weight, and target coordinate.  Freezing the other
factors destroys only their covariance; `(12)` records exactly that lost
information.

Combining `(4)` and `(10)` yields the usable gate

```text
E(L)<121 rho^2/2028

  => E(C)>=(11rho/sqrt(2028)-sqrt(E(L)))^2>0.       (13)
```

Consequently the THM-2365 successor marginal is nonconstant and its drift is
positive.

## 3. Sharp target-energy invoice

Let `S(s,t)` be the complete THM-2365 successor table and keep the row
`C_s=S(s,0)`.  The global mean minimizes squared distance on this row, so

```text
Var(S)
 =1/p^2 sum_(s,t)|S(s,t)-S_bar|^2
 >=E(C)/p.                                         (14)
```

The constant `1/p` is sharp: take every other row to be constant at the mean
of `C`.  THM-3674 now gives, at `p=13`,

```text
D_H   >=E(C)/2028,
E_dt  >=E(C)/26364.                                (15)
```

Here `E_dt` is the nonzero-deep, nonzero-target energy in THM-3674.  Under
the strict hypothesis in `(13)`, set

```text
Gamma=(11rho/sqrt(2028)-sqrt(E(L)))^2.
```

Then the entirely explicit consequence is

```text
Var(S)>=Gamma/13,
D_H>=Gamma/2028,
E_dt>=Gamma/26364.                                 (16)
```

Thus the remaining covering-row task is not “make a role harmonic”: that is
already forced.  It is the quantitative assertion that the other factors'
covariant leak cannot be exactly antiparallel with enough energy to erase
that harmonic.

## 4. The cancellation threshold is attained physically

Let

```text
D=(-1/14,1/14) mod 1,             w_0=d,
A_s(y)=d(y-s/13).                                  (17)
```

The frozen overlap profile is exactly

```text
R_0=1/7,
R_1=R_12=6/91,
R_s=0 otherwise.                                   (18)
```

Its nonzero-character energy is

```text
E(R)=2508/1399489
    >121/99372=121 rho^2/2028.                     (19)
```

Now covary the remainder by the same translation:

```text
w_s=A_s.
```

Then

```text
C_s=int A_s(y)^2dy=1/7             for every s,

E(C)=0,             L^o=-R^o,      E(L)=E(R).       (20)
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

int W_kappa>p_a^2.                                  (21)
```

Multiplying by the successor weight `(7)` only increases its mass.  Deleting
the redundant named danger factor and freezing the root-energy weight gives
the rigorously typed quadratic **role** floor

```text
E(R_kappa)>121 p_a^4/2028.                          (22)
```

The root-character action is not automatically the lawful target action.  If
one supplies a covariant quadratic family `w_(kappa,s)` whose base is this
frozen weight and whose `s` coordinate is the THM-2365 quotient-dual action,
then Sections 1--3 apply to its leak.  In that conditional lane, the concrete
sufficient test

```text
E(L_kappa)<=121 p_a^4/2028                          (23)
```

forces a nonconstant lawful **quadratic-energy** profile.  If
`p_a>rho_j/3`, the right side can be weakened to the source-only threshold

```text
121 rho_j^4/164268.                                 (24)
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
  bound E(L) below the role floor -> actual THM-2365 drift;

root-energy word:
  first construct the target/root covariant family,
  then bound E(L_kappa) below (23) -> quadratic drift,
  then prove a polarization/intertwiner sidecar -> linear target current.
                                                               (25)
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
gates.  The theorem does not bound the leak on a hypothetical scalar cover,
does not exclude THM-3670's thirteen-number recirculation locus, does not
identify quadratic and linear targets, and does not prove LRC(14).  **QED
pending independent hostile audit.**
