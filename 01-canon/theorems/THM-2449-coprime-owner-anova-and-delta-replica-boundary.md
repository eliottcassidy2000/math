---
id: THM-2449
title: "Coprime owner ANOVA and delta-replica boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. A rational
  C_7-by-C_13 lawful owner table has an all-or-all
  mixed spectrum: if one charged source-target coefficient vanishes,
  all 72 vanish and the table is additive. If its untwisted target
  column is the positive owner anchor a delta_0, the unique mixed-zero
  form is A_(0,s)=a+w_s and A_(ell,s)=w_s for all ell nonzero. Thus one
  rectangle defect forces all 72 colours. For rational finite-step
  source/word packets, each clock class has the exact form M+C/R, so
  persistent replica failure is equivalent to two finite additive
  tests. This sharpens the post-THM-2445/2442/2448 endpoint obstruction
  but does not prove a rectangle defect, semantic repair alignment,
  same-root cospan physicalization, a row exclusion, or LRC(14).
source: codex-2026-07-26-owner-delta-replica
depends_on:
  - THM-2340-owner-word-anova-target-landing
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2441-septimal-ancestry-event-period-collapse
  - THM-2442-delayed-word-septimal-source-completion
  - THM-2445-twenty-four-cell-graft-owner-conditioning
related:
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2448-right-endpoint-cospan-transition-atlas
script: 04-computation/lrc14_owner_delta_replica_anova_thm2449.py
output: 05-knowledge/results/lrc14_owner_delta_replica_anova_thm2449.out
script_sha256: 6ad02a29f2cfb9136e604f5785fb53853fc0ec10c88ca9b636a1ce738b666aff
output_sha256: 06d27a2247cb0941fb749c560ae9171f89efed319cb0bfc6b48e41a55a54007a
hash_basis: working-tree bytes (LF)
---

# THM-2449 -- one owner rectangle defect forces every mixed colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2445 reduces the surviving LRC(14) graft to `24` positive
partial-endpoint cells. THM-2442 restores a real delayed word on the
unique ghost. The other `23` cells retain literal blocker or
guard/unit-repair labels, but their semantic word alignment remains
open. THM-2448 expands every fixed current through the missing right
factors into a finite complete-mask/transition cospan; same-root
physicalization or lawful transition-to-repair recombination remains
open.

The operation-ready response is not one drift norm. Whenever a
continuation retains a lawful source phase and a first-target phase, it
has a rational

```text
C_7 x C_13 response table.
```

This theorem classifies the complete mixed-zero boundary of that table.
The hostile is unique after the owner anchor:

```text
one positive delta row
  + six identical nonnegative replicas.                         (1)
```

Consequently the next no-cancellation test is one rectangle, not all
`72` mixed Fourier colours.

## 1. The coprime response table

Let

```text
A=(A_(ell,s))_(ell in F_7, s in F_13) in Q^(7 x 13)
```

and define the normalized transform

```text
Ahat(kappa,b)
 =1/91 sum_(ell,s)
   A_(ell,s) zeta_7^(kappa ell) zeta_13^(b s).        (2)
```

The signs in (2) are immaterial. The **mixed locus** is

```text
kappa!=0,                         b!=0,
```

and contains `6*12=72` coefficients.

### Coprime Galois all-or-all law

Exactly one of the following occurs:

```text
Ahat(kappa,b)!=0 for every one of the 72 mixed pairs;

Ahat(kappa,b)=0  for every one of the 72 mixed pairs. (3)
```

Indeed,

```text
Q(zeta_7,zeta_13)=Q(zeta_91),

Gal(Q(zeta_91)/Q)
 =(Z/7Z)^* x (Z/13Z)^*.                              (4)
```

The group in (4) acts transitively on
`F_7^* x F_13^*`. Because every entry of `A` is rational, a Galois
automorphism sends any one value in (2) to any other mixed value.
Thus one zero forces all `72` zeros, and one nonzero forces all `72`
nonzeros.

Rationality is load-bearing. For arbitrary complex tables, one mixed
frequency can vanish without any of its non-Galois-conjugate peers
vanishing.

## 2. Mixed zero is exactly additive ANOVA

Define the row means, column means, and grand mean

```text
R_ell=1/13 sum_s A_(ell,s),

C_s=1/7 sum_ell A_(ell,s),

mu=1/91 sum_(ell,s)A_(ell,s).                       (5)
```

The interaction is

```text
I_(ell,s)=A_(ell,s)-R_ell-C_s+mu.                   (6)
```

As in THM-2340, `I` has zero row and column means. Its Fourier transform
is exactly the mixed block of (2), and two-dimensional Parseval gives

```text
sum_(kappa!=0,b!=0)|Ahat(kappa,b)|^2
 =1/91 sum_(ell,s)|I_(ell,s)|^2.                   (7)
```

Therefore the following are equivalent:

```text
one mixed coefficient vanishes;

every mixed coefficient vanishes;

I=0;

A_(ell,s)=f_ell+g_s for some rational f,g;

every rectangle
 A_(ell,s)-A_(ell,0)-A_(0,s)+A_(0,0)
 vanishes.                                          (8)
```

The first-to-second implication is the coprime Galois law; the
remaining equivalences are the exact ANOVA reconstruction.

## 3. The owner anchor has one sharp hostile

Assume now that the untwisted target column is the owner anchor

```text
A_(ell,0)=a 1_(ell=0),                 a>0.          (9)
```

Then the mixed-zero alternatives in (8) are equivalent to the existence
of one rational thirteen-vector `w` such that

```text
w_0=0,

A_(0,s)=a+w_s,

A_(ell,s)=w_s              for every ell!=0.        (10)
```

If `A` is nonnegative, then necessarily

```text
w_s>=0.                                             (11)
```

To prove this, write `A_(ell,s)=f_ell+g_s`. Equation (9) gives

```text
f_ell=-g_0                   for ell!=0,

f_0+g_0=a.
```

Put `w_s=g_s-g_0`. This gives (10), and the converse is immediate.
Nonnegativity gives (11) from any of the six replica rows.

Thus the decisive anchored rectangles are

```text
Delta_(ell,s)
 =A_(0,s)-A_(ell,s)-a,

ell!=0,                         s!=0.               (12)
```

One nonzero value in (12) implies

```text
Ahat(kappa,b)!=0
       for every kappa!=0 and every b!=0.           (13)
```

The test has collapsed from `72` algebraic coefficients to any one of
`72` rational rectangle differences. Failure of all tests is not an
untyped kernel: it is precisely the delta-plus-six-replica table (10).

The hostile is sharp even when the owner row is target-nonflat. Choose
any nonconstant nonnegative `w` with `w_0=0`. Then row zero is
`a+w`, the other six rows are `w`, the source column is a positive
delta, both pure axes can be active, and every mixed coefficient is
zero.

## 4. The exact lawful LRC table

Use the source label `j`, source-deleted present packet `U_(s,t)`,
deep probe

```text
Delta_r(x)=d(c_3x-r/13),
```

and fixed positive THM-2305 word of THM-2407/2442. Fix that word once,
and let `k_0` be large enough for both positive-overlap estimates in
THM-2407 equations (7b)--(7c). At

```text
R=13^k,                         k>=k_0,
```

put

```text
epsilon=R mod 7=(-1)^k,

d_(j,ell)(x)=d(c_jx-ell/7).                         (14)
```

Delete the source-safe factor from the terminal word, calling the
remaining word `T_sigma`, and make the full lawful source/word shift

```text
Q^epsilon_ell(y)
 =T_sigma(y)g(c_jy-epsilon ell/7).                  (15)
```

The word shift in (15) is mandatory: it is exactly THM-2409's
skew-diagonal law. Define

```text
H^R_ell(r,s,t)
 =integral_T
    U_(s,t)(x)d_(j,ell)(x)Delta_r(x)
    Q^epsilon_ell(Rx) dx,

A^R_(ell,s)=sum_r H^R_ell(r,s,0).                   (16)
```

Every entry is a nonnegative rational number. The deepest moving safe
factor is retained, so

```text
H^R_ell(t,s,t)=0.                                   (17)
```

At the untwisted target slice, scalar cover gives

```text
U_(0,0)g(c_jx)=0
```

almost everywhere. The seven translated source dangers partition one,
so `U_(0,0)` is supported on `d_(j,0)` and disjoint from
`d_(j,ell)` for `ell!=0`. The fixed word and the THM-2407
eventual-overlap estimate give positive owner overlap at every such
clock. Hence (16) has the exact anchor

```text
A^R_(ell,0)=a_R 1_(ell=0),             a_R>0.        (18)
```

In THM-2407's owner branch, its row

```text
A^R_(0,s)
```

is nonflat. Equations (10)--(13) now give the complete owner alternative:

```text
all 72 lawful source-target colours survive;

or

the actual delayed-word table is exactly delta plus six replicas. (19)
```

This is a statement about the full lawful present/word skew, not the
gauge-invalid operation of shifting only the present source factor.

### Deep/frequency lift

For completeness, retain the deep Fourier colours:

```text
B^R(kappa,alpha,b,tau)
 =1/(7*13^3) sum_(ell,r,s,t)
    H^R_ell(r,s,t)
    zeta_7^(kappa ell)
    zeta_13^(alpha r+b s+tau t),

J^R(kappa,alpha,b)=sum_tau B^R(kappa,alpha,b,tau).   (20)
```

Equations (2), (16), and the `t=0` projection give

```text
J^R(kappa,0,b)=Ahat^R(kappa,b)/13.                  (21)
```

Equation (17) gives

```text
sum_alpha J^R(kappa,alpha,b)=0.                     (22)
```

Therefore, in the first branch of (19), each prescribed
`kappa,b!=0` forces some `alpha!=0` and `tau` with

```text
B^R(kappa,alpha,b,tau)!=0.                          (23)
```

The exact-frequency expansion of THM-2365/2442 then supplies an atomic
term with nonzero source residue modulo seven, nonzero first-target
colour modulo thirteen, and

```text
m=alpha mod 13.
```

The centred deepest-danger coefficient vanishes at nonzero multiples of
seven, so every live term has

```text
gcd(m,91)=1.                                        (24)
```

The exact atom may depend on `(kappa,b)`. No common triangle or
all-coordinate-unit address is asserted.

## 5. Exact finite-clock `1/R` law

The infinite delayed-clock part of (19) has a stronger functional form
than a BV error bound.

Let

```text
E=disjoint_union_i [A_i/D,B_i/D),

D=13^K D_0,                      gcd(D_0,13)=1,

S=sum_i(B_i-A_i)=D measure(E),                       (25)
```

and let `Q` be any fixed measurable terminal set. For

```text
R=13^K N
```

put

```text
I_R(E,Q)=integral_T 1_E(x)1_Q({Rx}) dx.             (26)
```

If

```text
N'=N+D_0 t,                         t integer,

R'=13^K N'>0,                                       (27)
```

then

```text
R'[I_(R')(E,Q)-measure(E)measure(Q)]

 =R[I_R(E,Q)-measure(E)measure(Q)].                 (28)
```

To prove (28), disintegrate `y={Rx}`. The prefixes contributed by
`[A_i/D,B_i/D)` are the integers in

```text
[ceil(NA_i/D_0-y),ceil(NB_i/D_0-y)).
```

Under (27), the two endpoints shift by `tA_i,tB_i`, so the total prefix
count increases by `t(B_i-A_i)`. Hence

```text
R'I_(R')-RI_R=tS measure(Q)

              =(R'-R)measure(E)measure(Q),          (29)
```

which is (28). Rational finite-step weights follow by linearity. The
denominator of `Q` does not enter.

For `R_d=13^(K+d)`, the scaled covariance in (28) depends only on

```text
13^d mod D_0
```

and has period dividing

```text
ord_(D_0)(13).                                      (30)
```

This is the uncoloured scalar counterpart of THM-2418/2441's full
septimal ancestry classifier; its modulus can be smaller because the
seven prefix residues have been summed. For `D_0=1`, interpret the
period in (30) as one.

Apply (28) entrywise to the finite family in (16), using one common
source endpoint denominator. On each fixed clock/word-parity class
`rho`, there are exact rational matrices `M_rho,C_rho` such that

```text
A^R=M_rho+C_rho/R.                                  (31)
```

The rectangle functional kills no hidden terms:

```text
Delta_(ell,s)(R)
 =Delta_(ell,s)(M_rho)
  +Delta_(ell,s)(C_rho)/R.                          (32)
```

Consequently:

- if either coefficient on the right of (32) is nonzero, that rectangle
  can vanish at at most one clock in the class and is nonzero at every
  sufficiently large clock;
- replica failure persists at infinitely many clocks in one class if
  and only if **both** finite matrices `M_rho` and `C_rho` have zero
  interaction.

Thus the complete infinite owner-clock audit is finite:

```text
for each residue class rho,
  test the rectangles of M_rho and C_rho.           (33)
```

If the source phases introduce denominator seven, as in (14), then
`7|D_0`; the period in (30) is even and already retains the physical
reflection `epsilon=(-1)^k`. Otherwise the two parity families are
tested separately.

The “at most one clock” boundary is sharp. A rectangle of the form

```text
1-13/R
```

vanishes at `R=13` and at no other positive base-thirteen clock.

## 6. Post-THM-2445/2442 consequence and boundary

THM-2445's unique ghost has a literal delayed word after THM-2442, at
the same partial endpoint. THM-2448 makes the omitted-right-factor
expansion finite. The remaining frontier consists of semantic
word/repair alignment on the other `23` repair/blocker-labelled cells
and the physicalization or lawful repair conversion of its terminal
cospans.

The present theorem applies whenever one such continuation retains the
lawful table (16) and anchor (18). It says that owner-conditioned
cancellation has only one form:

```text
delta plus six exact replicas,
```

and that persistent realization of this form is equivalent to the
finite double-additivity test (33). A nonzero rectangle in either finite
matrix restores every mixed source-target colour, and then the
`91`-unit deep lift (24), at every sufficiently large clock in that
class. More sharply, there is at most one exceptional positive clock,
the possible root displayed in (32); the finitely many clocks before
`k_0` are direct finite cases.

THM-2448 now proves that every fixed marked current expands through the
omitted right factors into a finite complete-mask/transition cospan.
Its terminal pieces are fixed-frequency complex cross-endpoint
currents, however, not nonnegative rational tables of the form (16).
Consequently THM-2449 cannot be applied separately to a THM-2448
cospan piece merely because that piece retains a source colour. If a
same-root physicalization or lawful transition-to-repair recombination
preserves (16) and its anchor, the only remaining mixed cancellation is
exactly (10), and (33) makes its clock audit finite.

This does **not** prove that a selected THM-2445 typed cell has the
semantic THM-2305 word promised by its local label, that a first-failure
cell reaches its canonical THM-2379 repair endpoint, or that the
THM-2448 cospan can be physicalized while preserving the anchor or a
rectangle. No row is removed and LRC(14) remains open.

## 7. Exact companion

The dependency-free companion:

1. verifies the complete `72`-element CRT/Galois orbit;
2. exhausts all `4,096` binary nonnegative replica profiles and checks
   their anchored normal form and zero mixed interaction;
3. inserts a defect in each of the `72` decisive rectangles and checks
   all `72` mixed colours on every control;
4. checks the ANOVA/Galois dichotomy on `256` secondary signed exact
   matrices;
5. verifies (28) on `2,408` exact congruent-clock pairs, including all
   grid-cell sources through denominator seven and all single intervals
   on the `13` and `26` grids; and
6. realizes the sharp one-exception sequence at clocks
   `13,169,2197`.

Run

```bash
python3 04-computation/lrc14_owner_delta_replica_anova_thm2449.py
python3 -O 04-computation/lrc14_owner_delta_replica_anova_thm2449.py
```

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_owner_delta_replica_anova_thm2449.out
```

after LF normalization. Every load-bearing check remains active under
optimized Python.

## 8. Independent audit

An independent hostile audit reconstructed the coprime Galois orbit,
ANOVA and rectangle maps, the anchored normal form, the deep-frequency
lift, and the scalar covariance law. It found and repaired one
quantifier boundary: positivity in (18) begins at the THM-2407
eventual-overlap threshold `k_0`, and a nonzero matrix rectangle can
have the single exceptional clock allowed by (32). The audit also:

- obtained rank `72` independently for the mixed Fourier map, the
  rectangle map, and their stacked map;
- replayed the normal, optimized, and stored transcripts byte for byte;
- verified all frontmatter hashes and enumerated counts; and
- added `100` independent exact covariance controls.

No unresolved proof defect remains. QED.
