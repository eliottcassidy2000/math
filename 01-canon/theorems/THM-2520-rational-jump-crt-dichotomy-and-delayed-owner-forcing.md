---
id: THM-2520
title: "Rational-jump CRT dichotomy and delayed-owner forcing"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Once a Perron
  depth contains the full 13-primary part of a rational step function's
  endpoint grid, one nontrivial last-digit ladder vanishes if and only if
  all of them vanish, if and only if the Perron average is exactly constant.
  The equality branch is decided by one finite vector: the endpoint jump
  current aggregated modulo the prime-to-13 part of the grid.  Off that
  branch all twelve ladders are globally nonzero.  Any positive BV
  owner--word event, delayed sufficiently far beyond the collision, then
  gives all twelve weighted collision colours simultaneously; an explicit
  endpoint/variation bound supplies one common delay.  A pure-13 interval
  is a sharp zero-drift hostile.  The theorem does not orient antipodal
  cospans, rebase old source/deep sheets, prove nonzero mixed ANOVA drift on
  a live response table, exclude a scalar row, or prove LRC(14).
source: codex-2026-07-27-rational-jump-crt-dichotomy
depends_on:
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary
related:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
script: 04-computation/lrc14_rational_jump_crt_owner_forcing_thm2520.py
output: 05-knowledge/results/lrc14_rational_jump_crt_owner_forcing_thm2520.out
script_sha256: 059718af2c3c87ed0979136467d21264a01bb50339839c16190f38cf364acab1
output_sha256: dbaab5fce77bf471a7dc5967dbf25e76e8568bea2ef31836bf4366d357fc12d2
hash_basis: working-tree bytes (LF)
---

# THM-2520 -- the non-thirteen endpoint current decides delayed collision drift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2519 identifies last-digit collision drift with the weighted squared
norms of twelve Perron frequency ladders.  Its sharp hostile is measurable
at the marked future scale, so those norms can all vanish even when the bulk
square table is fully charged.  Two questions remain:

```text
when can a whole last-digit ladder vanish globally?

if it does not vanish globally, can the future owner be made to see it?   (1)
```

For rational step data both questions have exact answers.  Strip the
`13`-primary part from one common endpoint denominator and add the jumps in
each remaining residue class.  That finite vector is zero exactly on the
future-factor branch.  Otherwise CRT forces every last digit, and one more
owner delay supplies their common positive weight.

## 1. Endpoint grid, jumps, and Perron ladders

Let `F:T->R` be a periodic step function.  Choose a common endpoint grid

```text
D=13^K D_0,                    gcd(D_0,13)=1,                   (2)
```

so every discontinuity of `F` lies in `D^(-1)Z/Z`.  Minimality of `D` is
not required.  At `j/D`, use the cyclic jump convention

```text
Delta_j=F(j/D+)-F(j/D-),             j in Z/DZ.                 (3)
```

The total jump is zero.  Aggregate the jump current after forgetting the
`13`-primary sheet:

```text
C_r=sum_(j congruent r mod D_0)Delta_j,
                                      r in Z/D_0Z.              (4)
```

Fix

```text
M=13^m,                         m>=K,

h=P_MF,

(P_MF)(y)=1/M sum_(r=0)^(M-1)F((y+r)/M),

A=integral_T F.                                                (5)
```

Let `zeta=exp(2 pi i/13)`.  The twelve nontrivial last-digit projections are

```text
H_a(y)=1/13 sum_(u in F_13)zeta^(-au)h(y+u/13),

                                      a in F_13^*.              (6)
```

As in THM-2519, `H_a` contains exactly the Fourier frequencies of `h` which
are congruent to `a mod 13`.

## 2. The exact jump--CRT dichotomy

The following conditions are equivalent:

```text
H_(a_*)=0 almost everywhere for some a_* in F_13^*;            (7a)

H_a=0 almost everywhere for every a in F_13^*;                 (7b)

C_r=0 for every r in Z/D_0Z;                                  (7c)

P_MF=A almost everywhere.                                     (7d)
```

Consequently exactly one of the following occurs:

```text
constant branch:  h=A and all twelve last-digit ladders vanish;

charged branch:   H_a!=0 in L^2 for every a in F_13^*.         (8)
```

### Proof

Use the Fourier convention

```text
Fhat(n)=integral_T F(x)exp(-2 pi i n x)dx.                     (9)
```

Distributional integration by parts gives, for `n!=0`,

```text
2 pi i n Fhat(n)
 =sum_(j mod D)Delta_j exp(-2 pi i n j/D).                    (10)
```

Put `n=Mk`.  Since `m>=K`, grouping (10) by `j mod D_0` gives

```text
2 pi i Mk Fhat(Mk)
 =sum_(r mod D_0)C_r
    exp(-2 pi i 13^(m-K) k r/D_0).                            (11)
```

The Perron Fourier identity is

```text
hhat(k)=Fhat(Mk).                                              (12)
```

Suppose (7a) holds.  Then the right side of (11) vanishes for every
integer `k` satisfying

```text
k=a_* mod 13.                                                  (13)
```

As `k` ranges through (13), CRT says that `k mod D_0` ranges through all of
`Z/D_0Z`.  Multiplication by `13^(m-K)` is also a permutation modulo `D_0`.
Thus the `D_0`-point DFT of `(C_r)` vanishes at every character.  DFT
invertibility gives (7c).

If (7c) holds, (11)--(12) give

```text
hhat(k)=0                         for every k!=0.               (14)
```

Hence `h=A` in `L^2`, proving (7d), which immediately gives (7b) and (7a).
The reverse implications are included in this cycle.  This proves (7)--(8).

The depth condition `m>=K` is load-bearing.  Before the Perron depth has
passed the full `13`-primary endpoint grid, the residue reduction in (11)
does not exist and a pure-`13` future factor can remain nonconstant.

## 3. An endpoint-explicit energy floor

The dichotomy has a quantitative form.  Put

```text
S_C=sum_(r mod D_0)C_r^2.                                    (15)
```

Assume `S_C>0`.  For

```text
Q_t=sum_r C_r exp(-2 pi i tr/D_0),
```

finite Parseval gives

```text
sum_(t mod D_0)|Q_t|^2=D_0 S_C.                               (16)
```

Choose `t_0` with `|Q_(t_0)|^2>=S_C`.  For each `a!=0`, CRT supplies a
unique representative

```text
1<=k_a<=13D_0-1,

k_a=a mod 13,

13^(m-K)k_a=t_0 mod D_0.                                     (17)
```

Let

```text
E_a=integral_T |H_a|^2.                                      (18)
```

Equations (11)--(12), Parseval on the `a mod 13` ladder, and (17) imply the
uniform floor

```text
E_a>=|Fhat(Mk_a)|^2

   >=S_C/[4 pi^2 M^2(13D_0-1)^2]

                                      for every a!=0.          (19)
```

This floor is deliberately coarse but completely explicit in the endpoint
current.  No denominator-free floor is possible for arbitrary rational step
heights.

There is also a BV invoice.  Perron averaging contracts variation:

```text
Var(P_MF)<=Var(F)/M.                                          (20)
```

Indeed, each inverse branch traverses one of the `M` equal subintervals;
after the normalized sum, their variations total at most `Var(F)/M`.
Translation averaging in (6) therefore gives

```text
Var(H_a)<=Var(F)/M,

||H_a||_infinity<=||F||_infinity,

Var(|H_a|^2)
 <=2||F||_infinity Var(F)/M.                                  (21)
```

## 4. One delayed owner forces all twelve colours

Let `G:T->R_(>=0)` be BV with

```text
rho=integral_T G>0.                                           (22)
```

For `R>=1`, put

```text
W_R(y)=G(13^R y),

B^R_u=integral_T W_R(y)h(y)h(y+u/13)dy,

Bhat^R(a)=1/13 sum_(u in F_13)B^R_u zeta^(-au).               (23)
```

The weight is invariant under every last-digit translation.  THM-2519's
norm identity gives exactly

```text
Bhat^R(a)=integral_T W_R(y)|H_a(y)|^2dy>=0.                   (24)
```

The two-BV dilation estimate gives

```text
|Bhat^R(a)-rho E_a|
 <=Var(G)Var(|H_a|^2)/(12*13^R).                              (25)
```

Hence, on the charged branch, one common sufficiently large `R` satisfies

```text
Bhat^R(a)>0                       for every a in F_13^*.        (26)
```

Combining (19), (21), and (25), the following entirely explicit condition
is sufficient:

```text
3*13^R*rho*S_C

 >20 Var(G)||F||_infinity Var(F)
      M(13D_0-1)^2.                                          (27)
```

Indeed, (27) uses the elementary safe bound `pi^2<10` in the sharper
sufficient inequality

```text
13^R
 >(2 pi^2/3)
   Var(G)||F||_infinity Var(F)
   M(13D_0-1)^2/(rho S_C).                                   (28)
```

On the constant branch, no delay can help:

```text
B^R_u=rho A^2                         for every u and R.       (29)
```

Thus the exact eventual equality boundary is

```text
C=0
 iff for every BV G>=0 of positive mean and every R>=1,
     the delayed last-digit drift is zero;

C!=0
 implies every sufficiently delayed positive owner
         has all twelve last-digit colours.                   (30)
```

If `F` and `G` are rational step functions, every `B^R_u` is rational.
The prime-cyclotomic argument of THM-2519 then also gives, at each fixed
`R`,

```text
one Bhat^R(a_*)=0, a_*!=0
 iff every nontrivial Bhat^R(a)=0;

otherwise every nontrivial Bhat^R(a)>0.                       (31)
```

The delayed mixing theorem strengthens (31): when `C!=0`, it proves that
the positive alternative must eventually occur.

## 5. Lawful common-future cospan form

Return to THM-2519's collision notation

```text
N=13M,

d=u+13e,                         e in Z/MZ.                    (32)
```

Equation (23) is the exact reduction of the physical cospan average

```text
mathcal B^R_u
 =1/M sum_e integral_T
    G(13^(R-1)N x)
    F(x)F(x+(u+13e)/N)dx.                                    (33)
```

Both endpoints in every summand have the same `T^L` future, and the later
factor in (33) is evaluated at the same still-further future of both:

```text
G(13^(R-1)N(x+d/N))
 =G(13^(R-1)(Nx+d))
 =G(13^(R-1)Nx).                                              (34)
```

Thus delaying `G` does not freeze one packet role or choose an absolute
inverse branch.  For a THM-2478 owner--word block it is a lawful common
future owner and its terminal word.  The delay `R` can be taken in either
prescribed cofinal parity class after increasing the strict bound (27).

The construction remains a self-cospan.  It is antipodally even in `u`, so
the twelve positive colours occur in six equal opposite pairs.  It supplies
no orientation or tournament edge.

## 6. Sharp controls and equality mechanisms

### Pure-`13` hostile

Take

```text
F=1_[0,1/13),             D=13, K=1, D_0=1,

M=13.                                                            (35)
```

There is one positive and one negative endpoint jump, but after aggregation
modulo `D_0=1`,

```text
C_0=sum_j Delta_j=0.                                         (36)
```

Exactly one of the thirteen inverse samples lies in the interval, so

```text
P_13F=1/13.                                                   (37)
```

Every ladder and every delayed collision drift vanishes.  This is the
minimal future-factor mechanism behind THM-2519's larger
delta-plus-replicas hostile.  It also shows why the condition `m>=K` cannot
be weakened: at a shallower Perron depth the same indicator is nonconstant.

### Prime-to-`13` positive control

Take

```text
F=1_[0,1/2),              D=2, K=0, D_0=2,

M=1.                                                             (38)
```

The two jumps give

```text
(C_0,C_1)=(1,-1),                 S_C=2.                       (39)
```

Hence all twelve global ladders are nonzero, and every sufficiently delayed
positive owner sees all of them.  With `G=1`, no delay is needed.

The two controls isolate the exact coordinate: nonconstancy of `F` is not
enough; the part of its endpoint current surviving the prime-to-`13`
quotient is what matters.

## 7. Exact live-row reduction and nonclaims

At any fixed lawful THM-2449 clock, each response density is a rational step
function on a finite endpoint grid.  To test its global last-digit ladders:

```text
1. choose one common endpoint denominator D=13^K D_0;
2. enumerate only its actual boundary jumps Delta_j;
3. accumulate C_(j mod D_0);
4. test whether this sparse rational vector is zero.           (40)
```

There is no need to expand a `13^K`-cell tower.  When `C!=0`, equation (27)
installs a sufficiently delayed common future owner--word factor and forces
all twelve collision colours on that scalar response.  When `C=0`, the
Perron response is exactly constant and no such delay can create drift.

This reduces the analytic support problem to a finite endpoint-current
test.  It does not finish the table semantics:

- a drifting marked cell, zero-extended to the table, has every primitive
  table character, but this is not an intrinsic ANOVA interaction;
- different response cells can have different jump vectors, and signed
  mixed combinations can still cancel;
- the old deep label may be kept on an ancestry leg, but (34) does not rebase
  it at the future owner;
- a unit collision address eventually violates the old septimal/deep
  phase-bank preservation bounds of THM-2518;
- source-time and arrival-time atoms remain different projections;
- positivity of the root-colour norm does not orient `u` against `-u`; and
- no scalar row is excluded and LRC(14) remains open.

The next cheapest live computation is now precise: form the vectors (4) for
the owner-marked lawful response cells, then test their mixed table
combinations before attempting any further analytic estimate.

## 8. Exact referee

Run

```bash
python3 04-computation/lrc14_rational_jump_crt_owner_forcing_thm2520.py
python3 -O 04-computation/lrc14_rational_jump_crt_owner_forcing_thm2520.py
```

Both runs reproduce the stored transcript byte-for-byte.  The
dependency-free `Fraction` referee checks the cyclic jump aggregation,
direct Perron evaluation on rational grids, CRT coverage of every last
digit, variation contraction, the constant/charged dichotomy, the two sharp
controls, and exact finite conditional-variance identities.

An independent line audit rederived (10)--(14), including the CRT
representative range, the finite-DFT energy floor (19), Perron variation
contraction, the square-variation bound, the common-delay inequalities
(25)--(28), and the universal quantifiers in the zero branch.  It found no
unresolved proof defect.  **QED.**
