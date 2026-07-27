---
id: THM-2568
title: "Full-X transition annihilation and the refined-pair drift boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A Boolean danger-to-safe endpoint transition
  completed on one common lawful target twist has zero full-X current in
  every target character: endpoint Parseval recombines it to the pointwise
  product P_s(1-P_s)=0.  With independent left/right twists, the refined
  pair spectrum pushes to the coarse THM-2334 target by summing the line
  a+b=q; diagonal zero annihilates every such line.  Diagonal-shift drift is
  exactly the energy off a+b=0, but it detects only refined pair data and
  never defeats the zero coarse pushforward.  The complement hostile
  A(s,t)=1_(s!=t) has all twelve one-sided colours while living wholly on
  a+b=0.  An anchored duplicate-pair table over a fixed THM-2559 head does
  force a nonzero relative-offset colour, with phase-rotated floor
  9 rho/2028 for an ordinary role and 7 rho/2028 for the guard, but it
  freezes the old selector and duplicates its moving gates, so it is not a
  lawful present-endpoint orbit.  Escape requires retaining fixed X, a
  normal/jet weight, or an oriented reference before the diagonal
  pushforward.  No coarse target current, row exclusion, or LRC(14)
  conclusion follows.
source: common-endpoint-seam-2026-07-27-full-x-annihilation
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
related:
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2401-common-filter-endpoint-or-first-death-certificate
script: 04-computation/lrc14_full_x_transition_annihilation_thm2568.py
output: 05-knowledge/results/lrc14_full_x_transition_annihilation_thm2568.out
script_sha256: 3d29abe7c695281e97f68a1a6c6118df682ef1b1e4766c4a53da90cb555e9cb1
output_sha256: 82e667b072d4d36263c365f742eeabd24f719840a6779b5508060ce6f37d854a
hash_basis: working-tree bytes (LF)
---

# THM-2568 -- the full-X transition is exactly in the coarse-target kernel

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2563 inserts both paired blocker--graft dipoles on its moving endpoint,
but its old dangerous head is fixed.  The tempting completion is to co-shift
that head and then use the resulting left-minus-right colour.  At the full
physical-frequency level this cannot work: the old head has `k_a` dangerous
and the repaired endpoint has the same role safe, so a lawful common co-shift
makes the two endpoint layers pointwise orthogonal.

There are two distinct Fourier objects:

```text
independent left/right twists
  -> refined pair colours (a,b);

common lawful target twist
  -> coarse difference colour q=a+b.                         (1)
```

The map in (1) is a line sum.  The transition lies in its kernel, even when
the refined table has positive drift away from the line `a+b=0`.  This is the
endpoint version of the refined-to-coarse category boundary in MISTAKE-261.

## 1. Full-X transition annihilation

Put `p=13` and `zeta=exp(2 pi i/p)`.  Fix all finite colours other than one
target axis, and fix the deepest physical frequency offset `m`.  For every
common target twist `s in F_p`, let

```text
L_s,R_s in L^2(T)                                             (2)
```

be complete left and right endpoint layers.  Suppose one Boolean factor is
dangerous on the left and safe on the right.  Thus for some Boolean mask
`P_s`, co-shifted by the same lawful target action on both endpoint copies,

```text
L_s=P_s L_s,                  R_s=(1-P_s)R_s                 (3)
```

almost everywhere.  All other endpoint factors, a delayed word, and the
other target/deep finite twists may be absorbed in (2).

Use the physical Fourier convention

```text
f_hat(X)=integral_T f(x) exp(-2 pi i Xx) dx                  (4)
```

and define the fixed-frequency transition current

```text
J_X(s)=L_s_hat(X) conjugate(R_s_hat(X+m)).                   (5)
```

Cauchy--Schwarz makes the `X`-sum in (5) absolutely convergent.  Endpoint
Parseval gives

```text
sum_X J_X(s)
 =integral_T L_s(x)conjugate(R_s(x))exp(2 pi i m x) dx
 =0                                                        (6)
```

for every `s`, by (3).  Taking the normalized target transform

```text
Jhat_X(q)=1/p sum_s J_X(s)zeta^(q s)                        (7)
```

and interchanging the finite `s`-sum with the absolutely convergent
`X`-sum yields the exact all-character law

```text
sum_X Jhat_X(q)=0                         for every q in F_p. (8)
```

Equation (8) is not a zero-mode statement.  It annihilates every coarse
target character, pointwise in all other retained finite labels.  Individual
fixed-`X` currents in (5) may be nonzero; they must cancel in their complete
physical-frequency recombination.

This is precisely the THM-2452 off-diagonal-mask cancellation specialized to
the THM-2563 danger-to-safe transition and followed through the target DFT.

## 2. The refined-pair spectrum and its coarse pushforward

To expose the lost coordinate, shift the two endpoint copies independently
and put

```text
A(s,t)
 =sum_X L_s_hat(X)conjugate(R_t_hat(X+m))

 =integral_T L_s(x)conjugate(R_t(x))exp(2 pi i m x) dx.      (9)
```

Define the normalized two-dimensional transform

```text
Ahat(a,b)
 =1/p^2 sum_(s,t)A(s,t)zeta^(a s+b t).                     (10)
```

With the THM-2334 sign convention, `a` is the left dipole residue and `b`
is minus the right dipole residue.  The actual target is therefore

```text
q=a+b=eta.(u-v).                                           (11)
```

The common-twist table is the diagonal `C(s)=A(s,s)`.  Finite Fourier
inversion gives the exact pushforward

```text
Chat(q)
 =sum_(a+b=q)Ahat(a,b)
 =1/p sum_s A(s,s)zeta^(q s).                              (12)
```

Equation (3) makes `A(s,s)=0`.  Hence

```text
sum_(a+b=q)Ahat(a,b)=0                    for every q.       (13)
```

The refined pair space has dimension `169`, the coarse target has dimension
`13`, and the line-sum kernel has dimension

```text
169-13=156.                                                (14)
```

Thus a nonzero coefficient `Ahat(a,b)` with `a+b!=0` is not a surviving
THM-2334 target coefficient.  Its whole line still sums to zero by (13).
Calling such a refined pair a coarse target would repeat exactly the genus of
MISTAKE-261: reconstructing entries before a forgotten-fibre sum does not
control the sum.

## 3. Diagonal-shift drift is only a refined diagnostic

On the refined table define the simultaneous-shift projection

```text
(P_diag A)(s,t)=1/p sum_(g in F_p)A(s+g,t+g).                (15)
```

It is the orthogonal projection onto functions of `t-s`.  Normalized finite
Parseval gives

```text
D_pair(A)
 :=1/p^2 sum_(s,t)|A(s,t)-(P_diag A)(s,t)|^2

 =sum_(a+b!=0)|Ahat(a,b)|^2.                               (16)
```

Consequently

```text
D_pair(A)=0
 iff A(s,t)=G(t-s) for some G.                              (17)
```

If `A` is a nonnegative transition table, diagonal zero further gives
`G(0)=0`.  A cheap sufficient test for positive refined drift is variation
of a normalized row marginal

```text
r(s)=1/p sum_t A(s,t).                                     (18)
```

Indeed

```text
1/p sum_s |r(s)-mean(r)|^2
 =sum_(a!=0)|Ahat(a,0)|^2
 <=D_pair(A).                                              (19)
```

The same holds for the column marginal.  Equations (16)--(19) are useful
diagnostics, but even `D_pair(A)>0` does not change (13).  Refined drift is
not a lawful coarse target current.

## 4. The sharp complementary-mask hostile

The failure is already exact on the finite target square.  Put

```text
A(s,t)=1_(s!=t).                                           (20)
```

This is the overlap table of one cyclic singleton mask against its safe
complement under independent shifts.  It is nonnegative, has positive mass
off the diagonal, and satisfies (17) with

```text
G(d)=1_(d!=0).                                             (21)
```

Its normalized spectrum is

```text
Ahat(0,0)=(p-1)/p,

Ahat(a,-a)=-1/p                    for a!=0,

Ahat(a,b)=0                        when a+b!=0.              (22)
```

Thus its refined drift is zero and its coarse target is zero.  Nevertheless
the one-sided slice

```text
A(0,t)=1_(t!=0)                                            (23)
```

has Fourier coefficient `-1/p` in every nonzero colour.  This is the
two-endpoint completion of THM-2563's one-sided hostile: all twelve moving
endpoint colours can be present while every left-minus-right target
vanishes.

There are also diagonal-zero tables with `D_pair>0`; Section 7's exact
control records one.  Their refined coefficients off `a+b=0` cancel line by
line under (13), so positive refined drift is no escape.

## 5. A positive anchored duplicate-probe mode, and its first type failure

There is a tempting strengthening of THM-2563 which is mathematically true
but remains on the wrong carrier.  Write

```text
d_L(y)=1_(||y||<L/14),                u_L=1-d_L,

L in {1,2}.                                               (24)
```

Let `w>=0` have mass `rho>0` and satisfy

```text
support(w) subset {u_1(a x)=1} intersection {d_L(k x)=1}.  (25)
```

For displayed left-danger shift `r` and right-repair shift `s`, define

```text
A^#(r,s)
 =integral_T w(x)
    u_1(a x-r/p)d_L(k x+r/p)
    u_1(a x-s/p)u_L(k x+s/p) dx.                           (26)
```

The two displayed copies give

```text
A^#(r,r)=0                         for every r,              (27)

A^#(r,0)=0                         for every r.              (28)
```

At `r=0`, the first pair is redundant on (25).  THM-2379's translated-tooth
count and inclusion--exclusion give

```text
sum_s u_1(a x-s/p)u_L(k x+s/p)>=11-2L.                     (29)
```

Therefore

```text
sum_s A^#(0,s)>=(11-2L)rho.                                (30)
```

Choose `delta!=0` such that

```text
A^#(0,delta)>=(11-2L)rho/p                                 (31)
```

and take the fixed-relative-offset section

```text
F_delta(r)=A^#(r,r+delta).                                 (32)
```

Equations (28), (31) give

```text
F_delta(0)>=(11-2L)rho/p,

F_delta(-delta)=0.                                        (33)
```

For

```text
Fhat_delta(q)=1/p sum_r F_delta(r)zeta^(q r),               (34)
```

Fourier inversion at the forced zero shows

```text
sum_(q!=0)zeta^(q delta)Fhat_delta(q)=-Fhat_delta(0).       (35)
```

Since `F_delta>=0`, equations (31)--(35) force some `q!=0` with

```text
Re[-zeta^(q delta)Fhat_delta(q)]
 >=(11-2L)rho/[p^2(p-1)].                                  (36)
```

At `p=13`, this is

```text
ordinary k, L=1:             9rho/2028=3rho/676,

guard k, L=2:                7rho/2028.                     (37)
```

The mechanism is the proposed zero-column/positive-row obstruction: (28)
and (30) prove that `A^#` cannot have the relative-twist form (17), and the
forced-zero section makes one refined colour quantitative.

The first invalid implication is equally exact.  In the THM-2559
application, `w` is the fixed target-informed old head.  It already contains
the unshifted danger/safe facts in (25), and its singleton selector and slope
stratum were defined using that unshifted failure mask.  Equation (26)
multiplies `w` by an additional shifted left pair.  It therefore has the
form

```text
fixed old selector * unshifted gates * shifted duplicate gates,          (38)
```

not the THM-2365 present endpoint `F_s` in which every target-active factor
is replaced by its co-shifted copy.  No factorization

```text
w=w_0 u_1(a x)d_L(k x)
```

with `w_0` target-neutral or covariantly transported is proved.  Reading
(36) as the hidden left residue would repeat MISTAKE-266's frozen-moving-gate
error.  Moreover `delta!=0` in (32) is a relative target normal; the
canonical common endpoint is `delta=0`, where (27) makes the table zero.

Thus (36)--(37) are a proved positive **duplicate-probe / relative-offset
boundary**, not a THM-2334 target current or a completion of THM-2563.

## 6. Exact consequence for the live seam

The valid chain is

```text
THM-2563 one-sided moving-endpoint colour
  + independently displayed left duplicate probe
  -> quantitative refined relative-offset colour;                    (39)

lawfully co-shift the actual danger head and safe repair
  + full-X recombination
  -> zero in every coarse target character.                           (40)
```

Therefore the missing-left-residue seam cannot be closed by another
nonnegative full-X common-endpoint overlap.  The next object must stop before
the annihilating operation.  Three faithfully typed possibilities are:

1. construct a lawfully co-shifted complete left atom and retain one exact
   physical frequency `X` in (5) before summing;
2. retain a proved normal/jet weight `psi(X)` and show that
   `sum_X psi(X)Jhat_X(q)` is a physical current; or
3. supply an oriented spanning reference of the kind classified
   algebraically by THM-2383, then prove that the reference is carried by the
   same physical endpoint cospan.

All three require a covariant old-head/left-atom construction.  The current
target-informed selector has no such factorwise orbit.  A fixed-frequency or
jet coefficient must also retain the word, owner, deep leg, and future-root
typing before it can feed a semantic arrival theorem.

This result closes the `full-X completion` branch negatively and replaces it
with the typed fixed-`X`/normal problem.  It excludes no scalar row; the live
ledger remains `165`, and LRC(14) remains open.

## 7. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_full_x_transition_annihilation_thm2568.py
python3 -O 04-computation/lrc14_full_x_transition_annihilation_thm2568.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_full_x_transition_annihilation_thm2568.out.
```

The exact companion checks:

- all `13^3=2,197` basis instances of the line-sum identity (12), and the
  `169/13/156` refined/coarse/kernel dimensions;
- every coefficient of the complement hostile (20)--(23);
- a positive refined-drift diagonal-zero control whose coarse line sums still
  vanish;
- four nonzero fixed-frequency complement terms with exact zero full-X sum;
- all `69` ordinary and `161` guard anchored translated-tooth profile pairs,
  including the zero diagonal, zero column, sharp capacities `9/7`, and every
  phase identity behind (35); and
- both rational floors in (37).

The endpoint Parseval argument, absolute convergence, THM-2334 line-sum
typing, and MISTAKE-266 boundary are symbolic proofs above, not finite
extrapolations.

The independent root audit rederived the endpoint Parseval/all-`X` law, the
`a+b=q` refined-to-coarse line sum and its `156`-dimensional kernel, the
diagonal-shift projection, the complement hostile spectrum, and the anchored
duplicate-probe floor.  It also independently located the first invalid
implication at the frozen target-informed selector and confirmed that the
MISTAKE-266 rejection, rather than any coarse-current conclusion, is the
correct typed boundary.

**QED.**
