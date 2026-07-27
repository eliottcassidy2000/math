---
id: THM-2512
title: "Lawful interaction cut-bundle transplant and replica dichotomy"
status: >
  PROVED.  Transposing the doubly centred ANOVA interaction of any rational
  lawful THM-2449 owner table produces a row-zero F_13-by-F_7 defect to which
  the THM-2508 affine cut transform applies directly.  The anchored
  delta-plus-six-replicas branch is exactly the zero defect.  Off that branch,
  every one of the 5,184 primitive cut coefficients is nonzero and at least
  294 of the 504 toothpick components are nonzero; the construction is
  covariant under affine relabelling of both residue charts.  On each exact
  clock class it inherits the M+C/R law, so either the bundle is identically
  zero or at most one positive clock is exceptional.  This moves the static
  cut-bundle algebra onto a lawful live response table, but the centred table
  is signed, its diagonal sums are not Boolean events, the non-replica branch
  is not forced, and no common owner/arrival/deep ancestry current or LRC(14)
  row exclusion follows.  Independent audit requested.
source: codex-2026-07-27-lawful-interaction-cut-transplant
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
related:
  - THM-2340-owner-word-anova-target-landing
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
---

# THM-2512 -- the affine cut bundle lives on the lawful interaction table

**PROVED.  Independent audit requested.**

THM-2508 first arose on the punctured-stalk defect of THM-2436, a branch
already emptied by the septimal argument.  THM-2449 supplies a different
`C_7 x C_13` object on the still-live owner continuation.  Its ANOVA
interaction has exactly the two zero-margin laws needed to receive the cut
transform.  The bridge is simply the transpose

```text
d(s,ell)=I_(ell,s).                                             (1)
```

This gives a genuine algebraic transplant, but not yet a physical current:
`I` is signed and its toothpick diagonals need not be event masses.

## 1. From a lawful table to a doubly reduced defect

Let

```text
A=(A_(ell,s))_(ell in F_7,s in F_13) in Q^(7 x 13)              (2)
```

be a rational response table.  Write, as in THM-2449,

```text
R_ell=1/13 sum_s A_(ell,s),
C_s  =1/7  sum_ell A_(ell,s),
mu   =1/91 sum_(ell,s) A_(ell,s),

I_(ell,s)=A_(ell,s)-R_ell-C_s+mu.                              (3)
```

Define `d_A:F_13 x F_7 -> Q` by (1).  Both ANOVA margins vanish, so

```text
sum_ell d_A(s,ell)=0                 for every s,                (4)
sum_s   d_A(s,ell)=0                 for every ell.              (5)
```

Equation (4) is precisely the row-zero hypothesis of THM-2508.  Equation
(5) is extra information.  It kills the six-dimensional vertical kernel:
if `d_A(s,ell)=b(ell)` is independent of `s`, then (5) gives
`13b(ell)=0`, hence `d_A=0`.

Let `zeta=zeta_13`, `xi=zeta_7`, and use THM-2449's normalized transform

```text
Ahat(kappa,b)
 =1/91 sum_(ell,s) A_(ell,s) xi^(kappa ell)zeta^(b s).          (6)
```

For the unnormalized defect transform

```text
dtilde_A(alpha,gamma)
 =sum_(s,ell)d_A(s,ell)zeta^(-alpha s)xi^(-gamma ell),          (7)
```

the main effects in (3) vanish at mixed frequencies and therefore

```text
dtilde_A(alpha,gamma)=91 Ahat(-gamma,-alpha)

                         whenever alpha,gamma!=0.               (8)
```

This is the exact source-to-target dictionary; no dimension count or
approximation is being used.

## 2. The replica boundary is exactly the zero bundle

Assume the owner anchor

```text
A_(ell,0)=a 1_(ell=0),                    a>0.                   (9)
```

THM-2449 proves the following equivalence:

```text
I=0

iff there is w in Q^13, w_0=0, such that
    A_(0,s)=a+w_s,
    A_(ell,s)=w_s                         for ell!=0.             (10)
```

If `A` is nonnegative, then `w_s>=0`.  Thus (10) is exactly the sharp
delta-plus-six-replicas hostile, while

```text
d_A!=0                                                        (11)
```

is exactly the non-replica branch.

For `tau,alpha in F_13^*`, `a_0,beta in F_7^*`, and `c in F_7`, define
the THM-2508 cut bank

```text
R^A_(tau,a_0,c)(v)
 =sum_(ell in F_7)
   d_A(v-tau rep(a_0 ell+c),ell),                               (12)

Theta^A_(tau,a_0,c)(alpha)
 =sum_v R^A_(tau,a_0,c)(v)zeta^(-alpha v),                     (13)

Psi^A_(tau,a_0)(alpha,beta)
 =sum_c Theta^A_(tau,a_0,c)(alpha)xi^(-beta c).                (14)
```

The letter `a_0` labels a cut and is unrelated to the positive anchor mass
`a` in (9).  THM-2508's exact factorization and (8) give

```text
Psi^A_(tau,a_0)(alpha,beta)
 =K_(alpha tau,beta)dtilde_A(alpha,-beta a_0)

 =91 K_(alpha tau,beta)Ahat(beta a_0,-alpha),                  (15)

K_(u,beta)=sum_(j=0)^6(zeta^(-u)xi^(-beta))^j !=0.             (16)
```

Because `A` is rational, the coprime Galois group is transitive on its
`72` mixed characters.  Hence exactly one of the following happens:

```text
replica branch:
  d_A=0 and all primitive Psi coefficients vanish;

non-replica branch:
  Psi^A_(tau,a_0)(alpha,beta)!=0
  for every tau,a_0,alpha,beta!=0.                              (17)
```

The second line contains

```text
12*6*12*6=5,184                                                (18)
```

nonzero coefficients.  In particular, for any one prescribed primitive
quadruple,

```text
Psi^A_(tau,a_0)(alpha,beta)=0

iff I=0
iff A has the replica form (10).                               (19)
```

Thus a single primitive cut coefficient is a complete rational detector of
the replica boundary.  This statement is special to rational tables and to
the primitive mixed locus; it is false for arbitrary complex tables or at
`beta=0`.

There is also a pointwise component invoice.  If `d_A!=0`, equations
(5) and (11) show that it is not vertical.  THM-2508 therefore gives

```text
#{(tau,a_0,c):R^A_(tau,a_0,c)!=0}>=42*7=294.                   (20)
```

Every nonzero component has total sum zero, so it is a nonconstant rational
`13`-vector and all twelve of its nontrivial root colours are nonzero.
Consequently at least `294*12=3,528` entries in the raw `Theta` bank survive.
No uniform magnitude floor is claimed without a denominator or `L^1` bound
on the lawful table.

## 3. Affine chart covariance survives the transplant

Let

```text
g=(U,H;B,C) in AGL_1(F_13) x AGL_1(F_7)                        (21)
```

act on response tables by affine relabelling,

```text
(g.A)_(ell,s)
 =A_(B^(-1)(ell-C),U^(-1)(s-H)).                               (22)
```

Row means, column means, and the grand mean are merely permuted.  Therefore

```text
d_(g.A)(s,ell)
 =d_A(U^(-1)(s-H),B^(-1)(ell-C)),                              (23)
```

which is exactly the THM-2508 pullback action.  Its covariance laws become

```text
R^(g.A)_(tau,a_0,c)(v)
 =R^A_(U^(-1)tau,a_0 B,a_0 C+c)(U^(-1)(v-H)),                  (24)

Psi^(g.A)_(tau,a_0)(alpha,beta)
 =zeta^(-alpha H)xi^(beta a_0 C)
  Psi^A_(U^(-1)tau,a_0 B)(U alpha,beta).                       (25)
```

The anchor mark in (9) moves from `(0,0)` to `(C,H)`; after transporting
that mark, replica versus non-replica is invariant.  Equations (24)--(25)
are statements about residue-chart gauge.  They do not assert that every
affine relabelling is realized by one lawful temporal LRC operation.

This is the precise advance beyond using the `72` mixed coefficients alone:
the non-replica signal now sits in a finite covariant root-colour bundle,
rather than in a chart-dependent single toothpick projection.

## 4. Exact inheritance of the clock law

For the lawful THM-2449 table at clocks in one fixed word/parity class
`rho`, there are rational matrices `M_rho,C_rho` with

```text
A^R=M_rho+C_rho/R.                                             (26)
```

ANOVA centring, (12), and (14) are linear.  Hence

```text
d_(A^R)=d_(M_rho)+d_(C_rho)/R,                                (27)

Psi^(A^R)=Psi^(M_rho)+Psi^(C_rho)/R.                           (28)
```

There are two exact alternatives.

1. If `d_(M_rho)=d_(C_rho)=0`, every clock in the class is replica and the
   whole cut bundle vanishes.
2. Otherwise, at most one positive clock has `d_(A^R)=0`.  At every other
   sufficiently large lawful clock, all `5,184` primitive coefficients in
   (17) and the component invoice (20) hold.

Indeed, two distinct roots of (27) would force both coefficient matrices to
vanish.  This is the bundle form of THM-2449's finite double-additivity test.
The possible single exceptional clock is sharp already for a rectangle of
the form `1-13/R`.

## 5. What has and has not moved to the live branch

The proved chain is

```text
lawful rational owner response A^R
  -> doubly centred signed interaction I
  -> row-zero defect d_A(s,ell)=I_(ell,s)
  -> affine-covariant cut bundle
  -> either zero replicas or all 5,184 primitive cut coefficients.        (29)
```

Unlike the original THM-2508 application, (29) is attached to the lawful
THM-2449 response table rather than only to the already-empty high-septimal
punctured-stalk branch.  It closes the **algebraic transplant** of the cut
bundle to that table.

It does not close the LRC(14) proof, for four distinct reasons.

- THM-2449 does not force the non-replica alternative; the physically
  realizable uniform-offset hostile of THM-2456 remains compatible with the
  replica branch.
- ANOVA centring makes `I` signed.  A sum in (12) is a diagonal contraction
  of response masses, not the measure of one Boolean owner/arrival event.
- The output index `v` and cut phase `beta` have not been put on the same
  ancestry sheet as THM-2471's first-collision atoms or THM-2478's old deep
  probe.  THM-2449's deep witness may depend on the mixed colour.
- Affine chart covariance does not supply address/gain participation,
  terminal phase, owner-loop drift, or a nonzero integrated scalar current.

The cheapest next bridge test is therefore not another colour census.  It is
whether one signed toothpick contraction in (12), or a phase-neutral paired
version of it, can be realized on a common Boolean owner/deep ancestry fibre.
Until that happens, the exact LRC(14) ledger remains `165`. **QED.**
