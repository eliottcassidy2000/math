---
id: THM-2515
title: "Haar self-correlation disintegration and rational-shift recovery"
status: >
  PROVED + INDEPENDENTLY AUDITED.  For any finite rational-step response bank,
  its entrywise square
  moment is the exact average of its same-circle translation-autocorrelation
  tables.  Hence every nonzero linear coefficient of the square moment is
  already nonzero on one rational translation fibre, with no magnitude loss
  for that selected coefficient.
  On a rational C_7-by-C_13 table, one surviving mixed coefficient on that
  fibre Galois-saturates all 72 mixed colours and, after ANOVA centring, all
  5,184 THM-2508 primitive cut coefficients.  The product pair space therefore
  disintegrates into a finite difference-marked family, one fibre of which
  carries the signal.  Autocorrelation is antipodally even;
  the selected shift need not preserve a positive marked anchor, need not be
  a lawful 91-root/temporal shift, and still does not give a one-point Boolean
  owner/deep current or LRC(14).
source: codex-2026-07-27-self-correlation-disintegration
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
related:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2511-affine-cut-quadratic-root-service-and-pair-space-boundary
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - THM-2513-anchored-first-or-second-moment-spectrum-and-pair-space-boundary
---

# THM-2515 -- a square moment is an average of translation cospans

**PROVED + INDEPENDENTLY AUDITED.**

The entrywise square of a table of masses appears to live on two independent
copies of the underlying probability space.  On a compact group this pair
space has a canonical difference coordinate.  Disintegrating by that
difference rewrites the square as an average of one-parameter
self-correlations.

For rational step functions the disintegration is finite and exact.  It is
the parallel-class decomposition

```text
C_D x C_D = disjoint union_(k in C_D){(x,x+k):x in C_D}.       (1)
```

Thus a second-moment spectral signal already lives on one rational diagonal
fibre.  The new datum is the difference `k`; forgetting it recreates the
pair-space average.

## 1. Exact finite correlation disintegration

Let `J` be a finite label set, let `D>=1`, and let

```text
F_j:C_D -> Q,                         j in J.                   (2)
```

Use normalized counting measure and define the first moments

```text
A_j=1/D sum_(x in C_D)F_j(x),                                (3)
```

and the translation autocorrelation tables

```text
C_k(j)=1/D sum_(x in C_D)F_j(x)F_j(x+k),
                                      k in C_D.                (4)
```

Then, entry by entry,

```text
1/D sum_(k in C_D)C_k(j)=A_j^2.                              (5)
```

Indeed, the map `(x,k)->(x,y=x+k)` is a bijection of `C_D^2`, so

```text
1/D sum_k C_k(j)
 =1/D^2 sum_(x,y)F_j(x)F_j(y)
 =(1/D sum_x F_j(x))^2.                                      (6)
```

Consequently, for every linear functional

```text
Lambda:Q^J -> K                                               (7)
```

into a subfield `K` of `C`, equipped with the usual absolute value,

```text
Lambda(A^(circ 2))=1/D sum_k Lambda(C_k).                     (8)
```

In particular,

```text
Lambda(A^(circ 2))!=0

  => there is k in C_D with Lambda(C_k)!=0,                   (9)

max_k |Lambda(C_k)|>=|Lambda(A^(circ 2))|.                   (10)
```

No positivity, genericity, or asymptotic mixing is used.  If every `F_j` is
nonnegative, every `C_k(j)` is nonnegative as well.

For real-valued `F_j`, substitution `y=x-k` gives the antipodal law

```text
C_(-k)=C_k.                                                   (11)
```

Thus the disintegration retains the difference fibre but forgets its
orientation.  This is the exact first loss of the construction.

## 2. Rational step functions on the circle

Let `F_j:T=R/Z -> Q` be step functions, constant on the cells

```text
[m/D,(m+1)/D),                         m in C_D,                (12)
```

for one common denominator `D`.  Identify each function with its `D` cell
values in (2).  Then

```text
A_j=integral_T F_j(x)dx,

C_k(j)=integral_T F_j(x)F_j(x+k/D)dx.                         (13)
```

Equations (5)--(11) hold literally.  In particular, the selected fibre in
(9) is a rational translation `t=k/D`, and its table entries are rational.

For arbitrary `L^2` functions on a compact abelian probability group, the
same identity has the Haar form

```text
(integral F_j)^2
 =integral_t integral_x F_j(x)F_j(x+t) dx dt.                 (14)
```

This follows from Haar invariance and Fubini.  The finite rational version is
stronger for the present application because it produces an exact common-grid
shift without a limiting argument.

There is also a useful base-`13` refinement.  For `L^2` functions the map

```text
t -> Lambda(C_t)                                             (14a)
```

is continuous, because translation is strongly continuous in `L^2`.  If its
Haar integral is nonzero, its nonzero set contains a nonempty open interval.
The grids `{k/13^L}` become dense, so for arbitrarily large `L` one may choose

```text
t=k/13^L                                                     (14b)
```

with `Lambda(C_t)!=0`.  Rational step functions still give a rational table
at such a shift.  This choice guarantees exact nonvanishing; the finite
common-grid choice is the one that directly gives the magnitude invoice
(10).

Geometrically, (1) is a one-direction finite needle foliation of the product
torus: every ordered pair lies on exactly one slope-one diagonal, labelled
by its difference.  This is a decomposition, not a Kakeya covering theorem;
no second direction or small-measure assertion is being imported.

## 3. Primitive-colour recovery on `C_7 x C_13`

Take

```text
J=F_7 x F_13,

A=(A_(ell,s))_(ell,s),                                       (15)
```

and suppose one primitive mixed coefficient of its entrywise square is
nonzero.  With `xi=zeta_7`, `zeta=zeta_13`, assume for some
`kappa,b!=0` that

```text
Lambda_(kappa,b)(A^(circ 2))
 =1/91 sum_(ell,s)A_(ell,s)^2
       xi^(kappa ell)zeta^(b s)
 !=0.                                                         (16)
```

Apply (9) to this `Lambda`.  There is a rational grid shift `t=k/D` for
which the nonnegative rational correlation table

```text
C^t_(ell,s)
 =integral_T F_(ell,s)(x)F_(ell,s)(x+t)dx                    (17)
```

satisfies

```text
Lambda_(kappa,b)(C^t)!=0.                                    (18)
```

Because all entries of `C^t` are rational, the coprime Galois group

```text
Gal(Q(zeta_91)/Q)
 =(Z/7Z)^* x (Z/13Z)^*                                      (19)
```

is transitive on the `72` primitive mixed characters.  Therefore (18)
implies

```text
Lambda_(kappa',b')(C^t)!=0

for every kappa' in F_7^*, b' in F_13^*.                    (20)
```

Centre `C^t` by its two-way ANOVA interaction and transpose it:

```text
d_t(s,ell)=C^t_(ell,s)
            -rowmean_ell(C^t)-colmean_s(C^t)+grandmean(C^t).
                                                                    (21)
```

Then `sum_ell d_t(s,ell)=0`.  Its primitive transform is `91` times
the corresponding mixed transform in (20), with the sign convention of
THM-2512.  THM-2508's exact cut factorization now gives

```text
Psi^t_(tau,a)(alpha,beta)
 =91 K_(alpha tau,beta)
   Lambda_(beta a,-alpha)(C^t)
 !=0                                                          (22)
```

for every

```text
tau,alpha in F_13^*,               a,beta in F_7^*.           (23)
```

Thus all

```text
12*12*6*6=5,184                                               (24)
```

primitive cut coefficients survive on one rational translation fibre.  The
chosen mixed coefficient also retains the magnitude bound (10).  No claim is
made that the same `t` optimizes, or gives a common lower bound for, all its
Galois conjugates.

## 4. Application type and sharp boundaries

At any fixed finite THM-2449 clock, each response entry is the integral of a
nonnegative rational step density: sum its finitely many nonnegative
integrands defining `A^R` before integration.  Therefore any square-moment mixed
coefficient supplied on that table satisfies the hypotheses of Section 3.
The finite family admits one common grid denominator after taking the least
common multiple of all rational endpoints.  This makes `k/D` algebraically
available; it does not make that shift a lawful physical clock.
The conclusion replaces an unstructured independent pair by

```text
one base point x,
one translated point x+t,
one retained rational difference t.                           (25)
```

That is a real reduction in sidecar complexity, but it is not descent to one
ordinary LRC event.

Three losses are sharp.

1. **Orientation.**  Equation (11) identifies `t` and `-t`.  An
   autocorrelation fibre cannot by itself produce oriented owner-loop drift.
2. **Anchor mismatch.**  Even if a separately marked anchor autocorrelation
   is positive somewhere, the shift selected by (9) need not lie there.  On
   `C_4`, let

   ```text
   H=(1,0,0,0),
   U=(1,0,0,0),
   V=(3/5,4/5,0,0).
   ```

   The marked autocorrelation of `H` is positive only at shift zero.  For the
   linear functional `Lambda(C)=C(U)-C(V)`, the shift-zero value is zero,
   while the values at shifts `1,2,3` are `-3/25,0,-3/25`; their average is
   `-3/50=(mean U)^2-(mean V)^2!=0`.  Hence every signal-carrying shift misses
   the marked anchor.
3. **Physical shift type.**  The rational difference `t=k/D` need not be a
   canonical `1/91` root, a THM-2471 first-collision translation, a lawful
   target/deep action, or a delay already present among the fourteen speeds.
   Choosing the alternative base-`13` form `k/13^L` does not by itself supply
   any of those semantic identifications.
   Multiplying the two densities in (17) is a two-point cospan observable,
   not a one-point Boolean current.

The theorem therefore disintegrates pair space into a finite
translation-marked local system and restores the full primitive cut spectrum
on one fibre.  An independent audit rederived the double-counting
normalization, the Galois and cut-transform signs, the `5,184` count, and the
`C_4` hostile, and checked the exact distinction between common-grid and
base-`13` shift selection.  The next
physical question is whether its retained difference can be identified with
an existing collision/ancestry edge while preserving a positive owner mark.
It does not supply that identification, remove any of the `165` live rows,
or prove LRC(14). **QED.**
