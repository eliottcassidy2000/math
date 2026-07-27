---
id: THM-2530
title: "Anchored deep-Gram path cone and lossless skew target refinement"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In every lawful THM-2365 target cell, the Gram
  matrix of the thirteen translated deepest-comb indicators is positive
  semidefinite, symmetric, nonnegative, and has a zero row and column at the
  excluded target root.  Pointwise its comb mask is either one non-target
  root or two adjacent non-target roots, so the matrix lies in an exact
  23-ray anchored-path cone with lossless singleton/pair coordinates.  Every
  nontrivial anti-diagonal Fourier value is positive whenever the cell has
  positive mass, with a sharp cone-uniform lower bound.  Anticommuting this
  Gram matrix with any oriented cyclic Hilbert operator H_tau gives a skew
  ordered-root matrix whose twelve anti-diagonal modes are all nonzero.  The
  skewing is lossless, not merely detecting: the Cayley identity gives an
  explicit telescoping inverse from the anchored zero row.  Across target
  cells, the Gram drift is exactly thirteen times THM-2365 H-drift plus twice
  the adjacent-double-mask drift.  Thus the old circulant hostile either
  acquires a literal double-deep target mode or strengthens to invariance of
  the entire 23-coordinate mask law.  The pointwise object is a Boolean
  exterior-wedge lift, but singleton masks use unoccupied comparison roots;
  ties remain and no scalar row, uniform occupied pair, tournament, or
  LRC(14) proof follows.
source: codex-2026-07-27-anchored-deep-gram
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2526-affine-skew-orientation-gauge-boundary
related:
  - THM-2354-deep-shift-comb-cover-and-grouped-unit-current
  - THM-2515-haar-self-correlation-disintegration-and-rational-shift-recovery
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
  - THM-2529-deep-comb-adjacent-fibre-odd-consumer-and-zero-target-boundary
script: 04-computation/lrc14_anchored_deep_gram_skew_refinement_thm2530.py
output: 05-knowledge/results/lrc14_anchored_deep_gram_skew_refinement_thm2530.out
script_sha256: a5d4ff3a33e6f2b005033e3fc4de88a56f927bb1f746123311650c2a9fe5f192
output_sha256: 9a1e01ec07f39ef09d7fd96197deaf45f0da071de062339b2c2db8fd0a68eb98
hash_basis: working-tree bytes (LF)
---

# THM-2530 -- the excluded target root makes the deep Gram skewing lossless

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2365 leaves one exact residual law:

```text
H(r,s,t)=G(r-t).                                           (1)
```

That law sees one deepest-comb indicator at a time.  The next same-point
Boolean object is not another character sum but the full Gram matrix of all
thirteen translated deep indicators.  Its off-diagonal entries remember
whether two adjacent deep roots are simultaneously active.  The geometry is
much more rigid than an arbitrary positive-semidefinite matrix: after the
target root is removed, it is the moment cone of vertices and edges of a
twelve-vertex path.

There is also an exact reason a skew orientation becomes possible here.
THM-2526's cyclic Hilbert operator has a thirteen-dimensional Hankel kernel
when it anticommutes with an arbitrary matrix.  The forbidden target row is
zero, and that one anchored row intersects the Hankel kernel only at zero.
Thus the skew lift retains the whole Gram matrix.

## 1. Lawful anchored Gram matrix

Keep THM-2365's notation, put `p=13`, and fix one lawful target cell `(s,t)`.
Write

```text
F=F_(s,t),

delta_u(x)=Delta_(t+u)(x)
          =1_(||c x-(t+u)/13||<1/14),       u in F_13.     (2)
```

The deepest present-owner factor in `F` is the complement of `Delta_t`.
Consequently, up to the null endpoint set,

```text
0<=F<=1-Delta_t,
F delta_0=0.                                               (3)
```

Define

```text
K_(s,t)(u,v)
 =integral_T F_(s,t)(x) delta_u(x)delta_v(x) dx.            (4)
```

This is the Gram matrix of the real functions

```text
sqrt(F)delta_u.
```

Hence

```text
K>=0 entrywise,
K^*=K,
K is positive semidefinite,
K(0,v)=K(u,0)=0.                                           (5)
```

The last line is the **target anchor**.  It marks the target root as zero in
the relative root chart.  A single unlabelled matrix does not recover the
absolute residue `t`; the lawful family still carries its `(s,t)` cell label.

## 2. The exact 23-ray anchored-path cone

For a fixed `x`, let

```text
S_t(x)={u in F_13:delta_u(x)=1}.                            (6)
```

The thirteen centres are spaced by `1/13`, whereas the danger arc has length
`1/7`.  The nearest centre is at distance at most `1/26<1/14`, so the set in
(6) is nonempty almost everywhere.  Three centres cannot lie in one such arc
because

```text
1/7<2/13.
```

Thus `S_t(x)` is either one root or two cyclically adjacent roots.  On the
support of `F`, equation (3) removes root zero.  Using representatives
`1,...,12`, the only masks are therefore

```text
{j},                         1<=j<=12,
{j,j+1},                     1<=j<=11.                    (7)
```

Let their `F`-weighted masses be

```text
alpha_j=integral F 1_(S_t={j}),
beta_j =integral F 1_(S_t={j,j+1}).                        (8)
```

Then `alpha_j,beta_j>=0` and

```text
K
 =sum_(j=1)^12 alpha_j e_j e_j^T
  +sum_(j=1)^11 beta_j(e_j+e_(j+1))(e_j+e_(j+1))^T.        (9)
```

This representation is lossless.  With `beta_0=beta_12=0`,

```text
beta_j=K(j,j+1),

alpha_j=K(j,j)-K(j-1,j)-K(j,j+1).                         (10)
```

In particular, the Gram matrix is tridiagonal on the path

```text
1--2--...--12,
```

and (10) gives its exact cone inequalities.  If

```text
A=sum_j alpha_j,              B=sum_j beta_j,
```

then

```text
mu(F)=A+B,
tr K=A+2B,
1^T K 1=A+4B.                                            (11)
```

The endpoint irregularity is load-bearing: deleting root zero turns the
cyclic adjacency graph into a path and provides the boundary condition used
by the inverse below.

## 3. Every anti-diagonal Fourier mode is positive

Let `zeta=exp(2 pi i/13)` and use the unnormalised two-variable transform

```text
K_tilde(a,b)=sum_(u,v)K(u,v)zeta^(a u+b v).                (12)
```

Evaluating (9) on the anti-diagonal gives the exact two-parameter collapse

```text
K_tilde(a,-a)
 =A+(2+zeta^a+zeta^(-a))B
 =A+4 cos^2(pi a/13)B.                                    (13)
```

Thus, whenever `mu(F)>0`,

```text
K_tilde(a,-a)>0                    for every a in F_13.    (14)
```

More precisely,

```text
K_tilde(a,-a)
 >=min(1,4 cos^2(pi a/13))mu(F),                           (15)

K_tilde(a,-a)
 >=kappa_13 mu(F)                  for a!=0,

kappa_13=2-2 cos(pi/13)>0.                                 (16)
```

The last constant is sharp on the abstract cone: take a pure adjacent-pair
ray and `a=+/-6`.  Equation (13) is also a warning.  These positive modes see
only the total singleton and double masses `(A,B)`, not the locations of the
23 path rays.  Positivity alone is therefore not target drift.

## 4. The correct skew operation is the anticommutator

For `tau in F_13^*`, use THM-2526's pullback and cyclic Hilbert operator

```text
(T_tau f)_u=f_(u+tau),

H_tau=sum_(r=1)^6(T_(-2 tau r)-T_(2 tau r)).                (17)
```

They satisfy

```text
H_tau^*=-H_tau,
(I+T_tau)H_tau=T_tau-I.                                   (18)
```

The oriented slope `tau` is a retained sidecar.  In the live collision
carrier it may be chosen as THM-2526's guard-sheet slope
`tau_H=(-H)^(-1)`; in the target-coordinate chart, `tau=1` is the natural
positive co-shift direction.

The commutator of a skew and a symmetric matrix is symmetric.  The correct
skew operation is instead

```text
J_tau=H_tau K+K H_tau.                                    (19)
```

Indeed,

```text
J_tau^*=-J_tau.                                           (20)
```

If

```text
m_tau(a)=(zeta^(a tau)-1)/(zeta^(a tau)+1),                (21)
```

then the transform convention in (12) gives

```text
J_tilde_tau(a,b)
 =(m_tau(b)-m_tau(a))K_tilde(a,b).                         (22)
```

Since `m_tau(-a)=-m_tau(a)`, equations (13)--(14) imply

```text
J_tilde_tau(a,-a)
 =-2m_tau(a)K_tilde(a,-a)!=0             for every a!=0.  (23)
```

All twelve reflection-odd root modes survive.  A uniform quantitative form
is

```text
|J_tilde_tau(a,-a)|
 >=2 tan(pi/13)(2-2 cos(pi/13))mu(F)       for a!=0.       (24)
```

This uses only that multiplication by nonzero `tau` permutes the nonzero
residues.  It is deliberately a safe uniform bound, not the modewise optimum.

## 5. The skew lift is exactly invertible

Equation (23) proves nonvanishing but does not explain why the full matrix is
retained.  The sharper fact follows directly from the Cayley identity.  Put

```text
L_tau=(1/2)(I+T_tau)J_tau(I+T_tau).                        (25)
```

Using (18) on both sides of (19) gives

```text
L_tau=T_tau K T_tau-K.                                    (26)
```

In entries,

```text
L_tau(u,v)=K(u+tau,v-tau)-K(u,v).                          (27)
```

For each `u`, let `q_tau(u)` be the unique integer in `{0,...,12}` with

```text
u+q_tau(u)tau=0 mod 13.
```

Telescope (27) until its first coordinate reaches the anchored zero row.
Equation (5) yields the explicit inverse

```text
K(u,v)
 =-sum_(n=0)^(q_tau(u)-1)
      L_tau(u+n tau,v-n tau).                              (28)
```

Therefore

```text
K -> J_tau
```

is injective not only on the 23-ray Gram cone but on the entire
`156`-dimensional space of matrices with zero row.

The boundary is sharp.  Without the anchored row, (26) shows that the kernel
is exactly

```text
K(u,v)=g(u+v),                     g:F_13->C,               (29)
```

the thirteen-dimensional cyclic Hankel space.  Root zero evaluates every
value `g(v)`, so `K(0,*)=0` kills precisely this kernel.  This is why an
excluded target root can support a lossless ordered current even though an
unanchored affine skew sector cannot.

## 6. Pointwise Boolean exterior-wedge semantics

Let

```text
e(x)=(delta_u(x))_(u in F_13).                             (30)
```

Then

```text
K=integral F(x)e(x)e(x)^T dx,

J_tau
 =integral F(x)[(H_tau e(x))e(x)^T
                 -e(x)(H_tau e(x))^T] dx.                 (31)
```

Thus every point contributes a decomposable integral two-form.  Its positive
and negative coefficients can be separated on the finite rooted fibre
`(x,u,v)`, so (31) is a genuine signed Boolean exterior-wedge lift rather
than an analogy to one.

At target-coordinate slope `tau=1`, the two point-mask types are especially
transparent.

* If `e=e_j`, equation (31) is a complete unit oriented star centred at the
  unique occupied deep root `j`.
* If `e=e_j+e_(j+1)`, the Cayley telescope gives

  ```text
  H_1 e=e_j-e_(j+1),

  (H_1e)e^T-e(H_1e)^T
   =2(e_j e_(j+1)^T-e_(j+1)e_j^T).                       (32)
  ```

  This is a literal ordered pair of two simultaneously occupied adjacent
  deep roots.

For a general guard-selected `tau`, (31) remains lawful and lossless, but an
adjacent target-coordinate pair need not collapse to the single edge (32).
The comparison roots in `H_tau e` must then be retained.

There is no uniform two-occupied-root conclusion: all `beta_j` may vanish
while some `alpha_j` is positive, and equations (23)--(24) still hold.  Nor
is `J_tau` a tournament.  A pure adjacent-pair ray has only one nonzero
unordered edge, while a singleton ray has only a star; arbitrary mixtures
can create further ties and cancellations.  Thresholding signs would destroy
the amplitudes needed by (28).

## 7. Exact refinement of the THM-2365 zero branch

Restore all target cells and write `K_(s,t)` and `beta_(s,t),j`.  The old
tensor is exactly the diagonal of the new one:

```text
H(r,s,t)=K_(s,t)(r-t,r-t).                                (33)
```

Let

```text
K_bar=1/13^2 sum_(s,t)K_(s,t),

beta_bar_j=1/13^2 sum_(s,t)beta_(s,t),j,

E_K=1/13^2 sum_(s,t)||K_(s,t)-K_bar||_F^2.                 (34)
```

THM-2365's target-action projection fixes `u=r-t` and averages `(s,t)`.
Consequently its drift is

```text
D_H
 =1/13^3 sum_(s,t,u)
    |K_(s,t)(u,u)-K_bar(u,u)|^2.                           (35)
```

The path-cone decomposition has only diagonal and symmetric adjacent-edge
coordinates, so there is an exact orthogonal split:

```text
E_K
 =13D_H
  +2/13^2 sum_(s,t)sum_(j=1)^11
      |beta_(s,t),j-beta_bar_j|^2.                         (36)
```

Thus the stronger target-cell alternative is:

```text
D_H>0
  -> the existing THM-2365 nonzero-target/deep-colour landing;

D_H=0 but some beta_(s,t),j varies
  -> a nonzero target Fourier mode of the literal Boolean overlap
     integral F_(s,t) Delta_(t+j)Delta_(t+j+1);

D_H=0 and every beta_(s,t),j is constant
  -> every alpha_(s,t),j and beta_(s,t),j is constant,
     so the complete 23-coordinate anchored mask law is target-invariant.
                                                               (37)
```

The last implication uses (10): constant diagonals and constant `beta`
force constant `alpha`.  It is strictly stronger than the old residual law
(1), which controls only the twelve diagonal coordinates.

Because (19) and (28) are linear cell by cell, with

```text
J_bar_tau=H_tau K_bar+K_bar H_tau,
```

we also have

```text
J_(s,t),tau=J_bar_tau for every (s,t)
 iff K_(s,t)=K_bar for every (s,t).                         (38)
```

Thus the skew family loses none of the old `H`-drift and none of the new
adjacent-pair drift.

The new coefficient in the middle branch of (37) is a lawful same-point
higher moment: both deep indicators are honest Boolean factors and their
intersection has length `1/7-1/13=6/91` before other factors are imposed.
It is not automatically a coefficient of THM-2334's one-deep scalar row.
Turning that double-deep target mode, or the signed wedge (31), into one
source-to-arrival Boolean emission remains an additional semantic step.

THM-2529 gives the sharp contraction warning.  Its aggregate deep-mask
self-correlation retains only `(A,B)`, whereas (9) retains every location
`alpha_j,beta_j`.  More strongly, its consequence-level comb-compatible
hostile

```text
F_t=Delta_(t+1)d(13c x)                                  (38a)
```

has one fixed relative singleton coordinate, `alpha_1>0`, and `beta=0`.
Thus `K_(s,t)` and `J_(s,t),tau` are constant in the relative target chart:
all twelve root modes in (23) survive while every target-cell variation in
(36) vanishes.  This hostile is not asserted to be a canonical
`E_(s,t)W` owner packet.  It nevertheless proves at the exact comb-consequence
level that root all-mode support cannot substitute for target drift, even
before the Gram matrix is contracted.

## 8. Scope and exact companion

This theorem replaces the featureless circulant hostile by a recursive
boundary:

```text
one-deep target drift
  -> adjacent-double-deep target drift
  -> full anchored path-mask invariance.                   (39)
```

It also explains why the target anchor and the cyclic Hilbert orientation
fit: the former kills exactly the latter's Hankel anticommutator kernel.

It does **not** prove that the adjacent-pair drift in (37) is positive, force
a two-deep pair in a positive cell, identify the signed comparison roots with
a future owner/source edge, exclude a scalar speed row, or prove LRC(14).
The all-mode statement (23) fires even on the singleton hostile and must not
be promoted past those boundaries.

The dependency-free companion checks, with exact integer/rational arithmetic:

* all `23` point masks and their lossless path-cone coordinates;
* the anti-diagonal group-algebra identity;
* the singleton-star and adjacent occupied-edge wedges at `tau=1`;
* the Cayley and skew identities for every nonzero slope;
* the inverse (28) on all `12*156=1,872` zero-row basis cases;
* a cone-compatible zero-`H`-drift control with positive adjacent-pair drift;
  and
* the exact decomposition (36) on an independent nonuniform control.

Run

```bash
python3 04-computation/lrc14_anchored_deep_gram_skew_refinement_thm2530.py
python3 -O 04-computation/lrc14_anchored_deep_gram_skew_refinement_thm2530.py
```

Both executions agree exactly. **QED.**

An independent audit rederived the `23`-ray cone and inverse coordinates,
the sharp `kappa_13`, the two-variable Fourier signs, the Cayley sandwich and
Hankel kernel, the zero-row telescope, both pointwise wedge types, and the
normalisation in (36).  It also reproduced normal and optimized output
byte-for-byte.
