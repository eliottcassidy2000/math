---
id: THM-2543
title: "Augmentation--norm relative-phase local-system dichotomy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; independent hostile audit pending.
  Every live THM-2517 response row has a branch-free horizontal boundary
  local system.  If the seven lawful future-phase means are nonconstant, a
  four-time tensor with a literal Boolean atom refinement has all 72
  gauge-invariant relative-phase/target characters and 216 signed positive-boundary
  incidences.  If the means are constant, THM-2539's cyclic norm supplies
  its 216 source/target/root incidences.  This removes the zero-phase versus
  full-support split only at the labelled horizontal level.  It does not
  select the literal phase-zero slice in the augmentation branch, identify
  visible and inherited ancestry coordinates, create semantic arrival, or
  remove an LRC(14) row.
source: codex-2026-07-27-holotopy-transfer-norm
depends_on:
  - THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2539-diagonal-cubic-owner-clock-boundary-current
  - THM-2540-weighted-live-event-kakeya-flux-and-transverse-gain-boundary-refinement
related:
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
script: 04-computation/lrc14_augmentation_norm_relative_phase_thm2543.py
output: 05-knowledge/results/lrc14_augmentation_norm_relative_phase_thm2543.out
script_sha256: b7171af0af89c9e73ac7d5718478f61f3d55147ac6a7fb4d3cb7ba6d56f1a87a
output_sha256: b309d7dd7caeae94014c374cf3cf79f7f55149abe668aeef595f06cb00e3bbd5
hash_basis: working-tree bytes (LF)
---

# THM-2543 -- augmentation versus norm on the relative-phase local system

THM-2517 splits the seven lawful future owner--word phases into a zero phase
and a full-support phase branch.  The split is needed for its two physical
constructions, but it is not the intrinsic representation-theoretic split.
Over the rationals,

```text
Q[C_7]=Q direct-sum Q(zeta_7).                               (1)
```

The first summand is the constant phase profile and the second is the whole
six-dimensional augmentation representation.  This produces an exhaustive
horizontal dichotomy:

```text
nonconstant phase means -> opposite row/phase characters survive;
constant positive means -> the cyclic multiplicative norm survives.       (2)
```

In the first lane the opposite characters are invariant under simultaneous
translation of the visible response row and future phase.  In the second
lane THM-2539 supplies its gauge-safe row-to-clock scheduler.  Both lanes
admit three Galois-uniform physical-root slopes and hence 216 atomwise
positive-boundary incidences.

This is a holotopy result only in the horizontal sense.  The first lane
retains a lawfully shifted future phase but does not choose the literal
phase-zero owner; the second retains that owner at a row-dependent scheduler
clock.  Neither construction supplies the directed semantic/ancestry arrow
which THM-2542 isolates.

## 1. The four-time phase tensor

Fix one live THM-2449/2517 response row.  Use THM-2517's notation

```text
F_(ell,s)(x),                    A_(ell,s)=integral F_(ell,s),

G_gamma(x) rational Boolean,    q_gamma=integral G_gamma in Q_(>=0),

q_0>0,                  ell,gamma in F_7, s in F_13.         (3)
```

The full shift from `G_0` to `G_gamma` includes the mandatory terminal-word
skew.  Put

```text
C_(ell,s)=A_(ell,s)^3.                                     (4)
```

By THM-2517,

```text
Chat(kappa,b)!=0

for every kappa in F_7^*, b in F_13^*.                     (5)
```

For a sufficiently large even integer `L`, let `N=13^L` and define

```text
R^L_(ell,s)(x)
 =F_(ell,s)(x)F_(ell,s)(Nx)F_(ell,s)(N^2x),

V^L_(ell,gamma,s)(x)
 =R^L_(ell,s)(x)G_gamma(N^3x),

K^L_(ell,gamma,s)=integral V^L_(ell,gamma,s).               (6)
```

Higher-order mixing, with the same separated-frequency proof as THM-2517
(O9a)--(O9b), gives uniformly on the finite table

```text
K^L_(ell,gamma,s) -> q_gamma C_(ell,s).                     (7)
```

Indeed, after a common finite-band trigonometric approximation, a surviving
nonconstant term would give

```text
n_0+n_1N+n_2N^2+n_3N^3=0
```

with all `|n_i|` below the bandwidth.  The largest-index term strictly
dominates for sufficiently large `N`.  Telescoping the bounded-product error
and passing through the `L^1` approximations proves (7).  Since
`13=-1 mod 7`, even `L` gives `13^(jL)=1 mod 7` for every integer `j`; all
four epochs preserve the lawful septimal phase convention.

Use normalized transforms

```text
Khat^L(kappa,eta,b)
 =1/637 sum_(ell,gamma,s) K^L_(ell,gamma,s)
      zeta_7^(kappa ell+eta gamma)zeta_13^(bs),

Chat(kappa,b)
 =1/91 sum_(ell,s)C_(ell,s)zeta_7^(kappa ell)zeta_13^(bs),

qhat(eta)=1/7 sum_gamma q_gamma zeta_7^(eta gamma).          (8)
```

Then (7) is exactly the finite-tensor factorization

```text
Khat^L(kappa,eta,b) -> Chat(kappa,b)qhat(eta).               (9)
```

The normalization `637=7*7*13` is load-bearing; there is no extra factor
seven in (9).

## 2. The all-or-flat cyclotomic lemma

For a rational profile `q=(q_0,...,q_6)`, put

```text
Q_q(X)=sum_(gamma=0)^6 q_gamma X^gamma.                     (10)
```

If `qhat(eta)=0` for one `eta!=0`, then `Q_q` vanishes at a primitive
seventh root.  Its minimal polynomial over `Q` is

```text
Phi_7(X)=1+X+...+X^6.                                      (11)
```

Since both polynomials have degree at most six, (11) divides (10) only when
all seven coefficients of `q` are equal.  Conversely a constant profile has
every nontrivial Fourier coefficient zero.  Therefore

```text
q nonconstant iff qhat(eta)!=0 for every eta in F_7^*.      (12)
```

Use the lawful source-action convention

```text
(a.K)_(ell,gamma,s)=K_(ell-a,gamma-a,s).                   (13)
```

Its transform is multiplied by `zeta_7^(a(kappa+eta))`.  A character
`(kappa,eta)` is invariant under (13) exactly when

```text
kappa+eta=0.                                                (14)
```

Thus, if `q` is nonconstant, (5), (9), and (12) imply that for every
sufficiently large even `L`, simultaneously,

```text
Khat^L(kappa,-kappa,b)!=0

for all kappa in F_7^*, b in F_13^*.                       (15)
```

These are `72` gauge-invariant **relative row/future-phase** characters.
They are not absolute source residues and do not identify `ell` with the old
deep or inverse-ancestry coordinate.

There is a nonnegative quotient table behind the character notation.  Put

```text
B^L_(delta,s)=1/7 sum_gamma K^L_(delta+gamma,gamma,s).       (15a)
```

Its normalized `(delta,s)` transform is exactly
`Khat^L(kappa,-kappa,b)`.  Thus (15) lives on the diagonal-action quotient
and remembers `delta=ell-gamma`; it forgets both absolute labels.  In
particular, calling `kappa` an absolute source residue would reverse the
meaning of the construction.

## 3. Three physical-root slopes in the augmentation lane

Keep the same `L` in the nonconstant branch and define the function channel

```text
Z^L_(kappa,b)(x)
 =1/637 sum_(ell,gamma,s)V^L_(ell,gamma,s)(x)
   zeta_7^(kappa(ell-gamma))zeta_13^(bs).                   (16)
```

Equation (15) says `integral Z^L_(kappa,b)!=0`.  On the predecessor-root
disintegration `x=(z+u)/13`, every factor in (6) except the first response
factor is exactly root-invariant: its argument is multiplied by a positive
power of thirteen.  The first response factor obeys THM-2533's common guard
support.  Hence `Z^L` vanishes on three or four cyclically consecutive roots
on every generic predecessor fibre.

Let `E_a` be the physical predecessor-root projector and set

```text
mathcal A_(kappa,b)
 ={a in F_13:E_a Z^L_(kappa,b)!=0}.                         (17)
```

The integral gives `0 in mathcal A`.  If `#mathcal A<=3`, evaluate Fourier
reconstruction on the same number of consecutive forbidden roots.  The
square minor

```text
(zeta_13^(au))_(u forbidden,a in mathcal A)                 (18)
```

is a consecutive Vandermonde and is invertible.  Every active projector
would vanish, contradicting (15).  Thus

```text
0 in mathcal A_(kappa,b),          #mathcal A_(kappa,b)>=4. (19)
```

Choose a seed `(kappa_0,b_0)` and three distinct nonzero active roots
`a_1,a_2,a_3`; put `lambda_i=a_i/b_0`.  All functions in (6) are rational
step functions.  Independent cyclotomic Galois transport acts by

```text
(kappa,-kappa,b,a)
 ->(v_7kappa,-v_7kappa,v_13b,v_13a).                       (20)
```

It is transitive on the `72` relative channels and preserves the root/target
slopes.  Consequently

```text
E_(lambda_i b)Z^L_(kappa,b)!=0

for every kappa,b!=0 and i=1,2,3.                          (21)
```

The same three slopes work for all `72` channels of one fixed row and delay.
They need not be row-independent or stable as the delay changes.
The rationality of the phase bank and common guard support are both
load-bearing.  For example, a small positive profile
`1+epsilon cos(4 pi gamma/7)` is real and nonconstant but has a selected
primitive Fourier mode equal to zero; without the common consecutive guard
zeros, (18) supplies no root-spread conclusion.

## 4. Atomwise positive boundaries give 216 incidences

The character sum (16) is not Boolean.  Expand each of the three response
sums before transforming.  Its labelled atoms are

```text
v^L_(ell,gamma,r_0,r_1,r_2,s,t)(x)
 =h_(ell,r_0,s,t)(x)
  h_(ell,r_1,s,0)(Nx)
  h_(ell,r_2,s,0)(N^2x)
  G_gamma(N^3x).                                            (22)
```

Every (22) is a literal Boolean intersection and keeps the inherited deep
diagonal

```text
v^L_(ell,gamma,t,r_1,r_2,s,t)=0.                           (23)
```

On a predecessor fibre write an atom as `e_u(z)g(z)`, where `g` is the
root-independent product of the last three factors.  Apply THM-2537 **before
complex character contraction**:

```text
partial^+_tau v=g e_u(1-e_(u+tau)),

partial^-_tau v=g(1-e_u)e_(u+tau),

partial^-_tau v-partial^+_tau v=(P_tau-I)v.                 (24)
```

Sum the labelled atoms, then take the relative-phase/target transform.  For
`tau=lambda_i` and root colour `a=lambda_i b`, (21) and the nonzero multiplier

```text
zeta_13^(a tau)-1=zeta_13^(lambda_i^2b)-1                  (25)
```

show that the signed boundary current is nonzero.  Hence the nonconstant
branch has

```text
3*6*12=216                                                  (26)
```

nonzero relative-phase/target/root incidences.  They are contractions of an
ordered pair of positive Boolean packet systems.  At least one orientation
contributes to each incidence, but one orientation need not work uniformly.

Every labelled atom retains its lawful future phase `gamma` and terminal-word
skew.  Only the `gamma=0` atoms contain the literal phase-zero future owner.
Nonvanishing of the sum over `gamma` does not select those atoms.

## 5. The constant branch is exactly the norm lane

If `q` is constant, then `q_gamma=q_0>0` for every `gamma`; (12) makes all
relative augmentation characters vanish.  But full support holds, and the
multiplicative cyclic norm is positive:

```text
P=product_gamma q_gamma=q_0^7>0.                            (27)
```

THM-2539 applies without alteration.  Its seven separated Latin assignment
cells use every phase exactly once, its diagonal cell pairs visible row
`ell` with scheduler clock `c=ell`, and the literal phase-zero owner occurs at
slot `d=-ell`.  It gives all `72` visible source/target means, three
Galois-uniform physical-root slopes, and `216` atomwise positive-boundary
incidences.

Thus (12) and (27) are exhaustive:

```text
q nonconstant -> augmentation/relative-phase lane (15)--(26);

q constant    -> cyclic-norm scheduler lane of THM-2539.    (28)
```

If `q` is nonconstant but has full support, both lanes are available.  The
dichotomy chooses the cheaper augmentation lane; it does not assert that the
two labelled systems have a common root slope or common Boolean atom.

## 6. The literal owner slice and the missing 2-cell

There is an important simultaneous but non-equivariant control.  For every
phase profile, including profiles with zeros,

```text
K^L_(ell,0,s) -> q_0 C_(ell,s).                             (29)
```

Therefore the fixed `gamma=0` slice eventually has all `72` visible
source/target characters and contains the actual phase-zero owner at one
common future epoch.  But source translation sends that slice to
`gamma=a`, not back to itself.  It is an anchored local chart, not a section
of the diagonal quotient (13).

The augmentation current (15) is gauge invariant but sums all shifted
phases.  The fixed slice (29) has the desired owner label but is not gauge
invariant.  Current canon supplies no common root/atom selector identifying
their nonzero coefficients.  In the language of THM-2542, that selector is
precisely a candidate vertical 2-cell:

```text
anchored phase-zero chart
  <---- missing common selector/intertwiner ---->
diagonal relative-phase local system.                       (30)
```

Neither factorization (9), cyclotomic irreducibility, nor Galois transitivity
constructs (30).

## 7. Covariant endpoint transport is horizontal, not arrival

THM-2540 sharpens the boundary scope.  On each positive atom of (24), let
`K_tau=e(1-P_tau e)` be the occupied tail and

```text
H_tau=P_(-tau)K_tau
```

its translated empty head.  For every nonnegative, possibly root-visible
sidecar density `w`, it proves

```text
integral w K_tau=integral (P_(-tau)w)H_tau.                 (31)
```

Thus the lawful phase/word sidecar can be transported with the root edge
without being falsely declared invariant.  Equation (31) gives an exact
covariant edge-transport law on the labelled root packet.  It does not restore the
weighted Cayley skew pairing, turn the empty head into a later occupied
event, or identify this row/phase torsor with THM-2542's scheduler/root C91
mapping torus.

## 8. Sharp controls and exact stopping boundary

Three finite controls show why both lanes and the missing selector are real.

1. The abstract rational phase-profile control `q=(1,0,0,0,0,0,0)` has every
   nontrivial `qhat(eta)` nonzero, so the augmentation lane works, while the
   norm product is zero.  It is not asserted to be a separately realized live
   scalar-cover row.
2. `q=(1,1,1,1,1,1,1)` kills every augmentation character exactly, while
   its cyclic norm is positive.
3. The phase-zero slice of the first profile is nonzero, but translating the
   phase chart moves its support off `gamma=0`; literal owner selection is not
   a gauge-invariant operation.
4. Tensoring THM-2539's two-atom orientation hostile with any nonconstant
   rational phase profile shows that either fixed positive orientation,
   including the canonical selected-source orientation, may erase every
   relative mixed mean.  The ordered pair in (24) is sharp.

Finally, attach an auxiliary semantic-arrival indicator equal to zero on
every atom.  Equations (3)--(31), all `216` incidences, and all horizontal
endpoint transports remain unchanged.  This is the exact local-system hostile
to inferring a directed arrival from the branch-free boundary atlas.  It is a
finite typed hostile, not a claimed fourteen-speed scalar cover.

The theorem therefore closes

```text
arbitrary lawful future-phase profile
  -> augmentation or norm
  -> branch-free 216-incidence horizontal boundary system. (32)
```

It does not prove one Boolean character-contracted event, a common slope or
atom across the two lanes, a literal phase-zero owner inside the augmentation
coefficient, equality of visible source and inherited deep/ancestry labels, a
target-active semantic head, chronological arrival, scalar row exclusion, or
LRC(14).

## 9. Exact referee

Run

```bash
python3 04-computation/lrc14_augmentation_norm_relative_phase_thm2543.py
python3 -O 04-computation/lrc14_augmentation_norm_relative_phase_thm2543.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_augmentation_norm_relative_phase_thm2543.out
```

The dependency-free exact referee works in `F_547`.  It checks all `72`
cubic mixed modes; the all-or-flat lemma on `1,458` rational phase profiles;
`104,832` nonconstant relative diagonals and `144` flat zeros; `4,608` direct
normalized tensor factorizations and `4,608` normalized quotient-table checks;
`1,512` simultaneous-translation covariance cases;
the fixed-owner slice and a noncovariant control; all `4,901` consecutive
guard Vandermonde minors; a one-hot augmentation-positive/norm-zero bank with
all `216` positive-boundary incidences; `96` noninvariant endpoint-transport
checks; the flat norm-positive hostile; and an empty semantic vertical bank.

**QED (candidate; independent audit pending).**
