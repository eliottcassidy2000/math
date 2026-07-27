---
id: THM-2536
title: "Deep-comb start/terminal selector flow, target drift, and centered toothpick boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  On every THM-2530 target-anchored deep-comb
  cell, the positive start-selector masses gamma+ and terminal-selector
  masses gamma- are respectively alpha_j+beta_j and
  alpha_j+beta_(j-1).  Their difference is the discrete divergence of the
  adjacent-pair flow beta, and the anchored path endpoints make the two
  selector vectors jointly lossless for all 23 path-ray masses.  The
  opposite-selector divergence is nonzero exactly when an adjacent pair
  has positive mass and then has all twelve primitive root modes.  It is
  already zero-sum, but identifies the literal pair edge with a detour
  through the excluded target up to a target-triangle circulation.  Each
  positive rational selector vector has all twelve primitive root modes.
  Across the 169 lawful target cells, either selector has the same exact
  target-line cancellation and 1/13 nonzero-root-colour energy bound as
  THM-2365; its target drift vanishes exactly when that selector vector is
  cellwise constant.  The THM-2529 singleton hostile makes both selectors
  constant while retaining every root mode, so target variation is not
  uniform.  On the canonical typed non-cover word used by THM-2334's
  independently verified twist control, all 169 selector-cell totals are
  positive, take 75 exact rational values, and all 168 nontrivial target
  characters of that mass surface survive an exact cyclotomic certificate.
  This is a positive control, not a covering-row theorem.  The selectors
  are honest same-time one-deep-minus-adjacent-double-deep observables, not
  existing THM-2334 one-deep marked-current coefficients.  A THM-2508
  toothpick cannot realize the raw positive class vector, but it can realize
  either a uniform mean-centering or, more naturally, the zero-sum boundary
  of its genuine occupied-to-excluded-target star.  An explicit signed two-
  column lift proves this algebraically while the artificial C7 pair and
  divergence quotient still forget edge/triangle holonomy and do not create
  source/arrival semantics.
source: codex-2026-07-27-selector-target-current-bridge
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2530-anchored-deep-gram-cone-and-lossless-skew-target-refinement
  - THM-2531-prime-necklace-guard-boundary-selector
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2529-deep-comb-adjacent-fibre-odd-consumer-and-zero-target-boundary
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
script: 04-computation/lrc14_selector_target_drift_toothpick_boundary_thm2536.py
output: 05-knowledge/results/lrc14_selector_target_drift_toothpick_boundary_thm2536.out
script_sha256: 6864b37984da1900fb6fcbb952a63bae26e5ebe14d9b7d6f6d751500e5bd14d2
output_sha256: 79770220adf24c0b52b61a9711f35fb517479d78778a1f4e4fbc1859a6df6910
hash_basis: working-tree bytes (LF)
secondary_script: 04-computation/lrc14_selector_mass_drift_probe_codex_20260727.py
secondary_output: 05-knowledge/results/lrc14_selector_mass_drift_probe_codex_20260727.out
secondary_script_sha256: a5795aec23783af779bd855b756a26ca9bb630a48ba7768a7bf45df9176e97d7
secondary_output_sha256: 0177cf8634ca5e8321fc17c974d009f427eb7f15f80bdabed04361e579db13a4
---

# THM-2536 -- the two deep-comb selectors are a lossless boundary flow

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2531 finds a positive canonical occupied root in every target-anchored
deep-comb fibre.  On the actual one-run comb bank, that marker is only a
quadratic Boolean predicate.  This theorem determines exactly what the
resulting class masses remember, how their root colours differ from their
target characters, and why neither the marked-current algebra nor the
septimal toothpick can silently absorb them.

The central distinction is

```text
cellwise root Fourier support:       automatic in every positive rational cell;

target Fourier support:              equivalent to variation across target cells;

marked source-to-arrival current:    still requires a new semantic map.       (1)
```

## 1. Start and terminal selectors

Keep THM-2530's relative target chart.  In one lawful cell `(s,t)`, write

```text
e_a(x)=Delta_(t+a)(x),

K(a,b)=integral F e_a e_b,

K=sum_(j=1)^12 alpha_j e_j e_j^T
  +sum_(j=1)^11 beta_j(e_j+e_(j+1))(e_j+e_(j+1))^T.           (2)
```

As usual, put

```text
beta_0=beta_12=0.
```

There are two oriented boundary selectors:

```text
gamma^+_a
 =integral F e_a(1-e_(a-1))
 =K(a,a)-K(a-1,a),

gamma^-_a
 =integral F e_a(1-e_(a+1))
 =K(a,a)-K(a,a+1).                                      (3)
```

The first chooses the start of the unique occupied run in the positive
deep direction; the second chooses its terminal root.  Evaluating them on
the singleton and adjacent-pair rays in (2) gives

```text
gamma^+_j=alpha_j+beta_j,

gamma^-_j=alpha_j+beta_(j-1),                  1<=j<=12.      (4)
```

Thus both vectors are nonnegative and

```text
gamma^+_0=gamma^-_0=0,

sum_j gamma^+_j=sum_j gamma^-_j=mu(F).                       (5)
```

For THM-2531's arbitrary retained guard slope, its sign
`epsilon_tau in {+1,-1}` decides which of these two vectors occurs:

```text
gamma^tau=gamma^+          if epsilon_tau=+1,

gamma^tau=gamma^-          if epsilon_tau=-1.                 (6)
```

So the twelve guard slopes do not create twelve new mass vectors on this
bank.  They choose one of the two physical orientations of the same path.

## 2. The pair-flow divergence is lossless

Subtracting (4) gives the exact boundary-flow law

```text
gamma^+_j-gamma^-_j=beta_j-beta_(j-1).                        (7)
```

The anchored path, rather than a cycle, makes (7) invertible.  Starting
from `beta_0=0`, telescope:

```text
beta_j=sum_(i=1)^j(gamma^+_i-gamma^-_i),       1<=j<=11,

alpha_j=gamma^+_j-beta_j
       =gamma^-_j-beta_(j-1).                                  (8)
```

At `j=12`, equality of the two totals in (5) gives `beta_12=0`.
Consequently

```text
(gamma^+,gamma^-)
  <-> (alpha_1,...,alpha_12,beta_1,...,beta_11)
  <-> K                                                        (9)
```

is lossless.  This is a two-boundary alternative to THM-2530's skew-matrix
inverse.  Each selector separately forgets whether its chosen start or end
belongs to a singleton or a pair.  For example, a singleton `{j}` and an
adjacent pair `{j,j+1}` of the same mass have the same `gamma^+`, while
their `gamma^-` differ.  The opposite ambiguity holds after reversing the
path.

The signed opposite-selector difference is especially useful:

```text
delta_j:=gamma^+_j-gamma^-_j=beta_j-beta_(j-1),

sum_j delta_j=0,

delta=0 iff beta=0.                                           (9a)
```

The last equivalence is the anchored telescope from `beta_0=beta_12=0`.
Thus `delta` detects the existence of an adjacent double-deep ray in one
cell without retaining either selector's positive total mass.

Each vector also defines an honest occupied-to-excluded-target star:

```text
S^+_(a,0)= gamma^+_a,             S^+_(0,a)=-gamma^+_a,

S^-_(a,0)= gamma^-_a,             S^-_(0,a)=-gamma^-_a.        (10)
```

These are oriented stars with ties, not tournaments.  The endpoint `0` is
the excluded target root, not an occupied arrival.

There is an exact flow/homology invoice.  On the pair ray `{j,j+1}`,

```text
S^+-S^- = (j -> 0)+(0 -> j+1),

C_j       = (j -> j+1),                                      (10a)
```

where `C_j` is the literal occupied double-deep edge.  Both have vertex
divergence `e_j-e_(j+1)`, while

```text
(S^+-S^-)-C_j
 =(j -> 0)+(0 -> j+1)+(j+1 -> j)                             (10b)
```

is the directed triangle circulation through the excluded target.  Hence a
divergence-only or toothpick bridge identifies the literal pair edge with
its target detour only after forgetting this triangle holonomy.

## 3. Every positive rational cell has all root colours

Let `zeta=exp(2 pi i/13)` and define the unnormalised root transforms

```text
gamma_hat^epsilon(k)=sum_(a in F_13)gamma^epsilon_a zeta^(ka),

epsilon in {+,-}.                                             (11)
```

Assume the cell weights are rational and `mu(F)>0`.  If (11) vanished for
one `k!=0`, the rational coefficient polynomial of degree at most twelve
would vanish at a primitive thirteenth root.  Divisibility by

```text
Phi_13(X)=1+X+...+X^12
```

would force all thirteen coefficients to be equal.  The anchored zero in
(5) would then force every coefficient to vanish, contradicting positive
mass.  Hence

```text
gamma_hat^epsilon(k)!=0
             for every k!=0 and each positive selector cell.  (12)
```

This is a statement about the **relative root label inside one cell**.  It
does not yet mention a target character.

The pair divergence has a parallel saturation statement.  From (9a),

```text
delta_hat(k)=(1-zeta^k)sum_j beta_j zeta^(kj).                  (12a)
```

If the rational pair vector is nonzero, its anchored zeros prevent its
primitive transform from vanishing by the same `Phi_13` argument.  Therefore

```text
beta!=0
  -> delta_hat(k)!=0                    for every k!=0.         (12b)
```

This all-colour object is signed, but unlike a uniformly centered single
selector it is canonically the difference of two genuine opposite-slope
selector stars.

## 4. The selector target tensor and its exact drift dichotomy

Restore all target cells and define the absolute tensors

```text
Gamma^epsilon(r,s,t)=gamma^epsilon_(s,t),r-t.                  (13)
```

Use the normalized transform

```text
B^epsilon(k,b,h)
 =1/13^3 sum_(r,s,t)
    Gamma^epsilon(r,s,t)zeta^(kr+bs+ht).                       (14)
```

After setting `a=r-t`, equation (14) becomes

```text
B^epsilon(k,b,h)
 =1/13^3 sum_(a,s,t)
    gamma^epsilon_(s,t),a zeta^(ka+bs+(k+h)t).                 (15)
```

Therefore its target vector is exactly

```text
q=(b,k+h),                                                    (16)
```

the same target typing as THM-2365.  The selector colour `k` is not thereby
identified with the one-deep marked-current multiplier.

The anchor in (5) gives the complete line cancellation.  At fixed target
`(b,q_2)`, character orthogonality in `r-t` yields

```text
sum_k B^epsilon(k,b,q_2-k)=0.                                 (17)
```

The genuine star has a different and stronger zero-sum root object.  Its
vertex boundary is

```text
rho^epsilon_(s,t),a
 =gamma^epsilon_(s,t),a,                 a!=0,

rho^epsilon_(s,t),0=-mu_(s,t).                              (17a)
```

Thus `rho=div S`: every selected occupied root is a positive source and the
excluded target is the single negative sink.  It is rational, nonzero, and
has root total zero.  The `Phi_13` argument now gives

```text
rho_hat^epsilon_(s,t)(k)!=0              for every k!=0         (17b)
```

in every positive rational cell.

Let `B_rho` be the transform (14) with `gamma` replaced by `rho`, and put

```text
mu_hat(b,q)=1/13^2 sum_(s,t)mu_(s,t)zeta^(bs+qt).             (17c)
```

Root-total zero and the target sink give the two exact identities

```text
B_rho(0,b,q)=0,

sum_k B_rho(k,b,q-k)=-mu_hat(b,q).                            (17d)
```

The normalization in the second line is exact: summing over `k` contributes
`13` at relative root zero, cancelling one factor in (14).  Hence every
nonzero target character of the positive cell-mass surface directly forces
a `k!=0` star-boundary coefficient.  Unlike the zero line (17), no
Cauchy--Schwarz redistribution is needed.

Let

```text
gamma_bar^epsilon_a
 =1/13^2 sum_(s,t)gamma^epsilon_(s,t),a,

D_epsilon
 =1/13^3 sum_(a,s,t)
    |gamma^epsilon_(s,t),a-gamma_bar^epsilon_a|^2.             (18)
```

Finite Parseval and (15) give

```text
D_epsilon
 =sum_((b,k+h)!=0)|B^epsilon(k,b,h)|^2.                        (19)
```

On every nonzero target line, (17) and Cauchy--Schwarz give

```text
sum_(k!=0)|B^epsilon(k,b,q_2-k)|^2
 >=1/13 sum_k|B^epsilon(k,b,q_2-k)|^2.                         (20)
```

Consequently

```text
D_epsilon>0
  -> some coefficient has both a nonzero target and k!=0;

D_epsilon=0
  <-> gamma^epsilon_(s,t) is independent of (s,t).             (21)
```

In the zero branch, a positive rational constant cell vector still obeys
(12), but its complete spectrum is confined to

```text
b=0,                         h=-k.                             (22)
```

Thus all twelve root colours can coexist with no nonzero target character.

Combining (8) cell by cell gives the sharper two-selector alternative:

```text
D_+>0 or D_->0
  -> a positive selector target drift;

D_+=D_-=0
  <-> every alpha_j and beta_j is target-cell invariant
  <-> the complete 23-ray anchored path law is target-cell invariant.       (23)
```

One selector drift alone is not equivalent to THM-2530 Gram drift; the
singleton/pair ambiguity after (9) is the exact kernel.

For the adjacent-pair branch, put

```text
D_delta
 =1/13^3 sum_(a,s,t)|delta_(s,t),a-delta_bar_a|^2.             (23a)
```

The anchored inverse in (8) gives

```text
D_delta>0
  iff some beta_(s,t),j varies across target cells.            (23b)
```

Thus `delta` is precisely the divergence shadow of THM-2530's
adjacent-double-deep target-drift branch.  It may still have all root modes
and zero target drift when a nonzero `beta` profile is cellwise constant.

Likewise `gamma -> rho` in (17a) is injective, so the star-boundary family
varies exactly when its selector family varies.  The drift norms need not be
equal, but their zero/nonzero decisions are equivalent.

## 5. The target-invariant hostile is exact

THM-2529 supplies the comb-compatible rational control

```text
F_t(x)=Delta_(t+1)(x)d(13cx),

mu(F_t)=1/91.                                                 (24)
```

Its relative mask is always the singleton `{1}`.  Hence

```text
alpha_1=1/91,                  beta=0,

gamma^+=gamma^-=(1/91)delta_1                              (25)
```

Its pair divergence is zero and both star boundaries are the same constant
relative vector in every `(s,t)` cell:

```text
rho=(1/91)(delta_1-delta_0).
```

Equations (12) and (22) now give

```text
all twelve primitive root modes:       present;

both selector target drifts:           zero;

nonzero target characters:             absent.                 (26)
```

As in THM-2529, (24) is an exact nonnegative rational scalar-comb hostile
at the consequence level; it is not claimed to equal `E_(s,t)W` for one
canonical nine-coordinate owner packet.  It proves that target variation
does not follow from the anchored comb geometry alone.

### 5.1 Selector drift and THM-2365 diagonal drift are incomparable

Write the relative THM-2365 diagonal as

```text
h_j:=H(t+j,s,t)=K(j,j)
    =gamma^+_j+beta_(j-1)
    =gamma^-_j+beta_j.                                      (26a)
```

Thus the minimal sidecar missing from one start selector is the preceding
run-length-two mass.  The two selectors recover it by (8), but one selector,
its total mass, and the old diagonal have exact independent kernels.

Three two-cell controls lie entirely in the nonnegative anchored path cone:

```text
singleton {2}                 versus pair {2,3}:
  same mu and gamma^+, different h;

pair {2,3}                    versus singletons {2},{3}:
  same h, different mu and gamma^+;

pair {2,3}+singletons {6},{7}
  versus singletons {2},{3}+pair {6,7}:
  same h and same mu, different gamma^+.                       (26b)
```

Therefore positive selector drift neither follows from nor implies positive
THM-2365 `H`-drift, even after fixing cell mass.  Equation (23), which uses
both opposite selectors, is the lossless statement.

There is also a stronger lawful hostile than (24).  THM-2367 constructs a
full nine-factor bare-owner packet, on every valuation profile `(1,b,c)`,
whose peeled lower factors leave a target-cell-independent deep phase.  Put

```text
A=(13/30)(6/7)^3(1/7)=468/12005.                             (26c)
```

On the exact `182`-cell deep mesh, after excluding target root zero, its path
coordinates are

```text
alpha_j=A*(2/182),                   1<=j<=12,

beta_j =A*(12/182),                  1<=j<=11.                 (26d)
```

Hence, in every lawful target cell,

```text
mu=A*(156/182)=2808/84035,

gamma^+=A*(0,14,...,14,2)/182
       =(0,A/13,...,A/13,A/91),

h=A*(0,14,26,...,26,14)/182
  =(0,A/13,A/7,...,A/7,A/13).                                (26e)
```

The terminal selector is the reflected profile.  All of `gamma^+`,
`gamma^-`, `delta`, `rho`, and `h` are independent of `(s,t)`, although the
nonzero rational root profiles retain their primitive modes.  This is a
lawful full-packet, all-profile no-go for deriving selector or diagonal
target drift from packet form and valuation profile alone.  It is still a
typed non-cover packet with positive uncovered mass, not a scalar-cover row
or an LRC counterexample.

## 6. The selector is not already the THM-2334/2365 current

The start selector has the exact physical decomposition

```text
gamma^+_a
 =integral F Delta_(t+a)
  -integral F Delta_(t+a-1)Delta_(t+a),                       (27)
```

and `gamma^-` has the analogous forward adjacent overlap.  The first term
in (27) is THM-2365's one-deep tensor entry.  The second is THM-2530's
adjacent-double-deep higher moment.  Their difference is nonnegative
because the integrand is the literal Boolean selector

```text
F Delta_(t+a)(1-Delta_(t+a-1)).                               (28)
```

Therefore `gamma` is a lawful **same-time scalar mass observable**.  It is
not already a coefficient of THM-2334's marked relation current or
THM-2365's one-deep scalar row: those constructions contain one deep leg,
whereas (28) uses an additional adjacent deep complement.  Linear recovery
from the larger Gram matrix does not retroactively place that higher moment
in the smaller scalar algebra.

The typing boundary is equally important.  Both factors in (28) live at
the same predecessor horizon.  THM-2533 can multiply this root fibre by a
root-invariant terminal-word/future-owner factor without changing the
selector identities, so those data can coexist on one Boolean event.  The
excluded target root still has zero selector mass: neither (27) nor the star
(10) turns that empty predecessor into the later owner, a semantic arrival,
or a source-to-arrival edge.

There is a cheaper star which must not be confused with the selector star.
Since root zero is excluded, every occupied root has an intrinsic chord to
zero.  Integrating **all** of them gives

```text
H_a=K(a,a)=integral F e_a,

S_H=sum_(a!=0)H_a(a -> 0),

rho_H=div S_H=H-(sum_a H_a)delta_0.                           (28a)
```

This all-occupied-root star is already in THM-2365's one-deep algebra.  Its
boundary is rational, zero-sum, and has all primitive root modes in every
positive cell; its target-cell drift vanishes exactly when the old `H`-drift
vanishes (the two numerical norms need not be equal).  Thus
THM-2536 is not the first existence proof for an ordered target star.  The
selector cost buys something different: exactly one canonical occupied
root per fibre, start/terminal run data, and the lossless path-law recovery
(8)--(9).  On a two-root mask, `S_H` emits both target chords and forgets
which root was the start.  It also leaves the old zero-drift and semantic-
arrival boundaries completely unchanged.

## 7. Exact canonical 169-cell positive control

There is nevertheless a strong positive control on the exact same typed
row and delayed word used by the independently verified THM-2334
`169`-twist computation.  Take

```text
w=(1,14,27,40,53,66,13,13^3,2*13^5),

owner c_1=13,              target blockers c_2,c_3,

R=13^2=169,

Q=Q_(1,{c_2}).                                                 (29)
```

For every one of the `169` lawful quotient-dual coordinate shifts `ell`,
put

```text
F_ell=E_1^ell intersect T^(-2)Q,

mu_ell=mu(F_ell)=sum_a gamma^epsilon_(ell),a.                  (30)
```

The word is shift-neutral because `13|R`.  Exact rational interval
refinement gives

```text
#{ell:mu_ell>0}=169,

#{distinct mu_ell}=75,                                        (31)

min mu_ell=22841731/33256714010,

max mu_ell=17544684733/8962684425695,

mu_(0,0)=21376087/17907461390,

mu_(0,1)=595072523707/483984958987530.                         (32)
```

In particular, the selector family has positive target drift on this
control.  The exact statement is stronger.  Fourier transform the integer
overlap numerators in any quotient basis.  Under the ring map

```text
Z[zeta_13] -> F_79,                   zeta_13 -> 8,             (33)
```

all `168` nontrivial character sums have nonzero image.  Since `8` has
exact order thirteen, each nonzero image certifies that the corresponding
complex cyclotomic coefficient is nonzero.  Thus

```text
mu_hat(q)!=0                     for every q!=0                (34)
```

on this canonical typed control.  Equations (17d) and (34) say more
directly that every one of these `168` target lines has nonzero star-boundary
sum, while its `k=0` entry vanishes.  Hence every nontrivial target is
accompanied by at least one nonzero star-boundary root colour `k!=0`.

This does **not** prove a uniform theorem.  The row in (29) is the canon's
typed non-cover example, not a hypothetical covering row.  Equations
(24)--(26) remain the exact structural hostile.

## 8. Relation to the independently verified marked-current variation

The control in Section 7 and THM-2334's exact `169`-twist control use the
same:

```text
row and delayed word,
quotient-dual translation group,
169 lawful coordinate shifts,
target Fourier labels.                                        (35)
```

They do not use the same functional.  Section 7 evaluates the positive
zero-frequency masses

```text
mu_ell=hat(E_1^ell 1_(T^2x in Q))(0).                          (36)
```

THM-2334 evaluates a complex marked triangle of the form

```text
deep_phase(ell)
 *hat(E_1^ell 1_(T^2x in Q))(X)
 *conjugate(hat(E_1^ell)(Y)),              X,Y!=0.              (37)
```

Its independent referee proves nonconstancy of the `169` values in (37),
and hence at least one nonzero marked target aggregate.  Equation (34)
proves every nontrivial target character of the different mass surface
(36).  Neither statement implies the other without an intertwiner between
the zero-frequency selector mass and the marked `(X,Y)` triangle.

This is the precise connection: the target-action representation is shared;
the observable, phase, and consequence are not.

## 9. Star boundary and the exact two-tooth lift

Every THM-2508 toothpick output obeys

```text
sum_v R_(tau,a,c)(v)=0.                                      (38)
```

The raw selector instead has positive total `mu(F)`.  Therefore no
toothpick output can equal `gamma^+` or `gamma^-` in a positive cell.  A
uniform mean-centering is one possible zero-sum target:

```text
gamma^epsilon,circ_a
 =gamma^epsilon_a-mu(F)/13,

sum_a gamma^epsilon,circ_a=0.                                (39)
```

Centering changes only root colour zero, so every primitive root coefficient
in (12), and every target coefficient with `k!=0` in (20), survives.
It does, however, smear the target sink across all roots.  The more natural
object is the anchored star boundary `rho` from (17a), which retains one
positive source mass at every selected occupied root and the full negative
mass at the excluded target.

There is an atomwise exact two-tooth realization.  Fix a selector class
`a!=0` of mass `m=gamma_a`, and on the baseline septimal columns put

```text
d^a(h,r)=m 1_(h=a)(1_(r=0)-1_(r=1)).                          (40)
```

This is row-zero.  Choose the root edge slope

```text
tau_a=-a.                                                      (41)
```

Then the baseline THM-2508 toothpick is

```text
R_(tau_a)d^a(v)
 =d^a(v,0)+d^a(v+a,1)
 =m(1_(v=a)-1_(v=0)).                                        (42)
```

At the target output `v=0`, tooth `r=1` samples the source root `a` and
contributes the signed nonzero value `d^a(a,1)=-m`; at source output `v=a`,
tooth `r=0` contributes `+m`.  Thus every genuine occupied-to-excluded-
target star edge has a two-tooth boundary realization.  Summing over all
classes realizes `rho`, but across the twelve slope components `tau_a`, not
as one fixed toothpick.  The root-dependent edge slope is affine-covariant;
it is generally not the retained physical guard slope.

THM-2535 gives the compatible cut-clock description.  The edge slope makes
the target tooth the distinguished tooth `j=1`, and its C7 label transports
covariantly under the full cut chart.  More strongly, THM-2535 correlates
the selector with actual scheduler-clock cells and constructs a positive
full-clock-orbit incidence dipole.  It does not identify the signed table
`d^a` in (40) with the inherited THM-2449 source/deep residue, make one
incidence component a temporal transition, or turn root zero into an
arrival.  Thus the scheduler clock can be physical while the required old
source/deep and arrival semantics remain absent.

Concretely, for a positive marker class `E_a`, THM-2535 uses the seven
scheduler intersections

```text
X_kappa=E_a intersect P_kappa,             x_kappa=mu(X_kappa),

D_kappa=delta_(h=a)
 (x_kappa delta_(r=kappa)-x_(kappa+1)delta_(r=kappa+1)).       (42a)
```

In the transported chart `(tau,A,b)=(-a,1,-kappa)`,

```text
R_kappa=x_kappa delta_a-x_(kappa+1)delta_0,

sum_kappa R_kappa=M(delta_a-delta_0),       M=sum_kappa x_kappa>0. (42b)
```

This is a physically scheduler-clock-labelled version of the star boundary.
Its uncharted signed table satisfies

```text
sum_kappa D_kappa=0,                                           (42c)
```

so the nonzero dipole is chart holonomy around the full clock orbit, not one
common THM-2508/source-row defect or one observed temporal transition.

If one insists on one fixed `tau!=0`, solve the cyclic antiderivative

```text
g(v)-g(v-tau)=rho^epsilon_v.                                  (43)
```

The solution exists because the right side sums to zero.  On two of the
seven columns put

```text
d(h,0)=g(h),

d(h,1)=-g(h),

d(h,r)=0,                         r=2,...,6.                   (44)
```

Then `sum_r d(h,r)=0` and the baseline toothpick gives

```text
R_tau d(v)=d(v,0)+d(v-tau,1)=rho^epsilon_v.                    (45)
```

This one-slope form is an exact surjectivity statement, but it forgets the atomwise star-edge
decomposition.  The opposite-selector divergence `delta` also sums to zero,
so the same construction realizes it directly, with no mean subtraction.

All three lifts -- uniform mean-centering, anchored star boundary, and pair
divergence -- remain signed table constructions.  Their exact invoice is:

1. the ordered septimal columns `0,1` are cut-chart labels, not deep-comb
   ancestry atoms;
2. the atomwise version uses the edge slope instead of the generally
   different physical guard slope, while the one-slope version spends an
   antiderivative gauge;
3. a toothpick output sees only vertex divergence, so (10b)'s target-
   triangle circulation and all other cycle holonomy are forgotten; and
4. no signed target tooth turns the excluded root into a semantic arrival,
   source handoff, or future owner.

The target anchor and literal positive star remain recoverable from `rho`
only while root `0` and the star-edge sidecar are retained.  The bare
toothpick divergence alone does not retain them.  Applying (43)--(45)
separately in every target cell preserves the nonzero-`k` target Fourier
data; it still does not turn that data into THM-2334's marked current.

## 10. Scope and exact referee

The proved refinement is

```text
one oriented selector drift
  -> nonzero selector target/root coefficient;

both selector drifts zero
  <-> full 23-ray anchored path-law invariance;

positive selector star -> zero-sum anchored boundary
  -> exact atomwise two-tooth cut bank;

opposite selectors -> adjacent-pair divergence
  -> direct zero-sum toothpick lift;

both lifts forget target-triangle holonomy and do not create ancestry.       (46)
```

No theorem here forces either selector drift on a hypothetical covering row,
places the adjacent-double-deep factor in the existing one-deep scalar
algebra, identifies the excluded target with an arrival, excludes a scalar
speed row, or proves LRC(14).

Run

```bash
python3 04-computation/lrc14_selector_target_drift_toothpick_boundary_thm2536.py
python3 -O 04-computation/lrc14_selector_target_drift_toothpick_boundary_thm2536.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_selector_target_drift_toothpick_boundary_thm2536.out
```

byte-for-byte.  The exact companion checks all `23` path masks, the
start/terminal inverse, pair-divergence saturation, all eleven target-
triangle holonomy controls, the selector/star target-line identities, three
selector/diagonal incomparability controls, the singleton and THM-2367 full-
packet hostiles, all twelve atomwise star-edge, one-slope star-boundary,
uniform-centering, and direct-divergence lifts, the `169` canonical rational
cell masses, and the `168/168` finite-field cyclotomic target certificates.

The secondary exact probe independently checks the three path-cone
incomparability controls, the `182`-mesh deep-phase counts, the THM-2367
peeling constant, and its complete selector/diagonal profiles.

**QED.**

The independent audit rederived the anchored telescope recovering every
`alpha_j,beta_j`, the cyclotomic all-root-mode argument, the exact target-line
normalization in (17d), and the target-triangle circulation lost by vertex
divergence.  It separately checked the three minimal selector/diagonal
controls, the THM-2367 full-packet profile, and both toothpick formulas.
Normal and optimized executions of both exact scripts byte-match their
stored transcripts and all recorded hashes.
