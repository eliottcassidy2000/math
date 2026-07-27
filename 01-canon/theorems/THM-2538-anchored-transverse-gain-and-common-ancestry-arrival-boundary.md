---
id: THM-2538
title: "Anchored transverse gain and common-ancestry arrival boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY SPLIT-AUDITED.  Extending the
  THM-2533 mixed table through its
  lawful deep-target coordinate lifts each of the three selected gains to
  one full deep-target character at the same ordinary frequency.  Moving
  the deep target by c*tau together with root translation by tau preserves
  the THM-2530 excluded-root anchor exactly.  For a fixed nonzero root step,
  at least nine nonzero target shifts are transverse to all three lifted
  characters; atomwise exclusive-boundary splitting and four-scalar
  avoidance then give one rational nonnegative anchored table retaining all
  three frequencies and 216 Galois-saturated gain incidences.  An explicit
  finite anchored singleton table has all 216 incidences and atomwise
  disjoint transverse boundaries while its complete anchored Gram and both
  selector vectors are constant in every target cell.  Thus preserving the
  anchor still does not couple phase to deep-mask location.  More sharply,
  the empty head of a same-horizon boundary annihilates the event at that
  endpoint and cannot itself be an arrival.  If a genuine later Boolean
  target field is supplied on the same ancestry base, products formed
  before integration recover the complete source/arrival coupling matrix
  positively and exactly; predecessor and later fields may use different
  retained empty anchors.  Applied to THM-2537's selected head, this gives
  an exact cross-Kakeya form of its still-open positive target hit (56).
  Integrating the two endpoint banks separately loses a
  (m-1)(n-1)-dimensional mixed-Haar transportation kernel; sharp
  Frechet--Hall bounds and two- and three-state controls preserve every
  endpoint bank while reversing the arrival Haar sign or changing the
  candidate transport-loop type.  The live LRC carrier does not yet prove
  positive same-ancestry support and chronological typing for the formal
  later target-active field, so no semantic arrival, owner loop, row
  exclusion, or LRC(14) proof follows.
source: codex-2026-07-27-anchored-transverse-arrival
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2530-anchored-deep-gram-cone-and-lossless-skew-target-refinement
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2534-positive-kakeya-boundary-transform-and-crofton-reconstruction
  - THM-2536-deep-comb-selector-flow-target-drift-and-centered-toothpick-boundary
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2540-weighted-live-event-kakeya-flux-and-transverse-gain-boundary-refinement
related:
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
  - HYP-2991-lrc14-haar-zipper-cocycle
script: 04-computation/lrc14_anchored_transverse_arrival_boundary_thm2538.py
output: 05-knowledge/results/lrc14_anchored_transverse_arrival_boundary_thm2538.out
script_sha256: 2f9d0afad9d1a0d1b81501bf79d444f21e997c54cc16962cb3f439e3b1f75b55
output_sha256: 12dc0aea42982f7199db81cd4dbd6c2763e8569da23ed1ee6ad595648b5820fd
hash_basis: working-tree bytes (LF)
---

# THM-2538 -- the anchor can move covariantly; arrival must share ancestry

**PROVED + VERIFIED-EXACT.**

THM-2540's product-torus derivative creates atomwise exclusive transverse
target--root boundaries, but its uncorrected product shift moves the
THM-2530 deep target
anchor.  The first result below repairs that coordinate: the deep target must
move with its physical root colour.  This preserves the anchor and all three
THM-2533 gains.

The repair exposes a sharper obstruction.  Gain phase, deep-mask location,
and temporal ancestry are three different coordinates.  A finite hostile
keeps the first two while making the deep location completely constant.  A
second exact calculation shows what is still absent: a semantic transition
is a **cross-horizon coupling**, obtained by multiplying predecessor and
arrival events before integration.  Separate endpoint Gram or Kakeya banks
retain only its margins on the singleton subfamily.

The inheritance ledger is therefore

```text
closest positive mechanism:   THM-2540 transverse exclusive boundaries;
closest anchored mechanism:   THM-2530/2536 target stars and zero row;
canonical hostile:            THM-2529 target-constant singleton law;
repaired historical warning:  HYP-2991's fixed-margin mixed Haar switch;
new exact sidecar:             one later Boolean field on the same ancestry base.
                                                                    (1)
```

## 1. Lift the three gains through the deep-target character

Let `xi=zeta_7`, `zeta=zeta_13`, and let

```text
B_(ell,s,t,alpha):T->{0,1},

ell in F_7,          s,t in F_13                              (2)
```

be a finite rational-BV Boolean atomization of the full lawful mixed table.
The target coordinate `t` is the deepest target root.  At `t=0`, define

```text
W^0_(kappa,b)
 =1/91 sum_(ell,s,alpha)
    xi^(kappa ell)zeta^(bs)B_(ell,s,0,alpha),                 (3)
```

and on the full family define

```text
W_(kappa,b,q)
 =1/(91*13) sum_(ell,s,t,alpha)
    xi^(kappa ell)zeta^(bs+qt)B_(ell,s,t,alpha).              (4)
```

Finite Fourier inversion in `t` gives the exact evaluation identity

```text
W^0_(kappa,b)=sum_(q in F_13)W_(kappa,b,q).                  (5)
```

THM-2533 supplies one mixed channel `(kappa_0,b_0)`, three distinct gains

```text
lambda_i=a_i/b_0,                         i=1,2,3,             (6)
```

and positive integers `n_i=a_i mod 13` for which

```text
(W^0_(kappa_0,b_0))hat(n_i)!=0.                              (7)
```

Equation (5), at the same ordinary frequency, selects some

```text
q_i in F_13
```

such that

```text
(W_(kappa_0,b_0,q_i))hat(n_i)!=0.                            (8)
```

The three `q_i` need not agree.  No positivity or rationality of a
cyclotomic absolute-square energy is used here; (8) is only finite linear
expansion of the already-proved nonzero coefficients.

Here is the promised full-family typing.  Let `h_(ell,r,s,t)` be the lawful
THM-2449 Boolean packet integrand used in THM-2533, and let `g_own` be that
theorem's common sufficiently delayed Boolean owner--word factor.  In its
first branch take

```text
B_(ell,s,t,r)(x)
 =g_own(x)h_(ell,r,s,t)(x).                                  (8a)
```

In its delayed-square branch take the atom label `alpha=(r,r')` and

```text
B_(ell,s,t,r,r')(x)
 =g_own(x)h_(ell,r,s,t)(x)
             h_(ell,r',s,t)(13^Q x).                        (8b)
```

At `t=0`, summing `r`, or `(r,r')`, gives exactly THM-2533's first or
square channel.  Every atom is finite rational-BV and Boolean.  The lawful
present packet gives (10), while the delayed owner is target-neutral; in
the square branch both packet factors retain their displayed target label
rather than silently treating a target-dependent far factor as common.
Thus (2)--(8) are a construction from the canonical lawful family, not an
extra existence premise.

## 2. The deep-anchor-covariant product shift

Write the deepest danger indicators as

```text
Delta_t(x)=1_(||c x-t/13||<1/14),

c not=0 mod 13.                                               (9)
```

The lawful target exclusion is

```text
B_(ell,s,t,alpha)Delta_t=0.                                  (10)
```

For `tau!=0` and a target co-shift `h`, define

```text
(P^anc_(tau,h)B)_(ell,s,t,alpha)(x)
 =B_(ell,s+h,t+c tau,alpha)(x+tau/13).                       (11)
```

The target correction in (11) is forced by the exact identity

```text
Delta_(t+c tau)(x+tau/13)=Delta_t(x).                        (12)
```

Consequently `P^anc B` obeys the same excluded-root law (10).  This is the
coordinate missing from the uncorrected product shift; it is a transported
target anchor, not an invariant absolute target label.

At an ordinary frequency `n=a mod 13`, reindexing `s,t` and translating
`x` gives

```text
(P^anc_(tau,h)W_(kappa,b,q))hat(n)
 =zeta^(a tau-bh-qc tau)What_(kappa,b,q)(n).                 (13)
```

Thus the lifted character `(a,b,q)` has anchored characteristic direction

```text
h=b^(-1)(a-cq)tau.                                          (14)
```

For the three rows in (8), each equation (14) forbids at most one `h`.
Among the twelve nonzero target shifts, at least

```text
12-3=9                                                       (15)
```

are transverse to all three.  The three forbidden values may coincide or
equal zero, in which case more than nine choices survive.  Distinctness of
the original gains in (6) is not used for (15); it will make their Galois
orbits disjoint.

## 3. One anchored nonnegative pack retains all three gains

Fix one common transverse `(tau,h)` from (15).  Coordinatewise on the full
atom-indexed table put

```text
A=P^anc_(tau,h)B,

U=A(1-B),                V=B(1-A).                           (16)
```

Pointwise,

```text
A-B=U-V.                                                     (17)
```

Both `U,V` are literal Boolean exclusive-boundary orientations within each
atom.  Since
`A` and `B` both vanish on `Delta_t`, so do `U,V`.  Let `u_i,v_i` be their
three selected `(kappa_0,b_0,q_i,n_i)` complex linear functionals **after
summing the atom table**.  Equations (8), (13), and transversality give

```text
u_i-v_i!=0,                              i=1,2,3.             (18)
```

For fixed `i`, the equation

```text
u_i+r v_i=0                                                   (19)
```

excludes at most one scalar `r`.  If `v_i=0`, it excludes none because
`u_i!=0`; otherwise it has one complex solution.  Three rows cannot cover
four prescribed positive integers.  Hence some

```text
r in {1,2,3,4}                                                (20)
```

makes

```text
Z=U+rV>=0                                                     (21)
```

retain all three selected ordinary frequencies.  The table is rational and
still obeys the deep anchor.  It need not be Boolean, and exclusivity is
atomwise only; different atom labels may overlap after summation.

The seed ordinary coefficient implies the root projector is a nonzero
function:

```text
E_(a_i) Z_(kappa_0,b_0,q_i)!=0.                              (21a)
```

Rationality of the underlying atom table gives exact pointwise
source/target/deep-target/root Galois covariance

```text
(kappa,b,q,a)
 ->(v_7 kappa,v_13 b,v_13 q,v_13 a).                         (22)
```

Indeed the automorphism is applied to the finite character/projector sum,
not to arbitrary rational-endpoint ordinary Fourier phases.  It sends
`E_a Z_(kappa,b,q)` to the corresponding projector in (22), so nonvanishing
is preserved.  Because `kappa_0,b_0` are nonzero, each orbit has size
`6*12=72` even when `q_i=0`.  The three ratios `a_i/b_0` are distinct, so
their orbits are disjoint and have

```text
3*6*12=216                                                   (23)
```

incidences.  This is a projector statement.  The exact ordinary frequency
need not be common around an orbit; a grouped-jump lift may be applied in
each residue separately.  No order on a cyclotomic field and no rationality
of an absolute-square energy is asserted.

Thus the missing **deep anchor** in the transverse construction is
repairable without losing its gains.

## 4. Anchor plus all gains still does not force deep-mask drift

The repair is sharp.  There is a finite rational table satisfying all of the
algebra above while its anchored Gram is constant in every target cell.
For this hostile specialize the lawful deep multiplier in (9)--(14) to
`c=1`; no general-`c` claim is needed for sharpness.

Use a uniform root fibre `u in F_13`, three gains

```text
lambda_1=1,          lambda_2=2,          lambda_3=3,          (24)
```

and the ordered carry section

```text
[d]_7=(the representative of d in {0,...,12}) mod 7.         (25)
```

For atom type `i`, define

```text
B_(ell,s,t,i)(u)
 =1_(u=t+1) 1_(ell=[s-lambda_i t]_7).                        (26)
```

Every atom excludes root `t`.  On its support the complete relative deep
mask is the singleton `{1}`.  Summing the three atom types gives the same
mass and the same scalar multiple of `e_1e_1^T` in every `(s,t)` cell.

For `kappa,b!=0`, the source/target kernel is

```text
S_(kappa,b)
 =sum_(d=0)^12 xi^(kappa[d]_7)zeta^(bd)
 =sum_(d=0)^12 (xi^kappa zeta^b)^d.                          (27)
```

If `rho=xi^kappa zeta^b`, then

```text
rho^13=xi^(13kappa)=xi^(6kappa)!=1.                          (28)
```

Therefore the geometric sum in (27) is nonzero.  At deep-target character
`q=0`, atom `i` has the nonzero ordinary root mode

```text
a=lambda_i b.                                                (29)
```

This already gives all `216` primitive gain incidences.

Now take

```text
tau=1,                       h=4.                             (30)
```

For the three atom types, the section displacements `h-lambda_i tau` are
`3,2,1`.  The section in (25) changes under each of them at every cyclic
input: before wrap the change is `1,2,3 mod 7`, and after wrap it is
`2,3,4 mod 7`.  Thus `A=P^anc B` and `B` are disjoint within every atom.
Take the fixed positive pack

```text
Z=A+2B.                                                       (31)
```

At the gain in (29), its coefficient is the nonzero coefficient of `B`
times

```text
zeta^(b(lambda_i-4))+2!=0.                                  (32)
```

The inequality follows from absolute values over `C`.  Yet in every target
cell the packed table still has only the relative singleton `{1}`, with
constant Gram

```text
K_(s,t)=9/13 e_1e_1^T.                                      (33)
```

Accordingly its THM-2530 Gram drift, both THM-2536 selector drifts, and the
pair-flow divergence all vanish.  The table has a moving source/target
phase and a fixed deep location.  It is an exact finite hostile to

```text
deep anchor + three gains + positive exclusive-boundary pack
  => anchored Gram or selector target drift.                 (34)
```

This control satisfies the stated finite table, covariance, and anchor
axioms.  It is not asserted to be one physically realized THM-2449 owner
packet or a covering row.  Its mechanism identifies the missing coordinate:
the phase must be coupled to **deep-mask location**, not merely coexist with
an excluded-root anchor.

## 5. A same-horizon empty head cannot be its own arrival

The remaining semantic gap has an even sharper pointwise boundary.  For a
Boolean root field `e=(e_r)` define its directed Kakeya coordinate

```text
K^e_tau(r)=e_r(1-e_(r+tau)).                                 (35)
```

For distinct roots `r,s`,

```text
K^e_(s-r)(r)e_s
 =e_r(1-e_s)e_s
 =0.                                                         (36)
```

Thus the empty destination of the predecessor boundary literally
annihilates the same predicate at that horizon.  No relabelling, Fourier
mode, positive scalarization, or target anchor can turn that empty head into
an arrival of `e` itself.

This does not say arrival is impossible.  It says its Boolean event must be
a genuinely different field, normally at a later horizon or on a different
semantic factor.

## 6. One genuine later field gives the exact positive cross-Kakeya target hit

Let `(Omega,mu)` be one common ancestry base, let `F>=0` be a common weight,
and let

```text
e=(e_r)_(r in F_13),       a=(a_s)_(s in F_13)              (37)
```

be Boolean predecessor and later-arrival root fields on that same base.
Assume they retain (not necessarily the same) empty anchors

```text
e_(r_*)=0,                       a_(s_*)=0.                   (38)
```

The one boundary directed from an occupied coordinate to its own excluded
root gives

```text
K^e_(r_*-r)(r)=e_r,              r!=r_*,

K^a_(s_*-s)(s)=a_s,              s!=s_*.                    (39)
```

Hence the complete ordered source/arrival coupling is recovered entrywise
by the positive exact formula

```text
Gamma_(r,s)
 :=integral_Omega F e_r a_s dmu

 =integral_Omega F
    K^e_(r_*-r)(r)K^a_(s_*-s)(s)dmu,
                                      r!=r_*, s!=s_*.        (40)
```

Every bilinear arrival statistic, diagonal loop mass, mixed Haar mode, or
transport-cycle weight which is a linear functional of `Gamma` is therefore available
once the later field exists and the product in (40) is formed **before**
integration.  Formula (40) is a cross-Gram/cross-Kakeya bank; it is not
another self-correlation.

Different anchors matter semantically.  If the desired later arrival is at
the predecessor anchor `r_*`, imposing `a_(r_*)=0` would delete the event by
definition.  One must retain a different later empty anchor `s_*`, or use
the all-slope reconstruction below.  The two horizons need a common ancestry
base, not a common empty coordinate.

This types exactly against THM-2537's remaining positive hit.  On its mixed
predecessor locus define the categorical selected-head field and its positive
threshold layers by

```text
H_t(z)=1_(e(z) nonconstant)1_(t_tau(e(z))=t),

S_(j,t)(z)=g(z)1_(Psi_tau(e(z))>=j)H_t(z),
                                      1<=j<=98.              (40a)
```

Put `A_t(z)=A_tar(z,t)`.  Then THM-2537 equation (56) is identically

```text
integral g(z)Psi_tau(e(z))A_tar(z,t_tau(e(z)))dz
 =sum_(j=1)^98 sum_(t in F_13)
    integral S_(j,t)(z)A_t(z)dz.                             (40b)
```

Because `H` has either zero or one occupied coordinate, for every fixed
nonzero slope `nu`

```text
K^H_nu(t)=H_t.                                               (40c)
```

Thus every summand of (40b) is already a positive same-ancestry cross
product.  If the later field has a retained empty anchor, (39) turns it into
a product of two Kakeya coordinates.  If it has no globally fixed empty
anchor, retain an occupancy/anchor branch or use (41)--(42).  In particular,
THM-2537 has supplied the genuine selected wall and positive source layers;
the only absent theorem in (40b) is proved positive same-ancestry support
and chronological typing for its formally defined target-active field `A`.

There is a sharper all-slope form which uses the original predecessor
boundaries themselves.  For any root-indexed field `A` on the same base and
`m_e=sum_r e_r`, reindexing `s=r+nu` gives the pointwise needle identity

```text
sum_(nu!=0) sum_r K^e_nu(r)A_(r+nu)
 =sum_s (1-e_s)A_s sum_(r!=s)e_r
 =m_e sum_s (1-e_s)A_s.                                    (40d)
```

Define the selected target-active field

```text
A^sel_s=H_s A_s.                                            (40e)
```

The selected head is empty in `e`, so on every mixed fibre the right side
of (40d) for `A^sel` is exactly `m_e A_tar(z,t_tau(e(z)))`.
Consequently the open hit in (40b) has the positive all-needle expansion

```text
integral g Psi_tau(e) A_tar(z,t_tau(e(z))) dz
 =sum_(m=1)^12 1/m sum_(nu!=0) sum_r
 integral 1_(m_e=m) g Psi_tau(e)
   K^e_nu(r) H_(r+nu) A_(r+nu) dz.                          (40f)
```

Every term is nonnegative.  Formula (40f) removes slope choice as an
obstruction: the complete boundary bank automatically sends `m_e` needles
into any target-active selected empty head.  It does **not** create that
target-active event.  The remaining problem is exactly the existence and
alignment of `A_tar` on the selected head, not another boundary or Fourier
census.

The anchor is the cheapest inverse.  Without it, the complete all-slope
identity is

```text
sum_(tau!=0)K^e_tau(r)=(13-m_e)e_r,

m_e=sum_u e_u.                                               (41)
```

Thus on a branch with fixed occupancies `m_e=m` and `m_a=n`,

```text
Gamma_(r,s)
 =1/[(13-m)(13-n)]
   sum_(tau,sigma!=0)
    integral F K^e_tau(r)K^a_sigma(s).                        (42)
```

For the actual singleton/adjacent-pair comb, retaining the two occupancy
bits and summing the four branch-conditioned versions of (42) gives the
unweighted coupling.  Equation (40) avoids this branch ledger entirely.
If an occupancy is `13`, its field is identically one and its corresponding
factor is recovered from the separately retained branch weight or marginal
sidecar, not from its zero Kakeya bank.  With that sidecar the zero
denominator in (42) causes no loss; the boundary bank alone would lose the
constant-one branch.

The statement is conditional in exactly one semantic place: current canon
has the predecessor selector/star field `e` and the formal predicate
`A_tar`, but it has not proved positive same-ancestry support or
chronological arrival typing for that later Boolean field.

## 7. Separate endpoint banks lose the transportation kernel

The order of operations in (40) is load-bearing.  Consider the singleton
subfamily

```text
e_r=1_(X=r),                   a_s=1_(Y=s),                   (43)
```

with `X,Y` supported on the twelve anchored roots `1,...,12`.  At each
horizon, the Gram and all-slope Kakeya banks are exactly equivalent to the
one-point marginal:

```text
G^X_(r,r)=Pr(X=r),             G^X_(r,s)=0 for r!=s,

integral K^X_tau(r)=Pr(X=r),   tau!=0,                       (44)
```

and similarly for `Y`.  They say nothing about how the two endpoint atoms
are paired.

The joint matrix

```text
C_(r,s)=Pr(X=r,Y=s)                                             (45)
```

is a nonnegative coupling with the prescribed row and column margins.
Conversely every rational such coupling is realized by a finite rational
probability space.  If the positive source and arrival supports have sizes
`m,n`, the affine transportation fibre has dimension

```text
(m-1)(n-1).                                                   (46)
```

Its linear kernel is the space of matrices with zero row and column sums.
After choosing base indices `r_0,s_0`, a basis is

```text
E_(r,s)-E_(r,s_0)-E_(r_0,s)+E_(r_0,s_0),

r!=r_0,                         s!=s_0.                       (47)
```

These are precisely mixed two-dimensional Haar switches.  On all twelve
anchored roots the missing space has dimension

```text
11^2=121.                                                     (48)
```

The minimal hostile uses two equiprobable source states `1,2`:

```text
C_+=1/2 [[1,0],[0,1]],

C_-=1/2 [[0,1],[1,0]].                                      (49)
```

Both have the same endpoint Gram and every endpoint Kakeya mass.  But the
mixed Haar coordinate

```text
zeta_H=C_11-C_12-C_21+C_22                                  (50)
```

is `+1` on `C_+` and `-1` on `C_-`; their diagonal transport masses are `1` and
`0`.  This is the exact mathematical core anticipated, but not proved as a
live LRC theorem, by HYP-2991.

The one-coordinate obstruction is already sharp.  For two events `E,A` on
a base of total mass `M`, with masses `p,q`, the Frechet bounds are

```text
max(0,p+q-M)<=mu(E intersect A)<=min(p,q).                   (50a)
```

Both endpoints are attainable.  Consequently two separately positive
marginals do not force a positive hit unless their masses sum to more than
the entire base.  For categorical endpoint laws `(p_i)` and `(q_i)` of
total mass `M`, the diagonal-forbidden transportation graph and Hall's
theorem give the exact refinement

```text
minimum_C sum_i C_(i,i)
 =max_i (p_i+q_i-M)_+;                                      (50b)

there is a zero-diagonal coupling
 iff p_i+q_i<=M for every i.                                 (50c)
```

Only a common-base intersection estimate, or a structural restriction on
the transportation fibre, can improve these sharp marginal bounds.

There is one further arity boundary for an owner loop.  THM-2478 can place
a genuine future-owner handoff on an old source base.  Even after a later
arrival field has also been placed there, the three pairwise cross banks do
not in general determine the source--arrival--owner incidence.  The smallest
hostile is THM-2221's even/odd parity pair, read as laws of three Boolean
variables:

```text
mu_even =uniform{000,011,101,110},
mu_odd  =uniform{111,100,010,001}.                           (50d)
```

All one- and two-variable marginals agree, but

```text
E_even[XYZ]=0,                     E_odd[XYZ]=1/4,           (50e)
```

and their three-way parity characters have opposite signs.  Thus
source--arrival is a degree-two coupling problem, while a fully marked
owner-loop incidence can be degree three.  Once Boolean root fields
`e,a,o` genuinely share one ancestry base and retain respective empty
anchors, the exact repair is simply the three-fold positive bank

```text
Theta_(r,s,u)=integral F e_r a_s o_u
 =integral F K^e_(r_*-r)(r)
             K^a_(s_*-s)(s)K^o_(u_*-u)(u),                 (50f)
```

away from the three zero anchor coordinates; all-slope occupancy branches
give the corresponding anchor-free form.  This is an arity ledger, not a
claim that the live LRC construction has aligned all three fields.

There is a cycle-type version.  For three uniform states, let

```text
Y=pi(X),                   pi in S_3.                         (51)
```

Every permutation coupling has the same two endpoint banks.  The identity
has three fixed loops, a transposition has one, and a three-cycle has none.
Thus no separately integrated boundary/Gram quotient can determine an
underlying deterministic transport graph's loop type universally.  These
are abstract couplings; no semantic LRC owner loop has been constructed.

Equations (40) and (47) give the exact repair rank.  The cross-Kakeya bank
records `C` itself and hence every mixed-Haar coordinate.  A coupling
statistic can be inferred from separate endpoint banks only if its linear
functional annihilates the zero-margin kernel.

## 8. Exact gain and remaining obligation

The new boundary ledger is

```text
THM-2540 uncorrected transverse product shift
  -> loses the deep target anchor;

deep-target-covariant shift
  -> anchor + three gains + 216 incidences;

constant-singleton hostile
  -> those data still do not move deep-mask location;

same-horizon empty head
  -> annihilates as its own arrival;

genuine later target-active field on the same ancestry base
  -> complete positive source/arrival coupling by cross-Kakeya products;

separate endpoint integration
  -> loses exactly the mixed-Haar transportation kernel.     (52)
```

THM-2537 has already closed the same-predecessor phase-to-selected-positive-
wall seam; repeating that scalarisation cannot prove its equation (56).
This theorem refines THM-2536's triangle-holonomy invoice.  Vertex divergence and
separate endpoint margins forget cycle data for the same reason: they kill
the zero-margin coupling sector.  The next lawful proof object is not
another root colour, gain, marker, or target anchor.  It must construct a
positive same-ancestry realization and chronological typing of the formal
later Boolean target-active/owner field, or an equivalent polarized
cross-Gram field, and prove that its required mixed coupling functional is
nonzero while retaining the marked source/deep coordinate.

THM-2535's selector/scheduler intersections are genuine common-base events,
but their transported incidence cycle is not yet the inherited THM-2449
source/deep defect and its negative endpoint is still not an arrival.
THM-2533's owner is genuine and late, but current canon supplies no Boolean
root field `a_s` identifying that owner as the destination in (40).

No semantic source-to-arrival edge, owner loop, scalar-row exclusion, or
proof of LRC(14) is claimed.

## 9. Exact referee

Run

```bash
python3 04-computation/lrc14_anchored_transverse_arrival_boundary_thm2538.py
python3 -O 04-computation/lrc14_anchored_transverse_arrival_boundary_thm2538.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_anchored_transverse_arrival_boundary_thm2538.out
```

byte-for-byte.  The dependency-free referee verifies the deep-anchor
covariance, `4,826,809` full-character multiplier cases, every three-row
transverse and four-scalar avoidance pattern, the `216`-incidence orbit
count, all atoms and primitive modes of the constant-singleton hostile, the
complete same-horizon annihilation census, all retained-anchor and Crofton
pointwise reconstructions, all `106,496` all-slope empty-head hitting basis
cases, categorical selected-head and positive threshold-layer identities,
the two-state Haar switch, all `284` finite Frechet event pairs, `5,778`
exact Hall transportation problems, all six three-state permutation
couplings, the even/odd three-field parity hostile, and the rank-`23`
marginal map with its `121` mixed-Haar kernel.  In total it performs
`7,140,902` exact checks.

An independent split audit checked Sections 1--4 (lawful atoms, anchor
covariance, positive packing, `216` incidences, and the `c=1` hostile) and
Sections 5--8 (same-horizon annihilation, cross-Kakeya reconstruction,
transportation kernel, Hall bounds, and all stated non-consequences).  Both
halves passed without correction.

**QED.**
