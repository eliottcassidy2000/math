---
id: THM-2534
title: "Positive Kakeya boundary transform and Crofton reconstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every nonconstant
  Boolean mask on F_13, the all-slope field K_tau(r)=f_r(1-f_(r+tau)) is
  exactly the directed complete cut from the occupied roots to the empty
  roots.  Its matched opposite slopes are the disjoint positive and negative
  parts of the ordinary translation gradient (equivalently of
  (I+P_tau)C_tau), its divergence reconstructs the centered
  mask, and its outgoing/incoming packets positively reconstruct the mask;
  the two constant masks are the sole pointwise collision.  The scalar
  directional masses form a strict conditionally-negative Cayley variogram:
  their nontrivial Fourier coefficients are -|fhat(k)|^2 and their total is
  the exact Crofton variance n(13-n).  Each slope has at most six boundary
  roots and is reconstructed from at most seven Newton/Prony samples.  For a
  positive weighted mask law the boundary matrix is H_r-G_(r,s); its skew is
  the exact gradient H_r-H_s and its symmetric part is the Boolean cut
  semimetric.  Without an excluded root it determines (H,G) only up to one
  constant gauge, sharply.  With target root zero excluded it reconstructs
  H and the whole Gram matrix.  On the singleton/adjacent-pair deep path cone
  the twelve slopes collapse exactly to the all-occupied vector H and the
  two start/terminal selector vectors, which jointly recover all 23 ray
  masses.  A rational target anchor forces every nontrivial root colour in
  each fixed slope, but not general two-dimensional mixed modes: an exact
  anchored positive mixture has 132 such zeros.  The transform contains only
  same-time occupied-to-empty spatial boundaries.  It intrinsically supplies
  no C_7 clock, source charge, semantic arrival, future-owner identification,
  target-cell drift, row exclusion, or proof of LRC(14).
source: codex-2026-07-27-positive-kakeya-boundary-transform
depends_on:
  - THM-2530-anchored-deep-gram-cone-and-lossless-skew-target-refinement
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2439-cyclic-marker-replica-degree-and-homometric-gram-boundary
  - THM-2531-prime-necklace-guard-boundary-selector
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
  - THM-2536-deep-comb-selector-flow-target-drift-and-centered-toothpick-boundary
script: 04-computation/lrc14_positive_kakeya_boundary_transform_codex_20260727.py
output: 05-knowledge/results/lrc14_positive_kakeya_boundary_transform_codex_20260727.out
script_sha256: 45a037bb9e747bef2e77bf1311d0f417ab2db6086e6f2b1bb1ed8efe2a2dff9d
output_sha256: 45bfe3277dae41adced5115d281605395a734a38f961c819a3c43254a2e2f523
hash_basis: working-tree bytes (LF)
---

# THM-2534 -- the all-slope Boolean boundary is an oriented cut transform

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The word `Kakeya` here refers only to retaining every nonzero Cayley direction
on `F_13`.  No finite-field Kakeya theorem is invoked.  The underlying object
is simpler and more rigid:

```text
Boolean root mask
  -> directed complete cut occupied -> empty
  -> positive Jordan lift of every translation gradient `(I+P_tau)C_tau`
  -> strict cyclic variogram / Crofton energy
  -> target-anchored first-and-second-moment reconstruction.              (1)
```

This view removes a root-location and Gram-reconstruction problem.  It does
not turn a spatially empty root into a temporal arrival.

## 1. The boundary field is the complete directed cut

Put

```text
X=F_13,                       f:X->{0,1},
n=sum_r f_r,                  1<=n<=12.                       (2)
```

For every `tau in X`, including `tau=0`, define

```text
mathcal K_tau(r)=f_r(1-f_(r+tau)).                            (3)
```

The zero-slope field vanishes.  Change from tail/direction coordinates to
tail/head coordinates by

```text
M(r,s)=mathcal K_(s-r)(r)=f_r(1-f_s).                        (4)
```

Thus `M` is precisely the indicator of

```text
supp(f) x (X\supp(f)).                                       (5)
```

Pointwise it is the rank-one outer product `f(1-f)^T`, with zero diagonal.
As a binary relation, `r->s` exactly when `f_r>f_s`.  This is the strict part
of a two-level total preorder.  It has honest ties within the occupied class
and within the empty class, so it is **not** a tournament.  THM-2531's word
tournament is a tie-breaking refinement; (4) is the coarser intrinsic cut.

## 2. Matched slopes are the Jordan parts of the translation gradient

The edge opposite `(r,r+tau)` is represented by

```text
mathcal K_(-tau)(r+tau)=f_(r+tau)(1-f_r).                    (6)
```

The two terms are disjoint and satisfy

```text
mathcal K_tau(r) mathcal K_(-tau)(r+tau)=0,

mathcal K_tau(r)-mathcal K_(-tau)(r+tau)
 =f_r-f_(r+tau),

mathcal K_tau(r)+mathcal K_(-tau)(r+tau)
 =(f_r-f_(r+tau))^2.                                        (7)
```

Therefore the matched pair is exactly the positive/negative Jordan
decomposition of one discrete translation derivative.  With THM-2532's
convention `(P_tau f)_r=f_(r+tau)`, put

```text
B^+_tau=mathcal K_tau,
B^-_tau(r)=mathcal K_(-tau)(r+tau).
```

Then (7) and the Cayley law give the exact bridge

```text
B^-_tau-B^+_tau
 =(P_tau-I)f
 =(I+P_tau)C_tau f.                                        (7a)
```

Because `13` is odd, `I+P_tau` is invertible.  Thus the matched positive
walls recover `C_tau f`, but they are Jordan parts of the ordinary gradient,
not of `C_tau f` itself.  Summing all directions gives the complete Cayley
Laplacian:

```text
sum_(tau!=0)[mathcal K_tau(r)-mathcal K_(-tau)(r+tau)]
 =13f_r-n.                                                   (8)
```

There is also a wholly positive reconstruction.  Let

```text
O_r=sum_(tau!=0)mathcal K_tau(r),
I_r=sum_(tau!=0)mathcal K_tau(r-tau).                        (9)
```

The outgoing and incoming packets are

```text
O_r=(13-n)f_r,
I_r=n(1-f_r),
O_r I_r=0.                                                   (10)
```

Hence

```text
f_r=1 iff O_r>0,                    1-f_r=1 iff I_r>0.       (11)
```

The two constant masks both give the zero boundary field, and (11) proves
that they are the only pointwise collision.  A single matched pair of slopes
already gives every difference along one thirteen-cycle and hence recovers a
mixed mask; the full bank additionally retains every pairwise cut edge.

## 3. Affine covariance and exact converse

For `u!=0,b in X`, let

```text
g_x=f_(u^(-1)(x-b)).                                         (12)
```

Then direct substitution gives

```text
mathcal K^g_(u tau)(u r+b)=mathcal K^f_tau(r).                (13)
```

Thus the slope must be transported with the root chart.  In particular a
reflection `J_b:r->b-r` obeys

```text
mathcal K^(J_b f)_tau(r)=mathcal K^f_(-tau)(b-r).             (14)
```

Complement is a different involution and is the exact global converse:

```text
M_(1-f)(r,s)=M_f(s,r),

mathcal K^(1-f)_tau(r)=mathcal K^f_(-tau)(r+tau).             (15)
```

Thus root reflection, fixed-mask guard reversal, and semantic converse are
not conflated.  Complement reverses every actual cut arrow.

## 4. Crofton variance, phase loss, and the full Fourier lift

Let

```text
c_tau=sum_r f_r f_(r+tau),
b_tau=sum_r mathcal K_tau(r)=n-c_tau.                         (16)
```

Then `b_0=0`, `b_tau=b_(-tau)`, and for every `tau!=0`,

```text
b_tau=number of occupied runs in the tau-cycle
     =1/2 sum_r(f_r-f_(r+tau))^2 >0.                         (17)
```

Every ordered occupied/empty pair has one difference, so

```text
sum_tau b_tau=n(13-n)=13^2 Var_X(f),                         (18)

Var_X(f)=1/13 sum_r(f_r-n/13)^2.
```

This is the exact discrete Crofton variance formula.

Let `zeta=exp(2 pi i/13)` and use the unnormalized transform

```text
fhat(k)=sum_r f_r zeta^(-kr),
bhat(k)=sum_tau b_tau zeta^(-k tau).                         (19)
```

Wiener--Khinchin applied to (16) gives

```text
bhat(0)=n(13-n),
bhat(k)=-|fhat(k)|^2<0,                   k!=0.               (20)
```

Strictness follows without an estimate.  If one `fhat(k)` vanished, the
rational Boolean polynomial `sum_r f_r X^r` would vanish at a primitive
thirteenth root.  Divisibility by

```text
Phi_13(X)=1+X+...+X^12
```

would make all thirteen bits equal.  Thus `b` is a strict conditionally-
negative translation-invariant Cayley variogram.  Equivalently,

```text
b_tau=1/2 ||f-P_tau f||_2^2.                                 (21)
```

The scalar variogram retains the power spectrum but loses Fourier phase.
Write the complete field transform as

```text
Khat(a,b)
 =sum_(r,tau)mathcal K_tau(r)zeta^(-ar-btau).                (22)
```

With `s=r+tau`, (4) factors it exactly:

```text
Khat(a,b)=fhat(a-b) (13*1_(b=0)-fhat(b)).                    (23)
```

Every one of the `169` coefficients in (23) is nonzero for one mixed
Boolean mask.  Equation (20) is the positive anti-diagonal `a=0` of this
factorization.

The directional dispersion has a second exact form.  If

```text
E(A)=sum_tau c_tau^2
```

is the additive energy of `A=supp(f)`, then

```text
sum_(tau!=0)
 (b_tau-n(13-n)/12)^2
 =E(A)-n^2-(n^2-n)^2/12.                                    (24)
```

It vanishes exactly when every nonzero difference has one common
multiplicity.  The exhaustive finite census is

```text
n=1:  13 masks,      n=4: 52 masks,
n=9:  52 masks,      n=12:13 masks,                          (25)
```

forming `1+4+4+1=10` rotation necklaces.  The size-four and size-nine
layers are the `(13,4,1)` difference sets and their complements.

This scalar boundary is sharply phase-blind.  THM-2439's homometric pair

```text
A={0,1,3,9},                 B={1,2,5,7}                     (26)
```

has

```text
b_tau(A)=b_tau(B)=3                   for every tau!=0,       (27)
```

while `M_A!=M_B`.  The full directed field restores the mask; the Crofton
profile alone does not.

## 5. Run-terminal sparsity and the finite Prony invoice

Because every nonzero `tau` generates the thirteen-cycle, the support

```text
B_tau={r:f_r=1, f_(r+tau)=0}                                 (28)
```

is the set of terminal points of the occupied runs in that order.  Hence

```text
1<=|B_tau|=b_tau<=min(n,13-n)<=6.                            (29)
```

Put `x_r=zeta^r`.  The DC value gives `b=|B_tau|`, and the next `b` power
sums

```text
p_j=sum_(r in B_tau)x_r^j,                 1<=j<=b,           (30)
```

give the elementary symmetric functions by Newton's identities.  They
therefore reconstruct

```text
prod_(r in B_tau)(X-x_r)                                    (31)
```

and the complete boundary support.  Thus at most seven samples, including
DC, reconstruct one pointwise Boolean slope.  For an arbitrary weighted
`b`-atomic slope, where the atom weights are also unknown, the ordinary
Prony invoice is `2b` moments instead.  A single slope still does not recover
run lengths: the masks `{0}` and `{9,10,11,12,0}` have the same slope-one
boundary `{0}`.

## 6. Positive weighted laws and the one missing constant gauge

Now let `F>=0` weight a finite or measurable family of Boolean masks and
integrate every displayed pointwise quantity.  Put

```text
H_r=integral F f_r,
G_(r,s)=integral F f_r f_s,
M_(r,s)=integral F f_r(1-f_s).                               (32)
```

Then

```text
M_(r,s)=H_r-G_(r,s),

M_(r,s)-M_(s,r)=H_r-H_s,

M_(r,s)+M_(s,r)
 =integral F(f_r-f_s)^2.                                    (33)
```

Thus the skew part of the positive boundary matrix is the exact gradient
of `H`; every triangle circulation is zero.  The symmetric part is an `L^1`
cut semimetric.  In compact form,

```text
M_(r,s)=1/2[d_(r,s)+H_r-H_s],
d_(r,s)=integral F(f_r-f_s)^2.                               (34)
```

The divergence is

```text
D_r=sum_s(M_(r,s)-M_(s,r))=13H_r-sum_s H_s.                 (35)
```

Consequently an unanchored boundary law determines `H` only up to an
additive constant.  At the first/second-moment level the exact gauge is

```text
(H,G)->(H+c 1,G+c 1 1^T),                                   (36)
```

which leaves `M` fixed.  This is not a formal defect.  The uniform laws on
the thirteen singleton masks and on the thirteen co-singleton masks both
give

```text
M_(r,s)=1/13                    if r!=s,
M_(r,r)=0,                                                     (37)
```

but their occupancy vectors are respectively the constants `1/13` and
`12/13`.  Pointwise injectivity has been destroyed by averaging.

### The excluded target fixes the gauge and reconstructs the Gram matrix

Suppose now that

```text
f_0=0                         on every weighted fibre.        (38)
```

Then root zero is a common empty target and

```text
H_r=M_(r,0)=mathcal K_(-r)(r),
G_(r,s)=H_r-M_(r,s).                                         (39)
```

Thus the target-anchored all-slope boundary bank is algebraically equivalent
to the complete pair `(H,G)`.  In particular it contains the old
all-occupied-root star `H_a=G_(a,a)` as its target column.  It is a positive
replacement for storing `H-G` entry by entry; subtraction is used only when
recovering `G`.

The ordinary skew `M-M^T` is **not** THM-2530's lossless Hilbert
anticommutator.  It has the gradient form in (33) and remembers only `H`.
Losslessness here uses the full positive matrix together with the target
anchor.

Two spectral statements survive averaging, with different hypotheses.

1. If positive weight lies on nonconstant masks, then for every `b!=0`,

   ```text
   Khat(0,b)=-integral F |fhat(b)|^2<0.                       (40)
   ```

   Thus all twelve global slope characters survive without rationality or a
   target anchor.

2. If the weights are rational and (38) holds, then for each fixed
   `tau!=0` the rational vector `mathcal K_tau(r)` has coordinate zero at
   `r=0` and positive total.  It cannot be constant.  Irreducibility of
   `Phi_13` therefore forces all twelve of its nontrivial root Fourier
   coefficients to be nonzero.  This gives `12*12=144` fixed-slope/root
   coefficients.

Neither assertion forces all mixed two-dimensional modes after integration.
Take the uniform law on the `2^12-1` nonempty subsets of `{1,...,12}`.  It is
rational and target-anchored, and

```text
M_(r,0)=2048/4095                         for r!=0,
M_(r,s)=1024/4095                         for distinct r,s!=0,
M_(0,s)=M_(r,r)=0.                                           (41)
```

Its endpoint Fourier coefficient at `(k,l)` vanishes exactly when

```text
k!=0,               l!=0,               k+l!=0.             (42)
```

There are `12*11=132` such zeros.  Only the twelve mixed anti-diagonal
energy modes are protected in that quadrant.  This is the sharp distinction
between fixed-slope root support and a general joint slope/root transform.

## 7. On the deep path cone the bank collapses to three vectors

Use THM-2530's target-anchored singleton/adjacent-pair cone.  In one cell,
let

```text
alpha_j >=0,                         1<=j<=12,
beta_j >=0,                          1<=j<=11,
beta_0=beta_12=0                                                  (43)
```

be the masses of `{j}` and `{j,j+1}`.  Then

```text
H_j=alpha_j+beta_(j-1)+beta_j.                               (44)
```

The two natural boundary slices are exactly THM-2536's selector mass
vectors:

```text
gamma^+_j=M_(j,j-1)=mathcal K_(-1)(j)=alpha_j+beta_j,

gamma^-_j=M_(j,j+1)=mathcal K_1(j)=alpha_j+beta_(j-1).        (45)
```

If `tau` is neither `1` nor `-1`, no path ray occupies both `j` and
`j+tau`.  Therefore

```text
mathcal K_tau(j)=H_j,                    tau notin {1,-1}.    (46)
```

So the entire twelve-slope bank contains exactly three distinct root
vectors:

```text
(H,gamma^+,gamma^-).                                         (47)
```

They are jointly lossless for all `23` ray masses:

```text
beta_j=H_j-gamma^-_j,                         1<=j<=11,

alpha_j=H_j-beta_(j-1)-beta_j.                              (48)
```

Equivalently `beta_j=H_(j+1)-gamma^+_(j+1)`.  If

```text
A=sum_j alpha_j,                 B=sum_j beta_j,
N=sum_j H_j=A+2B,
T=sum_(r,s)M_(r,s)=12A+22B,                                (49)
```

then the total cell mass is `(T-10N)/2=A+B`.

Equation (39) identifies `H` with the all-occupied-root star already present
in the one-deep algebra of THM-2365.  Equations (45) identify the start and
terminal selector **coefficients** with two adjacent boundary slices.  The
THM-2531 selector star then relocates such a coefficient to the chord
`j->0`; the actual edge represented by (45) is `j->j-1` or `j->j+1`.
Their difference is the target-triangle circulation recorded in THM-2536.

Thus the all-slope transform removes the need for a degree-twelve marker or
a separate Gram inversion on this restricted cone.  On an arbitrary mask a
slope may contain several domain walls, and the boundary field alone does
not choose THM-2531's lexicographically distinguished one before averaging.

## 8. Exact semantic boundary and source-charge hostile

Every positive entry of `M` says

```text
root r is occupied now,               root s is empty now.    (50)
```

It does not say that an event travels from `r` to `s`, that `s` is occupied
later, or that an owner arrives there.  Its skew is the curl-free spatial
gradient (33), not a temporal or source current.

THM-2533 may multiply every root coordinate by the same root-invariant late
owner/terminal-word factor, so an attached future owner is compatible with
the identities above.  That external factor does not identify either wall
endpoint with the future owner or make the empty head a chronological
arrival.

The `C_7` obstruction in THM-2535 is untouched.  A pure septimal CRT
translation fixes the complete root mask and hence every entry of this
transform while moving the owner clock freely.  No function of the
all-slope bank can be a translation-covariant choice of a seven-epoch clock.

The old source-charge hostile becomes stronger, not weaker.  Take

```text
f=1_{1},

d(h,kappa)
 =(1_(h=0)-1_(h=2))(1_(kappa=0)-1_(kappa=1)).                (51)
```

Then

```text
mathcal K_tau(1)=1                         for every tau!=0,

dtilde(alpha,beta)
 =(1-zeta_13^(-2alpha))(1-zeta_7^(-beta))!=0
                         for all alpha,beta!=0,               (52)

d(1,kappa)=0                              for every kappa.    (53)
```

Thus every root direction is a positive boundary, the pointwise boundary
field has its complete spectrum, and the signed source defect has all `72`
primitive modes, yet every boundary tail lies on an uncharged source row.
All-slope geometry does not create the missing joint ancestry law.

The exact gain and loss ledger is therefore:

```text
REMOVED
  pointwise root-mask ambiguity (apart from the two constants);
  target-anchored first/second-moment and Gram reconstruction;
  root-colour loss in each rational anchored slope;
  separate H/gamma+/gamma- bookkeeping on the deep path cone.

REMAINING
  the averaging gauge without a common excluded root;
  arbitrary mixed slope/root Fourier cancellation;
  a canonical wall on general multi-run masks after quotient;
  the C_7 clock and inherited source/deep charge;
  target-cell variation, semantic arrival, intrinsic future-owner meaning,
  and row exclusion.
                                                                    (54)
```

No scalar speed row is excluded and LRC(14) remains open.

## 9. Exact companion

Run

```bash
python3 04-computation/lrc14_positive_kakeya_boundary_transform_codex_20260727.py
python3 -O 04-computation/lrc14_positive_kakeya_boundary_transform_codex_20260727.py
```

The dependency-free referee exhausts all `8,190` mixed Boolean masks.  It
checks `1,277,640` matched oriented-edge identities, `3,832,920` affine-
generator covariance identities, all `98,280` sparse Newton reconstructions,
and all `98,280` exact cyclotomic Crofton identities.  It also certifies the
directional-variance census, the fixed-slope and homometric phase hostiles,
the singleton/co-singleton gauge collision, the `132`-mode anchored hostile,
target-anchor reconstruction of `(H,G)`, the deep-path collapse and `23`-ray
inverse, and the all-slope zero-source-tooth hostile.  Normal and optimized
executions reproduce the stored transcript byte-for-byte. **QED.**

The independent audit reconstructed the endpoint matrix `M=f(1-f)^T`,
the Crofton/variogram transform, the anchored `(H,G)` inverse, and the
deep-cone three-vector collapse.  It also corrected the operator ledger:
the positive matched walls are the Jordan parts of `(P_tau-I)f`, which is
`(I+P_tau)C_tau f`, not of `C_tau f` itself.  Normal and optimized runs
byte-match the stored transcript and recorded hashes.
