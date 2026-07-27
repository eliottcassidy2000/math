---
id: THM-2556
title: "Reynolds duty curvature and the exact fibre-covariance mixed cell"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every finite scalar relation type, the all-unit pushforward splits
  exactly into a target-fibre duty multiplier plus a within-fibre covariance
  residual, and the seven-step translation square has an exact two-term
  curvature identity.  On the canonical THM-2541/2550 typed row the duty
  profile has all 168 nontrivial target characters; every nonzero gain leaves
  exactly 156 mixed-curvature characters.  A nonnegative nonunit-section
  hostile cancels that quotient curvature exactly.  The residual has an
  exact fibre-variance operator norm, but no lawful physical common-carrier
  intertwiner or sufficiently small physical covariance estimate is proved;
  no row is removed and LRC(14) remains OPEN.
source: codex-2026-07-27-holotopy
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2544-all-unit-target-projector-kernel-and-lawful-image-intersection-obstruction
  - THM-2548-seven-step-c91-transfer-and-full-norm-separation
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2541-canonical-typed-row-full-target-plane-support
  - THM-2550-canonical-typed-row-double-nondegeneracy
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
  - THM-2552-flat-common-base-gain-versus-unjoinable-wall-spectrum
script: 04-computation/lrc14_reynolds_duty_curvature_thm2556.py
output: 05-knowledge/results/lrc14_reynolds_duty_curvature_thm2556.out
script_sha256: 1ee1fcf3bf713b67fb9f5502cb9247517ce4a1b11da95ebf4bcdb3d95dfaa537
output_sha256: b14b861531eed46231189712bb2232942b70c9a6b52b37b94e79f7ca2e48f1c7
hash_basis: LF-normalized bytes
---

# THM-2556 -- the mixed square has an exact duty/covariance decomposition

THM-2551 proves that an untwisted product of the C91 transfer with the
all-unit projector preserves the projector kernel.  This theorem computes
the first non-product replacement: translate the relation-residue fibre and
ask how the all-unit mask changes.

The resulting curvature has two terms:

```text
quotient duty curvature + within-fibre covariance curvature.             (1)
```

The first is completely explicit and has full target spectrum on the
canonical row.  The second can cancel it exactly.  This identifies a concrete
norm inequality, rather than an unspecified gluing principle, as the next
possible positive bridge.

## 1. Reynolds decomposition on every scalar relation fibre

Work over `Q` (or `R,C`).  For any THM-2544 scalar type let

```text
Omega=K_91,                    Q:Omega->G=F_13^2,
F=|ker Q|,                     V=Q^Omega,
U:V->Q^G,                      Uc(q)=sum_(Qy=q)c(y),
S={y: every coordinate of y is a unit mod 91},
M c=1_S c,                     J=UM.                         (2)
```

Define the fibre-uniform lift and its projection by

```text
(La)(y)=a(Qy)/F,                Pi=LU.                       (3)
```

Then `UL=I`, `Pi^2=Pi`, and `Pi` projects onto the fibre-uniform currents.
Put

```text
nu(q)=|S intersection Q^(-1)(q)|/F,
M_nu a(q)=nu(q)a(q),
R=UM(I-Pi).                                                   (4)
```

Here `M_nu` is a positive Reynolds **multiplier**, not generally an
idempotent projector.  Splitting `c=Pi c+(I-Pi)c` gives the exact identity

```text
J=M_nu U+R.                                                  (5)
```

Equivalently, averaging the raw unit mask over neutral translations gives

```text
1/F sum_(h in ker Q) T_h M T_(-h)
  = multiplication by nu o Q.                              (6)
```

This is why Reynolds averaging is load-bearing.  A raw neutral translation
can cross a unit wall even though it has zero target gain; the averaged
multiplier commutes with every neutral translation by construction.

There is also an exact singular-value invoice.  Give both spaces their
counting `l2` norms and put

```text
rho_q=1_(S intersection Q^(-1)(q))-nu(q)1_(Q^(-1)(q)).      (6a)
```

The rows of `R` are the pairwise orthogonal vectors `rho_q`, so

```text
R R^*=diag_q(F nu(q)(1-nu(q))),
||R||_(2->2)^2=max_q F nu(q)(1-nu(q)).                      (6b)
```

Equality is attained by the normalized centred unit mask in any maximizing
fibre.  Thus the residual is controlled exactly by duty variance; small duty
alone is not the relevant quantity.

## 2. Exact mixed-square curvature

For `v in Omega`, write `g=Q(v)` and use

```text
(T_v c)(y)=c(y-v),             (tau_g a)(q)=a(q-g).          (7)
```

Then

```text
UT_v=tau_g U,                  T_vL=L tau_g.                 (8)
```

Define the seven-step translations

```text
D_v=sum_(j=0)^6 T_v^j,         C_g=sum_(j=0)^6 tau_g^j.      (9)
```

Using (5) and (8) once on each side gives the exact mixed-cell formula

```text
JD_v-C_gJ
 =[M_nu,C_g]U +(R D_v-C_gR).                               (10)
```

The first term is quotient duty curvature.  The second is the entire
within-fibre covariance debt.  No term is hidden.

On the uniform lift of the constant target vector,

```text
[M_nu,C_g]1(q)=7nu(q)-sum_(j=0)^6nu(q-jg).                  (11)
```

If a lawful common-carrier intertwiner identifies THM-2548's root transfer
with some `D_v`, then (10) is the desired mixed 2-cell on that carrier.
Canon has not made this identification.

## 3. Canonical-row duty profile

Now specialize only this section onward to the canonical typed row

```text
w=(1,14,27,40,53,66,13,2197,742586)                        (12)
```

and THM-2309's owner packet used by THM-2541/2550.  In canonical target
coordinates `q=(x,y)`, every point of the mod-thirteen packet fibre has the
form

```text
(-a1,-a2,-a3,-a4,-a5,
 a1+a2+a3+a4+a5,-a0,x-a1,y-a2).                            (13)
```

Consequently it is coordinatewise nonzero exactly when

```text
a0,...,a5!=0,
a1+...+a5!=0,
a1!=x,                      a2!=y.                          (14)
```

The standard nonzero-variable character count

```text
#{(b1,...,bk) in (F_13^*)^k:sum b_i=0}
 =(12^k+12(-1)^k)/13                                      (15)
```

and two-step inclusion--exclusion give the exact number of unit mod-thirteen
packet points:

```text
n(x,y)
 =2316060
  +210552(1_(x=0)+1_(y=0))
  +12 1_(x+y=0)
  +19128 1_(x=y=0).                                       (16)
```

Thus the histogram is

```text
2756304   once, at (0,0);
2526612   on the 24 punctured axis points;
2316072   on the 12 nonzero antidiagonal points;
2316060   on the remaining 132 points.                     (17)
```

The row has septimal support `s=8`, so THM-2544's exact count is

```text
N_7=6(6^8+6)/7=1439676.                                   (18)
```

Since

```text
F=|ker Q|=7^8 13^6=27825593350009,                         (19)
```

the full all-91-unit duty profile is

```text
nu(x,y)=N_7 n(x,y)/F.                                     (20)
```

Every value is positive, but positivity is not what drives the next result;
the nonconstancy and exact Fourier support do.

## 4. Full duty spectrum and 156-mode curvature

Use the unnormalized convention

```text
n_hat(u,v)=sum_(x,y)n(x,y)zeta^(-ux-vy).                   (21)
```

The four indicator terms in (16) give, for `(u,v)!=(0,0)`,

```text
n_hat(u,v)
 =19128
  +2737176(1_(u=0)+1_(v=0))
  +156 1_(u=v).                                            (22)
```

At zero, `n_hat(0,0)=396907776`.  In particular all `168` nontrivial
characters are positive and `nu` has trivial translation stabilizer.

For `g!=0`, the Fourier multiplier of `C_g` at `ell=(u,v)` is

```text
c_g(ell)=sum_(j=0)^6 zeta^(-j ell.g).                       (23)
```

Therefore (11) has multiplier

```text
(7-c_g(ell)) nu_hat(ell).                                  (24)
```

It vanishes exactly on the thirteen-character orthogonal line
`ell.g=0`.  If `ell.g!=0`, equality of the sum of seven unit complex numbers
with `7` would force every summand, hence `zeta^(-ell.g)`, to be one.  Thus
every nonzero gain has

```text
156 nonzero duty-curvature characters,
12 nonzero orthogonal characters killed,
1 constant character killed.                              (25)
```

The integer defect

```text
delta_g(q)=7n(q)-sum_(j=0)^6n(q-jg)                         (26)
```

has three exact squared energies:

```text
sum_q delta_g(q)^2 =
 24559042191264     for the 24 punctured-axis gains;
 49102678687104     for the 12 nonzero antidiagonal gains;
 49102698046752     for the remaining 132 gains.            (27)
```

Multiply (27) by `(N_7/F)^2` for the normalized duty profile.

## 5. Sharp cancellation hostiles

### Neutral raw-mask curvature is not target activity

On the toy carrier `Omega=F_13^2`, `Q(x,y)=x`, and
`S=(F_13^*)^2`, the neutral shift `(0,1)` has

```text
rank[M,T_(0,1)]=24.                                        (28)
```

It crosses unit walls but has zero target gain.  Reynolds averaging instead
gives `nu(0)=0`, `nu(x)=12/13` for `x!=0`, which commutes with neutral shifts.
Its one-step commutator with target shift one has rank two and all twelve
nontrivial target colours.  This proves that raw mask curvature is
semantically blind and that (6) cannot be omitted.

### The covariance residual can cancel the duty term exactly

In the same toy model, the nonnegative current

```text
c_host=1_(y=0)                                             (29)
```

has `Uc_host=1`, `Jc_host=0`, and is invariant under target translations.
Its full mixed curvature is zero even though the duty term on `U c_host` is
nonzero; the covariance term cancels it exactly.

The canonical version is THM-2544's

```text
c_host=sum_(q in G)delta_(v_q),                            (30)
```

where the nonunit section `q->v_q` is additive.  For any `g!=0`, put
`v_g=v_(g)-v_(0)` and

```text
H_g(k,y)=(k+1,y+v_g)                                      (31)
```

on `C_7 x Omega`.  This is one 91-cycle direction,
but the operators must be lifted before applying it to a current.  On
`Q^(C_7 x Omega)` put

```text
U~=I tensor U,       J~=I tensor J,       M_nu~=I tensor M_nu,
H~=S_clk tensor T_(v_g),                 C~=S_clk tensor tau_g, (31a)
```

and take the clock-uniform current `c~(k,y)=c_host(y)`.  Then `H~` has order
`91`, `H~c~=c~`, `U~c~=1`, and `J~c~=0`.  The fibrewise version of (10)
therefore has zero full mixed face while its quotient duty face (25) is
nonzero.  Positivity, full unrestricted target support, and C91 holonomy do
not control the covariance residual.  On `Omega` alone the corresponding
translation has order `13`; the clock lift in (31a) is load-bearing for the
literal C91 statement.

## 6. Exact conditional bridge

Suppose a lawful physical current `c` and a common-carrier intertwiner have
actually identified the relevant THM-2548 transfer with `D_v`, `Q(v)=g!=0`.
Then either of the following is sufficient for a nonzero all-unit mixed face:

1. `c` is fibre-uniform and its quotient duty term in (10) is nonzero; or
2. on the desired nonzero-target projection `pr_*`,

```text
||pr_*(R D_v-C_gR)c||_2
  < ||pr_*[M_nu,C_g]Uc||_2.                               (32)
```

The triangle inequality and (10) then force `pr_*(JD_v-C_gJ)c!=0`.
For the uniform constant lift on the canonical row, (25) supplies the full
156-mode signal explicitly.

The norm invoice (6b) makes the open hypothesis quantitative.  Since `Pi`
is the orthogonal fibre-average, it commutes with every `T_v`; moreover
`R Pi=0`, `||D_v||<=7`, and `||C_g||<=7`.  Hence

```text
||(R D_v-C_gR)c||_2
 <=14 sqrt(max_q F nu(q)(1-nu(q))) ||(I-Pi)c||_2.           (32a)
```

This is sharp at the one-row operator level but deliberately coarse after
the two seven-step sums.  It turns (32) into a concrete sufficient
near-fibre-uniformity test; it does not prove that a physical packet passes
that test.

The criterion is honest about what remains open.  The target-section
translation `v_g` and the algebraic 91-cycle exist, but no theorem identifies
them with the physical root `a`, ancestry sheet, semantic later edge, or a
lawful THM-2334 Abel current.  Nor is (32) known.  For arbitrary scalar types,
THM-2544 gives positive duty values when septimal support is at least two;
the full Fourier formula (22) is canonical-row-specific.

No row is excluded.  LRC(14) remains open.

## 7. Exact companion

The dependency-free companion reconstructs (16) both from its closed formula
and an independent `512`-subset affine inclusion--exclusion for every one of
the `169` target fibres.  It verifies the full cyclotomic Fourier spectrum by
exact exponent bins, all `28,392` gain/character curvature cases, the three
energies (27), every translation stabilizer, and the raw/Reynolds/hostile toy
boundaries.  Run

```bash
python3 04-computation/lrc14_reynolds_duty_curvature_thm2556.py
python3 -O 04-computation/lrc14_reynolds_duty_curvature_thm2556.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_reynolds_duty_curvature_thm2556.out
```

after LF normalization.  Every check raises explicitly under optimized
Python.

## 8. Independent hostile audit

An independent audit rederived the Reynolds decomposition, conjugacy laws,
mixed-cell identity, canonical affine count, full Fourier spectrum, all
`28,392` gain/character cases, and the three energy classes.  It caught and
repaired the order-`13` versus order-`91` typing boundary by requiring the
clock-fibre lift (31a).  It also independently derived the exact covariance
singular values (6b), checked the unique maximizing fibre and displayed norm
fraction, and verified the sharp nonunit-section cancellation.  Normal,
optimized, and stored transcripts agree exactly with `29,430` explicit
checks; the audit reproduced the LF hashes above.  No common physical
intertwiner or sufficiently small lawful covariance estimate is inferred.
