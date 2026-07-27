---
id: THM-2507
title: "Truncated Radon toothpick tomography and the nonaffine root boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For a
  rational row-zero array on F_p x {0,...,q-1}, q<p, the q-point
  toothpick pushforwards in any q-1 distinct nonzero slopes determine the
  array modulo the exact (q-1)-dimensional h-independent kernel. Every
  nonvertical array is visible in at least p-q+1 of the p-1 nonzero
  slopes, and this fraction is sharp in the abstract integral class. For
  the THM-2436 13 x 7 defect, at least seven of twelve slopes survive at
  almost every essential parent; every surviving pushforward has all twelve
  nonzero C_13 characters and coefficient floor 18^(-11). One fixed
  slope therefore survives on parent mass at least 1/6 or 1/4 in the two
  residual shapes. The fold is deliberately nonaffine: it uses an ordered
  representative cut on F_7, and THM-2436's carry translates the cut but
  does not sweep the slope. No standard thirteen-root, target, owner, or
  deep current is produced; physicalization requires the carry/cut sheet
  as an additional sidecar. No live row is excluded and LRC(14) remains
  open.
source: codex-2026-07-27-truncated-radon-toothpick
depends_on:
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
related:
  - THM-2368-owner-pivot-root-fibre-radon-invertibility
  - THM-2400-clean-parent-root-gauge-quotient-and-target-slope-boundary
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
script: 04-computation/lrc14_truncated_radon_escape_probe.py
output: 05-knowledge/results/lrc14_truncated_radon_escape_probe.out
script_sha256: 7dc0e33a090905e05babeb82efeeed954e43fff10b638514ae2e9c7f8cbb3e9b
output_sha256: fbf53503365231770e8a2b94c9b1727f2c513798764fcca393c7f38cb08b4393
hash_basis: working-tree bytes (LF)
---

# THM-2507 -- six toothpick directions detect the punctured stalk

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2506 proves a sharp negative statement: every affine homomorphic
pushforward of the punctured `C_13 x C_7` defect to a `13`-group vanishes.
The first escape is not a larger Fourier bank.  It is a smaller, explicitly
nonaffine family of seven-point fibres.

Choose the ordered representative section

```text
iota:F_7 -> F_13,             iota(r) in {0,1,...,6},            (1)
```

and, for `tau in F_13`, fold a defect along the toothpicks

```text
pi_tau(h,r)=h+tau iota(r),

R_tau(v)=sum_(r=0)^6 d(v-tau r,r).                              (2)
```

The fibres of `pi_tau` contain one point in every one of the seven rows.
They are finite Kakeya needles in a `13 x 7` strip, but not lines in the
group `C_13 x C_7`.  This failure of group-linearity is exactly why they can
see a defect that every affine `13`-primary quotient kills.

The result has three levels:

```text
abstract p x q strip:
  q-1 nonzero directions reconstruct modulo the vertical kernel;

THM-2436 essential 13 x 7 defect:
  at least 7 of the 12 nonzero directions survive pointwise;

physical scalar packet:
  OPEN -- retain the normalization carry, ordered cut, and sheet before
  interpreting a surviving fold as a target/owner/deep current.             (3)
```

## 1. The general truncated-Radon theorem

Let `p` be prime and let

```text
2<=q<p.
```

Write `I_q={0,...,q-1}`.  On the rational vector space

```text
V_(p,q)
 ={d:F_p x I_q -> Q:
     sum_(r=0)^(q-1)d(h,r)=0 for every h},                         (4)
```

define, for every `tau in F_p`,

```text
R_tau d(v)=sum_(r=0)^(q-1)d(v-tau r,r).                            (5)
```

For any set `S subset F_p^*` of `m` distinct nonzero slopes, the joint map

```text
R_S:V_(p,q) -> direct_sum_(tau in S) Q^(F_p)
```

has exact rank

```text
rank_Q R_S=(p-1)min(m,q-1).                                      (6)
```

In particular, as soon as `|S|>=q-1`, its kernel is exactly

```text
ker R_S
 ={d(h,r)=b(r):sum_r b(r)=0},

dim ker R_S=q-1.                                                  (7)
```

Thus any `q-1` distinct nonzero slopes detect whether `d` is nonvertical.
Equivalently, every nonvertical rational row-zero array satisfies

```text
#{tau in F_p^*:R_tau d!=0}>=p-q+1.                               (8)
```

### Proof

Fix a primitive `p`-th root `zeta`.  For `alpha in F_p`, put

```text
D_alpha(r)=sum_(h in F_p)d(h,r)zeta^(-alpha h),

P_alpha(X)=sum_(r=0)^(q-1)D_alpha(r)X^r.                          (9)
```

Fourier transform in the output variable of (5) gives the exact slice law

```text
widehat(R_tau d)(alpha)
 =P_alpha(zeta^(-alpha tau)).                                   (10)
```

The row-zero law in (4) says

```text
P_alpha(1)=0                                                    (11)
```

for every `alpha`.  At `alpha=0`, every output in (10) is zero and the
remaining `q-1` input dimensions are precisely the vertical arrays in (7).

For each `alpha!=0`, the allowed polynomials have degree at most `q-1`,
vanish at one, and form a space of dimension `q-1`.  Evaluation at `m`
distinct points

```text
zeta^(-alpha tau),               tau in S,
```

has rank `min(m,q-1)`: after dividing by `X-1`, this is an ordinary
Vandermonde evaluation on polynomials of degree at most `q-2`.  There are
`p-1` nonzero values of `alpha`, proving (6) after extending scalars; the
rank is unchanged from `Q` to `Q(zeta)`.

For `m=q-1`, all nonzero `alpha` slices vanish only when every column of
`d` is independent of `h`, proving (7).  If a nonvertical array had `q-1`
bad nonzero slopes, those slopes would put it in (7), a contradiction.
There are `p-1` nonzero slopes, so (8) follows. QED.

The rank increments in exact blocks of `p-1` until the `(q-1)`-st tooth,
then stops at the vertical kernel:

```text
0, p-1, 2(p-1), ..., (q-1)(p-1).                               (12)
```

This is the toothpick-ladder self-similarity.  Each new direction tests one
new evaluation of the same degree-`q-2` quotient polynomial in every
nontrivial horizontal character.

## 2. Galois saturation and sharpness

For rational `d`, one surviving output is automatically saturated in every
nonzero `C_p` colour.  Indeed `R_tau d` is rational and has total sum zero.
If one nonzero Fourier coefficient vanished, its rational coefficient
polynomial of degree at most `p-1` would be divisible by

```text
Phi_p(X)=1+X+...+X^(p-1).
```

Its entries would all be equal; total sum zero would then make the entire
output zero.  Therefore

```text
R_tau d!=0
  iff
widehat(R_tau d)(alpha)!=0 for every alpha in F_p^*.             (13)
```

The count (8) is sharp in the abstract integral row-zero class.  Choose any

```text
S subset F_p^*,                    |S|=q-2,
```

and put

```text
P(X)=(X-1)product_(tau in S)(X-zeta^(-tau))
    =sum_(r=0)^(q-1)D_rX^r.                                     (14)
```

Define the integral trace array

```text
d(h,r)=Tr_(Q(zeta)/Q)(D_r zeta^h).                               (15)
```

The `D_r` are cyclotomic integers, so (15) is integral.  Equation (14) at
`X=1` gives the row sums zero.  Moreover

```text
sum_h d(h,r)zeta^(-alpha h)=p sigma_alpha(D_r),
```

where `sigma_alpha(zeta)=zeta^alpha`.  Hence (10) vanishes exactly for

```text
tau in {0} union S.                                              (16)
```

There are exactly `p-q+1` good nonzero slopes.  This witness is a sharpness
control for the abstract theorem; it is not asserted to be one of the much
smaller THM-2436 cover defects or to obey their `L^1<=18` invoice.

The favorable strict inequality `q<p` is structural.  At `q=p` the lower
bound degenerates to one good direction.  Once `q>=p+1`, a nonzero multiple
of `X^p-1` fits inside the degree window and can vanish on every `p`-th root;
the truncated-direction detector is then no longer injective modulo the
vertical kernel.

## 3. The punctured `13 x 7` stalk

Return to the integral THM-2436 defect

```text
d_Y:F_13 x F_7 -> Z,

sum_r d_Y(h,r)=0.                                                (17)
```

THM-2436 proves

```text
d_Y is h-independent iff d_Y=0,                                 (18)
```

including one source, two distinct sources, and coincident labelled
sources.  It identifies `d_Y=0` with the flat locus `mathcal E`.  Apply
(8) with `(p,q)=(13,7)`.  At almost every essential parent

```text
Y in P minus mathcal E,
```

the twelve nonzero slopes obey

```text
#{tau in F_13^*:R_(Y,tau)!=0}>=7.                               (19)
```

By (13), every one of those good toothpicks has all twelve nonzero
`C_13` Fourier colours.  The affine slope `tau=0` is always zero:

```text
R_(Y,0)(v)=sum_r d_Y(v,r)=0.                                    (20)
```

Thus the first detector after the affine no-go is necessarily nonaffine.

There is a pointwise coefficient floor.  THM-2506 proves

```text
||d_Y||_1<=18.                                                   (21)
```

Pushforward cannot increase `L^1`, so every conjugate of a good
unnormalized coefficient has modulus at most `18`.  Such a coefficient is
a nonzero algebraic integer in the degree-twelve field `Q(zeta_13)`.  Its
field norm is a nonzero integer, hence for every `alpha!=0`,

```text
|widehat(R_(Y,tau))(alpha)|>=18^(-11).                           (22)
```

For the normalized `C_13` transform, divide the right side by `13`.

## 4. One slope can be fixed before the parent

Let

```text
E_tau={Y in P minus mathcal E:R_(Y,tau)!=0}.
```

Integrating (19) gives

```text
sum_(tau in F_13^*)mu(E_tau)>=7 mu(P minus mathcal E).            (23)
```

Consequently one fixed nonzero slope, selected before `Y`, satisfies

```text
mu(E_tau)>=(7/12)mu(P minus mathcal E).                           (24)
```

THM-2435's essential-parent floor is `(1+k_shape)/7`.  Therefore the two
typed residual shapes have the exact inherited fixed-slope floors

```text
k_shape=1:             mu(E_tau)>=1/6,

k_shape=2:             mu(E_tau)>=1/4.                           (25)
```

The same fixed slope carries every prescribed `alpha!=0` on its entire
locus, and (22) gives unnormalized fixed-colour energy floors

```text
k_shape=1:    integral_P |widehat(R_tau)(alpha)|^2 >=(1/6)18^(-22),

k_shape=2:    integral_P |widehat(R_tau)(alpha)|^2 >=(1/4)18^(-22). (26)
```

No union over the twelve `C_13` characters is paid.  The only remaining
selection is one of the twelve geometric toothpick slopes.

## 5. Why this evades the affine no-go but is not yet physical

The map `pi_tau` for `tau!=0` is not an affine homomorphism

```text
C_13 x C_7 -> C_13.
```

There is no nonzero homomorphism `C_7->C_13`.  More concretely, the section
`iota` has a seam.  Under the exact THM-2436 affine guard gauge, write the
CRT action as

```text
(h,r)->(a_13 h+c_13, a_7 r+c_7).                                (27)
```

For integer representatives `A,C` of `a_7,c_7`,

```text
iota(a_7r+c_7 mod 7)
 =A iota(r)+C-7 q_(A,C)(r),                                    (28)
```

where the wrap count `q_(A,C)(r)` depends on the row.  Pulling back the
toothpick phase therefore produces the extra cut term

```text
-7 tau q_(A,C)(r)                 mod 13.                        (29)
```

It is not a common root translation.  THM-2436's `kappa(Y)` is the
translation/carry part `c` of (27); it moves the cut and the absolute phase,
but it does not vary or sweep `tau`.

A literal sheet shear

```text
(h,r)->(h+sigma iota(r),r)                                     (30)
```

would send one toothpick slope to another.  But (30) depends on the chosen
representatives: crossing `r=6` changes its lift by `7 sigma`, which is
nonzero modulo thirteen for `sigma!=0`.  Thus it does not descend through
the `C_7` quotient and is not one of the affine root operations already
proved lawful.

Once the normalized root chart, the ordered `F_7` cut, the carry
`kappa(Y)`, and the full sheet `r` are retained, every `pi_tau` is a perfectly
well-defined finite set partition of the natural-extension stalk.  What is
**OPEN** is the next map:

```text
nonzero signed toothpick coefficient
  + carry/cut/sheet sidecar
  -> one Boolean owner-supported target/deep current on the physical circle.
                                                                    (31)
```

No current theorem proves that (30) preserves the scalar masks, the
THM-2365 target action, the THM-2471 owner/collision root, or the old deep
probe.  Calling `tau` a physical target slope would conflate it with the
factor slopes `eta_i/w_i` in THM-2400, which generally deform the packet and
are not this sheet shear.

## 6. Relation to the existing Radon and Kakeya views

The result is new at the exact interface left by THM-2506:

- THM-2368's Radon table is a cyclic correlation entirely on `C_13`; its
  pointwise all-mode theorem can still integrate to a circulant hostile.
- THM-2400 distinguishes harmless common-root gauge from unequal physical
  target slopes.  The toothpick parameter is neither one.
- THM-2506 proves that every affine `13`-primary pushforward kills the
  `C_13 x C_7` defect.  Equation (20) is that killed affine direction;
  equations (19) and (29) identify the exact nonaffine escape and its cost.
- THM-2504 warns that an oriented root chart is a consumed sidecar.  Here the
  ordered `F_7` cut and its wrap cocycle are consumed in addition.

This is a finite Kakeya statement in the precise sense that one seven-point
needle is tested in each slope.  It is not an import of Euclidean Kakeya
dimension theory.  Its proof is polynomial interpolation on a truncated
strip, and its destroyed information is explicit.

## 7. Exact gain and stopping boundary

The proved candidate chain is

```text
essential punctured-stalk defect
  -> at least seven of twelve nonaffine toothpick slopes
  -> all twelve C_13 colours on every good slope
  -> one fixed slope on parent mass 1/6 or 1/4.                    (32)
```

This improves the categorical conclusion of THM-2506: discarding `F_7`
affinely is fatal, but using it as an ordered sheet to bend the root chart is
injective modulo the already-classified vertical kernel.  It does not
transplant the now-empty high-septimal packet to the live `165` rows, make a
signed defect positive, identify an owner, move charge across time, descend
the deep sheet, exclude a scalar row, or prove LRC(14).

The highest-leverage next test is not another colour count.  It is whether
the carry/cut toothpick partition can be realized on the same Boolean
ancestry cospan as THM-2471/2478 without losing target charge or owner type.

## 8. Exact companion

Run

```text
python3 04-computation/lrc14_truncated_radon_escape_probe.py
python3 -O 04-computation/lrc14_truncated_radon_escape_probe.py
```

The dependency-free exact companion checks:

- rank `72` and kernel dimension six for the `13 x 7` row-zero map;
- all `924` six-slope banks and their common six-dimensional vertical
  kernel;
- a sharp integral cyclotomic-coefficient control with exactly five bad
  nonzero slopes;
- the explicit two-row THM-2506 defect, for which all twelve nonzero slopes
  survive;
- `C_13` all-colour saturation on every good output; and
- the representative-cut hostile (28)--(30).

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_truncated_radon_escape_probe.out
```

byte-for-byte.  An independent audit rederived the general rank formula,
the trace sharpness construction and its Fourier signs, Galois saturation,
the `18^(-11)` norm floor, the `1/6` and `1/4` fixed-slope invoices, and the
`q>=p+1` failure boundary.  It also verified that (28)--(30) expose a genuine
row-dependent cut and make no target, owner, or deep-current claim.  Normal
and optimized executions are byte-identical to the stored transcript, and
the source/output hashes match the metadata. QED.
