---
id: THM-2529
title: "Deep-comb adjacent-fibre odd consumer and the exact zero-target boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The thirteen shifted
  copies of the actual LRC danger comb contain exactly one or two roots on
  every generic fibre, and a two-root fibre consists of adjacent residues.
  Consequently the THM-2527 oriented cut coordinate collapses, in the
  natural deep-colour chart, to the pointwise identity
  Psi=2+5 d(13cx).  For every nonzero nonnegative same-integrand weight F,
  the resulting odd coordinate is exactly
  O(-4)=2 integral F+5 integral F d(13cx)>0.  The entire weighted deep
  correlation profile has an explicit two-mass normal form, all twelve
  nontrivial root Fourier modes are strictly positive, and its odd A_1 H_1
  image has all twelve modes nonzero.  This still need not create target
  charge: the comb-compatible selector F_t=Delta_(t+1)d(13cx) realizes the
  circulant diagonal law H(r,s,t)=1/91 when r=t+1 and zero otherwise, while
  its nonzero all-mode odd root bank is independent of (s,t).  This selector
  is an exact comb-compatible hostile, not a proved canonical owner packet.
  No target landing, semantic ordered edge, scalar-row exclusion, or LRC(14)
  proof is claimed.
source: codex-2026-07-27-deep-comb-odd-consumer
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2526-affine-skew-orientation-gauge-boundary
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
related:
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2524-translated-chi7-hamilton-polarization-inversion
  - THM-2525-unit-guard-collision-floor-and-word-owner-cross-cospan-collapse
script: 04-computation/lrc14_deep_comb_odd_consumer_thm2529.py
output: 05-knowledge/results/lrc14_deep_comb_odd_consumer_thm2529.out
script_sha256: ec8912b129fc85f528be735d27ca5d1ba12be10924ac5df307e5e52e55319c9a
output_sha256: 7d26ec6df530d52cc4cc72f6eea8af564a13e9fe0da8719b4375f11adeb29303
hash_basis: working-tree bytes (LF)
---

# THM-2529 -- the actual deep comb has only one adjacent edge

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2527 finds a positive Boolean coordinate of the oriented cyclic odd
collision transform for an arbitrary thirteen-bit root mask.  The actual
deep danger comb occupies a much smaller part of that cut polytope.  Every
one of its masks is either a singleton or one adjacent pair.  On these two
orbits the 98-valued general cut score takes only the values seven and two:

```text
singleton deep fibre:       Psi=7;
adjacent-pair deep fibre:   Psi=2.                          (1)
```

This gives a fixed positive odd coordinate for every nonnegative lawful
same-integrand weight, with no mixing threshold and no choice after
integration.  It also exposes a sharp semantic boundary.  Autocorrelation
forgets the location of the singleton or adjacent pair.  Therefore all root
modes and a positive odd coordinate can coexist with complete target
constancy.

## 1. The one-or-two adjacent-root law

Put `p=13` and use the centred danger interval

```text
d(y)=1_(||y||<1/14).                                       (2)
```

For a positive integer deep coefficient `c`, define the thirteen shifted
deep-comb bits

```text
Delta_r(x)=d(c x-r/13),                  r in F_13.          (3)
```

Away from their finite endpoint set, the exact successor identity is

```text
n(x):=sum_r Delta_r(x)=2-d(13c x).                         (4)
```

Indeed, the open arc in (2) has length `1/7`, so a translate of the
thirteen-grid meets it in one or two points.  The usual floor/ceiling count,
or the thirteen-shift comb identity, gives (4).  Since

```text
1/7<2/13,                                                   (5)
```

two occupied grid points must be consecutive in the cyclic residue order.
Thus the only generic root masks are

```text
{r},                         when d(13c x)=1;
{r,r+1},                     when d(13c x)=0,               (6)
```

up to reversing the printed orientation at a null boundary.  The natural
deep-colour chart in (3) orients increasing residue by the positive shift
`r -> r+1`.  A global reversal changes every odd-bank sign but none of the
nonvanishing or positive/converse-pair conclusions below.

## 2. Exact two-mass normal form

Let `F` be any nonnegative integrable function and define the local deep
autocorrelation

```text
c_u(x)=sum_r Delta_r(x)Delta_(r+u)(x).                      (7)
```

Use the branch normalization inherited from THM-2527:

```text
b_u=1/13 integral_T F(x)c_u(x) dx.                          (8)
```

Separate the weight carried by singleton and adjacent-pair fibres:

```text
alpha=integral_T F(x)d(13c x) dx,

beta =integral_T F(x)(1-d(13c x)) dx.                       (9)
```

Equations (6)--(7) give the whole profile, not merely one moment:

```text
b_0=(alpha+2 beta)/13,

b_1=b_(-1)=beta/13,

b_u=0,                         u notin {0,+1,-1}.            (10)
```

The mean of `b` may be removed without changing any operator below.  For
the normalized root transform

```text
b_hat(k)=1/13 sum_u b_u zeta_13^(-ku),                      (11)
```

equation (10) gives, for `k!=0`,

```text
b_hat(k)
 =1/169 [alpha+beta(2+zeta_13^k+zeta_13^(-k))]

 =1/169 [alpha+beta|1+zeta_13^k|^2].                        (12)
```

Therefore

```text
F nonzero and F>=0
  -> alpha+beta>0
  -> b_hat(k)>0                  for every k!=0.             (13)
```

This all-mode conclusion is geometric and does not require rationality.
When `F` is a rational step function, `b` is rational as well, so (13) also
agrees with the `Phi_13` irreducibility mechanism used elsewhere in the
canon.

## 3. The odd cut coordinate becomes `2+5d(13cx)`

Let `A_1` be the symmetric `chi_7` Hamilton operator of THM-2523/2524 and
let `H_1` be the skew cyclic Hilbert/tournament operator of THM-2526.  In
the natural deep-residue chart define

```text
O=13 A_1 H_1 b.                                            (14)
```

THM-2527's local fixed-coordinate score is

```text
Psi(e)=(A_1H_1 c(e))_(-4)

 =7c_0-12c_1+8c_2-6c_3+7c_4-6c_5+2c_6.                   (15)
```

For a singleton, only `c_0=1`, so `Psi=7`.  For one adjacent pair,

```text
c_0=2,                    c_1=c_(-1)=1,                    (16)
```

and every other entry vanishes, so `Psi=14-12=2`.  Combining this with
(4) proves the pointwise identity

```text
Psi(Delta(x))=2+5d(13c x).                                 (17)
```

The factor `1/13` in (8) cancels the leading `13` in (14).  Hence the same
integrand gives the exact positive consumer

```text
O(-4)
 =integral_T F(x)Psi(Delta(x)) dx

 =2 integral_T F(x)dx
    +5 integral_T F(x)d(13c x)dx

 =7alpha+2beta.                                            (18)
```

In particular,

```text
F>=0,  F nonzero
  -> O(-4)>=2 integral F>0,

O(4)=-O(-4)<0.                                             (19)
```

If `F` is Boolean, (18) is already a finite positive Boolean lift:

```text
O(-4)=2 measure(F)+5 measure(F intersection D_(13c)).       (20)
```

It uses only two Boolean atoms with integer multiplicities, rather than the
98 threshold layers needed for a general root mask.

Both `A_1` and `H_1` are invertible on the centred twelve-dimensional root
module.  Combining their nonzero Fourier multipliers with (13) gives

```text
O_hat(k)!=0                         for every k!=0.          (21)
```

Thus one fixed positive coordinate, its negative converse, and all twelve
primitive root modes occur simultaneously.

## 4. Application to the lawful target/deep tensor

Retain THM-2365's lawful target co-shifts.  Its nonnegative word-bearing
owner weight is

```text
F_(s,t)=E_(s,t)W>=0,                                      (22)
```

and the deep-coordinate complement inside `E_(s,t)` gives

```text
F_(s,t)(x)Delta_t(x)=0                                    (23)
```

almost everywhere.  For every fixed `(s,t)` with positive mass, equations
(8)--(21), with `F=F_(s,t)`, produce on that same integrand

```text
O_(s,t)(-4)
 =2 measure(F_(s,t))
   +5 measure(F_(s,t) intersection D_(13c))
 >0,                                                       (24)
```

together with all twelve root modes.  This retains the owner, delayed word,
lawful co-shift labels, and the fact that the root `t` is forbidden before
the final autocorrelation is taken.

Equation (24) is stronger than a post-integration sign choice.  But its
consumer is translation-invariant in the root label: the autocorrelation in
(7) knows whether the comb mask is a singleton or an adjacent pair and
forgets where that mask sits.  In particular, (23) is not visible in the
normal form (10).  This is the exact reason a positive root bank is not yet
a target landing.

## 5. Exact compatibility with the THM-2365 zero-target line

The loss is attained by a comb-compatible rational Boolean hostile.  For
each deep target label `t`, put

```text
F_t(x)=Delta_(t+1)(x)d(13c x).                            (25)
```

The factor `d(13cx)` restricts to singleton fibres.  Hence `F_t` selects
exactly those phases whose unique active deep root is `t+1`.  Haar
invariance and the thirteen equal root cells give

```text
measure(F_t)=1/91,

F_t Delta_t=0.                                             (26)
```

Make this control independent of the first target co-shift `s`.  Its
THM-2365 incidence tensor is then

```text
H(r,s,t)
 :=integral_T F_t(x)Delta_r(x)dx

 =1/91,                       r=t+1;

 =0,                          otherwise.                   (27)
```

Thus

```text
H(t,s,t)=0,

H(r,s,t)=G(r-t),

G(v)=1/91 1_(v=1).                                        (28)
```

Equation (28) is exactly THM-2365's circulant/zero-`H`-drift law.  It is a
nonnegative rational scalar-comb realization of that law and respects the
deep diagonal exclusion.  It is **not** claimed to equal
`E_(s,t)W` for one canonical nine-coordinate owner packet; its role is an
exact comb-compatible consequence-level hostile.

On (25), `alpha=1/91` and `beta=0`.  Therefore

```text
b=(1/1183)delta_0,

O(-4)=7/91=1/13,                                          (29)
```

and the root bank has all twelve primitive modes.  Yet the autocorrelation
of a singleton is `delta_0` regardless of the singleton's location.
Consequently the entire root bank in (29) is independent of `(s,t)`:

```text
root-mode support:                 all twelve nonzero modes;

target-shift variation:            zero;

nonzero target Fourier support:    none.                    (30)
```

This proves the sharp failed implication

```text
positive same-integrand odd root coordinate
  + all primitive root modes
  + retained forbidden deep label before contraction

does not imply

nonzero lawful target charge after autocorrelation.         (31)
```

The first missing map is now explicit: one must anchor the location of an
active deep root relative to `t`, or use an ordered cross-cospan before the
root-location quotient.  Another translation-invariant self-correlation
cannot recover it.

## 6. Relation to target-coshift cross-cospans

There is a parallel elementary diagnostic for the word/bare-owner repair
suggested by THM-2525.  Let `E_h` be a lawful nonzero target co-shift, let
`F_0=E_0W`, and use the root-invariance of the delayed word.  The relative
cross profile

```text
C_h(u)=integral F_0(x)E_h(x+u/13)dx                         (32)
```

has odd part

```text
C_h(u)-C_h(-u)
 =integral W(x)[
    E_0(x)E_h(x+u/13)-E_h(x)E_0(x+u/13)
   ]dx.                                                     (33)
```

Thus the first lawful ordered observable is a Boolean `2 x 2` plaquette
determinant.  Retaining `h` makes (33) legal, but `h!=0` alone does not make
the determinant nonzero.  Equation (31) is the deep-comb version of the
same boundary: a label retained before a symmetric contraction can still be
erased by that contraction.  An anchored Gram or cross-cospan must prove a
nonzero determinant, not merely exhibit its two labelled legs.

## 7. Exact companion and stopping boundary

Run

```bash
python3 04-computation/lrc14_deep_comb_odd_consumer_thm2529.py
python3 -O 04-computation/lrc14_deep_comb_odd_consumer_thm2529.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_deep_comb_odd_consumer_thm2529.out
```

byte-for-byte.  The dependency-free exact referee verifies:

- all `26` generic deep phase cells, split as thirteen singleton and
  thirteen adjacent-pair cells;
- the successor identity (4), adjacency, and the local values `7,2`;
- the two-mass profile (10), fixed consumer (18), oddness, and `96`
  primitive-mode controls;
- the complete `13^3` tensor (27), diagonal exclusion, and circulant law;
  and
- the coexistence of twelve nonzero root modes with zero target variation
  in the hostile.

The theorem gives a sharp positive Boolean consumer on every positive
lawful target/deep integrand.  It does not prove variation in `(s,t)`, retain
the active root's location through autocorrelation, identify one semantic
source--arrival edge, exclude any scalar row, or prove LRC(14).

## 8. Independent audit

An independent audit rederived the one-or-two adjacent-root law, the branch
normalization in (8)--(10), the normalized `1/169` Fourier factor and strict
`|1+zeta_13^k|^2` positivity, the leading-`13` cancellation in (18), and the
identity `7alpha+2beta=2 integral F+5alpha`.  It separately checked the
`1/91` singleton selector, the circulant tensor and zero-target scope in
(25)--(30), and the delayed-word change of variables in the cross-plaquette
identity (33).  Normal and optimized executions reproduce the stored output
byte-for-byte, and both recorded hashes match the working-tree artifacts.
**QED.**
