---
id: THM-2353
title: "Target plaquette holonomy and the cross-axis endpoint-slice scout"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For a nonzero 13 by 13 phase-free target response, the inverse
  deep-character boundary is equivalent to zero-free adjacent-plaquette
  flatness followed by two one-axis laws. The physical plaquette
  holonomy factors as the left endpoint holonomy times the conjugate
  right endpoint holonomy, and deep-character modulation preserves
  every plaquette minor up to a nonzero phase. Multiplicative
  plaquette curvature and THM-2340's additive ANOVA interaction are
  independent diagnostics. On a typed nine-factor c1-owner control,
  an exact finite beta=0 endpoint slice has all 169 target twists and
  all 169 cyclic adjacent plaquette minors nonzero. This is a retained
  finite slice only: the omitted infinite endpoint tail may cancel it,
  the control is not asserted to be a scalar row, and no canonical row
  exclusion or LRC(14) follows.
source: codex-2026-07-25-target-plaquette-holonomy
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2340-owner-word-anova-target-landing
  - THM-2343-deep-comb-affine-target-catalyst
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
related:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
script: 04-computation/lrc14_target_plaquette_holonomy_thm2353.py
output: 05-knowledge/results/lrc14_target_plaquette_holonomy_thm2353.out
script_sha256: e4dbdd977c2f2cb2bbc782465633cb233018e6c5055a4c88fc03456f4d371461
output_sha256: cbafa04c5ecd677de2aaa7875f6fb7b8812ef6667f0fd765d5ef06d42b41e2bf
hash_basis: working-tree bytes (LF)
---

# THM-2353 -- target plaquette holonomy and cross-axis mixing

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2343 reduces zero-only target landing to one inverse-character line,
and THM-2344 shows that arbitrary word content transported along one
aligned active coordinate can lie on that line. The first missing
coordinate is therefore genuinely two-dimensional. The smallest local
two-dimensional observable is a `2 x 2` determinant:

```text
Delta_A(s,t)
 =A(s,t)A(s+1,t+1)-A(s+1,t)A(s,t+1).              (1)
```

It measures multiplicative curvature. A bad inverse character has every
such determinant zero. More importantly, when no entry vanishes,
flatness reconstructs the whole matrix from its two axes, after which
the bad line is decided by two one-dimensional laws.

The canonical endpoint product gives this curvature a physical meaning:
its normalized holonomy is a mismatch between the two endpoint
holonomies. An exact finite slice of the first typed nine-factor control
exhibits the desired mismatch at every plaquette. That computation is a
scout, not a proof about the full endpoint sum: all omitted endpoint
harmonics remain capable of cancelling the retained slice.

## 1. The hierarchical inverse-character test

Write

```text
F=F_13,                 zeta=exp(2*pi*i/13),
K:F x F -> C,           m in F minus {0}.          (2)
```

In the LRC application, THM-2343 supplies

```text
K(0,0)!=0
```

and identifies the zero-only boundary as

```text
K(s,t)=c zeta^(-m t),       c!=0.                  (3)
```

Assume henceforth that `K(0,0)!=0`. All indices below are cyclic in `F`.
There are two useful levels of flatness.

First, with no nonvanishing hypothesis, every `2 x 2` minor of `K`
vanishes if and only if `K` has matrix rank one. Because `K(0,0)!=0`,
this is equivalently the reconstruction identity

```text
K(s,t)K(0,0)=K(s,0)K(0,t)       for every s,t.     (4)
```

Thus a nonzero arbitrary minor already excludes (3).

Second, suppose every entry of `K` is nonzero. Then it is enough to
check the `169` cyclic adjacent minors (1):

```text
Delta_K(s,t)=0 for every s,t

iff K(s,t)=K(s,0)K(0,t)/K(0,0) for every s,t.      (5)
```

Indeed, `Delta_K(s,t)=0` says that

```text
K(s+1,t)/K(s,t)
```

is independent of `t`. Multiplying these row ratios from row zero to
row `s` proves (5). The reverse implication is immediate.

Combining (3) and (5) gives the exact hierarchical criterion:

```text
K(s,t)=c zeta^(-m t), c!=0

iff
  K is nowhere zero,
  Delta_K(s,t)=0 for all s,t,
  K(s,0)=K(0,0) for all s,
  K(0,t+1)=zeta^(-m)K(0,t) for all t.              (6)
```

Consequently escape is certified at the first successful level:

```text
one zero entry;
or one nonzero plaquette minor;
or, after flatness, failure of either axis law.     (7)
```

This is stronger than using a global rank statistic. It names the first
local witness and, on the flat locus, leaves only two one-dimensional
profiles.

The axis stage is indispensable. Embed any nowhere-zero one-dimensional
array `b_t` as `K(s,t)=b_t`; every plaquette is flat, whether or not
`b_(t+1)/b_t=zeta^(-m)`. In particular, THM-2355's audited
perfect-autocorrelation hostile can be placed on one axis and survives
every curvature test while failing the last law in (6). Curvature is a
cheap escape certificate, not an equivalence by itself.

## 2. Physical endpoint holonomy

Retain THM-2334/2343's phase-free endpoint factorization

```text
K(s,t)=d L(s,t)conjugate(F(s,t)),                  (8)
```

where `d=delta_hat(m)!=0` is independent of `(s,t)`. At one plaquette,
abbreviate `(s,t),(s+1,t),(s,t+1),(s+1,t+1)` by
`00,10,01,11`. Direct expansion gives the exact four-twist identity

```text
Delta_K
 =d^2[
    L_00 L_11 conjugate(F_00 F_11)
   -L_10 L_01 conjugate(F_10 F_01)
 ].                                                (9)
```

Whenever the four values of an array `A` are nonzero, define its
normalized plaquette holonomy

```text
Omega_A
 =A_00 A_11/(A_10 A_01).                          (10)
```

Equation (8) then factors (10) as

```text
Omega_K=Omega_L conjugate(Omega_F).                (11)
```

Therefore the exact endpoint mismatch test is

```text
Omega_L conjugate(Omega_F)!=1.                    (12)
```

It is not generally the test `Omega_L!=Omega_F`: those are equivalent
only under an additional unit-modulus condition. This conjugate-
reciprocal distinction is the information lost by comparing the two
endpoints as unnormalised pictures.

The full target response is the deep modulation

```text
H(s,t)=zeta^(m t)K(s,t).
```

Both terms of a plaquette minor receive the same phase, so

```text
Delta_H(s,t)
 =zeta^(m(2t+1))Delta_K(s,t).                     (13)
```

Thus plaquette nonvanishing can be tested before or after restoring the
deep character. In particular a curved phase-free response cannot be
the inverse character in (3), and the full response cannot be constant.

## 3. Plaquette curvature and additive ANOVA are transverse

THM-2340's fork diagnostic is additive double-centring. Plaquette
holonomy is multiplicative. Neither implies the other.

For an arbitrary rank-one matrix

```text
A(s,t)=a_s b_t,
alpha=1/13 sum_s a_s,
beta =1/13 sum_t b_t,                             (14)
```

every plaquette minor vanishes, while its ANOVA interaction is

```text
I_A(s,t)=(a_s-alpha)(b_t-beta).                    (15)
```

The three nonconstant energy pieces are therefore

```text
row main effect   =|beta|^2 Var(a),
column main effect=|alpha|^2 Var(b),
fork interaction =Var(a)Var(b),                   (16)
```

with normalized variances. A flat rank-one matrix has nonzero fork
interaction whenever both one-axis factors are nonconstant.

Conversely, an additively separable matrix

```text
A(s,t)=f_s+g_t
```

has zero ANOVA interaction, but

```text
Delta_A(s,t)
 =(f_s-f_(s+1))(g_(t+1)-g_t).                     (17)
```

It can be curved at every plaquette. The exact companion uses rational
cyclic profiles with all `169` minors nonzero for (17), and a nonconstant
rank-one rational matrix with all `169` minors zero but nonzero
interaction.

Hence curvature certifies escape from the single bad line and therefore
some nonzero full target coefficient. It does not by itself say whether
that coefficient lies on THM-2340's pure-other, pure-deep, or fork
support locus. The additive ANOVA sidecar remains necessary for
word-matching.

## 4. A finite cross-axis slice of the typed nine-factor control

Use the `c_1`-owner typed control in labelled order

```text
w=(13,13^3,2*13^5,1,14,27,40,53,66),

X=13,
Y=X+w_2=742599,                                   (18)
```

with target axes the `c_2,c_3` coordinates `1,2` and speed-one guard
coordinate `3`; the deepest multiplier is `m=1`. On the typed control
used by THM-2331, choose the concrete THM-2309 owner packet which omits
`H`, grafts `c_2` to `q_1` at coordinate `4`, and grafts `c_3` to `q_2`
at coordinate `5`. Modulo thirteen, its two graft rows are

```text
e_H-e_(q_1)-e_(c_2),          e_H-e_(q_2)-e_(c_3). (19)
```

The packet has rank six. In the gauge `ell_H=0`, a quotient-dual basis is

```text
ell(s,t)
 =s(e_(c_2)-e_(q_1))+t(e_(c_3)-e_(q_2)).           (20)
```

Both displayed vectors annihilate every packet row, and together with
the speed gauge they have rank three. The speed vector annihilates the
packet because every packet row is an exact relation. Since the packet
rank is six, these three vectors are its full annihilator, and quotienting
by the speed gauge leaves the required two-dimensional dual. Thus (20),
rather than the raw ambient vectors `e_(c_2),e_(c_3)`, parametrizes the
canonical `169` target characters. The speed vector is the same strict
valuation-type control used in THM-2331/2344. It is **not** asserted to
satisfy the scalar-cover equations of a canonical row.

For `Z` equal to `X` or `Y`, retain the finite endpoint bank

```text
S_Z={
 u in Z^9:
 u_i in {-1,+1} for i!=3,
 u_3=Z-sum_(i!=3)w_i u_i,
 7 does not divide u_3
}.                                                 (21)
```

There is one candidate for each of the `2^8` non-guard sign vectors.
The guard filter in (21) removes every nonzero-mode Fourier zero of the
centred guard complement and deliberately leaves its zero mode outside
this chosen slice. Literal enumeration and an independent mod-seven
dynamic program give

```text
|S_X|=218,                   |S_Y|=220.            (22)
```

For the left marked endpoint, retain the terminal-word mode

```text
beta=0.
```

Its coefficient is a positive common scalar. Every non-guard base
factor is real and even, so its coefficients at `-1` and `+1` agree.
The direct centred interval/complement formula makes each such
coefficient nonzero. After removing all common nonzero endpoint scalars,
the only variable weight is

```text
q(h)=sin(2*pi*h/7)/(h sin(2*pi/7)).                (23)
```

Put `c=2cos(2*pi/7)`. For `h` nonzero modulo seven,

```text
h mod 7:       1,   2,       3,          4,     5,  6

q(h)*h:        1,   c,   c^2-1,   -(c^2-1),    -c, -1. (24)

c^3+c^2-2c-1=0.
```

Group each endpoint bank by its two target residues:

```text
U_Z(r_1,r_2)
 =sum_(
     u in S_Z,
     (u_1-u_4,u_2-u_5)=(r_1,r_2) mod 13
   ) q(u_3).                                       (25)
```

These are exactly the phases supplied by (20). All nine cells
`(r_1,r_2) in {0,2,-2}^2` occur. Their term counts, with rows and columns
ordered as `0,2,-2`, are

```text
X endpoint:
  [54 28 27]
  [26 14 14]
  [28 14 13]

Y endpoint:
  [56 27 27]
  [28 14 13]
  [28 13 14].                                      (26)
```

With `zeta=zeta_13`, define the two finite endpoint transforms and their
phase-free response by

```text
E_Z(s,t)
 =sum_(r in {0,2,-2}^2)U_Z(r)zeta^(s r_1+t r_2),

K_slice(s,t)=E_X(s,t)conjugate(E_Y(s,t)).          (27)
```

Up to one nonzero scalar, (27) is exactly the part of the physical
endpoint product obtained from (21) and `beta=0`, twisted by the lawful
chosen owner-packet target characters (20). The exact certificate proves

```text
K_slice(s,t)!=0                         for all 169 (s,t),

Delta_(K_slice)(s,t)!=0                 for all 169 cyclic
                                        adjacent plaquettes.         (28)
```

This is maximally curved at the local cyclic scale. It sharply differs
from THM-2344's arbitrary same-axis transported-word hostile, whose
response has all plaquettes flat.

## 5. Exact algebraic certificate

Every cell in (25), transform in (27), and minor in (28) lies in

```text
Q(c,zeta_13).
```

The two cyclotomic conductors are coprime, so the canonical
characteristic-zero basis

```text
{c^a zeta^b: 0<=a<=2, 0<=b<=11}                  (29)
```

has dimension `36`. The companion reduces its rational coefficients
modulo the prime

```text
P=1,000,003
```

in the quotient algebra

```text
F_P[c,zeta]/
 (c^3+c^2-2c-1, 1+zeta+...+zeta^12).              (30)
```

Every retained integer denominator is nonzero modulo `P`. If a
characteristic-zero element were zero, its canonical coefficient vector
would reduce to zero in (30). Thus a nonzero reduced vector is a valid
one-way certificate of characteristic-zero nonvanishing; no claim that
(30) is a field is needed.

The script obtains a nonzero reduced vector for every one of the `169`
entries and `169` minors. Before grouping, it independently evaluates
both endpoint transforms directly over all retained terms, for `338`
agreement checks. It also freezes the row-major modular coefficient
arrays with

```text
response digest:
b67c4454b273a2279bd092acebf8fc5a11cc7e1d8f886cddbe15f8f1f71a3b90

plaquette digest:
67b52827f88ab202dd4664e7ea0b98546d3e44092eaad6f8000b2fc485a08221. (31)
```

All load-bearing assertions use explicit exceptions, so optimized
Python retains them.

## 6. Scope and next decisive test

Equation (28) is a finite-slice statement, not a noncancellation theorem
for the canonical current. The full endpoints contain infinitely many
additional base harmonics and, on the left, nonzero terminal-word
harmonics. A determinant is nonlinear, so the full plaquette minor is
not the sum of the retained slice minors. The omitted tail can cancel
both entries and curvature:

```text
finite slice nonzero:             VERIFIED-EXACT;
full endpoint curvature:          OPEN;
word-matching ANOVA component:     OPEN;
all-91-unit visible component:     OPEN;
terminal phase retention:          OPEN;
canonical scalar row excluded:     NO;
LRC(14):                           OPEN.            (32)
```

The result nevertheless identifies a concrete next object. One should
retain four neighbouring THM-2334 twists together and control the
endpoint holonomy mismatch (12), rather than attempt to prove
nonconstancy at one twist at a time. The cheapest analytic upgrade would
bound the omitted endpoint tail strongly enough to preserve one
nonzero minor. The cheapest algebraic upgrade would isolate a
word-matching component of one curved four-twist packet and then attach
the all-`91`-unit and terminal-phase sidecars.

There is also an energy-side version of the same debt. THM-2355 shows
that complex pair-twist transports plus singleton magnitudes reconstruct
relative phases, whereas individual target-twist energies retain only an
autocorrelation. The THM-2356 candidate proposes a planar target--jet
chirp as a global replacement. Either service would preserve phase
information which the complex holonomy (11) needs; independent linear
twist magnitudes alone do not.

## 7. Exact companion

Reproduce with

```bash
python3 04-computation/lrc14_target_plaquette_holonomy_thm2353.py
python3 -O 04-computation/lrc14_target_plaquette_holonomy_thm2353.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_target_plaquette_holonomy_thm2353.out
```

byte-for-byte after LF normalization. Besides the finite-slice
certificate, the companion checks the rank-one reconstruction control,
the two ANOVA/curvature separation controls, inverse-character
flatness, deep-modulation covariance, the physical determinant identity,
and the conjugated endpoint-holonomy factorization.

## 8. Independent audit

The independent audit rederived the zero/curvature/axis hierarchy, the
endpoint determinant and conjugated-holonomy identities, deep-modulation
covariance, and the rank-one ANOVA formulas. It independently replayed
the characteristic-zero slice calculation as well as the ordinary and
optimized modular companion.

The audit found one material typing defect in the first candidate: raw
ambient `c_2,c_3` coordinate phases do not annihilate the grafted owner
packet. The repaired theorem now uses the quotient-dual differences in
(20). Independent mod-thirteen elimination verified packet rank six,
annihilator

```text
span(w,e_(c_2)-e_(q_1),e_(c_3)-e_(q_2)),
```

and the two-dimensional quotient basis after removing the speed gauge.
The repaired nine-cell computation again gave all `169` entries and all
`169` cyclic minors nonzero. The auditor also checked the guard zero-mode
exclusion, denominator reduction, non-field quotient logic, and every
omitted-tail and no-LRC-closure disclaimer.
