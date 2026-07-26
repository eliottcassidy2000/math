---
id: THM-2361
title: "Familywise fixed-colour cone and off-diagonal phase boundary"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. For any
  family of rational nonnegative root-fibre pairs supported in one open
  arc of length 1/7, one Galois automorphism depending only on the fixed
  primitive root colour places every cross-correlation in the same cone
  of width 4pi(ceil(N/7)-1)/N<pi, with an explicit positive real-part
  floor. The same holds on an arithmetic danger comb. Hence a full
  order-13 parameter family cannot be a nontrivial pure character.
  Off-diagonal colours instead acquire the exact absolute terminal
  phase exp(2pi i Delta x_s). At N=169 an exact two-interval hostile has
  diagonal and phase-compensated currents 1/2 but off-diagonal current
  zero. The compensation is a derived component observable not supplied
  by canon. No target landing, all-91-unit aggregate, scalar-row
  exclusion, profile decrement, or LRC(14) consequence follows.
source: codex-2026-07-25-familywise-fixed-colour-cone
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2334-relation-residue-current-and-character-twist-pushforward
related:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
script: 04-computation/familywise_fixed_colour_cone_thm2361.py
output: 05-knowledge/results/familywise_fixed_colour_cone_thm2361.out
script_sha256: b11d705b14cf9b997bb3035deb15f936092ad89acb2cc2d4a1d22af093ec8f07
output_sha256: 8df8d9f01900ce0dbcf9f4b749f24b139a2ef28ae76873da4fecba51704b2d0
hash_basis: working-tree bytes (LF)
---

# THM-2361 -- one root colour gives one familywise cone

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2323 proves noncancellation one fixed root colour at a time. The
Galois automorphism in that proof is independent of the target translate:
it straightens an entire family at once. This gives a genuine familywise
cone and rules out a tempting pure-character target law.

The boundary is exact. When the two root colours differ, the physical
correlation acquires a continuous terminal base phase. Discrete root
labels and diagonal cone confinement do not remove it.

## 1. A common cone for a fixed primitive colour

Let

```text
N>=2,

R=ceil(N/7)-1,

zeta=exp(2*pi*i/N).                                  (1)
```

Then

```text
R<N/7,

2*pi*R/N<2*pi/7<pi/2,

4*pi*R/N<4*pi/7<pi.                                 (2)
```

Let `t` range over any family. For each `t`, let `f_t,g_t` be
rational-valued step functions with rational breakpoints such that

```text
0<=f_t<=g_t,

support(f_t),support(g_t) subset I,

integral f_t>0,                                     (3)
```

where `I` is one common open circle arc of length `1/7`. Fix

```text
k in (Z/NZ)^*.
```

For `x_r(y)=(y+r)/N`, define

```text
M_(k,f_t)(y)
 =sum_(r mod N) f_t(x_r(y))zeta^(-kr),              (4)

C_k(t)
 =integral_0^1
   M_(k,f_t)(y)conjugate(M_(k,g_t)(y))dy.            (5)
```

As in THM-2323, expansion by root displacement gives

```text
C_k(t)
 =sum_(d=-R)^R c_(t,d)zeta^(kd),                    (6)

c_(t,d) in Q_(>=0),

c_(t,0)
 =N integral f_t g_t
 >=N integral f_t^2
 >0.                                                (7)
```

Use the single Galois automorphism

```text
sigma_(k^(-1)): zeta -> zeta^(k^(-1)).              (8)
```

It works simultaneously for every `t`:

```text
sigma_(k^(-1))(C_k(t))
 =sum_(d=-R)^R c_(t,d)zeta^d.                       (9)
```

Every summand in (9) lies in the closed sector

```text
|arg z|<=2*pi*R/N.
```

The positive diagonal makes the sum nonzero, and

```text
Re sigma_(k^(-1))(C_k(t))
 >=cos(2*pi*R/N)sum_d c_(t,d)
 >=cos(2*pi*R/N)c_(t,0)
 >0.                                               (10)
```

Thus the entire family lies in one cone of width

```text
4*pi*R/N<pi.                                       (11)
```

This is stronger than memberwise nonvanishing. In the common
Galois-conjugate embedding, sums with arbitrary nonnegative real weights
cannot cancel. Equivalently, in the original embedding, rational
nonnegative combinations of the `C_k(t)` cannot cancel, because the one
automorphism fixes those weights. No claim is made that it fixes arbitrary
irrational real weights.

## 2. Arithmetic-comb form

For a positive integer `a`, put

```text
D_a={x:||a x||<1/14}.                               (12)
```

Assume

```text
gcd(a,N)=1
```

and replace the common-arc condition in (3) by

```text
support(f_t),support(g_t) subset D_a.               (13)
```

If `c_(t,d)!=0`, the least signed residue `e_d` of `ad modulo N`
satisfies

```text
|e_d|<N/7.
```

The one automorphism

```text
zeta -> zeta^(a k^(-1))                             (14)
```

sends every surviving phase `zeta^(kd)` to `zeta^(e_d)`. Equations
(9)--(11) follow verbatim for the whole family. The automorphism depends
on the physical comb and frozen colour, not on `t`.

## 3. Pure order-thirteen character obstruction

Assume now that `13|N`, that the parameter set contains a complete
order-thirteen character fibre, and that

```text
C_k(t)=c chi(t),                 c!=0,               (15)
```

where `chi` has image all thirteenth roots of unity. Apply the same
automorphism (8), or (14) in the comb case. It permutes the thirteenth
roots, so the transformed values in (15) are a nonzero scaled and rotated
full set of thirteenth roots.

Their sum is zero. But (10)--(11) place every one of them in the same
open half-plane, where a sum of nonzero vectors cannot vanish. This
contradiction proves:

```text
one fixed-colour family cone
  => no nontrivial pure order-13 character law.      (16)
```

This does not prove that a nonzero target coefficient survives; a more
general family may remain inside the cone without being a character.

## 4. Exact off-diagonal identity

The cone argument is diagonal in root colour. Define the periodic gauges

```text
N_(k,f)(y)
 =exp(-2*pi*i k y/N)M_(k,f)(y).                     (17)
```

Let

```text
l=k+D
```

and, for a gauge offset `j in Z`, put

```text
C_(k,l)^(j)
 =integral_0^1
   exp(2*pi*i j y)
   N_(k,f)(y)conjugate(N_(l,g)(y))dy.               (18)
```

The gauge Fourier expansions give

```text
C_(k,l)^(j)
 =N^2 sum_h
   f_hat(k+Nh)
   conjugate(g_hat(l+N(h+j))).                      (19)
```

Put

```text
Delta=D+Nj.
```

Returning to root branches `x_r=(y+r)/N` gives the exact physical form

```text
C_(k,l)^(j)
 =sum_(r,s mod N) integral_0^1
   f(x_r)g(x_s)
   exp(2*pi*i Delta x_s)
   zeta^(k(s-r))dy.                                 (20)
```

The factor

```text
exp(2*pi*i Delta x_s)                               (21)
```

is the absolute terminal base phase. It depends continuously on the
terminal component position and is not controlled by the short
root-displacement support in Section 1.

If one multiplies each component contribution in (20) by the inverse of
(21), the result is exactly the diagonal fixed-colour correlation

```text
C_k(f,g).                                          (22)
```

This phase-compensated quantity is a valid derived component observable.
Current canon does not supply the componentwise inverse phase as a lawful
physical LRC probe.

## 5. Exact N=169 hostile

Take

```text
N=169,

k=1,

l=14,

D=13,

j=0,

epsilon=1/1352,                                    (23)
```

and let

```text
f=g
 =1_(
   (-1/52-epsilon,-1/52+epsilon)
   union
   (1/52-epsilon,1/52+epsilon)
  ).                                                (24)
```

Both colours are primitive modulo `169`, and the support lies strictly
inside `(-1/14,1/14)`. Multiplication by `169` sends the two interval
centres to

```text
3/4, 1/4 mod 1
```

and their half-width to `1/8`. The images are disjoint one-sheet fibres.
Therefore the diagonal is

```text
C_1=169 mu(f)=1/2.                                  (25)
```

For the off-diagonal current, (20) reduces to

```text
C_(1,14)^(0)
 =169 integral f(x)exp(2*pi*i 13x)dx.               (26)
```

The equal interval integrals have centre phases

```text
exp(-pi*i/2)=-i,

exp(pi*i/2)=i.
```

They cancel exactly:

```text
C_(1,14)^(0)=0.                                    (27)
```

After componentwise multiplication by the inverse terminal phase, the
current is again

```text
1/2.                                                (28)
```

Thus positive diagonal currents, primitive root colours, equal component
masses, one-sheet root fibres, and a common short support do not control
the off-diagonal phase. The failure is entirely the missing factor (21).

## 6. LRC application and stopping boundary

Here is the lawful reduction of the THM-2334 target translates to the
hypotheses above. Work modulo `13` with the THM-2309 owner-pivot packet.
If `j` is its owner blocker, the pivot row

```text
r_j congruent -w_(u_0)e_j mod 13
```

belongs to `L`. Hence every quotient-dual target twist

```text
ell in G^=L^perp/<w>
```

has

```text
ell_j=0.                                            (29)
```

This conclusion is representative-independent: `w_j=0 mod 13`, so
adding a multiple of `w` does not change (29). Thus all `169` lawful
twists leave the owner danger factor unshifted. The delayed word is also
unchanged because its clock `R` is divisible by `13`.

Let `E_ell` be the present owner packet with the remaining coordinates
translated by `ell/13`, and put

```text
E_(Q,ell)=E_ell intersection W.
```

Perron-normalize by the owner speed:

```text
g_ell=P_(|c_j|)1_(E_ell),

f_ell=P_(|c_j|)1_(E_(Q,ell)).                       (30)
```

Positivity of the Perron operator and the unshifted owner factor give

```text
0<=f_ell<=g_ell,

support(f_ell),support(g_ell) subset D_1.           (31)
```

These are rational step functions with rational breakpoints. Expose
`N=13`, fix any primitive colour `k`, and form the auxiliary diagonal
root correlation `C_k(ell)` from (5). If some `f_ell` is zero, then
`C_k(ell)=0`, already contradicting a law

```text
C_k(ell)=c chi(ell),             c!=0.              (32)
```

Under a hypothetical law (32), every `f_ell` therefore has positive
integral. For a nontrivial quotient character `chi`, restrict to an
affine target line on which it assumes all thirteen roots. Sections
1--3 rule out (32).

This is an auxiliary diagonal family made from each lawful target
translate. It is not THM-2334's off-diagonal marked target current
`H(ell)`, and the pure-character obstruction is not a target landing.

The marked edges used by THM-2327, THM-2349, and THM-2354 are
off-diagonal: their endpoint frequencies differ by a deepest-comb
multiple. Formula (20) therefore inserts the terminal phase (21).
Section 5 proves that no diagonal cone argument can simply discard it.

The theorem preserves

```text
fixed primitive root colour,
common shallow danger support,
whole finite target-translate family,
rational nonnegative correlation coefficients,
one familywise acute cone.
```

It does not preserve or select

```text
off-diagonal terminal phase,
relation target,
all-91-unit marked edge,
component phase-tree transport,
graph-channel pair ratio.
```

No scalar row is excluded, the ledger remains `165`, and LRC(14) remains
open. The exact missing service is a lawful component phase ratio capable
of compensating (21), or a different target observable which stays
diagonal in the frozen root colour.

## 7. Exact companion

The dependency-free companion uses integer and `Fraction` arithmetic. It:

- checks one common Galois exponent on `7,921,099` family/support cases
  for every modulus `2<=N<=200`;
- checks `10,473,536` arithmetic-comb straightening cases;
- verifies all `12` nontrivial order-thirteen character permutations;
- checks `42,000` off-diagonal gauge/frequency ledgers; and
- verifies the `N=169` hostile geometry, diagonal `1/2`, exact `Q(i)`
  cancellation, and compensated `1/2`.

Run

```bash
python3 04-computation/familywise_fixed_colour_cone_thm2361.py
python3 -O 04-computation/familywise_fixed_colour_cone_thm2361.py
```

Both transcripts must match

```text
05-knowledge/results/familywise_fixed_colour_cone_thm2361.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
