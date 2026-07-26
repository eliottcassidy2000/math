---
id: THM-2364
title: "Anchored corner forces a fully mixed deep/blocker-probe colour"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Let a
  positive set F avoid one deepest danger comb at shift zero and lie
  inside t=1 or 2 named blocker combs at shift zero. The fully mixed
  corner sum of its (t+1)-dimensional thirteen-shift table is exactly
  -13^(-(t+1)) times the integral over F of
  (2-d_deep^+) product_a(11+d_a^+). Thus some character nonzero in
  every coordinate has negative real part at most
  -11^t mu(F)/(13^(t+1)12^(t+1)), and the fully mixed Fourier energy
  is at least
  11^(2t)mu(F)^2/(13^(2(t+1))12^(t+1)). The finite transform is the
  joint Poisson-Abel boundary of a collapsed f-weighted
  deep/blocker-probe series, whose displayed live probe multipliers are
  separately coprime to 91; the literal word remains collapsed inside
  its indicator. THM-2354
  supplies such an F on all 165 rows. The exact pure-word coefficient
  floor is 2593 e_j/1195871040; the exact fork floor is
  28523 e_j/186555882240. On a pure word the blocker shift is the
  actual redundant named-factor shift; on a fork the two shifts are
  auxiliary duplicate-blocker coefficient functionals. No target
  quotient, relation address, endpoint phase, scalar-row exclusion,
  profile decrement, or LRC(14) consequence follows.
source: codex-2026-07-25-mixed-deep-blocker-corner
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2354-deep-shift-comb-cover-and-grouped-unit-current
related:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - THM-2362-thirteen-shift-successor-statistic-and-role-jet-floor
script: 04-computation/mixed_deep_blocker_corner_thm2364.py
output: 05-knowledge/results/mixed_deep_blocker_corner_thm2364.out
script_sha256: 9888d27fef54d9e1f6dc6d8cc9eeda541c0a3eac029e2da74949d64f5b5733df
output_sha256: 8d78a2ab25b94db1bf7a1ddd94893f312a25dd17b18abf3c5f6e567b03ce7bc2
hash_basis: working-tree bytes (LF)
---

# THM-2364 -- every named blocker gets one mixed probe colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2354 forces a nonzero colour after shifting the deepest danger comb.
The factor-role calculation behind THM-2362 forces nonzero colours after
shifting anchored blocker dangers. Separate certificates could live on
disjoint Fourier supports.

Retaining every shift coordinate before selecting a colour removes that
ambiguity. The deepest anchored zero and all named blocker anchored ones
force one coefficient which is nonzero in every coordinate. For a fork,
this gives a single deep-by-`a`-by-`b` coefficient rather than two
unrelated pair certificates.

## 1. The full blocker-probe shift table

On the circle put

```text
d_0(y)=1_(||y||<1/14),

D_(c,r)(x)=d_0(c x-r/13),       r in F_13.           (1)
```

Let `t in {1,2}`. Let `c,d_1,...,d_t` be positive integers and let `F`
be measurable with

```text
rho=mu(F)>0,

F intersection D_(c,0)=empty,

F subset intersection_(a=1)^t D_(d_a,0),            (2)
```

up to null boundaries. For

```text
r,s_1,...,s_t in F_13
```

define the physical overlap table

```text
H(r,s_1,...,s_t)
 =mu(
    F intersection D_(c,r)
      intersection intersection_a D_(d_a,s_a)
   ).                                               (3)
```

The first coordinate is the deepest comb. The remaining coordinates
insert one translated probe for every label in the nonempty blocker
word. The original word factors remain collapsed in `1_F`.

## 2. Pointwise corner factorization

For any function `V` on `F_13`, its normalized finite transform satisfies

```text
sum_(u!=0)Vhat(u)=V(0)-1/13 sum_r V(r).              (4)
```

Away from the pullbacks of the `26` base strict-open interval endpoints,
the exact thirteen-shift danger count is

```text
sum_r D_(q,r)(x)=2-D_(13q,0)(x).                    (5)
```

Thus (5) holds almost everywhere, which is the only sense used in all
integrals below. At an endpoint its two sides need not agree, but the
finite exceptional set has measure zero.

On `F` away from those null boundaries, the deepest zero in (2) makes
the right side of (4)

```text
-(2-D_(13c,0)(x))/13.                               (6)
```

For each named blocker, the anchored one in (2) instead gives

```text
1-(2-D_(13d_a,0)(x))/13
 =(11+D_(13d_a,0)(x))/13.                           (7)
```

Let

```text
Hhat(u_0,...,u_t)
 =13^(-(t+1))
   sum_(r,s_1,...,s_t)
   H(r,s_1,...,s_t)
   zeta^(u_0 r+sum_a u_a s_a),                     (8)

zeta=exp(2*pi*i/13).
```

Finite inversion may be performed pointwise almost everywhere inside the
integral in (3). Multiplying (6)--(7) and integrating yields the exact
fully mixed corner

```text
S_t
 :=sum_(u_0!=0,...,u_t!=0)Hhat(u_0,...,u_t)

 =-13^(-(t+1))
   integral_F
    (2-D_(13c,0))
    product_(a=1)^t(11+D_(13d_a,0)).                (9)
```

Every integrand in (9) is positive. More precisely,

```text
S_t<=-11^t rho/13^(t+1)<0.                          (10)
```

There are `12^(t+1)` fully mixed colours. Hence some tuple with every
`u_a!=0` satisfies the signed bound

```text
Re Hhat(u_0,...,u_t)
 <=-11^t rho/(13^(t+1)12^(t+1)).                   (11)
```

Cauchy--Schwarz gives

```text
sum_(all u_a!=0)|Hhat(u_0,...,u_t)|^2
 >=|S_t|^2/12^(t+1)
 >=11^(2t)rho^2/(13^(2(t+1))12^(t+1)).             (12)
```

Thus one coefficient, not merely one marginal in each direction, carries
the entire nonempty blocker word.

## 3. The exact one-blocker successor refinement

For `t=1`, write

```text
G_tot
 :=sum_r mu(F intersection D_(c,r))
 =integral_F(2-D_(13c,0)),

J
 :=integral_F
    (2-D_(13c,0))D_(13d_1,0).                      (13)
```

Then (9) becomes

```text
S_1=-(11G_tot+J)/169.                               (14)
```

If

```text
rho_3^+=mu(F intersection D_(13c,0)),

rho_1^+=mu(F intersection D_(13d_1,0)),
```

then

```text
G_tot=2rho-rho_3^+,

J>=rho_1^+.                                        (15)
```

Consequently

```text
S_1<=-(11rho+rho_1^+)/169.                          (16)
```

This is the exact synthesis of THM-2354's cover mass with THM-2362's
successor bit. Equality in the universal `t=1` floor holds exactly when

```text
D_(13c,0)=1,

D_(13d_1,0)=0
```

almost everywhere on `F`: the deepest cover then has multiplicity one
and the named-blocker successor overlap vanishes.

## 4. Equality boundary of the profile argument

The constants in (10)--(12) are sharp given only the abstract anchored
profile facts. They are not asserted sharp on the one-dimensional LRC
geometry.

Define one-coordinate profiles by

```text
a(0)=0,             a(r)=1/12 for r!=0,

b(0)=1,             b(s)=1/12 for s!=0,             (17)
```

and put

```text
H(r,s_1,...,s_t)
 =rho a(r) product_a b(s_a).                        (18)
```

The deepest profile has total one; each blocker profile has total two.
Every fully mixed coefficient in (18) is the same negative real number

```text
-11^t rho/(13^(t+1)12^(t+1)).                      (19)
```

Thus equality holds in (11)--(12). This is an abstract profile hostile,
not a claim that (18) is realized by a canonical LRC row.

## 5. Lawful joint Abel grouping

Assume now that `F` is a finite union of rational intervals, as in the
LRC application, and put `f=1_F`. This keeps the complete literal word
collapsed as one support indicator; only the added deepest and blocker
probe factors are expanded below. With

```text
h_hat(n)=integral_T h(x)exp(-2*pi*i n x)dx,
```

the danger coefficients are

```text
(d_0)_hat(0)=1/7,

(d_0)_hat(m)=sin(pi*m/7)/(pi*m)       for m!=0.      (20)
```

For `0<tau<1`, Poisson-smooth all `t+1` retained danger factors with the
same parameter before multiplying by `f`. At fixed `tau`, the
multi-Fourier expansion is absolutely convergent. For a fully mixed
colour tuple it is

```text
Hhat_tau(u_0,...,u_t)
 =sum_(m_a=u_a mod 13 for every a)
   tau^(sum_a |m_a|)
   f_hat(-m_0 c-sum_(a=1)^t m_a d_a)
   product_(a=0)^t (d_0)_hat(m_a).                  (21)
```

Poisson smoothing preserves `0<=d_(0,tau)<=1`. Bounded-product `L1`
convergence gives

```text
lim_(tau->1-)Hhat_tau(u_0,...,u_t)
 =Hhat(u_0,...,u_t).                                (22)
```

This defines the joint Abel grouping before any infinite rearrangement.
No absolute convergence or reordering of the undamped multi-series is
claimed.

Every selected residue `u_a` is nonzero modulo thirteen. Formula (20)
kills multiples of seven, so every live multiplier in (21) satisfies

```text
gcd(m_a,91)=1.                                      (23)
```

The coefficient forced by (11) is therefore a collapsed `f`-weighted
deep/blocker-probe `91`-unit Abel coefficient: two added
probe-multiplier groups for a pure word and three for a fork. Only these
displayed probe multipliers are proved separately `91`-unit. They are
not the original Bockstein word indices hidden inside `f`.

## 6. Uniform application to all 165 rows

Apply THM-2354 on any first-depth-one scalar row. It supplies a delayed
literal word

```text
F=E_j intersection T^(-k)Q_(j,sigma),

rho>=e_j eta/6,

eta=2593/90090,                                    (24)
```

where

```text
sigma in {{a},{b},{a,b}}
```

is nonempty. Put

```text
R=13^k,

t=|sigma|.
```

For every blocker label `a in sigma`, the terminal word contains
`D_(c_a,0)`, so its pullback gives

```text
F subset D_(R c_a,0).                              (25)
```

THM-2354 also gives

```text
F intersection D_(c_3,0)=empty.                   (26)
```

Use Sections 1--5 with

```text
c=c_3,

d_a=R c_a.                                         (27)
```

For a pure word, `t=1`, and some fully mixed pair has

```text
Re Hhat(u_0,u_1)
 <=-2593 e_j/1195871040.                            (28)
```

For a fork, `t=2`, and some deep-by-`a`-by-`b` colour has

```text
Re Hhat(u_0,u_1,u_2)
 <=-28523 e_j/186555882240.                         (29)
```

These statements hold at a finite coefficient-dependent clock on every
one of the `165` rows. They retain the literal word as a collapsed
support label, the selected owner and clock, the deepest probe colour,
one added probe colour for every named blocker, and the separate
`91`-unit property of every added multiplier group.

## 7. Pure words versus forks

The blocker coordinates have two different physical typings.

For a pure word `sigma={a}`, THM-2305's scalar cover gives

```text
A_0 intersection D_(c_j)^c intersection D_(c_b)^c
 subset D_(c_a).                                   (30)
```

Thus the displayed `D_(c_a)` factor is redundant. After pullback, the
profile `s -> H(r,s)` can be obtained by translating the actual named
blocker factor while all remaining word factors stay fixed. In this
branch, (28) is a genuine deep-colour by named-factor-colour response.

For a fork `sigma={a,b}`, deleting either danger factor enlarges the word:
the other danger already satisfies the scalar cover. Definition (3)
instead retains the literal fork and inserts translated duplicates of
both named blocker dangers. The auxiliary indices in (29) are not
THM-2337's `beta_a,b`, Bockstein-jet colours, or relation-address
coordinates. The coefficient is a lawful coefficient-derived
duplicate-factor functional, not an observable already supplied by
canon. A physical use requires explicit duplicate pair probes or an
equivalent graph-channel phase-ratio sidecar.

## 8. Scope

The theorem aligns the full factor-coloured word, not the relation target.
It preserves

```text
literal delayed word,
owner and clock,
deep probe-multiplier residue,
every added named-blocker probe residue,
separate 91-unit probe-multiplier conditions,
signed quantitative fully mixed coefficient.
```

It still forgets

```text
exact endpoint frequencies,
full relation address,
target quotient,
absolute terminal-component phase,
planar graph label,
bounded scalar-frequency visibility.
```

In particular, `R beta` remains target-neutral modulo thirteen, and the
relative deepest shift omits the endpoint phase. The theorem does not
produce a nonzero scalar-frequency survivor, realize THM-2356's missing
graph-channel phase ratio, exclude a scalar row, decrement a profile, or
prove LRC(14).

Its gain is precise: the collapsed support of the complete nonempty
blocker hyperedge now carries one signed, fully mixed, all-unit added-probe
coefficient on every row. It does not expand or identify the original
word indices. The next bridge must transport that probe coefficient into
the target/address channel without destroying the anchored overlap.

## 9. Exact companion

The dependency-free companion uses `Fraction` arithmetic to:

- exhaust the `26` open cells in the successor count (5);
- test `64` exact anchored two-coordinate profiles;
- test `64` exact one- and two-blocker successor-bit mixtures and their
  Cauchy lower-bound arithmetic;
- verify the exact one-blocker refinement (14)--(16);
- enumerate the `144` pure and `1,728` fork coefficients of the sharp
  factorized profiles and verify their actual mixed energies;
- enumerate `5,184` pure and `373,248` fork ordered `91`-unit residue
  tuples;
- realize both sharp abstract profiles (17)--(19); and
- reduce both all-row constants in (28)--(29) exactly.

Run

```bash
python3 04-computation/mixed_deep_blocker_corner_thm2364.py
python3 -O 04-computation/mixed_deep_blocker_corner_thm2364.py
```

Both transcripts must match

```text
05-knowledge/results/mixed_deep_blocker_corner_thm2364.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent hostile audit checked the almost-everywhere endpoint scope,
the exact one- and two-blocker corner formulas and constants, the
collapsed-word versus added-probe typing, pure/fork boundary, joint Abel
limit without undamped reordering, both equality profiles, all-row
quantifiers, normal and optimized transcripts, stored output, LF hashes,
and documentation routing. QED.
