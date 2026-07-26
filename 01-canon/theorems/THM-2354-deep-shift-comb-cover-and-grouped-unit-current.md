---
id: THM-2354
title: "Deep-shift comb cover and grouped unit-current energy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Let
  F be a positive-measure subset of E, with E disjoint from one
  c-fold danger comb. The overlap profile of F with the thirteen
  translates of that comb is zero at the original translate and has
  total mass at least measure(F). Its nonconstant Fourier energy is at
  least measure(F)^2/2028, so some nonzero mod-thirteen colour has
  magnitude at least measure(F)/156. That coefficient is an absolutely
  convergent collapsed deep-harmonic series and also the Abel boundary
  of the full mixed word/comb/owner current grouped by the deep
  multiplier modulo thirteen; every live multiplier in it is coprime
  to 91. Applied after THM-2349, every one of the 165 rows has such a
  literal-word grouped current, with magnitude at least e_j*eta/936.
  This is a factor-coloured relative deep-leg probe, not a character
  of the relation address or target quotient. It excludes no scalar
  row and does not prove LRC(14).
source: codex-2026-07-25-deep-shift-comb-cover
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
related:
  - THM-2326-vertexwise-septimally-primitive-c3-degree
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
script: 04-computation/lrc14_deep_shift_grouped_current_thm2354.py
output: 05-knowledge/results/lrc14_deep_shift_grouped_current_thm2354.out
script_sha256: 790ef9d1300de572603089c186d75eae0127a14dbe70dcf7a99b49fd7d82683e
output_sha256: ce64853010708e24194f5f320a0c7774f2e7336e8b90c67c21066612c42da687
hash_basis: working-tree bytes (LF)
---

# THM-2354 -- a relative deep shift forces a grouped unit current

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Selecting one Fourier triangle early retains an exact marked atom but loses
the phase of every alternative edge. There is a complementary operation:
sum the entire mixed Fourier identity before selecting the deep multiplier.
The sum is not mysterious. It is the overlap of the literal word set with
a translated deepest danger comb.

Thirteen translates cover the circle. The unshifted overlap vanishes by
owner exclusivity, while the sum of all overlaps is positive. Finite
Fourier inversion therefore forces a nonzero deep colour. The quantitative
form is

```text
literal word mass mu
  -> nonnegative 13-shift overlap profile G with G(0)=0
  -> nonzero-colour energy at least mu^2/2028
  -> one grouped deep colour of magnitude at least mu/156.         (1)
```

This is the acute grouped-current mechanism sought after THM-2344, but on
a deliberately factor-coloured coordinate. It preserves the deep
multiplier and literal word while forgetting the endpoint contribution to
the relation address. That loss is both the theorem's usefulness and its
sharp boundary.

## 1. The thirteen-translate comb cover

On the circle put

```text
d_0(y)=1_(||y||<1/14),

d_r(y)=d_0(y-r/13),                 r in F_13,

D_(c,r)(x)=d_r(c x),                c>=1.            (2)
```

Let `F,E` be measurable circle sets satisfying

```text
0<mu(F),
F subset E,
E intersection D_(c,0)=empty                      (3)
```

up to null boundaries. Define the real overlap profile

```text
G(r)=mu(F intersection D_(c,r)).                    (4)
```

Then

```text
G(r)>=0,
G(0)=0.                                             (5)
```

The thirteen arcs underlying `d_r` have centres spaced by `1/13` and
length `1/7`. Since

```text
1/13<1/7,
```

they cover the circle. Pullback by `x -> c x` gives the pointwise
almost-everywhere inequality

```text
sum_(r in F_13) D_(c,r)(x)>=1.                      (6)
```

Multiplying by `1_F` and integrating yields

```text
sum_r G(r)>=mu(F)>0.                                (7)
```

Thus `G` is nonzero and nonconstant.

The exact arrangement has denominator `182`. Each translate contains
`26` open cells; `26` cells are covered once and `156` twice. Every
translate has a two-cell unique core. These counts will give both a
one-sparse hostile and equality in the quantitative bound.

## 2. The sharp nonzero-colour energy

Let

```text
zeta=exp(2*pi*i/13)
```

and use the normalized transform

```text
C(a)=1/13 sum_(r in F_13) zeta^(a r) G(r).          (8)
```

Write

```text
G_bar=C(0)=1/13 sum_r G(r).
```

Normalized Parseval gives

```text
sum_(a!=0)|C(a)|^2
 =1/13 sum_r |G(r)-G_bar|^2.                        (9)
```

Because `G(0)=0`, Cauchy--Schwarz on the other twelve entries gives

```text
sum_(r!=0)G(r)^2
 >=(sum_r G(r))^2/12.
```

Substitution in (9), followed by (7), gives

```text
sum_(a!=0)|C(a)|^2
 >=G_bar^2/12
 >=mu(F)^2/(13^2*12)
 =mu(F)^2/2028.                                    (10)
```

There are twelve nonzero colours, so

```text
max_(a!=0)|C(a)|
 >=mu(F)/sqrt(12*2028)
 =mu(F)/156.                                       (11)
```

Both constants are sharp. Choose one equal cell from the unique core of
each translate `r=1,...,12` and let `F=E` be their union. Then

```text
mu(F)=6/91,
G(0)=0,
G(r)=1/182                  for r!=0,

C(a)=-1/2366                for a!=0.               (12)
```

Thus equality holds in (10)--(11).

## 3. The overlap transform is the Abel-grouped mixed current

Assume now that `F,E` are finite unions of rational intervals, as in the
LRC application. Write

```text
f=1_F,                 e=1_E,

h_hat(n)=integral_T h(x)exp(-2*pi*i*n*x)dx.         (13)
```

The centred base danger interval has

```text
(d_0)_hat(0)=1/7,

(d_0)_hat(m)=sin(pi*m/7)/(pi*m)       for m!=0.     (14)
```

Direct `L2` Parseval gives the exact functional form

```text
G(r)
 =sum_(m in Z)
    f_hat(-m c)(d_0)_hat(m)zeta^(-m r),             (14a)

C(a)
 =sum_(m=a mod 13)
    f_hat(-m c)(d_0)_hat(m).                        (14b)
```

Both series are absolutely convergent. Indeed, Cauchy--Schwarz bounds the
sum of

```text
|f_hat(-m c)(d_0)_hat(m)|
```

by the product of the two `l2` Fourier norms. Thus the grouped
deep-harmonic current itself needs no boundary convention.

There is nevertheless a load-bearing full mixed-current refinement which
retains the bare owner factor. Since `F subset E`,

```text
f e=f.
```

For `0<rho<1`, Poisson-smooth `f,e` in their physical frequency and
`d_r` in its base frequency before pulling back by `c`. Put

```text
G_rho(r)
 =integral_T f_rho(x)e_rho(x)d_(r,rho)(c x) dx.    (15)
```

All three factors lie in `[0,1]`, and their product converges almost
everywhere and in `L1` to

```text
1_F 1_E D_(c,r)=1_F D_(c,r).
```

Hence

```text
lim_(rho->1-)G_rho(r)=G(r).                         (16)
```

For fixed `rho`, the Fourier expansion is absolutely convergent and gives

```text
G_rho(r)
 =sum_(A,m in Z)
    rho^(|A|+|m|+|A+m c|)
    f_hat(A)(d_0)_hat(m)
    conjugate(e_hat(A+m c))
    zeta^(-m r).                                    (17)
```

Apply the transform (8) to (17). Character orthogonality selects exactly
one deep colour:

```text
C_rho(a)
 =sum_(A in Z; m=a mod 13)
    rho^(|A|+|m|+|A+m c|)
    f_hat(A)(d_0)_hat(m)
    conjugate(e_hat(A+m c)),                        (18)

lim_(rho->1-)C_rho(a)=C(a).                         (19)
```

Thus `C(a)` is not merely evidence that one monomial exists. It is both
the absolutely convergent collapsed series (14b) and the boundary of the
entire word/deep-comb/owner current grouped by the actual deep multiplier
modulo thirteen. Formally summing the endpoint index in (18) recovers
(14b), because

```text
sum_A f_hat(A)conjugate(e_hat(A+m c))
 =(f e)_hat(-m c)
 =f_hat(-m c).                                      (19a)
```

For each fixed `m`, the first sum is absolutely convergent by
Cauchy--Schwarz. The full undamped double series is not being reordered;
its canonical meaning remains the Abel limit in (19).

For `a!=0`, every `m` in (18) is a thirteen-unit. Equation (14) kills
exactly the nonzero multiples of seven. Therefore every term with nonzero
deep coefficient satisfies

```text
gcd(m,91)=1.                                       (20)
```

Combining (11), (18), and (20):

> Some nonzero deep colour has an Abel-grouped `91`-unit current of
> magnitude at least `mu(F)/156`.

The series at the boundary need not be absolutely convergent; (19) is an
Abel statement. No regrouping of an undamped conditional series is used.

## 4. Uniform application to all 165 rows

Use THM-2349 on any one of the `165` first-depth-one scalar rows. Let
`E_j` be the positive universal depth-one owner; on a repeated-first row
either shallow owner may be used. Put

```text
e_j=mu(E_j)>0,
eta=2593/90090.
```

At the delayed clock `k`, THM-2349 proves

```text
mu(E_j intersection T^(-k)R_j)>=e_j eta/2.         (21)
```

The static residual partition has three terminal words. Choose a largest
word `Q_(j,sigma)` and put

```text
F=E_(j,sigma,k)
 =E_j intersection T^(-k)Q_(j,sigma).
```

Then

```text
mu(F)>=e_j eta/6.                                   (22)
```

Owner exclusivity gives

```text
F subset E_j,
E_j intersection D_(c_3,0)=empty.                 (23)
```

Apply Sections 1--3 with `E=E_j` and `c=c_3`. For some `a!=0`,

```text
|C(a)|
 >=e_j eta/936
 =2593 e_j/84324240
 >0,                                               (24)
```

and `C(a)` is the Abel-grouped literal-word/deepest-comb/bare-owner
current (18), with every live deepest multiplier a `91`-unit.

This applies at a finite coefficient-dependent clock on all `165` rows.
It gives a grouped noncancellation statement which THM-2349's selected
triangle alone did not provide.

## 5. The exact target-quotient loss

Equation (18) groups by the deep harmonic `m mod 13`. It does **not**
group by the relation address

```text
r_full=u+R beta+m e_(c_3)-v.                        (25)
```

A target character sees

```text
ell.r_full
 =ell.u+ell.(R beta)+m ell_(c_3)-ell.v,             (26)
```

whereas the relative deep shift in (17) sees only `m`. The endpoint
phases `ell.u-ell.v` are absent, and `R beta` is target-neutral only after
the factor-coloured split is retained. Consequently:

```text
preserved:
  literal delayed word, full endpoint-frequency sum,
  deep multiplier colour, septimal nonvanishing, quantitative magnitude;

forgotten:
  exact X and m, exact relation address, target quotient,
  planar target--jet graph, bounded visibility, terminal-component phase.
                                                                  (27)
```

The loss is structural. In the canonical factorization, co-shifting the
present and bare deepest-coordinate complements together with the danger
comb would insert the missing endpoint character, but then the translated
complement and danger factors have product zero pointwise. The nonzero
signal comes precisely from **relative** deep-leg drift. A naive lawful
co-shift annihilates it.

Thus (24) does not settle THM-2334's target variance or THM-2356's planar
graph singleton boundary. A required next sidecar is an intertwiner which
couples this relative deep colour to the target/jet character without
destroying the overlap cone, or a coefficient-sensitive comparison showing
that the omitted endpoint phases cannot cancel (24).

## 6. Needle geometry and hostile controls

The thirteen arcs are a one-dimensional needle-cover object: the source is
one danger interval, the operation is translation by `F_13`, and the
preserved predicate is point coverage. There is no intrinsic pairwise
orientation, so a tournament would add no information. This is also not a
planar Kakeya theorem: direction has been quotiented out, and the useful
coordinate is the translate label together with the unique/double-cover
cell sidecar.

Two exact controls mark the boundary.

1. **One-sparse overlap.** Take one of the two unique cells of a nonzero
   translate and let `F=E` be that cell. Then `G` is supported at exactly
   one nonzero shift. All twelve nonzero Fourier colours survive, but the
   shift profile is one-sparse. Grouped deep-colour noncancellation alone
   therefore does not defeat the one-sparse planar-graph boundary.

2. **Sharp energy.** The twelve-cell construction in (12) realizes equality
   in (10)--(11). No better universal constant follows from only
   nonnegativity, `G(0)=0`, and the translate cover.

Neither control is asserted to be a complete nine-coordinate scalar cover.
They are exact hostiles for strengthening the abstract consequence.

## 7. Exact companion and status

The dependency-free companion freezes:

- the exact `182`-cell arrangement, including `26` singly and `156`
  doubly covered cells;
- the two unique cells per translate;
- the one-sparse and sharp-equality controls;
- the constants `2028`, `156`, and the LRC factor `936`; and
- all `84` nonzero mod-thirteen deep residues, split into `12` septimal
  zeros and `72` `91`-units.

Reproduce with

```bash
python3 04-computation/lrc14_deep_shift_grouped_current_thm2354.py
python3 -O 04-computation/lrc14_deep_shift_grouped_current_thm2354.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_deep_shift_grouped_current_thm2354.out
```

byte-for-byte after LF normalization.

No scalar row is excluded. The exact ledger remains `165`, and LRC(14)
remains open.

## 8. Independent audit

The independent audit rederived the `182`-cell cover, the sharp
`2028/156` energy constants, the Fourier sign selecting `m=a mod 13`,
the absolute convergence of the collapsed series, and the separate Abel
meaning of the raw endpoint double series. It checked the largest-word
floor `e_j eta/6`, the all-`165` quantifiers, the split into `12`
septimal zeros and `72` unit residues, and both exact transcripts.

It also audited the quotient boundary against THM-2356: `C(a)` is one
factor-colour aggregate, not a target/jet chirp table or a planar-graph
singleton location. The result is compatible with all final target mass
remaining at zero. This is why (24) is a proved grouped noncancellation
sidecar rather than an LRC profile exclusion.
