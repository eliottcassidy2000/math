---
id: THM-2571
title: "Deep-colour Cayley filling, primitive carry Bockstein, and norm-curvature split"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Every THM-2567 deep-colour augmentation-zero family is an exact rational
  boundary for THM-2532's Cayley derivative.  On the integral augmentation
  lattice its sole obstruction is the cyclic first moment modulo 13.  A
  lawful mass-one singleton realizes the nonzero obstruction with sharp
  denominator 13.  More decisively, the canonical K_11,9 carry tensor has,
  after one global primitive clearing, nonzero Bockstein in all 78 nonzero
  owner x target profiles.  The class factors as the 13-cyclotomic socle
  Omega times a unit Y in the septimal owner algebra.  It is target-flat,
  while the duty packet is plaquette-flat and can still have positive
  Galois norm.  These are exact coefficient/common-carrier statements, not
  a positive physical filling or a canonical semantic endpoint.  No scalar
  row is excluded and LRC(14) remains open.
source: root-holotopy-2026-07-27-cayley-bockstein
depends_on:
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2567-deep-coloured-duty-replica-cycle-and-augmentation-cancellation
related:
  - THM-2550-canonical-typed-row-double-nondegeneracy
  - THM-2562-canonical-duty-commutator-line-rank-and-anchor-rigidity
  - THM-2572-deep-augmentation-parseval-energy-and-nonlinear-holotopy-obstruction
script: 04-computation/lrc14_deep_colour_cayley_bockstein_thm2571.py
output: 05-knowledge/results/lrc14_deep_colour_cayley_bockstein_thm2571.out
script_sha256: 456c04fafac4242bd90041051928dac6f557e611d10f282e1202f283d5ce985d
output_sha256: 6fe025f7ebe43ae07fa054e22aa5a271c78731f5a1ee66e6c873ea91c617c9b7
hash_basis: LF-normalized bytes
---

# THM-2571 -- the coloured cycle is rationally fillable but integrally charged

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

THM-2567 exhibits a complete deep-colour cycle: every nonzero colour is
nonzero, but summing all thirteen colours gives zero.  THM-2532 supplies an
exact sawtooth inverse for the cyclic Cayley derivative.  Putting those two
objects on the same colour module gives a sharper answer than either alone:

```text
deep augmentation-zero cycle
        -> rational sawtooth filling
        -> one integral first-moment obstruction modulo 13
        -> a nonzero primitive class on the canonical K_11,9 carry tensor.
                                                                    (1)
```

Thus the ordinary rational homology vanishes.  The surviving class is
13-primary and integral.  It is different again from THM-2572's positive
quadratic energy: the latter can be nonzero when the first-moment class is
zero.  This rational/integral/nonlinear trichotomy is the precise holotopy
content of the theorem.

The common-carrier normalization is load-bearing.  A rational tensor is
cleared once, globally, to a primitive integral tensor.  Clearing each
Fourier profile separately would change the lattice, while inserting an
extra factor `13` would kill the obstruction by definition.

## 1. Rational Cayley filling of every deep cycle

Put `p=13`, let `E=Q(zeta_13)`, and write

```text
U_E={a=(a_m)_(m in F_13) in E^13 : sum_m a_m=0}.            (2)
```

For `tau in F_13^*`, let `P_tau` be cyclic translation in the colour
coordinate and let

```text
C_tau=sum_(d=1)^12 (-1)^(d+1) P_tau^d,                     (3)

B_tau=1/13 sum_(d=1)^12 (2d-13) P_tau^d.                   (4)
```

THM-2532 proves

```text
C_tau B_tau=B_tau C_tau=I-J/13.                            (5)
```

Therefore every `a in U_E` has the exact rational primitive

```text
a=C_tau(B_tau a).                                          (6)
```

Apply this coefficientwise to the THM-2567 target functions
`A_m(q)`.  Its deep augmentation law is

```text
sum_m A_m(q)=0                         for every q.          (7)
```

Consequently

```text
A(q)=C_tau B_tau A(q)                  for every q.          (8)
```

The target-duty operator `K_g` acts on `q`, whereas `B_tau` and `C_tau`
act on `m`.  They commute.  Hence the whole six-replica family is also an
exact rational boundary:

```text
K_g A=C_tau B_tau(K_g A).                                  (9)
```

There is no rational linear homology left in the deep-colour direction.
This does not make `B_tau` a physical positive filling: its sawtooth kernel
is signed and has denominator `13`.

## 2. The integral augmentation lattice and its Bockstein

Let

```text
O=Z[zeta_13],
Lambda_O={a in O^13 : sum_m a_m=0}.                         (10)
```

Define the cyclic first moment

```text
beta(a)=sum_(m=0)^12 m a_m                    in O/13O.     (11)
```

Then there is an exact sequence

```text
0 -> Lambda_O --C_tau--> Lambda_O --beta--> O/13O -> 0.    (12)
```

Equivalently,

```text
a has an integral Cayley primitive
<=> beta(a)=0
<=> B_tau a lies in Lambda_O.                              (13)
```

### Proof

THM-2532 computes the Smith form of `C_tau` on the integral augmentation
lattice as

```text
(1,1,1,1,1,1,1,1,1,1,1,13).                              (14)
```

Base extension from `Z` to its free `Z`-module `O` leaves that diagonal
form unchanged.  Thus the cokernel is `O/13O`.

For `x in Lambda_O`, cyclic translation preserves the first moment modulo
`13`:

```text
beta(P_tau^d x)=beta(x)-d tau sum_m x_m=beta(x).            (15)
```

The twelve alternating coefficients in (3) sum to zero, so
`beta(C_tau x)=0`.  The moment map is onto, because
`beta(e_1-e_0)=1`.  Its kernel and `C_tau Lambda_O` therefore have the same
index and are equal.  Equation (13) follows from the uniqueness of the
rational primitive on the augmentation subspace. QED.

The class is independent of the chosen origin in `F_13`: translating the
indices changes (11) by a multiple of `sum a_m=0`.  Multiplying the index by
a unit rescales `beta` by a unit, and reversal changes its sign.  Thus
vanishing is chart-gauge invariant, although a sign or oriented primitive
is not canonically chosen.

## 3. Physical coefficient formula and a sharp singleton

For a rational lawful table `H(r,s,t)`, use its global primitive clearing:
write `H=h/D` with `D>0`, integral `h`, and
`gcd(D,{h(r,s,t)})=1`, and put

```text
a_m(q)=D 13^3 A_m(q)
      =sum_(r,s,t) D H(r,s,t)
         zeta_13^[m(r-t)+q_1 s+q_2 t].                     (16)
```

The diagonal law `H(t,s,t)=0` gives `sum_m a_m(q)=0`.  Since

```text
sum_(m=0)^12 m zeta_13^(mu)=13/(zeta_13^u-1),   u!=0,      (17)
```

the integral obstruction has the exact carrier formula

```text
beta_q
 =13 sum_(r,s,t: r!=t) D H(r,s,t)
       zeta_13^(q_1 s+q_2 t)/(zeta_13^(r-t)-1)   mod 13O.  (18)
```

Every term in (18) is integral.  Formula (18) is a displacement-weighted
quantity; it retains the oriented difference `r-t` that ordinary deep
augmentation forgets.

The mass-one table

```text
H=delta_(r,s,t)=(1,1,0)                                   (19)
```

obeys the diagonal and both coordinate-zero laws used by THM-2563.  At
every target `q`,

```text
a_m(q)=zeta_13^(q_1+m),

beta_q=zeta_13^q_1 Omega !=0 mod 13O,                      (20)

Omega=sum_(m=0)^12 m zeta_13^m=13/(zeta_13-1).             (21)
```

Indeed `(zeta_13-1)Omega=13`; the cyclotomic valuation of `Omega` is
eleven, while that of `13` is twelve.  Thus `Omega` is nonzero modulo
`13O` and spans the one-dimensional socle of

```text
O/13O isomorphic to F_13[epsilon]/(epsilon^12),
epsilon=zeta_13-1.                                         (22)
```

The sawtooth primitive of (19) has denominator exactly `13`.  This proves
both nontriviality and sharpness of the denominator in (4).

## 4. Canonical K_11,9 carry class

The abstract singleton shows that the obstruction can occur.  The exact
old/future diagonal carrier underlying THM-2550(B) shows that it actually
occurs on the canonical typed row.

Use the fixed denominator and clock

```text
T=297836897838480,
R=13^6=4826809,
RT=1437601819018855810320.                                 (23)
```

The companion reconstructs, by exact interval sweeps, the rational tensor

```text
X_(ell,s,r),       ell in F_7, s,r in F_13,                (24)
```

formed by the THM-2550(B) old/future diagonal word, its seven source
phases, thirteen target shifts, thirteen deep-root shifts, and thirteen
base digits.  Exactly `1092=7*13*12` of its `1183` cells are positive:

```text
X_(ell,s,r)>0 <=> r!=0.                                    (25)
```

The integer numerators over denominator `RT` have global gcd exactly
`13`.  Divide once, globally.  The primitive tensor

```text
x_(ell,s,r)=D X_(ell,s,r) in Z,

D=RT/13=110584755309142754640                              (26)
```

has serialization digest

```text
a66ba96d31a33354468392b1dabc19865e6e925158efdab059fad9a98d4390f4.
                                                                    (27)
```

Let `M=Z[zeta_7,zeta_13]`.  For every nonzero owner colour
`kappa in F_7^*` and target colour `b in F_13`, define the deep cycle

```text
c_m^(kappa,b)
 =sum_(ell,s,r) x_(ell,s,r)
    zeta_7^(kappa ell) zeta_13^(b s+m r)       in M.        (28)
```

Equation (25) and root orthogonality give

```text
sum_m c_m^(kappa,b)=0.                                     (29)
```

Its first-moment class is nonzero for **all 78** pairs `(kappa,b)`:

```text
beta_(kappa,b)=sum_m m c_m^(kappa,b) !=0 in M/13M.          (30)
```

This is not a 78-case accident.  In `O/13O`, put

```text
Omega=(1,2,3,4,5,6,7,8,9,10,11,12)                        (31)
```

in the power basis `1,zeta_13,...,zeta_13^11`.  Then

```text
sum_m m zeta_13^(mr)=r^(-1) Omega,
zeta_13^a Omega=Omega.                                     (32)
```

The primitive carrier contracts to the single septimal polynomial

```text
Y(z)=6z+5z^2+z^3+12z^4+8z^5+7z^6        in F_13[z].        (33)
```

Exact reduction gives

```text
gcd(Y,Phi_7)=1,

Y^(-1)=3z^2+5z^3+8z^4+10z^5       mod (13,Phi_7).          (34)
```

Consequently

```text
beta_(kappa,b)=Omega Y(zeta_7^kappa),                       (35)
```

independently of `b`.  The six owner factors, in the standard septimal
power basis, are

```text
( 6,12,11, 7, 5,1), (12,11,5,7,4,6), (8,3,9,1,2,7),
( 5,10, 4,12,11,6), ( 1, 2,8,6,9,7), (7,1,2,6,8,12).       (36)
```

Every entry in (31) and (36) is nonzero, so every class in (30) has full
`6*12=72` product-basis support.  Before owner/target contraction, the
deep Bockstein is nonzero in `88/91` `(ell,s)` cells; the only zeros are

```text
(ell,s)=(0,0),(3,0),(4,0).                                 (37)
```

Owner contraction removes those exceptional directions because `Y` is a
unit.

The primitive convention in (26) is essential.  Multiplying every numerator
and the denominator by an extra `13` represents the same rational tensor but
makes (30) zero.  Such a rescaling is not an integral invariant of the
carrier.  Global primitive clearing is the canonical lattice; per-profile
clearing would likewise replace one common carrier by 78 unrelated lattices.

## 5. Why universality does not select the target

Equation (35) is simultaneously powerful and limiting.  The class survives
all six owner colours because `Y` is a unit, and it survives every target
modulation because `Omega` is the 13-cyclotomic socle.  Therefore it proves
an integral obstruction on each externally labelled target profile.  For
`b!=0`, that profile already has a nonzero first physical target coordinate.

But `beta_(kappa,b)` itself is independent of `b`.  It cannot distinguish
target phases, choose a semantic endpoint, or recover a prescribed target
from an unlabelled sum.  The same socle rigidity that makes the class
universal makes it target-flat.  The `b=0` profiles in particular do not
acquire a nonzero target merely from (30).

There is a second orientation boundary.  Algebraically, any `tau!=0` may be
used in (3).  On a physical `H` carrier, shifting the deep-colour index is
modulation by `zeta_13^[tau(r-t)]`; it is not a positive root translation.
Reversing `tau` reverses the primitive.  A common oriented physical sidecar
is still required before the signed filling itself can be used as a lawful
current.

## 6. Plaquette-flat duty and positive norm

Let

```text
a_m=A_m(0,0),                    m in F_13^*,               (38)
```

for a nonzero rational nonnegative lawful table.  THM-2567 proves every
`a_m` is nonzero.  For a gain `g!=0` and replica `j=1,...,6`, put

```text
R_(m;(g,j))=(K_g A_m)(jg)=-d_g a_m.                        (39)
```

The complete deep-colour x duty-replica matrix has rank one.  Every
`2 x 2` plaquette minor vanishes:

```text
R_(m;i)R_(n;j)-R_(m;j)R_(n;i)=0.                           (40)
```

Nevertheless each gain column has the positive nonlinear invariant

```text
prod_(m=1)^12 R_(m;(g,j))
 =d_g^12 prod_(m=1)^6 |a_m|^2 >0.                          (41)
```

The product is rational by Galois conjugacy.  Thus local rank-two curvature
is identically zero while the global Galois norm is positive.

The Bockstein and this norm/THM-2572 energy are independent diagnostics.
For the uniform off-diagonal hostile, after integral scaling,

```text
a_0=144,            a_m=-12 for m!=0.                      (42)
```

Its norm and Parseval energy are positive, but

```text
beta(a)=-12 sum_(m=1)^12 m=0 mod 13.                       (43)
```

The singleton (19) has nonzero `beta`.  Therefore neither adding more
proportional duty gains nor squaring the cycle recovers the missing integral
first moment.  Conversely, the integral class alone supplies no positive
Gram carrier.

## 7. Exact companion and proof boundary

Run

```bash
python3 04-computation/lrc14_deep_colour_cayley_bockstein_thm2571.py
python3 -O 04-computation/lrc14_deep_colour_cayley_bockstein_thm2571.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_deep_colour_cayley_bockstein_thm2571.out
```

byte-for-byte.  The standard-library exact referee uses only integers and
`Fraction` arithmetic.  It verifies:

- the Cayley inverse, determinant `13`, and mod-13 rank `11`;
- `4095` augmentation profiles, of which exactly `315` are integrally
  fillable and `3780` obstructed;
- the singleton denominator and `2028` physical displacement controls;
- `10098` rank-one plaquette controls and `4095` nonzero rational Galois
  profiles;
- the complete `1092`-cell canonical carry reconstruction, primitive gcd,
  valuation range, denominator, and digest;
- the `88/91` precontraction census, the closed `Omega*Y` factorization,
  explicit inverse (34), all `78/78` owner-target classes, and their full
  product-basis support.

There are `100860` explicit checks, none implemented with `assert`.

What is proved is an exact rational filling theorem, its complete integral
cokernel, a sharp mass-one obstruction, and a nonzero primitive obstruction
on one canonical typed-row carry carrier.  What is not proved is positivity
of the sawtooth primitive, a common THM-2562 Reynolds carrier, a selected
target phase, a semantic THM-2305 endpoint, a scalar relation, an owner/root
intertwiner, exclusion of a scalar-cover row, or LRC(14). **QED.**
