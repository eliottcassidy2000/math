---
id: THM-2508
title: "Affine cut-bundle covariance, primitive-mode intertwining, and the invariant-scalar no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The 42 ordered affine
  cuts of F_7, over all twelve nonzero F_13 toothpick slopes, form the exact
  gauge-stable completion of one THM-2507 detector.  Every physical CRT
  affine normalization acts by a permutation of slope/cut/output entries;
  in root Fourier coordinates it adds only the correct unit phase.  Fourier
  transform in the cut-translation label diagonalizes the construction:
  each primitive cut coefficient is a nonzero seven-term geometric factor
  times the corresponding THM-2506 mixed F_13 x F_7 coefficient.  Hence all
  5,184 primitive cut coefficients survive on every essential punctured
  stalk.  The zero cut character vanishes identically, and transitivity of
  the affine gauge makes the uniform sum the only invariant linear scalar;
  it too vanishes.  Thus the cut character is a necessary linear sidecar,
  while the full gauge-invariant scalar quadratic energy loses target phase.
  This closes the
  static affine carry/coherence seam, not the temporal ancestry transplant,
  owner/arrival/deep-current coupling, any live scalar row, or LRC(14).
source: codex-2026-07-27-affine-cut-bundle; mac-mini-2026-07-27 exact factorization
depends_on:
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
related:
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
script: 04-computation/lrc14_affine_cut_bundle_covariance_thm2508_referee.py
output: 05-knowledge/results/lrc14_affine_cut_bundle_covariance_thm2508_referee.out
script_sha256: 044025e9facb276fff6b966df38ad35abcca6f578e8845e14738d6e2d3cc76b0
output_sha256: f676e1e59f8ccedb495a69268ac9ad61be772b735c961fc18211f8bfd7bc9ef7
hash_basis: working-tree bytes (LF)
secondary_script: 04-computation/lrc14_affine_cut_bundle_probe.py
secondary_output: 05-knowledge/results/lrc14_affine_cut_bundle_probe.out
secondary_script_sha256: 99d7d7074cb18b1398e612e1d73eb620af055c66db15b71eade96ea84842824c
secondary_output_sha256: 9b701ce424cb87ed090db459501c22a827a0e6dee22054d1afa03129698d0fbb
---

# THM-2508 -- the ordered cut is a covariant primitive-mode bundle

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2507 finds the first nonaffine `13`-root observer of the punctured
`F_13 x F_7` stalk, but one ordered cut is moved by the physical affine
normalization.  The repair is finite and exact: retain the entire affine
torsor of cuts.  More strongly, Fourier transform in that retained cut label
recovers every primitive mixed character with a nonzero diagonal multiplier.

The theorem therefore separates two statements that had previously been
conflated:

```text
static affine carry coherence:             CLOSED by the cut bundle;

same-ancestry owner/arrival/deep current:   OPEN.                 (1)
```

## 1. The full cut bundle

Let `d:F_13 x F_7 -> Q` obey the row-zero law

```text
sum_(r in F_7)d(h,r)=0                     for every h.            (2)
```

Write `rep(r)` for the representative of `r` in `{0,...,6}`.  For

```text
tau in F_13^*,             a in F_7^*,             c in F_7,
```

define

```text
R_(tau,a,c)(v)
 =sum_(r in F_7)d(v-tau rep(ar+c),r).                           (3)
```

The pair `(a,c)` is an ordered affine cut, not a homomorphism
`F_7 -> F_13`.  There are

```text
6*7=42
```

cuts and `12*42=504` nonzero-slope components.

For each fixed cut `(a,c)`, equation (3) is THM-2507 after a permutation of
the seven columns.  Hence every nonvertical rational row-zero array is seen
by at least seven of its twelve slopes.  Pointwise,

```text
#{(tau,a,c):R_(tau,a,c)!=0}>=42*7=294.                         (4)
```

The common kernel of the full bank is still exactly the six-dimensional
space

```text
d(h,r)=b(r),                    sum_r b(r)=0.                    (5)
```

Thus adding the cut sidecar restores covariance without enlarging the
invisible kernel.

## 2. Exact affine covariance

Let

```text
g=(A,H;B,C) in AGL_1(F_13) x AGL_1(F_7),

A in F_13^*, H in F_13, B in F_7^*, C in F_7,
```

act by pullback

```text
d^g(h,r)=d(A^(-1)(h-H),B^(-1)(r-C)).                            (6)
```

Substitute `r=Bs+C` in (3).  The representative is deliberately retained
*after* the affine `F_7` argument; no false linearization across its seam is
made.  One obtains

```text
R^(d^g)_(tau,a,c)(v)
 =R^d_(A^(-1)tau,aB,aC+c)(A^(-1)(v-H)).                         (7)
```

This is the exact carry law.  The row-dependent wrap that spoiled one fixed
cut in THM-2507 has become a permutation of the cut label.

If `zeta=zeta_13` and

```text
Theta^d_(tau,a,c)(alpha)
 =sum_v R^d_(tau,a,c)(v)zeta^(-alpha v),                         (8)
```

then (7) gives

```text
Theta^(d^g)_(tau,a,c)(alpha)
 =zeta^(-alpha H)
  Theta^d_(A^(-1)tau,aB,aC+c)(A alpha).                          (9)
```

Thus root colours are permuted and acquire precisely one unit translation
phase; neither their magnitude nor their nonvanishing is lost.

The physical THM-2436 gauge

```text
s -> U s+K                       mod 91
```

is a special case of (6), with

```text
A=U mod 13, B=U mod 7, H=K mod 13, C=K mod 7.                  (10)
```

This includes reversal, the normalized AP multiplier, and the common
quotient-root carry `kappa(Y)`.

The `42` cuts are not an arbitrary overbank.  Under the right action

```text
(a,c).(B,C)=(aB,aC+c),                                         (11)
```

`AGL_1(F_7)` acts freely and transitively.  Therefore all `42` cuts are the
smallest gauge-stable family containing the baseline cut `(1,0)`.  Including
the independent nonzero `F_13` multiplier makes the `504` triples
`(tau,a,c)` one regular component orbit.

## 3. Fourier transform in the cut is the missing intertwiner

Put `xi=zeta_7` and retain THM-2506's primitive stalk transform

```text
dtilde(alpha,gamma)
 =sum_(h,r)d(h,r)zeta^(-alpha h)xi^(-gamma r).                  (12)
```

Fourier transform (8) in the cut translation:

```text
Psi_(tau,a)(alpha,beta)
 =sum_(c in F_7)Theta_(tau,a,c)(alpha)xi^(-beta c).              (13)
```

Set `s=ar+c`.  Then `c=s-ar`, so exact reindexing gives

```text
Psi_(tau,a)(alpha,beta)
 =K_(alpha tau,beta)dtilde(alpha,-beta a),                       (14)

K_(u,beta)
 =sum_(s=0)^6(zeta^(-u)xi^(-beta))^s.                           (15)
```

There is no approximation or averaging loss in (14).  If
`alpha,tau!=0`, put

```text
lambda=zeta^(-alpha tau)xi^(-beta).
```

Then

```text
lambda^7=zeta^(-7 alpha tau)!=1.
```

Consequently `lambda!=1` and

```text
K_(alpha tau,beta)=(1-lambda^7)/(1-lambda)!=0.                  (16)
```

For every fixed `tau,a!=0`, equations (13)--(16) give a diagonal
isomorphism between the `72` mixed modes and the
`alpha!=0,beta!=0` cut-character part of the bundle.  (The six
`alpha=0,beta!=0` coordinates belong to the vertical/pure-`F_7` kernel.)
In particular, THM-2506 says that on every essential
THM-2436 defect

```text
Psi_(tau,a)(alpha,beta)!=0

for every tau,a,alpha,beta!=0.                                  (17)
```

There are exactly

```text
12*6*12*6=5,184                                                (18)
```

such coefficients.  The bundle does not merely find one surviving root
colour: it transports the complete primitive `91`-module into a root-colour
bank with the `F_7` information retained as cut phase.

Equation (9) also diagonalizes compatibly.  Summing after
`c'=aC+c` gives

```text
Psi^(d^g)_(tau,a)(alpha,beta)
 =zeta^(-alpha H)xi^(beta aC)
  Psi^d_(A^(-1)tau,aB)(A alpha,beta).                            (19)
```

The product `(A alpha)(A^(-1)tau)=alpha tau` keeps the geometric kernel
fixed, while `-beta a` transforms to `-beta aB`, exactly as the stalk
character requires.

## 4. The cut character is necessary, not cosmetic

At `beta=0`, equation (14) becomes

```text
Psi_(tau,a)(alpha,0)
 =K_(alpha tau,0)dtilde(alpha,0)=0                              (20)
```

by the row-sum law.  Thus summing away the cut label recreates precisely the
pure `F_13` pushforward no-go of THM-2506.

There is also a representation-theoretic sharpness statement.  Let

```text
I={(tau,a,c,v):tau in F_13^*, a in F_7^*, c in F_7, v in F_13}.
```

The physical product affine group acts transitively on `I` through (7):
given two entries, `A,B,C,H` can be solved successively from their
`tau,a,c,v` coordinates.  Hence an invariant linear functional on the
permutation module `Q^I` is a multiple of the uniform sum.  But for every
component,

```text
sum_v R_(tau,a,c)(v)=sum_(h,r)d(h,r)=0.                          (21)
```

Therefore every gauge-invariant linear scalar readout vanishes on the whole
Radon image.  The full gauge-invariant scalar quadratic norm of the bundle is
positive, but it pairs each root phase with its conjugate and so destroys the
target charge.  This does not say that every local quadratic contraction is
constant or root-neutral.  A successful linear observer must remain
nontrivial equivariant data; the cut-character bundle supplies one such
observer.

There is a complementary single-map no-go.  Suppose a map
`pi:F_13 x F_7 -> F_13` is strictly translation-equivariant in the sense

```text
pi(x+g)-pi(x)=chi(g)                  for every x,g,              (21a)
```

for one function `chi` independent of `x`.  Taking `x=0` gives
`pi=pi(0)+chi`, and the cocycle law makes `chi` a group homomorphism.  It
kills the `F_7` factor and hence every row-zero defect.  The ordered cut is
the necessary symmetry-breaking datum for this construction; equations
(7),(19) show how to retain it coherently instead of quotienting it away.

## 5. Quantitative inheritance from THM-2507

For the THM-2436 defects, `||d||_1<=18`.  By (4) and THM-2507's Galois
saturation, at least

```text
294*12=3,528
```

entries `Theta_(tau,a,c)(alpha)` are nonzero, each with

```text
|Theta_(tau,a,c)(alpha)|>=18^(-11).                            (22)
```

Thus the unnormalized pointwise bundle energy satisfies

```text
sum_(tau,a,c,alpha!=0)|Theta_(tau,a,c)(alpha)|^2
 >=3528*18^(-22).                                               (23)
```

Integrating (4) over the essential-parent locus and averaging over the `504`
component labels in any chosen measurable bundle trivialization loses no
extra cut factor: some fixed `(tau,a,c)` is nonzero on at least `7/12` of
that locus.  Hence the inherited parent floors remain `1/6` and `1/4`.
This is pointwise nonvanishing in a covariant finite bundle, not yet one
integrated physical current.  THM-2507's finite-exact
eleven-of-twelve sharpening applies to its normalized baseline cut and, by
(7), to every transported affine image with the transported cut.  It is not
asserted for all `42` cuts of each fixed normalized atlas defect.

## 6. Exact gain and remaining physical boundary

The proved chain is now

```text
essential punctured-stalk defect
  -> all 72 primitive 91-colours                         [THM-2506]
  -> finite nonaffine toothpick detection                [THM-2507]
  -> exact affine carry covariance of the whole bank
  -> lossless diagonal primitive-mode intertwiner        [THM-2508]. (24)
```

This closes the **static normalization/carry coherence** question left by
THM-2506/2507.  The carry is a permutation cocycle, not an information loss,
provided the cut torsor is retained.

It does **not**:

- turn one `R_(tau,a,c)` into a standard homomorphic `13`-root current;
- put the cut phase on the same ancestry sheet as THM-2471's source/arrival
  atoms or THM-2478's old deep probe;
- transfer this already-empty high-septimal packet to any of the live `165`
  low-septimal rows;
- prove a nonzero pairing with the THM-2334 relation-address action, a
  terminal-component phase, or owner-loop drift; or
- exclude a scalar row or prove LRC(14).

The next object should therefore be a **bilinear pairing of local systems**:
retain the cut-character bundle on one leg and construct its dual on the
owner/endpoint/deep ancestry cospan.  Another invariant scalar collapse or
another colour count cannot perform that transport, by (20)--(21).

## 7. Exact companion and audit

Run

```text
python3 04-computation/lrc14_affine_cut_bundle_covariance_thm2508_referee.py
python3 -O 04-computation/lrc14_affine_cut_bundle_covariance_thm2508_referee.py
```

The dependency-free exact referee works in `F_547`, containing primitive
seventh and thirteenth roots.  It verifies (14) coefficientwise on all
`471,744` basis rows, proves all `5,184` indexed kernel evaluations nonzero,
checks
all `5,184` primitive cut coefficients on the explicit THM-2506 two-row
hostile, verifies `864` zero-cut-character vanishings, and checks `596,232`
raw translation-covariance components.  Normal and optimized executions
reproduce

```text
05-knowledge/results/lrc14_affine_cut_bundle_covariance_thm2508_referee.out
```

byte-for-byte.  Independent audits rederived the full multiplier/translation
law (7), Fourier phase (9), signs in (14), nonvanishing proof (16), covariance
(19), regular cut torsor, and invariant-linear-scalar boundary.  A separate
companion, `lrc14_affine_cut_bundle_probe.py`, checks the regular `6,552`-entry
affine action, all `25,041,744` source-index covariance identities, `26,364`
DFT phase identities, the full factorization bank, the invariant-linear
collapse, and the quadratic-energy boundary, including the `kappa=2` hostile.
Its normal and optimized transcripts also match its stored output.  The scope
distinction between static mixed-bundle covariance and temporal physical
transplantation is explicit. **QED.**
