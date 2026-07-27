---
id: THM-2607
title: "Constant-six rail-boundary holonomy invoice"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On THM-2600's literal common-x q=0 carrier, 138 of the 162 positive
  constant-six middle rails have unit primitive Bockstein.  Their complete
  seven-clock choice polytope has 900 selectors; 891 have nonzero oriented
  rail boundary.  If n is the number of theta-one rails, the retained
  physical endpoints give c=w-v=7 on those n edges and zero otherwise, so
  the cyclic correction is 7n and exactly cancels THM-2542's constant-a
  class when a=-n mod 13.  The forward orientation realizes every marker
  a=6,...,12 on 57 displacement/marker/count lanes; the published THM-2600
  selector has sums 10 at s=6,11 and 7 at s=8, matching a=6,6,12.
  Each matched combined cochain is exact and has exactly thirteen vertex
  gauges.  Equivalently, the two rail edges are the C91 lifts 78 and 85
  of one C7 clock edge; the combined lifted loop closes exactly on the
  same marker lanes.  This is an abstract root/rail local-system intertwiner only:
  the natural chart rotation x->x+m/91 moves future digit six into the
  pairs {4,5},{2,3},{0,1},{11,12},{9,10},{7,8}, giving zero adjacent-chart
  overlap with the entire depth-five arrival/future toothpick.  Chronological
  transport does not repair this: every positive Koopman iterate kills the
  source C13 deck exactly and rebases C91 to an alternating C7 owner phase;
  the same future point has positive rail antecedents with c=0 and c=7.  No semantic
  vertical arrow, physical chart/deck identification, row exclusion, or
  LRC(14) conclusion is proved.
source: wild-holotopy-2026-07-28-constant-six-rail-boundary
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2599-rootwise-opposite-shift-paired-slice-law
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
related:
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2605-inverse-root-dipole-connection-and-mixed-square-invoice
script: 04-computation/lrc14_constant_six_rail_holonomy_thm2607.py
output: 05-knowledge/results/lrc14_constant_six_rail_holonomy_thm2607.out
script_sha256: 92e73022d65aa44974a286079ab9759a60a5de7e93f156b2dfe931f1527772dd
output_sha256: aacae8074f5707bec9c3e18b7de677517f5fa2864effec0069e9a81a635410b1
hash_basis: LF-normalized bytes
---

# THM-2607 -- constant-six rail-boundary holonomy invoice

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2602 shows that a root chosen separately at every chart contributes
only a telescoping vertex coboundary.  THM-2600 contains a different object:
each selected middle rail retains **two** physical root labels, the arrival
digit `v` and the rescaled deepest root `w=7t`.  Their oriented difference is
an edge cochain rather than a vertex selector.

That cochain has nonzero cyclic sum on the actual common-`x`, unit-Bockstein
carrier.  Abstractly it cancels seven of the twelve nonzero THM-2542 marker
classes.  The physical chart identification needed to interpret it as the
missing semantic two-cell is still absent, and the natural candidate for
that identification has exactly zero overlap.

## 1. The two constant-six rail endpoints

Work on THM-2600's canonical typed row, at its fixed target section

```text
q=0.                                                       (1)
```

The two middle rails have

```text
theta-zero: (v,t,w)=(6,12,6),
theta-one:  (v,t,w)=(6, 0,0),              w=7t mod 13.    (2)
```

Orient a rail from arrival to deepest root and define

```text
c=w-v in F_13.                                             (3)
```

Then

```text
c=0  on the theta-zero rail,
c=7  on the theta-one rail.                               (4)
```

Unlike a chartwise root choice, (3) retains both endpoints of one positive
THM-2584 rail cell.  The assertion that it is a **THM-2542 deck correction**
will remain conditional below; its status as a physical rail boundary inside
THM-2600 is unconditional.

## 2. The complete q=0 unit-edge choice polytope

The exact companion rebuilds only the `q=0` face of THM-2600's full
common-`x` bank.  It uses the same global primitive content

```text
g0=4,244,240                                                 (5)
```

and never reprimitivizes an edge.  Exactly

```text
138/162                                                     (6)
```

positive middle rails have unit septimal Bockstein at `q=0`.  For a fixed
collision displacement `s` and clock `ell`, write the admissible set as the
set of `t` values whose slice is both positive and unit.  The complete table
compresses to

| `s` | admissible `t` by clock `ell=0,...,6` | attainable `n=#{ell:t_ell=0}` |
|---:|:---|:---|
| `1,3,4,9,10,12` | `{0,12}` at every clock | `0,...,7` |
| `2,7` | `{12}` at every clock | `0` |
| `5` | `{0,12}` except `{12}` at `ell=5` | `0,...,6` |
| `6,11` | `{0}` at every clock | `7` |
| `8` | `{0}` at `ell=2`, `{0,12}` otherwise | `1,...,7` |

Thus the fixed `q=0` carrier has

```text
900 total seven-clock unit selectors,
891 selectors with n>0,
57 distinct (s,marker,n) count-level lanes.                (7)
```

For any one selector, (4) gives the cyclic correction

```text
C=sum_ell c_ell=7n mod 13.                                 (8)
```

THM-2542's constant marker `a` has class `7a`.  Hence the forward-oriented
rail pays its exact invoice precisely when

```text
C=-7a,
a=-n mod 13.                                               (9)
```

Because `1<=n<=7`, the forward orientation realizes exactly

```text
a in {6,7,8,9,10,11,12}.                                  (10)
```

The numbers of the 891 selectors paying the seven marker classes
`a=6,...,12` are respectively

```text
9, 49, 147, 245, 245, 147, 49.                            (11)
```

These counts are selectors, not masses, and no uniform lower bound is
inferred from them.

If the physical semantic arrow were legitimately oriented in the reverse
direction, `c` would be replaced by `-c`; then the same choice polytope would
realize marker classes `a=1,...,7`.  This is an algebraic orientation
boundary, not permission to reverse a directed semantic arrival.  No reverse
orientation is used in the forward theorem.

## 3. The published selector already has nonzero boundary

THM-2600 freezes one explicit selector by preferring `t=12` whenever it is a
unit.  Its theta-one counts are

```text
n_s=7 for s=6,11,
n_8=1,
n_s=0 otherwise.                                          (12)
```

Therefore

```text
C_s=10 for s=6,11,
C_8=7,
C_s=0 otherwise.                                          (13)
```

The nonzero values satisfy

```text
C_6=C_11=10=-7*6 mod 13,
C_8=7=-7*12 mod 13.                                       (14)
```

Thus the prescribed selector alone pays the abstract invoice on

```text
(s,a)=(6,6),(11,6),(8,12).                                (15)
```

For `s=6,11`, every edge has `a+c_ell=6+7=0`; cancellation is pointwise
around the seven-cycle.  The `s=8,a=12` lane is genuinely global:

```text
(a+c_ell)_(ell=0)^6=(12,12,6,12,12,12,12),
sum_ell(a+c_ell)=0 mod 13.                                 (16)
```

## 4. Exact Cech trivialization

Let the THM-2542 root-deck one-cochain be

```text
g_ell=a.                                                   (17)
```

For any matched selector in (9), the combined cochain is

```text
b_ell=g_ell+c_ell=a+c_ell,
sum_ell b_ell=7a+C=0.                                     (18)
```

On the cycle graph `C_7`, an `F_13`-valued one-cochain is a coboundary if
and only if its cyclic sum is zero.  Hence there is a vertex gauge `h` with

```text
h_ell-h_(ell-1)=-(a+c_ell),
a+c_ell+h_ell-h_(ell-1)=0.                                (19)
```

Choosing one base value determines the other six, and all thirteen base
values work.  Therefore every matched rail selector has exactly

```text
13 vertex gauges.                                         (20)
```

The companion checks all `891*13=11,583` gauges.  This proves more than the
necessary sum invoice: after granting the root/rail identification and the
orientation (3), the tensor product of the constant-`a` root local system
with the selected rail local system is abstractly trivial.

### 4.1 The `C_91` mapping-cylinder normal form

The same calculation has a canonical lifted-loop form.  Use the exact
sequence

```text
0 -> F_13 --iota--> Z/91Z --pi--> F_7 -> 0,

iota(r)=14r mod 91,                 pi(z)=z mod 7.           (20a)
```

Thus `iota(F_13)` is exactly the kernel of the clock projection `pi`.
The unique lift of one positively oriented clock edge whose kernel/root
correction is `c` is

```text
delta(c)=78+iota(c) mod 91,

pi(delta(c))=1,                    delta(c)=c mod 13.        (20b)
```

The two physical constant-six rails are therefore

```text
theta-zero: delta(0)=78,
theta-one:  delta(7)=85.                                  (20c)
```

For a selector with `n` theta-one edges,

```text
sum_ell delta(c_ell)=7n mod 91.                            (20d)
```

The constant marker contributes `iota(a)` on each of seven edges, hence

```text
sum_ell [delta(c_ell)+iota(a)]
 =7(n+a) mod 91.                                          (20e)
```

The lifted loop closes exactly when `a=-n mod 13`, which is (9).  Once its
quotient clock vertex is fixed, a lift has thirteen possible kernel-fibre
origins.  These are precisely the thirteen additive vertex gauges in (20),
after the same rail/deck identification has been granted.  This is a normal
form for the conditional abstract intertwiner, not an additional physical
identification.

It also explains the chronological failure before any interval geometry.
Multiplication by thirteen on `Z/91Z` has kernel `iota(F_13)`, and

```text
13 delta(c)=13 mod 91                  for every c in F_13,

13 iota(r)=0 mod 91.                                      (20f)
```

In particular it sends both `78` and `85` to the same edge `13`, whose
clock image is `-1 mod 7`.  Positive time reverses the quotient-clock
orientation and annihilates exactly the kernel component that paid the
marker invoice.  Section 5 shows that the corresponding physical overlap
also vanishes on the retained constant-six carrier.

## 5. The natural physical chart edge is exactly zero

The grant in Section 4 is not currently supplied by the LRC geometry.  The
most natural chart transport is the source-owner rotation

```text
x -> x+m/91,                  m=1,...,6.                   (21)
```

At THM-2600's response clock `R=13^6`, it changes the delayed phase by

```text
R*m/91=m*13^5/7=-m/7 mod 1.                               (22)
```

The future digit-six cell has length `1/13`.  Exact denominator-91 interval
arithmetic gives the possible new digit pairs

```text
m=1: {4,5},
m=2: {2,3},
m=3: {0,1},
m=4: {11,12},
m=5: {9,10},
m=6: {7,8}.                                                (23)
```

None contains digit six.  Hence every selected fixed-six slice is disjoint
from every nontrivial naturally rotated copy.

The obstruction is stronger than the fixed selector.  THM-2584's arrival
vertices are only

```text
{0,6,12}.                                                  (24)
```

Among (23), the only possible arrival/future matches are digit zero at
`m=3` and digit twelve at `m=4`.  But THM-2600 proves that delayed-word
digits zero and twelve are empty in every owner phase.  Therefore the
natural `C_91` rotation has zero adjacent-chart overlap with the **entire**
depth-five arrival-to-future diagonal toothpick, not only with the chosen
`q=0` unit selector.

This is the sharp physical stopping boundary.  A future edge must retain a
nontrivial digit carry and change/broaden the terminal word, or use a
different non-rotation correspondence.  THM-2542's BV mixing can visit the
seven positive rail supports at widely separated times, but that itinerary
is not (21) and does not identify the chart local system.

### 5.1 The chronological Koopman edge erases the root deck

THM-2604 and THM-2599 now supply genuine later same-root role events on one
orbit.  They do not transport the rail class.  Indeed, for every `N>=1`,

```text
T^N(x+q/13)=T^N(x),                                      (24a)

T^N(x+m/91)
 =T^N(x)+m 13^(N-1)/7
 =T^N(x)+(-1)^(N-1)m/7 mod 1.                            (24b)
```

Thus positive Koopman time kills the entire source `C_13` deck action and
rebases a source `C_91` chart translation to the coarse, parity-alternating
`C_7` owner phase.  Requiring

```text
floor(13x)=floor(13T^N x)=6                              (24c)
```

matches two root **values**; it is not an equivariant transport of the root
torsor.

There is a sharp rail-level hostile.  At `(s,ell)=(1,0)`, both the
theta-zero and theta-one `q=0` rail cells are positive open sets.  Call
positive open subintervals of them `J_0,J_7`, according as their correction
is `c=0,7`.  For every sufficiently large `N`, each expanding image is the
whole circle:

```text
T^N(J_0)=T^N(J_7)=R/Z.                                  (24d)
```

Consequently every future point `y` has antecedents `x_0 in J_0` and
`x_7 in J_7` with

```text
T^N x_0=T^N x_7=y,              c(x_0)=0, c(x_7)=7.      (24e)
```

So the rail correction is not a function of the chronological future.
The exact cylinder specification in THM-2599 can instead realize arbitrary
positive source and future rail cells, which makes the two corrections
independent rather than covariantly identified.  A lawful positive result
must therefore retain an inverse-branch/ancestry section carrying the rail
endpoint labels, not merely a later same-root occurrence.

## 6. Exact gain and scope

The proved and conditional parts are deliberately separated:

```text
PROVED physical data:
  one common-x positive/unit q=0 rail bank;
  both endpoint labels v,w on every chosen edge;
  the exact edge cochain c=w-v and its complete choice polytope;

PROVED abstract consequence:
  if c is the THM-2542 deck correction, every matched lane has
  trivial combined H^1 class and exactly thirteen gauges;

PROVED hostile:
  neither the natural C91 rotation nor the bare chronological Koopman edge
  can realize that identification;

NOT SUPPLIED:
  rail-root/deck-root intertwiner, chart-clock identification,
  semantic vertical arrow, or one naturally composable physical K bank. (25)
```

The unit Bockstein is a globally primitive aggregate on the same
pre-marginal carrier as the positive rail slice.  It does not select the
hidden positive Perron sheet, so no sheetwise charge is inferred.  The
collision displacement `s`, marker root `a`, target section `q`, rail roots
`v,w`, and two clock labels retain their separate types.

No typed row is excluded and LRC(14) remains open.

## 7. Exact companion

Run

```text
python 04-computation/lrc14_constant_six_rail_holonomy_thm2607.py
python -O 04-computation/lrc14_constant_six_rail_holonomy_thm2607.py
```

The companion reconstructs the `q=0` common-`x` bank directly from the
THM-2592/2600 exact interval machinery.  It checks all 4,752 inherited global
content divisions, all 138 unit edges, the complete 900-selector polytope,
all 891 nonzero invoices and 11,583 gauges, the published selector, both
orientation spectra, every denominator-91 digit intersection in (23), and
the exact C13-erasure/C91-to-C7 Koopman identities (24a)--(24b).
All decisions are exact and every check raises under optimized Python.

Normal and optimized executions byte-match the stored transcript after LF
normalization.

QED (candidate; independent audit pending).
