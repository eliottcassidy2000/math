---
id: THM-2630
title: "Old-wall affine clutching and successor-sector no-go"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the old THM-2587 speed-two wall, if t is the absolute deep digit,
  tau=t+epsilon is a chosen incident translated-gate edge, and
  theta=t-2v, then the THM-2586/2600 physical root and THM-2616 framed
  target label obey q'=h=v=7(tau-epsilon-theta) mod 13.  The THM-2586
  priority carriers have a positive middle wall in 80 of 84 cells, with
  four sharp pure-low cells.  THM-2600's 69 theta-zero and 15 theta-one
  q=0 carriers have positive middle wall in all 84 cells even after the
  complete common-x pullback; omitting theta fails exactly on the 15
  theta-one cells.  This does not identify tau with THM-2600's later probe
  r in Delta_r=d(2*13^5 x-r/13).  Splitting every later tooth into its
  left/right edge on THM-2616's raw q=h diagonal gives 3,792 signatures,
  of which 3,476 are nonfunctional; 158 of 162 rails and every one of the
  84 base cells retain multiple h values at fixed (theta,r,epsilon).
  Retaining the complete future-owner clock leaves 9,768 of 10,480 fully
  labelled signatures nonfunctional.  Exhausting all 156 affine bijections
  r=alpha*h+beta on the 888 visible clock/edge signatures leaves at least
  five h values everywhere.  In particular THM-2629's coefficient-optimal
  r=-h-1 graph has physical support histogram
  7^132,8^143,9^266,10^257,11^90 and is not a successor checksum.  Even the
  later half-bit and middle-wall state leave seven successor sectors.  The
  missing datum is a physical successor-sector/carry sidecar, not another
  binary wall, owner-clock label, or affine coefficient graph.  No row
  exclusion or LRC(14) conclusion follows.
source: deep-energy-audit-2026-07-28-wall-successor-square
depends_on:
  - THM-2586-depth-five-arrival-to-future-root-diagonal
  - THM-2587-deep-root-coshift-incidence-wall-and-theta-selector-no-go
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
  - THM-2613-canonical-root-diagonal-opposite-shift-section
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
related:
  - THM-2614-punctured-target-root-cosupport-and-principal-deck-no-go
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2629-fixed-deep-affine-graph-spectrum-and-puncture-cancellation-boundary
script: 04-computation/lrc14_old_wall_successor_sector_thm2630.py
output: 05-knowledge/results/lrc14_old_wall_successor_sector_thm2630.out
script_sha256: 06f4eeddb498e03b1cfcd94af3d9ffc57a62144e79ec4649329b888d8810f701
output_sha256: cafeefc7d993cc12da0158464110e9b72d3201259e6d437803a1cbaff8a455b0
hash_basis: LF-normalized bytes
---

# THM-2630 -- the old wall clutches; the later probe does not

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

Four nearby constructions print residues in `F_13`:

```text
THM-2587: old absolute deep digit t and translated-gate label tau;
THM-2600: later translated deep-probe label r;
THM-2613: physical inverse-root digit h;
THM-2616: framed next-target label q'=h.                         (1)
```

The first, third, and fourth labels do form an exact affine clutching square
on the old wall.  The second label lives five base-thirteen scales later.
One left/right wall bit does not bridge that time change.  The distinction is
not merely formal: an exact refinement of the proved common carrier is
nonfunctional in every base cell.

## 1. The old wall gives an affine clutching formula

Use THM-2584/2587's route-five physical coordinate `x`.  Away from null cell
endpoints, write

```text
2x=(t+z)/13 mod 1,

t=floor(26x) mod 13,                   0<=z<1.             (2)
```

For the old translated danger gate with label `tau`, THM-2587 proves

```text
tau=t+epsilon,

epsilon=0 permitted iff z<13/14,
epsilon=1 permitted iff z>1/14.                              (3)
```

Thus low wall points force `epsilon=0`, high wall points force
`epsilon=1`, and the open middle wall permits both.  At the strict seams the
rule is unambiguous:

```text
z=1/14:   epsilon=0 only,
z=13/14:  epsilon=1 only.                                  (4)
```

Let `v` be the arrival root and retain the physical rail bit

```text
theta=t-2v in {0,1}.                                       (5)
```

Since `2^(-1)=7 mod 13`, equations (3)--(5) give

```text
v=7(t-theta)=7(tau-epsilon-theta).                         (6)
```

THM-2586's later physical root equals `v` on its frozen carrier.  On
THM-2600's constant-six carrier, the delayed physical digit is also
`h=v=6`.  Finally, for the canonical ordinary role `q_1=14=1 mod 13`,
THM-2613 gives local shift `s=-h`, while THM-2585's convention is `q'=-s`.
Consequently the old-wall square commutes as

```text
q'=h=v=7(tau-epsilon-theta) mod 13.                        (7)
```

This is affine and equivariant in the fixed physical gauge.  Normalize
`tau_bar=7tau`.  Under diagonal translation `v -> v+a`, one has

```text
t -> t+2a,       tau_bar -> tau_bar+a,       h -> h+a,

epsilon,theta fixed.                                      (8)
```

Hence (7) is an equivariant clutching formula, not numerical coincidence.
On a middle wall, the overlap-edge bit `epsilon` is necessary and sufficient
to invert the two-valued incidence relation.

## 2. The middle ambiguity is physical on the frozen carriers

The THM-2586 priority selector uses `(v,t)=(0,0)` except in its three
fallback cells, where it uses `(6,12)`.  Refining those exact positive
carriers by the old wall gives the support patterns

```text
(low,middle,high) positive:

(1,1,0): 44 cells,
(0,1,0): 36 cells,
(1,0,0):  4 cells.                                        (9)
```

Thus the overlap-edge bit is genuinely needed in `80/84` frozen carriers.
The four sharp pure-low exceptions are

```text
(s,ell)=(7,0),(7,1),(7,2),(7,3).                          (10)
```

On those four cells `epsilon=0` is already forced.

THM-2600 uses a different selector: arrival `v=6`, with

```text
(t,theta)=(12,0) on 69 cells,
(t,theta)=( 0,1) on 15 cells.                             (11)
```

Its exact old-wall patterns are

```text
theta=0: 39 middle+high, 30 middle-only;
theta=1:  8 low+middle,   7 middle-only.                  (12)
```

In particular, the middle ambiguity occurs in every one of the `84` rail
cells.  This remains true on a literal full common-`x` pullback, not merely
before unrelated factors are inserted.  Refine THM-2600's selected `q=0`,
future-digit-six carrier simultaneously by

```text
the old middle wall,
the full present F_(ell5,0) packet,
the actual later Delta_r probe,
the delayed Q^+_(ell5,6) word.                            (13)
```

Every base cell has positive mass.  The number of positive `(ell5,r)` fine
entries per cell has exact histogram

```text
12:12 cells,   24:23 cells,   36:37 cells,   48:12 cells. (14)
```

Equation (14) licenses the pointwise old-wall formula on the complete
carrier.  It does **not** say that a Bockstein remains a unit after the wall
restriction, and no hidden Perron sheet is charged by aggregate unitness.

The theta correction is sharp.  The tempting theta-zero formula

```text
h=7(tau-epsilon)                                          (15)
```

works on the `69` theta-zero cells and fails on exactly all `15` theta-one
cells, by the constant error `7`.  Equation (7) repairs all `84`.

## 3. THM-2600's later probe is not the old gate

The old label in (3) belongs to

```text
d_1(2x-tau/13).                                           (16)
```

THM-2600's displayed deep-probe label `r` instead belongs to

```text
Delta_r(x)=d_1(2*13^5 x-r/13).                            (17)
```

The factor `13^5` is load-bearing.  Put

```text
y={13^5 x},             j=floor(13y)=floor(13^6x) mod 13,

a=floor(13{2y}),                    z={13{2y}}.            (18)
```

Splitting the later translated tooth in (17) at its centre gives another
THM-2587-shaped local wall:

```text
r=a+epsilon,

epsilon=0: right half of the tooth,
epsilon=1: left half of the tooth.                         (19)
```

But THM-2616's physical successor digit is one digit deeper.  Write

```text
u={13y}={13^6x},             kappa=floor(2u),

z={2u},              h=floor(13u).                        (20)
```

The exact predecessor-carry identity is

```text
a=2j+kappa mod 13,

r=2j+kappa+epsilon mod 13.                                (20a)
```

Equivalently,

```text
a=floor(13*(2*13^5x)) mod 13
 =floor(2*13^6x) mod 13
 =2 floor(13^6x)+floor(2{13^6x}) mod 13.                  (20b)
```

Thus the tempting successor formula `r=2h+kappa+epsilon` is false: `j` is
the predecessor/carry residue, whereas `h` is the next digit of `u`.

Then the exact recovery formula is

```text
h=floor(13(kappa+z)/2).                                   (21)
```

The later edge bit in (19) recovers `a=r-epsilon`.  Supplying `kappa` then
recovers the predecessor carry `j=7(a-kappa)`, but still does not recover the
successor-sector position of `z` in (21).  The old rail bit `theta` is at the
route-five wall and cannot replace the fresh later-scale branch `kappa`.

## 4. A strict two-point hostile and the seven-sector minimum

The failure is visible away from every seam.  Take local phases

```text
y_1=11/260,                 y_2=3/52.                     (22)
```

Both have

```text
r=1,       epsilon=0,       kappa=1,

1/14<z<13/14,                                                (23)
```

with exact middle coordinates

```text
z(y_1)=1/10,                 z(y_2)=1/2.                   (24)
```

Nevertheless

```text
h(y_1)=7,                    h(y_2)=9.                     (25)
```

Thus even `(r,epsilon,kappa,middle-wall-state)` is not a successor map.
The failure is information-theoretically sharp.  At fixed
`(r,epsilon)=(1,0)` and fixed open middle wall, the two values of `kappa`
support exactly

```text
kappa=0: h in {0,1,2,3,4,5,6},
kappa=1: h in {6,7,8,9,10,11,12}.                         (26)
```

Every displayed sector has a strict rational interior witness.  Therefore,
after `kappa` is supplied, a finite deterministic sidecar recovering `h`
still needs at least seven states on either half.  The sufficient minimal
refinement is the seven-sector partition of `z` in (21); after the two halves
are glued, it is equivalently the full successor digit `h`.

## 5. Exact no-go on THM-2616's common carrier

The ambient hostile survives the canonical packet with nearly maximal
support.  Start from THM-2616's raw positive `q=h` diagonal **before** the
seven-clock unit test.  For every one of its `162` rails, every future clock,
every `h`, and every nonzero later probe `r`, split (17) into the two halves
in (19) before delayed integration.  All `61,248` nonempty fine tooth
partition candidates over nonempty pre-probe overlaps split exactly into left
plus right; a candidate half may itself be empty.

For a base rail cell retain the support signature

```text
(s,ell4,theta,r,epsilon) -> {positive h values}.           (27)
```

After the future clock and hidden positive sheets are united, the exact
census is

```text
present signatures:       3,792,
nonfunctional signatures: 3,476,
functional singletons:      316,
ambiguous rails:         158/162,
base cells with ambiguity: 84/84.                         (28)
```

Every nonfunctional signature has support size ten or eleven.  The complete
ambiguous histogram by old rail bit is

```text
theta=0: size 10 -> 484,     size 11 -> 1,254;
theta=1: size 10 -> 286,     size 11 -> 1,452.             (29)
```

The `316` singleton signatures do not assemble a map: every base cell also
contains a nonfunctional signature.  Equations (28)--(29) are support of the
lawful raw diagonal carrier.  They do not assign an aggregate unit
Bockstein to one half-tooth sheet.

Retaining the future-owner clock `ell_5` rather than taking the union in
(27) gives the strictly finer signatures

```text
(s,ell_4,theta,ell_5,r,epsilon) -> {positive h values}.   (29a)
```

Their exact census is

```text
support size          1     7    8     9     10    11
signatures          712  1452  484  2904   2948  1980,

total=10,480,       nonfunctional=9,768,       ambiguous base cells=84/84.
                                                               (29b)
```

Thus the already present seven-valued owner clock is not the missing
seven-sector coordinate.  Some fully labelled signatures still carry all
eleven allowed future digits.

There is also a direct comparison with THM-2629's affine optimum.  For each
of the `888` visible tuples `(s,ell_4,theta,ell_5,epsilon)` and every affine
bijection

```text
r=alpha h+beta,                    alpha !=0,             (29c)
```

retain the `h` for which the corresponding half-tooth coefficient is
positive, deleting the unique `h` at which `r=0`.  Exhausting all `156`
choices in (29c) gives at least five surviving `h` values in **every**
visible tuple.  If each affine map is scored by its largest fibre, the exact
histogram is

```text
largest fibre          9    10    11
affine bijections     74    74     8.                    (29d)
```

For the unique coefficient-level unit optimum from THM-2629,
`r=-h-1`, the complete physical fibre histogram is

```text
size                   7    8    9    10   11
visible tuples        132  143  266   257   90.           (29e)
```

Consequently the opposite graph is not a hidden chronological checksum on
this carrier.  It belongs to the worst `11`-fibre class in (29d), even
though it is uniquely best for THM-2629's nonlinear seven-clock unit score.
This cleanly separates spectral cancellation from physical successor
selection.

## 6. Holotopy square and exact stopping boundary

The completed and broken squares are

```text
old absolute t -- chosen old wall edge epsilon --> old tau
       |                                            |
       | h=7(t-theta)                               | equation (7)
       v                                            v
physical h ---------------------------------------> q'=h


later absolute a -- chosen later wall edge epsilon --> THM-2600 r
       |                                                |
       | successor sector z,kappa missing               | no map
       v                                                v
physical h -------------------------------------------> q'=h.       (30)
```

The first square is an affine clutching in one fixed gauge.  The second is
not a square until a physical successor-sector sidecar is carried on the
same ancestry.  Merely adding the absent future support sheets, taking a
safe/danger cospan, or selecting a fixed affine graph at coefficient level
does not provide that sidecar or its adjacent-clock deck action.

The exact connection ledger is:

| item | content |
|---|---|
| old positive control | `q'=h=7(tau-epsilon-theta)` |
| sharp old ambiguity | one edge bit on the open middle wall |
| first failed implication | old `tau` and later `r` have the same residue alphabet, but different speeds/times |
| later hostile | fixed `(theta,r,epsilon)` supports up to eleven `h` values on the canonical diagonal |
| failed refinements | later edge bit, future-owner clock, and every affine bijection `r=alpha h+beta` |
| minimal local repair | fresh `kappa` plus the seven-sector position of `z`, equivalently `h` |
| still missing after repair | one principal adjacent-clock action and semantic endpoint/current transport |

The result is confined to the canonical typed row.  It does not select one
unit hidden sheet, retain the old relation residue, identify the two owner
clocks, exclude a scalar row, or prove LRC(14).  The ledger remains `165`.

## 7. Exact companion

Run

```bash
python 04-computation/lrc14_old_wall_successor_sector_thm2630.py
python -O 04-computation/lrc14_old_wall_successor_sector_thm2630.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_old_wall_successor_sector_thm2630.out.
```

The companion rebuilds the old wall independently from the exact packet,
checks the strict seams, all affine gauges, both frozen-carrier censuses, and
the full common-`x` middle restriction.  Four deterministic shards then
check all `61,248` candidate additive partitions over nonempty pre-probe
overlaps on the raw `q=h` carrier and reproduce (28)--(29b).  They retain the
future clock, exhaust all `156` affine bijections, and check (29c)--(29e).
A separate rational construction verifies (22)--(26).  Every logical
decision is an explicit optimized-mode guard.

QED (candidate; independent hostile audit pending).
