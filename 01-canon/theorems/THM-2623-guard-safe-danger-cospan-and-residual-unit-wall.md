---
id: THM-2623
title: "Guard-safe/danger cospan, half-tooth repair, and carry-loss wall"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Complementing only the delayed guard factor gives two disjoint Boolean
  endpoint sectors whose physical future-digit supports are respectively
  {1,...,11} and {0,1,11,12}; their labelled union is all F13 and retains
  global primitive content 26.  The unsplit q=h Bockstein atlas is full on
  76/84 base cells, with eight explicit singleton unit failures.  Splitting
  the actual high-speed deep tooth into its left/right halves preserves
  content 26 and repairs exactly four failures, leaving only
  ((2,2),11), ((6,5),2), ((7,2),11), ((11,5),2).  It does not create a
  private deep-root row: every nonzero refined row has support size 11 or 12.
  A further future half-digit split still gives zero private unit roots on
  all 84 cells because it forgets floor(Rx) mod 13.  Thus the cospan fills
  support and partially repairs unitness, but supplies neither a positive
  root decoder nor a principal C13 transition, row exclusion, or LRC(14).
source: carry-transition-cell-2026-07-28-guard-cospan
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
related:
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2571-deep-colour-cayley-filling-bockstein-and-norm-curvature-split
  - THM-2587-deep-root-coshift-incidence-wall-and-theta-selector-no-go
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2614-punctured-target-root-cosupport-and-principal-deck-no-go
  - THM-2624-two-clock-signed-target-root-tomography-and-carrier-holotopy-boundary
  - THM-2629-fixed-deep-affine-graph-spectrum-and-puncture-cancellation-boundary
  - THM-2630-high-speed-half-tooth-ambiguity-and-carry-coordinate-no-go
script: 04-computation/lrc14_guard_cospan_half_tooth_thm2623.py
output: 05-knowledge/results/lrc14_guard_cospan_half_tooth_thm2623.out
script_sha256: edac7583aae8f20d9de105ab50bd497b74437a01fafdfc4ce78b1f6ad38405b6
output_sha256: bd02c6ee318a128ad40b401d8d5970f9e8ba99df9e416af1a598862697e3450b
secondary_script: 04-computation/lrc14_successor_halfcell_carry_no_go_thm2623.py
secondary_output: 05-knowledge/results/lrc14_successor_halfcell_carry_no_go_thm2623.out
secondary_script_sha256: 39fbff2a846c30e4d93dd881b4789a1af9175260d852b7e613e43dce62d43c6a
secondary_output_sha256: a5e84fccb64dabc28cc4a611501e384dc97ea57a11c3fb11fd8e11af25fb1a66
hash_basis: LF-normalized bytes
---

# THM-2623 -- the guard complement fills the digit puncture, not the carry torsor

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2616 found a positive cross-time diagonal `q=h`, but its delayed word
has no physical digits `0,12`.  The cheapest lawful completion is not an
invented sheet: complement the one delayed guard factor and retain both
Boolean sectors with their labels.  This fills the future support exactly.
It also exposes two sharper boundaries.  First, the nonlinear unit test does
not become full merely because the support does.  Second, even the natural
left/right deep-tooth refinement forgets the dilation carry needed to make a
private root response.

## 1. The disjoint delayed guard cospan

Use THM-2616's canonical typed row, `162` middle rails, present target section
`q`, present owner clock `ell_5`, nonzero deep probe `r`, and delayed physical
digit `h`.  Write `Q^s_(ell,h)` for its original delayed endpoint sector.  It
contains

```text
D_(c2), all five ordinary safeties, c3-safety, c1-clock safety,
and guard safety.                                             (1)
```

Define `Q^d_(ell,h)` by replacing only guard safety in (1) by guard danger.
The two sectors are disjoint, and their pointwise sum before the unit test is
the lawful word with the guard factor omitted:

```text
Q^s + Q^d = Q^(omit H).                                     (2)
```

No endpoint factor other than the guard is changed.  In particular, this is
a physical Boolean cospan, not a formal coefficient completion.

The exact delayed supports are

```text
supp_h Q^s = {1,2,3,4,5,6,7,8,9,10,11},
supp_h Q^d = {0,1,11,12},
supp_h Q^s union supp_h Q^d = F_13.                        (3)
```

The safe sector is the THM-2616 carrier.  Across all digits its positive
pre-route census and content are

```text
safe:       649,968,       gcd 26;
danger:     244,992,       gcd 86,814=26*3,339;
labelled joint global gcd: 26.                              (4)
```

Thus the physical route-two content remains `13*26=338`.  One must retain
the sector label: unitness in
`F_13[z]/(Phi_7)` is nonlinear, so summing the two coefficient tables before
testing a slice need not preserve either sector's unit.

## 2. The full-support cospan still has eight unit holes

On the numerical cross-time diagonal `q=h`, normalize both sectors by the one
labelled global content `26`.  The safe and danger sectors have respectively

```text
unit slices:          1,483 and 533;
positive rail/q pairs: 1,703 and 613.                      (5)
```

Every nonzero fixed `(rail,q,ell_5,h=q)` row in either sector is positive on
all twelve nonzero deep probes:

```text
safe:   4,156 rows of support size 12;
danger: 1,526 rows of support size 12.                     (6)
```

Consequently no nonnegative combination of the unrefined rows can be a
nonzero singleton root response: every selected positive summand brings all
twelve roots.

For each base cell `(s,ell_4)`, allow either middle rail and either guard
sector.  The unit-section set has size thirteen on `76` cells and size twelve
on `8`.  The exact missing pairs are

```text
((2,2),11), ((6,5),2),  ((6,6),10), ((7,2),11),
((8,2),12),((11,2),12),((11,5),2),  ((11,6),10).           (7)
```

Thus filling the two physical digit punctures is strictly weaker than filling
the unit atlas.

## 3. The actual high-speed half-tooth bit repairs exactly four holes

For the actual deep probe

```text
Delta_r(x)=1_(||c3*x-r/13||<1/14),
c3=2*13^5,                                                (8)
```

split its strict danger tooth at its centre:

```text
Delta_(r,L)=1_(-1/14<c3*x-r/13<0),
Delta_(r,R)=1_(0<c3*x-r/13<1/14).                         (9)
```

The equality is a.e. and all endpoints have measure zero.  This is a
THM-2587-style incidence-edge bit on the **actual high-speed probe**.  It is
not THM-2587's older low-scale `(tau,theta)` wall, and no identification of
their labels is made.

Split (before primitive reduction) every safe/danger fine coefficient by
(9).  All `353,808` sector/tooth partition identities hold exactly.  The
refined full positive census is

```text
1,661,176,              global gcd 26.                    (10)
```

The four labelled unit counts are

```text
                 L       R
guard safe     1,442   1,450
guard danger     547     530,                              (11)
```

while both safe halves retain `1,703` positive rail/q pairs and both danger
halves retain `613`.

The refined atlas is full on `80/84` base cells.  The bit genuinely repairs
four entries of (7), rather than merely relabelling them.  Its exact residual
wall is

```text
((2,2),11), ((6,5),2), ((7,2),11), ((11,5),2).            (12)
```

The common unit sections are

```text
{0,1,3,4,5,6,7,8,9,10,12}.                               (13)
```

Nevertheless the refined response rows are still maximally nonprivate:

```text
safe L,R:   size 11 on 3,800 rows; size 12 on 356 rows each;
danger L:   size 11 on 1,150 rows; size 12 on 376 rows;
danger R:   size 11 on 1,130 rows; size 12 on 396 rows.    (14)
```

In particular, adding the edge bit does not evade the positive-decoder
obstruction.

## 4. Fixed affine-graph cross-check

As a comparison with THM-2629, externally insert the physical factor

```text
r(h)=-h-1 mod 13                                         (15)
```

before summing over roots.  This is a lawful conditioned subcarrier, not a
decoder extracted from (6) or (14).  On the guard/half-tooth cospan its unit
section set has size twelve on `80` cells and eleven on the four cells in
(12).  The new endpoint `h=0` is unit on all `84` cells, necessarily through
the guard-danger sector and `r=12`.  The endpoint `h=12` is zero on all `84`
cells because (15) gives the structurally absent probe `r=0`.  Apart from
that universal puncture, the only failures are exactly (12).

This is a positive punctured-graph control.  Since choosing (15) is external,
it does not contradict the no-decoder statement.

## 5. Why the tempting successor half-digit is not the missing sidecar

It is tempting to split the future digit further by

```text
h=floor(13{R*x}),
kappa=floor(26{R*x})-2h in {0,1},                         (16)
```

and infer a singleton law such as `r=2h+kappa` or
`r=2h+kappa+1`.  That inference is false on the actual speeds.  Since

```text
R=13^6,          c3=2*13^5=2R/13,                        (17)
```

write `R*x=n+y`, with `n=floor(R*x)` and `y={R*x}`.  The
absolute high-speed deep digit is

```text
a=floor(13{c3*x})
  =2n+floor(2y) mod 13.                                   (18)
```

The left/right tooth then has `r=a+1` or `r=a`, respectively.  Formula (16)
retains information inside `y` but forgets

```text
c=n mod 13.                                               (19)
```

This is the missing affine carry coordinate.

The second exact companion inserts all four `(L/R,kappa)` lanes into the
entire labelled cospan.  Its global content remains `26`, and it has
`2,986,852` positive fine entries.  Yet every nonempty clock row and every
nonempty seven-clock slice still has support size `11` or `12`.  There are
`20,778` direct failures of the false singleton formula, and the union of all
four lanes has

```text
number of private unit roots = 0             on all 84 cells. (20)
```

Thus even a four-chart cover supplies no positive private-root matching.
The minimal corrected refinement is to retain (19), or an equivalent
same-packet affine carry sidecar, before the tooth split.  THM-2555 proves the
general carry mechanism, but not positivity of (19) on this selected
THM-2616 carrier.  THM-2571 retains old/future digits and an integral
Bockstein on the same typed row, not `floor(Rx) mod 13` on this cospan.
THM-2579 is an integral torsor statement.  None supplies (19) here.

THM-2630 independently isolates the same high-speed carry law on the original
safe carrier; it is related evidence, not a dependency while under audit.

## 6. Holotopy interpretation and exact boundary

The connection contract is now exact:

| item | retained or lost content |
|---|---|
| source | two disjoint physical delayed guard sectors on one THM-2616 carrier |
| target | present absolute section `q` and delayed physical digit `h` |
| map | the cross-time diagonal `q=h`, refined by guard and half-tooth labels |
| preserved | all thirteen future supports, content `26`, common physical `x`, both clocks |
| repaired | four of eight nonlinear unit holes |
| still lost | four unit cells, private deep root, `floor(Rx) mod 13`, inter-clock transition |
| needed sidecar | the carry (19) on the same positive packet, plus a lawful clock/action gluing |

Topologically, the guard complement closes the support puncture, and the
half-tooth bit refines the cover.  It does not trivialize the affine carry
torsor on overlaps.  The surviving obstruction is therefore not another
missing Boolean endpoint; it is transition data.

No coefficient-level unit identifies a hidden Perron sheet.  No sector or
chart here is a semantic THM-2305 terminal endpoint.  No adjacent principal
`C_13` action, owner/root provenance, row exclusion, or LRC(14) proof follows.
The scalar ledger remains `165`.

## 7. Exact companions

Run

```bash
python3 04-computation/lrc14_guard_cospan_half_tooth_thm2623.py
python3 -O 04-computation/lrc14_guard_cospan_half_tooth_thm2623.py
python3 04-computation/lrc14_successor_halfcell_carry_no_go_thm2623.py
python3 -O 04-computation/lrc14_successor_halfcell_carry_no_go_thm2623.py
```

Normal and optimized executions must byte-match their stored transcripts.
The primary companion reconstructs all four THM-2616 rail shards, verifies
the safe/danger Boolean partition and all high-speed half-tooth partitions,
recomputes both global contents, all unit atlases, all support histograms, and
the fixed affine-graph control.  The secondary companion independently
refines every delayed digit into its two physical half-cells, checks the four
lane spectra, exhibits direct failures of the false singleton formula, and
proves the all-`84` no-private census.  Every logical check is an explicit
optimized-mode guard.

QED (candidate; independent hostile audit pending).
