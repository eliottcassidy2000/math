---
id: THM-2794
title: "Rail-eight fixed-configuration all-component semantic-atom no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the fixed THM-2672 configuration
  (rail,sector,edge,kappa,h)=(8,0,1,1,6), source carry 12 omits label one
  and has positive twelve-label facet mass exactly on clocks 1,2,3.  The
  matching THM-2749 fully marked rail-eight source carrier is also positive
  exactly on those clocks.  Nevertheless the facet base retains
  F_(ell,7), whereas the marked base retains its complement, so their strict
  open pre-prefix bases are disjoint on every clock.  Exhaustive
  delayed-prefix refinement gives
  165291,356973,356973 positive facet components, respectively, and zero
  positive overlap for all 879237 components.  There are 503 closure-only
  contacts, all oriented from a facet right endpoint to a marked-carrier
  left endpoint.  This excludes one carry in one fixed configuration, not
  another configuration, the full rail-eight union, a row, or LRC(14).
source: lrc-a12-carry-bridge/full-rail8-component-intersection-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2785-a12-minuscule-carry-fourier-character-and-common-atom-boundary
related:
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
script: 04-computation/lrc14_rail8_fixed_configuration_all_component_no_go_thm2794.py
output: 05-knowledge/results/lrc14_rail8_fixed_configuration_all_component_no_go_thm2794.out
script_sha256: fe836495e56e3fc731b3c6f59dea943ff6d1fda68bad0dcb7bd38863e7c99c31
output_sha256: d980fac00326a613cd867bfaf3600a4d5560226f88125ed666838d61568e3dae
hash_basis: LF-normalized bytes
---

# THM-2794 -- one full fixed rail-eight facet misses the semantic carrier

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2785 proves that one selected component of the rail-eight slope-seven
facet misses THM-2749's semantic carrier.  That could have been an artefact
of the lexicographic section.  It is not.  On the same source carry and
fixed configuration, the two positive open bases are already disjoint
before the delayed word is imposed.  Consequently every one of the
`879,237` delayed components misses.

The result is a complete fixed-configuration no-go, not a global rail-eight
no-go.  Its sharp boundary is visible: the two closures touch at `503`
delayed seams.  Changing a half-edge, configuration, or strictness convention
can therefore change the answer.

## 1. Exact fixed-configuration universe

Use

```text
p=13,                    R=13^6,                    S=13^5.           (1)
```

Fix the THM-2672 configuration

```text
e=(rail,sector,edge,kappa,h)=(8,0,1,1,6),                             (2)
```

whose rail metadata is `(1,4,12)`.  Fix source carry

```text
c=12.                                                                (3)
```

Its private root is `12`, its unique missing target label is

```text
m(12)=2(6-12)=1                         in F_13,                      (4)
```

and its active target labels are

```text
A=F_13\{1}={0,2,3,...,12}.                                           (5)
```

For each delayed clock `ell in F_7`, define the **fixed-facet base**
`B_ell` by intersecting, in source coordinates:

1. the twelve rail-eight supports pulled back by the slope-seven
   translations `7 delta/R`, `delta in A`;
2. the twelve matching relative-present supports at clock `ell`; and
3. the edge-one private half with root `12`.

No delayed carry prefix has yet been applied.  Since `-h=7 mod 13`,
the label-zero anchor in item 2 gives

```text
B_ell subset F_(ell,7).                                               (6)
```

Thus `B_ell` retains the rail, sector, edge, `kappa`, `h`, source
carry/root typing, and all twelve target labels, but not yet the
predecessor-carry word cylinder.

Let `D_ell` be the strict open subset obtained from `B_ell` by imposing the
clock-`ell`, sector-zero, `h=6`, `kappa=1`, carry-`12` delayed prefix.
Its connected open interval components are the delayed components counted
below.

Independently, let `M_ell` be THM-2749's fully marked rail-eight **source
base** before its delayed prefix.  It retains:

```text
source and pulled-back target rail profiles with equal weight;
both relative-present complements;
the source/root-12 and pulled-back target root halves;
S_(ell,0,4)=E3 intersect F_(ell,0,4);
the pulled-back target copy of that semantic section.
```

The THM-2749 relative-overlap construction instead contains

```text
M_ell subset F_(ell,7)^c.                                             (7)
```

The THM-2749 source coefficient then applies its semantic delayed prefix to
`M_ell`.  The comparison is clock-matched and source-coordinate matched,
but the two banks have opposite anchor-present polarity.  It is not a
comparison of unrelated coefficient marginals.

## 2. Both banks are positive on the same three clocks

Exact delayed integration of the fixed facet gives

```text
(mu_D(ell))_(ell=0)^6
 =(
 0,
 1509311690628117483066960,
 3259498009127699308489520,
 3259498009127699308489520,
 0,0,0).                                                            (8)
```

The sum is

```text
8028307708883516100046000,                                           (9)
```

the carry-`12` positive numerator already present in THM-2785.

The matching fully marked THM-2749 source vector is

```text
(mu_M(ell))_(ell=0)^6
 =(
 0,
 339633525654239542165440,
 750593782703678965571520,
 719200126392878704654080,
 0,0,0).                                                           (10)
```

Thus both exact positive-clock supports are

```text
{1,2,3}.                                                            (11)
```

This is the load-bearing positive control.  The no-go below is not caused
by putting one side on an inactive clock.

## 3. The strict bases are already disjoint

The exact positive weighted-piece counts of `B_ell` and `M_ell` are

| clock `ell` | pieces of `B_ell` | pieces of `M_ell` |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 239 | 239 |
| 2 | 516 | 526 |
| 3 | 516 | 504 |
| 4 | 0 | 0 |
| 5 | 0 | 0 |
| 6 | 0 | 0 |

The literal polarity `(6)--(7)` already gives

```text
B_ell intersect M_ell=empty                       for every ell.     (12)
```

Here `empty` means the strict open intersection used by the physical
carriers.  Their closures do meet.  On clocks `1,2,3`, respectively, there
are

```text
189, 416, 406                                                       (13)
```

shared base endpoints.  Every one has the same orientation:

```text
right endpoint of a fixed-facet base interval
  = left endpoint of a marked semantic interval.                    (14)
```

There is no contact of the opposite orientation.  Hence `(12)` is a
strict-seam separation rather than a large metric gap.

Because `D_ell` is obtained by intersecting `B_ell` with another predicate,
`(12)` already proves

```text
D_ell intersect M_ell=empty.                                        (15)
```

The delayed component enumeration below is an exact independent refinement
and a count of how much positive facet structure is excluded.

## 4. Exhaustive delayed-component refinement

Write the delayed prefix on the coordinate `{Rx}` as a union of strict
integer intervals inside period `13T`, with `T` the canonical interval
denominator.  For every weighted interval `(a,b)` of `B_ell`, the companion
intersects

```text
R(a,b)
```

with every integer translate of the carry-`12` delayed prefix.  All
resulting intervals are sorted and checked to be pairwise nonoverlapping.
Their exact counts are

| clock `ell` | positive delayed components | positive overlap with `M_ell` | closure seams |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 1 | 165,291 | 0 | 94 |
| 2 | 356,973 | 0 | 207 |
| 3 | 356,973 | 0 | 202 |
| 4 | 0 | 0 | 0 |
| 5 | 0 | 0 | 0 |
| 6 | 0 | 0 | 0 |

Therefore

```text
total positive delayed components =879237,
total positive overlap             =0,
total closure-only seams           =503.                            (16)
```

The `503` contacts again all have the orientation `(14)`.  Since every set
is strict open, a closure contact contributes neither a component
intersection nor positive coefficient mass.

Equations `(12)` and `(16)` are two exact verification routes:

```text
base route:       intersect first, then observe the empty base;
component route:  apply every delayed prefix interval, enumerate all
                  components, and sweep them against the marked carrier. (17)
```

Both give the same zero conclusion.

## 5. Consequence and precise stopping boundary

No clock-matched component in the carry-`12` facet of the single fixed
configuration `(2)` is a common atom with THM-2749's fully marked
rail-eight source carrier.  In particular, replacing THM-2785's
lexicographically first component by another delayed component in the same
fixed facet cannot repair the attachment.

This sharpens the next move:

```text
a successful attachment must leave at least one of
  fixed configuration e,
  source carry 12,
  twelve-label facet A,
  exact strict seam orientation,
  or THM-2749 semantic source base.                                 (18)
```

A configuration/edge switch can help only if it changes or removes the
anchor polarity `(6)--(7)`.  The closure contacts show where such a repair
could open, but `(13)--(14)` do not prove that any lawful switch does open a
positive overlap.

## 6. Compatibility with THM-2782 and THM-2785

There is no contradiction with either theorem.

- THM-2785's all-twelve carry colours are ordinary Fourier coefficients of
  the fixed-facet component indicators.  Fourier integration retains their
  character ratios while forgetting semantic membership.  Equation `(12)`
  is the exact lost sidecar.
- THM-2782 constructs positive fully marked semantic arm packets and
  explicitly withholds a full physical root-deck/current attachment.  It
  does not claim that those packets lie in the particular THM-2672 facet
  `(2)--(5)`.

The source-to-target contract is now:

| item | content |
|---|---|
| source | all delayed components of one carry-`12` fixed facet |
| target | THM-2749 fully marked rail-eight source carrier |
| proposed map | literal same-clock physical intersection |
| preserved | physical coordinate, rail, clock, source carry/root, target-label facet |
| first failure | anchor polarity `F_(ell,7)` versus `F_(ell,7)^c` |
| surviving boundary | `503` one-sided closure seams |
| needed next sidecar | a lawful configuration/edge switch opening a seam |

## 7. Scope

The theorem does **not** prove:

- that another source carry in configuration `(2)` is disjoint;
- that another rail, edge, sector, `kappa`, `h`, or base cell is disjoint;
- that the full union of THM-2672 configurations misses the semantic
  carrier;
- that a closure seam cannot open after a lawful switch;
- that THM-2782's semantic arms are zero;
- a row exclusion or LRC(14).

This is exactly the fixed-configuration versus union boundary already
required by THM-2672 and THM-2687.

## 8. Exact companion

Run

```bash
python 04-computation/lrc14_rail8_fixed_configuration_all_component_no_go_thm2794.py
python -O 04-computation/lrc14_rail8_fixed_configuration_all_component_no_go_thm2794.py
```

Both modes byte-match the stored transcript.  The companion uses exact
integer/rational interval arithmetic, explicit exception gates, no floating
point, and no truth-bearing Python assertions.  It pins the audited
THM-2785 and THM-2749 helpers; checks the primitive content, metadata, root,
carry, missing label, and active labels; rebuilds both positive vectors;
checks every base and marked piece; enumerates all `879,237` delayed
components; and independently verifies `(12)` before prefix refinement.

An immutable independent hostile audit rebuilt the universe and both
enumeration routes.  It reproduced the positive clocks, every piece and
component count, all `503` delayed seams and `1,011` base seams with their
one-sided orientation, normal/optimized/stored replay, dependency hashes,
and LF hashes.  The audit also exposed `(6)--(7)` as the actual first
obstruction; that typing repair is incorporated above.

QED.
