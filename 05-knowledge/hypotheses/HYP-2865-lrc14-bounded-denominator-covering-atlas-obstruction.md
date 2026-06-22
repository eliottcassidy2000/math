---
id: HYP-2865
title: LRC14 bounded-denominator covering atlas -- fixed caps are impossible; scaled residue atlases remain useful
status: CORRECTED / fixed-denominator version refuted by THM-566; scaled-atlas version open
source: codex-2026-06-22-S93
depends_on:
  - THM-566
  - THM-523
  - THM-524
  - THM-366
related:
  - HYP-2052
  - HYP-2864
  - THM-565
  - OPEN-Q-108
results:
  - 05-knowledge/results/lrc14_covering_bounded_denominator_obstruction_codex_s93.out
---

# HYP-2865 -- Bounded-Denominator Covering Atlas Obstruction

## Corrected Claim

The attractive statement

```text
every LRC14 covering 13-set has a witness tau=a/D with D <= B0
```

is false for every fixed `B0`, even after the THM-523 reduction to covering
sets.

The correct residue-atlas target must be scaled or conditioned: it may be true
for bounded speed boxes, non-divisor-loaded rows, or with `D` bounded by the
first modulus not killed by a speed.  A fixed absolute denominator cap cannot
close LRC(14).

## Refutation Mechanism

THM-566 gives the exact obstruction.  For any `B`, set

```text
S_B = {1,2,...,11,13,84*lcm(1,...,B)}.
```

This is primitive and covering.  For every `D <= B`, the last speed is divisible
by `D`, so every `a/D` puts that runner at the observer.  Therefore no
denominator `D <= B` can witness loneliness.

This is the covering-hard-core version of the older HYP-2052 no-finite-sieve
warning.  The new point is that "restrict to covering sets" does not rescue a
uniform denominator cap.

## Positive Local Evidence

The fixed cap is false, but the scout explains why it looked plausible.

Named rows:

```text
{1,...,11,13,84}:       tau=17/41
{1,...,11,13,84*6}:     tau=22/53
{1,...,12,182}:         tau=2/27
{1,2,3,4,5,7,...,13,84}: tau=4/23
```

Tower scout for `S_m={1,...,11,13,84m}`, `m<=5000`, found zero failures with
search bound `D<=100` and max least denominator `67`.  The first row exceeding
`41` is already `m=6`, with witness `22/53`.  The first rows at the displayed
high denominators are:

```text
D=41: m=1
D=53: m=6
D=55: m=53
D=65: m=416
D=67: m=1681
```

AP-repair rows `({1,...,13}\{drop}) union {14m}`, `m<=500`, had zero failures
with `D<=100`, max denominator `55`.  A 300-row random covering-obligation bank
up to speeds `10^6` had zero failures with `D<=120`, max denominator `32`.

So bounded denominators are a strong diagnostic for typical and bounded boxes,
not a theorem-level global cap.

## Synthesis With Current Routes

- **THM-523:** covering-set reduction remains the right first split, but the
  residual can still be divisor-loaded against every fixed denominator.
- **THM-524:** binding-pair switches are the exact global max-min language; a
  scaled denominator atlas should be phrased through binding crossings or first
  unblocked moduli, not arbitrary fixed `D`.
- **HYP-2864:** the same arithmetic appears as a quotient/gcd branch.  When a
  sheet or denominator argument fails, it exposes divisor loading or a bounded
  quotient, not a mysterious analytic obstruction.
- **THM-565:** finite-ruler sampling is still needed for divisor-loaded towers;
  fixed rational probes cannot replace the Node-1/Node-3 floor route.
- **Selberg/minorant tools:** still useful for per-row finite certificates, but
  cannot be compressed into a universal denominator list.

## Tournament Analysis And Assumption Challenge

Vertices considered:

```text
runners, denominators D, numerator residues a mod D, forbidden residue arcs,
binding-pair crossings, sheet quotients b/g, safe components, Fourier modes,
and proof obligations.
```

Chosen quotient: proof carriers.  Pairwise observable:

```text
(rules out false route, preserves LRC witness predicate,
 finite-check value, formalization readiness)
```

Switch/gauge: lexicographic comparison of the observable.  Hamiltonian path
from the script:

```text
covering_denominator_no_go
> scaled_denominator_or_speed_bound
> THM523_covering_reduction
> THM524_binding_pair_switches
> HYP2864_sheet_gcd_quotient
> THM565_three_gap_floor
> fixed_B_residue_atlas
> raw_runner_vertices
```

Fingerprint: transitive tournament, score histogram `{0,1,2,3,4,5,6,7}`, no
directed 3-cycles, singleton SCCs, one Hamiltonian path.

Challenged assumption: tournament vertices should be runners or fixed residue
sections.  That quotient hides the decisive divisor-loading obstruction.  The
proof-carrier quotient preserves the LRC predicate "there exists a lonely
witness" and destroys only the irrelevant identity of individual runners.

## Next Scaled Target

Replace the false fixed-cap statement with:

```text
q_min(S) is controlled by the first modulus D for which no speed of S is
divisible by D, plus a bounded local repair term.
```

This is the HYP-2052+ direction specialized to LRC14 covering sets and should
be tested against the tower `S_m`, THM-524 binding denominators, and the
HYP-2864 gcd-quotient fallback.
