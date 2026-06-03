---
source: codex-2026-06-03-S571
status: reflection + finite Burnside/Fourier audit
tags: [LRC, time, Burnside, Fourier, reset, representation, n14]
---

# Time Burnside and Fourier Quotients

The user framed time in LRC as the real object: not just a snapshot of runners,
but the whole period before reset, or the proof that there is no exact reset.
That is the right instinct.  The correction from HYP-2080 is that for primitive
integer speeds the raw reset time is always one lap, so the classifying object
is not reset length.  It is the folding of the time orbit.

S571 turns that into a finite Burnside/Fourier diagnostic.

## The Burnside Trap

If we discretise time to `Z/N` and let `X` be the lonely slots, the full cyclic
group `Z/N` acts on ambient time slots.  But it usually does not act on `X`:
shifting a lonely time can land on a non-lonely time.

So the tempting formula

```text
|X/G| = (1/|G|) sum_g |X^g|
```

is not a legal orbit count for the full group.  In practice it degenerates to
`|X|/N`, a density.

The legal replacement is:

```text
K = stabilizer of the whole lonely time word.
```

Then `X` is `K`-invariant, and Burnside gives `|X/K|=|X|/|K|`.

For the audited primitive n=14 rows, `K` is trivial.  That is itself useful:
the binary lonely word has no nontrivial time-shift symmetry on the chosen
grid.  Any meaningful quotient symmetry must live in a richer labelled object.

## Fourier Dual

The representation side is better behaved immediately.  The irreps of `Z/N`
are frequencies.  The lonely indicator `1_X` has a Fourier transform, and its
energy distribution says which clocks are active.

On the S571 grid `N=38640`, the dominant frequencies line up with the existing
clock story:

```text
near_AP_apex:          6, 12
S562_packet_n14:      21, 42
S562_packet_n14_lift: 42, 84
no_small_pinch_proxy:  2, 24
random_low_resonance: 34, 2
```

So the dual version of "which family does this speed set belong to?" is:

```text
which frequencies carry the lonely-time word?
```

The `n`-clock, pair-sum pinch clocks, dyadic packet gears, and generic
low-resonance clocks show up as frequency families.

The incoming HYP-2083 antipodal/summand-unit bridge should be folded into this
same clock atlas: unit-visible shells modulo `2n-1` are labelled event orbits,
while nonunit holes are the composite-clock branch the binary lonely word
forgets.

## What This Means For n=14

The finite-time model suggests a layered proof object:

```text
raw time word       -> lonely density and wall/non-wall split
stabilizer K        -> legal Burnside quotient
Fourier spectrum    -> clock family / gear support
owner-labelled word -> no-return handoff and endpoint core
```

The binary word is not enough for a proof because its stabilizer is usually
trivial.  The next layer should label time events by:

```text
runner owner, pair-sum denominator, G component, wall endpoint, core endpoint.
```

That is where Burnside can become structural rather than merely diagnostic.

## Handoff

Use S571 as the Fourier front end for S570:

```text
for a speed set:
  find primal witness or endpoint core;
  build the labelled time-event word;
  compute stabilizer subgroups and Fourier energy by labels;
  run Tournament Analysis on active frequencies or owner events.
```

If a purported n=14 counterexample has to be both binary-stabilizer trivial and
owner-event highly constrained, that tension may be exploitable: the time word
looks generic, but the endpoint core must act resonantly.

**Artifacts:** `04-computation/lrc_time_burnside_fourier_s571.py`,
`05-knowledge/results/lrc_time_burnside_fourier_s571.out`, HYP-2085; related:
HYP-2083.
