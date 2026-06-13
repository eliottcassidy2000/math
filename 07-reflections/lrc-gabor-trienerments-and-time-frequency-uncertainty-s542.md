---
source: codex-2026-06-01-S542
status: computational probe plus synthesis with t(r)ienerments
tags:
  - lonely-runner
  - gabor
  - time-frequency
  - trienerment
  - symbolic-dynamics
  - tournament-analysis
  - uncertainty
---

# LRC Gabor T(r)ienerments And Time-Frequency Uncertainty

HYP-2027 left a specific handoff: build the `(sector, harmonic)` Gabor
tournament and test the uncertainty restriction.  The t(r)ienerment writeup
changes the right answer: this object should not be a tournament first.  Equal
or phase-opposed Gabor atoms naturally produce a third state, and the
t(r)ienerment language already knows how to treat that state as bidirectional
for Hamiltonian paths.

## Construction

For a denominator `n`, an open LRC chamber has sector occupancy vector

```text
c = (c_0, ..., c_{n-1}), sum c_k = n-1.
```

The Gabor vertices are:

```text
(sector k, harmonic m), 0 <= k < n, 1 <= m < n.
```

The local atom key is exact:

```text
(amplitude, phase) = (c_k, m*k mod n).
```

The relation uses the cyclic phase gauge:

```text
larger amplitude wins;
equal amplitude uses the regular cyclic phase tournament;
equal phase ties;
for even n, antipodal phase also ties.
```

This last line is the t(r)ienerment move.  Antipodal phase is not a nuisance
to break by an arbitrary lexicographic order.  It is an actual symmetry
ambiguity, so it belongs in the third state.

## Tournament Analysis Declaration

```text
vertices:
  Gabor cells (sector, harmonic), not runners or arcs

pairwise observable:
  amplitude and phase residue of the local atom c_k zeta_n^(mk)

switch/gauge:
  amplitude priority followed by cyclic phase majority

tie Hamiltonian path:
  lexicographic (phase, sector, harmonic) inside tie classes

fingerprints:
  tie counts/rates, ternary B1, score histograms, strict 3-cycles,
  SCC counts, and Hamiltonian path counts for <=12 vertices
```

The script is:

```text
04-computation/lrc_gabor_trienerment_s542.py
05-knowledge/results/lrc_gabor_trienerment_s542.out
```

## Main Output

The full Fourier support uncertainty product holds in every bounded sample:

```text
S * (F_nonzero + 1) >= n.
```

Counting the DC coefficient matters.  If the zero harmonic is omitted, a
single-sector occupancy has `S=1` and `F_nonzero=n-1`, so the product is
`n-1`; that is a non-DC Gabor support diagnostic, not the full uncertainty
theorem.

Bounded scan summary:

```text
n=3, max_speed=10: 154 open cells, 6 sector vectors, 1 good vector
n=4, max_speed=10: 1018 open cells, 20 sector vectors, 4 good vectors
n=5, max_speed=8: 976 open cells, 56 sector vectors, 12 good vectors
n=6, max_speed=7: 438 open cells, 105 sector vectors, 15 good vectors
n=7, max_speed=7: 188 open cells, 87 sector vectors, 10 good vectors
```

Fingerprint separation:

```text
n=3: classes=4,  pure good / pure bad / mixed = 0 / 3 / 1
n=4: classes=10, pure good / pure bad / mixed = 1 / 6 / 3
n=5: classes=12, pure good / pure bad / mixed = 2 / 8 / 2
n=6: classes=53, pure good / pure bad / mixed = 1 / 43 / 9
n=7: classes=16, pure good / pure bad / mixed = 0 / 13 / 3
```

This is the key negative: unanchored Gabor fingerprints are not pure target
certificates.

## What Survives

The Gabor t(r)ienerment still gives useful signals.

Good states usually sit at lower full support product:

```text
n=5 good/bad average S*F_full: 10.417 / 14.091
n=6 good/bad average S*F_full: 18.800 / 21.200
n=7 good/bad average S*F_full: 30.100 / 30.273
```

Tie density is generally lower in good states, except for a small `n=4`
near-symmetry:

```text
n=5 good/bad average tie ppm: 64912 / 73206
n=6 good/bad average tie ppm: 137931 / 147229
n=7 good/bad average tie ppm: 48781 / 55975
```

Cyclic phase produces real cycles for odd denominators:

```text
n=5 good/bad average strict 3-cycles: 78 / 57
n=7 good/bad average strict 3-cycles: 469 / 691
```

Even denominators instead emphasize antipodal ties.  That is exactly what we
should expect from the phase circle and from the repo's earlier even-denominator
endpoint debt/tie themes.

## Why The Linear Gauge Was Rejected

The first version ordered equal amplitudes by the integer phase residue.  That
made every relation a total preorder, so strict directed 3-cycles were zero in
every sample.  It was computationally tidy and conceptually wrong: phase is
circular.  The cyclic gauge turns the Gabor layer back into an actual
phase-object and lets odd cycles appear.

This is a useful methodology lesson for the t(r)ienerment program.  A tie
state is meaningful only if the binary states also respect the geometry that
created the tie.

## Assumption Challenge

This session considered:

```text
runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes,
Gabor cells, matroid circuits, proof obligations.
```

The selected vertices are Gabor cells.  They preserve the LRC target only
through the underlying sector vector:

```text
c_0 = c_{n-1} = 0.
```

They destroy:

```text
runner identity,
wall order,
gate owner,
endpoint debt,
and too much observer anchoring after fingerprint compression.
```

The challenged assumption was:

```text
time-frequency support alone might separate the LRC target.
```

The mixed classes say no.  Gabor is a strong coordinate system, but it is not a
standalone proof certificate.

## Synthesis With The Recent Threads

HYP-2022 says sector occupancy is DFT-dual to the Fourier/resonance picture.
HYP-2028 says every raw sector vector is existentially realizable, so sector
existence is not the obstruction.  The recent HYP-2029 threads give the two
sides this probe needed: the symbolic target-shift view says the fixed clock is
a compactified word needing consistency labels, while the t(r)ienerment/Gabor
view poses exactly this `(sector,harmonic)` ternary lift.

HYP-2031 slots between them:

```text
sector occupancy -> Gabor t(r)ienerment -> uncertainty/tie fingerprint.
```

It adds a serious phase/tie layer, but it still has to be attached to the
fixed-clock event word or observer endpoint labels.

## Next Concrete Move

Do not merely enlarge `max_speed`.  The next useful script should anchor the
Gabor layer:

```text
vertices = (observer-windowed sector, harmonic)
or
vertices = (wall event, harmonic)
```

Then measure:

```text
mixed fingerprint classes,
tie-axis B1,
strict phase cycles,
event-word bad SCC leakage,
and AP wall behavior.
```

If observer-windowed or event-Gabor cells reduce mixed classes, the Gabor angle
becomes a proof microscope rather than a global descriptor.
