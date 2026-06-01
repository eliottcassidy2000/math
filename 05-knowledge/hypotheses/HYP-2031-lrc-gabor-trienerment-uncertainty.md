---
id: HYP-2031
status: OPEN
source: codex-2026-06-01-S542
related:
  - HYP-2022
  - HYP-2024
  - HYP-2027
  - HYP-2028
  - HYP-2029
---

# HYP-2031: Gabor-cell t(r)ienerments expose a ternary uncertainty axis, but need observer anchoring

**Claim.** The time-frequency/Gabor vertex set proposed in HYP-2027 should be
treated as a t(r)ienerment, not forced into a tournament.  For denominator
`n` and sector occupancy vector `c`, take vertices

```text
(sector k, harmonic m), 0 <= k < n, 1 <= m < n.
```

Attach to each vertex the local Gabor atom key

```text
(amplitude, phase) = (c_k, m*k mod n).
```

The natural ternary relation is:

```text
larger amplitude wins;
equal amplitude is compared by cyclic phase majority;
equal phase, and antipodal phase for even n, is a tie.
```

Thus the third state is not noise.  It records exact Gabor degeneracy or
phase-opposition ambiguity, matching the t(r)ienerment model where a tie is a
bidirectional edge and contributes to the ternary Krawtchouk axis.

**Evidence.** `lrc_gabor_trienerment_s542.py` scans bounded primitive LRC
clocks for `n=3..7`.  It computes the Gabor t(r)ienerment fingerprint of each
sector vector: tie count/rate, ternary `B1 = 2M - 3*ties`, score histogram,
strict directed 3-cycles, SCC counts, and Hamiltonian path counts for
`n<=4`.

The full Fourier support uncertainty check holds in every seen vector:

```text
S * (F_nonzero + 1) >= n.
```

This is the correct support product because the DC coefficient is always
present for a nonzero occupancy vector.  The non-DC Gabor product is still a
useful diagnostic, but it is not the uncertainty theorem by itself.

Bounded fingerprint separation:

```text
n=3: classes=4,  pure good / pure bad / mixed = 0 / 3 / 1
n=4: classes=10, pure good / pure bad / mixed = 1 / 6 / 3
n=5: classes=12, pure good / pure bad / mixed = 2 / 8 / 2
n=6: classes=53, pure good / pure bad / mixed = 1 / 43 / 9
n=7: classes=16, pure good / pure bad / mixed = 0 / 13 / 3
```

So the unanchored Gabor fingerprint is not a target certificate.  It carries
real structure, but mixed good/bad classes remain.

**Main signal.** Good states tend to have lower full support product and lower
or comparable tie density than bad states, while cyclic phase creates genuine
directed 3-cycles for odd denominators.  For example:

```text
n=5 good: S*F_full avg 10.417, tie rate ppm avg 64912, strict 3-cycles avg 78
n=5 bad:  S*F_full avg 14.091, tie rate ppm avg 73206, strict 3-cycles avg 57

n=7 good: S*F_full avg 30.100, tie rate ppm avg 48781, strict 3-cycles avg 469
n=7 bad:  S*F_full avg 30.273, tie rate ppm avg 55975, strict 3-cycles avg 691
```

The cyclic gauge matters.  A linear phase order collapses the relation to a
total preorder and kills strict cycles everywhere.  Phase is circular, so the
phase comparator must be circular too.

**Connection to t(r)ienerments.** The Gabor quotient supplies a concrete LRC
source for the ternary state:

```text
forward/backward = amplitude or cyclic phase preference;
tie = equal atom, or antipodal atom when the phase comparison is undecidable.
```

This aligns with the trienerment writeup's ternary tiling picture and its
first Krawtchouk coordinate.  The measured `B1 = 2M - 3*ties` is exactly the
tie-axis statistic for this Gabor relation.

It also executes the Gabor-cell construction posed by the incoming LRC
t(r)ienerment/Gabor thread: the earlier work made ties equal nearness in
runner space; this probe moves the ternary relation to the joint
`(sector,harmonic)` phase space.

**Assumption challenge.** This session explicitly rejects the default that
tournament vertices must be runners or arcs.

Considered vertex sets:

```text
runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes,
Gabor cells, matroid circuits, proof obligations.
```

Chosen quotient:

```text
vertices = Gabor cells (sector, harmonic).
```

Predicate preserved:

```text
the underlying sector vector still records the LRC target c_0=c_{n-1}=0.
```

Information destroyed:

```text
runner identity, fixed-clock wall order, gate ownership, endpoint debt,
and enough observer anchoring that the compressed fingerprints become mixed.
```

Challenged assumption:

```text
time-frequency support alone might separate the LRC target.
```

The mixed classes refute that as a standalone proof route.  The Gabor relation
must be anchored to the observer or coupled to the event language of HYP-2029.

**Predictions.**

1. Observer-anchored Gabor cells, such as cells weighted by distance to sectors
   `0` and `n-1`, will reduce mixed fibers more than increasing the speed box.
2. Event-Gabor cells `(wall event, harmonic)` should connect the symbolic
   target-shift HYP-2029's decorated arithmetic word to HYP-2031's phase-cycle
   t(r)ienerment.
3. The tie-axis statistic `B1` should be useful in a vector of gauges, not as a
   scalar certificate.
4. Even denominators should keep extra antipodal ties; this is a feature,
   not a degeneracy to break.
5. A proof-shaped quotient should combine sector support uncertainty, cyclic
   phase cycles, and observer endpoint labels.

**Next tests.**

1. Add observer-windowed Gabor weights and retest mixed fingerprint classes.
2. Compare cyclic phase t(r)ienerments with event-word tournaments from
   HYP-2029 on the same fixed clocks.
3. Compute anchored isomorphism classes preserving the two observer-adjacent
   sectors instead of using only coarse score/fingerprint data.
4. Extend the scan to hard rows such as `n=14` AP/scalar-puncture families,
   where antipodal phase ties should be pronounced.

**Files.** `04-computation/lrc_gabor_trienerment_s542.py`;
`05-knowledge/results/lrc_gabor_trienerment_s542.out`;
`07-reflections/lrc-gabor-trienerments-and-time-frequency-uncertainty-s542.md`.
