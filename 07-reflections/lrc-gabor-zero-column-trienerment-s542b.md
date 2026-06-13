---
source: codex-2026-06-01-S542
status: finite Gabor trienerment construction plus bounded audit
tags:
  - lonely-runner
  - gabor
  - trienerment
  - sector-vectors
  - tournament-analysis
  - time-frequency
  - fixed-clock-menu
---

# LRC Gabor Trienerment: The Zero-Column Target

The web search made one distinction useful.  "Gabor angle" can mean the
literal orientation/phase parameters of Gabor patches, but the LRC-relevant
version is the finite time-frequency one: Gabor systems on cyclic groups use
translations and modulations on a finite phase space, and finite Gabor frame
theory studies the geometry, coherence, and uncertainty of those atoms.  That
is the same kind of object as our sector/harmonic LRC picture.

External anchors used for this interpretation:

```text
Pfander, "Gabor frames in finite dimensions":
  finite cyclic/finite Abelian Gabor frames, coherence, geometry.
  https://www.mathematik.uni-marburg.de/~pfander/pubs/FiniteGaborFramesChapter.pdf

Feichtinger-Kozek-Luef, "Gabor analysis over finite Abelian groups":
  finite time-frequency planes, projective representations, symplectic Poisson.
  https://arxiv.org/abs/math/0703228

Grafakos-Sansing, "Gabor frames and directional time-frequency analysis":
  directionally sensitive Gabor/Radon atoms, useful for the "angle" reading.
  https://grafakos.missouri.edu/preprints/GrafakosSansing9.pdf

PsychBench Gabor element docs:
  literal Gabor patch orientation and phase parameters; useful but less central.
  https://www.psychbench.org/elementdocs/gabor
```

Upstream had already posed exactly the right next object:

```text
vertices = (sector, harmonic) atoms.
```

This session makes a concrete finite version and asks whether it actually
restricts the realizable trienerment classes.

## The Map

Start with a sector-vector:

```text
c = (c_0,...,c_{n-1}),  sum c_i=n-1.
```

Use the smallest LRC-visible finite Gabor window: two adjacent sectors.  For
window `a` and harmonic `m`, define:

```text
G_c(a,m) = c_a zeta^(ma) + c_{a+1} zeta^(m(a+1)).
```

The observer danger window is `a=n-1`, covering sectors `n-1` and `0`.  Hence:

```text
open LRC sector target c_0=c_{n-1}=0
iff
G_c(n-1,m)=0 for all m.
```

This is the cleanest Gabor translation I found.  The target is not "large
Fourier coefficient" or "low frame bound" first.  It is a marked zero column
in a finite time-frequency grid.

## The Trienerment

The Gabor trienerment has vertices `(a,m)`.

The pairwise observable is the complex phase of `G_c(a,m)`.

The switch is:

```text
tie     if either coefficient is zero/unresolved or phases lie in one
        atom-resolution cell;
direct  otherwise by half-turn phase order.
```

The tie Hamiltonian path is lexicographic atom order after the fixed
`(sector, harmonic)` scaffold.  I canonicalize under the dihedral sector
scaffold, with harmonic sign reversal under reflection, not full unmarked
`S_{n^2}` isomorphism.  Full unmarked isomorphism would throw away the observer
column and therefore destroy the LRC predicate.

Fingerprints recorded:

```text
scaffold class count,
zero observer atoms,
tie count,
directed triangle count,
out-score and tie-degree histograms.
```

## Computation

Artifact:

```text
04-computation/lrc_gabor_zero_column_s542b.py
05-knowledge/results/lrc_gabor_zero_column_s542b.out
```

Global sector-vector image:

```text
n=3: raw=6   raw_dihedral_orbits=2   Gabor scaffold-classes=5
n=4: raw=20  raw_dihedral_orbits=4   Gabor scaffold-classes=18
n=5: raw=70  raw_dihedral_orbits=10  Gabor scaffold-classes=44
n=6: raw=252 raw_dihedral_orbits=26  Gabor scaffold-classes=178
```

Good image:

```text
n=3: good_raw=1   good_Gabor_classes=1
n=4: good_raw=4   good_Gabor_classes=4
n=5: good_raw=15  good_Gabor_classes=11
n=6: good_raw=56  good_Gabor_classes=37
```

Observer zero-atom histograms:

```text
n=3: ((0,5),   (3,1))
n=4: ((0,14),  (1,2),  (4,4))
n=5: ((0,55),  (5,15))
n=6: ((0,172), (1,24), (6,56))
```

The single-zero cases at even `n` are harmonic cancellation in a nonempty
two-sector window.  That is a good warning: the target is not "some Gabor
coefficient vanishes."  It is "the marked observer column vanishes in every
harmonic."

## Fixed-Clock Menus

Bounded primitive scans:

```text
n=4, B<=10:
  speed_sets=109, open_good=108, wall_only=1, missing=0
  Gabor classes seen=18/18, good seen=4/4

n=5, B<=8:
  speed_sets=69, open_good=67, wall_only=2, missing=0
  Gabor classes seen=44/44, good seen=11/11

n=6, B<=7:
  speed_sets=21, open_good=20, wall_only=1, missing=0
  Gabor classes seen=135/178, good seen=22/37
```

The AP row is wall-only in all three cases.  So the Gabor object inherits the
same compactification debt as the sector, section, and symbolic views: open
zero-column classes catch generic witnesses, but AP needs a boundary version.

## What This Adds To Trienerments

The runner trienerment from the upstream note says:

```text
tie = near pair;
observer lonely = observer tie-degree 0.
```

The Gabor trienerment reverses the local tie story:

```text
observer lonely = observer-window Gabor atoms are zero,
                  hence phase-unresolved,
                  hence a zero/tie block appears in Gabor space.
```

Numerically, the Gabor tie count rises at observer-tie-free chambers:

```text
n=4: observer-tie-free avg Gabor ties=73.3, observer-tied=58.4
n=5: observer-tie-free avg Gabor ties=148.4, observer-tied=97.9
n=6: observer-tie-free avg Gabor ties=304.3, observer-tied=252.3
```

That is the conceptual payoff.  LRC is local tie deletion in real pair-space
but local tie creation in time-frequency atom-space.  The Gabor angle is not
just a metaphor for uncertainty; it is a dual way to mark the same event.

## Relation To Sector-Vector Realizability

HYP-2028 says all sector-vectors are existentially realizable.  Therefore all
finite Gabor trienerments induced by sector-vectors are existentially
realizable too.  The lift does not create a global obstruction.

The useful restriction is:

```text
fixed primitive clock V
  -> finite menu of sector-vectors
  -> finite menu of Gabor trienerment classes
  -> must hit the marked zero-column family or a wall compactification.
```

This agrees with the broader lesson from HYP-2023/HYP-2024/HYP-2029:
compression helps only when the observer target and arithmetic memory survive.

After this work was drafted, upstream HYP-2031 landed an unanchored
Gabor-cell trienerment.  Its negative signal is important: coarse
time-frequency support and phase fingerprints have mixed good/bad classes.
The zero-column construction here should be read as the anchored answer to
that: the Gabor lift becomes target-visible only when the observer danger
window is built into the atom definition.

## Assumption Challenge

Considered vertex sets:

```text
runners,
runner pairs,
fixed sectors,
section boundaries,
wall events,
Fourier modes,
Gabor atoms (sector, harmonic),
proof obligations.
```

Chosen vertex set:

```text
Gabor atoms (two-sector window a, harmonic m).
```

Predicate preserved:

```text
open LRC target c_0=c_{n-1}=0
as a marked all-zero observer column.
```

Information destroyed:

```text
runner identities,
positions inside sectors,
exact wall timing,
and full continuous Fourier phase of individual runners.
```

Challenged assumption:

```text
the Gabor route should be phrased as generic frame density or scalar
frequency-floor first.
```

Those may be useful later, but the exact finite object is the zero-column
target.  It is simpler and better anchored.

## Next

1. Replace floating phase comparison with exact cyclotomic comparison.
2. Add compactified wall atoms: half-weight endpoint windows or wall-owner
   labels should make AP rows hit a closed zero-column analogue.
3. Test `n=7` global Gabor classes after optimizing canonicalization.
4. Compare Gabor zero-column classes with HYP-2029 decorated target words.
5. Use the finite Heisenberg/metaplectic scaffold to decide which harmonic
   relabelings are legitimate isomorphisms for the LRC window.
