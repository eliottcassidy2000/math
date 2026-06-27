---
id: HYP-3109
title: Lee-Yang root-curve extremality and the zero-real ear map for the LRC miss-count PGF
status: EVIDENCE / exact anchored scout plus proof-program hypotheses; not a proof
source: codex-2026-06-27
script: 04-computation/lrc_lee_yang_root_curve_earmap_codex_20260627.py
result: 05-knowledge/results/lrc_lee_yang_root_curve_earmap_codex_20260627.out
reflections:
  - lee-yang-extremality-root-curve-earmap-codex-20260627
related:
  - HYP-3103
  - HYP-3104
  - HYP-3105
  - HYP-3106
  - HYP-3108
  - HYP-3085
  - HYP-3043
  - HYP-2879
  - THM-534
  - THM-556
  - THM-576
  - THM-577
  - OPEN-Q-108
---

# HYP-3109: Lee-Yang Root-Curve Extremality And The Zero-Real Ear Map

## Claim

The next Lee-Yang formulation of the LRC14 extremality problem should use the
whole miss-count PGF

```text
G_N(z) = sum_t q_t z^t,  q_t = P(# empty inner sectors = t),
```

not only the single value `G_N(0)=p0` or finitely many factorial moments.  The
proof-facing object is the root locus plus a controlled-forgetting sidecar:

```text
lee_yang_packet =
  (root_zero_locus,
   zero-free segment clearance,
   nearest-root radius,
   apex-7 angle error,
   quartic potential profile,
   root-collision wall,
   zero-real ear component)
```

The anchored scout supports the strong version of HYP-3103: the `#real=0`
root stratum is not just a label.  In the bounded anchored residue bank it is
one connected, articulation-free one-swap component containing `consec_8`, and
`consec_8` is the unique `p0` and `L_yK8` leader.

Post-rebase integration note: HYP-3108 now reserves the broader
Lee-Yang/Savitch/Bravais/ear-lattice synthesis lane.  This HYP supplies the
completed root-curve and zero-real ear-map scout for that lane to import,
while incoming HYP-3107 remains the Lean proof-frontier ledger.

## Exact Scout Readout

Artifact:

```text
04-computation/lrc_lee_yang_root_curve_earmap_codex_20260627.py
05-knowledge/results/lrc_lee_yang_root_curve_earmap_codex_20260627.out
```

Named rows:

```text
consec_8:
  p0=0.32721, L_yK8=3.58231, #real=0
  min |z| = 1.48858
  dist(roots, [-1,0]) = 0.91189
  upper angles = [55.3, 111.0, 154.2]
  mean apex-7 angle error = 4.042 deg
  quartic fit: double-well, b=-0.10232, lambda=+0.00566

break_spread:
  p0=0.04528, #real=2, min |z|=0.12125
  dist(roots, [-1,0]) = 0

random_spread:
  p0=0.01152, #real=4, min |z|=0.04724
  dist(roots, [-1,0]) = 0
```

Anchored bounded bank `0 union A`, `|A|=7`, `A subset {1,...,13}`:

```text
count = C(13,7) = 1716

#real=0: n=290,  max_p0=0.32721, mean_p0=0.10551, max_L_yK8=3.58231
#real=2: n=1426, max_p0=0.20833, mean_p0=0.05732, max_L_yK8=2.37381

corr(p0, -#real)             = +0.483
corr(p0, min root modulus)   = +0.899
corr(p0, segment clearance)  = +0.681
corr(p0, -apex7 angle error) = +0.485
```

The zero-real component:

```text
vertices = 290
edges = 1048
components = 1
cycle_rank = 759
articulation_points = 0
consec_component_size = 290
boundary edges from #real=0 to #real=2 = 10084
```

This makes the root stratum an actual graph object.  If its one-swap edges are
bidirected, the component is strongly connected and therefore admits a directed
ear decomposition.  The theorem content is not that tautology; it is whether a
small seed around the consecutive/AP root packet plus root-collision ears
generates the full zero-real component before crossing into the `#real=2`
wall.

## Phi4 / Quartic Potential Signal

The prompt's `exp(-lambda S^4 - b S^2)dS` suggests a compact signal.  For each
row, symmetrize the miss-count law around `N=3` and fit

```text
V(s) = -log r_s ~= c + b s^2 + lambda s^4.
```

For `consec_8`, this gives a clean double-well diagnostic:

```text
b=-0.10232, lambda=+0.00566, well_radius ~= 3.005, r2=0.9783.
```

This is not a proof model and should not be over-read.  In the consecutive
ladder `k=9..13`, the naive quartic fit becomes unstable, which is itself a
signal: the effective potential needs a moving center, a sextic term, or an
unsymmetrized field.  The useful retained sidecar is therefore
`quartic_phi4_profile`, not the scalar `b` alone.

## Bold Hypotheses

1. **Lee-Yang extremality theorem.**  Among anchored bounded LRC sector rows,
   the `p0` maximizer is characterized by a root-locus condition stronger than
   `#real=0`: no zero on `[-1,0]`, large root product at the origin, and an
   apex-7-constrained upper-half-plane arc.  `#real=0` is necessary but not
   sufficient; nearest-root radius alone is also not sufficient.

2. **Root-collision wall.**  The transition from extremal/correlated packets to
   spread packets is a collision of a conjugate pair onto the real segment
   `[-1,0]`.  The discriminant/resultant of `G_N` should detect the same
   finite wall as the k=8 break and the one-swap nonglobal traps in HYP-3104.

3. **Zero-real ear generation.**  The 290-row anchored `#real=0` component has
   an ear decomposition whose seed should be the AP/consec root packet.  The
   ears should correspond to controlled root deformations; crossing an ear
   endpoint into the boundary is the birth of a real zero.  This turns the
   prompt's strong-connectivity/ear theorem into a concrete LRC signal.

4. **Odd-ear and factor-critical analogy.**  The right factor-critical object is
   not the runner graph.  It is the root-collision boundary graph: odd ears may
   mark parity-protected ways a conjugate root pair can hit the real segment
   and return.  Test by building the boundary graph between `#real=0` and
   `#real=2` rows and checking whether minimal collision cycles are odd.

5. **Nested-ear and series-parallel analogy.**  If the zero-real component has a
   nested-ear decomposition after quotienting by AP symmetries, it may explain
   why one-swap exchange is non-transitive yet finite: the obstruction is a
   nested local wall, not a global scalar descent.

6. **Two-map synthesis.**  The gold is in the HYP-3043 lens map and the
   HYP-3104/HYP-3106 signal/functor maps: a proof must retain the root sidecar
   precisely when a scalar moment quotient would otherwise forget the
   zero-collision coordinate.

## New Signals To Carry Forward

- `root_zero_locus`: full multiset of PGF roots.
- `lee_yang_segment_clearance`: minimum distance from roots to `[-1,0]`.
- `nearest_root_modulus`: coverage-pole clearance, useful but not decisive.
- `apex7_angle_error`: deviation of upper roots from `360/7` multiples.
- `quartic_phi4_profile`: `(b, lambda, well_radius, fit_status)`.
- `root_collision_discriminant`: discriminant/resultant wall for real-zero birth.
- `zero_real_ear_component`: component id, ear distance from AP, cycle rank.
- `root_boundary_degree`: number of one-swap moves from `#real=0` to `#real=2`.
- `odd_ear_collision_parity`: parity of minimal wall-crossing ears.
- `nested_ear_width`: width/depth of an ear decomposition after quotienting.

## Tournament Analysis

Vertices are signals/proof obligations, not runners.  Considered alternatives
were runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, cover arcs, Fourier modes, PGF roots, quartic potential
wells, root-collision events, ear graph components, and proof obligations.

The selected vertex set is:

```text
whole_pgf_curve
root_zero_locus
quartic_phi4_fit
zero_stratum_ear_graph
root_collision_discriminant
nearest_zero_clearance
real_root_count
factorial_moments
single_value_p0
```

Pairwise observable: which signal retains more proof-relevant coordinates
under predicate preservation, whole-curve retention, root geometry, phase-wall
detection, local graph structure, controlled-forgetting safety,
computability, and hypothesis yield.  Tie path follows the displayed order in
the script.

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles = 0
Hamiltonian_path_count = 1
priority path =
  root_zero_locus > whole_pgf_curve > zero_stratum_ear_graph
  > root_collision_discriminant > quartic_phi4_fit
  > nearest_zero_clearance > real_root_count
  > factorial_moments > single_value_p0
```

The quotient preserves the proof question "which coordinate must be retained
next?" and destroys row identity.  That is intentional for this exploratory
stage.

## Caveats

- The root solver is a self-contained Durand-Kerner implementation to avoid
  external numerical dependencies in this checkout.  Degree is only six, but
  future certification should replace numerical roots by interval or
  Sturm/discriminant certificates.
- The anchored bank is a bounded residue model, not the full LRC14 universe.
- The quartic fit is a diagnostic inspired by the phi4 density, not an
  asserted continuum limit.
- The zero-real ear component shows the right graph shape; the actual proof
  must give labelled ears that preserve the LRC predicate and identify the
  destroyed coordinates.
