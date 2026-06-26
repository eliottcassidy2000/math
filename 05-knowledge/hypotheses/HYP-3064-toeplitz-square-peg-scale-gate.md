---
id: HYP-3064
title: Toeplitz square-peg scale gate and four-witness carrier for LRC14
status: SYNTHESIS
source: codex-2026-06-26-S229
tangent: T1146
technique: LTI-211
tournament_technique: LTT-109
script: 04-computation/lrc14_toeplitz_square_peg_carrier_s229.py
result: 05-knowledge/results/lrc14_toeplitz_square_peg_carrier_s229.out
forum_posts:
  - poke-forum/posts/20260626-lrc14-toeplitz-square-peg-scale-gate/post.md
external_anchors:
  - https://arxiv.org/abs/1611.07441
  - https://arxiv.org/abs/1001.0186
  - https://arxiv.org/abs/2203.02613
  - https://arxiv.org/abs/2404.05179
  - https://arxiv.org/abs/2407.07798
  - https://arxiv.org/html/2407.20412v1
related:
  - HYP-3063
  - HYP-3062
  - HYP-3061
  - HYP-3060
  - HYP-3059
  - HYP-3058
  - HYP-3057
  - HYP-3056
  - HYP-3054
  - HYP-3053
  - HYP-3037
  - HYP-2997
  - HYP-2974
  - HYP-2963
  - THM-572
---

# HYP-3064: Toeplitz Square-Peg Scale Gate

## Claim

Toeplitz's square-peg conjecture should be imported into the LRC14 proof stack
as a **four-witness configuration carrier with an explicit positive-scale
gate**, not as a general geometric analogy.

The square-peg problem asks whether every continuous Jordan curve contains the
four vertices of a non-degenerate square.  The general continuous case remains
open, while many regular, generic, nearby, Lipschitz-graph, rectangular, and
periodic variants have positive results.  The recurring proof danger is that
compactness or approximation can produce a limiting square of side length zero.

The LRC14 analogue is direct: a quotient or certificate may show an approximate
safe witness but lose strictness when the limit collapses to a boundary/AP-GW
zero-open atom.  Therefore future final-assembly packets should retain a
positive-scale sidecar before trusting any four-witness, rectangle, Toeplitz,
homological, or Floer-style obstruction.

## Evidence

The exact scout `04-computation/lrc14_toeplitz_square_peg_carrier_s229.py`
records the square constraints in diagonal form:

```text
ambient_variables=4 plane points = 8 real coordinates
constraint_count=4 equalities plus one open inequality
square_family_dimension=4 (center 2 + diagonal vector 2)
midpoint_balance: p0+p2=p1+p3
equal_diagonal_radius: |p0-p2|^2=|p1-p3|^2
quarter_turn_orthogonality: (p0-p2) dot (p1-p3)=0
positive_scale_gate: p0 != p2 and p1 != p3
```

It also records the finite symmetry side:

```text
D4_group_size=8
all_pair_partitions=3
opposite_pair_partition_orbit_size=1
cyclic_order_orbit_size=8
all_labelled_orders=24
diagonal_pairing_is_not_a_raw_four_point_choice=True
```

The square is not just four points.  It is two antipodal pairs sharing a
midpoint, equal radius, quarter-turn orthogonality, cyclic order, and positive
scale.  Forgetting any of those fields turns a theorem-facing configuration
into a scalar shadow.

## Current External Status

Sources checked on 2026-06-26:

- Tao's integration approach proves new two-Lipschitz-graph cases and isolates
  the small-square compactness problem: <https://arxiv.org/abs/1611.07441>.
- Matschke's survey/results frame generic and open classes:
  <https://arxiv.org/abs/1001.0186>.
- Chambers proves square pegs for curves close to a `C^2` curve under a
  curvature-scale condition: <https://arxiv.org/abs/2203.02613>.
- Greene-Lobb Floer work gives rectangular/square-peg lanes for smooth or
  rectifiable settings: <https://arxiv.org/abs/2404.05179> and
  <https://arxiv.org/abs/2407.07798>.
- Hugelmeyer resolves Tao's periodic square-peg variant by a Floer-homology
  torus-intersection argument: <https://arxiv.org/html/2407.20412v1>.

These are not LRC inputs.  They are proof-design anchors for which coordinates
must not be forgotten.

## LRC14 Translation

The square constraints become packet fields:

```text
toeplitz_square_scale_gate
midpoint_balance_residue
diagonal_equal_radius_residue
quarter_turn_residue
ordered_quad_collapse_mode
d4_orbit_word
toeplitz_psd_bridge_degree
```

Suggested reading:

- `midpoint_balance_residue` is the AP/GW endpoint-credit analogue.
- `diagonal_equal_radius_residue` is the equal-slack / danger-radius analogue.
- `quarter_turn_residue` is the Haar rectangle-zeta / hourglass analogue.
- `ordered_quad_collapse_mode` records whether a four-witness object is
  collapsing to a point, boundary pair, diagonal pair, AP/GW atom, K33 state
  lift, or true strict-open witness.
- `toeplitz_square_scale_gate` is the strictness condition: the witness must
  have positive safe mass or positive geometric scale.
- `toeplitz_psd_bridge_degree` joins this Toeplitz-square thread to
  HYP-2974's Fourier-Toeplitz PSD dual certificates.

Thus HYP-3064 belongs immediately after the HYP-3063 Moser-fibbinary
partial-cube sidecar, the HYP-3062 Roth-Minkowski fence, the HYP-3061
geometry-regime audit, and the HYP-3060 finalizer: Moser says not to use
sequence/cube/simplex motifs without native-transition, Theta, and gated
sidecars; Roth-Minkowski says not to use volume/height pressure without lattice
and exceptional-approximant sidecars; geometry-regime says not to confuse axes;
Desargues/Beal says not to call a rectangle-free residual structureless; and
Toeplitz-square says not to call a four-witness residual proved until the
scale/strictness gate survives.

## Tournament Analysis

Vertices are proof carriers and sidecar gates, not runners and not curve
points.  Candidate vertices considered and rejected as first-level objects:
curve points, square corners, raw diagonal pairs, raw D4 orbits, runners, and
raw Toeplitz moments.  The useful vertex set is proof obligations because it
preserves strict LRC witness status and named residual handoff.

Pairwise observable:

```text
predicate retention
noncollapse protection
reconstructibility
dual certificate value
route handoff
```

The carrier tournament is transitive:

```text
labelled_packet_sheaf >
toeplitz_square_configuration_space >
positive_scale_gate >
midpoint_balance_residue >
diagonal_equal_radius_residue >
quarter_turn_orthogonality_residue >
cyclic_order_D4_orbit >
floer_spectral_invariant_lane >
integration_sign_pattern_lane >
fourier_toeplitz_PSD_dual_bridge >
raw_square_peg_analogy
```

Fingerprint:

```text
score_hist={0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1}
directed_3cycles=0
scc_sizes=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
hamiltonian_path_count=1
```

## Next Pull

Add the Toeplitz-square fields after HYP-3060's Desargues/Beal finalizer fields
in the HYP-2963 packet ledger:

```text
desargues_girth6_residue
beal_common_owner_gate
toeplitz_square_scale_gate
midpoint_balance_residue
diagonal_equal_radius_residue
quarter_turn_residue
ordered_quad_collapse_mode
d4_orbit_word
toeplitz_psd_bridge_degree
```

A final residual should then be routed by the first nonzero or nonclosed gate:

- zero scale: boundary/AP-GW/degenerate witness debt, not a strict proof;
- midpoint/equal-radius/quarter-turn failure: Haar/rectangle/hourglass repair;
- D4/cyclic-order ambiguity: observer-cut orbit or value-origin repair;
- PSD bridge failure: Fejer/Toeplitz interval certificate;
- all gates survive: named family descent or THM-572/F7 state-lift debt.
