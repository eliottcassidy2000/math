---
id: HYP-3253
title: LRC14 contact-holonomy curvature sheaf
status: EVIDENCE / exact bounded-bank scout plus proof-frame synthesis; not an LRC14 proof
source: codex-2026-06-28
tangent: T1347
technique: LTI-347
tournament_technique: LTT-247
script: 04-computation/lrc14_contact_holonomy_curvature_codex_20260628.py
result: 05-knowledge/results/lrc14_contact_holonomy_curvature_codex_20260628.out
reflection: 07-reflections/lrc14-contact-holonomy-curvature-sheaf-codex-20260628.md
related:
  - HYP-3252
  - HYP-3251
  - HYP-3250
  - HYP-3249
  - HYP-3248
  - HYP-3247
  - HYP-3246
  - HYP-3245
  - HYP-3244
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3239
  - HYP-3228
  - HYP-3204
  - OPEN-Q-108
---

# HYP-3253: LRC14 Contact-Holonomy Curvature Sheaf

## Claim

The current LRC14 frames can be organized as one certificate sheaf:

```text
global section / index     = HYP-3246 index-degree frame
calibrated boundary        = HYP-3246/HYP-3245 unit equioscillation
local quotient curvature   = HYP-3247 shell-lag commutator
connection repair          = zeta_7 contact holonomy
frontier split             = HYP-3250 finite tight-locus + uniform margin
```

The new bounded-bank evidence is that HYP-3247's ordered contact sidecar has a
compact cyclotomic connection coordinate.  Over the exact
`ordinary lag profile + residue histogram mod 7` fibers, the first contact
moment

```text
contact_holonomy(E) = sum_{j in contact_support(E)} zeta_7^j
```

repairs every HYP-3228 shell-magic ambiguity.  It is not a global terminal
quotient: empty and full contact supports have the same zero holonomy.  Its
proof-facing role is narrower and more useful: it is the holonomy coordinate
for the lag/residue-to-shell quotient square.

## Exact Scout

The script enumerates the same exact anchored bounded k=8 bank as HYP-3247:

```text
E = {0} union A, A subset {1,...,14}, |A|=7
rows = 3432
```

The base quotient is:

```text
Q(E) = (ordinary support-autocorrelation, residue histogram mod 7).
```

Curvature means that HYP-3228 shell magic is nonconstant on a `Q`-fiber.
Exact result:

```text
lag+residue mixed_shell_fibers = 62
mixed_rows = 124
mixed_fiber_size_hist = {2: 62}
```

Repair table:

```text
support_size             mixed_shell_fibers = 62
position_sum_mod7        mixed_shell_fibers =  2
position_power_sums      mixed_shell_fibers =  2
minmax_position          mixed_shell_fibers = 14
ordered_gap_values       mixed_shell_fibers =  9
gap_multiset             mixed_shell_fibers = 62
contact_support          mixed_shell_fibers =  0
zeta7_contact_holonomy   mixed_shell_fibers =  0
```

So the first cyclotomic moment is as effective as the full ordered support for
the exact HYP-3247 interface, while several plausible scalar repairs fail.

The first scalar-position failure is:

```text
E=(0,1,3,5,9,10,11,13)
  magic=126097/90090
  contact_support=(1,2,3,6)
  holonomy=(-1,0,0,0,-1,-1)

E=(0,2,3,4,8,10,12,13)
  magic=5989/5460
  contact_support=(0,3,4,5)
  holonomy=(1,0,0,1,1,1)
```

They have the same scalar position sum modulo `7`, but different holonomy and
different shell magic.

## Proof-Frontier Use

This reframes the proof route from "add one more sidecar" to "show curvature
is exact after the right connection is installed."  The theorem-facing target:

```text
For every legal primitive packet, the lag/residue quotient either
  has zero shell curvature,
  has curvature killed by zeta_7 contact holonomy,
  lifts to an endpoint arrangement cell / finite chamber,
  or emits named residual debt.
```

That target is compatible with HYP-3246: a nonzero global index can be viewed
as holonomy of the certificate sheaf, while the local contact holonomy prevents
lag-space scalarization from erasing the endpoint coordinate that carries the
index.  It also absorbs HYP-3249's P2 warning: the naive odd map may produce a
runner at the observer rather than a lonely point, so the proof map must keep
cover-hole / endpoint-cell data rather than only a signed nearest-runner
coordinate.  HYP-3253 supplies a local version of that guardrail: the
lag/residue quotient becomes legal only after the endpoint/contact holonomy is
kept or an endpoint chamber names the missing coordinate.  It also tightens
HYP-3245/HYP-3247: the shell-lag commutator is not just "positional"; on the
bounded interface it is cyclotomic positional.

Post-rebase integration with HYP-3250/HYP-3251/HYP-3252 makes the role more
specific.  HYP-3250 says the proof frontier should split into finite
tight-locus rigidity plus a uniform-margin floor.  HYP-3251 and HYP-3252 say
the index/equioscillation count describes the AP saddle but does not by itself
perform the S-dependent proof.  Contact holonomy is therefore not a new global
obstruction; it is a local sidecar for deciding whether a lag/residue packet
belongs to the tight endpoint chamber, descends to the floor/margin side, or
has named residual curvature debt.  The new reflection
`the-index-theorem-describes-the-floor-proves.md` points the next bridge toward
`Q(sqrt(-7))`; HYP-3253's `zeta_7` holonomy is the finite local coordinate
that should be tested before trying to reorganize the floor in that field.

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
oriented-matroid topes/cocircuits, endpoint cells, and proof obligations.

Chosen quotient vertices: lag/residue fibers lifted by a `zeta_7` contact
holonomy.

Preserved LRC predicate: the HYP-3228 shell-magic value as it feeds the
HYP-3245/HYP-3246 proof route.

Destroyed information: the base quotient forgets endpoint placement of
contact defects; holonomy alone forgets the empty/full distinction in rare
global fibers.  Therefore holonomy is a connection coordinate, not permission
to discard endpoint cells.

Challenged assumption: HYP-3247's ordered endpoint repair must remain a raw
support set.  Exact evidence says the first cyclotomic contact moment is enough
over the lag/residue connection.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
raw_runner_vertices
lag_plus_residue_histogram
gap_multiset
position_power_sums
ordered_contact_support
zeta7_contact_holonomy
endpoint_arrangement_cell
index_degree_sheaf
```

Pairwise observable: which carrier kills lag/residue-to-shell curvature while
retaining the LRC predicate.  Switch/gauge: `A -> B` iff `B` preserves strictly
more curvature/index payload with a named destroyed coordinate.

Bounded-bank fingerprint:

```text
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1]
edge_flips = 0
hamiltonian_path_count = 1
tie_hamiltonian_path =
  index_degree_sheaf
  -> endpoint_arrangement_cell
  -> zeta7_contact_holonomy
  -> ordered_contact_support
  -> position_power_sums
  -> gap_multiset
  -> lag_plus_residue_histogram
  -> raw_runner_vertices
```

## Status

This is evidence and a proof frame, not an LRC14 proof.  The bounded exact
result gives a compact local curvature repair.  The open theorem is to prove
that the same connection either kills all primitive residual curvature,
descends through endpoint/finite-chamber carriers, or produces a named
residual class compatible with the HYP-3250 finite-tight-locus/uniform-margin
split, the HYP-3251/HYP-3252 verdict that the index is descriptive rather than
S-dependent, and HYP-3249's forcing-gap guardrail.
