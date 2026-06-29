---
id: HYP-3494
title: LRC14 history-mined missed carriers for random031 seam cobordism
status: SYNTHESIS / repo-history carrier mining plus local PGF; not an LRC14 proof
source: codex-2026-06-29 history-mining pass after HYP-3486, HYP-3485, HYP-3490 Lean firewall ledger, HYP-3493 seam-sheaf scaffold, HYP-3510/HYP-3511 seam packet, and HYP-3480 Lean ledger
tangent: T1454
technique: LTI-454
tournament_technique: LTT-354
script: 04-computation/lrc14_history_missed_carriers_codex_20260629.py
result: 05-knowledge/results/lrc14_history_missed_carriers_codex_20260629.out
reflection: 07-reflections/lrc14-history-mined-missed-carriers-codex-20260629.md
related:
  - HYP-3511
  - HYP-3510
  - HYP-3490
  - HYP-3493
  - HYP-3486
  - HYP-3485
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3451
  - HYP-3437
  - HYP-3428
  - HYP-3422
  - HYP-3402
  - HYP-3300
  - HYP-3243
  - HYP-3140
  - HYP-3034
  - HYP-3023
  - THM-523
  - OPEN-Q-108
---

# HYP-3494: LRC14 History-Mined Missed Carriers

## Claim

The repo-history pass suggests a sharper proof interface for
`random_covering_031`:

```text
private-label firewall side       phase-flow seam-complement side
        HYP-3490          <->              HYP-3486
```

The bridge should be treated as a seam cobordism with sidecars, not as another
wall, projection-current, or scalar phase-count problem.  The older missed
carriers that now look most useful are:

After the HYP-3490 Lean update, the firewall side is also a formal ledger
interface: `TournamentH7.LRCPrivateLabelFirewall` records the finite
private-row count, the unique random031 private/hard overlap, and the dispatch
arithmetic.  HYP-3494 should therefore treat the remaining random031 work as
the phase/topology/PGF side needed to meet that formal firewall ledger.

The incoming HYP-3493 seam-sheaf scaffold is the first right-boundary
execution of this plan: it keeps the same `79` legal components, finds no
mixed/debt stalks, and isolates the single pure-bypass stalk as
owner-boundary debt rather than mirror-closure or projection-current debt.
Thus HYP-3494's quotient-price matrix should import HYP-3493 as the first
relative-H1/stalk sidecar and then focus the hard proof on owner-boundary
persistence.

The later HYP-3510/HYP-3511 seam packet further tightens the right boundary.
HYP-3510 gives a coarse branch-ordered seam-complement carrier with one
phase component after hard-seam deletion, while HYP-3511 brackets the `40`
free-hole cells as `10` ordinary-bracketed packets plus `2` ordinary-bracketed
doublets.  So HYP-3494 should now price three local terminal lemmas rather
than a vague free-hole debt: rank-2 routing, free-hole bracketing, and pure
bypass owner-boundary persistence.

```text
two-adic branch loss ledger
observability matrix / Morse descent
endpoint-owner current
oriented topes / cocircuits
overlap-tax Menger cuts
relative H1 owner persistence
magnitude-cocycle zipper
component-cover conductance
local fiber-PGF moments
```

The first new concrete object is the local random031 escape-sheet polynomial.
On occupied `u=2t mod 1` fibers after deleting the hard seam,

```text
P_escape(y) = 24 + 226 y + 8 y^2
mean_escape_sheets = 121/129.
```

This is the HYP-3140 sheet-count idea localized to the two-sheet
seam-complement cylinder.  It is stricter than the raw `282` witness count and
more proof-facing than only saying `242` cells hit rank-2 gates.

## Exact Readout

Script:

```text
04-computation/lrc14_history_missed_carriers_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_history_missed_carriers_codex_20260629.out
```

The script imports HYP-3486 and recomputes the fiber data:

```text
witness_cells=282
occupied_fibers=258
terminal_sheet_hist={1:234, 2:24}
rank2_escape_pgf={0:24, 1:226, 2:8}
mean_rank2_escape_sheets_per_occupied_fiber=121/129
class_signature_top=[
  (ordinary, 200),
  (free_hole, 22),
  (free_hole+ordinary, 14),
  (bypass, 12),
  (ordinary+ordinary, 8),
  (free_hole+free_hole, 2)
]
legal_vs_vertical_components=79 legal, 69 vertical-glued, 7 vertical-mixed
```

Interpretation:

```text
24 fibers have no rank-2 exit sheets, but are free-hole terminal packets;
226 fibers have exactly one rank-2 exit sheet;
8 fibers have two rank-2 exit sheets.
```

So the HYP-3486 trichotomy can be restated as a finite coefficient theorem:
rank-2 exit sheets plus HYP-3511-bracketed free-hole sheets cover the occupied
seam-complement phase fibers, while the pure bypass is one of the one-exit
components rather than a second wall.  HYP-3510 supplies the coarse connected
phase-carrier witness that this split is a refinement of a single
seam-complement transport object, not a disconnected set of holes.

## Missed-Carrier Ranking

The history-mining script scored current anchors separately, then ranked older
carriers by random031 predicate retention, sidecar cost, and executable
follow-up clarity.  The top older carriers were:

```text
HYP-3480  formal singleton-current contrast
HYP-3428  two-adic descent loss ledger
HYP-3300  observability matrix and Morse descent
HYP-3402  endpoint-owner current / tropical wall
HYP-3243  oriented topes/cocircuits atlas
HYP-3437  overlap-tax Menger cut
HYP-3034  owner-essential relative H1
HYP-3422  two-adic relocation identity
HYP-3023  magnitude-cocycle zipper guard
HYP-3451  component-cover conductance router
HYP-3140  fiber-PGF first-moment certificate
```

The ordering should not be read as theorem confidence.  It is a work queue:
which old carriers become newly relevant after HYP-3486 proves that vertical
half-turn is an address projection, not a legal topology quotient.

## New Experiments

1. **Local PGF escape moment.**  Prove that `24 + 226y + 8y^2` is the legal
   random031 seam-complement escape polynomial after hard-seam deletion, with
   the zero-exit coefficient delegated to the HYP-3511 free-hole bracket lemma.
2. **Relative H1 exit.**  Build
   `H1(seam complement, rank2 exits union free-hole boundary)` and test
   whether the pure bypass is the only nontrivial hard-component flow class
   after HYP-3510 connectedness and HYP-3511 bracketing are attached.
3. **Quotient-price matrix.**  Instantiate HYP-3300 with columns:
   `u_index`, branch, mirror mate, cell class, endpoint rank, owner word,
   private firewall, HYP-3511 bracket type, HYP-3510 carrier component, and
   vertical-halfturn sidecar.
4. **Owner-current cobordism.**  Prove seam-only owners
   `(45,147,169,173)` are boundary charge while bypass owners `(23,93,113)`
   are flow charge, or emit named owner/two-adic/SPEC debt.
5. **Formal dispatch join.**  Extend the HYP-3480 Lean dispatch ledger so the
   private-firewall rows split into six singleton-current terminals plus one
   random031 seam-cobordism terminal.

## Proof Pull

The current random031 terminal theorem should not try to rediscover a
projection current.  HYP-3490 says that carrier is firewalled.  The proof
interface should instead be:

```text
private-label firewall
  + relative seam-sheaf stalk audit
  + branch-ordered seam-complement connectivity
  + free-hole bracket atlas
  + legal mirror-run fiber graph
  + local escape-sheet PGF
  + relative H1 / owner-current sidecar
  + quotient-price matrix
  => rank-2 route discharge
     or free-hole bracket discharge
     or pure-bypass owner-boundary discharge
     or named owner/two-adic/SPEC debt.
```

This also explains why `n+2` versus `n*2` felt slippery.  The `n*2` side gives
the address coordinate `u=2t`; the `n+2` side gives the owner-boundary seam.
They are not competing recursions.  They are the two ends of the seam
cobordism, and HYP-3428 is the loss ledger that prevents collapsing one into
the other.

## Tournament Analysis

Vertices are history-mined proof carriers, not runners or raw gates:

```text
formal_singleton_current_contrast
two_adic_descent_loss_ledger
observability_morse_matrix
endpoint_owner_current
oriented_tope_cocircuit_atlas
overlap_tax_menger_cut
relative_h1_owner_persistence
two_adic_relocation_identity
magnitude_cocycle_zipper
component_cover_conductance
fiber_pgf_moment
```

Pairwise observable:

```text
random031 predicate retention
+ quotient sidecar cost
+ executable follow-up clarity
+ compatibility with HYP-3486 vertical guardrail
```

Fingerprint from the script:

```text
score_hist={81:1,94:1,99:1,105:1,106:1,108:1,114:1,115:1,116:1,119:1,148:1}
directed_3cycles=0
hamiltonian_path=
  HYP-3480 -> HYP-3428 -> HYP-3300 -> HYP-3402 -> HYP-3243 ->
  HYP-3437 -> HYP-3034 -> HYP-3422 -> HYP-3023 -> HYP-3451 -> HYP-3140
```

## Assumption Challenge

Candidate vertices considered: runners, gaps, hard gates, witness cells,
u-fibers, branch sheets, owner labels, dead islands, punctures, relative
cycles, component cuts, PGF coefficients, observability columns, Lean dispatch
packets, and proof obligations.

Chosen vertices: proof carriers plus the local fiber polynomial.

Preserved predicate:

```text
random031 terminal discharge after deleting the hard seam.
```

Destroyed coordinates are priced explicitly: vertical half-turn sheet identity,
free-hole status, bypass purity, endpoint-rank exit, private-label firewall,
owner-boundary charge, and formal dispatch status.
