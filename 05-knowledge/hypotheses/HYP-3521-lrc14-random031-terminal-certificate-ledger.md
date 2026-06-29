---
id: HYP-3521
title: LRC14 random031 terminal certificate ledger
status: EVIDENCE / finite terminal-dispatch ledger; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3511, HYP-3510, HYP-3494, HYP-3493, HYP-3490, HYP-3486, HYP-3485, HYP-3484, HYP-3483, HYP-3482, HYP-3481, HYP-3480, HYP-3477, HYP-3460, and HYP-3455
tangent: T1521
technique: LTI-521
tournament_technique: LTT-421
script: 04-computation/lrc14_random031_terminal_certificate_ledger_codex_20260629.py
result: 05-knowledge/results/lrc14_random031_terminal_certificate_ledger_codex_20260629.out
reflection: 07-reflections/lrc14-random031-terminal-certificate-ledger-codex-20260629.md
related:
  - HYP-3511
  - HYP-3510
  - HYP-3494
  - HYP-3493
  - HYP-3490
  - HYP-3486
  - HYP-3485
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3477
  - HYP-3460
  - HYP-3455
  - THM-523
  - OPEN-Q-108
---

# HYP-3521: LRC14 Random031 Terminal Certificate Ledger

## Claim

The remaining `random_covering_031` packet should be formalized as a terminal
certificate ledger, not as another search for a projection current.

The new ledger joins four recent facts:

```text
HYP-3486: legal mirror-run fiber graph
HYP-3511: free-hole bracket / doublet classification
HYP-3510: coarse seam-complement connected phase carrier
HYP-3490: private-label firewall blocking projection-current deletion
```

The useful theorem target is the finite packet identity:

```text
282 phase cells
= 230 ordinary endpoint-rank-2 route cells
+  40 bracketed free-hole cells
+  12 pure bypass owner-boundary cells.
```

The older `242` count remains correct but should be named carefully:

```text
242 gate-routed cells = 230 ordinary route cells + 12 bypass route cells.
```

That naming split matters.  The bypass cells also hit endpoint-rank-2 gates,
but they are the only gate-routed hard-component flow and carry the remaining
owner-boundary obligation.

## Exact Readout

Computed by:

```text
04-computation/lrc14_random031_terminal_certificate_ledger_codex_20260629.py
05-knowledge/results/lrc14_random031_terminal_certificate_ledger_codex_20260629.out
```

Cell partition:

```text
witness_cells=282
cell_class_counts={bypass:12, free_hole:40, ordinary:230}
gate_routed_cells=242
ordinary_rank2_hit_cells=230
bypass_rank2_hit_cells=12
free_hole_cells=40
endpoint_rank_hist_on_hit_cells={2:242}
```

Component and certificate partition:

```text
legal_mirror_components=79
component_type_hist={bypass:1, free_hole:14, ordinary:64}
component_terminal_hist={
  ordinary_rank2_route_component:64,
  free_hole_single_bracket_packet:10,
  free_hole_doublet_packet:4,
  pure_bypass_owner_boundary_component:1
}
cell_terminal_hist={
  ordinary_rank2_route_component:230,
  free_hole_single_bracket_packet:26,
  free_hole_doublet_packet:14,
  pure_bypass_owner_boundary_component:12
}
terminal_certificate_count_after_doublet_collapse=77
```

Thus:

```text
79 legal components
= 64 ordinary components
+ 14 free-hole components
+  1 bypass component

77 terminal certificates
= 64 ordinary route certificates
+ 10 free-hole single certificates
+  2 free-hole doublet certificates
+  1 bypass owner-boundary certificate.
```

The free-hole part collapses exactly as HYP-3511 predicted:

```text
free_single_cells=26
free_doublet_cells=14
doublet_cell_summary={(2,3):4, (8,13):10}
```

## Fiber PGF Sidecar

The HYP-3494 fiber PGF becomes more interpretable after the certificate split:

```text
occupied_fibers=258
rank2_escape_pgf={0:24, 1:226, 2:8}
mean_rank2_escape_sheets=121/129
zero_exit_signature_hist={
  (free_hole):22,
  (free_hole, free_hole):2
}
one_exit_signature_top=[
  ((ordinary),200),
  ((free_hole,ordinary),14),
  ((bypass),12)
]
```

So the zero-exit coefficient is exactly pure free-hole fibers.  The mixed
`free_hole+ordinary` fibers are not failures: each carries one endpoint-rank-2
sheet and one bracketed no-gate sheet.  This gives a cleaner PGF theorem
target than a raw mean or count.

## Pure Bypass Boundary Packet

The unique bypass certificate is:

```text
component_index=1
size=12
hit_components=(43,54)
owners=(23,93,113)
seam_debt=(45,147,169,173)
```

This is the only terminal certificate that is both gate-routed and not an
ordinary route component.  It should be treated as the owner-boundary lemma:
prove this packet discharges the seam boundary, or name the exact
owner-current/two-adic/signed-SPEC sidecar still missing.

## Proof Pull

The remaining random031 terminal theorem can now be split into five finite
lemmas:

1. **Ordinary-route lemma.** Every ordinary legal mirror-run component has
   endpoint-rank-2 hits and exits the seam complement.
2. **Free-hole bracket lemma.** All `40` no-gate cells are terminal local
   beads: `10` ordinary-bracketed singles plus `2` same-branch doublets.
3. **Bypass owner-boundary lemma.** The only nonordinary gate-routed
   hard-component packet is the `12`-cell bypass with owners `(23,93,113)` and
   seam debt `(45,147,169,173)`.
4. **Firewall compatibility lemma.** HYP-3490 forbids replacing the previous
   three lemmas by adjacent-label projection-current deletion.
5. **Quotient guardrail lemma.** Vertical half-turn is an address projection;
   it is not allowed to identify ordinary, free-hole, and bypass terminal
   classes.

If these five statements are proved, the random031 part of the current LRC14
frontier has a compact terminal packet theorem.

## New Angles Exposed

### 1. Count Ambiguity Becomes A Proof Split

The phrase "`242` rank-2 routed cells" is convenient but hides the only
remaining hard-component owner-boundary packet.  HYP-3521 replaces it with:

```text
230 ordinary rank-2 route cells
12 bypass rank-2 route cells
```

This prevents the bypass from disappearing inside the easy endpoint-rank
lemma.

### 2. Free-Hole Debt Is Not All Zero-Exit Debt

Only `26` of the `40` free-hole cells lie over zero-exit fibers.  The remaining
`14` live in mixed `free_hole+ordinary` fibers with one ordinary escape sheet.
This suggests a small "one escaped sheet brackets one missing sheet" lemma.

### 3. Doublet Collapse Is A Formalization Shortcut

The HYP-3511 half-open packets reduce `14` free-hole components to `12`
free-hole certificates.  That is exactly the kind of compression a Lean
terminal ledger should use: prove the doublet rule once, instantiate it twice.

### 4. Bypass Is The Only Owner-Boundary Certificate

After HYP-3521, owner-boundary persistence is no longer distributed over the
whole random031 packet.  It is concentrated in one certificate with flow owners
`(23,93,113)` and seam-only debt `(45,147,169,173)`.

## Tournament Analysis

Vertices are terminal proof obligations, not runners, raw arcs, or scalar cell
counts:

```text
terminal_certificate_ledger
ordinary_rank2_route_lemma
free_hole_bracket_lemma
pure_bypass_owner_boundary_lemma
private_firewall_negative_certificate
fiber_pgf_zero_exit_split
vertical_halfturn_guardrail
raw_count_shadow
```

Pairwise observable: which carrier preserves the terminal random031 discharge
predicate with the least forgotten sidecar.  The switch orients toward higher
retained certificate payload; ties follow the Hamiltonian path above.

Fingerprint:

```text
score_hist={10:1,66:1,73:1,79:1,84:1,88:1,91:1,100:1}
directed_3cycles=0
sccs=8 singleton SCCs
hamiltonian_path=
  terminal_certificate_ledger
  -> ordinary_rank2_route_lemma
  -> free_hole_bracket_lemma
  -> pure_bypass_owner_boundary_lemma
  -> private_firewall_negative_certificate
  -> fiber_pgf_zero_exit_split
  -> vertical_halfturn_guardrail
  -> raw_count_shadow
```

## Assumption Challenge

Candidate vertices considered: runners, gaps, hard gates, witness cells,
`u`-fibers, branch sheets, free-hole packets, doublet clusters, owner words,
endpoint ranks, vertical pairs, and proof obligations.

Chosen vertices are terminal proof obligations.  This quotient preserves the
random031 terminal predicate because it retains cell class, endpoint rank,
free-hole bracket type, bypass owner debt, private-firewall status, and the
vertical guardrail.
