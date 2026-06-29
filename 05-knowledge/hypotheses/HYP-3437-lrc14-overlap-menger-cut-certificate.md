---
id: HYP-3437
title: LRC14 overlap-tax Menger-cut certificate
status: EVIDENCE / exact interval-atom cut certificate stress; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3436/HYP-3435/HYP-3434/HYP-3429/HYP-3425
tangent: T1398
technique: LTI-398
tournament_technique: LTT-298
script: 04-computation/lrc14_overlap_menger_cut_certificate_codex_20260628.py
result: 05-knowledge/results/lrc14_overlap_menger_cut_certificate_codex_20260628.out
reflection: 07-reflections/lrc14-overlap-menger-cut-certificate-codex-20260628.md
related:
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3418
  - HYP-3415
  - THM-523
  - OPEN-Q-108
---

# HYP-3437: LRC14 Overlap-Tax Menger-Cut Certificate

## Claim

HYP-3434 identified the exact one-branch compression identity

```text
branch0 = naive_slack + overlap_tax,
naive_slack = |E_safe| - sum_o |E_safe cap B0_o|,
overlap_tax = sum_atoms max(multiplicity(atom)-1,0) length(atom).
```

HYP-3437 turns the overlap tax into a finite graph/cut certificate.  For each
row, atomize every `E_safe` component by all branch-0 odd-bad interval
endpoints.  Each atom records which odd blockers cover it.  A rescue cut is a
small odd-blocker subset whose own overlap tax is strictly larger than the
negative naive deficit.

The proof-facing theorem target is:

```text
Every primitive covering row either has naive_slack >= 0,
or has a bounded endpoint-labelled odd-blocker cut core
whose overlap tax beats -naive_slack.
```

This is one-branch and complementary to HYP-3436: HYP-3436 is reserved for the
two-color bad-core extractor `E_safe cap B0_odd cap B1_odd`, while HYP-3437
supplies the graph/cut lower bound needed to make HYP-3434's overlap-tax rescue
structural instead of scalar.

## Exact Readout

Script:

```text
04-computation/lrc14_overlap_menger_cut_certificate_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_overlap_menger_cut_certificate_codex_20260628.out
```

The audit uses the `150`-row HYP-3429/HYP-3434 endpoint-spine bank.

```text
rows_audited=150
branch0_positive=150/150
naive_slack_nonnegative=91/150
naive_slack_negative=59/150
negative_rows_with_rescue_cut=59/59
all_row_rescue_rank_hist={0:91, 2:7, 4:2, 5:48, 6:2}
negative_row_rescue_rank_hist={2:7, 4:2, 5:48, 6:2}
```

The maximum minimum rescue rank is `6`, occurring on the canonical tight row
and its duplicate:

```text
covering_AP_with_84 = (1,2,3,4,5,6,7,8,9,10,11,13,84)
branch0 = 563/105105
naive_slack = -18586/315315
overlap_tax = 4055/63063
tax/deficit = 1.090875
rescue_subset = (3,5,7,9,11,13)
rescue_margin = 563/105105
```

The worst scalar deficit row is not close to failure after the cut sidecar is
kept:

```text
random_covering_12
naive_slack = -941318888891/12839701214340
tax/deficit = 2.739021
rescue_subset = (27,47,63,119,137)
```

The smallest rescue margin in the bank is still positive:

```text
random_spine_02
rescue_rank = 4
rescue_subset = (17,29,87,91)
rescue_margin = 724590768379/1008622450475640
```

For the canonical `84m` tower in the audited range, `m=1` is the only rank-`6`
case; all listed larger `m` use the stable core `(5,7,9,11,13)` with rank `5`.
That aligns with HYP-3431: the canonical tower is the tight corridor-fence base
case, not generic high-rank chaos.

## Structural Interpretation

The one-branch odd cover has a simple atom algebra:

```text
multiplicity 0 atoms = branch0 survivor mass,
multiplicity 1 atoms = scalar single-cover mass,
multiplicity >=2 atoms = overlap-tax curvature.
```

Thus HYP-3434's "overlap tax" is not diffuse harmonic noise.  In every
negative-slack row scanned, it has a small labelled cut core.  The hard proof
obligation is now bounded and local:

```text
prove bounded cut rank / endpoint word,
or emit named two-adic, owner-current, exact-period, state-lift, or SPEC debt.
```

## Creative Route Translation

Schwarz-Christoffel enters as a faithful finite metaphor: interval endpoints
are prevertices of the branch arrangement.  The proof must keep those
prevertices and their endpoint-owner labels; it cannot replace them with a
continuous area scalar.

Bring-radical language marks monodromy debt: if a row requires a high-rank cut
core, the branch choice is the object to retain, not a reason to scalarize.

Barban-Davenport-Halberstam and Mertens can enter only after this exact finite
object exists, as mean-square or reciprocal-tail estimates over residue classes
of labelled cut cores with endpoint exceptions removed.

The AP six-touch equioscillation is the equality boundary/prevertex model.
The covering floor remains driven by the two-adic branch atoms, not by the
off-path 7-adic census alone.

## Candidate Lemma

For every primitive covering row `S=O union 2E`, after atomizing
`E_safe(1/14)` by branch-0 odd bad endpoints, either:

1. `|E_safe| - sum_o |E_safe cap B0_o| >= 0`;
2. a bounded odd-blocker subset `K` has
   `overlap_tax(K) > sum_o |E_safe cap B0_o| - |E_safe|`;
3. the row routes to HYP-3436's two-color bad-core extractor, HYP-3431's
   corridor-fence theorem, HYP-3429 endpoint-spine rank, HYP-3428 loss ledger,
   HYP-3427 wall word, owner-current labels, exact-period/state-lift debt, or
   HYP-3129 signed-SPEC.

The evidence suggests "bounded" may be `6` for the current bank, with the only
rank-`6` stress being the canonical AP-with-84 tight row.  A proof should try
to lower that bound for noncanonical rows or prove that the rank-`6` canonical
case is exactly the corridor-fence family already covered by HYP-3431.

## Tournament Analysis

Vertices are proof carriers and cut interfaces, not runners, arcs, or raw
rows.

```text
pairwise_observable =
  retained LRC predicate + overlap-tax payload + scalar-forgetting debt
switch_gauge =
  higher weighted proof-facing score; ties by declared carrier order
score_hist = {21:1, 49:1, 53:1, 54:1, 61:1, 65:1, 66:2}
directed_3cycles = 0
hamiltonian_path_count = 1
hamiltonian_path =
  atomic_interval_arrangement
  -> overlap_tax_menger_core
  -> endpoint_spine_cut_lift
  -> two_adic_branch_induction
  -> bdh_mean_square_sidecar
  -> schwarz_christoffel_prevertex_model
  -> bring_monodromy_branch_guard
  -> raw_harmonic_scalar
```

Assumption challenge: runners, odd blockers, even-half gates, `E_safe`
components, atomic cells, endpoint labels, wall-crossing events, residues,
cover arcs, Fourier modes, matroid circuits, graph cuts, and proof obligations
were considered.  The chosen vertices are atomic cells plus proof-carrier cuts.
This preserves the one-branch LRC predicate exactly because branch0 and
overlap tax are reconstructed from atom multiplicities; it destroys raw runner
order and most continuous geometry.
