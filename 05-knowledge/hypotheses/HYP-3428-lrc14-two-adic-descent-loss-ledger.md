---
id: HYP-3428
title: LRC14 two-adic descent loss ledger
status: SYNTHESIS / exact AP-collar descent scout; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3418 and HYP-3417; renumbered after mainline HYP-3419-HYP-3427
tangent: T1389
technique: LTI-389
tournament_technique: LTT-289
script: 04-computation/lrc14_two_adic_descent_loss_ledger_codex_20260628.py
result: 05-knowledge/results/lrc14_two_adic_descent_loss_ledger_codex_20260628.out
reflection: 07-reflections/lrc14-two-adic-descent-loss-ledger-codex-20260628.md
related:
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3414
  - HYP-3410
  - HYP-3407
  - HYP-3406
  - HYP-3403
  - HYP-3311
  - HYP-3310
  - HYP-3265
  - HYP-3260
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3428: LRC14 Two-Adic Descent Loss Ledger

## Claim

HYP-3418 made the remaining LRC14 floor route 2-adic: even speeds, not the
apex-7 census, pull the witness away from the odd-speed maximum at `t=1/2`.
HYP-3428 sharpens the next proof obligation:

```text
A legal 2-adic descent quotient must carry a loss ledger recording
odd blockers, the halved even-child packet under u=2t, and any labelled
owner-current/even-hinge sidecar needed to keep theorem exits separated.
```

The exact scout is intentionally on the HYP-3410/HYP-3417 AP-collar
mixed-fiber substrate.  These rows are q-witness/non-covering proxies, not
arbitrary covering rows.  Their value is diagnostic: they show exactly why the
raw odd/coprime-to-14 witness is an illegal quotient and what information a
covering-floor version must retain.

After fetching the live mainline, this ledger should be read downstream of the
new covering-floor sequence HYP-3421 through HYP-3425.  In particular,
HYP-3422 already begins the genuine covering-packet lift by splitting
`S = O union 2E` and relocating from the dead half-witness to two odd branches,
while HYP-3425 sharpens that into a two-color Helly/interval-piercing target.
HYP-3428 contributes the complementary AP-collar loss audit: it names exactly
which descent information is unsafe to forget before invoking those covering
relocation lemmas.

After the second fetch, mainline added HYP-3426 one-branch mirror support and
HYP-3427 two-branch wall signatures.  Those are now the sharper proof carriers:
mirror first, wall words second, loss ledger third.  HYP-3428 is therefore the
exception router for mirror/wall failures, not the main interval-piercing
engine.

## Exact Readout

Script:

```text
04-computation/lrc14_two_adic_descent_loss_ledger_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_two_adic_descent_loss_ledger_codex_20260628.out
```

On the current substrate:

```text
rows_total=24
reconstructed=23
opaque_named_debt=1  # P10+GW
q_status_hist={'q-witness:14': 23}
half_witness_failure_rows=23/23
rows_with_even_binders=22/23
rows_where_even_child_carries_M_at_u=2t*=22/23
```

The fiber ledger is:

```text
height_leak_12_family:
  loss_hist={'even_child_with_odd_coupling': 2, 'opaque_named_row_debt': 1}
  selected owner-current cut=('5:g1',)

persistent_owner_leak_26_40_54_family:
  loss_hist={'even_child_with_odd_coupling': 11}
  selected owner-current cut=('1:g1',)

height_persistent_owner_leak_10_20_drop_add_family:
  loss_hist={'even_child_with_odd_coupling': 9,
             'odd_binder_after_even_shift': 1}
  selected owner-current cut=('2:g2', '11:g1', '13:g1')
```

Thus the naive half witness fails on every reconstructed AP-collar row, while
the even child still carries the exact maximin bottleneck in all but the single
odd-binder-after-even-shift case.  The one exceptional row is not noise; it is
the first named loss class the covering-floor descent must be able to discharge.

## Proof Target

The next covering-floor theorem should try to prove:

```text
For every covering-floor packet, either the halved even child inherits a
certified floor with controlled odd blockers, or the packet emits a bounded
owner-current/even-hinge certificate or named off-grid/state-lift debt.
```

Equivalently, in the newer HYP-3422/HYP-3425/HYP-3426/HYP-3427 language, prove
that the even-safe child interval set is not swallowed by the one-branch odd
bad core, then certify survivor windows by bounded wall words; record an
owner/current label whenever an AP-collar loss class reappears.

This is weaker and more usable than "odd subset lonely implies whole set
lonely."  It is also stricter than a v2 profile: the HYP-3410 fibers already
show mixed v2/survival signatures, while HYP-3417's hardest current contains
one even-cover label plus two binding labels:

```text
{2:g2, 11:g1, 13:g1}.
```

That is the current local avatar of the global `14 = 2 * 7` split: the proof
route is 2-adic, but the terminal sidecar still has to remember labelled
binding owners.

## Tournament Analysis

Tournament vertices are proof carriers and descent ledgers, not runners.  The
pairwise observable is retained floor/descent payload minus forgotten descent
debt.  The exact fingerprint is:

```text
vertices=7
score_hist={-12: 1, 6: 1, 16: 1, 24: 1, 34: 1, 36: 1, 55: 1}
directed_3cycles=0
hamiltonian_path_count=1
priority_hamiltonian_path:
  two_adic_descent_loss_ledger
  owner_current_even_hinge
  even_child_induction_packet
  odd_half_witness_failure_gate
  q_witness_noncovering_exit
  apex7_census_offpath_guard
  raw_coprime_to_14_reduction
```

## Assumption Challenge

Considered vertices included runners, odd speeds, even speeds, v2 layers,
owner labels, halved child packets, q-witness exits, off-grid witnesses, and
proof obligations.  The chosen quotient preserves the covering-floor predicate
to be proved and destroys raw row order.  The challenged assumption is that
the coprime-to-14 or odd-speed subproblem can supply a witness without
recording how even speeds move the optimum away from `t=1/2`.
