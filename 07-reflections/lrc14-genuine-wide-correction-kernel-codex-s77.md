# LRC14 Genuine-Wide Correction Kernel (codex S77)

## What changed

The HYP-2805 adjacent-doublet correction is now split from the generalized-doublet proof obligations.

The exact kpswf9 script `lrc14_genuine_wide_true_maximizer_kpswf9.py` was rerun because the stored output was NUL-interleaved.  The clean result confirms the reported adjacent-doublet table:

- k=8: `348487/1681680 < 2243/5880`
- k=9: `321/980 < 1979/4004`, with non-primitive bounded base
- k=10: `265/588 < 55/91`, with non-primitive bounded base
- k=11: `5617/10780 < 66/91`
- k=12: `347/588 < 6/7`

The smallest reported cap margin is the k=10 row:

```text
55/91 - 265/588 = 783/5096 < 4/25.
```

So the actual cap predicate remains true on the reported table, but the stronger robust `0.16` margin predicate is false.  This is now formalized in `TournamentH7.LRCGenuineWideCorrection`.

## What this does not prove

This Lean module is only an arithmetic import boundary for the completed adjacent-doublet correction table.  It does not prove:

- the HYP-2807 generalized-doublet reduction,
- the frozen-room inequality `Phi_frozen(B,g)<cap`,
- the HYP-2808 Mordell-Tornheim R-tail bound,
- the finite-window closure for all `(B,g,M)`.

The newly pulled `lrc14_genwide_finite_window_claudeopus_0622.py` currently uses a naive exact `p0_fast` loop over all bounded bases, gaps, and positions.  A direct run stayed CPU-bound for several minutes before any k=10 row output, so it should not yet be cited as a completed certificate.  It needs a THM-563-style endpoint tiling/reuse engine, stronger pruning before exact integration, or a smaller certified atlas.

## Tournament carrier

The useful tournament vertices here are proof carriers rather than runners:

```text
adjacent_exact_table
> Lean_arithmetic_kernel
> primitive(FULL E)_guardrail
> generalized_doublet_window_runner
> frozen_room_bound
> Tornheim_R_tail
> wide_gK8_direct_route
> raw_runner_tournament
```

This preserves the cap-vs-robust-margin distinction and the base-primitivity guardrail.  It intentionally destroys the row geometry, so it cannot replace the generalized-doublet finite-window proof.
