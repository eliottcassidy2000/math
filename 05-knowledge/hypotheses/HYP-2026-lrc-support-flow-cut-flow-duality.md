---
id: HYP-2026
status: OPEN
source: codex-2026-06-01-S540
related:
  - HYP-1794
  - HYP-1802
  - HYP-1811
  - HYP-1812
  - HYP-2001
  - HYP-2007
  - HYP-2016
  - HYP-2017
  - HYP-2022
  - HYP-2023
  - HYP-2024
  - HYP-2025
  - THM-379
  - THM-380
---

# HYP-2026: LRC is trapped between support-flow cancellation and a zero-cut cover flow

**Claim.** HYP-2025's nowhere-zero flow dictionary should be expanded from the
full-support resonance to the entire LRC measure.  The exact lonely-time
measure is a support-decomposed flow enumerator:

```text
|SAFE| = sum_{m : sum_i m_i v_i = 0} prod_i ghat_n(m_i)
       = sum over supports S of NZ integer flows on the speed sub-dipole G_v[S].
```

On the dual cover side, each danger interval is a directed arc on the time
circle.  The coverage count on the endpoint arrangement is a nonnegative flow
on the directed cell cycle:

```text
positive lonely interval  <=>  a zero-flow cut in the cover-flow cycle.
wall-only AP extremal     <=>  nowhere-zero cover flow on open cells,
                               but empty strict endpoint-protection core.
```

Thus a genuine counterexample would have to satisfy two simultaneous directed
graph conditions:

```text
1. the support-flow enumerator cancels the positive main term completely;
2. the cover-flow side has a nonpeeling nowhere-zero endpoint-protection core,
   with best margin still below 1/n.
```

This is much more rigid than merely having many full-support resonances.

**Evidence.** `lrc_flow_cut_support_s540.py` audits AP wall cases and open
sample cases at `n=4,5,6,7`, plus the `n=14` AP wall.  It computes:

1. exact open safe measure and closed best margin;
2. the danger-arc cover-flow profile on the endpoint cell cycle;
3. the strict endpoint-protection core after peeling unprotected interval
   endpoints;
4. truncated Fourier support-flow ledgers by support size;
5. modular `Z_{n*}` nowhere-zero sub-dipole profiles;
6. a support-layer tournament fingerprint.

The repeated split is:

```text
AP wall cases:
  open |SAFE| = 0, best_margin = 1/n, min cover-flow = 1,
  zero open cells = 0, strict endpoint-protection core = 0.

Open cases:
  open |SAFE| > 0, min cover-flow = 0, zero cells > 0,
  strict endpoint-protection core = 0.
```

Selected output:

```text
n4 AP:        open SAFE 0,     best 1/4,  min_coverage 1, zero_cells 0, core 0
n4 (1,2,4):  open SAFE 1/8,   best 5/16, min_coverage 0, zero_cells 2, core 0
n6 AP:        open SAFE 0,     best 1/6,  min_coverage 1, zero_cells 0, core 0
n6 Z3 bridge: open SAFE 5/54,  best 7/36, min_coverage 0, zero_cells 8, core 0
n14 AP:       open SAFE 0,     best 1/14, min_coverage 1, zero_cells 0, core 0
```

The Fourier side also separates the issue.  For `n=14` AP with truncation
`M=3`, support layers already show large alternating flow cancellation:

```text
main term +0.13480057
nonzero support contribution -0.05988926  (truncated)
largest layers: support 4..8 alternate + - + - +
```

The truncation has not converged to zero, but it reveals the demanded shape:
the AP wall is not killed by a single full-support resonance.  It is a
coordinated alternating cancellation across many sub-dipole NZ-flow layers.

**Interpretation.** There are now two graph-dual ways to try to prove LRC:

```text
flow enumerator route:
  prove the sub-dipole NZ-flow layers cannot cancel the positive main term
  unless the cover is the regular-polygon wall class.

cut-flow route:
  prove any nowhere-zero open-cell cover flow has an empty or peelable
  endpoint-protection core, hence leaves either a closed wall witness or a
  positive zero-flow cut.
```

HYP-2025 says full-support parity is a nowhere-zero flow existence question.
HYP-2026 says the real proof target is not the existence of flows, but the
failure of those flows to assemble into both a total cancellation and a
nonpeeling cover core.

**Predictions.**

1. For AP `1..n-1`, the exact support-flow ledger is a finite/infinite
   alternating cancellation whose cover-flow core always peels to empty; this
   is the graph-flow signature of a wall extremal rather than a counterexample.
2. Any positive-measure LRC witness appears as a zero edge in the directed
   cover-cell cycle before the endpoint-protection core can become nonempty.
3. A true counterexample would force a nonempty SCC in the endpoint-protection
   digraph from THM-379/THM-380 and simultaneously a complete cancellation in
   the support-flow enumerator.
4. Hard rows such as `n=14` and `n=18` should be attacked by intersecting the
   support-flow cancellation ledger with endpoint/carry SCC peeling, not by
   full-support parity alone.

**Next tests.**

1. Replace the truncated Fourier support ledger with an exact residue-class
   or interval-kernel summation where possible.
2. Build the endpoint-protection digraph explicitly and report SCCs, feedback
   arcs, and cut certificates instead of only the peeled core size.
3. Run the cut-flow core on the hard `n=14` and `n=18` ladder families from
   HYP-2001/HYP-1992.
4. Search for artificial full-cover interval systems with nonempty cores, then
   test which arithmetic constraints of LRC forbid them.

**Files.** `04-computation/lrc_flow_cut_support_s540.py`;
`05-knowledge/results/lrc_flow_cut_support_s540.out`;
`07-reflections/lrc-support-flow-cut-flow-duality-s540.md`.
