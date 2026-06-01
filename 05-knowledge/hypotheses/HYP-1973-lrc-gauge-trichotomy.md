---
id: HYP-1973
status: OPEN
source: codex-2026-06-01-S506
related:
  - HYP-1967
  - HYP-1968
  - HYP-1969
  - THM-373
  - THM-380
---

# HYP-1973: LRC Tournament Analysis needs ranker, entropy, and pressure gauges

## Statement

No single arc-assignment rule should be expected to carry the whole LRC proof.
Useful Tournament Analysis of LRC splits into three gauge types:

```text
ranker gauges:
  scalar loneliness coordinates; H collapses.

phase/entropy gauges:
  cyclic global shape; H is a meaningful loneliness/spread metric.

pressure gauges:
  blocker-debt dependencies; score/SCC/source-sink structure is meaningful.
```

The current best clean H-loneliness gauge is the half-turn phase tournament.
The best endpoint-aware H perturbation is `safe_phase_gate`.  The best
pressure-shape candidate from S506 is `pressure_k2`.

## Evidence

`04-computation/lrc_tournament_gauge_lab_s506.py` audits twelve gauges on the
hard `n=14` and `n=18` LRC recursive rows.

The half-turn phase gauge has broad H range and sees the recursive plateau:

```text
n14 scales 7->14->28: H_ratio 2.0831 -> 1.6752 -> 1.6752
n18 scales 9->18->36: H_ratio 2.3981 -> 2.0965 -> 2.0965
```

The scalar rankers all collapse:

```text
distinct H values = 1,
cyclic triples = 0,
largest SCC = 1.
```

The endpoint-aware `safe_phase_gate` keeps strong spread correlation but can
collapse H too aggressively:

```text
H_ratio range 0.0204..2.3675,
spread correlation 0.922.
```

The pressure gauges have tiny H but useful dependency shape:

```text
pressure_k2: cyclic triples 0..40, largest SCC 1..18.
```

## Prediction

Future LRC frontier scans should report the composite metric

```text
(
  H_ratio(phase_half),
  H_ratio(safe_phase_gate),
  score_width(pressure_k2),
  largest_SCC(pressure_k2),
  endpoint_debt,
  gap * endpoint_debt
).
```

A real counterexample-shaped row should break all three layers at once:

```text
non-plateau phase H,
endpoint-safe hybrid alarm,
and nontrivial pressure SCC with surviving endpoint core.
```

Rows that only shrink scalar gap while phase H plateaus and endpoint debt grows
are denominator-depth recursions, not new global loneliness shapes.

## Proof Program

1. Prove that every scalar vertex-score gauge yields a transitive tournament
   after tie breaking, hence H cannot be the right statistic there.
2. Prove the S505 phase plateau for hard first-even ladders.
3. Add endpoint labels to `phase_half` rather than overriding too many arcs.
4. Use `pressure_k2` transitive reductions and SCCs as the pressure realization
   layer in THM-380.
5. Treat `safe_phase_gate` as an alarm coordinate: useful when it disagrees
   sharply with phase H.

## See Also

`07-reflections/lrc-tournament-gauge-lab-s506.md`,
`05-knowledge/results/lrc_tournament_gauge_lab_s506.out`,
HYP-1967, HYP-1968, HYP-1969, THM-373, THM-380.
