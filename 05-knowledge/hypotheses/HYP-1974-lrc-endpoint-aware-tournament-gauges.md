---
id: HYP-1974
status: OPEN
source: codex-2026-06-01-S507
related:
  - HYP-1931
  - HYP-1940
  - HYP-1951
  - HYP-1967
  - HYP-1968
  - HYP-1971
  - HYP-1972
---

# HYP-1974: LRC loneliness metrics need endpoint-aware tournament gauges

## Statement

Hamiltonian-path count `H` is a useful loneliness metric only after the arc
gauge has the right threshold in it.  The pure half-turn circular tournament is
a sharp `1/2`-gap meter, but it is too coarse for the LRC threshold
`1/(k+1)`.  For LRC Tournament Analysis, use a two-layer gauge family:

1. **Status gauges** for witness detection.  Unsafe/safe dominance orients an
   edge by whether one runner is inside the forbidden origin interval
   `||s_i t|| < 1/(k+1)`, with half-turn ties inside equal status blocks.
   Its score sequence makes "zero unsafe runners" a visible tournament event.
2. **Shape gauges** for non-tautological geometry.  Theta-close and
   antipodal-open switches flip a fixed Hamiltonian base path using pair
   separation near `2/(k+1)` or near antipodal separation.  These do not directly
   read the target variable, but they preserve pair geometry that can separate
   safe corridors from generic unsafe phases.
3. **Boundary/pressure gauges** for proof certificates.  Endpoint-wall pressure
   and LRC slack bands are weaker witness classifiers, but they mark boundary
   pressure near `||s_i t|| = 1/(k+1)` and should feed endpoint-core and
   pressure-DAG scripts.

Thus the useful metric is not "H alone"; it is `(gauge, H, score sequence,
cycle/SCC fingerprints)` with the gauge declared.  For LRC, `H+score sequence`
under endpoint-aware gauges is the basic loneliness fingerprint.

## Evidence

`04-computation/lrc_loneliness_tournament_gauges_s507.py` tests ten arc gauges
on six small speed families over the grid `t=a/840`.

The declared Tournament Analysis object is:

```text
vertices             = moving runners in the speed set
pairwise observable  = phase separation, origin distance, endpoint slack,
                       endpoint-wall nearness, origin side, or deletion relief
switch/gauge         = named arc rule
tie Hamiltonian path = increasing speed/index order
target               = (k+1) * min_i ||s_i t||
```

Stored output:

```text
05-knowledge/results/lrc_loneliness_tournament_gauges_s507.out
```

Key fingerprints:

```text
unsafe_dominance:
  H+score pure safe/unsafe fraction 0.996
  safe isolation 0.845, unsafe isolation 0.866
  perfect safe isolation on consecutive-4/5/6/7 and prime-6

safe_dominance:
  same aggregate performance as unsafe_dominance, with reversed status order

theta_close_switch:
  best numeric monotone feature = score_width, Spearman rho -0.266
  mixed-7 safe isolation 0.767 with H+score span 0.552

antipodal_open_switch:
  best numeric monotone feature = scc_count, Spearman rho -0.207
  mixed-7 safe isolation 0.433 with H+score span 0.480

half_turn:
  safe isolation only 0.028 across all cases
```

The half-turn result matches HYP-1968: it is a spread and `1/2`-threshold meter,
not an LRC witness detector.  The status-dominance gauges are nearly tautological
but valuable as a sanity-check witness meter.  The theta-close and antipodal-open
gauges are the promising shape-only candidates because they extract LRC-scale
structure without directly asking whether every runner is safe.

This complements HYP-1972: the metric-vector session audits a broader set of
marked-origin and pressure criteria, while this benchmark isolates which
completed tournament fingerprints can classify or separate safe times on small
clock grids.

## Predictions

1. Future LRC scripts should report an endpoint-aware tournament fingerprint
   beside any scalar gap: at minimum status-dominance score sequence,
   theta-close `H`, theta-close score width, antipodal-open SCC count, and
   endpoint-wall pressure layers.
2. If a hard row has a genuine lonely corridor, theta-close and antipodal-open
   movies should show stable score/H signatures across the corridor, not just
   isolated time samples.
3. Endpoint-wall pressure should align with boundary witnesses and endpoint
   peel layers, even though it is poor at classifying safe times by itself.
4. A counterexample-shaped row should make status dominance report persistent
   nonzero unsafe count while shape gauges still look deceptively corridor-like;
   this mismatch is a useful warning signal.
5. The right LRC "loneliness meter" is likely a small vector of tournament
   fingerprints, not a single scalar `H`.

## Next Tests

- Add these gauges to the n=14/n=18 corridor scripts and compare witness
  intervals against theta-close and antipodal-open movie signatures.
- Run the same gauge menu on product-sieve surviving tuples from the current
  fixed-case proof frontier.
- Test whether endpoint-wall pressure sources overlap HYP-1961/HYP-1968
  pressure-DAG source/sink and private-endpoint layers.
- Replace single-time classification with corridor persistence: Hamming
  diameter, stable score sequence runs, and edge-flip counts across each lonely
  interval.

## See Also

- `04-computation/lrc_loneliness_tournament_gauges_s507.py`
- `05-knowledge/results/lrc_loneliness_tournament_gauges_s507.out`
- `07-reflections/lrc-loneliness-tournament-gauges-s507.md`
- HYP-1931, HYP-1940, HYP-1951, HYP-1967, HYP-1968, HYP-1971, HYP-1972
