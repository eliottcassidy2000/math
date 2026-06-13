---
id: HYP-2082
status: SUPPORTED by exact bounded audit - global theorem open
source: codex-2026-06-03-S570
related:
  - HYP-2081
  - HYP-2080
  - HYP-2075
  - HYP-2059
  - THM-369
---

# HYP-2082: the right exact LRC checker is witness-or-core - pair pinches recover the maximin, the n-clock isolates tightness, and endpoint cores are the dual obstruction ledger

**Grounded by exact bounded audit (`lrc_witness_or_core_s570.py`):**
- For every audited primitive box
  - `k=3, max_speed=20` (`997` rows),
  - `k=4, max_speed=16` (`1745` rows),
  - `k=5, max_speed=13` (`1281` rows),
  - `k=6, max_speed=11` (`462` rows),
  the exact audited maximin `M(S)` is recovered by the cheap primal catalogue: pair-sum pinch or antipode.
- Route split is simple: the resonance-maximal tight rows route through the `n`-clock, while the generic rows route through pair-sum pinch. In the audited boxes the route histogram is exactly `{n_clock, pair_sum}` with no surviving fallback/core cases.
- Sample `n=14` families fit the same picture: AP and `V*` certify at the `n`-clock threshold `1/14`; Fibonacci, Sidon, and random rows certify at pair-sum times strictly above threshold.
- Every audited primitive box peels to an empty endpoint-protection core (`nonempty_core=0`), so the core behaves as a dual obstruction ledger rather than a routine witness source.
- Bonus bounded margin probe: against the speculative second floor `2/(2n-1)`, the only audited rows below that value were already `n`-clock-tight rows with `M(S)=1/n` (counts `1,2,2,1` in the four scanned boxes). No off-stratum pair-sum row in those boxes fell below the second-gap floor.

**Hypothesis / proof shape:**
- For primitive integer LRC rows, the exact maximin should always be attained by a pair-pinch time `t=m/(v_a+v_b)` or by an antipode.
- Tightness `M(S)=1/n` should be a special arithmetic stratum detected by the `n`-clock / CRT gear picture.
- Strict failure `M(S)<1/n` should force a nonempty labelled endpoint-protection core, giving a finite dual certificate instead of an analytic "almost witness".

**Interpretation:** the continuous-time question is probably best treated as a primal-dual finite checker:

```text
find a cheap primal witness,
or export a finite protected core.
```

This matches the repo's current synthesis:
pair-sum pinch collapses time, the `n`-clock handles tight equality, and the
endpoint core is the residual obstruction object.

**Honest scope:** the bounded audits strongly support the architecture but do
not prove the global theorem. The open step is to show that the empty-core
behavior persists globally, or that any surviving core contradicts the pinch /
gear / margin structure.

**See:** `07-reflections/lrc-1d-periodic-maximin-breakthrough-route-s570.md`, `04-computation/lrc_witness_or_core_s570.py` (+.out); HYP-2081 (which clocks matter), HYP-2080 (orbit/resonance), HYP-2075 (pair-sum completeness), HYP-2059 (pair-pinch exactness), THM-369 (divisibility witness sieve).
