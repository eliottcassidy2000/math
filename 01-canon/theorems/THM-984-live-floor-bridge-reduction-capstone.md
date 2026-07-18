---
id: THM-984
title: Live-floor sampling bridge and trapped-core reduction target
status: PAPER CLAIM / COMPUTATIONAL REFEREE. The direct `liveCount > 0` consumer is now Lean; the general endpoint-counting inequality, exact modulus arithmetic, T_s composition, and exhaustive trapped-core reduction are not yet formalized.
source: kind-pasteur-2026-07-17-S128 (cont.47; owner: finish the math, then all the formalization)
depends_on:
  - THM-979 (the proposed B₅-route modulus)
  - LRCQ25Obstruction and LRCLiveCountLonely (the exact live-count consumer)
  - THM-928/932/933 (cascade/gluing), the assembly branches (LRC14GrandAssembly)
related:
  - the fee ledger (death-star THM-938..945) and witness machinery (THM-936) — the two live routes to the remaining μ₀-floors
  - OPEN-Q inf-L (the classical nucleus, now in its final cut)
script: 04-computation/live_floor_bridge_referee_kps_S128c47.py -> 05-knowledge/results/live_floor_bridge_referee_kps_S128c47.out
---

# THM-984 — live-floor bridge and reduction target

## Audited bridge statement

For a signed integer family use

```text
E(V) = 2 * sum_i |v_i|.
```

If the good set in one period is represented by finitely many interval
components with at most `E(V)` boundary points, the intended discrepancy bound
is

```text
liveCount(v,q) >= q * μ₀(v) - E(V).
```

Consequently, when `μ₀(v)>0`, the exact integer choice
`q₀=ceil(2E(V)/μ₀(v))` would give `liveCount(v,q₀)>0`.  The repository's
referee checks two positive-measure examples, but computes `q₀` through a
floating-point conversion; it is evidence, not the general proof or exact
ceiling formalization.

`LRCLiveCountLonely` now closes the consumer side in Lean:
`0 < liveCount v q` directly gives `∃ t, Lonely 14 v t`.  No
`CoverageCapped`, depth-six race, or call to `lonely_of_census` is needed once
a live multiplier exists.

## Correction to the original control claim

The earlier text said that the tight family `(1,...,13)` has zero live points
at every modulus because its good measure is zero.  This is false: at `q=14`
the exact histogram has six live multipliers and `B5=5`.  The referee only
observes zero live points at `q=2003` and `q=8009`.  Measure zero does not imply
that every rational grid misses the safe set; it merely prevents this
positive-measure bridge from selecting a modulus.

## Relation to the B5 route

| branch | route | modulus | conditional on |
|---|---|---|---|
| dissociated | B₅ (THM-979 proposal) | ≈ 51,900·Σ|v| | the unformalized T_s assembly |
| μ₀ > 0 | live (this proposal) | ceil(2E/μ₀) | endpoint discrepancy + a μ₀-floor |

These are complementary large-modulus routes.  They do not activate the sharp
`q<=98` seven-stalk theorem: thirteen distinct positive magnitudes already
give `E>=182`, hence the live-route modulus is at least `364` when `μ₀<=1`.

## Honest remaining target

Positive safe measure on every genuinely trapped primitive core remains a
natural sufficient target.  The current Lean corpus proves the grand-assembly
reduction to `SafeMeasureFloorPrimitive`; it does not yet prove that all of the
additional prose cuts in the proposed trapped-core list form an exhaustive
formal partition.  Nor does it prove the required positive measure.

## Named next

- Formalize the finite-component endpoint discrepancy with `2*sum|v_i|` and
  exact ceiling arithmetic.
- Formalize the T_s/B5 sampling route rather than treating its referee as a
  kernel theorem.
- Prove the exhaustive trapped-core partition and the positive-measure floor.
- Bridge either large-modulus route to the remaining colored/depth-six-seven
  payment, or supply a genuinely adaptive small modulus.
