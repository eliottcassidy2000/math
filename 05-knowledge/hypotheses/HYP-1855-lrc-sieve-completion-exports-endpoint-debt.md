---
id: HYP-1855
status: EXPLORATORY
source: codex-2026-05-31-S388
related:
  - THM-366
  - HYP-1844
  - HYP-1850
  - HYP-1851
---

# HYP-1855: Sieve completion exports endpoint debt

## Statement

The small-denominator sieve from THM-366 is only the first gate in the LRC
obstruction.  When an initial-segment-like speed set is modified to satisfy a
missing denominator by replacing a dropped speed `d` with `lcm(n,d)`, the
immediate rational witness disappears, but the obstruction reappears as
positive gap and endpoint debt.

In slogan form:

```text
small-denominator debt can be paid only by exporting it to endpoint layers.
```

The expected proof shape is a conservation law.  Each added lcm gate pays the
local modulus invoice, but because it is a sparse high-frequency interval
family, it creates or doubles boundary obligations on descendant denominators.

## Evidence

S388 compares naive `add n` swaps with `lcm(n,d)` repairs for `n=14,15,18`.
Every lcm-repaired transfer row in the table is still positive-gap.  Focused
endpoint audits of the largest repairs are terminal-core empty:

```text
n=14 drop 13 add 182:
  sieve-complete, gap/th=0.065934, unprotected=48, coreE=0

n=15 drop 14 add 210:
  sieve-complete, gap/th=0.061905, unprotected=24, coreE=0

n=18 drop 17 add 306:
  sieve-complete, gap/th=0.052288, unprotected=64, coreE=0
```

The seven-ladder is already sieve-complete at `n=14`, but still has

```text
gap/th=0.005411, unprotected=84, coreE=0.
```

This separates two regimes:

```text
missing denominator -> explicit a/m witness;
sieve complete      -> debt has moved to endpoint incidence.
```

## Predictions

1. In near-initial families, the map
   `drop d, add lcm(n,d)` never creates full measure at `n=14`; it opens or
   preserves a positive gap and leaves terminal endpoint core empty.
2. The endpoint debt created by lcm repairs is monotone in the number of
   primitive low-denominator obligations that the replacement speed absorbs.
3. A disproof search should treat sieve completion as a constraint, not as
   progress: after the sieve is satisfied, rank/private-pivot and
   critical-radius diagnostics dominate.
4. The same debt-export pattern should appear at the `n=18` battlefield,
   with `18=2*3^2` creating a richer denominator tree than `14=2*7`.

## Test Plan

1. Compute endpoint debt transfer matrices for all `drop d -> lcm(n,d)` rows
   at `n=14,15,18`.
2. Split debt by endpoint denominator and by first peel layer.
3. Compare the debt matrix with HYP-1844's quotient-layer export law.
4. Search for a counterexample to the pattern: sieve-complete, full-measure,
   and nonempty endpoint core.

## Sources

- THM-366
- `04-computation/lrc_small_denominator_sieve_s388.py`
- `05-knowledge/results/lrc_small_denominator_sieve_s388.out`
- HYP-1844
- HYP-1850
- HYP-1851
