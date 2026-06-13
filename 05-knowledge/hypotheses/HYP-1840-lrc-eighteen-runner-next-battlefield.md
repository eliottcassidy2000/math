---
id: HYP-1840
status: EXPLORATORY
source: codex-2026-05-31-S383
related:
  - HYP-1831
  - HYP-1837
  - HYP-1838
  - HYP-1839
---

# HYP-1840: Eighteen runners are the next best LRC battlefield above fourteen

## Statement

Among the first denominators above `n=14`, the best next proof/disproof target
is `n=18`: seventeen moving speeds against the stationary runner, with
threshold `1/18`.

The claim is not that `n=18` has the largest raw pressure.  Larger composite
denominators like `n=20,22,24` already look louder by raw endpoint/gap metrics.
The claim is that `n=18` is the best balance of:

```text
small enough for exact certificates,
large enough to show new recursive behavior,
arithmetically rich enough to expose the proof mechanism.
```

## Evidence

`lrc_above14_runner_favorites_s383.py` audits largest-proper-divisor ladders
for `15 <= n <= 24`.

The `n=18` row:

```text
n=18 = 2*3^2
lpd=9
best skip=8
gap/th=0.005682
unprotected endpoints=176
pressure/n^2=95.604938
phi(n)/(n-1)=0.352941
product-sum comparison m(18)=24 from seed (2,3,4) plus 15 ones
```

This makes `n=18` the first denominator above fourteen with all of the
following visible at once:

1. a tiny largest-proper-divisor quotient-ladder gap;
2. a very low unit-skeleton density, second only to `n=24` in the scanned
   window;
3. an even gate plus a square `3`-torsion layer;
4. a genuinely multi-factor product-sum comparison `(2,3,4)`;
5. exact endpoint pressure that is still small enough to audit directly.

The triage also separates neighbors:

```text
n=16: clean 2-adic proof laboratory;
n=18: favorite next battlefield;
n=21: pretty 3-by-7 transfer target;
n=24: divisor-rich stress test, likely after the mechanism is known.
```

## Predictions

1. The `n=18` proof branch should split first by the `18`-gate unit skeleton,
   then by the lpd layer `9`, then by `3`-torsion descendant endpoints.
2. A no-`18`-multiple branch should be dismissed by the same unit-skeleton
   obstruction as HYP-1839.
3. The lpd ladder at `d=9` should have a small endpoint-debt certificate
   analogous to, but cleaner than, the `n=14` seven-ladder debt.
4. If a counterexample search is ever viable above fourteen, `n=18` should
   show the first credible endpoint-protection cycle before larger denominators
   like `20,22,24`.

## Sources

- `04-computation/lrc_above14_runner_favorites_s383.py`
- `05-knowledge/results/lrc_above14_runner_favorites_s383.out`
- HYP-1831
- HYP-1837
- HYP-1838
- HYP-1839
