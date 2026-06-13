---
id: HYP-1970
status: OPEN
source: oracle-2026-06-01-S26; codex-2026-06-01-S26b
related:
  - HYP-1951
  - HYP-1967
  - THM-373
  - THM-374
  - THM-369
---

# HYP-1970: H is a half-turn loneliness meter, not a scalar max-gap meter

## Statement

For a half-turn circular tournament built from `n` points on the circle, the
Hamiltonian-path count `H` is a loneliness meter only at half-turn resolution:

```text
H = 1  iff  the point set lies in an open semicircle
       iff  max circular gap > 1/2
```

in the nondegenerate case.  With the repo's base-path tie completion, the
boundary `max_gap = 1/2` belongs to the same transitive-side reading, so the
computational probe records the condition as:

```text
H = 1  iff  max circular gap >= 1/2.
```

Above that transitive endpoint, `H` is not a function of the largest gap and is
not pointwise monotone in the largest gap.  It should instead be treated as a
finite circular-tournament cell invariant: a spread-class meter that records
cyclic arrangement information lost by max-gap, score, and 3-cycle summaries.

## Evidence

`04-computation/H_loneliness_meter_s26.py` samples half-turn circular tournaments
for `n=5..9` and computes exact `H`, max/min gaps, two-neighbor safe counts at
threshold `1/n`, score variance, 3-cycles, and gap entropy.

Sample results:

- `H=1 iff max_gap>=1/2` had zero mismatches for `n=5..9` under the
  tie-completed half-turn convention.
- Rounded max-gap buckets often map to multiple `H` values: `26/74` at `n=5`,
  `30/73` at `n=6`, `33/70` at `n=7`, `34/63` at `n=8`, and `36/57` at `n=9`.
- Pointwise monotonicity fails.  For `n=7`, the sample contains `H=105` with
  `max_gap=0.1812` and higher `H=151` with larger `max_gap=0.4881`.
- `H` separates score/3-cycle collisions: at `n=6`,
  `(scores=(2,2,2,3,3,3), c3=8)` splits into `H=41` and `H=45`; at `n=7`,
  `(scores=(2,2,3,3,3,4,4), c3=12)` splits into `H=123` and `H=137`.
- The two-neighbor safe count at `1/n` is not recoverable from `H`; by `n=8,9`,
  its sample correlation with `H` is close to zero.

The n=14 LRC overlay gives the practical warning.  At selected hard rows, the
half-turn clock remains high-spread:

```text
n14 initial       H=24104937, H/H0=1.000, max_gap=1/14
n14 row-parent    H=22168229, H/H0=0.920, max_gap=481/3696
n14 gate          H=17826951, H/H0=0.740, max_gap=481/3696
n14 double-gate   H=17826951, H/H0=0.740, max_gap=481/3696
```

These rows are not half-turn-bunched, but the LRC question is anchored at the
observer and threshold `1/14`.  This is why HYP-1967 requires two clocks.

## Consequences

The slogan "H is a loneliness meter" should always be expanded:

1. `H=1` is the exact open-semicircle detector.
2. `H>1` is a circular-tournament class feature correlated with spread.
3. Scalar loneliness for LRC requires extra local data: two-neighbor gaps,
   endpoint protection, and the anchored `1/n` clock.

## Next Tests

1. Prove the 3-cycle geometry formula:
   `c3(T)=binom(n,3)-sum_v binom(outdeg(v),2)` counts triples not contained in
   any semicircle.
2. Characterize the reach-word or arc-arrangement feature separating `H=41`
   from `H=45` and `H=123` from `H=137`.
3. Build the `alpha=1/k` comparator family and identify which exact reading
   replaces the `H=1 <-> 1/2-gap` theorem at LRC threshold `1/n`.

## Sources

- `04-computation/H_loneliness_meter_s26.py`
- `05-knowledge/results/H_loneliness_meter_s26.out`
- `07-reflections/H-as-a-loneliness-meter-s26.md`
- `05-knowledge/hypotheses/HYP-1951-runner-tournament-clock-circular-menu.md`
- `05-knowledge/hypotheses/HYP-1967-lrc-two-clock-corridor.md`
- THM-373
- THM-374
