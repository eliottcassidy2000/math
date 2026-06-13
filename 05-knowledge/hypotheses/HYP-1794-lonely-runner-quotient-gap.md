---
id: HYP-1794
status: EXPLORATORY
source: codex-2026-05-30 S356
related:
  - HYP-1781
  - HYP-1783
  - HYP-1785
  - HYP-1802
  - THM-349
  - THM-355
  - THM-357
---

# HYP-1794: Lonely Runner as a Quotient-Gap Certificate

## Statement

The reduced lonely runner conjecture can be attacked as a finite
quotient-gap theorem.

For a primitive integer speed set

```text
V = {v_1, ..., v_k},
```

pull the open forbidden arcs

```text
dist(v_i t, Z) < 1/(k+1)
```

back to the time circle `R/Z`.  Every forbidden boundary lies in the finite
quotient

```text
Q(V) = (k+1) * lcm(V).
```

The conjectural witness is therefore either:

1. a positive-length gap between occupied forbidden fibers; or
2. a boundary residue where all runners are exactly outside the open forbidden
   arcs.

The hypothesis is that every minimal counterexample candidate has a detectable
transport residue in this boundary quotient before one needs continuous
Diophantine approximation.

## S357 Update

The positive-gap/boundary-residue split is a finite-open-cover trichotomy, not
the hard part by itself.  Since the forbidden set is a finite union of open
intervals, any nonempty complement is either a positive interval or a finite
set of boundary points.

Thus the sharpened target is HYP-1802: a genuine counterexample would be a
full open cover of the circle, equivalently full forbidden measure plus every
forbidden endpoint strictly protected by another forbidden interval.

THM-357 now proves this trichotomy and endpoint-protection certificate
equivalence.  HYP-1794's remaining value is the quotient-gap language:
positive cells and unprotected endpoints are the two surviving finite residue
types.

## Why This Belongs To Residue Calculus

The map

```text
t in R/Z  |->  (v_1 t, ..., v_k t) in (R/Z)^k
```

projects one time parameter into a torus.  The forbidden set is a union of
coordinate slabs, and the desired lonely time is the residue left outside their
pullback.  This is the same shape as:

- good-cut interval gas: interval unions on a path leave bucket gaps;
- THM-355 quotient gaps: empty fibers force zero transport rows and columns;
- deletion-residue rank: deleting a high-loss vertex leaves a small surviving
  obstruction;
- Erdős-Straus base-42 notes: residue classes are covered by finite identity
  families.

The useful object is not the scalar speed set.  It is the finite boundary
quotient plus the residue that survives the forbidden interval cover.

## Evidence

`04-computation/lonely_runner_residue_probe_s356.py` computes the pulled-back
forbidden interval union exactly over `Fraction`.

The S356 sample table contains 16 deterministic speed sets:

- Initial segments `{1,...,k}` for `k=2..9` have no positive gap, but have
  boundary witnesses such as `t=1/(k+1)`.  These are tight quotient-boundary
  cases.
- Lacunary, arithmetic, prime, CRT-mixed, and random speed sets all have
  positive gaps in the sample.
- The mixed CRT sample with eight speeds has the tightest positive gap in the
  table:

```text
V=(5,8,13,21,34,55,89,144)
threshold=1/9
max_gap=22/7209
max_gap/threshold=0.027466
boundary_modulus=9814044240
```

This supports a two-regime picture:

```text
tight boundary residue  vs.  positive quotient gap
```

## Predictions

1. Known tight families should have no positive gap but many boundary witnesses
   in `Q(V)`.
2. Hard random or CRT-correlated families should first appear as very small
   positive gaps, not as absence of boundary witnesses.
3. A minimal counterexample search should be more efficient if stratified by
   quotient features:

```text
boundary_modulus,
component_count,
positive_gap_count,
boundary_witness_count,
max_gap/threshold,
boundary residue pair.
```

4. The row/column gap vocabulary of THM-355 should have an LRC analogue: once
   the forbidden interval quotient is discretized at `Q(V)`, every surviving
   boundary residue is a silent target column for the forbidden transport.

## Test Plan

1. Extend the probe from handpicked samples to bounded primitive speed sets.
2. Rank sets by `max_gap/threshold` and isolate the smallest positive gaps.
3. Add a boundary-only mode that records tight witnesses when positive gaps
   vanish.
4. Compare the hardest speed sets against known LRC proof reductions and
   minimal-counterexample bounds.
5. Build the endpoint-protection graph of HYP-1802 and try to prove that a
   full-measure forbidden union must have an unprotected endpoint.

## Sources

- `04-computation/lonely_runner_residue_probe_s356.py`
- `05-knowledge/results/lonely_runner_residue_probe_s356.out`
- `07-reflections/famous-problem-residue-bridges.md`
- HYP-1781, HYP-1783, HYP-1785.
