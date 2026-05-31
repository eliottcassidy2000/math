---
id: HYP-1857
status: OPEN
source: codex-2026-05-31-S389
related:
  - HYP-1838
  - HYP-1841
  - HYP-1842
  - HYP-1844
  - HYP-1850
  - HYP-1851
  - HYP-1853
---

# HYP-1857: n=16 LRC gates obey a Cayley-Dickson dyadic debt law

## Statement

At Lonely Runner denominator `n=16`, a counterexample must contain a speed
divisible by `16`, but every such dyadic gate pays for closing the old unit
endpoint skeleton by exporting endpoint debt to deeper dyadic layers.

The hoped-for proof theorem is:

```text
Every primitive n=16 all-protected endpoint system with a 16-gate
either has a positive gap, or its labelled endpoint-protection cycle
has a private dyadic leaf.
```

Equivalently, the dangerous branch should not be "has a 16-multiple" alone.
It should be "has a 16-multiple and also defeats the dyadic debt norm plus
private-leaf descent."  The S389 data suggest that second condition is the
impossible one.

## Cayley-Dickson analogy

The Cayley-Dickson tower doubles dimensions:

```text
1, 2, 4, 8, 16, ...
```

At dimension `16`, the sedenion layer exposes zero-divisor behavior: a formal
doubling operation still exists, but the old normed-division intuition has
broken.  The LRC `n=16` gate behaves similarly.

It protects the old unit endpoints `a/16`, but not by making the system
globally safe for a counterexample.  It creates a new dyadic kernel at
denominators `32,64,128,256,...`.  The proof attempt should therefore track a
norm-like ledger:

```text
visible gap shrinks, but endpoint debt grows.
```

## Evidence

The unit-gate branch is an actual lemma.  For odd `a`, the unit endpoint
`t=a/16` is strictly protected by speed `v` only if

```text
v*a == 0 mod 16.
```

Since `a` is odd, this is equivalent to `16 | v`.  Therefore no speed divisible
by `16` implies an odd unit endpoint remains a boundary witness, so no open
cover counterexample can live in the no-gate branch.

`lrc_n16_cayley_dickson_attempt_s389.py` then audits the gated branch.  The
exact dyadic ladders are all positive-gap and endpoint-core-empty:

```text
d=2  ladder: gap/th=1/33,  unprotected=34,  debt_norm=34/33
d=4  ladder: gap/th=1/66,  unprotected=68,  debt_norm=34/33
d=8  ladder: gap/th=1/132, unprotected=140, debt_norm=35/33
d=16 ladder: gap/th=1/264, unprotected=280, debt_norm=35/33
```

The visible gap halves under doubling, but unprotected endpoints double
up to the small `34/35` transition.  The first unprotected endpoint layer also
moves one dyadic level deeper:

```text
d=2  -> v2(endpoint denominator)=5
d=4  -> v2(endpoint denominator)=6
d=8  -> v2(endpoint denominator)=7
d=16 -> v2(endpoint denominator)=8
```

The normalized scalar quotient gives a second obstruction.  At `n=16`, the
exact finite cell system has `1152` alpha patterns and `18432` candidates.
Every scalar ramp is full, but puncturing a scalar ramp is expensive:

```text
radius-1 best defect: 128 missed cells, half-turn at coordinate 15
radius-2 best defect: 160 missed cells, half-turns at coordinates 10 and 15
```

Fixed high-speed gate samples with two or three `16`-multiples also remained
positive-gap and endpoint-core-empty, with hundreds of unprotected endpoints.

## Predictions

1. The full `n=16` proof should split cleanly into no-`16`-gate and
   has-`16`-gate branches.
2. In the gated branch, any sequence of repairs that shrinks the visible gap
   must increase a dyadic endpoint-debt functional.
3. A true all-protected `n=16` endpoint core, if abstractly drawn, should fail
   integer realizability or have a private dyadic leaf in the labelled incidence
   matrix.
4. The `n=16` scalar-puncture moat should be provable symbolically from the
   half-turn action on the last coordinates, analogously to the `n=14` moat
   but with pure 2-adic folds rather than `2*7` CRT residue.
5. `n=16` should be the clean laboratory for proving the debt-export law
   (HYP-1844) before returning to the messier `n=14` and `n=18` cases.

## Sources

- `04-computation/lrc_n16_cayley_dickson_attempt_s389.py`
- `05-knowledge/results/lrc_n16_cayley_dickson_attempt_s389.out`
- `07-reflections/lonely-runner-n16-cayley-dickson-attempt-s389.md`
