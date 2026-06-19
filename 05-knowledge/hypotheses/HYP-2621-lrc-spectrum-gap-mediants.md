---
id: HYP-2621
title: LRC spectrum gap mediants and AP-defect constant ladder
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S16
depends_on:
  - THM-524
related:
  - HYP-2052
  - OPEN-Q-108
---

# HYP-2621 - LRC Spectrum Gap Mediants and AP-Defect Constant Ladder

## Claim

Let

```text
M(S) = max_t min_{v in S} ||v t||
```

for a primitive `k`-speed set, with AP floor `1/(k+1)`.  The doubled-top
family

```text
D_k = {1,2,...,k-1,2k}
```

realizes the explicit second-mediant value

```text
M(D_k)=2/(2k+1),
M(D_k)-1/(k+1)=1/((k+1)(2k+1)).
```

Thus the G2/lift-depth denominator forced by this family is exactly
`(k+1)(2k+1)`.

The new observation is that the doubled-top value is not the only low
AP-defect constant.  In the one-defect family

```text
A_{k,r} = {1,2,...,k} \ {k-1} union {r(k-1)},
```

the `r=3` branch has a residue split in the tested range:

```text
k == 1 mod 30        -> M(A_{k,3}) = 1/(k+1)       (AP-tight again)
k == 7,13,19,25 mod30 -> M(A_{k,3}) = 3/(3k+2)    (third-mediant)
```

On the third-mediant branch,

```text
M(A_{k,3}) - 1/(k+1) = 1/((k+1)(3k+2)).
```

This improves the upper bound for the spectral gap constant on an infinite
tested residue-class family, but it still has the same `Theta(1/k^2)` scale.
No computation here finds an `o(1/k^2)` dip.

## Computation

Script:

- `04-computation/lrc_spectrum_gap_mediants_codex_s16.py`
- output: `05-knowledge/results/lrc_spectrum_gap_mediants_codex_s16.out`

The script uses the THM-524 exact candidate envelope:

```text
t = (2a+1)/(2v_i),     t = m/(v_i+v_j),     t = m/|v_i-v_j|.
```

It evaluates `M(S)` with exact rational arithmetic.

## Doubled-Top Family

Verified exactly for `k=2..30`:

```text
k=2   M=2/5    depth=15
k=3   M=2/7    depth=28
k=4   M=2/9    depth=45
...
k=30  M=2/61   depth=1891
```

The depth column is

```text
1 / (2/(2k+1) - 1/(k+1)) = (k+1)(2k+1).
```

This confirms the user's explicit family and pins the G2 lift-depth obstruction
for that route.

## Third-Mediant Branch

For

```text
A_{k,3} = {1,...,k}\{k-1} union {3(k-1)}
```

the stored table through `k=103` gives:

```text
k=7,13,19,25 mod30:  M=3/(3k+2)
k=1 mod30:           M=1/(k+1)
```

Concrete rows:

```text
k=7   S=(1,2,3,4,5,7,18)                  M=3/23   depth=184
k=13  S=(1,...,11,13,36)                  M=3/41   depth=574
k=31  S=(1,...,29,31,90)                  M=1/32   AP-tight
k=37  S=(1,...,35,37,108)                 M=3/113  depth=4294
k=91  S=(1,...,89,91,270)                 M=1/92   AP-tight
```

The `k=31,61,91` tight returns are the guardrail: the branch is not simply
`k == 1 mod 6`; the mod-`30` address matters.

## Multiplier Ladder Probe

The same one-defect shape was probed for `r=2..5`, `k<=60`:

```text
r=2: 0 below-family hits; many AP-tight returns.
r=3: the third-mediant branch above.
r=4: one below-family hit at k=31, M=4/127, t=72/127.
r=5: 0 below-family hits.
```

This suggests a real AP-defect constant ladder, but not a monotone or dense one.
The immediate next computation should scan `A_{k,r}` by residue classes in
`(k,r)` rather than sweeping arbitrary boxes.

## Bounded Small-k Scans

Finite primitive boxes still match the mediant-tight statement for the smallest
sizes:

```text
k=2, B=40: best above floor = 2/5
k=3, B=36: best above floor = 2/7
k=4, B=24: best above floor = 2/9
k=5, B=20: best above floor = 2/11
```

These are finite scouts, not proofs of `sigma_2(k)`.

## Tournament Analysis

The script explicitly challenges the default vertex set.  It considered:

```text
runners,
residue classes mod 2k+1,
pair-crossing denominators,
AP-defect families,
finite scan boxes,
proof obligations.
```

It uses proof-obligation/family routes as vertices, because that quotient
preserves the predicate `M(S)-1/(k+1)` and its reciprocal lift-depth while
discarding individual witness-time geometry except for the argmax address.

The proof-route tournament is transitive:

```text
infinite_lower_bound_for_g(k)
> third_mediant_family
> doubled_top_family
> one_defect_scan
> two_defect_scan
> bounded_exhaustive_box
> raw_runner_vertices
```

Score histogram: `[6,5,4,3,2,1,0]`; directed `3`-cycles: `0`.

## Status

Partially confirmed by exact rational computation.  The doubled-top identity is
fully verified in the stored range and has a clear elementary witness.  The
third-mediant and multiplier-ladder patterns are exact in the stored range but
still need a symbolic residue-class proof.

The lower-bound problem remains open:

```text
Does g(k) ever dip below c/k^2 along an unbounded sequence for every c>0?
```

The present evidence says "not yet"; it improves constants but not order.
