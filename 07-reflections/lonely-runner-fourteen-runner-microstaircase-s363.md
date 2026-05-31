---
source: codex-2026-05-31-S363
status: research sprint
tags:
  - lonely-runner
  - fourteen-runners
  - micro-staircase
  - finite-checking
  - ansatz
---

# Lonely Runner Fourteen-Runner Micro-Staircase Sprint

## Frontier

The public finite frontier now reaches `k<=12` in the reduced Wills
formulation: thirteen total runners including the stationary runner.  The next
case is therefore:

```text
k = 13 moving speeds,  n = k+1 = 14,  threshold = 1/14.
```

The recent `k in {10,11,12}` proof has two ingredients:

1. use ansatz primes `p` to force many prime divisors of any counterexample;
2. for `k+1` prime, avoid the huge `c=k+1` lift of the tight class
   `(1,...,k)` with a polynomial/floor-vector argument.

For `k=13`, the obstacle is immediate: `k+1=14` is composite.  The field
argument does not directly apply.

## What Broke

I tested the direct analogue: can every lifted tight-class vector in
`(Z/14Z)^13` with a zero coordinate be resolved by times of the form

```text
t = s/14 + r/14?
```

No.  The vector

```text
v = (7,4,9,6,7,8,5,0,1,12,13,12,7)
```

blocks every pair `(s,r) mod 14`.  In the proof language, for every candidate
there is some coordinate landing in residue `0` or `13`, so the candidate
interval still touches the forbidden edge.

This is good news in disguise: it says the next proof should not waste time
trying to copy the prime polynomial lemma word for word.

## What Worked

Actual prime ansatz times are finer:

```text
t = s/14 + r/p.
```

The floor vector

```text
R_alpha(i) = floor(14 * {i alpha})
```

has many more cells than the coarse `alpha=r/14` samples.  The obstruction
above is resolved by a tiny staircase cell:

```text
R_alpha = (0,0,0,1,1,1,2,2,2,3,3,3,4)
alpha in [2/91, 1/42), width 1/546.
```

With `s=1`, the obstructing vector becomes

```text
(7,4,9,7,8,9,7,2,3,1,2,1,11),
```

which avoids `0` and `13` in every coordinate.

So the composite case seems to want a micro-staircase lemma:

```text
coarse prime-field cells fail,
fine floor-vector cells may still resolve all tight lifts.
```

## Verifier Data

I cloned the public `vzsky/13-lonely-runners` code to a temp directory and ran
its experimental `LrcVerifier<13>` path on primes through `101`.

Every tested prime finished with empty final set.  More interestingly, after
`I(13,p,1)` was computed, the first `c=2` squeeze already killed everything.
The large cost was the initial cover search:

```text
p=89:  |I(13,p,1)| = 11804159, find_cover = 91.812971s
p=101: |I(13,p,1)| = 12697411, find_cover = 113.117640s
```

This matches the paper's closing warning: progress toward `k=13` needs a
better understanding of tuples without ansatz witnesses at level `l=1`.

## Candidate Route

The next proof attempt should split into two independent lemmas.

### Lemma A: Tight-Class Micro-Staircase

For every non-gcd lifted tight-class vector `v in (Z/14Z)^13`, find a cell of
the `R_alpha` arrangement and an `s mod 14` such that

```text
s v_i + R_alpha(i) notin {0,13} mod 14
```

for all `i`.

If the resolving cell has width at least `1/M`, then all primes `p>M` have a
usable grid point.  Smaller primes can be checked directly.

### Lemma B: Initial Bad-Set Compression

For `I(13,p,1)`, the published verifier's data suggests an unexpectedly strong
claim:

```text
all bad sets that survive l=1 have a c=2 witness obstruction that is
destroyed immediately after lifting.
```

A proof-grade version would be a pruning theorem inside the DFS:
as soon as a partial cover has the right micro-staircase witness, it can be
discarded before full depth `13`.

## Why This Feels Alive

The obstruction vector is not a failure of the ansatz; it is a failure of using
too coarse a shadow of the ansatz.  The prime grid already contains the missing
cells.  The proof task is to name and control those cells before the search
explodes.

That is a much sharper problem than "prove fourteen lonely runners":

```text
classify floor-vector cells for n=14
rank tight-lift residues by resolving cell width
turn the resulting certificate into cover-search pruning
```

## Artifacts

- `04-computation/lonely_runner_k13_microstaircase_s363.py`
- `05-knowledge/results/lonely_runner_k13_microstaircase_s363.out`
- `05-knowledge/hypotheses/HYP-1817-lrc-k13-microstaircase.md`
