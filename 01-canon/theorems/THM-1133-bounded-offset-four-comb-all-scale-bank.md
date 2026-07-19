---
id: THM-1133
title: All-scale sharp four-comb theorem for every offset shape of span at most thirty
status: PROVED / COMPUTER-ASSISTED — THM-1123 bottom bank, 3,539,936-row exact finite complement, and THM-1129 rectangle tails cover every legal K for all 495 cores and all 4,060 normalized offset shapes of span at most 30
source: codex-2026-07-18-S73 bounded-offset continuation
depends_on: [THM-1123, THM-1129]
related: [THM-1126, THM-1127, THM-1128, MISTAKE-164]
verification:
  - 04-computation/lrc14_r5_bounded_offset_all_scale_complement_exact_codex_S73.cpp
  - 04-computation/lrc14_r5_bounded_offset_complement_fraction_replay_codex_S73.py
  - 05-knowledge/results/lrc14_r5_bounded_offset_all_scale_complement_exact_codex_S73.out
  - 05-knowledge/results/lrc14_r5_bounded_offset_complement_fraction_replay_codex_S73.out
---

# THM-1133 — all-scale bounded-offset four-comb theorem

Let `P` be any eight-element subset of `{1,...,12}` and let

```text
A={0<a<b<c<=30},             K>13 max(P).
```

Put

```text
E(P,A,K)={t in R/Z : ||pt||>=1/14 for p in P,
                       ||(K+alpha)t||>=1/14 for alpha in A}.
```

Then `E(P,A,K)` contains a component of length `L` satisfying

```text
7(K+c)L>1.                                                (1)
```

This quantifies over all

```text
C(12,8)*C(30,3)=495*4060=2,009,700                       (2)
```

core/offset rays and every legal base `K`.  It is the complete bounded-offset
stratum of the open `r=5` sharp component problem.

## 1. Three exact ranges

Write

```text
lo(P)=13 max(P)+1.
```

The proof partitions each ray without gaps.

### Bottom

THM-1123 proves (1) whenever

```text
lo(P)<=K,              K+c<=lo(P)+39.                  (3)
```

These are exactly the first `40-c` rows on the ray.

### Tail

THM-1129 assigns each pair `(P,A)` an exact one-sided polygon-rectangle
certificate.  If its certified time width is `d(P,A)`, then

```text
K>=ceil((8/7)/d(P,A))                                  (4)
```

places a complete `1/7` vertical-arc preimage inside the rectangle.  The
resulting safe component has length

```text
1/(7K)>1/(7(K+c)),                                     (5)
```

which proves (1).  The worst individual threshold in (4) is `832`, but
`1,802,872` of the `2,009,700` rays already satisfy (4) at `K=lo(P)`.

### Exact finite complement

The only remaining rows satisfy

```text
lo(P)+40-c <= K < ceil((8/7)/d(P,A)).                  (6)
```

The THM-1129 Fraction atlas counts them exactly:

```text
raw legal rows below the individual tails       6,040,056
rows already in the bottom range (3)            2,500,120
rows in the complement (6)                      3,539,936. (7)
```

The C++ referee independently reconstructs every `d(P,A)` from the exact
core endpoints and offset-center collision cells.  It reproduces every count
in (7), performs exact endpoint subtraction on all `3,539,936` complement
rows, and finds

```text
failures of 7(K+c)L>1: 0.                              (8)
```

Equations (3), (6), and (4) are consecutive ranges, so (8) completes the
all-scale theorem.

## 2. Exact residual extremal

The smallest sharp metric inside the newly checked complement is attained at

```text
P={1,2,4,5,7,8,9,11},
A=(0,4,6,8),
K=186,
(k1,k2,k3,k4)=(186,190,192,194),
L=86/61845,
7k4L=16684/8835=1.888... .                              (9)
```

Its two longest intervals are the reflected pair

```text
[267/2660,265/2604], [2339/2604,2393/2660].            (10)
```

For this pair the direct rectangle replay gives

```text
d(P,A)=1/308,        individual tail start=352,         (11)
```

with a witness cell beginning at `71/154`.  Its legal start is `144`, and
its first row outside THM-1123 is `176`, so `K=186` genuinely belongs to the
finite complement rather than either neighboring certificate.

The analytic tail metrics from (5) approach `1` from above as `K` grows.
Thus (9) is the residual finite-bank extremal, not a positive uniform margin
for the all-scale theorem.

## 3. Exact arithmetic and independent guards

The primary success criterion uses only integer-rational endpoints and
signed `__int128` cross-products.  It emits a deterministic row digest

```text
15941569285517864349.                                   (12)
```

As a shared-boundary guard it independently recovers THM-1123's finite
extremal

```text
P={1,2,4,5,7,9,11,12}, (158,160,162,164),
L=41/25920, 7k4L=11767/6480.                            (13)
```

Optimized, unoptimized, and ASan/UBSan C++ executions are byte-identical to
the frozen primary output; the sanitizer run reports no error and Clang static
analysis is clean.

The Python replay uses a structurally different engine: it builds the
complete simultaneous breakpoint arrangement for all twelve speeds and
classifies cells by exact Fraction midpoints.  It reproduces (9)--(11), both
intervals in (10), the shared row (13), and the exact ledger identity in (7).
It explicitly leaves the full `3,539,936`-row quantifier to the primary.
Normal and optimized Python outputs are byte-identical.

Frozen SHA-256 values are

```text
C++ source     bc4b3d622e9b7b5c1906024e6f62b35ba4499edbacc37af23a801e5702e83af6
C++ output     ed26470874e6e57e19bdbcc93c8269848fd44ab35b21ce00863878d6c6afb8a4
Python source  1d207325f5e376641ea2a8c3aea531cff31bd66c6e69347477b7e6687ef9e25e
Python output  889d1bdd2ee772283170600df041575b16b96688351a9fededf40153ded7e233
```

## 4. Kakeya/tournament carrier

The all-scale proof is one compact Kakeya object split into three regimes:
finite bottom endpoint word, labelled polygon rectangle, and the exact cells
between them.  Runner, comb, core-cell, section-boundary, wall-event, residue,
and proof-obligation vertices were all challenged.

Exact boundary order gives a transitive tournament: no directed cycles,
singleton SCCs, and one sorted Hamiltonian path after ties are coalesced.  It
is not faithful.  It destroys cyclic gap size, endpoint owner, wall slope,
safe side, and removal stage.  The predicate-preserving carrier is

```text
(exact endpoint coordinate, owner, stage, adjacent safe side, wall slope). (14)
```

The metric sidecars in (14), not the naked tournament, carry (1).

## 5. Frontier after the theorem

THM-1133 closes every four-comb ray with offset span at most `30`.  THM-1128
closes an arbitrary-shape centred cone, with midpoint corollary

```text
K>=max(1272,26(c+1)).                                   (15)
```

Uniform `r=5` is therefore reduced to offset shapes with `c>30` lying outside
the THM-1128 cone and outside the other analytic gates.  No finite-shape
translation or monotonicity claim is made beyond the proved ranges above.
