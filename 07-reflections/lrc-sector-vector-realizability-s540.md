---
source: codex-2026-06-01-S540
status: computation plus theorem-shaped clarification
tags:
  - lonely-runner
  - sector-vectors
  - realizability
  - occupancy
  - forced-hitting
  - tournament-analysis
---

# Realizable Sector-Vectors

HYP-2022 made the sectors of the `1/n` grid into tournament vertices.  The
next natural question is: which sector-vectors are actually realizable?

This session's main answer is deflationary, but useful:

```text
Every raw sector-vector is existentially realizable.
```

So the word "realizable" has to be handled carefully.  There are at least
three different notions:

```text
existentially realizable:
  some primitive speed set and some open time exhibit c;

clock-reachable:
  a fixed speed set V exhibits c somewhere along its orbit;

forced:
  every primitive speed set exhibits c, or every primitive speed set hits a
  target family of c's.
```

LRC is about the third notion, not the first.  It asks whether every primitive
clock hits the good face:

```text
c_0 = c_{n-1} = 0.
```

## Computation

Artifact:

```text
04-computation/lrc_sector_vector_realizability_s540.py
05-knowledge/results/lrc_sector_vector_realizability_s540.out
```

The script enumerates all compositions of `n-1` into `n` parts, constructs
explicit primitive witnesses for every vector, and performs bounded exact
open-cell searches.

Raw counts:

```text
n=3: all vectors=6,   good vectors=1
n=4: all vectors=20,  good vectors=4
n=5: all vectors=70,  good vectors=15
n=6: all vectors=252, good vectors=56
n=7: all vectors=924, good vectors=210
```

These are:

```text
all vectors  = C(2n-2,n-1)
good vectors = C(2n-4,n-3)
```

Observer-clasp reflection orbits and free dihedral orbits:

```text
n=3: clasp=4,   dihedral=2
n=4: clasp=10,  dihedral=4
n=5: clasp=38,  dihedral=10
n=6: clasp=126, dihedral=26
n=7: clasp=472, dihedral=76
```

Constructive witness audit:

```text
n=3: 6/6 realized,     max witness speed <= 12
n=4: 20/20 realized,   max witness speed <= 21
n=5: 70/70 realized,   max witness speed <= 32
n=6: 252/252 realized, max witness speed <= 45
n=7: 924/924 realized, max witness speed <= 60
```

The max speeds follow `(n-1)(n+3)` in this construction.

Bounded exact searches:

```text
n=3, B<=10: vectors seen=6/6,     good seen=1/1
n=4, B<=12: vectors seen=20/20,   good seen=4/4
n=5, B<=14: vectors seen=70/70,   good seen=15/15
n=6, B<=16: vectors seen=252/252, good seen=56/56
n=7, B<=14: vectors seen=897/924, good seen=202/210
```

So low-complexity clocks already see all raw vectors through `n=6`, and most
through `n=7`.  This reinforces the main point: the globally possible
sector-vector set is not restrictive.

## Constructive Proof

Let `c=(c_0,...,c_{n-1})` be any composition of `n-1`.  Choose `L` larger than
every `c_k`, and put `q=nL`.  Sector `k` contains the open integer residues:

```text
kL + 1, kL + 2, ..., (k+1)L - 1.
```

Choose `c_k` distinct residues from this interval for each `k`, with total
gcd `1`; this is possible for large enough `L` because each interval has many
choices, and the total count is finite.  Use those residues as speeds and set:

```text
t = 1/q.
```

Then every speed `r` lands at `r/q`, strictly inside its assigned sector.  The
speeds are distinct, primitive, and the sector-vector is exactly `c`.

Thus:

```text
global sector-vector existence = the full composition simplex.
```

## Forced Vectors Are The Wrong Ones

The bounded searches also compute the intersection of all vectors seen by all
primitive speed sets in the box.  For `n>=5`, the forced intersection is the
same four bad boundary fans:

```text
(n-1,0,0,...,0)
(n-2,1,0,...,0)
(0,...,0,1,n-2)
(0,...,0,n-1)
```

For example:

```text
n=5:
  (4,0,0,0,0), (3,1,0,0,0), (0,0,0,1,3), (0,0,0,0,4)

n=6:
  (5,0,0,0,0,0), (4,1,0,0,0,0),
  (0,0,0,0,1,4), (0,0,0,0,0,5)

n=7:
  (6,0,0,0,0,0,0), (5,1,0,0,0,0,0),
  (0,0,0,0,0,1,5), (0,0,0,0,0,0,6)
```

These are forced because near `t=0` all runners lie in sector `0`, then the
fastest runner leaks into sector `1`; near `t=1` the symmetric sector
`n-1` versions occur.  They are not lonely.  The forced intersection has:

```text
good=0
```

This is a crucial correction to the naive sector-vector program:

```text
There is no single good sector-vector forced by every clock.
LRC is forced hitting of a good face, not forced hitting of a point.
```

## What HYP-2022's Restriction Really Means

HYP-2022's sector-tournament restriction cannot mean "only some
sector-vectors exist globally."  They all exist.

The restriction lives in one of four sharper places:

1. **Fixed-clock menus.**  For a fixed speed set, the path through sector
   vectors is highly constrained.
2. **Low-complexity stratification.**  Some vectors require larger max speed
   or more arithmetic complexity to realize.
3. **Quotient restriction.**  Sector-tournament classes collapse the full
   vector simplex into a small image, and that image is far from all
   tournaments.
4. **Forced target hitting.**  LRC asks every primitive clock to hit the good
   face `c_0=c_{n-1}=0`, with boundary compactification for AP tightness.

The sector-vector simplex is the ambient lattice.  LRC is not about whether the
good face exists; it is about whether every arithmetic line crosses it.

Concurrent upstream work sharpened this interpretation while this session was
closing.  HYP-2024's section/boundary functors are pure because they keep the
observer target visible inside this full sector-vector simplex.  HYP-2026's
support-flow/zero-cut language is the dual menu statement: the question is
whether each fixed clock creates a zero-cover cut, not whether zero-cover sector
vectors exist somewhere in the global universe.

## Implications

### 1. The Good Face Is A Face, Not A Class

The good sector-vectors form a face of dimension `n-3`:

```text
c_0=c_{n-1}=0,  sum interior c_i = n-1.
```

Counting good vectors gives:

```text
C(2n-4,n-3).
```

That is already a large target family.  Any proof language that tries to name
one favorite good vector is probably too rigid.  The target should be an ideal,
face, or absorbing set.

### 2. The Boundary Fans Are Universal But Bad

Every clock begins and ends in the bad observer sectors.  This suggests a
flow/hitting theorem:

```text
every primitive sector-vector path starts in the bad left boundary fan,
ends in the bad right boundary fan,
and must cross the good face or only touch its compactified AP boundary.
```

That is a more geometric statement than "some sector vector is forced."

### 3. The Sector Program Should Study Menus, Not Existence

The next invariant should be the menu:

```text
M(V) = { sector vectors reached by V on open cells }.
```

Useful fingerprints:

```text
|M(V)|,
|M(V) ∩ Good|,
conv(M(V)) ∩ Good,
first-hit face,
missing-good-vector profile,
boundary-only good profile,
sector-tournament image of M(V).
```

The bounded search already shows AP-like tight rows are the only uncertified
open-cell misses in small boxes.  The menu perspective should separate:

```text
open good hit
boundary-only good hit
no good hit in sampled grid
```

### 4. Tournament Mapping Should Use Good-Face Ideals

For sector tournaments, the target should not be a single isomorphism class.
It should be:

```text
classes whose color/anchor data certifies c_0=c_{n-1}=0
```

or, in HYP-2023 language:

```text
anchored two-hole classes.
```

This is exactly the lesson `purity = compression + the right anchor`.

## HYP-2028

**All sector-vectors are existentially realizable.  The LRC sector-vector
problem is forced hitting of the observer-empty face by fixed primitive clock
menus, not classification of the globally realizable vectors.**

Predictions:

1. For every `n`, every composition of `n-1` into `n` sector counts has a
   primitive distinct-speed witness at an open time.
2. For `n>=5`, the only sector-vectors forced by every primitive speed set are
   the four bad boundary fans above.
3. The AP initial segment is the unique family whose menu hits the good face
   only on the compactified boundary.
4. The correct tournament target for sectors is a good-face ideal in anchored
   sector/hole classes, not a singleton class.

The slogan:

```text
global realizability is full;
forced realizability is the conjecture.
```
