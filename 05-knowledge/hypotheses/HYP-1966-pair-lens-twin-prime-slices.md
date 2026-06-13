---
id: HYP-1966
status: OPEN
source: codex-2026-06-01-S495
related:
  - HYP-1931
  - HYP-1941
  - HYP-1953
  - HYP-1962
  - HYP-1963
  - HYP-1964
  - HYP-1965
  - HYP-1960
  - HYP-1961
---

# HYP-1966: twin primes are the base ray of the prime-pair plane

## Statement

Every unordered pair of natural numbers should first be written in midpoint
coordinates:

```text
(a,b) = (m-h, m+h).
```

This turns a surprising amount of the additive, tournament, and LRC material
into different slices of one pair plane.

- Twin primes are the fixed half-gap ray `h=1`, with `m` varying.
- Prime pairs of gap `2h` are fixed half-gap rays, with a local
  Hardy-Littlewood chord factor.
- Goldbach for `N` is the fixed midpoint row `m=N/2`, with `h` varying.
- Zeckendorf carry debt is a normal-form cost attached to a pair edge.
- Tournament Analysis is the complete choice of a binary relation on every
  unordered pair.

The core claim is methodological: before lifting to triples, hypergraphs, or
global certificates, compute the pair movie in coordinates

```text
(midpoint, half-gap, local residue product, normal-form carry debt).
```

Twin primes are the cleanest test case because the pair is forced:

```text
(m-1, m+1).
```

There is no summand choice left.  The whole problem becomes whether infinitely
many midpoints survive all local pair sieves.

## Evidence

`04-computation/pair_lens_twin_primes_s495.py` builds the pair-plane probe and
stores its output in `05-knowledge/results/pair_lens_twin_primes_s495.out`.

### Twin-prime ray

The script counts twin pairs up to `200000` and compares them to the finite
Euler-product approximation

```text
2*C_2*x/log(x)^2.
```

Selected output:

```text
limit      twins    2*C2*x/log^2x     ratio
10000        205             155.64    1.317
100000      1224             996.11    1.229
200000      2160            1772.39    1.219
```

This is not offered as new analytic evidence.  It is a calibration: the ray
`h=1` is the unboosted baseline because `h` has no odd prime factors.

### Pair chords

For gap `2h`, the Hardy-Littlewood prime-pair correction relative to twins is

```text
prod_{p|h, p odd} (p-1)/(p-2).
```

The midpoint explanation is local.  For each odd prime `q`, the pair
`m-h, m+h` forbids midpoint residues

```text
m = +h mod q,  m = -h mod q.
```

If `q` divides `h`, the two forbidden residues collapse to one, so the slice is
less obstructed and receives the factor `(q-1)/(q-2)`.

The S495 counts reproduce this qualitative chord pattern:

```text
h=1  gap=2   chord=1     ratio/twin=1.000
h=3  gap=6   chord=2     ratio/twin=1.988
h=5  gap=10  chord=4/3   ratio/twin=1.325
h=15 gap=30  chord=8/3   ratio/twin=2.679
```

### Midpoint spine

For `h=1`, the midpoint sieve forbids `m=+/-1 mod q` for each odd `q`.  The
first nontrivial case is `q=3`, so every twin pair beyond `(3,5)` has midpoint
divisible by `3`.  Parity forces the midpoint to be even.  Thus the baseline
ray lives on the familiar spine:

```text
(6k-1, 6k+1).
```

In pair language this is not a modular curiosity; it is the first local
projection of the pair ray.

### Goldbach row

Goldbach is perpendicular to twin primes in the same plane.  It fixes
`m=N/2` and scans half-gaps `h`.

Selected S495 rows:

```text
N=1000   m=500   Goldbach pairs=28   h range=9..497
N=5000   m=2500  Goldbach pairs=76   h range=117..2493
N=10000  m=5000  Goldbach pairs=127  h range=81..4941
```

Thus Goldbach asks for coverage of every large even midpoint row; twin primes
ask for persistence along the single smallest half-gap ray.

### Zeckendorf carry debt

S495 overlays Zeckendorf supports on twin pairs.  Among `705` twin pairs up to
`50000`, no zero-carry pair was found in this scan.  The first carry-score bins
begin at `4,7,8,9,...`.

The absence is useful: prime-pair survival and Fibonacci normal-form
compatibility are independent filters.  A pair can be locally prime-compatible
without being carry-compatible in a chosen additive normal form.

### Polygonal pair layer

Two `k`-gonal atoms cover only part of the target interval `1..300`:

```text
k=3: pair-covered 180, misses 120
k=4: pair-covered 113, misses 187
k=8: pair-covered  51, misses 249
```

Fermat polygonal is therefore the warning label for pair-first thinking:
pairs are the first relational layer, not always the final proof layer.
Higher arity repairs the missed pair layer.

## Consequences

1. Prime-pair studies should record both coordinates.  Fixed-gap rays and
   fixed-sum rows are dual probes of the same object.

2. Hardy-Littlewood singular factors are pair-local residue products.  In the
   pair plane, the factor is literally the collapse or separation of forbidden
   midpoint residues.

3. Bounded-gap theorems show that at least one bounded half-gap ray, or a
   bounded family of admissible rays, survives infinitely often.  They do not
   identify the `h=1` ray, but they are exactly the right kind of theorem for
   a pair-plane attack.

4. Tournament Analysis should orient candidate pairs by debt-reduction
   metrics: lower local obstruction product, lower carry debt, stronger
   coverage, better pressure-peel behavior, or greater endpoint export.

5. For LRC, the analogous pair object is not only nearest distance from the
   stationary runner.  It is the full pair movie of runner-runner distances,
   two-neighbor handoffs, endpoint protection pairs, and pressure dependencies.

## External anchors

- Hardy and Littlewood's 1923 prime-sums paper is the classical source for the
  singular-series viewpoint: <https://link.springer.com/article/10.1007/BF02403921>.
- Zhang proved the first finite bound for prime gaps in Annals of Mathematics:
  <https://annals.math.princeton.edu/2014/179-3/p07>.
- Maynard's refinement proved bounded `m`-prime gaps for every fixed `m` and
  improved the one-gap bound: <https://annals.math.princeton.edu/2015/181-1/p07>.
- The Polymath8 retrospective records the route from Zhang's bound through
  later improvements to `H_1 <= 246`: <https://arxiv.org/abs/1409.8361>.

## Next probes

- Build a prime-pair tournament on half-gaps `h`, orienting `h1 -> h2` when
  the `h1` ray has lower residual debt after normalizing by its
  Hardy-Littlewood chord.
- For Goldbach rows, orient prime pairs by a lexicographic debt vector
  `(local product, balance, Zeckendorf carry, gap position)` and track whether
  low-debt pairs form stable tournament kernels as `N` grows.
- For LRC, copy the same pair-plane schema: midpoint becomes time/phase
  midpoint, half-gap becomes circular separation, local product becomes
  p-adic endpoint debt, and carry becomes Zeckendorf or denominator-depth debt.

## See Also

- `04-computation/pair_lens_twin_primes_s495.py`
- `05-knowledge/results/pair_lens_twin_primes_s495.out`
- `07-reflections/pair-lens-twin-primes-s495.md`
- HYP-1953
- HYP-1962
- HYP-1963
- HYP-1964
- HYP-1965
- HYP-1931
- HYP-1941
