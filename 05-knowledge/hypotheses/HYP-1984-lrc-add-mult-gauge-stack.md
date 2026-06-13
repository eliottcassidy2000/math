---
id: HYP-1984
status: OPEN
source: codex-2026-06-01-S513
related:
  - HYP-1963
  - HYP-1965
  - HYP-1966
  - HYP-1972
  - HYP-1975
  - HYP-1976
  - HYP-1977
  - HYP-1978
  - HYP-1979
  - HYP-1980
  - HYP-1981
  - HYP-1982
  - HYP-1983
---

# HYP-1984: LRC loneliness is an add/multiply gauge stack over the odd-core grid

## Statement

The LRC metric that connects `H`, score sequences, A000568, Goldbach/Lemoine
fibers, and product-sum equations is not a single arc rule.  It is a stack of
gauges over the natural-number address

```text
N = 2^h * odd_core.
```

The horizontal coordinate is the odd-core chain

```text
odd_core -> odd_core + 2,
```

and the vertical coordinate is dyadic doubling

```text
N -> 2N.
```

On that grid:

```text
even N:  additive fiber = Goldbach pairs p + q = N,
odd N:   additive fiber = Levy/Lemoine pairs p + 2q = N,
product/multiply labels = divisors, dyadic height, largest proper divisor,
product-sum collisions = P - S = D and (a-1)(b-1) = N-1,
LRC labels = endpoint debt phi(N), threshold pressure, and marked fibers.
```

Addition is the flattening operation: as an unlabeled shadow it becomes the
total order, but its fibers remember representation abundance.  Multiplication
is the branching operation: divisibility, odd core, and dyadic height form the
sparse skeleton.  Product-sum equations are collision equations between these
two operations.  LRC needs the collision stack because loneliness is controlled
by where additive smoothing, multiplicative branch depth, endpoint debt, and
marked tournament quotient data disagree.

## Evidence

`lrc_add_mult_gauge_stack_s513.py` builds denominator tournaments on

```text
10 12 14 15 16 18 20 21 22 24 26 27 28 30
```

using explicit Tournament Analysis declarations.

Scalar operation gauges are ledgers, not loneliness metrics.  Every scalar
denominator gauge completed to a transitive tournament:

```text
additive_prime_fiber          H=1, c3=0
additive_scarcity             H=1, c3=0
multiplicative_branch_depth   H=1, c3=0
product_sum_collision         H=1, c3=0
lrc_endpoint_debt             H=1, c3=0
a000568_odd_survival          H=1, c3=0
scalar_grid_stack             H=1, c3=0
```

The nontrivial shape appears only when criteria are combined pairwise by
majority vote:

```text
add_mult_majority:
  H=242407, c3=42, SCCs=(13,1)

loneliness_pressure_majority:
  H=509, c3=15, SCCs=(9,1,1,1,1,1)
```

This supports the stack thesis: individual arithmetic labels are static, but
their conflicts create a useful Tournament Analysis fingerprint.

The edge-flip table shows which coordinates are almost the same and which are
nearly perpendicular:

```text
additive_prime_fiber vs lrc_endpoint_debt:      0.10 flip rate
additive_prime_fiber vs a000568_odd_survival:   0.77 flip rate
multiplicative_branch_depth vs product_sum:     0.69 flip rate
scalar_grid_stack vs loneliness_pressure:       0.16 flip rate
```

So endpoint debt currently shadows additive abundance on the selected
denominators, while A000568 odd survival and product-sum collision see very
different axes.

The operation grid also explains why `14` and `18` are natural hard-row
neighbors:

```text
14 = 2 * 7  lives on first-even row h=1, odd core 7.
18 = 2 * 9  lives on first-even row h=1, odd core 9.
16 = 2^4    is the pure dyadic column.
```

They are adjacent first-even children in the `x+2` horizontal direction, while
dyadic ladders move vertically by `x*2`.

## Predictions

1. Corridor-level LRC scans should report three layers at once:

```text
runner dynamic gauges
+ pair-cell danger/operation labels
+ denominator add/multiply/A000568 labels.
```

2. Static arithmetic labels will separate denominator families better than
time cells inside one family; dynamic endpoint and danger gauges will separate
time cells better than static labels.

3. A counterexample-shaped row must align several bad coordinates at once:

```text
low additive representation abundance,
high endpoint debt,
unhelpful product-sum collision,
marked A000568 projection trap,
and a nonpeeling pressure SCC.
```

4. The best scalar `H` reading is not the proof object.  `H` should be read
as a height inside a chosen gauge layer, while score sequence, directed
triangles, SCCs, edge flips, and marked fibers tell which layer is carrying
the obstruction.

5. Product-sum equations should become useful in LRC when attached to
pair-cells or denominator transitions, not as standalone number identities.

## Proof Program

1. For each denominator `N`, attach the operation address
   `(odd_core(N), v2(N))` plus additive, multiplicative, product-sum, and
   endpoint-debt fibers.
2. For each LRC corridor, compute the runner and pair-cell gauge vectors from
   HYP-1975/HYP-1976 at every half-turn and endpoint wall.
3. Build a majority or sheaf-style tournament whose arcs compare full stack
   labels rather than scalar ranks.
4. Search for invariant features of hard rows: stable SCC cores, stable edge
   flip rates, and denominator transitions that preserve `gap*endpoint_debt`
   while changing dyadic height.
5. Translate any surviving stack obstruction into a marked A000568 quotient
   predicate from HYP-1977/HYP-1978.

## Sources

- `04-computation/lrc_add_mult_gauge_stack_s513.py`
- `05-knowledge/results/lrc_add_mult_gauge_stack_s513.out`
- `07-reflections/lrc-add-mult-gauge-stack-s513.md`
- HYP-1975
- HYP-1976
- HYP-1977
- HYP-1978
