# Variable: good-cut bucket polynomial

**Symbol:** `B_N(x)`, with connected-run counts `c_L`
**Type:** ordinary generating polynomial
**Defined in:** THM-349

## Definition

Let `N=n-1` be the number of legal cuts in the fixed-base staircase. Let

```text
B_N(x) = sum_g b(N,g) x^g
```

where `b(N,g)` is the number of tilings whose good-cut set has size `g`.

Let `c_L` be the number of interval subsets covering a connected run of `L`
cuts, using only intervals of length at least `2`.

## Equations

Connected-run inclusion-exclusion:

```text
c_L = sum_{A subset [L]} (-1)^|A| 2^{sum_R binom(|R|,2)}
```

where `R` ranges over the contiguous runs of `[L] \ A`.

Bucket recurrence:

```text
B_N(x) = B_{N-1}(x) + sum_{L=2}^N c_L x^L B_{N-L-1}(x)
```

with `B_{-1}(x)=B_0(x)=1`.

For a specific good-cut set with connected run lengths `L_1,...,L_r`, the
number of tilings realizing it is

```text
prod_i c_{L_i}.
```

## Values

```text
c_2=1
c_3=5
c_4=50
c_5=903
c_6=30773
c_7=2032504
c_8=264271477
```

Bucket polynomials:

```text
B_5 = 1 + 4x^2 + 15x^3 + 101x^4 + 903x^5
B_6 = 1 + 5x^2 + 20x^3 + 153x^4 + 1816x^5 + 30773x^6
B_7 = 1 + 6x^2 + 25x^3 + 206x^4 + 2739x^5 + 61671x^6 + 2032504x^7
```

## Relationships

- Refines `g(τ)`: `b(N,g)` is the fiber size of the good-cut count over the
  full tiling cube.
- Recovers THM-336 small-bucket formulas:
  `[x^2]B_{n-1}=n-2`, `[x^3]B_{n-1}=5(n-3)`, and
  `[x^4]B_{n-1}=50(n-4)+binom(n-4,2)`.
- Gives a cheap engineering checksum for tiling-cube quotient computations.

## Sources

- Theorem: `01-canon/theorems/THM-349-good-cut-interval-union-recurrence.md`
- Script: `04-computation/goodcut_interval_union_s15.py`
- Results: `05-knowledge/results/goodcut_interval_union_s15.out`
