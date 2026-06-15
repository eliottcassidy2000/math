# HYP-2523: A000568 delegates to odd-divisor profiles, but Python is partition-bound

**Status:** CONFIRMED as an exact structural reformulation and NEGATIVE as a
standalone Python speedup (`a000568_divisor_profile_codex.py`, codex-2026-06-15).

## Claim

For tournament Burnside counting, the surviving combinatorics after the odd-cycle
filter can be rewritten almost entirely in number-theoretic language:

```text
gcd(a,b) = Σ_{d|a,b} φ(d),
```

so the exponent of `2` in the A000568 Burnside term can be accumulated from the
odd-divisor profile

```text
S_d = #{chosen odd cycle blocks with length divisible by d}
```

instead of explicit pairwise `gcd` loops across distinct part sizes.

However, in the current Python regime this does **not** produce a net wall-clock
speedup, because odd-part partitions already have very small distinct-size support.
The true bottleneck is partition enumeration itself.

## Exact reformulation

For an odd-part cycle type with multiplicities `m_k`, the usual Burnside exponent
is

```text
e(λ) = Σ_k m_k (k-1)/2 + Σ_k C(m_k,2) k + Σ_{r<s} m_r m_s gcd(r,s).
```

Using `gcd(r,s)=Σ_{d|r,s} φ(d)`, if we build the partition in descending odd part
sizes and maintain

```text
S_d = Σ_{d|k already used} m_k,
```

then adding `m` copies of a new odd part `k` contributes

```text
self  = m (k-1)/2 + C(m,2) k,
cross = m Σ_{d|k} φ(d) S_d.
```

So the cross-size interaction is fully delegated to totients and divisibility.

## Exact evidence

`04-computation/a000568_divisor_profile_codex.py` verifies the divisor-profile
engine through `n=14` against the OEIS values:

```text
1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168,
48542114686912,28401423719122304.
```

It also matches an internal pairwise-gcd baseline exactly through all benchmarked
sizes.

## Benchmark result

The arithmetic reformulation is exact, but the Python timings are slightly worse:

```text
n   pairwise(s)  divisor-profile(s)
20    0.0005         0.0006
30    0.0025         0.0030
40    0.0080         0.0086
50    0.0253         0.0322
```

This is not a refutation of the reformulation. It identifies the wrong target for
optimization.

## Why the speedup does not appear in Python

The script measures the odd-part partition geometry:

```text
n=50:  odd partitions = 3658,   avg distinct sizes = 4.053
n=100: odd partitions = 444793, avg distinct sizes = 5.666
```

So the baseline pairwise work per partition is already tiny:

```text
n=100: avg pair gcd slots  = 13.752
       avg odd-div slots   = 12.378
```

The arithmetic recurrence does not shrink the dominant combinatorial layer,
namely the number of odd partitions. It only rewrites the per-partition local
interaction.

## Interpretation

This is the honest boundary of "delegate to pure number theory":

- **yes:** odd-cycle filter, totients, divisor profiles, multiplicative strata;
- **no:** that alone does not defeat the partition-count growth.

So the right next step is not another Python micro-optimization. It is to feed
the divisor-profile recurrence into the **compiled** A000568 engines, where the
arithmetic simplification can matter, or to find a second compression that reduces
the partition state space itself.

## Consequences for future work

1. The divisor-profile recurrence is the right arithmetic core for CRT / GMP /
   C implementations of A000568.
2. A genuine asymptotic speedup needs a second quotient:
   profile memoization, divisor-stratum recurrence, or another compression of
   odd-part partitions.
3. The FKN / transitive-corner picture suggests where such a quotient may come
   from conceptually: many labeled perturbations collapse onto a much smaller
   interaction profile, but the right quotient has not yet been found for the
   full Burnside sum.
