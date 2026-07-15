# THM-204: Grand Fourier Level Formula -- corrected boundary notice

**Status:** REFUTED AS AN ALL-n FORMULA at `n=8`; the exhaustive `n<=7`
computation remains valid.  The exact replacement is proved in THM-848.
**Source:** opus-2026-03-14-S89; corrected codex-2026-07-15-S15

## Historical statement

This theorem originally proposed, for the Walsh-Hadamard decomposition of the
Hamiltonian path count `H(T)`,

```text
E_(2k)/E_0 = 2(n-2k)^k/(n)_(2k),                        (old)
```

where `E_(2k)` is the sum of squared Walsh coefficients at level `2k`,
`E_0=(n!/2^(n-1))^2`, and `(n)_(2k)=n!/(n-2k)!`.  Exhaustive FWHT checks
established `(old)` for every level through `n=7`; they did not establish an
all-size law.

## Exact correction

Put `m=n-2k` and `Q(w)=(1+w)/(1-w)`.  THM-848 proves from the weighted
even-path-forest EGF that

```text
E_(2k)/E_0 = [w^k]Q(w)^m/(n)_(2k)
            = 2g_k(m)/(n)_(2k),                         (correct)
```

where `g_k` is THM-201/217's matching polynomial.  Explicitly,

```text
[w^k]Q(w)^m
 = sum_(j=0)^min(m,k) binom(m,j)binom(m+k-j-1,k-j).
```

The old numerator `2m^k` agrees with the exact coefficient for `k=1,2`, and
for every `k` when `m=1`.  Those two coincidences cover every level occurring
at `n<=7`.  At the first new interior case,

```text
n=8, k=3, m=2:
[w^3]Q(w)^2=12,       whereas       2m^3=16,
E_6/E_0=1/1680,       not           1/1260.
```

The corrected `n=8` variance ratio is

```text
3/14 + 2/105 + 1/1680 = 131/560,
```

not the historical prediction `59/252`.  This agrees independently with
THM-589's exact succession formula `W(8)/8!-1=49752/40320-1=131/560`.

## What survives

1. Odd Walsh levels vanish.
2. The level-two and level-four formulas in the original theorem are exact.
3. For odd `n`, the top level has `m=1`, so its original formula is exact.
4. The exhaustive `n<=7` data and the asymptotic statement that level two
   dominates are unaffected.
5. Every use of `2(n-2k)^k` at an interior level `k>=3,m>=2` must be replaced
   by `[w^k]Q(w)^m=2g_k(m)`.
