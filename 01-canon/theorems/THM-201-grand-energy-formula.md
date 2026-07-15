# THM-201: Grand Fourier Energy Formula

**Status:** CORRECTED (kind-pasteur-S112; exact coefficient EGF proved in
THM-848).  The historical `m^k` numerator is exact only for `k<=2` or `m=1`.
**Session:** kind-pasteur-2026-03-14-S105--S107, corrected
kind-pasteur-2026-03-15-S112 and codex-2026-07-15-S15

## Statement

For the Hamiltonian path count `H(T)` on tournaments with `n` vertices, let
`E_(2k)` be its Walsh energy at level `2k`, let
`E_0=(n!/2^(n-1))^2`, and put `m=n-2k`.  Then

```text
E_(2k)/E_0 = 2g_k(m)/(n)_(2k)
            = [w^k]((1+w)/(1-w))^m/(n)_(2k),
```

for `1<=k<=floor((n-1)/2)`.  Equivalently,

```text
2g_k(m)
 = sum_(j=0)^min(m,k) binom(m,j)binom(m+k-j-1,k-j).
```

THM-217 gives the matching interpretation of `g_k`; THM-848 gives a direct
weighted even-path-forest EGF proof.  The first coefficient polynomials are

```text
2g_1(m)=2m,
2g_2(m)=2m^2,
2g_3(m)=2m(2m^2+1)/3.
```

The original proposal `2m^k/(n)_(2k)` is therefore exact for `k=1,2` and on
the boundary `m=1`, but it is not even the correct leading coefficient in
general: for fixed `k`,

```text
2g_k(m) ~ 2^k m^k/k!.
```

## Corollaries

**Variance:**

```text
Var(H)/Mean(H)^2
 = sum_(k=1)^floor((n-1)/2)
     [w^k]((1+w)/(1-w))^(n-2k)/(n)_(2k).
```

**Concentration:** `n Var(H)/Mean(H)^2 -> 2`; in particular the ratio is
`O(1/n)`.

**Spectral purification:** `E_2/Var(H) -> 1`.  Level two contributes exactly
`2(n-2)/(n(n-1))` after normalization by `E_0`.

**Fixed-level asymptotics:** for fixed `k`,
`E_(2k)/E_0 ~ 2^k/(k! n^k)`.

## Small-size verification and first refutation

| n | k | Exact `E_(2k)/E_0` | Note |
|---|---|---|---|
| 3 | 1 | 1/3 | exact |
| 4 | 1 | 1/3 | exact |
| 5 | 1,2 | 3/10, 1/60 | exact |
| 6 | 1,2 | 4/15, 1/45 | exact |
| 7 | 1,2,3 | 5/21, 3/140, 1/2520 | exact |
| 8 | 1,2,3 | 3/14, 2/105, 1/1680 | first interior `k=3` correction |

At `n=8,k=3`, the exact numerator is `[w^3]Q(w)^2=12`, while the old
numerator was `2(2)^3=16`.  Hence

```text
Var(H)/Mean(H)^2=131/560,
```

not `59/252`.  THM-589 independently obtains `131/560` from the exact
succession count `W(8)=49752`.

## Historical proof heuristic

The original spectator-freedom argument replaced the exact weighted matching
count `g_k(m)` by `m^k`.  This happens to preserve `k<=2` and `m=1`, explaining
all checks then available.  For `k>=3,m>=2`, component collisions change the
weights and the exact matching/forest polynomial must be retained.
