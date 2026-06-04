# Computational speedups for LRC and unit distances (S628)

Engineering deliverable: turning the structural understanding (THM-411/412/415/418, HYP-2230/2245)
into measured algorithmic speedups. All correctness-checked against brute references; all code in
`04-computation/lrc_fast_s628.py` and `04-computation/lrc_unit_speedup_bench_s628.py`, benchmarks in
`05-knowledge/results/lrc_fast_s628.out`, `lrc_unit_speedup_bench_s628.out`.

## Headline: tight-instance enumeration — **100–660× faster, identical output**

| n | range | brute | fast | speedup | same result |
|---|---|---|---|---|---|
| 5 | 14 | 8.58 s | 0.02 s | **×388** | yes (2 tight) |
| 6 | 16 | 58.65 s | 0.09 s | **×661** | yes (1 tight) |

This is the exact task that kept timing out in S621–S623. Two structural levers:

1. **`gap_shells` filter (the big one, ×56 on its own).** The loneliness gap `M(S)` is, for tight and
   near-tight configs, attained on the **`(ℤ/m)*` witness orbit at shells `m ≤ 2n−1`** (THM-411/415).
   So instead of the brute critical-time set (all corners `(2k+1)/(2v_i)` and crossings
   `k/(v_i±v_j)`, `O(n²·v_max)` candidates), evaluate only `t = a/m`, `m ≤ 2n−1`, `a ∈ (ℤ/m)*` —
   `O(Σ_{m≤2n-1} φ(m)·n)` candidates, integer-only. **Exact when the optimum is on the orbit;
   a guaranteed lower bound otherwise** (174/400 exact, 400/400 `≤` brute). As a *filter*: if
   `gap_shells(S) > target` the config is provably loose and rejected with **no brute call** — and
   that rejects the vast majority instantly. Brute runs only on the high-shell residual.
2. **`v_max ≤ 2n−1` pruning** (THM-411 finiteness): the search box for tight instances is bounded;
   skip any candidate with a speed above `2n−1`. (Empirically safe: `same=True`.)

## LRC primitives (`lrc_fast_s628.py`)

- `gap_shells(S, n, Mmax=2n-1)` — fast exact-or-lower-bound gap on the shell witnesses. **×56** vs
  `gap_brute`. Use as the default; fall back to `gap_brute` only when an exact value off the orbit is
  needed.
- `good_mult_fast(S, p)` — O(n) good multiplier at a prime shell: the band `{0,±1}` forbids exactly
  `a = ±v_i⁻¹` (THM-418), so the bad set is built directly instead of scanning all `p−1` multipliers.
  *Asymptotic* win (large `p`, or near-tight configs where good multipliers are rare); at the small
  `p ≈ 2n` we usually run, brute short-circuits and the two are comparable — use `good_mult_fast`
  when `p` is large or a guaranteed-`O(n)` worst case matters.
- `is_loose_fast(S, n)` — O(n) lower bound `M ≥ 2/p` (prime-shell dodge, THM-418). Fast certificate
  of *near*-looseness; note `2/p ≲ 1/(2n)` so it does not by itself certify `M > 1/n` (the factor-2
  gap is the ramified-shell frontier, HYP-2240).

## Unit-distance primitives (`lrc_unit_speedup_bench_s628.py`)

- `tri_count_fast(ij)` — O(k) unit-distance count on the triangular (= Eisenstein CM) lattice via
  the neighbor rule `di²+di·dj+dj²=1` (6 neighbors), vs O(k²) all-pairs. **×3.7** at k=22, grows
  with k. Use integer lattice coordinates, not floats.
- **CM-rotation restriction** (HYP-2230): the rotations that create extra unit distances are the
  bounded-height modulus-1 elements `α/ᾱ` of `ℚ(ω)`; for `|height| ≤ 4` there are **24** of them.
  A rotation-overlay search ranges over these 24 angles, not a continuum — making the
  "lattice + rotations" optimum search (the n=22 frontier) a finite, small problem.

## When to use what

- Enumerating tight / collapse-family / extremal configs → `gap_shells` filter + `v_max` prune
  (100–660×).
- Deciding looseness with a fast lower bound → `is_loose_fast` (O(n)).
- Exact `M` off the witness orbit (the residual, high-shell optima) → `gap_brute` (fallback).
- Unit distances on a lattice → `tri_count_fast` (O(k)); searching rotated overlays → the 24 CM
  rotations.

The throughline: **structure shrinks the search.** Every speedup replaces an exhaustive geometric
search with evaluation on the arithmetic object our theory says the answer lives on — the `(ℤ/m)*`
shell orbit for LRC, the CM-rotation set for unit distances.
