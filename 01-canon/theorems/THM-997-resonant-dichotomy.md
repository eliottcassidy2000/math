# THM-997 — The resonant dichotomy: units live, non-units deep (boxeph-2026-07-17-S81)

**Status:** PROVED (elementary; verified N=3..22 — `lrc_prefix_census_boxeph_S81.py`).
The `q = N` slice of the uniform live law [[THM-996-resonance-confinement-census-law-deathstar-S56]] (death-star; its Part II is this live law), stated as a
perfect partition. Generalizes the N=14 observation (death-star: at q=14, live = the
6 units, deep = {p=7} and the evens) to a clean ring-theoretic law at every N.

## Statement

Fix `N ≥ 2`, family `V = {1,…,N-1}`, threshold `1/N`, and take the **base resonant
modulus `q = N`**. Every multiplier `p ∈ {1,…,N-1}` is exactly one of:

- **live** — iff `p ∈ (ℤ/N)^×` (a unit): then `{ip mod N} = {1,…,N-1}`, so
  `min_i ‖ip/N‖ = 1/N` (the equality instants). There are `φ(N)` of these.
- **deep (collision)** — iff `p` is a **zero-divisor** (`gcd(p,N) = g > 1`): the runner
  `i = N/g ∈ {1,…,N-1}` satisfies `N ∣ ip`, so runner `i` sits **exactly on the
  origin**. There are `N − 1 − φ(N)` of these.

> At `q = N`, `{1,…,N-1} = (ℤ/N)^× ⊔ (zero-divisors)` = `live ⊔ deep`, with **nothing
> in between**: `#live = φ(N)`, `#deep = N − 1 − φ(N)`.

## Proof

`q = N ⟹ c = ⌈N/N⌉ = 1`, so `p` live `⟺ ip mod N ∈ [1, N-1]` for all `i` `⟺ N ∤ ip`
for all `i ∈ {1,…,N-1}` `⟺ gcd(p,N) = 1`. If `gcd(p,N) = g > 1`, then `i = N/g` is a
valid speed and `ip = p(N/g) = (p/g)N ≡ 0`, a collision (deep). `live ∩ deep = ∅` since
a live `p` has no runner within `1/N`, hence none at 0. The two counts sum to `N-1`. ∎

## The arithmetic fabric (N-dependence)

| regime | `#deep(q=N) = N−1−φ(N)` | reading |
|---|---|---|
| `N` prime | `0` | **every** multiplier live — LRC(prime) equality case is maximally lonely at resonance |
| `N` prime power `p^k` | `p^{k-1}−1` | one collision ray `p^{k-1}` and its multiples |
| `N` highly composite | large | many collision rays; e.g. N=14: 7 deep, N=12: 7, N=18: 11 |

The "difficulty" of the LRC(N) equality-case census is measured by the **zero-divisor
deficit `N−1−φ(N)`**: prime N has a trivial census (all live), composite N loses exactly
the non-units to origin collisions. N=14's `7 = 13 − 6` deep multipliers are the
zero-divisors of `ℤ/14` (evens + 7).

## Significance

- Recasts the census as **`ℤ/N`-arithmetic**: live = unit group, deep = zero-divisors.
  The live/deep race that `lonely_of_census` adjudicates is `φ(N)` vs `N−1−φ(N)`.
- Explains uniformly why death-star's `q=14` census closes by `decide`: the partition is
  forced, not incidental.
- Feeds the "generic wagner-circle theorem" — the per-family circle list is the family's
  resonance-denominator profile; the tight family's is maximally simple.

Related: [[THM-996-resonance-confinement-census-law-deathstar-S56]] (death-star; its Part II is this live law), [[THM-998-farey-circle-deep-law]], [[THM-991]].
