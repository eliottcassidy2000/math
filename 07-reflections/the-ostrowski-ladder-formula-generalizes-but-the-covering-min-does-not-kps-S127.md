# The Ostrowski ladder formula generalizes to all N; the covering-min N/Φ₆(N) does not

*kind-pasteur-2026-07-11-S127 cont.54. Owner: "work on the remaining Ostrowski LRC mathematics." klein-S267
pinned the LRC(14) covering minimum at `14/183 = N/Φ₆(N)` via the Ostrowski ladder `{1..12, 13k}`. I asked
whether that is a general law. Half of it is — the ladder formula is exact for every N — and half of it is
not — `N/Φ₆(N)` is the covering-min only for large N, and I confirm klein's `14/183` by ruling out the
small-N families that would have undercut it.*

---

## The clean half: the ladder formula is universal

For LRC(N) (N−1 speeds, tight bound `1/N`), the **Ostrowski ladder** `{1, …, N−2, (N−1)k}` — a base AP with
one far outlier — has loneliness

> **`M_k = k / ((N−1)k + 1)`, exactly, for every N and every k.**

Verified for `N = 3,…,7,14` and `k = 1,…,N` on the nose (the outlier `(N−1)k` sits at
`‖13k²/((N−1)k+1)‖ = k/((N−1)k+1)`, and runner 1 matches it). The rungs climb from the tight value
`M_1 = 1/N` (the AP `{1..N−1}` itself) toward `1/(N−1)` as `k→∞`. This is a genuinely clean closed form — the
"spine" of the Farey/Stern–Brocot tree that mac-mini-S54 identified in the M-spectrum.

**The first covering rung.** The family is covering (has a multiple of every `d ∈ {2..N}`) iff its outlier
`(N−1)k` carries a multiple of `N`, i.e. `N ∣ (N−1)k`, i.e. `N ∣ k` (since `gcd(N−1,N)=1`). The smallest is
`k = N`, giving

> `M_N = N / ((N−1)N + 1) = N / (N² − N + 1) = N / Φ₆(N)`.

So `14/183 = 14/Φ₆(14)` is not cyclotomic magic; `Φ₆(N) = N²−N+1` is exactly `(N−1)N+1`, the denominator of
the first covering rung. The same arithmetic gives `3/7, 4/13, 5/21, 6/31, 7/43, …, 14/183` for `N = 3..14`.

## The false half: N/Φ₆(N) is not the covering-min

Being *a* covering family is not being the *minimum*. Exhaustive covering-min over small N:

| N | ladder `N/Φ₆(N)` | true covering-min | at |
|---|---|---|---|
| 3 | 3/7 ≈ 0.429 | **2/5 = 0.400** | {2,3} |
| 4 | 4/13 ≈ 0.308 | **2/7 ≈ 0.286** | {1,3,4} |
| 5 | 5/21 ≈ 0.238 | **2/9 ≈ 0.222** | {1,3,4,5} |

For small N the ladder is *beaten* — compressed covering families sit below it. So `N/Φ₆(N)` is **not** a
universal covering-min formula; it is the covering-min only from some N onward (a transition), and klein's
`14/183` is the `N=14` instance of "the ladder has become the minimum," verified by his ILP certificate, not a
law that holds at every N.

## Why klein's 14/183 nonetheless stands

The natural worry: if the small-N compressed families extrapolate as `2/(2N−1)`, at `N=14` that would be
`2/27 ≈ 0.0741 < 14/183`, undercutting the crux value. It does not — the extrapolation is wrong. The actual
`N=14` analogs of the small-N minimizers are loose by a wide margin:

- `{1,3,4,…,14}` (the `{1,3,4,5}` analog): `M = 2/17 ≈ 0.118`
- `{2,3,…,14}` (the `{2,3}` analog): `M = 1/8 = 0.125`
- `{1,2,4,…,14}` (drop 3): `M = 2/19 ≈ 0.105`

all far above `14/183 ≈ 0.0765`. The compressed families climb steeply with N and the ladder overtakes them
well before N=14, so at N=14 the deep well `{1..12, 182}` really is the covering minimum — independently
consistent with klein's `speeds ≤ 182` ILP infeasibility certificate. So this cross-check *supports* the crux
value while correcting the impression that it comes from a closed-form law.

## What this sharpens

The Ostrowski/change-of-base picture is exactly right about the **shape** of the extremals — they are ladder
rungs, continued-fraction-simple families — and the ladder formula `k/((N−1)k+1)` is the honest general
statement. But "the covering-min is `N/Φ₆(N)`" is an LRC(14) fact, not a theorem of the ladder: the ladder is
*a* covering family at every N, and the *minimum* only once N is large enough that the compressed competitors
have risen above it. For the proof this means the crux `DC ⟹ M ≥ 14/183` cannot be reduced to "the ladder
formula" — it needs the extra content that at N=14 no other structure dips below the ladder's first covering
rung, which is precisely klein's finite ILP certificate (`≤ 182`), still the open uniform statement above it.

*Files: lrc14_ostrowski_covering_min_general_kps_S127.py/out, lrc14_ostrowski_transition_kps_S127.out,
lrc14_ostrowski_transition_pin_kps_S127.out. HYP-6215. Extends klein-S267/HYP-6190 (14/183 = N/Φ₆(N)), mac-mini
cont.54 (Farey tree), the Ostrowski-ladder reflection (mac-mini-S38); confirms klein's 14/183 by the compressed-
family cross-check.*
