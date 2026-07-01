# A BRAVE attack on the tight-lowness ("tight ⟹ S is the AP minus O(1) elements") via covering-preservation: the chain "tight ⟹ covering (a multiple of every q, THM-523) ⟹ swapping q out needs a multiple mq that PATCHES the hole H_q (Jacobsthal-gated) ⟹ the gate is closed for LARGE PRIMES (they can't be swapped, census-confirmed) ⟹ few k swappable ⟹ lowness" is STRUCTURALLY right and gave a clean new prediction — the doubling family k=n−2 → 2(n−2) is exactly the runners-n ≡ 2 mod 6 case (= repo's speeds ≡ 1 mod 6), NEWLY verified tight at n=20 (18→36), extending the repo's n=14 — BUT my exact gate formula "Jacobsthal run ≥ n−k+2" is REFUTED (the n=8,14,20 doublings all have run J=3 < n−k+2=4, yet are tight), so the true gate is J(k) ≥ 3-ish (a run of 3 nonunits, ⇔ 6|k for k=n−2), the repo's HYP-2893; the full lowness reduces to bounding the count of Jacobsthal-swappable k — the repo's open crux (HYP-3740), which I confirm and sharpen but do not close

*opus-2026-06-30. Owner: throw patch-tuning at the lowness bravely. Result: the structure holds and extends
(n=20), one sub-claim (large primes forced) is clean, my exact Jacobsthal formula was wrong (corrected), and
the full lowness reduces to the repo's Jacobsthal counting.*

## The brave chain
1. **Tight ⟹ covering** — a multiple of every `q∈{2,…,n−1}` (q-witness at `t=1/q`, THM-523). Rigorous.
2. **Swapping `q` out ⟹ a multiple `mq` (m≥2) must be in S** (to still cover `q`), and it must **patch the
   hole `H_q`** left by removing runner `q` (patch-tuning). For `q∈(n/2,n−1]` the only candidates are `q`
   itself or `2q` (`3q>2n−2`).
3. **The patch `2q` is Jacobsthal-gated** — `2q` covers `H_q`'s inner half (`‖qt‖<1/(2n)`); the outer residual
   must be covered by the surviving AP runners, which needs a run of consecutive nonunits mod `q`.
4. **Large primes are forced.** For `q` prime, the nonunit run mod `q` is minimal (`J(q)=1`: only `0`), so `2q`
   cannot cover the residual ⇒ `q` cannot be swapped ⇒ **`q∈S` forced**. Census-confirmed: no swapped-out `k`
   is a large prime (`>n/2`); the swaps are `k=n−2` (composite) or the small prime `k=2` (cross-type).
5. **Few `k` swappable ⟹ lowness** — S is the AP minus the (few) swappable elements.

## What HELD (and a new result)
- **Large primes `(n/2,n−1]` are forced** — clean, census-verified through n=16.
- **The doubling family is `n_runners ≡ 2 mod 6`** (`k=n−2 → 2(n−2)`): verified tight at **n=8, 14, and now
  n=20** (`18→36`, `M=1/20` exactly) — a NEW extension of the repo's n=14. This equals the repo's
  "speeds ≡ 1 mod 6" law (`n_speeds = n_runners−1`), confirming HYP-2893 and adding a data point.
- **The single-swap tight locus is sparse**: `k=n−2` doublings (only when `6|(n−2)`) plus rare small
  cross-types (`n=5,6`), else AP-only — the census I reproduced (n=6..16).

## What was WRONG (honest correction)
> My proposed gate **"Jacobsthal run `J(k) ≥ n−k+2`" is REFUTED.** The n=8,14,20 doublings all have `J(n−2)=3`
> while `n−k+2 = 4` — so `J < n−k+2`, yet they are tight. The required run is **not** `n−k+2`. The true gate is
> `J(k) ≥ 3`-ish (a run of 3 consecutive nonunits mod `k`, `⇔ 6|k` for `k=n−2`) — the repo's HYP-2893
> "`[interval]` is a run of `v`-nonunits". My distance-dependent formula overcounted; the real condition is a
> near-constant run length gated by `k`'s small-prime structure.

## Where it leaves the lowness
- **The lowness reduces to: bound the number of Jacobsthal-swappable `k`.** Steps 1–2 are rigorous; step 4
  (primes forced) is clean; the residual is **counting the composite `k` with a sufficient nonunit run** — a
  number-theoretic problem = the repo's HYP-2893 (Jacobsthal) + HYP-3740 (lowness). I confirm and sharpen it
  (primes out, doubling = `n≡2 mod 6`) but do NOT close it.
- **The honest gap:** (a) making "large primes forced" fully rigorous (the residual-coverage failure for
  prime `q`, not just census); (b) the exact required-run length (repo's Jacobsthal interval); (c)
  multi-swap and cross-type coverage (the `k=2→9` type escapes the doubling picture); (d) bounding the count
  to `O(1)` or `O(polylog)`.

## Status
- **Rigorous (opus):** tight ⟹ covering (mult of every `q`); for `q∈(n/2,n−1]` only `q,2q` cover.
- **Verified/new (opus):** doubling family `n≡2 mod 6`, extended to **n=20** (`18→36` tight); large primes not
  swapped (census); single-swap locus sparse (n=6..16).
- **Corrected:** the gate is NOT `J≥n−k+2` (doublings have `J=3<4`); it is `J(k)≥3`-ish (repo HYP-2893).
- **Open (the crux, not closed):** bound the swappable-`k` count ⇒ full lowness (HYP-3740); rigorous
  "primes forced"; cross-type & multi-swap coverage.

Related: the-tight-skeleton-is-COMBINATORIAL-COVERING (the covering frame), PATCH-TUNING-… (the hole/patch),
covering-rigidity-…-dead-end (why not moments), THM-523 (q-witness), HYP-2893 (Jacobsthal accelerations),
HYP-3740 (lowness lemma, the crux), HYP-3750 (near-AP classification); OPEN-Q-108.
