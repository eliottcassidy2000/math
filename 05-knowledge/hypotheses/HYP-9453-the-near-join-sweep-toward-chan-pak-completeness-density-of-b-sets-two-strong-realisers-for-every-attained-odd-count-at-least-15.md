---
id: HYP-9453
title: "Toward the completeness of the Hamiltonian-path spectrum (Chan-Pak): (1) the near-join family T' = T -> M_k -> z with the arc z -> t0 reversed, H(T') = h(T) + 2^k B(T, t0), realises every odd number except 1, 7, 11, 13, 21, 47, 67 once the base tournaments are large enough (with bases of order <= 8 it realises every odd number up to 22,227 except those seven); (2) for every attained odd h the set {B(T, t) : H(T) = h} is dense enough that h + 2B sweeps (3h, ...), which would make completeness an induction on windows; (3) every attained odd h >= 15 below the top band of its order is the count of at least two strong tournaments (the strong singletons are 3, 5, 11, 13 only)"
status: >
  OPEN (klein-2026-11-02). (1) and (3) are FINITE-EXACT through base order 8 and order 9 respectively; (2) is the
  missing lemma named in the reflection: B(T, t0) counts ordered partitions of T - t0 into a path and a path from t0,
  B >= h(T), and no structural control of the B-sets is known; the pair (h, B) is not a closed state under iteration at
  t0 (counterexample pairs with equal (h, B, k) and B(T', t0) = 45 vs 47), while B(T', z) = H(T') always. The
  completeness statement itself is Pi^0_1 (bounded witnesses of order <= floor(3 log_3 h); rotation_owners_spectrum_20261020e)
  and reduces to the odd primes != 7 together with the 7p lane (THM-4094 Theorem 6.1); the repository's digital bank and
  strong-ear banks already realise every odd number <= 99,999,999 except 7, 21 (square_sum_spectrum_overnight_20261031d;
  independently audited to 999,999 in session124_million_path_spectrum_20261031), so the near-join family is a third,
  weaker construction whose interest is its one-parameter formula, not its reach.
session: klein-2026-11-02 (codex/session-ramsey-20261030)
source: 05-knowledge/results/spectrum_analogies_nearjoin_20261102/README.md
related:
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-4861-the-five-fermat-primes-are-hamiltonian-path-counts-of-tournaments-with-least-orders-3-4-6-8-11-a-vertex-transitive-tournament-has-a-prime-count-only-for-c3-and-fermat-primes-are-free-vertices-of-the-prime-square-sum-graph.md
  - 07-reflections/partners-windows-and-bands-analogies-for-the-completeness-of-the-hamiltonian-path-spectrum-klein-20261102.md
---

# HYP-9453 — the near-join sweep toward Chan–Pak completeness

## 1. The near-join family misses exactly seven odd numbers below its reach (FINITE-EXACT)

`H(T′) = h(T) + 2^k B(T, t₀)` for `T′ = T → M_k → z` with `z → t₀` reversed (PROVED in the reflection, §2). With bases of order `≤ 7` the family realises every odd number `≤ 1,437` except `1, 7, 11, 13, 21, 47, 67`; with bases of order `≤ 8` every odd number `≤ 22,227` except the same seven. **Claim:** for every `X` there is an order `n(X)` such that bases of order `≤ n(X)` realise every odd number `≤ X` other than these seven. Since 11, 13, 47, 67 are attained by strong tournaments and 7, 21 are impossible, this claim implies Chan–Pak completeness; the digital bank already proves the analogous statement through `10^8`, so the claim's value is structural (one formula) rather than numerical.

## 2. The density of `B`-sets (OPEN; the lemma that would turn §1 into an induction)

`B(T, t₀) = Σ_S h(T[S]) · h_{t₀}(T[S^c ∪ t₀]) ≥ h(T)`. If for every attained odd `h` the set `𝔅(h) = {B(T, t) : H(T) = h}` contained every integer in `[h, c·h]` for some `c > 3/2`, then `h + 2B` with `h ≤ v/3` would realise every odd `v` in `(3h, (1 + 2c)h]`, and an induction over dyadic windows would prove completeness from a finite base (the same shape as HYP-9029's tiling law for strong intervals). Computed `𝔅(h)` for `h ≤ 31` (bases `≤ 7`) have 14–42 members, contain `h`, and are not intervals; what is needed is a construction that moves `B` by one while fixing `h`.

## 3. Two strong realisers (OPEN beyond order 9)

Every attained odd `h ∈ [15, 2881]` is the count of at least two non-isomorphic strong tournaments of order `≤ 9`; the strong singletons through order 9 are exactly `3, 5, 11, 13` and band values above 2,881. **Claim:** for every attained odd `h ≥ 15` the strong fibre has at least two members. Its role: the rigid values 11 and 13 are the only ones the near-join family can never reach below 47, and rigidity is what a "leaf" is in the square-sum graph of THM-4860.

## kind-pasteur-2026-11-09-S20 note (THM-5109, HYP-9656)

- **B a power of two.** When B = 2^j the near-join is a pure cone: h + 2^K for every K > j. Such "dyadic ports" exist
  for 69 of the 148 attained odd values <= 300 at order <= 8.
- **Cone chains.** A cone step lands on a dyadic port again only if h + 3^a = 2^t (THM-5109 (D)). By Ellison's bound
  this happens only at the Pillai near-misses 1, 5, 13, and no two dyadic cone steps chain (THM-5109 (E)).
- **Bounded odd parts.** Closures that keep the odd part of B bounded are finite (HYP-9656). So the response intervals
  sought here must come from B with unbounded odd part.
