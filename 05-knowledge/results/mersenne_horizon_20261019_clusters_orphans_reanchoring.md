# The uncovered Mersenne branch at the 50000-exponent horizon: clusters, orphans, the mean-cluster law, and anchor re-anchorings

Session mac-mini-2026-10-19-mersenne (2026-10-19). Owner prompt (2026-10-14, continued): keep working on the uncovered
Mersenne branch; blocks that succeed together where single steps fail (the two-prime repair of 561); source-specific
guards and grounded children toward universal coverage; the Wilson triplet {5, 13, 563}, ternary trees and `(n−1)! mod n`;
and, today, pull in the five days of new work and extend its angles. **Status: THM-4712 PROVED (statements 2–3; 4 at
THM-4581's sketch level) + FINITE-EXACT (exact horizon to K = 50000); THM-4713 PROVED (corollary of THM-4640) + FINITE-EXACT;
HYP-9332, HYP-9333 OPEN (numerical). Universal Collatz, universal Mersenne coverage and ROOT for any new exponent remain OPEN.
Pre-audit.**

## 0. Inheritance: what the five days changed, and what this session adds

* **Same prompt, already answered in part.** THM-4640 (repair algebra at −1: collisions = increasing {2,3}-representations,
  primitive collisions = ladders, Gersonides carries, switch coverage → 1/2) and THM-4641 (Wilson/{2,4,8}/Fermat dictionary)
  by mac-mini-2026-10-14-wilson; codex-tiling's Mersenne packets of 2026-10-14c (three short paired repairs with exponent
  phases 79 mod 128, 199 mod 256, 323 mod 512; the fixed-head classifier `2·3^(E−D−1) C = R`; the finite-anchor closure
  bound: zero-defect receipts ground no exponent above the largest supplied terminal rank); codex's successor packets
  (561's unique two-prime child 561·401·3041, the obstruction repeating at 15201 = 27·563). The planned repair algebra of
  this session was therefore dropped; THM-4713 keeps only the one corollary THM-4640 did not state.
* **Also new since 2026-10-14 and used here:** THM-4651 (Giuga's sum at p² is the Wilson quotient; Korselt/Giuga parity
  divide), THM-4654 (safe Mahler codes: no two consecutive odd steps), the 2026-10-19 bank/reset packets (expanding bank
  followed by reset bank has a uniform descent head once the multiplier contracts; the Mersenne run is the extreme
  expanding word, and codex's joint-guard §2 shows its native two-run is too short to pay).
* **What was not done by anyone:** the exact horizon of the Mersenne line stood at K = 12800 (THM-4556 (vi), HYP-9242,
  barriers note §4.2, "next step 3: exact orbits up to K ~ 10^5 would test HYP-9243's exponent"). This session pushes it to
  K = 50000 and proves the mean-cluster law it was testing.
* **Board.** Anchor: the Mersenne line's cluster structure (THM-4712). Niche: anchor collisions as trunk re-anchorings
  (THM-4713). Wildcard: the Terras-transition reading of Korselt and the Mahler-orbit reading of the run (§4, dictionary
  rows, KNOWN/CITED).

## 1. The exact horizon to K = 50000 (THM-4712 (6))

`horizon.c` (GMP, 6 threads, ≈ 15 min) computes for every `2 ≤ K ≤ 50000` the odd count, the Terras count, the last odd
value before 1, the peak bit length and the payment time of `2^K − 1`. The first three columns reproduce
`runcompress_20261007/mersenne_sigma_12800.txt` exactly (`diff` empty): an independent re-derivation of the old table.

| window of odd K | orphans / odd K | fraction·√K | mean cluster size (min in window) | 2/mean·√K | members below K / √σ_T(x_K) |
|---|---|---|---|---|---|
| 3–400 | 35/199 | 2.50 | 12.4 | 2.29 | 0.44 |
| 401–800 | 13/200 | 1.59 | 32.0 | 1.53 | 0.26 |
| 801–1600 | 13/400 | 1.13 | 65.7 | 1.06 | 0.44 |
| 1601–3200 | 17/800 | 1.04 | 90.4 | 1.09 | 0.47 |
| 3201–6400 | 27/1600 | 1.17 | 118.7 | 1.17 | 0.41 |
| 6401–12800 | 30/3200 | 0.92 | 220.2 | 0.89 | 0.51 |
| 12801–19200 | 25/3200 | 0.99 | 247.2 | 1.02 | 0.47 |
| 19201–25600 | 18/3200 | 0.84 | 368.2 | 0.81 | 0.56 |
| 25601–38400 | 31/6400 | 0.87 | 416.5 | 0.86 | 0.53 |
| 38401–50000 | 21/5800 | 0.76 | 534.7 | 0.79 | 0.62 |

Facts (all FINITE-EXACT, `analyze_horizon.out`, `analyze_horizon_extra.out`):

* 231 clusters (level sets of `(σ_T − K, odd)`); the odd count alone separates them at every cutoff (22, 36, 57, 103,
  136, 177, 231 distinct odd counts for K ≤ 100, 400, 1200, 6000, 12800, 25000, 50000).
* Reset pairs (THM-4556 (ii)) hold without exception; the least shift to a smaller mate is odd in all 24,769 cases
  (THM-4556 (iii)); it is 1 in 95.52% of them, and the counts for odd `D = 3, …, 41` fit `D^(−2.01)`; the largest is 243.
* 23 clusters are intervals; the hole fraction is 0.147; at most three cluster spans cover any exponent; the largest
  cluster has 1276 members (`40713 ≤ K ≤ 42002`, 14 holes), `6.3 √K_min`.
* Terminal channels: 5 (46,186 exponents, 215 clusters), 85 (1,912 / 8), 341 (1,900 / 8); nothing else; no mixed cluster.
* Payment: `(pay − K)/K = 2.8213 ± 0.0714` (K ≥ 100) against the drift value
  `(log_2 3 − 1)/(1 − (log_2 3)/2) = 2.8188`; peak bits `1.5851 K` (run end `log_2 3`).
* Orphans by two-run length `J`: rates 0.0056, 0.0085, 0.0113, 0.0078 for `J = 1, 2, 3, ≥ 4` (K ∈ [1001, 50000]): the
  uncovered branch is not concentrated in the deep reset-2 shells.
* The 95 orphans in (12800, 50000] are listed in `analyze_horizon_extra.out` (first: 12963, 13079, 13129; last: 48445, 49601).

## 2. The mean-cluster law (THM-4712 (2)–(5))

In the Haar model of THM-4581 6(d) the children `(x+1)/3^D − 1` of a post-run source `x` are driven by one source, so the
chains are dependent; but the **mean** number absorbed by Terras time `T` is `Σ_D P(τ_D ≤ T)` regardless. THM-4581's proof
of 3(i) shows that the level `k` of a chain is a martingale with increments in `{−1, 0, 1}`; Azuma's maximal inequality
then gives `P(τ_D ≤ T) ≤ exp(−D²/(2T))` and the clean bound `E[N_T] < √(πT/2)`. THM-4581 (4)'s union bound with its
`|k_0|` first-passage term gives `E[N_T] ≥ c √(T/log T)` at the sketch level. So the mean number of deletion children merged
by time `T` is `√T` up to a root-log — the exponent of HYP-9243 for the mean, in the Haar model. The SRW skeleton predicts
the envelope `√(T/π) = 0.564 √T`; on the integers (members below K, `T = σ_T(x_K)`) the ratio runs 0.41–0.62 from K ≈ 800
(table above). The `E[J] ≈ 10` excursion deficit of THM-4581 (7) is not visible at this scale.

**Blocks that succeed together.** A cluster is a block of exponents with one common future. The single step `D = 1`
repairs 95.5% of the odd exponents that have any repair; longer blocks follow the `D^(−2)` law; orphans (0.36% in the top
window) are the residue. When the `D = 1` chain of `M_K` is absorbed, the whole cluster of `K − 1` joins at that moment —
this is the barriers note's simultaneous absorption of 1916 deletions, read at the horizon. None of it grounds anything:
certificates preserve `σ_T(x_K)` (barriers Prop. 4.1) and zero-defect receipts cannot ground exponents above the largest
supplied terminal rank (codex, 2026-10-14c). Coverage → 1 (THM-4581 6(d), THM-4712 (4)) and grounding power 0 coexist.

## 3. Anchor collisions are trunk re-anchorings (THM-4713)

THM-4640 (1) listed, as a FINITE-EXACT fact, that `2^a − 3` (the numerator of the single-letter head) has a second
canonical representation exactly for `a ∈ {7, 13, 17, 21, 28, 31}` (`a ≤ 31`). The reason: a head `h` of total `A` shares
the numerator iff `f_h(−1) = (4^k − 1)/3 − 2^(1−A)` with `a = A − 1 + 2k` — it lands one unit of its own denominator below
a predecessor of 1. For `k = 1` the numerator of `h` is a Mersenne number `2^(A−1) − 1`: the heads are

    (4,1,1) ↦ 31,   (5,1,2,3,1) ↦ 2047,   (4,3,2,2,1,1,3) ↦ 2^15 − 1,   (6,1,3,1,1,1,2,2,3) ↦ 2^19 − 1,
    (4,1,3,6,2,3,3,2,2,3,1) ↦ 2^29 − 1,

and `a = 28` is the one `k = 2` case, `(5,2,3,1,1,1,4,1,1,2,1,1,2)` landing below 5. The list is complete to `a ≤ 50`
(`anchor_reanchoring.py`, about four minutes): the only further exception is `a = 49`, with the 23-letter head
`(4,4,2,1,4,2,1,3,2,2,1,1,1,1,2,1,2,2,1,4,1,5,1)` of numerator `2^47 − 1` (`k = 1`, distance 22); an `a ≤ 64` run did not
finish under machine load. Each row is the uniform switch
`(a+1, c) ~ (h, c + 2k)` of distance `|h| − 1` (THM-4555 (iv)); the first row is the Catalan carry `9 = 8 + 1`
(THM-4640 (5)) and the block of THM-4601 (iii). Among `2^a − 3`, `3 ≤ a ≤ 400`, only `125 = 5^3` is a perfect power.

## 4. Dictionary rows (KNOWN / CITED / NUMEROLOGY; no theorem numbers)

* **Korselt as a transition (KNOWN).** For odd `n`, the Terras map `T(x) = (3x+1)/2` is an affine bijection of `Z/n` with
  fixed point `−1`, conjugate to multiplication by `3/2`; `T^(n−1) = id` on `Z/n` iff `n | 3^(n−1) − 2^(n−1)`. The composite
  such `n` are OEIS A073631 (65, 133, 529, 793, 1105, 1649, 1729, …; 339 below 2·10^6, `threehalves_pseudoprimes.out`);
  every Carmichael number coprime to 6 is one (R. Israel's comment there; 52 of the 55 below 2·10^6), and 561, 62745,
  656601 are excluded by their prime 3 — the multiplier of the transition. The smallest member is 65 = 5·13, the product of
  the two small Wilson primes (`ord_5(3/2) = 2`, `ord_13(3/2) = 4`): a coincidence, recorded as such. The non-squarefree
  members below 2·10^6 are the nine multiples of 23² (23 is a (3/2)-Wieferich prime), as A073631 conjectures.
* **The run as a Mahler orbit (elementary).** For `ξ = 2^E`, `ξ (3/2)^n` is an integer for `n ≤ E` and
  `{ξ (3/2)^(E+1)} = 1/2` exactly: the Mersenne source is a Z-number candidate for E steps that fails on the boundary at
  step E+1; afterwards `{ξ(3/2)^n} = (3^n mod 2^(n−E))/2^(n−E)`, the low bits of `3^n`. THM-4654's safe codes forbid two
  consecutive odd steps; the Mersenne run is the opposite extreme, E − 1 consecutive odd steps.
* **Wilson primes on the Mersenne line (NUMEROLOGY).** The orbit of `M_3 = 7` ends `… → 13 → 5 → 1` (the two small Wilson
  primes consecutively), and `U(563) = 845 = 5·13²`; so all three Wilson primes lie on one Collatz trajectory
  (`563 → 845 → 317 → 119 → 179 → 269 → 101 → 19 → 29 → 11 → 17 → 13 → 5 → 1`). THM-4641 already records `5 = 2^3 − 3`,
  `13 = 2^4 − 3` as the anchor reciprocals of letters 3 and 4; nothing structural is claimed.
* **Channels and THM-4517 (HYP-9333).** The three terminal channels of the Mersenne line and their shares against the
  trunk-entry basin densities `e_2, e_4, e_5` of THM-4517.

## 5. What would change the picture

* A fourth terminal channel, or an exponent covered by four cluster spans, or the orphan constant turning upward.
* A certificate type with nonzero odd-count defect on the Mersenne line (a child whose odd count differs from the source's
  by more than the deletion distance): none exists among Mersenne-to-Mersenne merges, which are all deletions, so it must be
  a non-Mersenne child — the only route to grounding (barriers note, next step 2).

## 6. Reproduction

    cc -O2 -o horizon horizon.c -I/opt/homebrew/include -L/opt/homebrew/lib -lgmp -lpthread
    ./horizon 50000 6 horizon50000 && cat horizon50000_t*.txt | sort -n > mersenne_horizon_50000.txt
    python3 analyze_horizon.py mersenne_horizon_50000.txt; python3 analyze_horizon_extra.py
    python3 anchor_reanchoring.py 44; python3 threehalves_pseudoprimes.py

All in `04-computation/experiments/mersenne_horizon_20261019/`; the data file is committed gzipped (`.txt.gz`, 443 KB).
