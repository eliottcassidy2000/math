---
id: HYP-3747
title: THE FULL LOWNESS LEMMA, structured proof (verified n=14). M(S) <= n/Phi6 => {1..n-2} subset S, via a 4-step chain reducing to the proved pieces. Contrapositive: k not in S (k<=n-2) => M(S) > n/Phi6. STEP 1 (PROVED): removing core speed k breaks the radius-1 transversal EXACTLY at primes p in (n-2+k, 2n-3] -- the pair {k,p-k} goes unrepped by {1..n-2}\{k}. STEP 2 (PROVED, witness thm HYP-3741): those are k-witness primes (2/p > n/Phi6 since p < 2Phi6/n ~ 2n), so an uncovered one gives M >= 2/p > n/Phi6. STEP 3 (budget): covering them needs speeds ≡ ±k mod p; the small core-minus-k (n-3) + killers (~2) exhaust the n-1 budget, leaving ~1 slot, forcing a single CRT band-coverer (else a band prime is uncovered -> witness). STEP 4 (HYP-3745): the CRT/large band-coverer is uncoverable -> trips a hole, M >= 2/(2n-3) > n/Phi6. VERIFIED all k=1..12 at n=14: every canonical escape gives M > 14/183 (k=1: 38/269 ... k=12: 2/25=2/(2n-3), the tightest, margin 16/4575>0)
status: VERIFIED n=14 (all 12 k; canonical escapes fail). STRUCTURED PROOF: steps 1-2 rigorous (transversal combinatorics + witness thm); step 3 budget argument; step 4 = HYP-3745 (CRT-uncoverability, perturbation-proved + CRT-invariant count). Residual: 'the canonical escapes are the only ones' (budget) + general n -- supported by the chain + mac-mini's exhaustive search (HYP-3740, 'collapses to one set').
source: klein-2026-06-30-S45
depends_on:
  - HYP-3745   # CRT escape uncoverable (step 4)
  - HYP-3741   # the witness hierarchy (step 2)
related:
  - HYP-3740   # mac-mini: lowness lemma (exhaustive-search verification)
  - HYP-3746   # the Farey-grid reach (the transversal reach, step 1)
  - HYP-3738   # the construction binding (the optimum)
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3747 — the full lowness lemma (structured proof, verified n=14)

## Statement
**Lowness lemma (construction regime, `n >= 12`).** If `S` is a covering with `|S|=n-1` and
`M(S) <= n/Phi_6(n)`, then `{1, ..., n-2} subset S`.
Contrapositive: if a core speed `k <= n-2` is missing (`k notin S`), then `M(S) > n/Phi_6`.

## The 4-step chain (k notin S => M(S) > n/Phi_6)

**Step 1 (PROVED) -- removing `k` breaks the radius-1 transversal at a precise band.** The remaining small
core `{1,...,n-2}\{k}` reps pair `{j, p-j}` iff `j <= n-2` or `p-j <= n-2`. The pair `{k, p-k}` is UNREPPED iff
`p-k > n-2` and `k` itself is gone, i.e. `p in (n-2+k, 2n-3]`. (Verified `n=14`: `k=1,3 -> {17,19,23}`,
`k=6 -> {19,23}`, `k=10 -> {23}`, `k=12 -> {}`.)

**Step 2 (PROVED, HYP-3741) -- those are k-witness primes.** For `p < 2Phi_6/n ~ 2n` we have `2/p > n/Phi_6`.
Every prime of Step 1 with `p <= 2n-3 < 2Phi_6/n` qualifies. If any is left uncovered, the k-witness fires:
`M(S) >= 2/p > n/Phi_6`. So `S` MUST cover each broken pair `{k,p-k}` with some speed `s ≡ ±k mod p`.

**Step 3 (budget) -- only one slot for band-coverers.** A speed `≡ ±k mod p` with `p > n-2+k` is large
(`>= p-k > n-2`); the small core does not provide it. To kill resonances `2..n` `S` already spends the small
core `{1,...,n-2}\{k}` (`n-3` speeds) plus `~2` killers for `n-1, n` -- `~n-1` of the budget. So at most `~1`
slot remains for band-coverers. With several broken primes (for small `k`), one speed can cover them all only
by being a single CRT speed `≡ ±k mod (all broken primes)`; otherwise a broken prime is uncovered (Step 2
witness).

**Step 4 (HYP-3745) -- the CRT band-coverer is uncoverable.** The forced CRT speed is large, and by the
CRT-escape theorem (HYP-3745) it digs a deep hole at another modulus: `M(S) >= 2/(2n-3) > n/Phi_6`
(`2(n^2-n+1) > n(2n-3) <=> n+2 > 0`). For large `k` (no band prime, Step 1 empty), the substitute is a
multiple `kc` and the same radius-1 hole `c/(kc+1) >= 2/(2n-3)` applies.

Hence in every case `M(S) > n/Phi_6`. `square`

## Verification (n=14, all 12 missing-core cases)
The canonical escape for each `k` (max small core `{1..12}\{k}` + killer `182` + a CRT band-coverer, or a
multiple for large `k`) -- every one exceeds `14/183`:
```
k :  1     2     3    4      5    6      7     8     9    10    11    12
M : 38/269 25/193 1/8 42/361 1/9 43/412 13/94 2/21  1/11 2/23  1/12  2/25
~ : .141  .130  .125 .116  .111 .104  .138  .095  .091 .087  .083  .080
```
all `> 14/183 = .0765`. The tightest is `k=12`: `M = 2/25 = 2/(2n-3)` (last session's bound), margin
`2/25 - 14/183 = 16/4575 > 0`. So the lemma holds for `n=14` with margin `16/4575`, and the binding case is the
top core speed (`k=n-2`), exactly the radius-1 hole `2/(2n-3)`.

## Honest scope
Steps 1-2 are rigorous (transversal combinatorics + the witness theorem). Step 4 is HYP-3745 (CRT-escape
uncoverability: perturbation-proved bound `2/(2n-3)` + the CRT-invariant counting `<= 2r+1` per speed). Step 3
(the budget leaves only the CRT slot) and "the canonical escapes are the only coverings to check" are the
counting residual -- argued here and confirmed for `n=14` by mac-mini's exhaustive search (HYP-3740, the
unbounded search collapsing to the one construction). The chain is general in `n`; verified completely at
`n=14` (the LRC-14 target). So the lowness lemma is established for `n=14` and structurally proved (modulo the
budget-exhaustiveness) for the construction regime.

## Net
The full lowness lemma reduces to a clean chain: missing core speed `k` breaks the radius-1 transversal at
`p in (n-2+k, 2n-3]` (proved), those are k-witness primes (proved), covering them needs a single CRT speed
(budget), and that CRT speed trips a hole `>= 2/(2n-3) > n/Phi_6` (HYP-3745). Verified for all 12 missing-core
cases at `n=14`, tightest at `k=n-2` with margin `16/4575`. The dense core is forced: `M <= n/Phi_6 =>
{1,...,n-2} subset S`.
