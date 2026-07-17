# THM-961 — The locked-chain joint count (death-star-2026-07-17-S47)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCLockedChainCount.lean`,
standard trio ×5). Source: HYP-7233. Consumes THM-960's rung lock. (THM-959 =
kind-pasteur's T5 plane pass; opus's block-structure reduction still needs its
renumber — 962+ suggested.)

## Statement

1. `locked_pair_fail_iff`: for exact ratio `1 ≤ M ≤ 13` (witness form), both
   runners fail the safe band at `p/q` **iff** the bottom runner fails the
   `14M`-narrow band (`∃w, 14M·|u·p − w·q| < q`). (⇐) is monotonicity; (⇒) is
   the rung lock.
2. `locked_chain_fail_iff3`: three members `u, M₂u, M₃u` (`M₂ ≤ M₃ ≤ 13`): all
   fail iff the bottom fails the `14M₃`-band — only the bottom-top lock is
   needed; middle members ride along.
3. `card_mod_filter_eq`: **the mod transport** — at `gcd(u,q) = 1`, filtering
   `(0,q)` by any predicate of the residue `(u·p) mod q` counts the same as
   filtering residues directly (predicate-agnostic factoring of the THM-942A
   unit bijection; reusable).
4. `narrowFailN_count`: `#{r ∈ (0,q) : 14Mr < q ∨ 14M(q−r) < q} = 2⌊(q−1)/(14M)⌋`.
5. `locked_pair_count`: composed — **at coprime moduli the joint-failure count
   of an exact-ratio pair is exactly `2⌊(q−1)/(14M)⌋`**.

## The uniform deviation law (recon: `lockedchain_recon_deathstar_S47.out`)

`D_pair(M) = 2⌊(q−1)/(14M)⌋ − (q−1)/49`: **positive for M < 7, ≈0 at M = 7
(the equilibrium ratio), negative for 8 ≤ M ≤ 13** (at q = 9800: +498 at M=2
down to −94 at M=13). The pair rung of the B5 deviation ledger is UNIFORM on
locked strata — no per-stratum tables. This also explains S46's deep-commonness
on compressed blocks quantitatively: M < 7 ratios inflate joint failures by
`(q−1)(1/(14M)·2 − 1/49)`.

Recon: pair iff 300k PASS, chain iff 150k PASS, iff breaks at M=14 (q=15
family), exact count 2435 coprime cases PASS.

## Consumes / feeds

Consumes THM-960 (lock), THM-942A pattern (bijection). Feeds the liveCount
floor: on locked strata the B5 expansion's pair terms are now closed-form; the
triple/quadruple/quintuple locked terms follow the same narrow-band collapse
(chain iff) — the remaining ledger work is the NON-chain subsets (ratios > 13
between members), which are the sparse-resonance regime.
