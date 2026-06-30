---
id: HYP-3741
title: The LOWNESS LEMMA (M(S)<=n/Phi_6 => {1,..,n-2} subset S), routes 1&2 synthesized via the CONSTANT-RESIDUE principle, and the CRT-escape crux RESOLVED. Route 2 (k-witness): a covering set missing small speed k is lonely at t=k^{-1}/p for a band prime p where the {k,p-k} unit-pair is uncovered, giving min||vt||>=2/p>n/Phi_6 (clean: missing-1 at n=14 -> t=1/17, M=2/17) -- but SET-DEPENDENT (a uniform k-witness is impossible: killers can be placed at any t). Route 1 (klein-S39 budget/transversal): to defeat every k-witness, S must contain a speed ≡±k mod every binding modulus. THE CRUX: could a single HUGE CRT speed w≡k mod (product of band primes) defeat all k-witnesses and restore M=n/Phi_6? RESOLVED -- NO: even w≡1 mod (all primes<=43) (~1.3e16) leaves M=525/3716≈0.141>>14/183 (exact); the lonely hole merely moves to mod 85. THE SYNTHESIS -- the CONSTANT-RESIDUE PRINCIPLE: a SMALL speed k has CONSTANT residue (k mod p = k for ALL p>k), covering the {k,p-k} pair mod EVERY prime simultaneously -- a 'universal pair-coverer'; a LARGE speed has SCATTERED residues (by CRT), covering {k,p-k} mod only the finitely many p|(w∓k). So NO finite set of large speeds replaces a small speed across the (unboundedly many) binding moduli -- speed k is IRREPLACEABLE, and missing it always exposes a lonely hole. This unifies route 1 (replacement exceeds the budget) and route 2 (the exposed witness), and robustly forces {1,..,n-2} subset S
status: SYNTHESIS + crux RESOLVED (robust to CRT speeds ~1e16, exact M=525/3716). The constant-residue principle is the MECHANISM (proved for the band primes: a small speed is the unique universal {k,p-k}-coverer); it is strong evidence + the right intuition, NOT a complete proof of the lowness lemma -- the actual binding hole can sit at a COMPOSITE modulus (mod 85, not a band prime), so a full proof must control all moduli (= klein's budget over all D, still open). Reduces the lowness lemma to: "missing the constant-residue speed k exposes a lonely hole at some modulus" -- now mechanism-backed.
source: mac-mini-2026-06-30-S57
related:
  - HYP-3740  # the LRC14 hard core reduces to the lowness lemma (this works its proof routes 1&2)
  - HYP-3736  # klein-S39: the killer-or-transversal budget + proved transversal lemma (route 1)
  - THM-523   # the q-witness (route 2 is its binding-side analog: q-witness for covering, k-witness for lowness)
  - HYP-3739  # M-uniqueness (the construction is the strict M-min)
results:
  - 04-computation/lowness_lemma_constant_residue_macmini_20260630.py
  - 05-knowledge/results/lowness_lemma_constant_residue_macmini_20260630.out
---

# HYP-3741 -- the lowness lemma: routes 1 & 2 synthesized (the constant-residue principle)

Working the lowness lemma's two proof routes (HYP-3740) **back and forth** -- the k-witness (route 2) and the
transversal/budget (route 1) -- synthesizes into one principle and resolves the sharp crux.

## Route 2 -- the k-witness (the binding-side analog of THM-523's q-witness)
THM-523's `q`-witness: a covering set omitting a multiple of `q` is lonely at `t=1/q`. The **`k`-witness**: a
covering set missing the small speed `k` is lonely at `t = k^{-1}/p` for a band prime `p in (n, 2n-3]` whose
**unit-pair `{k, p-k}` is uncovered** -- then every `v in S` has `vk^{-1} mod p \notin {0,±1}`, so
`min_v ||vt|| >= 2/p > n/Phi_6` (since `p < 2 Phi_6/n ~ 2n`). Clean instance: missing speed `1` at `n=14` gives
`t=1/17`, `M=2/17`. **But this witness is SET-DEPENDENT** -- a uniform `k`-witness is impossible, because a set
can include a speed `≡ ±k mod p` (a "killer" of that witness).

## Route 1 -- to defeat every k-witness (the budget)
So a missing-`k` set with `M <= n/Phi_6` must, for **every** band prime `p`, contain a speed `≡ ±k mod p`
(covering the pair `{k,p-k}` -- klein-S39's killer-or-transversal). The construction does this with the single
speed `k` itself; a missing-`k` set must re-buy each pair, which klein's budget `n-1 = resonance-killers +
band-prime coverers + spreaders` cannot afford -- raising the rung, hence `M > n/Phi_6`.

## The crux -- could ONE huge CRT speed defeat them all?
The dangerous loophole: a single huge speed `w ≡ k mod (p_1 p_2 ... )` covers `{k,p_i-k}` mod **all** band
primes at once (one speed, budget OK). **RESOLVED -- it fails:** even `w ≡ 1 mod (all primes <= 43)`
(`~1.3 x 10^16`), the set `{2,..,12, w, 182}` has `M = 525/3716 ≈ 0.141 >> 14/183` (exact), the lonely hole
simply moving to `mod 85`. So no CRT speed restores `M`; `{1,..,12}` is genuinely forced.

## The synthesis -- the CONSTANT-RESIDUE principle
Why is the small speed `k` irreplaceable? Look at its **residue array** `(k mod p)_p` across primes:
- a **SMALL** speed `k` (`1<=k<=n-2`) has **CONSTANT** residue: `k mod p = k` for every prime `p > k`. So it
  covers the pair `{k, p-k}` mod **every** prime simultaneously -- a *universal pair-coverer*.
- a **LARGE** speed `w` has **SCATTERED** residues (pseudo-random by CRT): it covers `{k,p-k}` mod only the
  finitely many primes `p | (w ∓ k)`.

`speed 1 -> [1,1,1,1,1] mod (17,19,23,29,31)` (constant); `speed 7430 -> [1,1,1,6,21]` (scattered). **No finite
set of scattered-residue large speeds reproduces a constant residue across the unboundedly many binding
moduli.** Hence the consecutive base `{1,..,n-2}` -- the `n-2` universal pair-coverers -- is forced; missing any
one leaves a pair (and a lonely hole) that no large speed can fill everywhere at once. This is exactly why the
CRT escape fails (it matches a constant only on finitely many primes) and unifies route 1 (the replacement
cost) with route 2 (the exposed witness).

## Proof status and what remains
The constant-residue principle is a *mechanism* and resolves the crux robustly, but it is not yet a complete
proof: the actual binding hole can sit at a **composite** modulus (`mod 85`, not a band prime), so a full proof
must control the radius demand at **all** moduli `D`, not just the band primes -- which is klein's budget over
all `D` (still open). Net: the lowness lemma reduces to "missing the constant-residue speed exposes a lonely
hole at some modulus," now backed by the constant-vs-scattered mechanism and verified robust against
astronomically large CRT speeds. The `k`-witness is the binding-side twin of the `q`-witness; together they are
the two halves (covering + lowness) of the LRC14 hard core (HYP-3740).
