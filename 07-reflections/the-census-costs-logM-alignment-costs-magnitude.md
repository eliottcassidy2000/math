# The census costs log M; alignment costs magnitude

*kind-pasteur-2026-07-03-S34. The owner asked me to crack the open crux — the compressed
covering family, where there is no far runner to peel and the rational census denominator
grows. I could not prove LRC(14) (that is the crux, and the crux is the conjecture). But the
disagreement in the fleet resolved cleanly, and the mechanism underneath it is sharp enough
to state, and one piece of it is rigorous.*

## The disagreement, and who was right

opus (S55) briefly held that compressed covering families are census-able at **bounded** `q`
(17–19), so the crux dissolves. mac-mini (S27) retracted the mirror claim and asserted the
opposite: the census denominator **grows**, `q ~ log M`, and the crux is scale-independent
CRT capacity. Both had earlier been fooled by weak adversaries (MISTAKE-095/096/097/098) —
constructions that block one modulus per runner and so give a false ceiling.

The decisive computation builds the strongest *compressed* adversary directly: 13 speeds all
in `[N, 2N]` (ratio `≤ 2`, genuinely no far runner), divisibility-blocking `{2, …, Q}` for
`Q` as large as 13 runners of magnitude `N` allow. Its census denominator:

| N | 10³ | 10⁶ | 10⁹ | 10¹² | 10¹⁵ |
|---|-----|-----|-----|------|------|
| q_min | 57 | 127 | 181 | 265 | 335 |

`q_min ≈ 8.7·ln M`, growing without bound, ratio `≤ 2` throughout. **mac-mini was right.**
opus's "bounded 17–19" was an under-magnitude artifact; it is already retracted in
MISTAKE-097. There is no fixed-`q` census.

## The one rigorous statement: q* ≤ 13 ln M

The census denominator equals the **first free modulus** `q*` — the smallest `q` dividing no
runner — because any `q` that divides a runner `v_i` is doomed (`v_i·a ≡ 0 (mod q)` for every
`a`, so that runner is at an integer at every `t = a/q`). And `q*` is bounded, provably:

> Every prime `p < q*` divides some runner (else `p` would itself be a free modulus below
> `q*`). Distinct primes dividing runners multiply into the lcm of the runners, so
> `∏_{p < q*} p  |  lcm(v_1, …, v_13)  ≤  ∏ v_i  ≤  M¹³`. Taking logs,
> `θ(q*) = Σ_{p<q*} ln p ≤ 13 ln M`. By the prime-number theorem `θ(x) ~ x`, so
> **`q* ≤ (1+o(1))·13 ln M = O(log M)`**.

Verified in every row (`θ(q*)/ln M ≤ 8.7 < 13`, with room). This is clean and unconditional:
the first free modulus — and hence the census denominator, *if* the witness lives there — is
`O(log M)`. The census is never a fixed check, but it is always a *logarithmic* one.

## The mechanism: divisibility is 1×, covering is 13×

Why does the witness live at the first free modulus? Because the two ways to block a modulus
`q` cost wildly different amounts of the adversary's finite CRT capacity (`13 ln M` total, the
log-magnitude of 13 runners):

- **Divisibility-block** `q`: put one runner `≡ 0 (mod q)`. Cost `= ln q` (one prime factor).
- **Covering-block** `q` (free `q`, no runner at 0): arrange the 13 danger sets
  `{a : ‖v_i a / q‖ < 1/14}` — dilated intervals of size `~q/7` each, total `13/7 ≈ 1.86 > 1`
  — to cover `ℤ/q`. This is possible (measure permits it, and greedy search finds covers for
  every `q` tested up to 129), **but it constrains all 13 residues `v_i mod q` at once**.
  Cost `~ 13 ln q`.

Divisibility is thirteen times cheaper. So the optimal adversary spends its capacity
divisibility-blocking `{2, …, Q}` with `Q ~ log M`, and the first free modulus `q*` is left
with *un-arranged* residues — where the danger sets do **not** cover (generic residues leave
`~(6/7)¹³ ≈ 13.5%` of `ℤ/q*` safe), so a witness sits there. `q_min = q*`.

## Alignment costs magnitude — why the two-sided architecture closes

The adversary *can* pay the 13× and cover the free primes: a hybrid that covers 6 free primes
pushes `q_min` from 127 to 148. But covering `q` by residue-alignment forces each runner into
a prescribed class `mod q` — which, by CRT, requires a **cofactor divisible-adjusted through
`∏ (covered primes)`**. That product spreads the magnitudes: in the measured hybrid the ratio
`max/min` blew from `≤ 2` to `410–2777`. And this is not a bug of the construction — it is the
content:

> **To cover a free prime, the adversary must multiply the CRT modulus into the runners. That
> either (a) raises the whole family's magnitude to a new `M'`, keeping `q_min = Θ(log M')` —
> self-consistent, no super-logarithmic blow-up; or (b) makes one runner `13×` larger than the
> rest, i.e. DOMINANT — which is exactly the peelable case (`far_peel_lonely_of_cite`, opus's
> `covering_lonely_of_dominant_or_compressed`, HYP-4054).**

So the only move that defeats the log-`M` census is the same move that either preserves the
log-`M` bound or triggers the far-peel. **Alignment and compression are mutually exclusive;
alignment and peelability are the same phenomenon.** This is the double-alignment impossibility
(mac-mini's kernel) seen from the cost side: you cannot both keep the family compressed *and*
buy your way past the first free modulus — the currency of alignment is magnitude, and
magnitude is what the peel consumes.

That is why the two-sided architecture is closed as a *structure*: dominant → peel; compressed
→ census at the first free modulus `≤ 13 ln M`, where alignment is unaffordable.

## What is genuinely still open

The reduction is now crisp and, on one axis, rigorous:

> **LRC(14) ⟺ for every covering 13-family of magnitude `M`, the danger sets fail to cover
> `ℤ/q*` at the first free modulus `q* ≤ 13 ln M`** (equivalently, a lonely `a/q*` exists) —
> and by the MSS bound `M ≤ 91¹²`, this is a finite check over `q ≤ ~450`.

The rigorous half is `q* ≤ 13 ln M`. The open half is *the witness exists there*: that the
13 danger sets, with the residues the adversary is forced into, cannot cover `ℤ/q*`. The
cost argument says the adversary is out of capacity at `q*` and so cannot arrange a cover —
but "out of capacity ⟹ generic residues ⟹ no cover" is a genericity step, and turning it
into a proof is exactly LRC(14). The character-sum lower bounds that would decide it are the
known-insufficient direction for `n = 14`. So the crux does not fall here — but it is now a
single, magnitude-independent covering statement at a `O(log M)`-bounded modulus, with the
escape route (alignment) proven to cost the one resource (magnitude) that hands the family to
the peel.

---
*Linked: [[the-rational-irrational-duality]] (S31–33, the region-measure / peel side),
[[the-bounded-denominator-route]] (S28, corrected S29 — the census side, now with an explicit
`O(log M)` denominator). Scripts: `lrc14_compressed_crux_kps_S34.py` (growth),
`lrc14_free_modulus_bound_kps_S34.py` (the rigorous `q* ≤ 13 ln M`),
`lrc14_hybrid_primes_kps_S34.py` (alignment-costs-magnitude). Resolves the opus/mac-mini
disagreement (mac-mini right; MISTAKE-097). HYP-4040 (mac-mini CRT capacity), HYP-4054 (opus
dominant/compressed dichotomy), HYP-4053 (kps far-peel threshold).*
