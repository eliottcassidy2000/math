# The wide residual is the log-census, and the deep well was never the gap

*mac-mini-2026-07-03-S26. Working the honest remaining crux of LRC(14) — the compressed / wide-far case.*

## Two corrections that reshaped the picture

**The deep well is census-able.** The canonical "hard" family `{1,…,12,182}` (`182 = 13·14`) is lonely at
`t = 3/40`: `min-dist = 3/40 = 0.0750 > 1/14` (runner `1` is the tightest; there is no runner `13`). It does
NOT need the Eisenstein rational `14/183`. `M = 14/183` is the *maximum* min-distance, not the only lonely
point — loneliness (`≥ 1/14`) holds on a positive-measure set that includes `3/40`. The "lonely only at
`14/183`" reading came from measuring the `p = 1` covering sieve (`t = 1/q`), which does fail; the general
`p/q` census (kps `lonely14_of_ratio`, here `p = 3`) reaches it. So the beautiful Eisenstein resonance
(`14` a primitive 6th root mod `183`) is real mathematics but **not** the census gap.

**The gap is the wide-compressed lcm-blockers.** Hunting for the family that actually defeats a fixed census,
the answer is not a tight AP or a resonance — it is the **lcm-product band-blocker**: `13` runners all `~N`
(compressed: no ratio-`13` dominant), each divisible by the `lcm` of a *group* of band moduli, so that every
small `q` divides some runner (residue `0` = danger). A weak adversary (one modulus per runner) looks
bounded (`q ≤ 56`); the real adversary packs `lcm(15..21) ≈ 2·10⁸ ≤ N` into one runner and blocks seven
moduli at once. Measured: witness `q = 61, 71, 83, 101, 109, 121, 128` at `N = 10⁴,…,10¹²` — slope `≈ 9` per
decade, i.e. **`q ~ 3.6·ln(max-speed)`, unbounded.** This is HYP-4040 (no uniform band) reaching inside the
*compressed* class, not just the dominant lcm family.

## The crux, stated exactly

> Every compressed covering `gcd=1` family with `max|v| = M` is lonely at some `a/q` with
> `q ≤ Q(M) = O(log M · log log M)` — and no smaller uniform bound exists.

Combined with the far-peel (the dominant-far case, step 1 = THM-609 now formalized) + the citation +
non-covering, this **log-census** closes LRC(14). It is the honest content of the whole wide-residual /
two-sided-split program.

## Why it is hard (the witness, not the modulus)

The tempting half is easy: a **free modulus** exists below `O(log M · log log M)`, because `13` speeds `≤ M`
carry at most `13·(log M / log 17)` primes `> 14`, while `{15,…,Q}` has `~Q/ln Q` primes — so for
`Q = O(log M · log log M)` some prime divides no speed (elementary prime counting).

But a free modulus does **not** hand you the witness. At a free prime `q`, the danger set has `~q/7` residues,
so the "bad" numerators `{a : v·a ∈ danger}` number `~13·q/7 > q` across the `13` speeds — the crude union
bound **fails to leave a good `a`.** The good `a` exists (the census finds it) only because those bad sets
*overlap* in a structured way. So the compressed closure is genuinely the census `∃ a ∈ (0,q), ∀ v, va ∉ danger`
at a `q ~ log M`, and its positivity is a real arithmetic fact, not a counting corollary. That is exactly where
the difficulty lives, and it is why every route — bounded-denominator census (kps), renormalization tower
(opus HYP-3901), measure floor (klein) — is attacking the same wall from a different side.

## What this rules in and out

- **Rules out**: a fixed-`q` census closing everything; the Eisenstein `183` being the obstruction; the deep
  well being special.
- **Rules in**: the closure is a *magnitude-indexed* census `q ≤ Q(max-speed)`, `Q = Θ(log)`; the extremal
  families are the lcm-product blockers; the open kernel is the *witness existence* at the log-scale free
  modulus (equivalently, the overlap structure of the `±k` bad-numerator sets).

The mathematics keeps pointing at the same place: controlling all residues at a logarithmically-growing
denominator, where the naive union bound is exactly critical (`13·(2/14) = 13/7 > 1`). That `13/7 > 1` is the
same over-one margin that makes LRC hard in the first place — here it reappears as the reason the free modulus
is not enough. The wall has one face.
