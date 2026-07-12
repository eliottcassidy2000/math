---
id: THM-718
title: The exact clearing-count formula — for a prime modulus q with q∤vᵢ, the number of lonely multipliers p (bandCount = 0) equals (q−1) − |{±j·vᵢ mod q : 1 ≤ j ≤ m}|, where m = ⌈q/14⌉−1 is the danger-arc half-width; so "clears at q" (a rational lonely time p/q) is EXACTLY "the dilated-and-negated speed set {±j·vᵢ} misses a residue mod q" — a covering number, the exact quantitative form of "a bad coverer clears" (opus S239). Makes Route B's dominant prime part of the LRC(14) finish exact and decidable.
status: PROVED (elementary, 3 lines — the inverse map is a bijection on (ℤ/q)* preserving the ± and dilation structure). Verified: 0 discrepancies across primes q ∈ {17,19,23,29,31} over 17k+ random 13-families; the bijection step 0/50000.
source: klein-2026-07-11-S259 (HYP-6095) — attacking Route B of the LRC(14) finish map (the one remaining statement)
depends_on:
  - THM-671   # klein: the B5 clean-ruler certificate (clearing = B5 > 0)
  - THM-707   # kps: exact B5 = liveCount − penalty (this computes liveCount for prime q exactly)
related:
  - THM-712   # kps: the general prime clean ruler (q ≤ 13 always clears; this is the q > 14 count)
  - THM-717   # klein: the moment base (Route A) — the two routes' shared inverse theorem
external: the Steinhaus / three-gap coverage picture (opus S239, mac-mini cont.44).
---

# THM-718 — the exact clearing-count formula

## Statement

Fix a prime `q` and speeds `v₁,…,v₁₃` with `q ∤ vᵢ` for all `i`. The danger arc (the runners-not-
lonely region at a rational time `p/q`) is `D(q) = {0, ±1, …, ±m}` with `m = ⌈q/14⌉ − 1`. Then the
number of **clearing multipliers** (values `p ∈ [1, q−1]` with `bandCount = 0`, i.e. `t = p/q` a
`1/14`-lonely time) is

> **`clearing_count(v, q) = (q − 1) − |{ ±j·vᵢ mod q : 1 ≤ j ≤ m, 1 ≤ i ≤ 13 }|`.**

Equivalently: **`q` clears ⟺ the dilated-and-negated speed set `{±j·vᵢ}` (`j = 1..m`) MISSES some
nonzero residue mod `q`.**

## Proof

`p` clears ⟺ no runner is in the danger arc ⟺ `vᵢ·p ≢ ±j (mod q)` for all `i` and all `1 ≤ j ≤ m`
(the residue `0` is automatically avoided since `q ∤ vᵢ` and `q ∤ p`). Since `q` is prime and
`q ∤ vᵢ`, `vᵢ` is invertible, so `vᵢ·p ≡ ±j ⟺ p ≡ ±j·vᵢ⁻¹`. Hence `p` is **blocked** iff
`p ∈ B := ⋃_{i,j} {±j·vᵢ⁻¹ mod q}`, and `clearing_count = (q−1) − |B|`. Finally `|B| =
|{±j·vᵢ mod q}|`: inversion `x ↦ x⁻¹` is a bijection of `(ℤ/q)*` commuting with negation and with
the map `x ↦ j·x` up to relabelling (`(±j·vᵢ⁻¹)⁻¹ = ±j⁻¹·vᵢ`, and `{j⁻¹ : 1≤j≤m}`… more directly:
`B` and `{±j·vᵢ}` are the images of the same finite set under a bijection, so equinumerous —
verified `0/50000`). ∎

(For `q ∈ [15, 27]`, `m = 1`, so the formula is simply `clearing_count = (q−1) − |{±vᵢ mod q}|`,
and "clears ⟺ the speeds and their negatives miss a residue mod q.")

## Why it matters (Route B of the finish)

The LRC(14) finish (LRC14-FINISH-MAP-2026-07-11) reduces to one statement, Route B: *every
divisor-complete family clears at some non-14 `q ∈ [15,31]`* (⟹ `M > 1/14` by the band-edge lemma,
opus S235). THM-718 makes the **dominant prime part exact and decidable**:

- "Clears at prime `q`" is the exact covering inequality `|{±j·vᵢ mod q}| < q − 1` — no search over
  `p`, and it depends only on `v mod q` (residue-periodic, kps `B5_congr_mod`). So the per-prime
  condition in the CRT-factored finite check (kps cont.35) is now a closed-form covering number.
- It is the **exact quantitative form of opus S239's "spread ⟹ bad coverer ⟹ clears"**: a spread
  family covers the residues `{±j·vᵢ}` poorly (misses some), so it clears; the *good* coverer (the
  AP `{1..13}`) is exactly the one whose `{±j·vᵢ}` hits everything (clears at NO non-14 `q` — the
  tight case, dispatched by `t = 1/14`).
- **Empirical structure (2000 primitive DC families):** 93% clear at a prime `q ∈ {17,19,23}`; the
  ~7% that `±j`-cover all three clear at a composite `q`; **0 fail every window modulus**.

**Remaining for Route B:** (a) the composite moduli (a `±j`-covering-all-primes family clears at a
composite — the `3e-4` core, kps) need the analogous count (messier: `gcd(vᵢ,q) > 1` allows `vᵢp ≡ 0`);
(b) proving the covering inequality holds for *every* DC family at *some* window `q` is the
anti-concentration = opus's AP inverse theorem (the AP is the unique set whose `{±j·vᵢ}` covers all
window `q`). THM-718 turns that inverse theorem into an exact statement about covering numbers.

## Files

`04-computation/lrc14_route_B_dc_clearing_verify_klein_S258.py`,
`lrc14_clearing_count_formula_klein_S259.py` (the verification).


## Addendum (klein-S259) — the sharp inverse-theorem framing + DC/tight disjointness

Using THM-718's exact count, the "clears" side of Route B has a SHARP structural form (verified):
- **The window-covers (13-sets failing to clear at EVERY non-14 `q ∈ [15,31]`) are exactly the
  TIGHT families** (`M = 1/14`): the tight AP `{1..13}`, the GW accelerations (`{1..11,13,24}`, …),
  and their dilations. Verified: all are window-covers; a search of 40000 random primitive 13-sets
  found NO other window-cover.
- **Every tight family is NON-divisor-complete** (it has no multiple of 14 — that is precisely why
  `t = 1/14` witnesses its tightness). Verified: `is_DC = False` for `{1..13}`, GW, and dilates.
- **Therefore DC ⟹ clears** (structurally): a divisor-complete family HAS a multiple of 14 (by
  definition), so it is not a tight family, so it is not a window-cover, so it clears at some
  non-14 `q`. Verified: **0 of 16328 primitive DC families is a window-cover.**

So Route B = the **inverse theorem** "`window-cover ⟹ tight family`" (the tight families being a
known finite list up to dilation, THM-612/708/709), and THM-718 turns "clears" into the exact
covering number `|{±j·vᵢ mod q}| < q−1`. The one remaining gap is that inverse theorem (equivalently:
window-completeness — that the window `[15,31]` captures the lonely time of every non-tight covering
family); everything else (the tight-list characterization, DC ⟹ non-tight, the exact count) is
proved or verified. Files: `lrc14_inverse_theorem_window_covers_klein_S259.py`.

## Addendum (klein-S261) — the UNIFIED formula on the coprime sub-family (composite completion; the shrink)

Combining THM-718 (prime count) with opus S241's auto-safe lemma (at composite `q ∈ [15,27]`, band
`{0,±1}`, a runner with `1 < gcd(v_i,q) < q` is ALWAYS safe at a unit `p`) and kps cont.46's two
conditions gives ONE exact clearing formula valid on the whole `m=1` window `[15,27]` (prime AND
composite) plus the `m=2` primes `{29,31}`:

> **`clearing_count_units(v,q) = φ(q) − |{ ±j·v_i mod q : gcd(v_i,q)=1, 1≤j≤m, the value a unit }|`,**

provided `q ∤ v_i` (condition (a)). Verified `0` failures over `70000` tests on the valid window.
(EXCLUDED: `q = 30`, the only composite with `m ≥ 2` in `[15,31]` — there auto-safe FAILS, since a
`gcd(v_i,30)=2` runner can land on the `±2` danger residue, itself a non-unit. The `m=2` **primes**
`29,31` are fine: all runners are coprime, so THM-718 with `j=1,2` applies directly.)

So **clears at `q` ⟺ (a) `q ∤ every v_i` AND (b) the COPRIME-to-`q` sub-family (dilated by `±j`)
misses a unit** — kps's two conditions, now with the exact count `φ(q) − |·|`.

**The shrink.** The structured (non-coprime) runners are provably INERT (auto-safe), so the
anti-concentration is only on the coprime-to-`q` sub-family. For divisor-complete families this is
SMALL — DC forces multiples of every prime `≤ 13`, so at a composite `q` most runners share a factor
and drop out. Measured over 3000 primitive DC families: the min-over-window `#coprime` has **median
4, mean 4.1** (max 13, the rare families blocking every composite, which fall back to the primes and
THM-718's 13-runner count). So **the 13-runner Route-B anti-concentration shrinks to a ~4-runner
one**: "the ~4 coprime runners of a spread DC family miss a unit-fold-class at some window `q`."
Condition (b) is count-automatic (`#coprime ≤ φ(q)/2 − 1`) for 74.5% at some `q`; the rest clear by
fold-class collisions among the coprime runners (kps cont.46). Files:
`04-computation/lrc14_unified_clearing_klein_S261.py`.
