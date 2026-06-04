---
id: THM-412
title: The twisted-shell dodge — a constructive looseness criterion on the 2n-1 shells
status: PROVED (the criterion); coverage of C′(14) is CONJECTURE (HYP-2197)
source: claudebox-2026-06-03-S622
depends_on:
  - THM-398   # the reduction LRC(n) <= C′; Lemma A (1-clock) and Lemma B (dominance dodge)
  - THM-401   # pair-sum sieve modulus 2n-1 (= the shell ceiling here)
related:
  - HYP-2197  # coverage conjecture: twisted-shell ∪ B′ = all multiple-of-n configs
  - THM-407   # CONVERGENT (opus-S599g, same prompt): the dual shell-FOLDING view —
              # G=⟨2,-1⟩ folds 13 shells→3 gcd-strata mod 27=3³ (WHY n=14 is hard).
              # THM-412 is the constructive DODGE (HOW to escape); THM-407 the orbit space.
  - THM-411   # tight sets witnessed on (ℤ/m)* — the same multiplier/perspective object
  - HYP-2130  # rigidity = orbit-type (the multiplier a = the "perspective"/twist)
---

# THM-412 — the twisted-shell dodge

## Statement (constructive looseness criterion)

Convention as in [[THM-398]]: `n` runners, speed set `S` of `n-1` distinct positive integers,
`M(S) = max_t min_{v∈S} ‖v t‖`, gap level `1/n`; `S` is **loose** if `M(S) > 1/n`.

> **THM-412.** Let `m ≥ 2` and `a` with `gcd(a,m)=1`. If for every `v ∈ S`
> ```
>     dist(a v mod m, 0) > m/n          (equivalently  n·min(r, m-r) > m,  r = a v mod m),
> ```
> i.e. every residue `a v (mod m)` avoids the **central danger band** `(m/n, (n-1)m/n)`,
> then `t = a/m` witnesses `M(S) > 1/n`, so `S` is loose.

The witness `t = a/m` is a **twisted clock**: the modulus `m` chooses a *shell*, the multiplier
`a ∈ (ℤ/m)*` is the *twist* (a "perspective", cf. [[HYP-2130]]) that rotates every runner off
the band.

## Proof

`‖v · a/m‖ = dist(av/m, ℤ) = min(r, m-r)/m` with `r = av mod m`. The hypothesis gives
`min(r,m-r)/m > 1/n` for every `v`, so `min_{v∈S} ‖v a/m‖ > 1/n`. The point `t=a/m` thus lies in
the open safe set `G(S) = {t : ‖vt‖ > 1/n ∀v}`, which is therefore nonempty (and open), so
`M(S) > 1/n`. ∎

## Scope and the shell ceiling

- **`m ≤ n-1`, `a` arbitrary:** the band excludes only `0`, so the criterion is "no `v ≡ 0 (mod m)`"
  — the twist is irrelevant; this is the strict-looseness **1-clock** (Lemma A of THM-398, sharpened
  to `>`). One clock per `m ∈ {2,…,n-1}`.
- **`n+1 ≤ m ≤ 2n-1`:** here `1/m < 1/n ≤ 2/m`, so the criterion needs **band-distance ≥ 2**
  (`av mod m ∉ {0, ±1}`) and the multiplier `a` is **essential**. The ceiling is exactly
  `2/m > 1/n ⟺ m ≤ 2n-1` — the **pair-sum-sieve modulus** [[THM-401]] (`2n-1 = 27 = 3³` at `n=14`).
  These are the **flow shells**.
- The achievable witness values `d/m` (band-distance `d`) form the realized "second-gap" ladder
  on the shells.

## Why it matters (the C′(14) frontier)

`M(S) ≥ 1/n` for configs with no multiple of `n` is Lemma A; the whole of LRC(14) is the residual
**C′(14)** — multiple-of-14 configs are loose (THM-398). THM-398's measure tool (Criterion B′,
dominance) leaves a `3.2%` "all-short" residual at `n=14`. **The twisted-shell dodge closes it**
empirically: over 6000 multiple-of-14 configs, `twisted-shell (m ≤ 27) ∪ B′` leaves **0 residual**,
and the dodge alone subsumes B′ on all near-tight configs (the `846` B′-failures are all
twisted-shell-covered; the only configs needing B′ are far-from-tight dominant ones, `m ≳ 29`).

At the critical shell `m = 2n-1 = 27 = 3³` the dodge has a clean form: a witness exists iff some
**unit ±-pair `{u,-u} ⊂ (ℤ/27)*`** (a **twisted-involution orbit**) is disjoint from `S mod 27`
(then `a = u^{-1}` works). The inner 3-adic shells (multiples of 3) are the degenerate "cross"
where the `±` involution has fixed points (the apex). See [[HYP-2197]].

## Verification

`04-computation/lrc_n14_flowshells_s622.py`, `lrc_n14_twisted_shell_s622.py`,
`lrc_n14_ceiling27_s622.py`, `lrc_n14_residual_empty_s622.py`
→ `05-knowledge/results/lrc_n14_*_s622.out`. Exact over ℚ; coverage over 6000+40000 configs.
