---
id: THM-616
title: One tightener is useless at every scale — M(mU ∪ {w}) = M(U) for any m≥1 and any w with m∤w (when M(U)≤1/4; always ≥ min(M(U),1/4)). A single non-m-divisible runner cannot reduce the loneliness margin of an m-divisible family, because its orbit-max over {t+j/m} is ≥ ½−gcd(w,m)/(2m) ≥ 1/4 ≥ M(U), so it never binds at the argmax. Corollary: f=1 confinement for ALL m (M(mU∪{w}) ≥ 1/13 > 1/14 for e=12) — generalizes mac-mini's Lemma C (m=2 only) with a one-line orbit-max argument and a quantitative margin.
status: PROVED (elementary orbit-max; verified exactly for m=2..7, 0 violations over 200+ families).
source: opus-2026-07-04-S67
depends_on:
  - LRCUpTo13   # M(U) ≥ 1/(e+1): the e=12 corollary uses M(U) ≥ 1/13
related:
  - THM-612     # Lemma C (mac-mini, m=2 f=1) — this generalizes it to all m
  - THM-615     # the folding identity (m=2 f=2); this is the f=1 companion, clean because Ψ₁≥1/4 (no extremity)
  - HYP-4062    # tight-locus rigidity — confinement (this closes its f=1 slice for all m)
---

# THM-616 — one tightener is useless at every scale

**Setup.** `U` a finite set of positive integers, `m ≥ 1`, `w` a positive integer with `m ∤ w`. The
"tightener" `w` is added to the `m`-dilated family `mU = {mu : u∈U}`. Write `g(t) = min_{u∈U} ‖mu·t‖ =
min_u ‖u·(mt)‖`; it is `(1/m)`-periodic and `max_t g = M(U)`.

## Theorem (PROVED)
> **`M(mU ∪ {w}) = M(U)`** whenever `M(U) ≤ 1/4`; in general `M(mU ∪ {w}) ≥ min(M(U), 1/4)`.

### Proof
`M(mU∪{w}) = max_t min(g(t), ‖wt‖)`. Fix `t` and look at the `m`-orbit `{t + j/m : j=0..m-1}`: `g` is
constant `= g(t)` on it (`(1/m)`-periodic), so
`max_j min(g, ‖w(t+j/m)‖) = min(g(t), Φ(t))`, `Φ(t) := max_j ‖w(t+j/m)‖ = max_j ‖wt + jw/m‖`.
As `t` ranges over `[0,1)` the orbits tile it, so `M(mU∪{w}) = max_t min(g(t), Φ(t))`.

Now `{jw/m mod 1}` is the cyclic subgroup `⟨w/m⟩ ≤ ℝ/ℤ` of order `m/d` and spacing `d/m`, `d = gcd(w,m)`.
Since `m ∤ w`, `d` is a proper divisor of `m`, so `d ≤ m/2` and the spacing `d/m ≤ 1/2`. The coset
`wt + ⟨w/m⟩` therefore has a point within `d/(2m) ≤ 1/4` of `1/2`, giving
> `Φ(t) ≥ 1/2 − d/(2m) ≥ 1/4`   for every `t`.

At a `g`-argmax `t₀` (`g(t₀)=M(U)`): `min(g(t₀), Φ(t₀)) = min(M(U), Φ(t₀)) ≥ min(M(U), 1/4)`. Hence
`M(mU∪{w}) ≥ min(M(U),1/4)`, `= M(U)` when `M(U) ≤ 1/4`. And `M(mU∪{w}) ≤ M(mU) = M(U)` (adding a runner
never increases the margin). ∎ *(Verified exactly, m=2..7, 200+ families, 0 violations.)*

## Corollary (f=1 confinement, ALL m)
A tight family with hiding denominator `q* = 14m` and exactly **one** tightener has `e = 12` `m`-divisible
runners `E = mU` (`U` 12 runners) plus one `w`, `m∤w`. By LRC≤13, `M(U) ≥ 1/13`, so
> `M(mU ∪ {w}) ≥ min(M(U), 1/4) ≥ 1/13 > 1/14`.

So the family is **not** tight — no `q*=14m` (`m≥2`) tight family has a single tightener. This
**generalizes mac-mini's Lemma C** (THM-612, proved only for `m=2`) to *all* `m`, replaces its shift
obstruction with a one-line orbit-max argument, and adds the quantitative margin `1/13`.

## Significance and scope
- **The mechanism, sharply:** an `m`-divisible family lives on the coarse `(1/m)`-grid; a single
  non-`m`-divisible runner always has an orbit point ≥ `1/4` from the integers, so it can never be the
  binding runner at the family's best time. *Each tightener alone is harmless.*
- **Why f≥2 is different (honest):** two tighteners `w₁,w₂` can be jointly *extremal* — the folding
  quantity `Ψ = max(min(a,b), ½−max(a,b))` (THM-615) can vanish, unlike the single-tightener
  `Ψ₁ = max(‖wt‖, ½−‖wt‖) ≥ 1/4` here. So confinement's residual is precisely the **joint** reduction of
  two-or-more small tighteners — each individually useless, together able to reach `1/14` only on the
  near-AP arithmetic corner (THM-615 Lemma 3 + the tight-locus rigidity, HYP-4062). This theorem removes
  the entire `f=1` slice at every scale.
