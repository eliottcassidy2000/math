---
id: THM-403
title: The AP's floor witnesses are exactly the primitive n-th roots of unity (the (ℤ/n)*-orbit)
status: PROVED
source: opus-2026-06-03-S592
depends_on:
  - THM-369   # the 1/n clock witness
related:
  - HYP-2124  # static rigidity (S584)
  - HYP-2130  # exp/cyclotomic (S588)
---

# THM-403 — the cyclotomic witness orbit of the AP

## Statement

Let `AP = {1,2,…,n-1}` (the arithmetic-progression speed set, `δ=1/n`). Then:

**(a)** `t = j/n` is a witness (`‖v·j/n‖ ≥ 1/n` for all `v ∈ AP`) **iff `gcd(j,n)=1`.**

**(b)** Hence the AP's floor-witness set on the `n`-clock is *exactly* the `φ(n)`
**primitive `n`-th roots of unity** `{ e^{2πi j/n} : j ∈ (ℤ/n)^* }` — a single
`(ℤ/n)^*`-orbit, on which the multiplicative group acts (by `t ↦ u·t`) **simply
transitively**. In particular `M(AP) = 1/n`, attained on this finite (measure-zero,
rigid) orbit.

## Proof

**(a)** For `v ∈ {1,…,n-1}`, `‖v j/n‖ = (\,\min(r, n-r)\,)/n` where `r = vj \bmod n`; this
is `≥ 1/n` iff `r ≠ 0`, i.e. `n ∤ vj`. So `t=j/n` is a witness iff `n ∤ vj` for every
`v ∈ {1,…,n-1}`.
- If `gcd(j,n)=1`: `n ∤ vj` for `1 ≤ v ≤ n-1` (else `n | vj` with `gcd(j,n)=1` forces
  `n | v`, impossible). Witness. ✓
- If `d = gcd(j,n) > 1`: take `v = n/d ∈ {1,…,n-1}`; then `vj = (n/d)j = n·(j/d)` is a
  multiple of `n`, so `r=0` and `t` is *not* a witness. ✗

Thus witness `⟺ gcd(j,n)=1`. ∎

**(b)** `{ j : gcd(j,n)=1 } = (ℤ/n)^*`, and `e^{2πi j/n}` is a primitive `n`-th root iff
`gcd(j,n)=1`. The AP is **unit-invariant** (`u·AP ≡ AP \bmod n` for every unit `u`, since
multiplication by a unit permutes `{1,…,n-1}` mod `n`), so if `t` is a witness then for a
unit `u`, `‖v·(u t)‖` ranges over `{‖(uv)·t‖} = {‖w t‖}` (same multiset), so `ut` is a
witness too. The action `t ↦ ut` carries `j/n ↦ uj/n`, and `(ℤ/n)^*` acts on `{ j/n :
gcd(j,n)=1 }` simply transitively (regular representation). The value at any such `t` is
`1/n` (the minimum is achieved by `v ≡ ±j^{-1}` mod `n`), so `M(AP)=1/n`. ∎

## Significance

This is the bedrock of the cyclotomic/rigidity/exp arc:
- **Static rigidity** (HYP-2124): the witness set is a single rigid `(ℤ/n)^*`-orbit
  (measure zero, the clock points) — symmetry pins `M=1/n`.
- **exp / helix** (HYP-2130): the witnesses are literally the *primitive roots of unity*,
  so the AP's loneliness is a cyclotomic statement (`Φ_n`, Galois group `(ℤ/n)^*`).
- **The break** is at non-units `gcd(j,n)>1` (a runner `n/gcd` lands on the observer) —
  the composite/2-adic seam (THM-398 C′; the `2q` apex; `n=14` ramifies at `2`).

## Verification

`lrc_rigidity_scaling_s584.py` / `lrc_hyperoperation_helix_s588.py`: the AP clock
witnesses equal `(ℤ/n)^*` for `n=6,7,8,12,14` exactly.

**Artifacts:** see S584/S588. Builds on THM-369. (Convergent with oracle-S581o's
χ-separates-Paley result.)
