---
id: THM-433
name: unit-vector-count-of-the-Moser-ladder-bridge-lattice-is-12-plus-rE(t)
status: PROVED (elementary, complete) + VERIFIED exact-integer t=2..31
date: 2026-06-07
session: monad-explorer-2026-06-07-S5
depends_on:
  - THM-432   # the Moser lattice = bridge ring; defines L_t = Z[zeta6, w_t]
  - HYP-2298  # the ladder w_t = ((2t-1)+i sqrt(4t-1))/(2t); the OPEN "characterize 6+6k" this CLOSES
resolves:
  - "HYP-2298 ADDENDUM open question (INDEX line ~8716): characterize k = #transverse-orbits"
---

# THM-433: the Moser-ladder bridge lattice has exactly `12 + r_E(t)` unit vectors

## Statement

Fix an integer `t ≥ 2` and let `d = 4t − 1`. Assume `d` is **not** three times a
perfect square (equivalently `ω_t ∉ ℚ(√−3)`), so that the **bridge lattice**

```
   L_t  =  ℤ[ζ₆] ⊕ ℤ[ζ₆]·ω_t  ⊂  ℂ,
   ζ₆ = (1 + i√3)/2   (|ζ₆| = 1, the triangular generator),
   ω_t = ((2t−1) + i√(4t−1)) / (2t)   (|ω_t| = 1, the t-th Moser angle, cos = (2t−1)/2t)
```

is a genuine rank-4 ℤ-module (a Moser-type lattice; `t=1` is the rosette ℤ[ζ₆],
`t=3` is the Moser lattice of THM-432). Then the number of **unit vectors** of
`L_t` — lattice points `z ∈ L_t` with `|z| = 1` — is **exactly**

```
        # units(L_t)  =  12  +  r_E(t),
```

where `r_E(t) = #{ α ∈ ℤ[ζ₆] : N(α) = α·ᾱ = t }` is the **Eisenstein
representation count** of `t`. Explicitly `r_E(t) = 6·(d_{1}(t) − d_{2}(t))`,
six times the excess of divisors of `t` that are `≡ 1 (mod 3)` over those
`≡ 2 (mod 3)` (the classical Loeschian / `x²+xy+y²` count).

Writing the count as `6 + 6k` (every unit vector lies in a free `⟨ζ₆⟩`-orbit of
size 6): **`k = 1 + r_E(t)/6`** = one plus the number of `ζ₆`-orbits of
norm-`t` Eisenstein integers.

### The three families of unit vectors
| family | vectors | count | condition |
|---|---|---|---|
| triangular rosette | `ζ₆^j` | 6 | always |
| `ω_t` rosette | `ζ₆^j · ω_t` | 6 | always |
| **transverse** | `α·(1 − ω_t)`, `N(α)=t` | `r_E(t)` | depends on `t` |

The transverse vectors exist because **`|1 − ω_t|² = 1/t` exactly**
(`1 − ω_t = (1 − i√(4t−1))/(2t)`, so `|1−ω_t|² = (1 + (4t−1))/(4t²) = 4t/4t² = 1/t`),
hence `|α(1−ω_t)|² = N(α)/t = 1 ⟺ N(α) = t`.

## Why it matters (the sharpening)

The count depends on `t` **only through the splitting of `t` in `ℚ(√−3)`** — the
*triangular* (first) CM field. The second glued CM direction `√−(4t−1)` plays **no
role in the count**. This corrects/sharpens the HYP-2298 guess that `k` is governed
by splitting in the *biquadratic* field `ℚ(√−3, √(4t−1))`: it is a `ℚ(√−3)`-only
invariant. Concretely:
- `r_E(t)=0` ⟺ `t` is **not** Loeschian (some prime `≡2 mod 3` divides `t` to an
  odd power), e.g. `t = 2,5,6,8,10,11,14,…` → **12 units** (`k=1`).
- `r_E(t)=6` ⟺ `t` Loeschian with a single `ζ₆`-orbit, e.g. `t = 3,4,9,12,16,25,27`
  → **18 units** (`k=2`). (The Moser lattice `t=3` sits here — its "18" is generic.)
- `r_E(t)=12` ⟺ `t` carries a **split** prime `≡1 mod 3`, e.g. `t = 13,21,28,31`
  → **24 units** (`k=3`). (`t=21` — the campaign number — is a 24-unit lattice.)
- In general `# units` takes every value in `12 + 6·{0,1,2,3,…}`; the max degree of
  a Moser-ladder unit-distance graph is unbounded along the ladder.

## Proof (complete, elementary)

**Setup.** `ℤ[ζ₆]` is the ring of Eisenstein integers `E` (a PID); complex
conjugation is its nontrivial automorphism, `N(α)=αᾱ`. Multiplication by `ζ₆`
preserves `L_t` (since `ζ₆²=ζ₆−1` and `ζ₆·(ζ₆ω_t)=(ζ₆−1)ω_t`), so the set `U` of
unit vectors is closed under the order-6 group `⟨ζ₆⟩`, which acts **freely** on `U`
(no `|z|=1` vector is fixed by a nontrivial rotation of order dividing 6 except via
`ζ₆^j z = z ⟹ z=0`). Hence `|U| = 6·#(orbits)`.

Write `z = α + β·ω_t` with `α, β ∈ E` (this is the `E`-module decomposition of
`L_t`).

**Step 1 — rationality forces `ᾱβ ∈ ℤ`.** Compute
`|z|² = N(α) + N(β) + 2 Re(γ ω_t)` with `γ = ᾱβ ∈ E`. Writing `γ = p + qζ₆`
(`p,q ∈ ℤ`), `Re(γ) = p + q/2`, `Im(γ) = q√3/2`, and
`Re(ω_t) = (2t−1)/(2t)`, `Im(ω_t) = √(4t−1)/(2t)`, one gets
```
   |z|²  =  [ N(α) + N(β) + (p + q/2)(2t−1)/t ]  −  (q/(2t))·√(3(4t−1)).
```
`N(α), N(β) ∈ ℤ`. Since `d=4t−1` is not `3·square`, `√(3(4t−1))` is irrational, so
`|z|² ∈ ℚ` forces **`q = 0`**, i.e. `γ = ᾱβ = p ∈ ℤ`.

**Step 2 — collinearity.** `ᾱβ ∈ ℤ ⊂ ℝ` means `β/α = ᾱβ/N(α) = p/N(α) ∈ ℚ`. So
(unless `α=0` or `β=0`) `α` and `β` lie on one ray `ℚ·α₀` through the origin, where
`α₀ ∈ E` is the primitive Eisenstein point on that ray.

**Step 3 — the two rosettes.** If `β = 0`: `|z|=|α|=1 ⟺ N(α)=1`, the 6 units
`ζ₆^j` (one orbit). If `α = 0`: `z = βω_t`, `|z|=|β|=1`, the 6 vectors `ζ₆^j ω_t`
(one orbit). Total `12`, always present.

**Step 4 — transverse vectors and the rigid factorization.** Suppose `α,β ≠ 0`.
By Step 2 write `λ = β/α = −u/w ∈ ℚ` in lowest terms (`w ≥ 1`). Since
`f(λ) := |1+λω_t|² = 1 + λ² + λ(2t−1)/t` has minimum `d/(4t²) > 0` at
`λ = −(2t−1)/(2t)`, and `|z|² = N(α)·f(λ) = 1` needs `f(λ) ≤ 1`, we have
`λ ∈ (−(2t−1)/t, 0)`, so `u > 0`. As `gcd(u,w)=1` and `−(u/w)α = β ∈ E`, we get
`w | α`, say `α = w·α₀`; then
```
   z = α(1 + λω_t) = wα₀(1 − (u/w)ω_t) = α₀·(w − u·ω_t),
   |z|² = N(α₀)·|w − u ω_t|² = N(α₀)·[ t(u²+w²) − (2t−1)uw ] / t.
```
The bracket is the binary form of discriminant `−(4t−1)`, and it **factors**:
```
   Q(u,w) := t(u²+w²) − (2t−1)uw  =  t(u − w)²  +  uw.
```
Thus `|z|²=1 ⟺ N(α₀)·[ t(u−w)² + uw ] = t`. Both factors are positive integers,
so `t(u−w)² + uw ≤ t`. But `uw ≥ 1`, and if `u ≠ w` then `t(u−w)² ≥ t`, giving
`t(u−w)²+uw ≥ t+1 > t` — impossible. Hence **`u = w`**, and `gcd(u,w)=1 ⟹ u=w=1`,
i.e. `λ = −1`. Then `N(α₀)·1 = t`, so `α := α₀` has `N(α)=t` and
`z = α(1 − ω_t)`.

Conversely every `α` with `N(α)=t` gives a unit vector `z=α(1−ω_t)` (Step 0 identity
`|1−ω_t|²=1/t`); distinct `α` give distinct `z`; and `ζ₆·(α(1−ω_t))=(ζ₆α)(1−ω_t)`
with `N(ζ₆α)=t`, so these `r_E(t)` vectors split into `r_E(t)/6` free `⟨ζ₆⟩`-orbits.

**Total:** `12 + r_E(t)`. ∎

## Verification

`04-computation/moser_ladder_unitcount_formula_monad_s5.py` (exact integer
arithmetic over the value-basis `{1, √3, √d, √3d}`, box-stable `R`): for every
non-degenerate `t = 2,…,31`, the exact unit count equals `12 + r_E(t)` AND the
both-nonzero count equals `r_E(t)` — 30/30, zero failures. Output:
`05-knowledge/results/moser_ladder_unitcount_formula_monad_s5.out`.

## Honest scope

- PROVED in full (the `Q(u,w)=t(u−w)²+uw` rigidity needs no computation); the
  computation is an independent exact check.
- Hypothesis `d ≠ 3·square` is necessary (else `ω_t ∈ ℚ(√−3)` and `L_t` is rank 2;
  these are `t = 1, 7, 19, 37, … = (3m²+1)/4`, `m` odd).
- This counts unit **vectors** (max possible vertex degree), not the densest
  unit-distance graph the lattice supports — that is a separate (harder) question.
