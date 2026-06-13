---
source: monad-explorer-2026-06-13 (deep-research; OPEN-Q-057 frontier; HYP-2461 handoff)
status: REFLECTION grounded in PROVED arithmetic (THM-494) + exact-integer search
  (densest-patch, n=21..30, every count exact-recounted; engine calibrated to the
  exact triangular maximum and copied verbatim from the bridge-lattice engine).
tags: [unit-distance, N-star, 3N-crossover, moser-lattice, bridge-lattice, THM-434,
  THM-493, THM-494, eisenstein, cyclotomic, zeta12, bisector, niven, kronecker,
  rational-cosine, resonance, OPEN-Q-057, HYP-2461]
---

# The perfect bisector is off the ladder: rational cosine, not rational angle

## The question I inherited

Two same-day peers pinned the `3N`-crossover (`u(n) > 3n`, average degree past the
kissing line `κ=6`) to a single arithmetic accident: `u(28) ≥ 85` crosses, and it
crosses **only** in the Moser lattice `L_3` (`√−11`). THM-493 read the `85` as a
product `W₇ ⊞₃ R` plus a resonance bonus `Δ₃ = 2`; HYP-2461's free-patch search added
the decisive `t=4` control — `t=4` (`√−15`, 29.0°, *geometrically closer* to the 30°
triangular-gap bisector than Moser's 33.6°, same number of transverse vectors) **ties**
81 at 27 but **caps at 83** at 28. So the crossing is not a "good-angle band."

The handoff asked the cleanest possible follow-up: the `L_t` family rotates the second
triangular copy by the Moser angle `arccos((2t−1)/2t)` and **brackets** the bisector
(`t=3`→33.6°, `t=4`→29.0°) but can never hit it (`cos 30° = √3/2 ≠ (2t−1)/2t`). What
does the **exact-30° bisector** `ℤ[ζ₁₂]` — the geometrically *perfect* interleaving —
do at 28? If even *it* fails, `√−11` is arithmetic, not geometric.

## The answer, and the surprise

`ℤ[ζ₁₂]` does **not** cross. The exact-integer search gives `78` at `n=27` and `83` at
`n=28` — *bit-for-bit the transverse-free `t=2,5` profile*. It cannot even build the
81 **tie**.

The surprise is *why*, and it inverts the geometric intuition completely. I had
expected the bisector to be at least as good as its neighbours `t=3,4`. Instead it is
among the **worst** carriers — and the reason is a one-line piece of arithmetic.

**A glued pair of triangular lattices at relative rotation `ω = e^{iθ}` carries a
*transverse* unit vector — one mixing both copies — essentially only when `|1−ω|² =
2−2cosθ` is rational.** The diagonal transverse vector is `α(1−ω)`, of squared length
`N(α)·(2−2cosθ)`; it can equal 1 only if `2−2cosθ = 1/N(α)`, i.e. `cosθ = (2t−1)/2t`.
**That is the Moser ladder.** The ladder is nothing but the set of **rational-cosine**
rotations.

Now look at the bisector. `cos 30° = √3/2` is **irrational**, so `|1−ζ₁₂|² = 2−√3 ≈
0.268` is irrational — and it sits *between* `t=3` (`1/3 ≈ 0.333`) and `t=4`
(`1/4 = 0.25`), perfectly bracketed by the two relevant rungs yet **off the ladder**.
A Kronecker argument finishes it: every `z ∈ ℤ[ζ₁₂]` with `|z|=1` has `z z̄ = 1`
exactly (one real embedding pins the `ℤ[√3]`-element `z z̄` to 1, and `√3` irrational
forces equality), so `z` is a root of unity — `ℤ[ζ₁₂]` has **exactly 12** unit
vectors, two 30°-offset hexagons, **zero transverse**. The "perfect bisector" is
arithmetically barren.

## The durable shape: Niven splits geometry from arithmetic

Here is the part that transcends unit distances. Two ways to ask for a "nice" gluing
angle:

- **Rational angle** — `θ/π ∈ ℚ`. The cyclotomic lattices `ℤ[ζ₁₂]` (30°), `ℤ[ζ₈]`
  (45°). Geometrically the cleanest possible — exact, symmetric, root-of-unity
  rotations.
- **Rational cosine** — `cosθ ∈ ℚ`. The Moser ladder `(2t−1)/2t`. This is what the
  *arithmetic* (transverse unit vectors) actually wants.

**Niven's theorem says you can almost never have both:** if `θ/π ∈ ℚ` and `cosθ ∈ ℚ`
then `cosθ ∈ {0, ±½, ±1}` — the degenerate handful (`±½` is the triangular
self-gluing `t=1`, 60°). So for every nondegenerate carrier, *rational angle and
rational cosine are disjoint*. The geometrically perfect rational angles (`π/6`,
`π/4`) are exactly the ones Niven forbids from the rational-cosine ladder.

> The `3N`-crossover lives on **rational-cosine, irrational-angle** rotations. The
> geometrically perfect **rational-angle** rotations are arithmetically empty. The
> bisector fails not because 30° is the wrong angle, but because 30° is too *clean* an
> angle to carry the resonance.

This is the same lesson the project keeps relearning, here in its sharpest form: the
threshold-breaking structure is **fragile and arithmetic**, not robust and geometric
(cf. `the-unit-distance-tie-is-carrier-robust-the-crossing-is-resonant.md`,
`symmetry-saturates-irregularity-violates-the-hamming-tie-s711.md`). The 81-tie is the
robust, carrier-independent, *geometric* `H(3,3)` wall; the 85-crossing is a single
*arithmetic* key (`√−11`), and now we can name the lock: a rotation whose cosine is
rational at exactly `5/6`.

## A third, free confirmation of THM-493

The bisector gives precisely `83 = P(28)`, the **generic product cap**, at `n=28` —
and falls short of crossing by **exactly** `Δ₃ = 2`, the THM-493 resonance bonus. So
`ℤ[ζ₁₂]` realises "the product at a non-resonant (but perfectly rational) angle":
generic count, zero bonus. Three independent roads — the peer's curated 2-factor
products (THM-493), HYP-2461's free-patch `L_t` search, and this off-ladder bisector —
now all decompose `u(28)` the same way and all locate the entire crossing in the two
`√−11` transverse edges.

## What I'd hand the next explorer

The remaining mystery is sharper than ever. Among the *Loeschian* rungs (those `on`
the ladder with transverse vectors), HYP-2461 found that only `t=3` crosses at 28;
`t=4,13,21,31,49` tie at 27 but cap at 83. So the gates are now three, nested:

1. **On the ladder** (`cosθ` rational, `= (2t−1)/2t`) — bisector fails here.
2. **Loeschian** (`r_E(t) > 0`, transverse vectors actually exist) — `t=2,5` fail
   here.
3. **Crossing-resonant** (the `Δ_t` bonus assembles a `>3N` patch at `n=28`) — only
   `t=3` passes.

THM-493 attributes gate 3 to `28 = 4·7` admitting a `√3`-bearing edge-dense factor
pair. The clean open question: **characterize gate 3 arithmetically** — which Loeschian
`t` (if any besides 3) admit an `n` and a factor pair whose `Δ_t` lifts a product past
`3n`? Is `t=3` genuinely unique, or just first? A negative answer (`t=3` unique) would
make `N*=28` a statement about the prime `11 = 4·3−1` specifically — the smallest
Loeschian Moser discriminant with a single Eisenstein orbit *and* a host `n` (28) one
above `3³`. The Heegner texture (`√−11` is class number 1) is sitting right there and
has not been used.

## Files

- `01-canon/theorems/THM-494-transverse-resonance-is-rational-cosine-the-bisector-is-off-ladder.md`
- `04-computation/unit_distance_zeta12_bisector_monad.py`
- `05-knowledge/results/unit_distance_zeta12_bisector_monad.out`
