---
source: claude-2026-06-06-S684
status: computed (AG_n spectrum) + operational J0 plane bound + the spectral-floor/gadget-growth synthesis
tags: [alternating-group-graph, cayley-graph, nonabelian, distance-graph, spectral, hoffman, eisenstein, prime-3, chromatic, hadwiger-nelson, metagraph, de-grey, unification]
---

# The alternating group graph: a scale-invariant distance Cayley graph, and the spectral-floor / gadget-growth split

Extending the distance-Cayley-graph unification (HYP-2265): I ran the LRC
spectral machinery on the plane (operational test), then studied the alternating
group graph `AG_n` as a *nonabelian* distance Cayley graph and how it varies with
`n`. The result sharpens the whole unification into a clean dichotomy.

## Operational test: LRC machinery → χ(ℝ²)

The LRC spectral (Hoffman) bound `χ ≥ 1 − λ_max/λ_min`, with `λ` = the
connection-set Fourier transform, transcribed to the plane with `λ = J₀(2π r)`
(Bessel = the unit circle's transform), gives `λ_min = −0.4028`, independence
ratio `m₁ ≤ 0.287`, hence `χ(ℝ²) ≥ 3.48`. Literal code reuse of the LRC
distance-graph bound — so HYP-2265's unification is **operational**, not just a
framing.

## AG_n as a nonabelian distance Cayley graph

`AG_n` = Cayley graph of the alternating group `A_n` with the 3-cycle generators
`{(1,2,i),(1,i,2): i=3..n}`, degree `2(n−2)`. It is a distance Cayley graph on a
**nonabelian** group — its eigenvalues are the 3-cycle connection set evaluated on
`A_n`'s irreducible representations (representation-theory "Fourier" in place of
abelian characters). Computed for `n=4..7`:

- `λ₂ = deg − 2` ⟹ **spectral gap = 2, constant in `n`** (a poor expander);
- `λ_min = −deg/2 = −(n−2)` ⟹ ratio `λ_min/deg = −1/2`, **constant**;
- so the Hoffman bound `χ(AG_n) ≥ 1 − λ_max/λ_min = 3`, **constant for all `n`**;
- but the **true `χ` grows** (greedy upper bounds `3, 4, 6, 7` for `n=4..7`).

**Why `λ_min = −deg/2`:** the 3-cycle generators act on the irrep where each
generator is the cube root of unity `ω = e^{2πi/3}` (the **Eisenstein / prime-3**
root); there `2cos(2π/3) = −1` per generator pair, summed over the `n−2` pairs
gives `−(n−2) = −deg/2`. So the spectral floor is a **prime-3 phenomenon** — the
same `π/3`/Eisenstein root behind `Cl₂(π/3)` (HYP-707), the hexagonal 7-coloring,
and the LRC 3-shell.

## The dichotomy: spectral floor vs gadget growth

`AG_n`'s spectrum is **scale-invariant** (gap 2, ratio `−1/2`), so the
spectral/Fourier distance-Cayley bound is a **constant floor** (`χ ≥ 3`) — it is
**blind** to the chromatic growth. This is *exactly* the repo's tournament
**metagraph** phenomenon (`χ = n−1` while the Hoffman bound stays ~3,
chromatic-number-synthesis). The lesson, now sharp:

> The spectral / Fourier bound is the **shared floor** of the whole
> distance-Cayley family — it captures bounds where the connection-set transform's
> negativity *scales* (LRC distance graphs; `χ(ℝ²) ≥ 3.48`). But chromatic
> **growth** that is combinatorial — `χ(AG_n)`, `χ(metagraph)=n−1`, `χ(ℝ²)=5` —
> is **invisible to the spectrum** and requires **finite rigid gadgets**.

And the gadget method is *also* shared across the family:
- **de Grey's graph** forces `χ(ℝ²) ≥ 5` (a finite unit-distance gadget);
- **LRC tight configs** (AP, `V*`) pin the loneliness bound;
- the **pancyclic odd-cycle count** is what actually closed `H=21` (HYP-2258) —
  a "gadget" (the long cycles) the spectral/profile view missed.

So the unification has two layers: a **spectral floor** (Fourier transform of the
connection set, universal across LRC/HN/AG_n/metagraph) and a **gadget ceiling**
(finite rigid certificates, where the real growth lives). `AG_n` is the cleanest
specimen: its floor is *exactly constant*, so 100% of its `χ`-growth is gadget.

## Where AG_n sits in the web

- **Distance-Cayley frame:** nonabelian instance — extends LRC (`ℝ/ℤ`), HN (`ℝ²`),
  unit-distance (`ℤ[ζ₆]`) to `A_n`.
- **Prime-3 / Eisenstein:** `λ_min = −deg/2` via `ω` — the same root as
  everywhere.
- **Tournaments:** the 3-cycle generators ↔ tournament 3-cycles (OCF); the
  chromatic floor `3` is the triangle `id→g→g²`; `A_n` = even permutations = even
  Hamiltonian-path orderings (Rédei parity — the even/odd Ham-path split).

## New idea / next

1. **The floor-vs-gadget invariant.** Define, for a distance Cayley graph, the
   ratio `χ_true / χ_Hoffman` — `AG_n` and the metagraph have it `→ ∞` (gadget-
   dominated), LRC's AP has it `≈ 1` (spectrally tight). Classify the family by
   this ratio; the "hard" problems are the gadget-dominated ones.
2. **AG_n as a tournament-parity carrier.** Since `A_n` = even Ham-path orderings,
   is the signed Ham-path count (`Σ sign(π)` over a tournament's Ham paths)
   computed by an `AG_n`-flow? That would put the Rédei *parity* on the
   alternating group graph directly.
3. **Gadget transfer:** the H=21 gadget was "count the long odd cycles." Its
   AG_n analogue: the long cycles in `AG_n` (its girth is 3, but longer cycles
   carry the chromatic growth) — do they give a combinatorial `χ(AG_n)` lower
   bound beating Hoffman, the same way the pancyclic count beat the profile?
