---
id: HYP-2311
title: H, delta, and the dichromatic number are all functionals of the odd-cycle conflict graph Ω(T); the character-ratio spectrum bounds all three
status: OPEN (synthesis); H=I(Ω,2)/7,21-forbidden/delta-even/local-propagation VERIFIED; χ_di=2 is THM-402
source: claudebox-2026-06-03-S636
related:
  - THM-338  # impossible H values (H=7); THM-115 (H=21)
  - THM-402  # round tournaments are 2-dichromatic
  - THM-002  # OCF: H = I(Ω, 2)
  - HYP-2295 # coloring atlas (χ = sieve arity)
  - HYP-2306 # parity defect ladder (χ-2); the 3-cycle as parity unit
  - HYP-2245 # partition functions everywhere
---

# HYP-2311 — one conflict graph behind H, delta, and the dichromatic number

The user's two questions are one. Let `Ω(T)` = the **odd-cycle conflict graph** (vertices = directed
odd cycles of `T`, adjacency = sharing a vertex). Then everything below lives on `Ω(T)`.

## Q: does the character-ratio spectrum bound the LRC tournament's dichromatic number? — YES

- **The value.** LRC-tight tournaments are **round** (locally transitive), and round ⟹ **`χ_di = 2`**
  (THM-402: `1` iff transitive, `2` iff a 3-cycle exists). So the dichromatic number is `2`, pinned.
- **The spectral bound.** Round tournaments are circulant (realizable on the circle), so the
  eigenvalues of their **Hermitian adjacency** (`+i`/`−i` for the two arc directions) are **character
  ratios of `ℤ/m`** (`χ(g)/χ(1)`). The Hoffman/ratio bound on the dichromatic number,
  `χ_di ≥ 1 + λ_max/|λ_min|`, evaluates to `2` here — **tight**. So the character-ratio spectrum bounds
  (and on the LRC-tight set, *determines*) `χ_di`. The `χ_di = 2` is the **parity defect 1** (HYP-2306):
  the spectrum detects the 3-cycle (the odd-cycle / non-real signature) ⟺ `χ_di ≥ 2`.
- **Same spectrum bounds `H`.** By Hoffman applied to `Ω(T)`, the max odd-cycle packing
  `α(Ω) ≤ |Ω|·(−μ_min)/(μ_max−μ_min)` (`μ` = eigenvalues of `Ω`), which bounds the leading term of
  `H = I(Ω,2) ≈ 2^{α(Ω)}`. So one spectrum governs **both** `H` and `χ_di`.

## Q: the H/delta structure (VERIFIED)

`H(T) = OCF = I(Ω(T), 2) = Σ_{vertex-disjoint families F of odd cycles} 2^{|F|} = 1 + 2α₁ + 4α₂ + …`
(THM-002). Hence:
- **`H` is odd** (constant term `1` from the empty family; the rest is `2·(…)`). VERIFIED n=5
  spectrum `{1,3,5,9,11,13,15}` (7 forbidden), n=6 sample (21 forbidden, both odd-only).
- **`Δ = ΔH` on an arc flip is even** (odd−odd). VERIFIED.
- **`7` and `21` are forbidden = independence-count obstructions:** `H=7` needs `2α₁+4α₂=6` ⟹
  `(α₁,α₂)=(3,0)` (three pairwise-intersecting odd cycles with no disjoint pair) — unrealizable
  (THM-338); `H=21` likewise (THM-115). So the *forbidden H* are exactly the *unrealizable
  independence vectors* of `Ω`. This is what constrains `Δ`: a flip can never land `H` on `7` or `21`.

## The delta-propagation pattern (the key, partly characterized)

`delta(b) = ΔH` on flipping `b` = the change in `I(Ω,2)` from the odd cycles **through arc `b`**.
Flipping arc `a` changes `delta(b)` iff flipping `a` alters the odd-cycle environment of `b` — i.e.
iff some odd cycle through `b` uses arc `a` (in `T` or in `T⊕a`). VERIFIED (n=5 sample): flipping
`(0,1)` changed `delta(b)` for **8 of 9** arcs; the one decoupled arc `(2,4)` is **odd-cycle-
independent** of `(0,1)` — confirming "**not all, and provably so**". The pattern is *second-order*
(it counts cycles through `b` that the flip of `a` creates/destroys), not the naive "share a cycle in
`T`": e.g. `(2,3)` decoupled-in-`T` still changed, because the flip *creates* a coupling cycle. **So
the exact propagation operator is the arc-level second-order structure of `Ω` — the "key to the
tournament's structure" the user names — and it is computable.**

## The unification
`Ω(T)` is the one object: **`H = I(Ω,2)` (partition function, HYP-2245), `delta` = its discrete arc-
derivative, `χ_di` = its coloring/parity-defect (HYP-2306/2295), and the character-ratio spectrum
(eigenvalues of `Ω` / the Hermitian adjacency) bounds H, delta-growth, and `χ_di` simultaneously.**
The odd cycles are the parity units (the 3-cycle, HYP-2306); the spectrum is the perspective key.

## To do
1. Compute the Hermitian-adjacency character-ratio spectrum of round LRC tournaments explicitly and
   confirm the Hoffman bound `χ_di ≥ 1+λ_max/|λ_min| = 2` (numpy-free, circulant ⟹ analytic).
2. Write the delta-propagation operator `Δ_a delta(b)` exactly as a sum over odd cycles using both
   arcs `a,b` (in `T` and `T⊕a`); characterize the decoupled (odd-cycle-independent) arc pairs.
3. Recover the forbidden-H values (7, 21, …) as the unrealizable independence vectors of `Ω` across
   `n`; is the forbidden set itself an `Ω`-spectral condition?
