---
source: opus-2026-06-06-S699 (de Grey field-set question)
status: HONEST CORRECTION + refinement + cross-problem unification. de Grey's χ=5 graph does NOT provably force ≥3 distinct √−d — the field COUNT is the wrong invariant (Minkowski sums preserve the field; by Hilbert 90 the unit-norm group |β|=1 of ONE imaginary-quadratic field is INFINITE rank, so independent rotations pack within fixed fields; de Grey plausibly lives in ℚ(√−3,√−11)). Corrects S699m's Heegner-field-count conjecture. The RIGHT invariant is the ROTATION RANK = # multiplicatively-independent non-torsion unit-norm rotations; χ = 3 + rank (cyclotomic floor 0, Moser 1, de Grey ≥2). Shared object across HN/UD/LRC = unit-norm elements: LRC = torsion (cyclotomic clock), UD = ball-count (Sawin CM), HN = rank (de Grey).
tags: [hadwiger-nelson, de-grey, unit-norm, hilbert-90, rotation-rank, cyclotomic-floor, torsion, CM, sawin, LRC, unit-distance, heegner, correction]
---

# The unit-norm group, the rotation rank, and the cyclotomic floor

**Prompt (user):** does de Grey's graph provably force ≥3 distinct √−d? How does this relate to LRC
and unit distance?

## 1. Honest answer: NO — and the field count is the wrong invariant

> **de Grey's χ=5 graph does not provably force ≥3 distinct imaginary-quadratic fields.** Two
> reasons, both structural:
> - **Minkowski sums preserve the field.** de Grey's graph is a Minkowski sum of rotated copies of a
>   spindle-type gadget; sums of points in a field `K` stay in `K`. Only the *rotations* can enlarge
>   `K`, and they enlarge it once (to include the rotation), not per copy.
> - **The unit-norm group of one imaginary-quadratic field is INFINITE rank** (Hilbert 90:
>   `β = γ/γ̄` has `|β|=1` for every `γ`). **Verified:** `8` distinct *non-torsion* unit-norm
>   rotations (`cos = 5/6, −0.467, −0.10, 0.185, …`, none a root of unity) **all inside the single
>   field `ℚ(√−11)`**. So you can force *many independent rotations without any new field.*

de Grey plausibly lives in `ℚ(√−3, √−11)` (the Moser-spindle field), 2 imaginary-quadratic fields.
**This corrects my own S699m conjecture** ("χ = 2 + #Heegner fields"): the field *count* is bounded
and is not what grows with χ.

## 2. The right invariant: the rotation rank (χ = 3 + rank)

> **The growing invariant is the ROTATION RANK** `= ` the number of *multiplicatively-independent
> non-torsion unit-norm elements* (`|β|=1`, not roots of unity) the graph forces. The torsion
> (roots of unity) is the lattice; each independent non-torsion rotation is one "lattice escape."

- `χ = 3`: **rank 0** — only roots of unity (the cyclotomic lattice `ℤ[ζ₆]`); `χ ≤ 3` there (S687).
- `χ = 4`: **rank 1** — the Moser spindle forces *one* non-torsion rotation `(5+√−11)/6 ∈ ℚ(√−11)`,
  `cos = 5/6`, Niven-irrational, minimal poly `3z²−5z+3` (verified). Provable: `χ≥4 ⟹ rank ≥1`
  (S687: `χ≤3` on rank-0).
- `χ = 5`: de Grey forces **rank ≥ 2** *if* `χ = 3 + rank`. The sharp open question is exactly this:
  **does de Grey force a *second* multiplicatively-independent non-torsion rotation, or does it use
  only powers of the single Moser rotation (rank 1)?** That is the right reformulation of "≥3
  fields" — and it lives *within* `ℚ(√−3,√−11)`, not in new fields.

> **Refined conjecture: `χ(ℝ²) = 3 + (max non-torsion unit-norm rotation rank of a finite UD
> graph)`.** Narrowing `{5,6,7}` ⟺ bounding this rank (rank ≤ 2 ⟹ χ = 5; rank can reach 3 ⟹ χ ≥ 6).

## 3. The cross-problem unification: unit-norm elements are the shared object

All three problems are about **unit-norm elements `|β|=1` of imaginary-quadratic / CM fields**, read
three ways:

| problem | what it uses | reading of the unit-norm group |
|---|---|---|
| **LRC** | witnesses = `n`-th roots of unity (THM-403) | the **torsion** (roots of unity) = the **cyclotomic floor**; hardness = the *torsion arithmetic* (shells mod `2n−1`, prime-3) — rank 0 |
| **unit distance** | Sawin's `n^{1.014}` | the **count** of unit-norm elements in a ball (CM field + class-field tower) = the non-torsion group's **density** |
| **Hadwiger–Nelson** | de Grey / Moser | the **rank** of forced non-torsion rotations; `χ = 3 + rank` |

> **The escape from the cyclotomic floor = using non-torsion unit-norm elements.** LRC is the
> *torsion* member (the clock = roots of unity; its difficulty is the pure torsion arithmetic). UD
> and HN are the *non-torsion escape* members: UD measures the **density** of unit-norm elements
> (Sawin's CM construction beats the kissing-6 lattice floor), HN measures their **independent rank**
> (de Grey beats the χ=3 lattice floor). Same object, three functionals: *torsion clock* (LRC),
> *ball-count* (UD), *independent-rank* (HN).

This also re-reads THM-411/S703 (unit-distance density quantized by a roots-of-unity *rosette*): the
rosette **is** the torsion subgroup; the lattice density is the torsion count; Sawin's escape adds
non-torsion. And it re-reads S699g's "χ > χ_f = Vitali wall": the fractional/measure bound sees only
the *torsion average*; the combinatorial jump is the *non-torsion rank* — invisible to measure.

## 4. The LRC ↔ HN bridge made concrete

> The LRC worry-set lives at the **torsion** (clock points = roots of unity, THM-403); its
> hardness is that the torsion is *maximally tight* (the cyclotomic shell arithmetic, `n=14` =
> prime-3). HN's hardness is that the **non-torsion rank** grows (de Grey). The two are the
> *torsion* and *non-torsion* faces of the unit-norm group. **A unit-distance graph forcing high
> non-torsion rank is the HN analog of the LRC worry-set being torsion-tight** — and the Moser
> rotation's incommensurability (S699g's 33.56°) is the *non-torsion* witness, the exact analog of
> the LRC's irrational-speed resonances.

## 5. Honest status

- **Verified:** the Moser rotation `(5+√−11)/6` is a non-torsion unit-norm element of `ℚ(√−11)`
  (`|β|=1`, `cos=5/6`, not a root of unity); `ℚ(√−11)` has many independent non-torsion unit-norm
  rotations (Hilbert 90, infinite rank).
- **Established (rigorous):** Minkowski sums preserve the field; the unit-norm group is infinite
  rank ⟹ the field *count* is the wrong invariant; `χ≥4 ⟹ rank ≥1` (from S687).
- **Honest correction:** S699m's "χ = 2 + #Heegner fields" is too naive — replaced by `χ = 3 + rotation rank` (rank within fixed fields).
- **Open (the sharp reformulation):** does de Grey force rank ≥ 2 (a *second* independent
  non-torsion rotation)? — the correct version of "≥3 fields." I cannot resolve it without de
  Grey's explicit vertex set; it is the right finite/arithmetic target.
- **New (the unification):** unit-norm elements as the shared object; LRC=torsion, UD=count,
  HN=rank; the cyclotomic-floor escape = non-torsion; the re-readings of THM-411 and the Vitali wall.

**Artifacts:** `04-computation/degrey_rotation_rank_unitnorm_s699o.py` (+`.out`). Builds on S687
(field tower), S699m (Heegner — corrected here), S703/THM-411 (rosette), S699g/n (spectral / Vitali
wall), THM-403/404 (cyclotomic witnesses / doubling), Sawin (CM/class-field tower), Hilbert 90. New:
**HYP-2279**.
