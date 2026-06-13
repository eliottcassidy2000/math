---
source: opus-2026-06-06-S699 (HN deep dive, extending S685/S687/S703)
status: SYNTHESIS — resolves the sharp question (the LRC ±2π/3 is the WITNESS/F(T,x) face, not the depth PGF) and extends the repo's HN cluster (S685 Niven=Dehn, S687 field tower, S703 density rosette) with: (1) the z⁷−z KEYSTONE unification — the 7-coloring = ℤ[ζ₆]/(2+ζ₆)≅𝔽₇, and the HN upper bound 7 = the forbidden tournament value 7 = Φ₃(2) = N(2+ζ₆) = |PG(2,2)| (the SAME 7); (2) the HEEGNER tower conjecture (χ=3↦√−3, χ=4↦√−11, χ=5↦√−19; class number 1 = the rigidity to force a color); (3) the Niven=Dehn=field-escape=Cl₂(π/3) coherence. Many HN statements collected.
tags: [hadwiger-nelson, z7-z, eisenstein, 7-coloring, F_7, Phi3, PG(2,2), heegner, field-tower, niven, dehn, Cl2, lee-yang, witnesses, FTA, chromatic-plane]
---

# Hadwiger–Nelson: the z⁷−z keystone, the "same 7", and the Heegner tower

**Prompt (user):** focus the sharp open question; make as many Hadwiger–Nelson statements as you
can; explore the repo heavily, including seemingly-unrelated tangents; find the buried needle.

The repo already holds a deep HN cluster — **S685** (Niven = Dehn-triviality = lattice escape),
**S687** (the FTA field-tower: χ≤3 on `ℤ[ζ₆]`, χ=4 via `√−11`, `z⁷−z` keystone), **S703/THM-411**
(unit-distance density = a roots-of-unity rosette), and my **S599u/g/h** (`Cl₂(π/3)`, the
spectral/Cayley unification, the alternating group). This session resolves the sharp question and
extends the cluster.

## A. The sharp question, resolved

> **The LRC depth-PGF roots are NOT at `±2π/3`.** Computed (`…s699m.py`): the AP worry-set's
> `P(z)=Σp_k z^k` has nonzero-root angles `{±73°,180°}` (n=5), `{±54.5°,±146°}` (n=6) — not `±120°`.

So the depth PGF was the *wrong* polynomial. The LRC's cyclotomic `±2π/3` (Eisenstein) structure
lives in **two other places**, both genuinely cyclotomic:
- the **witness clock points** = the `n`-th roots of unity (THM-403): for `3∣n`, `t=1/3` (angle
  `2π/3`) *is* a worry-set lonely time;
- the **worry-set round tournament's `F(T,x)`** (the forward-edge polynomial) — its Lee–Yang zeros
  cluster at `±2π/3`, all on `|z|=1` for the regular `H=9` tournament (opus-S44).
Via THM-402 (worry-set = round tournament), the LRC's Eisenstein angle is the **witness/tournament
face**, and it matches HN's `W₆` (the `z⁶−1` hexagon, roots at `0,±60°,±120°,180°` — *including*
`±2π/3`). *The ±2π/3 is real and cyclotomic in LRC; it just isn't the depth-PGF's roots.*

## B. The z⁷−z keystone and the "same 7" (extending S687)

`z⁷−z = z(z⁶−1)` has two readings (S687): over `ℂ`, roots `= {0}∪{6th roots}` = the **wheel `W₆`**
(center + hexagon), the `χ=3` gadget (the triangular lattice's non-2-colorability); over `𝔽₇`,
roots `= 𝔽₇` (Fermat), the **7-coloring palette**. I sharpen the arithmetic side:

> **The hexagonal 7-coloring is coloring by `ℤ[ζ₆]/(2+ζ₆) ≅ 𝔽₇`.** Verified: `7 = N(2+ζ₆) =
> a²+ab+b²` at `(2,1)`, so `(2+ζ₆)` is an **Eisenstein prime of norm 7**; the Isbell 7-coloring's
> period sublattice has index 7 = `(2+ζ₆)ℤ[ζ₆]`, so the 7 colors are the residues mod that prime.

> **The "same 7."** `7 = N(2+ζ₆) = Φ₃(2) = |PG(2,2)| = ` the HN upper bound `χ(ℝ²)≤7` = the first
> **forbidden tournament Hamiltonian-path value** (HYP-2180). One number, five guises, all via
> `z⁷−z`: the geometric `W₆` (χ=3 obstruction), the arithmetic `𝔽₇` (the 7-coloring upper bound),
> the cyclotomic `Φ₃(2)` (the tournament phantom volume), the projective `PG(2,2)`, the Eisenstein
> prime `N(2+ζ₆)`. The HN *upper bound* and the tournament *forbidden value* are literally the same
> `7`.

## C. The Heegner field tower (sharpening S687's conjecture)

S687: `χ≤3` on the cyclotomic `ℚ(√−3)`; `χ=4` (Moser spindle) needs `ℚ(√−11)` (the rotation
`cosθ=5/6`, `e^{iθ}` a root of `3z²−5z+3`, disc `−11`, Mahler measure 3 — *non-cyclotomic*). I
note both discriminants are **Heegner** (class number 1: `d∈{1,2,3,7,11,19,43,67,163}`):

> **[CONJECTURE, sharpening S687] each chromatic step adjoins a new class-number-1 (Heegner)
> imaginary-quadratic rotation field:** `χ=3 ↦ √−3`, `χ=4 ↦ √−11`, and `χ(ℝ²)=5 ↦ √−19` (the next
> Heegner). **Class number 1 is the rigidity** — unique factorization of the rotation, so the
> rotation group it generates cannot collapse — needed to *force* a new color. de Grey's `χ≥5` graph
> would then be realized in `ℚ(√−3,√−11,√−19)` (or a tower of Heegner fields), and `χ(ℝ²)` = the
> number of independent Heegner rotations one unit-distance graph forces.

## D. Niven = Dehn = field escape = Cl₂(π/3) (coherence with S685/S599u/v)

The chromatic jump is one phenomenon under four names:
- **Niven** (S687): `cosθ=5/6` is rational but `θ/π` is irrational (Niven's theorem) — the rotation
  is *not a root of unity*.
- **Dehn / equidecomposability** (S599v, S685): a non-Niven angle has *nontrivial Dehn invariant* —
  the tetrahedron of dihedral angle `θ` is not scissors-congruent to a cube; `Cl₂(θ)` is its volume.
- **Field escape** (S687): the rotation leaves `ℤ[ζ₆]` for `ℚ(√−11)`.
- **`Cl₂(π/3)`** (S599u): the Eisenstein/HN constant `1.0149` is the Dehn invariant / volume of the
  ideal regular tetrahedron (dihedral `π/3`).

> **The HN chromatic obstruction = the Dehn-invariant nontriviality of the rotation = leaving the
> cyclotomic field = Niven failure.** `χ=3` is the Dehn-trivial (cyclotomic, root-of-unity, `Cl₂`-
> flat) floor; each higher color needs a Dehn-nontrivial (non-cyclotomic, Heegner) rotation.

## E. Many HN statements (the collection)

1. `χ(ℤ[ζ₆]) ≤ 3` — every unit-distance graph on the Eisenstein lattice is 3-colorable (S687).
2. `χ(ℚ²)=2` (rational plane bipartite); the chromatic number climbs with the number field.
3. `χ=4` requires a non-cyclotomic (Niven-failing, Dehn-nontrivial) rotation; smallest is `√−11`
   (Moser spindle) (S687).
4. `z⁷−z` keystone: `ℂ`-roots = `W₆` (χ=3 gadget, incl. `±2π/3`); `𝔽₇`-roots = the 7-coloring.
5. 7-coloring = `ℤ[ζ₆]/(2+ζ₆) ≅ 𝔽₇` (color = residue mod the norm-7 Eisenstein prime) **[new]**.
6. HN upper bound `7` = forbidden tournament value `7` = `Φ₃(2)` = `N(2+ζ₆)` = `|PG(2,2)|` (same 7) **[new unification]**.
7. Spectral/Hoffman: triangular lattice `χ≥3` (tight); plane Bessel-`J₀` bound `χ≥3.48` (S699g).
8. `χ_f(ℝ²)=1/m₁≈4.36`; the integrality gap `χ>χ_f` *is* the LRC Vitali wall (S699g).
9. Density quantization (THM-411/S703): unit-distance density is a roots-of-unity rosette; the
   triangular `6`-rosette ↔ the Eisenstein `π/3`.
10. Heegner tower conjecture (C above): `χ=k ↦ √−d_k`, `d_k` Heegner; `χ=5 ↦ √−19` **[new]**.
11. The LRC worry-set ↔ round tournament ↔ `F(T,x)` Lee–Yang at `±2π/3` ↔ the `W₆` (THM-402/403,
    S44) — the LRC and HN share the Eisenstein angle through the witnesses **[A above]**.
12. `A₅ ≅` icosahedral: the unit-distance Cayley graph on `A₅` is a *spherical* HN problem (S699h).

## F. Honest status

- **Verified:** depth-PGF roots not at `±120°` (sharp question); `7=N(2+ζ₆)` Eisenstein prime;
  `z⁷−z` hexagon roots include `±2π/3`; the Heegner list; the Moser `√−11`/`3z²−5z+3` (S687).
- **New (mine):** the 7-coloring `= ℤ[ζ₆]/(2+ζ₆) ≅ 𝔽₇` reading; the five-guise "same 7"
  unification (HN-bound = forbidden-value = `Φ₃(2)` = `N(2+ζ₆)` = `|PG(2,2)|`); the sharp-question
  resolution (`±2π/3` = witnesses/`F(T,x)`, not the PGF); the Heegner-tower sharpening (`χ=5↦√−19`,
  class-number-1 = rigidity); the Niven=Dehn=`Cl₂(π/3)`=field-escape coherence.
- **Conjectural:** the Heegner tower (`χ(ℝ²)` = # Heegner rotations forced) — speculative, extends
  S687's field-tower conjecture; the honest verified content is the single-field realizations.
- **Built on (credit):** S685 (Niven=Dehn), S687 (field tower, `z⁷−z`), S703/THM-411 (density
  rosette), S599u (`Cl₂(π/3)`), S699g/h (spectral/Cayley/icosahedral), opus-S44 (`F(T,x)` Lee–Yang),
  THM-402/403 (round/cyclotomic witnesses), HYP-2180 (forbidden value `Φ₃(2)`).

**Artifacts:** `04-computation/hn_z7z_keystone_and_sharp_question_s699m.py` (+`.out`). New: **HYP-2277**.
