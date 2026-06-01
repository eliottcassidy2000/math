---
source: oracle-2026-06-01-S529b
status: synthesis + rigorous framing (complements concurrent oracle-S529)
tags: [LRC, regular-polygon, dihedral, diagonals, inside-debt, resonance, n3-proof, S526]
---

# The outside vs the hidden inside: a polygon's diagonals are its LRC difficulty, and the triangle has none

**Prompt (user):** a tournament is a binary relation on a simplex, but also a
regular polygon (dihedral D_n); how does the OUTSIDE relate to the hidden INSIDE
arcs, and how does that apply to LRC?

**Note on convergence.** A concurrent session, `oracle-2026-06-01-S529`
(`04-computation/lrc_polygon_inside_outside_s529.py`), built the full computation
for this prompt: the skip-shell split (sides = outside, diagonals = inside), the
regular-polygon tournament with Gauss sum `|Σ χ(k)ω^k| = √p` (the QR/Paley
character that balances the inside), and — crucially — the **LRC covering sum
graded by resonance order**, defining the **"inside debt" = Σ_{r≥3} (resonance
orders)** and finding it is `0` at n=3 and "born at n=4." It also checked the
order-2 term against my S526 Legendre closed form (match). This note does not
duplicate that; it adds the one rigorous statement those numerics point at.

## The outside/inside dictionary

Place the `n` vertices as a regular polygon. Then:
- **OUTSIDE** = the `n` sides (skip-1 chords) = the boundary `n`-cycle = the cyclic
  order / ranking (the base path & cut/score space of the tiling model, S524).
- **INSIDE** = the diagonals (skip ≥ 2 chords) = the genuinely cyclic content.
- **D_n** acts: rotation = cyclic shift of the gaps; reflection = reverse them =
  complement `T ↔ T^op`.

Two facts make this sharp:
1. **The inside is enslaved to the outside.** For points on a circle the every
   diagonal is a *determined function* of the outside gap sequence. The tournament
   has `n` gap degrees of freedom, not `C(n,2)`. The realizable image is the
   **dihedral necklace `A000016`** (round tournaments, S523) — a vanishing slice of
   `A000568`. ("Arcs are not independent," S524, is exactly this enslavement.)
2. **In the LRC harmonic measure** (S526) `|SAFE| = Σ_{Σ k_i s_i = 0} Π g_{k_i}`,
   the grading by resonance order = grading by chord depth:
   - order 0 = the **outside** mean-field `(1-2/n)^{n-1}` (opus-S524 independence);
   - order 2 = the **sides** (pairwise) — a single character sum;
   - order ≥ 3 = the **inside diagonals** = the "inside debt" = multi-way correlations.

## The clean geometric theorem (the WHY)

> **A polygon has `n(n-3)/2` diagonals, which is `0` iff `n = 3`.** The triangle has
> no inside. An LRC "inside-debt" term (a `≥3`-term resonance `Σ k_i s_i = 0` with
> `≥3` nonzero `k_i`) requires `≥3` runners, i.e. `n ≥ 4`. So the inside debt
> vanishes identically for `n = 3` and is first possible at `n = 4`.

Therefore:
- **n = 3 (the triangle): no diagonals ⇒ no inside debt ⇒ |SAFE| is the outside +
  sides only ⇒ the single Legendre sum of S526 ⇒ PROVED.** The only "inside" a
  triangle could have is its one cyclic class (the 3-cycle), and that is already
  carried by the order-2 term.
- **n = 4 (the square): the first diagonal appears — the diameter.** The diameter
  is the antipodal chord = the dihedral reflection axis = the half-turn WALL
  (S522: even-gon = wall). Its arc is the order-≥3 resonance. The inside debt is
  born exactly here, and so is the open case.

So the user's image resolves precisely: **the OUTSIDE (sides, ranking, mean-field)
is always harmless — it gives the positive `(1-2/n)^{n-1}` that never vanishes. The
HIDDEN INSIDE (diagonals, ≥3-term resonances) is the entire LRC difficulty.** The
conjecture is the statement that the inside debt never overpays the outside credit.
The reason small cases are easy and `n ≥ 4` is hard is not arithmetic accident: it
is that **the triangle has no diagonals and the square does.**

## Where this points

- The first open case, n=4, is the *first diagonal* = the diameter = the dihedral
  reflection. Its inside-debt is a single new 3-term character sum (n=4 uses only
  odd harmonics). Evaluating/bounding it is the exact next step (HYP-2004 / the
  concurrent S529 "inside debt"): the analogue of the S526 Legendre evaluation for
  the square's diameter.
- The Gauss-sum `√p` the concurrent session found is the inside-balance of the
  *Paley* (QR) regular tournament — yet S522 showed the LRC-reachable regular class
  is the *rotational* `R_m`, not Paley. So the inside has two regular fillings
  (rotational vs QR); LRC reaches the rotational one. Worth reconciling.

## Anchors
Concurrent `04-computation/lrc_polygon_inside_outside_s529.py` (+ `.out`) — the
graded covering sum and inside-debt. This note:
`05-knowledge/results/lrc_polygon_inside_outside_triangle_s529b.out` (diagonal
counts). Builds on S522 (roots of unity), S523 (round=A000016), S524 (arcs not
independent), S526 (n=3 proof), HYP-2004.
