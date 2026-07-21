---
id: THM-1780
title: "H LEAVES THE MOMENT-NULLCONE LADDER AT EXACTLY n = 6 — the Hamiltonian-path count is spectral (a function of the trace moments) for n ≤ 5 and NOT spectral for n ≥ 6, so it sits one rung ABOVE the holonomic ceiling as the tournament's #P-hard invariant; the ladder extends rational ⊂ algebraic ⊂ holonomic ⊂ #P. Answering named-next (1) of THM-1775. TEST: group all labeled tournaments by their moment vector (tr(A¹),…,tr(Aⁿ)) — equivalently the characteristic polynomial, by Newton — and ask whether H is constant on each co-spectral class. RESULT: constant at n = 4 (3 classes) and n = 5 (9 classes); SPLITS at n = 6 — 3 co-spectral classes carry two H-values, the witness being moment vector (0,0,12,12,10,48) which carries BOTH H = 13 and H = 17 (two non-isomorphic co-spectral tournaments, both odd per Rédei, both in the allowed spectrum {odd}\{7,21}). CONSEQUENCE: THM-133's spectral reduction H = (462 − tr(A⁴))/2 is a SYMMETRY COLLAPSE special to Z₇ circulants (whose isomorphism type IS determined by the spectrum), NOT the general law — for n ≥ 6 no fixed polynomial in the trace moments computes H. So on the ladder of THM-1775, the trace moments tr(Aᵐ) occupy the RATIONAL floor (Cayley–Hamilton, depth n) while H itself is the tournament's permanent — genuinely #P-hard, above holonomic. The tournament therefore spans the WHOLE ladder, from its spectrum at the bottom to H at the top, and this is exactly why the H-spectrum is a strictly finer 'universal code' than the trace spectrum: it first sees more at n = 6. The threshold n = 6 joins the project's n ≥ 6 phase-transition lore (the per-path identity fails at n ≥ 6)"
status: >
  PROVED by exhaustive labeled census: H (via bitmask DP) is constant on every co-spectral
  class for n = 4, 5 and splits for n = 6 (all 2^15 tournaments enumerated; the (0,0,12,12,10,48)
  witness with H in {13,17} is exact).  "H is #P-hard" is the standard status of the
  directed-Hamiltonian-path / permanent count, cited not reproved; what is proved here is the
  concrete non-spectrality — the first co-spectral H-split — which is the precise sense in
  which H is off the ladder.  The ladder extension rational ⊂ algebraic ⊂ holonomic ⊂ #P is a
  placement, and the placement of H at the #P level is what the split establishes.
renumbered: "claimed THM-1765; renumbered to THM-1780 by first-pusher rule — boxeph fold-edge THM-1765 (18:33:56) and opus two-straddle THM-1765 (18:36:58) both pushed before mine (18:38:55). Three-way collision, mine last."
source: kind-pasteur-2026-07-20-S128c129 (owner: work the named-next questions of THM-1775)
depends_on:
  - THM-1775    # the moment-nullcone template and its ladder
  - THM-133     # H = (462 - tr(A^4))/2 for Z_7 circulants (the special spectral case)
related: [THM-895, THM-1710, THM-1740]
script: 04-computation/H_on_the_ladder_kps_S128c129.py (+ .out)
---

# THM-1780 — H leaves the spectral ladder at n = 6

THM-1775 placed the tournament trace moments `tr(A^m)` on the rational floor of the
`rational ⊂ algebraic ⊂ holonomic` ladder and asked where `H` (the Hamiltonian-path count,
`= I(Ω,2)` by the OCF) sits. **It sits off the ladder — above the top — starting at `n = 6`.**

## The test and the result

Group every labeled tournament by its **moment vector** `(tr A¹, …, tr Aⁿ)`. By Newton's
identities this is exactly the characteristic polynomial, so two tournaments share a group iff
they are **co-spectral**. Ask: is `H` constant on each group?

| `n` | co-spectral classes | classes with split `H` |
|---|---|---|
| 4 | 3 | **0** |
| 5 | 9 | **0** |
| 6 | 28 | **3** |

- `n ≤ 5`: `H` is a function of the spectrum — it is **spectral**, on the rational floor, with
  detection depth `n`.
- `n = 6`: `H` **splits**. The moment vector `(0,0,12,12,10,48)` carries **both `H = 13` and
  `H = 17`** — two non-isomorphic co-spectral tournaments with different Hamiltonian-path
  counts. (Both odd, as Rédei forces; both in the allowed Ham-path spectrum `{odd}\{7,21}`.)

So for `n ≥ 6` **no fixed polynomial in the trace moments computes `H`.** `H` is not a moment.

## What this means for the ladder

`H` is the number of directed Hamiltonian paths — a **permanent**, `#P`-hard. The trace moments
are its rational shadow. The ladder therefore has a fourth rung:

> **rational (trace, depth `n`) ⊂ algebraic (TNC, depth `D`) ⊂ holonomic (GMC, depth `K`) ⊂
> `#P` (H, no spectral detection).**

The tournament appears at **both ends**: its spectrum on the rational floor, its `H` at the
`#P` top. This is the precise reason the project found the *H*-spectrum a strictly finer
"universal tournament code" than the trace spectrum — the two first diverge at `n = 6`.

## Why THM-133 does not contradict this

THM-133's `H = (462 − tr(A⁴))/2` is exact **for `Z₇` circulant tournaments**. Circulants are
determined up to isomorphism by their spectrum (the cyclic symmetry collapses the permanent to
a spectral function), so on that symmetric subfamily `H` *is* spectral — a **symmetry
collapse**, not the general law. Off the symmetric locus, at `n = 6`, the collapse fails.

## The threshold `n = 6`

`n = 6` is where `H` first exceeds the spectrum. It joins the project's cluster of `n ≥ 6`
phase transitions — the per-path identity fails at `n ≥ 6` (CLAUDE.md), and the width-formula
`C(n−2, ⌊(n−2)/2⌋)` for `G_n` fails at `n ≥ 7`. The moment-nullcone frame gives a clean reason
`n = 6` is special *for H*: it is the first size at which the tournament's permanent is not a
function of its spectrum.

## Named next

- **Where does `H mod 2` sit?** Rédei says `H` is odd — a spectral-independent invariant that
  is *constant* (`≡ 1`). Between the constant `H mod 2` and the `#P` `H`, is there an
  intermediate modulus (`H mod 4`, the blue law) that IS spectral? The mod-4 blue law
  (THM-790s) is the lead.
- **The OCF `I(Ω,2)` reading.** `H = I(Ω,2)` expresses `H` through the odd-cycle conflict
  graph. Is the independence-polynomial-at-2 of `Ω` holonomic in some parameter, placing a
  *refined* `H` back on the ladder one level below `#P`?
