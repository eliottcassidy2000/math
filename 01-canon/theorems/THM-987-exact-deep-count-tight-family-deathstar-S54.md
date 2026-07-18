# THM-987 — The exact deep count and the tight family through the funnel (death-star-2026-07-17-S54)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCDeepCountExact.lean`,
standard trio ×3; `decide` only on the q = 14 census gates). Source: HYP-7295.
The capstones on THM-982.

## Statement

1. `deep_count_exact`: for EVERY q ≥ 2, the canonical family's deep count in
   closed form:
   **#deep(q) = 2B + (B + 1 − (q + B) % 2)**, B = ⌊(q−1)/84⌋
   — circle I gives the arcs [1,B] ∪ [q−B, q−1]; circle II the
   parity-corrected central arc ⌈(q−B)/2⌉…⌊(q+B)/2⌋; the circles are
   provably disjoint for all q. All constants literal ⟹ the whole count is
   omega-native. (Recon: exact on all q ∈ [2,2000).)
2. `canonical_lonely`: **∃t, Lonely 14 canonical** — THE TIGHT FAMILY
   (1,…,13), the equality case of LRC(14), PROVED LONELY THROUGH THE B5
   FUNNEL at the resonant modulus q = 14: one deep multiplier (p = 7, the six
   evens at the half-integer — the mirror-fixed point of circle II) against
   six live multipliers (p ∈ {1,3,5,9,11,13}, the equality witnesses);
   THM-951's `lonely_of_census` fires by kernel decide.

## Why this matters

The S42 census pipeline + the S52–53 two-circle theorem close a loop: the
funnel, built for dissociated strata, ALSO handles the maximally-resonant
tight family — at the resonant modulus the S49 scope obstruction dissolves
(deep is capped at exactly 6 and live points are the equality instants). The
adaptive-q moral in final form: the funnel closes every family at ITS OWN
modulus. Lean craft note: omega atomizes syntactically-distinct spellings of
identical div-terms ((q+B)/2 vs (B+q)/2) — align spellings across statement
and proof.
