# The Octal Lens on the H-Spectrum: residue 1 mod 8 is gap-free

**Session:** mac-mini-2026-06-16-S1
**Builds on (does not repeat):** T830 / HYP-2557 (kind-pasteur boat session, same human dispatch, earlier today), THM-485 (two temperatures), THM-486 (Pisano π(24)), THM-224 (golden exceptional points), the Pascal-slope-d Pisot tower, THM-462/463 (H-spectrum gap-freeness).

## The ground that was already mapped

The human re-raised a constellation — Fibonacci, triangular, square, prime, Heronian
triangles, the octal odd-square fact, Fib∩Tri = {1,3,21,55}, the 2×6-marble Burnside,
α(p) | p−(5/p). T830 already records all of it: 21 = Φ₃(1)·Φ₃(2) = 3·7 as the meeting
point of Fib∩Tri, T₆ = C(7,2), the (2n−1)(4n−1) = T_{4n−2} family, and the forbidden
phantom; odd² = 8T_k+1; the Euler-criterion/Paley reading of the Fibonacci entry point;
the Heron 4-product / four-square / 4 = 2² motif; the 5-inradius-2-triangles ↔ 5-Platonic
**count coincidence**. I re-verified every numeric claim independently (third confirmation;
`fib_tri_square_prime_heron_mac-mini-2026-06-16-S1.py`). Nothing there needed redoing.

The spine worth restating once, because it is *structural*, not a pun: the project's
tile-count **m = C(n−1,2) = T_{n−2}** is a triangular number — the dimension of the tiling
hypercube Q_m, the exponent in 2^{T_{n−2}}. The staircase is not only triangle-*shaped*;
its area-in-tiles **is** a triangular number. So the integers' own triangular-number
arithmetic is, a priori, arithmetic about Q_m.

## The new move: the octal identity as a lens, not a fact

T830 lists odd² = 8·T_k + 1 as a fact about integers. The new step is to *aim* it. Two
things are true at once and had never been put together:

1. Odd squares are **exactly** the residue-1 class mod 8 (8·T_k + 1), and the "octal peel"
   floor/8 recovers the triangular number T_k.
2. H(T) = #directed Hamiltonian paths is **always odd** (Rédei), so H mod 8 ∈ {1,3,5,7} —
   the same odd-residue world odd squares inhabit. The "octal peel" floor(H/8) is the
   direct analogue of T_k.

So: *does the forbidden-H frontier see the residue-1 class differently?* This is a real,
cheap question, and it lands on an open project problem (THM-462/463, the H-impossibility
spectrum). The answer at n ≤ 6 (exhaustive, Held-Karp exact;
`H_mod8_octal_probe_mac-mini-2026-06-16-S1.py`):

```
realized H-spectrum (union n≤6): 1 3 5 9 11 13 15 17 19 23 25 27 29 31 33 37 41 43 45
odd integers in [1,45] NEVER realized:  7   21   35   39
   7  ≡ 7 mod 8        21 ≡ 5 mod 8        35 ≡ 3 mod 8       39 ≡ 7 mod 8
residue 1 mod 8 in [1,45]:  1 9 17 25 33 41  — ALL realized (zero gaps)
realized perfect squares:   1 9 25          — all residue 1, all present
```

**The finding:** every gap in the H-spectrum (at n ≤ 6) avoids residue 1 mod 8. The
odd-square residue class is gap-free. This is exactly consistent with — and sharpens —
T830's Φ₃ observation: the forbidden phantoms sit Φ₃(2) = 7 ≡ 7 and 21 = Φ₃(1)Φ₃(2) ≡ 5,
both **off** the odd-square residue. The octal identity isn't decoration; it predicts where
the forbidden values are *allowed* to be.

## What is structural, what is suggestive

- **STRUCTURAL.** m = T_{n−2}; H odd (Rédei) ⟹ H mod 8 ∈ {1,3,5,7}; the n ≤ 6 gap set
  {7,21,35,39} and its residue distribution (independently reproduces the documented
  universal-forbidden 7 and 21). The 8 = 2³ is the parity cube / Bott period already
  flagged (T232).
- **HONEST CAVEAT.** Four gaps over n ≤ 6 is a small sample. 7 and 21 are documented as
  universally forbidden; 35 = 5·7 and 39 = 3·13 may well be realized at n ≥ 7 (they are
  "not yet realized," not proven forbidden). The conjecture "residue 1 mod 8 is gap-free in
  the H-spectrum" (HYP-2558) must be tested at n = 7, 8 before it earns weight. If it
  breaks, that is itself the answer: the octal lens sees only the smallest phantoms.
- **SUGGESTIVE (flagged in T830, not load-bearing).** 5 Heron triangles ↔ 5 Platonic
  solids (count coincidence); m Fibonacci at n ∈ {3,4,8,12}; the Klein-four Burnside hit.

## The reflection proper

The recurring lesson of this project is that a pattern true at n = 3,4,5 breaks at 6 or 7;
T830 itself is a catalogue of formulas that hit {1,3,21,55} and then diverge. So the
*right* posture toward the human's constellation is not "these all connect" but "which of
these is a lens that survives contact with the tournament data?" The octal identity does —
weakly, at n ≤ 6, but in the one place it can be made to *predict* (the gap residues) it is
not contradicted. The triangle keeps being the board; the surprise is that its area-count
is triangular, so the integers' parity-cube arithmetic (odd = 1 + 2³·triangular) reaches
into the Hamiltonian-path spectrum's residues. Whether that reach extends past n = 6 is the
first thing the next session should check — it costs one Held-Karp sweep at n = 7.
