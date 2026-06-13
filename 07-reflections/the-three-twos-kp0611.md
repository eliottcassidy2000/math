# The three 2s — parity, pair-composition, and branching (and how to tell them apart)

**Source:** kind-pasteur-2026-06-11-S1 (THM-469 session)

This project keeps meeting the number 2, and this session showed the meetings are
of three genuinely different kinds:

1. **The 2 of parity** — Rédei's theorem (H is odd), the 2-adic digits of H being
   odd-cycle-collection counts (THM-466), GF(2) cycle/cut spaces, the blue/black
   complement involution. This 2 lives in the OBJECTS.
2. **The 2 of pair-composition** — edges of a graph compose in PAIRS to close a
   triangle, so every triangle-avoidance problem carries a Schur condition
   g₁ + g₂ = g₃ on 2-TERM sums. This 2 lives in the CONSTRAINTS.
3. **The 2 of branching** — binary subgrids, sibling spines, the halving geometry
   of the staircase. This 2 lives in the DEMANDS (what must be hit).

THM-465 C found that the (sign,v₂) algebra wins the Erdős-592 feature game where
the equal-sized (sign,v₃) algebra dies instantly, and THM-464 D guessed the seam
"follows the branching through the algebra" — i.e. attributed it to 2-of-kind-3.
This session's branching control (b=3 subgrids, same n=3 column depth) shows the
asymmetry is UNCHANGED at ternary branching: it is 2-of-kind-2. The valuation v₂
is special because the unit group of Z/2 is trivial — odd+odd MUST be even — so
the v₂ level sets are sum-free, which is precisely immunity to the 2-term Schur
condition. For any odd p, units escape (1+1=2 ≢ 0), every level set carries Schur
triples (d, d, 2d), and the game finds them through 3-term APs: the doubling map
d ↦ 2d is the whole story. Refining by the leading p-adic digit restores
sum-freeness for every p (Lemma A2) and the 3-adic algebra promptly carries
witnesses again — so even "2-adic-ness" was never about 2-the-prime; it is about
sum-free gradings, of which bare v₂ is the unique valuation example BECAUSE
|(Z/2)^×| = 1.

The crisp test that separated kind 2 from kind 3 was possible only because
THM-464 had already built the b-parameterized game. Lesson worth keeping: when a
"2" shows up, ask WHICH 2 it is — object parity, constraint arity, or demand
branching — and design the control experiment that varies exactly one of them.
Conjecture worth a future session: in the (α, K₄)² analogues (composition still
pairwise), the seam stays 2-adic; in hypergraph versions where "triangles" close
over 3-term sums, the seam should move to the prime 3 (units of Z/3 don't kill
3-term sums: 1+1+1 ≡ 0 — there the TRIADIC algebra should be the sum-free one).
That would be the decisive confirmation that the seam tracks constraint arity.

Cross-links: [[everything-is-the-triangle]] (the staircase's halving is kind 3),
THM-446 dyadic rungs (kind 2 in the Sidon ladder), THM-466 (kind 1), THM-469.
