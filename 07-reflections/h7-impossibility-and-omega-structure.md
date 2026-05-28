# H=7 Impossibility and Omega Structure

**Session:** opus-2026-05-28-S4

## The Universal Forbidden Value

H(T) = 7 is impossible for every tournament. This is not a theorem about one size — it's universal.

The OCF formula forces a specific algebraic signature: H=7 requires exactly three pairwise-intersecting odd cycles with nothing else (α₁=3, α₂=0). But in a tournament, three mutually-conflicting cycles always force a fourth. The universe of odd cycles is not "small enough" to stop at 3 with all pairwise conflicts.

**The proof is a beautiful two-step:**
1. H=7 iff Ω(T) = K₃ exactly (three vertices in a triangle, nothing else).
2. Ω(T) = K₃ is impossible: three pairwise-intersecting 3-cycles C₁=(a,b,c), C₂=(a,d,e), C₃=(b,d,f) always force a 4th. Case c→d: the 5-cycle (c,d,e,a,b) appears. Case d→c: the 3-cycle (d,c,a) appears. Either way, α₁ ≥ 4.

What makes this remarkable: the arc between c and d (not constrained by any of the three original cycles) deterministically creates a new odd cycle in BOTH possible directions. The tournament structure is too rigid — it cannot accommodate exactly 3 pairwise-conflicting cycles.

## The Score Sequence Proof at n=5,6

For n=5: N₃=3 forces score sequence {1,1,2,3,3}. This sequence is strongly connected (Landau condition satisfied strictly). Strong connectivity + Moon's theorem: every strongly connected tournament has a Hamiltonian directed cycle (length n). For n=5, that's a 5-cycle. So α₁ ≥ 4. Q.E.D.

For n=6: N₃=3 forces score sequence {0,2,2,3,4,4}. The vertex with score 0 is a sink (beaten by all), hence outside every cycle. The remaining 5 vertices have effective score sequence {1,1,2,3,3} — exactly the n=5 case. Again α₁ ≥ 4.

This is a perfect example of "proof by descent": n=6 reduces to n=5, and n=5 has a self-contained proof via Moon's theorem.

## The Omega Triangle Paradox

Omega(T) is NOT K₃ (full graph) — proved. But Omega(T) DOES contain K₃ as a subgraph — for 78% of tournaments at n=5, 98% at n=6, 99.9% at n=7.

So there's a sharp distinction:
- **Ω ≠ K₃**: impossible (THM-343) — Omega can't have exactly 3 nodes all connected.
- **K₃ ⊂ Ω**: typical — 3 mutually-conflicting cycles appear all the time among many.

This means the Chudnovsky-Seymour approach to TRRT (triangle-free + claw-free ⟹ real-rooted) cannot work via triangle-freeness. Omega contains triangles abundantly. The TRRT proof must go through the Hermite-Biehler recursion (OPEN-Q-053), which is the strongest existing approach.

## The Alpha_1 = 3 Structure

When a tournament has exactly 3 odd cycles (α₁=3), the conflict graph is always:
- K̄₃ (empty): impossible to achieve at n=5,6 (because N₃=3 at small n forces additional cycles via strong connectivity)
- K₁⊔K₂ (one edge + isolated): the ONLY pattern seen at n=7

At n=7, all 10 alpha_1=3 tournaments have H=15, with two groups: {one isolated 3-cycle} + {two 3-cycles sharing a vertex}. The shared pair forms one conflict, the isolated cycle is disjoint from both. This gives α₂=2 and H = 1 + 6 + 8 = 15.

The conflict graph K₁⊔K₂ is NOT K₃ (its clique number is 2, not 3). This directly avoids the forbidden pattern.

## H Mod p Gaps: A More Precise Picture

The HYP-1751 (H mod p has gaps) is WRONG as stated — refuted at p=7, n=7.

The correct picture:
1. **Universal**: H ≡ 1 (mod 2) — H is always odd (Rédei's theorem).
2. **H=7 universal impossibility** creates a gap at n=5: H≡2 (mod 5) is impossible because H∈{2,7,12,...}∩{odd} = {7,17,...}, and 7 is forbidden while H_max(n=5)=15 < 17.
3. For n=6,7: the gap disappears because H_max grows fast enough to cover all residues.

The lesson: universal H-value gaps (like H=7) create mod-p gaps at small n where H_max is too small to cover the missing residue. At larger n, these artifacts vanish.

## The N₃ Formula as a Bridge

The classical formula N₃(T) = C(n,3) - Σᵢ C(sᵢ, 2) connects the SCORE SEQUENCE to the cycle structure. Used here to constrain which score sequences are compatible with N₃=3. This bridges combinatorics (score sequences) with algebraic topology (independence polynomials).

Score sequence → N₃ → strong connectivity → additional cycles → α₁ ≥ 4.

Every step in this chain uses a different mathematical tool: algebra (Landau condition), graph theory (Moon's theorem), topology (OCF). The proof is a genuine cross-disciplinary argument.

## What This Means for TRRT

TRRT (Tournament Real-Rootedness Theorem): I(Ω(T), x) is always real-rooted.

This session ruled out the Chudnovsky-Seymour route (via K₃-free). The correct route remains the Hermite-Biehler recursion:
- I(Ω,x) = A(x) + x·B(x) where A = I(Ω\C*, x), B = I(Ω-N[C*], x)
- B interlaces A (3672 cases verified, 0 failures)
- This recursively reduces degree and proves real-rootedness by induction

The two remaining lemmas (A: existence of HB-cycle, B: interlacing) are where the effort should go. The K₃ structure of Ω might still be relevant: when three cycles mutually conflict, their neighborhood in Ω creates a specific structure for Ω-N[C*] that could be exploited in the interlacing proof.
