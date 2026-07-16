---
id: THM-913
title: THE PARALLEL-CLASS BOOK DRAWING — for odd n, the 2-page book drawing of K_n on a cyclic spine with pages = a CONTIGUOUS SPLIT of the parallel-class circle achieves Guy's Z(n) exactly (verified n = 5,7,9,11,13: min over all 2^{n−1} class-colorings = Z(n) = A000241(n), attained at classes {1..(n−1)/2} vs the rest) — hence, by the Ábrego–Aichholzer–Fernández-Merchant–Ramos–Salazar theorem (2-page crossing number = Z(n)), the construction is OPTIMAL. Structure: (i) parallel chords (a+b ≡ c+d mod n) NEVER cross on the cyclic spine (within-class crossings = 0, all n); (ii) the class-crossing matrix is CIRCULANT (spine rotation shifts sums by 2, a class n-cycle since gcd(2,n) = 1); (iii) Σ_{s<t} X[s][t] = C(n,4) (every 4-subset has exactly one interleaved pairing, never same-class) ⟹ 2-page optimality ⟺ MAX-CUT of the class circulant = C(n,4) − Z(n)
status: **PROVED for all odd n** (upgraded death-star-S28): (L1) the class-crossing profile ξ(d) = (d−1)(n−1−d)/2 — refereed exactly all odd n ≤ 31, all d (separation count: exactly-one-endpoint-in-the-symmetric-arc bookkeeping; recursion ξ(d+1)−ξ(d) = m−d); (L2) THE MAIN IDENTITY, three lines: F(contiguous split) = Σ_{e=1}^{m−1} e(M−2e)(M−e)/2 (M = 2m−1) = [M²S₁ − 3M·S₂ + 2S₃]/2, and S₂ = m(m−1)M/6 cancels the M²-terms EXACTLY, leaving m²(m−1)²/4 = Z(2m+1) (refereed m ≤ 60 + independent arc bookkeeping m ≤ 24); (L3) optimality over ALL 2-page drawings by AAFRS 2012 (2-page cr = Z(n) lower bound) — coloring-optimality NOT needed for the theorem. REFINEMENT (the class-coloring max-cut): ŵ(k) = (n−1)/2 − n/(4sin²(πk/n)) IDENTIFIED (exact at n = 9, 13), so by Parseval the size-terms drop and coloring-optimality ⟺ the contiguous arc maximizes the cycle GREEN'S FORM G(A) = Σ_{k≠0}|χ̂_A(k)|²/(4sin²(πk/n)) ⟺ the arc minimizes Σ_{pairs}d(n−d) — verified exhaustively all odd n ≤ 19 (2^{18} colorings at n = 19); the general-n arc-Green maximality is the one named remaining lemma (affects only coloring-universality, not the theorem)
source: death-star-2026-07-16-S27 (owner: work the cyclic 2-page book drawing, spine = Z_9); grew from THM-906(II)'s parallel-chord law
depends_on: [THM-906(II) (parallel chords = sum classes), opus THM-900 (Guy/square-triangular weave), T-crossing-numbers]
external: Ábrego–Aichholzer–Fernández-Merchant–Ramos–Salazar 2012 (2-page cr(K_n) = Z(n)); Guy's conjecture; A000241
verification: 04-computation/glue_run_z9_book_deathstar_S27.py -> 05-knowledge/results/glue_run_z9_book_deathstar_S27.out
---

# THM-913 — the parallel-class book drawing

Place V(K_n) = Z_n on the spine in cyclic order (n odd). The n(n−1)/2 chords partition into
n PARALLEL CLASSES by midpoint sum s ≡ a+b (mod n), each of (n−1)/2 pairwise-noncrossing
chords (equal-sum chords are nested or disjoint — never interleaved). Assigning whole
classes to pages makes the page problem a 2-coloring of the class circle Z_n with circulant
crossing costs X[t−s]; the total cross-class mass is C(n,4) exactly, so

> **2-page optimality (Z(n)) ⟺ the contiguous split is a MAX-CUT of the class circulant,
> with cut value C(n,4) − Z(n).**

**PROOF (all odd n, S28).** ξ(d) = (d−1)(n−1−d)/2 (L1; refereed n ≤ 31). The contiguous
split's same-page mass is F = Σ_d [(m−d)⁺ + (m+1−d)⁺]ξ(d) = Σ_{e=1}^{m−1} e(M−2e)(M−e)/2
with M = 2m−1; expanding, the M²·S₁ and −3M·S₂ Faulhaber terms cancel exactly
(S₂ = m(m−1)M/6), leaving F = m²(m−1)²/4 = Z(2m+1). With the AAFRS lower bound the
drawing is optimal among ALL 2-page drawings. ∎  (Referee: m ≤ 60 exact; arc bookkeeping
independent check m ≤ 24; exhaustive coloring minimum = Z(n) for all odd n ≤ 19.)
The spectral refinement: ŵ(k) = (n−1)/2 − n/(4sin²(πk/n)), so among fixed-size class sets
the coloring problem is the cycle Green's-form maximization — the arc-maximality lemma
(⟺ the arc minimizes Σ_{pairs} d(n−d)) is the one open refinement, verified to n = 19.

LRC face: sum classes mod n ARE the repo's kernel-quadruple/parallel-chord structure
(THM-906(II)); the optimal book drawing is the round-robin 1-factorization split into two
contiguous half-rounds — Guy's conjectured optimum realized by exactly the additive
structure that makes near-APs extremal in the LRC covering program. The two extremal
worlds (crossing minimization, additive-energy maximization) share one object: the
parallel-class circle.
