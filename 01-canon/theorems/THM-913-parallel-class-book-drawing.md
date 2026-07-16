---
id: THM-913
title: THE PARALLEL-CLASS BOOK DRAWING — for odd n, the 2-page book drawing of K_n on a cyclic spine with pages = a CONTIGUOUS SPLIT of the parallel-class circle achieves Guy's Z(n) exactly (verified n = 5,7,9,11,13: min over all 2^{n−1} class-colorings = Z(n) = A000241(n), attained at classes {1..(n−1)/2} vs the rest) — hence, by the Ábrego–Aichholzer–Fernández-Merchant–Ramos–Salazar theorem (2-page crossing number = Z(n)), the construction is OPTIMAL. Structure: (i) parallel chords (a+b ≡ c+d mod n) NEVER cross on the cyclic spine (within-class crossings = 0, all n); (ii) the class-crossing matrix is CIRCULANT (spine rotation shifts sums by 2, a class n-cycle since gcd(2,n) = 1); (iii) Σ_{s<t} X[s][t] = C(n,4) (every 4-subset has exactly one interleaved pairing, never same-class) ⟹ 2-page optimality ⟺ MAX-CUT of the class circulant = C(n,4) − Z(n)
status: CONSTRUCTION + EXACT VERIFICATION n = 5..13 (all 2^{n−1} colorings enumerated; min = Z(n) at the contiguous split every time); (i) and (iii) are one-line lemmas (parallel chords have equal midpoint-axis ⟹ nested or disjoint, never interleaved; the C(n,4) count is the classic one-interleaving-per-4-subset); (ii) verified exactly (row profiles e.g. n=9: [0,0,3,5,6,6,5,3,0]); the general-n proof of contiguous-split optimality = a circulant max-cut statement (named)
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

Verified exactly for n = 5..13: the contiguous split {1..(n−1)/2} attains Z(n) =
A000241(n) = 1, 9, 36, 100, 225 — the 2-page crossing number (AAFRS 2012), so the
construction is optimal. The circulant profiles (n = 9: [0,0,3,5,6,6,5,3,0]) are the
crossing spectra between sum classes at circulant distance t — symmetric near-quadratic
ramps; their closed form and the general-n max-cut proof are the named follow-ups.

LRC face: sum classes mod n ARE the repo's kernel-quadruple/parallel-chord structure
(THM-906(II)); the optimal book drawing is the round-robin 1-factorization split into two
contiguous half-rounds — Guy's conjectured optimum realized by exactly the additive
structure that makes near-APs extremal in the LRC covering program. The two extremal
worlds (crossing minimization, additive-energy maximization) share one object: the
parallel-class circle.
