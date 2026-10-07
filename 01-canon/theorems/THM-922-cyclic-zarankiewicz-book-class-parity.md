---
id: THM-922
title: THE CYCLIC ZARANKIEWICZ BOOK + THE CLASS-PARITY INCLUSION-EXCLUSION — (I) on the Z_{2m} spine with parts = parities, K_{m,m}'s edges are exactly the ODD-sum parallel classes (m classes of m pairwise-noncrossing chords); the contiguous class-split achieves Zarankiewicz's Z(m,m) = ⌊m/2⌋²⌊(m−1)/2⌋² EXACTLY (verified by full coloring enumeration m = 3..7: min = 1, 4, 16, 36, 81) — ordinary-optimal through m = 7 (m ≤ 6 by Kleitman; m = 7 by Woodall); (II) the parity split of K_{2m}'s class circle is the GEOMETRIC INCLUSION-EXCLUSION between Guy and Zarankiewicz: odd-sum classes = K_{m,m}, even-sum classes = K_m ⊔ K_m (the two parts' internal edges) — the classic counting bridge cr(K_{2m}) vs cr(K_{m,m}) + 2cr(K_m) realized as a decomposition of ONE drawing; (III) EVEN-n completeness: the full K_n class-coloring minimum equals Guy's Z(n) at n = 6..14 as well (unequal class sizes m−1/m, within-class crossings still zero) — THM-913's construction is optimal for ALL n (the even-n closed-form identity = the named algebra follow-up); (IV) the ledger: on odd n, each 4-subset has exactly ONE crossing pairing and at most one parallel pairing and THEY ARE NEVER EQUAL (overlap 0, verified) — the crossing/additive dichotomy is exact
status: (I)/(III) EXACT ENUMERATION (all class colorings; m ≤ 7 resp. n ≤ 14); ordinary optimality is unconditional for m ≤ 6 by Kleitman and for m = 7 by Woodall (1993); (II) the decomposition is definitional once seen (sum parity ⟺ endpoint parities differ/agree) — the identities C(2m,4)-split verified; (IV) proved (equal-sum ⟹ nested-or-disjoint) + verified; general-m proof of (I) via the odd-class ξ-profile = the named follow-up (same Faulhaber shape as THM-913's)
source: death-star-2026-07-16-S29 (owner: crossing numbers complete vs bipartite, odd vs even, inclusion vs exclusion — the recursive thread exploration)
depends_on: [THM-913 + LEM-030 (the odd-n machinery), THM-906(II), mac-mini THM-920 (crossing-energy frame), opus HYP-7100 stub (the parallel-class grammar)]
external: Zarankiewicz's conjecture (min(m,n) ≤ 6 by Kleitman 1970; K_{7,7} by Woodall 1993; general case open); Guy's conjecture; AAFRS 2012
verification: 04-computation/arc_green_polarization_deathstar_S29.py -> 05-knowledge/results/arc_green_polarization_deathstar_S29.out
---

# THM-922 — the cyclic Zarankiewicz book and the parity inclusion-exclusion

Put Z_{2m} on the spine; color vertices by parity. An edge's sum-class is odd iff its
endpoints lie in DIFFERENT parts — so the bipartite K_{m,m} is exactly the union of the m
odd parallel classes (each m pairwise-noncrossing chords), and the two parts' internal
K_m's are the even classes. One spine, one class circle, and the parity split performs the
inclusion-exclusion between the complete and bipartite crossing worlds geometrically.

Enumerating ALL 2-page class-colorings of the odd classes: the minimum equals
Z(m,m) = ⌊m/2⌋²⌊(m−1)/2⌋² at every m = 3..7 (1, 4, 16, 36, 81), attained at contiguous
splits — the Zarankiewicz drawing's crossing count, from the same three-line-species
machinery as THM-913. For m ≤ 6 this is optimal by Kleitman; at 7×7 it is optimal by
Woodall's 1993 cyclic-order computation. The even-n complete case also closes
computationally (min = Z(n), n ≤ 14) despite unequal class sizes.

The dichotomy underneath (IV): on odd n a 4-subset's three pairings contain exactly one
crossing and at most one parallel pairing, never the same one — additive coincidence and
geometric crossing are mutually exclusive, which is why sum classes are page-friendly
bundles in every variant.

## Update 2026-10-06 (mac-mini-2026-10-06-oaimath2; `05-knowledge/results/oai2_openai_math_second_reading_20261006.md`, section 1)

**(III) for every `n` (PROVED).** The DDS drawing (THM-913, a classical construction) has exactly `Z(n)` crossings for every `n`, even included. This is classical, and checked here for `n <= 120`. With Ábrego–Aichholzer–Fernández-Merchant–Ramos–Salazar (2012: `ν_2(K_n) = Z(n)`), the class-coloring minimum is `Z(n)` for every `n`, unconditionally. The enumeration above had covered `n <= 14`.

**(I) for every `m`.**
* *The upper bound is PROVED.* The contiguous split of the odd sum classes on `Z_(2m)` has exactly `d_m^2 = Z(m,m)` crossings for every `m`, by a telescoping count:

      Σ_(κ=1)^K (2κ(m−κ) − m)(m − 2κ) = K^2 (m−K−1)^2,   K = floor((m−1)/2).

  It is the DDS drawing of `K_(2m)` restricted to its even–odd edges.
* *The lower bound is CONDITIONAL on openai/math #165* (`cr(K_(m,n)) = Z(m,n)`). With it, the class-coloring minimum is `Z(m,m)` for every `m`, and the 2-page Zarankiewicz conjecture `ν_2(K_(m,n)) = Z(m,n)` of de Klerk–Pasechnik–Salazar (2014) holds for all `m, n`.

Script: `04-computation/experiments/oai2_20261006_crossing_books.py` (+ `.out`).
