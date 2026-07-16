---
id: THM-882
title: THE FLAT = 2×GOOD LAW — the clock theorem's flat set is measure-exactly-twice the deep-well good set, with full mechanism. (I) cl(G) ⊂ F; both are 12 intervals, one per Farey-12 cell (a/i, b/j) with i+j = 13; per-cell |F| = 1/(14·min(i,j)·|i−j|) vs |G| = 1/(14ij), ratio max(i,j)/|i−j|; λ(F) = 2λ(G) = 2H₁₂/91 = 6617/97020. (II) The overhang F∖G has measure EXACTLY λ(G) — flat-but-not-lonely times. (III) The mechanism is NOT a pointwise 2-to-1 map (per-cell ratios are 7, 8/3, 9/5, 10/7, 11/9, 12/11 — non-uniform); it is the even/odd harmonic split of the unit group (Z/13)^* under halving u = 2v, riding THM-819's primitive harmonic law. (IV) GENERAL LAW: for every ODD N, λ(flat of S₂^{(N)}) = 2·m({1..N−1}; 1/(N+1)); FAILS at every even N (no even units), exact ratios recorded. In LRC terms n = N+1: the 2× law is an EVEN-n phenomenon — n = 14 qualifies
status: PROVED (elementary; every step machine-verified in exact ℚ: 12/12 per-cell interval predictions, 12/12 two-gap word structures, 34/34 three-gap αβ-adjacencies, 2× law exact at N = 3..17 odd, failure exact at N = 4..16 even, unit identities to N = 25; Part 7 of the script independently re-derives and re-verifies THM-819's prime criterion through n = 28)
source: boxeph-2026-07-16-S23; resolves the "6617 identity" backlogged by mac-mini-S113 and named by kind-pasteur-S128 cont.28 (their candidate mechanism "±/2-fold S₂-kink symmetry" is replaced by the unit-group halving)
depends_on:
  - LEM-020 + MISTAKE-152   # klein: S₂ ≥ 6/7 with flat bottom = the polytope P (gaps ≤ 1/7, adjacent 2-sums ≥ 1/7)
  - THM-878                 # the clock theorem (flat census 6617/97020; D(q) = 0 ⟺ q ∈ {7,13,14})
  - THM-819                 # the primitive harmonic law m({1..k};1/(k+2)) = (2/((k+1)(k+2)))·Σ_{units} 1/u — the good-set side of the ledger below
  - THM-853(II) / THM-805   # the corridor law / harmonic law giving m({1..12}; 1/14) = H₁₂/91 = 6617/194040
related: [THM-874 (Möbius-log² grammar), THM-873 (Ramanujan expansion), HYP-6885/THM-819 (the good-set closed form; my per-site derivation re-proves it en route), THM-826 (effective arcs)]
script: 04-computation/lrc14_flat_2x_law_boxeph_S23.py -> 05-knowledge/results/lrc14_flat_2x_law_boxeph_S23.out
---

# THM-882 — the flat = 2×good law

**Objects.** For N ≥ 3 let S₂^{(N)}(t) = Σ_{m=1}^{N−1} (N−m)·(κ − ‖mt‖)₊ with κ = 2/(N+1)
— the pair-overlap energy of the orbit {t, 2t, …, Nt} at radius κ, since N−m = #{(u,v) :
1 ≤ u < v ≤ N, v−u = m} and ‖(v−u)t‖ = d(ut, vt). LEM-020's floor argument gives
S₂^{(N)} ≥ Nκ − 1 = (N−1)/(N+1), with flat set F_N = the orbit slice of the polytope
P = {all gaps ≤ κ, all adjacent 2-sums ≥ κ}. Let G_N = {t : ‖mt‖ ≥ κ/2, m ≤ N−1}, the
good set of the (N−1)-core at threshold 1/(N+1). N = 13 is the LRC(14) instance:
κ = 1/7, floor 6/7, F = THM-878's flat set, G = the deep-well good set at 1/14.

**Theorem.**
1. **(Cell census)** F_N and G_N both consist of exactly one interval in each Farey_{N−1}
   cell (a/i, b/j) with i + j = N (one per ordered coprime pair, φ-many... exactly
   #{i ∈ (Z/N)* } cells), and F_N carries zero measure in every cell with i+j > N.
2. **(Per-cell anatomy)** In the cell with left denominator i and right denominator j,
   writing α = it − a and β = b − jt (so jα + iβ = 1 on the cell):
   F∩cell = {α ∈ [κ/2, (1−iκ)/(j−i)]} for i < j (mirror for i > j),
   G∩cell = {α ≥ κ/2, β ≥ κ/2}; hence cl(G) ⊆ F with
   |F∩cell| = κ/(2·min(i,j)·|i−j|), |G∩cell| = κ/(2ij), ratio max(i,j)/|i−j|.
3. **(The 2× law)** For N odd: λ(F_N) = 2·λ(G_N) = (2κ/1)·Σ_{i∈(Z/N)*, i<N/2} 1/(i(N−i)).
   At N = 13: λ(F) = 2H₁₂/91 = 6617/97020, λ(G) = H₁₂/91·(1/2)… = 6617/194040.
   For N even the law FAILS (exact ratios 3/2, 5/4, 35/22, 21/16, 77/46, 1287/952, 4719/2846
   at N = 4..16).
4. **(Overhang)** λ(F_N ∖ G_N) = λ(G_N) for N odd: the flat set = the good set plus an
   equal-measure overhang — the times where the minority orbit gap has dipped below κ/2
   (loneliness fails) while every adjacent 2-sum still clears κ (pair energy stays minimal).
5. **(Good-set closed form — THM-819 re-proved per-site)** m({1..N−1}; 1/(N+1)) =
   (2/(N+1))·Σ_{i∈(Z/N)*, i<N/2} 1/(i(N−i)) — algebraically identical to THM-819's
   primitive harmonic law (2/(N(N+1)))·Σ_{units} 1/u, whose prime criterion ("the mod-6
   pattern was primality's shadow") my Part 7 independently re-verifies through n = 28,
   including the exact value m({1..24}; 1/26) = 15987925/1546777232 ≠ H₂₄/325 at the
   first odd-composite-coprime-to-6 station n = 26. Nothing here is new relative to
   THM-819; it is recorded as the per-site (Farey-cell) route to the same law, which is
   the route the flat side (1)–(4) needs.

**Proof.**

*(a) Successor rule.* Fix a cell (a/i, b/j) of Farey_{N−1} (so bi − aj = 1, i+j ≥ N, and
each ordered coprime pair (i,j) with i+j = N occurs exactly once). For t inside, the orbit
{mt}, m = 1..N, has right-neighbor structure: adding i moves right by α, adding j moves
left by β. Re-indexing m ∈ {0,…,M}, M = N−1: the gap after point m is α if m ∈ [0, M−i],
γ = α+β if m ∈ (M−i, j), β if m ∈ [j, M] — multiplicities N−i, i+j−N, N−j. Two-gap ⟺
i+j = N.

*(b) Three-gap cells carry no flat mass.* The walk order is m → m+i (α-steps), m → m−j
(β-steps), m → m+i−j (γ-steps). Take m = max(0, j−i): then m ≤ M−i (using j ≤ M) and
m+i ≥ j, so an α-gap is immediately followed by a β-gap: **αβ-adjacency exists in every
cell**. In a cell with a γ-gap, flatness forces γ = α+β ≤ κ and (αβ-adjacency) α+β ≥ κ,
so α+β = κ; since α+β = (i−j)t + const with i ≠ j, this is a single t — measure zero. ∎

*(c) Two-gap word.* For i+j = N the walk step is uniformly +i mod N (m−j ≡ m+i), a full
cycle since gcd(i, N) = 1; the gap word is the rotation coding w_r = [ r·i mod N ∈ [0, j−1] ]
(α ⟺ residue < j). Consecutive letters sit at residues x, x+i mod N. ββ-adjacency ⟺
∃ x ∈ [j, N−1] with x+i−N ∈ [j, N−1]∖… ⟺ ∃ m ∈ [2j, N−1] ⟺ j < i; likewise αα ⟺ j > i:
**the minority letter is never self-adjacent; the majority letter always is.** ∎

*(d) Per-cell solve (i < j WLOG).* Flatness on the cell = {α ≤ κ, β ≤ κ (all gaps ≤ κ);
2α ≥ κ (αα occurs); α+β ≥ κ (αβ occurs); no ββ constraint}. On the line jα + iβ = 1:
β ≤ κ ⟺ α ≥ (1−iκ)/j, which is ≤ κ/2 for all i ≥ 1 (equality at i = 1); α+β ≥ κ ⟺
α ≤ (1−iκ)/(j−i) (the slope dβ/dα = −j/i has modulus > 1, so α+β decreases in α); and
(1−iκ)/(j−i) ≤ κ automatically. So F∩cell = {α ∈ [κ/2, (1−iκ)/(j−i)]}, of α-length
[2(1−iκ) − κ(j−i)]·/(2(j−i)) = (2 − κ(i+j))/(2(j−i)) = κ/(2(j−i)) using κ = 2/(N+1),
i+j = N; t-length κ/(2i(j−i)). G∩cell = {α ≥ κ/2, β ≥ κ/2} has t-length
(1/(ij))(1 − (κ/2)(i+j)·(N+1)/…) = κ/(2ij). Both are anchored at the same left endpoint
α = κ/2, and (1−iκ)/(j−i) > the G-endpoint (1−iκ/2·…): cl(G∩cell) ⊂ F∩cell. ∎

*(e) Aggregate.* λ(F) = 2·Σ_{i<N/2, gcd(i,N)=1} κ/(2i(N−2i)) (the mirror t ↦ 1−t pairs
the (i,j) and (j,i) cells). Partial fractions: 1/(i(N−2i)) = (1/N)(1/i) + (2/N)(1/(N−2i)).
As i runs over units < N/2, the map i ↦ N−2i is a bijection onto the ODD units of (Z/N)*
(N odd; inverse u ↦ (N−u)/2), so Σ 1/(N−2i) = Σ_{u odd unit} 1/u. Meanwhile
λ(G) = Σ κ/(2i(N−i)) summed the same way = (κ/2)(1/N)·Σ_{all units u} 1/u ×… and the even
units are exactly u = 2v with v a unit < N/2, giving Σ_{units<N/2} 1/v = 2Σ_{even units} 1/u.
Assembling: Σ 1/(i(N−2i)) = (1/N)[Σ_{units<N/2} 1/i + 2Σ_{odd units} 1/u]
= (1/N)[2Σ_{even} + 2Σ_{odd}] = (2/N)Σ_{all units} 1/u = 2·Σ 1/(i(N−i)). For N even there
are no even units and the identity fails (LHS strictly smaller by the same bookkeeping). ∎

*(f) Statement (5).* λ(G_N) per-cell summation gives the closed form for every N; the
unit inversion 1/(i(N−i)) = (1/N)(1/i + 1/(N−i)) turns it into THM-819's form. When N is
prime every i < N/2 is a unit and Σ 1/(i(N−i)) = H_{N−1}/N, recovering H_{n−2}/C(n,2);
composite N drops the non-unit terms — THM-819's prime criterion. ∎

## Reading

- **Three clocks, three roles.** Within each flat interval: the left endpoint is the
  ‖it‖ = 1/14 boundary (the 14-clock: 1/14 and 13/14 are the extreme endpoints; the other
  four primitive 14ths are ISOLATED flat points inside i+j = 14 cells — measure zero, so
  the a.e. two-gap law is untouched); each flat interval contains exactly one primitive
  13th in its interior (the 13-clock); and the two central cells (6,7)/(7,6) have flat
  intervals reaching the Farey endpoints 1/7-class fractions (the 7-clock). THM-878's
  {7, 13, 14} are the three stations of every flat interval.
- **The overhang is the adversary's exact headroom.** F∖G (flat-but-not-lonely, measure
  = λ(G)) is where a covering family can sit at minimal pair-energy without conceding a
  1/14-lonely time: the covering adversary's playground has measure exactly twice the
  reward. Any pair-energy (second-moment) instrument is blind to the distinction — the
  quantitative face of LEM-020's "minimal certifying moment order = 3".
- **The mechanism is mod-N halving, not a kink symmetry.** The per-cell ratios
  max(i,j)/|i−j| are non-uniform (7, 8/3, 9/5, 10/7, 11/9, 12/11), so no pointwise
  measure-doubling map exists; the doubling is the invertibility of 2 in (Z/N)*. Note
  −2 ≡ 11 is a primitive root mod 13: the site map i ↦ N−2i is a primitive-root walk,
  and (i,j) ↦ (i, j−i) is the Stern–Brocot parent step — the overhang trades each Farey
  pair for its mediant-parent pair.
- **Even-n exclusivity.** For LRC(n) the 2× law needs n even. The parity of n enters
  through the unit group of Z/(n−1) — a different clock from the mod-14 tight-locus μ₆
  (opus-S320): the flat geometry runs on the N = 13 clock, the tight locus on the n = 14
  clock; the pair (13, 14) = (walk modulus, threshold modulus) is the coprime gear pair
  183 = 13·14 + 1 rides on.

## Evidence log
- [x] N = 13: 12/12 cell predictions exact; word structure 12/12 + 34/34; measures exact
- [x] N = 3..17: 2× law ⟺ N odd, exact ratios both ways
- [x] unit identities N ≤ 25; harmonic-law criterion n ≤ 28 incl. the n = 26 counterexample
- [ ] Lean: the per-cell solve is decide-shaped (rational linear algebra); the unit
      bijection is Finset bookkeeping — a candidate for the certificate batch
