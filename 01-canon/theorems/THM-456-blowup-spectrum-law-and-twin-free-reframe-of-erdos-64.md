# THM-456 — The Blowup Spectrum Law and the Twin-Free/Turán-Corridor Reframe of Erdős–Gyárfás (Problem 64)

**Status:** Interval lemma + path law PROVED (constructive) + verified exhaustively (995
connected atlas graphs, 0 failures). Exact spectrum law VERIFIED 995/995 (and the circumference
formula c = 2s verified independently by the adversarial verifier on all 144 nontrivial cases,
0 failures). Census results EXHAUSTIVE n ≤ 12. Erdős–Gyárfás itself remains OPEN — this is a
structural reframe + complete small-case census, not a proof.
**Source:** kind-pasteur-2026-06-09-S2 (branch II + verifier; `verify_II_blowup_interval_kps2.out`).
**Related:** THM-445/446 (E-G statement + Sidon ladder), THM-447/454 (doubling, twin-lift = the
tournament parity mixer), HYP-2358 (resolved into this), MISTAKE-068.

## (1) The interval lemma (PROVED)

If a graph G has a k-cycle (k ≥ 3), the twin blowup G[K₂] has cycles of EVERY length in
[k, 2k]: double t of the k cycle-vertices via twin detours u⁰→u¹, giving length k+t for all
0 ≤ t ≤ k. Since [k, 2k] always contains 2^⌈lg k⌉ ≥ 4, and δ(G[K₂]) = 2δ(G)+1 ≥ 3 whenever
δ(G) ≥ 1 — and in fact a SINGLE EDGE uv already plants the twin 4-cycle u⁰u¹v¹v⁰ —
**blowups can never be Erdős–Gyárfás counterexamples**, robustly (a power of 2 in every dyadic
window, not just one).

## (2) The full spectrum law (path law PROVED; exact law VERIFIED 995/995)

- PATH LAW: a path on r vertices in G gives cycles of lengths 2r and 2r−1 in G[K₂];
  so spec(G[K₂]) ⊇ [3, 2p(G)] (p = longest-path order). Verified 995/995.
- EXACT LAW: spec(G[K₂]) is the GAP-FREE interval [3, c(G[K₂])], and
  **c(G[K₂]) = 2·s(G)** where s(G) = max order of a subgraph carrying edge multiplicities {1,2}
  with all multi-degrees in {2,4} (paths, cycles, thetas, figure-eights, sun/jellyfish walks —
  closed stutter-walks visiting each vertex ≤ 2 times).
- c = 2p for 939/995 graphs; 56 beaters, smallest = the NET graph (C₃ + 3 pendants):
  c = 12 > 2p = 10 via the sun walk c₁p₁c₁c₂p₂c₂c₃p₃c₃. (Tie at n=6: C₄+2 pendants also beats.)
- 896/995 blowups are Hamiltonian.

## (3) The twin-free / Turán-corridor reframe (the structural heart)

An E-G counterexample contains no C₄ (a 4-cycle is the smallest power of 2), i.e. **no two
vertices with two common neighbors** — in particular any twin pair (true or false) in a δ≥3
graph creates a C₄, so counterexamples are automatically twin-free. Quantitatively:
C₄-freeness caps e ≤ ex(n; C₄) while δ ≥ 3 forces e ≥ ⌈3n/2⌉. The **Turán corridor**
ex(n;C₄) − ⌈3n/2⌉ has margins:
```
n:       4   5   6   7   8   9   10   11   12
margin: −2  −2  −2  −2  −1  −1   +1   +1   +3
```
**The corridor is CLOSED for n ≤ 9** (a C₄ is forced by density alone: E-G trivially true)
and opens at n = 10 — E-G counterexamples live within O(1) edges of C₄-Turán-extremality.
(ex(n;C₄) values re-derived internally for n ≤ 12, matching Clapham–Flockhart–Sheehan 1989.)

## (4) Exhaustive census n ≤ 12 (0 counterexamples, with the structural reason)

Connected δ≥3 graphs: 1, 3, 19, 150, 2589 at n=4..8 (all contain C₄ for n ≤ 9).
ALL C₄-free δ≥3 graphs: 5 at n=10 (incl. Petersen; 3 of the 5 have girth 3 — C₄-free ≠ girth 5),
9 at n=11, 57 at n=12. **Every one of the 71 has dyadic profile EXACTLY {8}** — uniformly
killed by a forced C₈, none even needs C₁₆. (First C₁₆-only kills appear at 24 vertices:
Markström's four cubic graphs, one planar.) Consistent with Royle–Markström: counterexamples
need ≥ 17 vertices, cubic ones ≥ 30.

## (5) Local move taxonomy (604,445 instances verified, 0 failures)

- M1 (+1): an off-cycle vertex dominating an EDGE of the cycle (two consecutive cycle
  neighbors) forces length k+1 — this, not a bare triangle, is the correct +1 move.
- M1' (fan): off-cycle vertex with cycle neighbors at cycle-distance d forces d+2 and k−d+2.
- M2 (+2): off-cycle edge xy dominating consecutive cycle vertices. M3: chord split.
- **M5 (sunlet interval, the weakest interval-forcer found): t vertex-disjoint edge-dominators
  on a k-cycle force the FULL interval [k, k+t]** — the blowup is the extreme case t = k.
- Bondy–Vince (JGT 27 (1998) 11–15) verified on all 173 δ≥3 graphs n≤7; exactly 2 graphs
  achieve only difference 2 (K_{3,3} — bipartite spectra all-even): "1 or 2" is sharp.
- Theta law: spec(θ(a,b,c)) = {a+b, a+c, b+c} — the lacunary extreme — but subdivisions have
  degree-2 vertices: precisely why E-G needs δ ≥ 3 and why the hard core is NOT subdivisions
  but near-C₄-extremal graphs.

## (6) Curiosity

The unique cubic 10-vertex graph avoiding C₈ has cycle spectrum {3,4,5} (g6: IsTX@?B?w) —
circumference 5 on 10 vertices.

## Scripts

`blowup_interval_lemma_kps2.py`, `twin_local_moves_kps2.py`, `mindeg3_spectrum_census_kps2.py`,
`verify_II_blowup_interval_kps2.py` (+ .out). Citations web-verified (Bondy–Vince; Wikipedia
E-G; Markström Congressus Numerantium 171 (2004); Clapham–Flockhart–Sheehan JGT 13 (1989);
Carr arXiv:2605.22844 (2026 preprint, unrefereed); Hu–Shen Discrete Math 2024 P10-free).
