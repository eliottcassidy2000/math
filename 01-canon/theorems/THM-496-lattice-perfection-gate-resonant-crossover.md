---
id: THM-496
name: The lattice-perfection gate — the two-factor resonant-product crossover is at
      n=28 because 9 is the FIRST lattice-imperfect size (Harborth(9)=16 < u(9)=18);
      27=3³ is doubly obstructed (chord-free 3-factor ∧ imperfect 9-factor) and the
      RESONANT cap at 27 is 75, not the 81 generic tie; 28=4·7 is lattice-perfect AND
      chord-bearing
status: PROVED (bounds) + VERIFIED exact-integer (full connected-patch enumeration
        k≤9; exhaustive 2-factor resonant maxima n=24..27)
date: 2026-06-13
session: monad-explorer-2026-06-13 (deep-research; OPEN-Q-057 frontier)
depends_on:
  - THM-493   # resonant-product decomposition U = e(G)|H|+|G|e(H)+Δ_t
  - THM-495   # (concurrent) chord-spectrum bottleneck: Δ_t>0 ⟹ t∈ChordSpec(small factor)
  - THM-432   # generic-angle product cap; products tie 3N at 27,30
  - THM-437   # cube K₃^□3 angle-rigid at 81
  - THM-431   # u(n) exact table (AMP); N*∈[25,28]
external:
  - "Harborth (1974): max unit edges of an n-point triangular-lattice patch
     = ⌊3n−√(12n−3)⌋ — here RE-VERIFIED exhaustively for n≤9 (all connected patches)."
  - "Alexeev–Mixon–Parshall arXiv:2412.11914: exact u(n), n≤21 (u(9)=18, u(10)=20)."
complements:
  - "THM-495 (same session, parallel monad-explorer): the chord-spectrum gate answers
     WHICH t crosses. THM-496 answers WHY n=28 not n=27, on the orthogonal
     lattice-perfection axis, and CORRECTS the resonant cap at 27 (75, not the 81 tie)."
---

# THM-496: the lattice-perfection gate for the resonant-product crossover

THM-495 (the chord-spectrum gate, concurrent) proves the resonance norm `t` is read
off the small factor's chord spectrum, and that `27 = 3³` gets zero bonus because the
triangle `K₃` is chord-free. This theorem adds the **second, orthogonal gate** —
*lattice-perfection of the factor sizes* — which (i) explains why the crossover is at
`n = 28` and not earlier, and (ii) **corrects** the resonant cap at `n = 27`: a
*resonant* product there reaches only **75**, not the `81` generic tie. The two gates
together give the complete account of the `27→28` jump.

## Definitions

For a finite triangular-lattice patch `X ⊂ ℤ[ζ₆]` (`N(x+yζ₆)=x²+xy+y²`), `e(X)` =
number of unit edges. Let
- **`H(k)` = Harborth(k)** = the max of `e(X)` over all `k`-point patches `X ⊂ ℤ[ζ₆]`
  = `⌊3k − √(12k−3)⌋` (Harborth 1974);
- **`u(k)`** = the planar unit-distance max (over *all* `k`-point sets, AMP exact `k≤21`).
- Size `k` is **lattice-perfect** iff `H(k) = u(k)` (the triangular lattice attains the
  global optimum on `k` points).

## Statement

**(A) [VERIFIED exhaustively] Lattice-perfection holds for `k ≤ 8`; `k = 9` is the
FIRST lattice-imperfect size.**
```
   k       2   3   4   5   6   7   8    9   10  11  12  13  14
   H(k)    1   3   5   7   9  12  14   16   19  21  24  26  29
   u(k)    1   3   5   7   9  12  14   18   20  23  27  30  33
   perfect ✓   ✓   ✓   ✓   ✓   ✓   ✓    ✗    ✗   ✗   ✗   ✗   ✗
```
`H(k)=u(k)` for `2≤k≤8`; at `k=9`, `u(9)=18 > 16=H(9)` — a deficit of **2**. (`H(k)`
for `k≤9` is RE-VERIFIED here by enumerating **all** connected triangular-lattice
patches up to size 9 — 3, 11, 44, 186, 814, 3652, 16689, 77359 of them — and taking
the max edge count; matches the formula exactly. No reliance on the 19-hex heuristic.)

**(B) [PROVED] The resonant cap at `n = 27` is `75`, not the `81` tie.** The only
nontrivial two-factor split is `27 = 3·9`. For a *resonant* product `G ⊞_t H`
(`|G|=3, |H|=9`, both `⊂ ℤ[ζ₆]`),
```
   U = e(G)·9 + 3·e(H) + Δ_t,    e(H) ≤ H(9) = 16,    Δ_t ≤ 8·(3 − e(G)),
   ⟹  U ≤ 9·e(G) + 48 + 8(3−e(G)) = e(G) + 72 ≤ 75.
```
The bonus bound: `G` has `3` pairs, `e(G)` of them unit (norm 1 ≠ t), so at most
`3−e(G)` have norm `t`; hence `Σ_{N(α)=t} m_α(G) ≤ 2(3−e(G))`, and `m_α(H) ≤ |H|−1=8`
(displacement-`α` pairs form ≤ `|H|−1` chains), giving `Δ_t ≤ ½·8·2(3−e(G))`.
*Exact maximum (full enumeration) = 75*, attained by `K₃ ⊞ (3×3 patch)` with
`Δ_t = 0` (`K₃` chord-free) — so **resonance strictly HURTS at 27** (`75 < 81`). The
`81`-tie is the **generic / off-lattice** product `K₃^□3` (or `K₃□G₉` with the
*generic* `u(9)=18` factor), which carries **no** resonance. (With only the planar
bound `e(H)≤u(9)=18`, the same computation gives `U ≤ e(G)+78 ≤ 81`: no resonant
product beats `81` at 27 in any case.)

**(C) [PROVED/VERIFIED] `n = 28` is the FIRST two-factor resonant crossing, gated by
lattice-perfection ∧ chord-bearing.** A resonant product beats `3n` iff it has a
factorization meeting all three:
  (i)  **lattice-perfect** sizes (all factors `≤ 8`) — so the Cartesian part equals the
       generic cap `e(G)|H|+|G|e(H)`, with no Harborth penalty;
  (ii) a **chord-bearing** factor (size `≥ 4`; sizes 2,3 are chord-free, THM-495) — so
       `Δ_t > 0`; and
  (iii) `Δ_t > gap(n) := 3n − P_gen(n)`, where `P_gen(n)=max_{ab=n} u(a)b+a·u(b)`.
Exact data (`gap` and exhaustive resonant max `U*`):
```
   n    P_gen  gap   LP-chord factrzn   exhaustive resonant U*   3n   verdict
   24    66     6    4×6, 3×8           68                       72   no  (Δ<gap)
   25    70     5    5×5                72                       75   no  (Δ<gap)
   26    73     5    —  (2·13: 13 imperfect)   ≤65               78   no  (no LP factrzn)
   27    81     0    —  (3·9: 9 imperfect, 3 chord-free)  75      81   no  (B): resonance hurts
   28    83     1    4×7                85                       84   CROSS (Δ₃=2>1)
```
`n=24,25` satisfy (i)+(ii) but fail (iii) (gap 6,5 ≫ realizable bonus 2); `n=26,27`
fail (i)+(ii) (no lattice-perfect chord-bearing factorization — `13,9` are imperfect);
`n=28=4·7` is the FIRST to satisfy all three. The explicit crosser is THM-493's
`W₇ ⊞₃ R`: both sizes `4,7` lattice-perfect (`H(4)=u(4)=5`, `H(7)=u(7)=12`), the
rhombus `R` carries a `√3` chord, `gap(28)=1 < Δ₃=2`, so `U = 83 + 2 = 85 > 84`.

## Why 9 is imperfect — and why that propagates to 27 (the deep link)

`u(9)=18` is realized by `K₃ □ K₃` (the Erdős product of two unit triangles) **only at
a generic relative angle**; forcing the two triangle directions onto one triangular
lattice (`60°`) collapses the graph to the `3×3` rhombic patch with just `16` edges.
So `u(9) > H(9)` is itself a *"the product needs an irrational angle"* phenomenon — the
**smallest instance of the integrality/irreducibility premium** (cf. THM-433's
`N*`-is-non-product, the `χ > χ_f` Vitali wall). And the conjectured `27`-optimum is the
generic cube `K₃^□3 = K₃ □ (K₃ □ K₃) = K₃ □ G₉`: it is `K₃` times the *generic*
`9`-optimum. So the **same off-lattice product structure** that makes `9` imperfect
makes `27` tie at `81`, and the lattice (resonant) route — which could add a `√t` bonus
— cannot replicate it, capping at `75`. Lattice-imperfection at `9` *propagates
multiplicatively* to `27`.

## Relation to canon

- **Complements THM-495 (no conflict, no court case).** THM-495: *which `t`* (small
  factor's chord spectrum, `t=3` forced at 28). THM-496: *why `n=28` not `27`* (the
  lattice-perfection axis) + the **corrected resonant cap 75 at 27**. THM-495 (C) says
  "the best two-factor product can only tie 81 at 27"; precisely, the *generic* product
  ties 81 while the *resonant* product caps at 75 — resonance hurts at 27. Together the
  gate is: lattice-perfect factorization ∧ chord-bearing factor ∧ bonus > generic gap.
- **Sharpens THM-493:** its curated-search caps `68,72,…,75` at `n=24,25,…,27` are now
  EXACT MAXIMA over the full two-factor connected-patch family (not search lower bounds).
- **Re-derives & extends THM-437/THM-432** for the product family from the
  lattice-perfection side.
- **Scope / honesty:** this is the **two-factor triangular-lattice product family**
  (a lower-bound construction lens), not arbitrary configs. It does NOT prove
  `u(27)=81` or `N*=28` (AMP's upper bound at 27 is still 90). It proves: *no two-factor
  product (generic ≤81, resonant ≤75) beats 81 at 27*, and pins the `27→28` crossover to
  the first lattice-imperfect size. Strong structural support for `N*=28` (HYP-2299).

**Artifacts:** `04-computation/lattice_perfection_gate_monad.py`,
`05-knowledge/results/lattice_perfection_gate_monad.out`. New hypothesis **HYP-2467**
(the lattice-imperfection-propagation conjecture). Reflection
`07-reflections/the-lattice-perfection-gate-nine-is-the-first-imperfect-size.md`.
