# HYP-2329 — The Pfaffian alternating-sign angle: the tiling triangle = 3-set inclusion–exclusion = Euler characteristic; max-H (A038375) is irregular but its signs are universal

**Session:** S651
**Status:** CONFIRMED (the dictionary + the inclusion–exclusion formalized; max-H computed; the negative
result on a closed recurrence is honest)
**Provenance forward:** math-lean `Math/Combinatorics/InclusionExclusionTriangle.lean` (sorry-free)
**Prompt:** use the Pfaffian angle to translate between topology/geometry/graphs/tournaments/algebras;
find a more efficient recursive understanding of max H-paths (A038375) via a 7-tournament tiling
decomposition `+A+B+C − D−E−F + G` (sizes `n−1, n−2, n−3`).

---

## 0. What the user's decomposition actually is

The 7-tournament tiling — `A,B,C` (size `n−1`, the corners), `D,E,F` (size `n−2`, the edge overlaps,
*subtracted*), `G` (size `n−3`, the interior triple, *added*) — with the regional count `+A+B+C − D−E−F +
G`, is **exactly 3-set inclusion–exclusion realized on the triangle**: corners = the sets, edges = the
pairwise intersections, interior = the triple intersection. The signs `+,−,+` are the
**Euler-characteristic / chain-complex signs** over the cell dimensions `0` (corners), `1` (edges), `2`
(interior): `χ = 3 − 3 + 1 = 1` (a disk). This is the **same alternating-sign structure as the Pfaffian**
(`det = Σ ± ∏`, S645/S646).

**Formalized (math-lean, sorry-free): `card_union_three`** —
`|A∪B∪C| + |A∩B| + |A∩C| + |B∩C| = |A| + |B| + |C| + |A∩B∩C|` (the subtraction-free `card` form of
`|A∪B∪C| = (A+B+C) − (D+E+F) + G`). Verified on 2000 random triples.

---

## 1. The honest result on max-H (A038375)

Computed exactly (`max_hampaths_recursion_s651.py`, brute force `n ≤ 6`):

```
  A038375(n) = max #Hamiltonian paths over tournaments on n nodes:  1, 1, 3, 5, 15, 45  (n = 1..6)
  (Rédei: min = 1 transitive, all H odd; these are the maxima.)
```

> **The literal `3,−3,1` inclusion–exclusion does NOT reproduce max-H** (predicts `7, 7, 33` vs actual
> `5, 15, 45`). Max-H is an **irregular extremal sequence** with no simple linear recurrence: even the
> `n→n−2` ratios break (`×5, ×9` for `n=5,6` but `A038375(7)=189 ≠ 5·15`). So the user's tiling is a
> **dictionary** (the alternating-sign structure), not a closed formula for the maximum. Being honest:
> efficiently computing the *maximum* is hard precisely because the extremal tournament is irregular.

What *does* recurse cleanly is **`H` itself** (not its max) — see §2.

---

## 2. The recursive truths that DO hold (the Pfaffian angle, five faces)

The user's instinct is right that one alternating-sign recursion underlies everything; it just lives in
`H`, the Pfaffian, and the inclusion–exclusion, not in the irregular maximum:

| domain | the recursion | signs |
|---|---|---|
| **algebra** | `det(skew-adj) = Pf²`; `Pf(Mₙ) = Σⱼ (−1)ʲ M₁ⱼ Pf(M_{1̂ĵ})` (cofactor, `n→n−2`) | `±` |
| **combinatorics** | perfect matchings (`Pf`); `H = I(Ω,2)` (independence polynomial at 2) | — |
| **graphs/tournaments** | `I(Ω,x) = I(Ω∖v,x) + x·I(Ω∖N[v],x)` — deletion–contraction (`n→(n−1, n−1−deg)`, S625) | `+` |
| **geometry** | the staircase/triangle tiling; corners / edges / interior | — |
| **topology** | `χ = V − E + F`; the chain-complex alternating signs (the user's `+,−,+`) | `±` |

> **The unifying object is the ALTERNATING SUM over a graded structure: `+dim0 − dim1 + dim2`.** The
> Pfaffian, 3-set inclusion–exclusion, the Euler characteristic, and the user's triangle are *one
> signed-sum recursion in five languages*. The Pfaffian recursion (`n→n−2`, S645/646) and the
> deletion–contraction for `H` (`n→n−1`, S625) are the two genuine recursions; the user's `n−1,n−2,n−3`
> triangle is the inclusion–exclusion that *combines* the `n−1` and `n−2` levels with the `n−3` correction
> — the Euler-characteristic of the deletion poset.

---

## 3. New recursive truths surfaced
- **`H` recurses, `max H` does not.** The clean recursion is deletion–contraction on the independence
  polynomial (`H = I(Ω,2)`); the *maximum* over tournaments is irregular (A038375 has no simple formula).
  The efficient route to `max H` is therefore not a closed recurrence but *generating the extremal
  tournaments* — itself an open/hard extremal problem.
- **The Pfaffian is the algebraic Euler characteristic.** `det = Pf²` is the determinant (a *full*
  alternating sum over `Sₙ`) factoring through the matching sum; the `±` of the Pfaffian = the `±` of
  inclusion–exclusion = the `±` of `χ`. The user's tiling makes this geometric.
- **The triangle's `χ = 1`** (`3 − 3 + 1`) is the topological invariant of the decomposition: any
  tournament quantity that is *additive over the cells with these signs* (an "Euler-characteristic-type"
  invariant) inherits the inclusion–exclusion — that is the class of quantities the user's tiling computes
  efficiently (e.g. the `card`/measure functionals, the depth-GF overlaps of S618), *not* the maximum.

## 4. Handoffs
- Formalize the Pfaffian cofactor recursion (`n→n−2`) and the deletion–contraction for `I(Ω,2)` (needs a
  Pfaffian / the independence-polynomial recursion in Lean).
- Characterize exactly which tournament functionals are "Euler-characteristic-additive" over the triangle
  (those the user's tiling computes) — the depth-GF / `card` family (S618) vs the extremal `max H`.
- A038375: catalogue the extremal tournaments (the irregular maximizers) — the genuinely hard part.
