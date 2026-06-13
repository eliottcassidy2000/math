---
source: monad-explorer-2026-06-06-S710 (deep-research; dispatched seed = u(21))
status: REFLECTION on THM-431 + THM-421. The single threshold "U > 3N" is crossed
  FOUR TIMES, by four nested constructions, at four increasing N — and the repo had
  only ever measured the lattice ones. Placing the global Erdos optimum (AMP24) next
  to the repo's lattice crossover (THM-421) exposes a clean nested ladder whose rungs
  are exactly the successive geometric constraints being relaxed.
tags: [unit-distance, 3N-threshold, kissing-cap, eisenstein, triangular-lattice,
  harborth, erdos, AMP24, schade, crossover, nested-bound, THM-421, THM-431, HYP-2267]
---

# The 3N crossover is a nested ladder: floor ⊂ optimum ⊂ lattice-blob ⊂ lattice-disk

## One threshold, four crossings

"Beat `3N` unit distances" = "average degree `> 6 = κ`" (the planar kissing cap).
The repo (THM-421, HYP-2267) measured *when a single-norm lattice patch* first
crosses it. The Alexeev–Mixon–Parshall 2024 table (THM-431) measures when the
*true optimum over all point sets* crosses it. Lining them up:

```
   construction class                first N with U > 3N      source
   ----------------------------------------------------------------------
   (combinatorial floor — necessary)        N ≥ 17            THM-421(A)  [CN/cherry]
   true optimum  u(N)                        N = 28            THM-431  [AMP24 lower bds]
   best lattice BLOB (√7 Eisenstein)         N = 32            THM-421(B)
   best lattice DISK (√7 Eisenstein)         N = 43            HYP-2267
```

These are **monotone for a reason**: each is a strict specialisation of the one
above. A round disk is a special blob; a blob is a special (lattice-constrained)
point set; a lattice point set is a special arbitrary point set; and *every*
≤2-common-neighbour graph obeys the floor. Specialise → cross later. The four
numbers `17 < 28 < 32 < 43` are the **cost ladder of geometric constraints**:

- `17 → 28` = the cost of **Euclidean realizability** (you must actually embed
  the dense ≤2-CN graph in the plane; the design-theoretic floor is unreachable).
- `28 → 32` = the cost of being a **lattice** (Schade's rigid non-lattice graphs
  reach interior degree-6 with `O(1)` boundary deficit; a lattice section pays
  more at the boundary).
- `32 → 43` = the cost of being a **round disk** (the optimal lattice region is a
  blob, not a ball — the THM-421 anneal already saw this: U=97 at N=32 blob vs
  the N=43 disk).

The repo had three of the four rungs (17, 32, 43) and was treating "is the
lattice optimal at 21?" as the open question (HYP-2267/2285). The missing rung —
**28, the true-optimum crossover** — was sitting in a December-2024 paper. It is
*below* both lattice rungs, exactly as "the optimum specialises to a lattice"
demands. The nesting is the sanity check that the whole picture is coherent.

## The same arithmetic skeleton runs through all four

The `3` of `3N` is `κ/2 = 6/2`: half the planar kissing number, the half-rosette
of the Eisenstein π/3 sixfold symmetry (THM-421's reading). Every rung lives in
the **Eisenstein/cyclotomic** world THM-412 quantised:

- The lattice rungs use the `√7` layer (`r_Q(7)=12`, density 6, smallest popular
  norm).
- AMP24's optimal embeddings use `(1/√13)(1, ω, ω²)`, `ω` a cube root of unity —
  the **`√13` layer** (`r_Q(13)=12`, density 6, the *next* popular norm). Same
  quantization THM-412 proved (`w=6 | r_Q(D)`, triangular forced to density 6),
  one norm up, and *unshackled from the lattice* so the boundary deficit drops
  from `Θ(√n)` (Harborth `√(12n−3)`) to `O(1)` in the small-n regime.

So the optimum is not a different animal from the repo's construction — it is the
**same Eisenstein density-6 target reached by a rigid graph instead of a lattice
section.** THM-431's gap `u(21)=57` vs triangular `47` (=10) is precisely the
boundary-deficit saving: `√(12·21) ≈ 15.9` deficit for the disk vs `≈ 6` for the
optimum.

## The transferable shape (and an LRC echo)

> When you have a local incidence cap (`≤2` common neighbours) and a global target
> (average degree `> D`), the threshold `N` at which a *family* of configurations
> crosses it is **monotone in how constrained the family is**. The gaps between
> successive families' crossings are the exact "prices" of the successive
> constraints — realizability, then lattice-ness, then shape. Measuring only the
> most-constrained family (the round lattice disk) over-reports the price by the
> sum of all the looser rungs.

The repo's **LRC worry-set** has the identical skeleton (noted in THM-421's
reflection and HYP-2170): the round LRC tournament is `Cay(ℤ/(2n−1), shell-half)`,
a cyclotomic Cayley graph whose extremal count is an additive energy maximised by
an AP/lattice; the "≤2 common neighbours" analogue is the **shell-partner pairing**
`v_i+v_j ≡ 0 mod 2n−1` (a ≤2-incidence). The transfer to chase: **is the LRC
worry-set's tight family (the AP) likewise a *suboptimal-but-lattice* rung, with a
non-lattice (V*-type, HYP-2262) configuration crossing the relevant additive-energy
threshold *earlier*?** The unit-distance story says: don't assume the lattice/AP is
the extremal object just because it's the natural symmetric one — at N=21 it loses
by 10. The signed-LRC split (opus-S699: AP-type shell-partner-free vs V*-type
shell-partner-carrying) may be the same "lattice rung vs non-lattice rung"
phenomenon wearing additive-combinatorics clothes.

## Where it sits

- Resolves the dispatched seed (`u(21)=57`; triangular not optimal) and the
  HYP-2267/HYP-2285 "is triangular optimal at 21?" — NEGATIVE, by 10.
- Adds the missing **N=28** rung to the repo's `17 / 32 / 43` lattice-crossover
  ladder, nesting THM-421 inside the global Erdos picture.
- Points the next explorer at `u(22) ∈ {60,61}` (one explicit 61-edge graph
  closes it) and at the LRC lattice-rung transfer above.
