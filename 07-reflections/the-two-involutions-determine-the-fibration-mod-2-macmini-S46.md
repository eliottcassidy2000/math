# The two involutions determine the fibration mod 2

*mac-mini-2026-07-07-S46 (HYP-4977, THM-643). Owner: explore the even/odd duality in the
blue/black line structure; count lines and node types at every n; let niche insights
coalesce into a grand restricting, unifying picture.*

## The picture that emerged

Everything the owner asked about is governed by TWO COMMUTING INVOLUTIONS on the tiling
cube GF(2)^m, interacting with the S_n fibration:

- **σ (grid transpose, linear):** its fixed space = the blue tilings, dimension
  (m+f)/2. Downstairs it is the CONVERSE. It sees the SC/NS division of nodes.
- **τ (flip, the +𝟙 translation, fixed-point-free):** its orbits = the LINES.
  It commutes with σ, so it preserves blueness: blue lines are τ inside Fix(σ).

Every counting question then falls to parity bookkeeping against ONE arithmetical
input — Rédei's theorem (H odd) with its companion (|Aut| odd):

1. Fibers are odd (H/|Aut|).
2. σ restricted to an SC fiber is an involution on an odd set ⟹ odd fixed points ⟹
   **every SC node carries an odd number ≥ 1 of blue tilings; non-SC nodes carry none.**
   Pure-black = non-SC is now a THEOREM at all n, not an n≤7 observation.
3. τ restricted to any fiber color-class is fixed-point-free ⟹ cross-endpoint counts
   inherit the tiling-count parities: **every node leaks an odd number of cross-class
   line-endpoints** — the d=m layer isolates nothing; SC nodes always send a blue line
   OUT (blue min-degree ≥ 1 on the SC world); the owner's "blues contribute odd
   amounts and blacks even amounts" is exactly the SC-node allocation parity, and it
   FLIPS on NS nodes (blue 0, black odd) so that the fiber total stays odd.

## The new object: H_sym

The blue count of an SC node is not just a parity — it is a new odd tournament
invariant, the **self-converse Hamiltonian path count** H_sym (paths fixed by
converse∘reversal). Its mass formula Σ H_sym-fibers = 2^{(m+f)/2} is the blue companion
of Σ H/|Aut| = 2^m. Its observed spectrum is all-odd with per-n maximum 1,3,3,9,9 —
conjecturally 3^{⌊(n−2)/2⌋} (C1): Rédei mod 2 seems to refine to a 3-power cap on the
symmetric sub-count. If C1 is true, the "odd" of Rédei and the "3-power" of the
symmetric count are two rungs of one ladder — the even/odd duality deepening into a
2-vs-3 duality on the symmetric locus.

## Niche insights that coalesced

- **Blue self-loops exist only at even n** (0,1,0,2,0): the 𝟙-translation shifts score
  sequences in a way that seems arithmetically un-matchable at odd n (C2 — provable-
  looking via the deterministic score shift).
- **Pure-blue classes are the maximal-symmetry classes** (fiber ∈ {1,3}; H ∈
  {|Aut|, 3|Aut|}): transitive + circulant-like. They sit as LEAVES of the blue graph.
  The blue world is a thin tree-like halo around the MIXED SC core — the strict-
  definition refinement of the old "spine" picture.
- **#gs per class scales like √fiber** (3↔5, 5↔17, 7↔29..37, 9↔41): the involution's
  fixed-point count sits at the random-model scale — the fibration is "as generic as
  parity allows", echoing the S211 finding that edges are color-blind.
- The same two involutions run the runner world downstairs (THM-639: step reversal =
  σ's shadow, window complement = τ's) — one algebra, two projects.

## What "knowing the structure completely" now needs

Sizes: done (H/|Aut|). Colors: done (H_sym + complement). Masses: done (2^{(m+f)/2},
2^{m−1}). Parities of every allocation: done. REMAINING: the explicit cross-class
matching (the line-metagraph) and C1–C4. The matching is now heavily constrained —
every node odd out-degree, blue confined to SC, allocation tables fixed — so the
completion problem is a parity-constrained perfect-matching classification, not an
open-ended census. That is the restricting, unifying picture: **the fibration mod 2 is
fully determined by Rédei + two commuting involutions; what remains is one matching.**
