---
id: HYP-3809
title: CONJECTURE ATLAS for the merged-metagraph bucket/parity structure (extends kps S12/S675 blue-black work). Core new results: (1) the PARITY THEOREM REDUCES to a sigma-orbit count -- sigma (=complement-fold=transpose, the half-tiling mirror) PRESERVES merged nodes, so a node is a union of sigma-orbits (fixed=grid-sym, else pairs), tiling_count = #gridsym + 2*#pairs, hence parity = #gridsym mod 2 (SC has ODD #gridsym, NS has 0). (2) IDENTITY #grid-sym tilings = 2^floor((n-1)^2/4) = 2^(half-tiling quarter-square size) => the BLUE structure lives on the half-tiling. (3) SC classes = A051337 (self-converse tournaments) = #pure_blue+#mixed; V_merged=(A000568+A051337)/2. (4) All-orders even blue-count multiplicities are REFUTED at n=8, but the repaired ODD-order theorem is PROVED (2026-10-20): auxiliary hollow-frame toggle cubes preserve blue count and reverse cyclic-triangle parity; these are not blue-metagraph automorphisms. (5) STRUCTURE=CONSTRAINT is UNDER-DETERMINED: valid color+eligibility+degree-preserving 2-swaps exist (n=5,6), so the owner's parity/eligibility constraints are NECESSARY but NOT SUFFICIENT -- the metagraph's 'parity skeleton', with the tournament-iso structure selecting the actual realization.
status: MIXED -- several CONFIRMED (n<=6 exhaustive: grid-sym=2^floor((n-1)^2/4); SC=A051337; sigma preserves nodes; parity reduction; historical blue-count even-multiplicity checks in the stated finite range; 2-swap existence) + a CONJECTURE ATLAS (many small conjectures, some open). Extends kps S12/S675 (tripartite structure, tiling_count=degree, black-even, grid-sym<=>tau-fixed<=>SC, n=6 self-loop correction -- kps has PRIORITY on these; this session confirms + extends). The original note did not prove general parity; THM-643-blue-black-line-parity-structure supplies that theorem. The all-orders twin-pairing extrapolation is REFUTED, while the odd-order equal-shore theorem and free auxiliary toggle cubes are PROVED (2026-10-20 correction below).
source: mac-mini-2026-07-01-S84
related:
  - HYP-3808   # S83 the parity decomposition (black=even graph, blue=odd) -- this reduces/extends it
  - THM-549    # half-tiling: fold = complement = transpose; quarter-square floor((n-1)^2/4)
  - HYP-3798   # S81 minimal-flip kappa -- the n=6 gauge break = the sea-onset (same threshold)
  - HYP-3799   # S82 even graphs -- the black-line even graph; the even-graph lift target
results:
  - 04-computation/metagraph_bucket_conjecture_atlas_macmini_20260701.py
  - 05-knowledge/results/metagraph_bucket_conjecture_atlas_macmini_20260701.out
  - 04-computation/merged_metagraph_blue_black_lines_macmini_20260701.py
---

# HYP-3809 -- conjecture atlas for the merged-metagraph buckets

## Current correction and repaired odd-order theorem — 2026-10-20

The all-order extrapolation is **REFUTED**: at `n=8`, exactly **three**
self-converse classes have blue tiling count `1`. Thus there cannot be a
fixed-point-free class involution preserving blue count at **every** order.
The [exact order-eight spectrum](../results/merged_metagraph_black_quotient_bridges_20261020/README.md#two-further-hostile-controls-at-n8)
remains a counterexample; the theorem below does not erase the even-order
failure.

The repaired **odd-order** statement is now **PROVED for every odd `n>=3`**.
There is a class-level involution preserving blue tiling count and reversing
cyclic-triangle parity. More strongly, self-converse classes partition into
free auxiliary toggle cubes `(C2)^r`, constant in blue count and checkerboard
in triangle parity. For `n=2k+1`, the cube dimension satisfies `r>=1` and
`r=k (mod 2)`. Consequently:

- The blue-count spectra on the even- and odd-triangle sides are equal at
  every odd order, not only in the historical `n=3,5,7` data.
- Every blue-count multiplicity is even at odd order.
- At `n=1 (mod 4)`, `n>=5`, each full multiplicity is divisible by four,
  and each shore multiplicity is even.

The [hollow-frame theorem and proof](../results/odd_blue_spectrum_hollow_cubes_20261020/ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md)
delete only the matching edges of an anti-involution. An odd-order frame
automorphism group acts on these pairs, whose orbits have odd size.
Toggling one orbit preserves the centralizer and the fixed-reversal path
set, hence the blue count, while reversing triangle parity. Conjugacy of
anti-involutions makes the construction descend to unmarked classes.

This is **not** the old unproved `flip.sigma` candidate and is not a blue-edge
preserving map. It need not preserve purity, total fibre mass, `H`, or the
full automorphism group. At `n=7`, `C3 -> 1 -> C3` has `(B,mass)=(1,1)`
and is pure-blue, but its switched partner has `(B,mass)=(1,45)` and is
mixed; the two are not adjacent by an actual blue complement line.
The cubes are auxiliary class structures, not asserted subgraphs of the
original metagraph.

C4 and Atlas item 3 below retain the original conjecture and finite evidence.
Their universal version is refuted, and their odd-order refinement is now
closed by this different mechanism. General SC-odd-blue parity is supplied
by [THM-643, blue/black line parity](../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md).
The count/colour-preserving rewiring obstruction survives; no complete
metagraph realizability or unique-arrangement theorem follows.



Extends kps's blue/black metagraph work (`metagraph-blue-black-even-odd-s675.md`,
`the-blueblack-line-pairing-is-a-degree-tiling-count-realization-kps.md`) and my HYP-3808. Owner's asks:
generate multitudes of conjectures; test "structure = constraint"; understand which tilings share a node
(half-tiling symmetry). Setup: merged node = merged iso class; line `= {tiling t, flip(t)}` (`flip` =
complement-tiling), blue iff `isGridSym(t)`; tiling_count(node) = its colored line-degree.

## The symmetry frame: `<flip, sigma> = Z2 x Z2` (which tilings share a node)
Two involutions on the `2^m` tilings: **`flip`** (complement all tile bits -> the LINES; no fixed points)
and **`sigma`** (the anti-diagonal / `y=x` mirror = complement `T -> T^op` = transpose, THM-549 -> the
MERGE and grid-symmetry; fixed points = grid-sym tilings = the half-tiling). VERIFIED: **`sigma` PRESERVES
merged nodes** (`t` and `sigma(t)` always share a node), so a node is a union of `sigma`-orbits: singletons
(grid-sym) and pairs. Hence `tiling_count = #gridsym + 2*#pairs`, and **`tiling_count parity = #gridsym mod
2`.** This is the clean mechanism behind the parity theorem.

## CONFIRMED (n<=6 exhaustive)
- **C1 (parity theorem, kps + reduced here):** SC merged node `<=>` ODD tiling count; NS `<=>` EVEN.
  Reduces to: **an SC node has an ODD number of grid-sym tilings; an NS node has ZERO.** (NS-zero is clear;
  SC-odd is the crux.)
- **C2 (grid-sym identity):** `#grid-sym tilings = 2^{floor((n-1)^2/4)}` = `2^{half-tiling size}` (`4,16,64`
  for n=4,5,6). `#blue lines = 2^{floor((n-1)^2/4) - 1}`. So the BLUE structure lives entirely on the
  half-tiling (the `sigma`-fundamental domain, THM-549).
- **C3 (SC count):** `#SC classes = A051337` (self-converse tournaments) `= 2,2,8,12` `= #pure_blue +
  #mixed`; `#pure_black = (A000568 - A051337)/2 = #NS-pairs`; `V_merged = (A000568 + A051337)/2`.
- **C4 (twin-pairing, NEW):** the multiset `{#grid-sym tilings per SC node}` has ALL-EVEN multiplicities
  (`n=5: {1:4,3:4}`, `n=6: {1:2,3:2,5:2,7:4,9:2}`) -- a fixed-point-free involution twins SC nodes by
  grid-sym count. (The tiling-count multiset is NOT all-even -- the twinning is on grid-sym count.)
- **C5 (structure=constraint UNDER-DETERMINED):** valid color+eligibility+degree-preserving 2-swaps EXIST
  (n=5,6). So the parity/eligibility/degree constraints are **necessary but not sufficient** to pin the
  metagraph -- they are its "parity skeleton"; the tournament-iso structure selects the realization.
- **C6 (total lines):** `#lines = 2^{C(n-1,2)-1}`; `flip.sigma` has 0 fixed points (n=4,5,6).

## THE CONJECTURE ATLAS (open, ranked)
1. **Blue = half-tiling metagraph (kps target):** the blue subgraph on SC nodes is (isomorphic to / a cover
   of) the merged metagraph of the folded half-tiling; C2 (`#gridsym = 2^{half}`) is the first evidence.
2. **SC-odd-grid-sym (the parity crux):** every SC class contains an odd number of grid-sym tilings.
   Likely a unique "central" grid-sym representative (on a principal diagonal) forces the parity.
3. **Twin involution (C4):** identify the fixed-point-free involution on SC nodes preserving grid-sym count
   -- candidate `flip.sigma` or a residual `Z2` on the half-tiling; does it also pair pure_blue with mixed?
4. **Sea-onset criterion:** NS-NS black lines + pure-black self-loops onset at `n=6`; conjecture the onset
   is EXACTLY the minimal-flip gauge break (HYP-3798 `kappa` jump at n=6) and precedes the `n=7` diameter
   jump -- one threshold, three faces (metagraph sea, `kappa`-redundancy, mod-3 three-channel law).
5. **Category-count sequences:** `#pure_blue = 2,1,3,2` (n=3..6), `#mixed = 0,1,5,10`, `#pure_black =
   0,1,2,22`; find closed forms/recurrences (Mode-A/B). `#pure_black = (A000568-A051337)/2` is closed.
6. **Sea/self-loop growth:** `#sea = 0,0,0,290`; `#pb-selfloops = 0,0,0,24`; extend to n=7, find the law.
7. **Black Eulerian cycle structure `<->` H-spectrum:** the even black graph decomposes into cycles; do the
   cycle lengths / count encode the H-values (tiling counts)?
8. **Odd/even shape parity:** odd `n` half-tiling = square `k^2`, even = pronic `k(k-1)` (THM-549); does the
   pure_blue vs mixed split follow this parity (square -> more pure-blue)?
9. **Even-graph lift (S82):** the black-line subgraph is an even-DEGREE graph; do blue vs black lines have
   distinct footprints under the tournament -> even-graph (cycle-space) projection?
10. **Forbidden line-degrees:** are `7, 21` (forbidden H, apex-prime) also forbidden as node line-degrees?
11. **`flip.sigma` fixed-point-free (C6):** prove no tiling equals the complement-tiling of its transpose.
12. **Handshake corollary:** `#SC = #pure_blue + #mixed` is EVEN (blue all-odd-degree handshake) -- so
    `A051337(n)` is even for all `n>=3` (a parity statement about self-converse tournament counts).

## Concrete next targets
- **Prove C1 via C2/the crux:** show every SC class has an odd number of grid-sym (`sigma`-fixed) tilings
  (a unique central one + flip-pairs). This proves the whole parity theorem.
- **Prove/identify the C4 twin involution** and the blue=half-tiling recursion (Atlas 1,3).
- **Characterize the sea onset** (Atlas 4) and unify with `kappa` (HYP-3798) and the n=7 transitions.

## Honest scope
C1-C6 VERIFIED exhaustively n<=6 (C1 also in kps S12; C2/C3/C4/C5/C6 confirmed here). The atlas items are
CONJECTURES. The owner's "structure = constraint" is REFINED: the constraints are a necessary parity
skeleton, not a full determination (C5). Priority note: kps holds the core tripartite / degree = tiling-
count / black-even / n=6-self-loop results; this session confirms and extends with the sigma-orbit parity
reduction, the grid-sym=2^half identity, SC=A051337, the twin-pairing, the under-determination, and the atlas.
