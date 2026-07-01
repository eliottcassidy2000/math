---
id: HYP-3798
title: THE MINIMUM-FREE-ARCS INVARIANT kappa(n) = 1 + C(n-2,2) = A000124(n-3) (lazy-caterer numbers). In the arc-hypercube Q_{C(n,2)} colored by tournament iso class, kappa(n) is the minimum dimension of an axis-aligned subcube (choose k FREE arcs, FIX the other C(n,2)-k to specific orientations) whose 2^k tournaments realize ALL A000568(n) iso classes. EXACT n<=6: kappa = 1,2,4,7 (n=3..6); formula 1+C(n-2,2) = m(n)-(n-3) = m(n-1)+1 predicts 11,16 for n=7,8. The naive Hamiltonian-path tiling (m=C(n-1,2) free tiles) is REDUNDANT by exactly n-3 arcs. Optimal free-arc SHAPES are clique-packings: n=4 a perfect MATCHING {01,23}, n=5 a TRIANGLE+edge {012,34}, n=6 TWO TRIANGLES + bridge {012,345,03} -- C(n-2,2)+1 edges. Companion invariant: beta(T) = min tile-flips to express T's class = the MINIMUM FEEDBACK ARC SET (beta_path=beta_order always -- the median order is always a Hamiltonian path); covering radius R(n)=max beta = 1,1,3,4 (n=3..6).
status: CONFIRMED n<=6 (exhaustive search: kappa=1,2,4,7 EXACT, k=6 provably fails at n=6; beta=MFAS and R(n)=1,1,3,4 exact; optimal shapes + explicit configuration rule verified). CONJECTURAL n>=7 (formula predicts kappa(7)=11 but a randomized n=7 search did NOT reach full coverage -- upper bound unverified; info-floor ceil(log2 A000568(7))=9 is the lower bound). This is a NEW pair of tournament invariants (a coding-theory covering problem on the S_n-orbit-colored arc-hypercube), the direct answer to the owner's "n=4: 3 tiles naively but 2 free arcs with a configuration rule."
source: mac-mini-2026-07-01-S81
related:
  - THM-591    # metagraph observer/rotation (AGL(1) on arcs); the arc-hypercube is the ambient space
results:
  - 04-computation/min_flip_tiling_beta_macmini_20260701.py
  - 04-computation/min_free_arcs_transversal_subcube_macmini_20260701.py
  - 04-computation/kappa_shapes_and_n7_macmini_20260701.py
  - 05-knowledge/results/min_flip_tiling_beta_macmini_20260701.out
  - 05-knowledge/results/min_free_arcs_transversal_subcube_macmini_20260701.out
  - 05-knowledge/results/kappa_config_rule_macmini_20260701.out
---

# HYP-3798 -- the minimum-free-arcs invariant kappa(n) and its companion beta = MFAS

The owner's seed: for n=4 tournaments, all 4 iso classes are reached by flipping 3 tiles in the naive
Hamiltonian-path tiling, "but they may also be expressed by 2 arcs only, provided a configuration rule on
the fixed arcs." Study how the minimizing shape changes with n.

## The two invariants
Color the arc-hypercube `Q_{C(n,2)}` (one bit per arc orientation) by tournament isomorphism class.
- **beta(T)** = min number of tile-flips to express `T`'s class (over the best base Hamiltonian path)
  `= C(n,2) - max_{Ham-path order} #forward-arcs = the MINIMUM FEEDBACK ARC SET`. Proved along the way:
  `beta_path = beta_order` always (a minimum-feedback-arc order can never have a backward *consecutive*
  arc -- an adjacent swap would improve it -- so it is a Hamiltonian path). Covering radius `R(n) = max_T
  beta = 1, 1, 3, 4` for `n=3..6` (`~ n^2/4 - c n^{3/2}`, the max-MFAS asymptotic).
- **kappa(n)** = min dimension `k` of an axis-aligned subcube -- choose `k` FREE arcs, FIX the other
  `C(n,2)-k` -- whose `2^k` tournaments hit ALL `A000568(n)` iso classes. This is the owner's object.

## The result: kappa(n) = 1 + C(n-2,2) = A000124(n-3)  [lazy-caterer numbers]
Exhaustive search (exact, `n<=6`): `kappa = 1, 2, 4, 7`. The formula `kappa(n) = 1 + C(n-2,2)` fits and
gives `11, 16` for `n=7, 8`. Equivalent forms:
- `kappa(n) = m(n) - (n-3)` where `m(n)=C(n-1,2)` is the naive tile count: **the Hamiltonian-path tiling is
  redundant by exactly `n-3` arcs.**
- `kappa(n) = m(n-1) + 1` (a size-`n-1` tile count plus one).
- `kappa(n) = A000124(n-3)` -- the central polygonal / "lazy caterer" numbers `1,2,4,7,11,16,...`.

Info-theoretic floor `ceil(log2 A000568(n)) = 1,2,4,6,9` for `n=3..7`: `kappa` MEETS the floor for `n<=5`
and first EXCEEDS it at `n=6` (`7 > 6`) -- the subcube can be information-optimal only up to `n=5`.

## The optimal SHAPE and the configuration rule
The minimizing free-arc set is a dense clique-packing with `C(n-2,2)+1` edges (unique up to free-graph
iso for `n<=5`):
- `n=4`: a perfect **matching** `{01, 23}` (free-graph degseq `(1,1,1,1)`), fixed rule e.g. `0->2, 3->0,
  2->1, 3->1`. This is exactly the owner's "2 arcs, configuration rule."
- `n=5`: a **triangle + disjoint edge** `{012, 34}` (degseq `(1,1,2,2,2)`, 1 triangle), fixed rule
  `0->3,0->4,1->3,4->1,3->2,4->2`.
- `n=6`: **two triangles + a bridge** `{012, 345, 03}` (2 triangles).
The triangle is the natural "generator" (its 8 orientations give the cyclic-vs-transitive iso bit); the
count of triangles in the optimal shape grows with `n`.

## Related things to study (open)
1. **Prove `kappa(n) = 1 + C(n-2,2)`**: upper bound = a clique-packing construction with the right fixing
   (a transitive fixing FAILS at `n=7`: 225/456 -- the fixing is subtle); lower bound beyond the info-floor.
2. **Why `n-3` redundancy** in the Hamiltonian-path tiling -- which `n-3` tiles are always fixable?
3. **The info-floor gap** `kappa - ceil(log2 A000568)` (`0,0,0,1,2,...`): its growth `~ n log n - O(n)`.
4. **Coding-theory framing**: `beta` = covering radius, `kappa` = min transversal-subcube dimension, of the
   `S_n`-orbit coloring of `Q_{C(n,2)}` -- both are covering-code parameters (repo engineering mandate).
5. **The even-graph dual** `E_n` (A002854): the same two invariants for even graphs.
6. **`n=7` verification**: is `kappa(7)=11`? (needs a smarter search or a construction; currently open).

## Honest scope
`kappa=1,2,4,7` and `R=1,1,3,4` are EXACT (`n<=6`, exhaustive). `beta=MFAS` and `beta_path=beta_order` are
proved. The formula `kappa(n)=1+C(n-2,2)` and `kappa(7)=11` are CONJECTURES (fit `n<=6`; `n>=7` unverified,
random search inconclusive). A new invariant + a clean conjectural formula, not yet a theorem for all `n`.
