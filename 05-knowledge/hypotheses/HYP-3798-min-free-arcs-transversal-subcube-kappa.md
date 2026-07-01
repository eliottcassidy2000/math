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
Two equivalent views of the optimal configuration:

**(a) Free arbitrary arcs** -- a dense clique-packing with `C(n-2,2)+1` edges (unique up to free-graph iso
for `n<=5`): `n=4` a perfect **matching** `{01,23}`; `n=5` a **triangle+edge** `{012,34}`; `n=6` **two
triangles+bridge** `{012,345,03}`.

**(b) The clean rule INSIDE the tiling model** (`kappa_tiling(n) = kappa(n)`, verified `n<=6`): keep the
Hamiltonian path fixed AND freeze the `n-3` tiles of the **skip-2 diagonal** `{(i, i+3) : 0<=i<=n-4}` (the
line one step inside the staircase hypotenuse), freeing the other `1+C(n-2,2)` tiles. The frozen tiles carry
no isomorphism information once the rest of the board is set. So the owner's "configuration rule on the
fixed arcs" = *fix the Hamiltonian path and the skip-2 diagonal.* (Exact `n<=6`; at `n=7` the forward-fixed
skip-2 diagonal covers **454/456** classes -- a near-miss, so the clean diagonal rule is exact only through
`n=6` and `kappa(7)=11` needs a slightly different `11`-arc fixing, still open.)

## Related things to study (open)
1. **Prove `kappa(n) = 1 + C(n-2,2)`**: upper bound = a clique-packing construction with the right fixing
   (a transitive fixing FAILS at `n=7`: 225/456 -- the fixing is subtle); lower bound beyond the info-floor.
2. **Why `n-3` redundancy** in the Hamiltonian-path tiling -- which `n-3` tiles are always fixable?
3. **The info-floor gap** `kappa - ceil(log2 A000568)` (`0,0,0,1,2,...`): its growth `~ n log n - O(n)`.
4. **Coding-theory framing**: `beta` = covering radius, `kappa` = min transversal-subcube dimension, of the
   `S_n`-orbit coloring of `Q_{C(n,2)}` -- both are covering-code parameters (repo engineering mandate).
5. **The even-graph dual** `E_n` (A002854): the same two invariants for even graphs.
6. **`n=7` verification**: is `kappa(7)=11`? (needs a smarter search or a construction; currently open).
   **[opus-2026-07-01-S15 RESOLUTION — HYP-3805]: very likely `kappa(7)=12`, NOT 11 — the lazy-caterer breaks.**
   Best 11-config = fix spans {1,3} of a linear order (free span-2 + spans≥4) → 454/456 (better than the
   clique-packing's 225); it misses exactly the |Aut|=1 class `(1,2,2,3,4,4,5)` and the **Paley heptagon**
   `(3^7)` |Aut|=21. For that free-set, all 184 bases hosting both missing classes were enumerated — **none**
   reaches 456. `kappa(7)≤12` proven (fix span1 + span3−(0,3)). The obstruction = argmax|Aut| = the Paley/QR
   tournament on Z/7 = the LRC atoms (HYP-3802); mechanism = n!/|Aut| few-reps. So the formula `1+C(n-2,2)`
   holds only `n≤6`.

## Honest scope
`kappa=1,2,4,7` and `R=1,1,3,4` are EXACT (`n<=6`, exhaustive). `beta=MFAS` and `beta_path=beta_order` are
proved. The formula `kappa(n)=1+C(n-2,2)` and `kappa(7)=11` are CONJECTURES (fit `n<=6`; `n>=7` unverified,
random search inconclusive). A new invariant + a clean conjectural formula, not yet a theorem for all `n`.


> **CORRECTION (mac-mini-S90, HYP-3819):** the formula kappa(n)=1+C(n-2,2) gives kappa(7)=11 but the TRUE flip-rank kappa(7)=12 (opus/klein) -- the lazy-caterer formula BREAKS at n=7 (Paley-heptagon obstruction). The formula is exact only for n<=6.
