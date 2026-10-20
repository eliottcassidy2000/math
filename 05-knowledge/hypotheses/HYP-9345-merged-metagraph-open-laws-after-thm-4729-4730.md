---
id: HYP-9345
title: "Open laws of the merged tournament metagraph after THM-4729/4730: the closure conjecture (H > n HC for every strong tournament of order >= 6), connectivity, the blue self-line classes, the black self-line law, the determination problem, and the strong-tiling sequence"
status: OPEN (pre-audit); each item carries its finite evidence and the first untested case
source: mac-mini-2026-10-20-metagraph
depends_on:
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
related:
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-793-mode-b-line-tower.md
  - 01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md
---

# HYP-9345 -- what the merged metagraph still hides

Setting and names: THM-4729 (closed/open tilings, closure number `cl(C) = n HC/|Aut|`,
fully closed / fully open nodes, first-cut count `a(n)`, `NS(n)`), THM-4730 (points of
symmetry = pure-blue nodes; strong ones are `K1, C3, T5`).

## (a) The closure conjecture

**Claim.** Every strongly connected tournament of order `n >= 6` has a Hamiltonian path
that does not close to a Hamiltonian cycle: `H(T) > n HC(T)`. Equivalently (THM-4729 L2)
no node of the metagraph with `n >= 6` is fully closed.

Evidence: exhaustive over all strong classes for `n <= 9` (fully closed classes: `C3` at
`n = 3`, `Q5 = a -> C3 -> b -> a` at `n = 5`, none at `4, 6, 7, 8` and none among the
178,133 strong classes at `n = 9`, `fully_closed_n9.out`). First untested case `n = 10`. A proof would need an open
Hamiltonian path in every strong tournament of order `>= 6`; the two small exceptions
have the pattern "one vertex beats a rotational core, one vertex is beaten by it".

## (b) Connectivity

**Claim.** For every `n >= 3` the merged metagraph is connected, and so is the blue graph
(SC nodes, non-loop blue lines). Evidence: `n = 3..8` exact (THM-4729 census). THM-4729
L3 shows every component contains strong nodes; a proof must connect the strong world
through complement lines and class identifications.

## (c) Blue self-line classes

Counts `1, 2, 4, 24` at `n = 4, 6, 8, 10` (klein-S169/S173; zero at odd `n`, THM-023/648).
**Claims.** (c1) Every class carries at most one blue self-line (true for `n <= 10`).
(c2) Every blue self-line class is strong and near-regular, scores
`(n/2-1)^{n/2} (n/2)^{n/2}` (true for `n <= 10`; only the two path ends are forced to
`n/2` and `n/2 - 1` by the score argument). (c3) A mechanism for the counts; the `n = 12`
value is the first discriminator (the `2^{n/2-2}` and `(n/2-1)!` guesses are both dead,
MISTAKE-124).

## (d) Black self-lines

THM-793's conjecture `black self-lines (same class) = SC(n)/2` at odd `n` (`4 = 8/2`,
`44 = 88/2` at `n = 5, 7`) is still untested at `n = 9` (`2^28` tilings; a score-filtered
census is feasible in C). Merged black loops (same node, converse sheets included):
`0, 0, 4, 24, 114, 618` for `n = 3..8`, all even (THM-796 (4.7)); no formula is known.

## (e) Rigid strong tournaments of composite order >= 63

THM-2477's lemma ML. THM-4730 shows the count law for pure-blue nodes is equivalent to
ML beyond `n = 125`. A rigid strong tournament is vertex-transitive with `Aut` regular on
Hamiltonian paths, hence circulant with respect to every Hamiltonian cycle and with
`S n U(Z_n) = M` an index-two multiplier group not containing `-1`; the prime case is
`C3` only (the Paley tournament has `(p-3)/4 >= 1` non-progression Hamiltonian paths
through the arc `0 -> 1`, for `p >= 7`).

## (f) The strong-tiling sequence

`a(n) = 1, 0, 1, 5, 50, 903, 30773, 2032504` (strong tilings = labelled strong
tournaments on `[n]` containing a fixed Hamiltonian path) and `NS(n) = 2^{C(n-1,2)} - a(n)`
are not in OEIS. Questions: an asymptotic `NS(n) ~ 2^{C(n-2,2)} (1 + o(1))` (the `k = 1`
term dominates; `k = n-1` contributes `a(n-1)`), the full expansion, and whether
`a(n)` is the inverse-Euler-type transform of a known strong-tournament sequence.

## (g) The determination problem

For which `n` do the known laws (fibre sizes `H/|Aut|`, blue parities THM-643, closure
numbers and the closed/open bipartition THM-4729, mirror pairing THM-796, score
complement THM-646, cyclic-triangle flux THM-785) determine the merged metagraph up to
isomorphism? THM-4729 proves `n <= 4`; at `n = 5` the blue system is one of several
perfect matchings between closed and open symmetric tilings compatible with the degrees,
and it is open whether THM-646/785 single it out. A finite check at `n = 5, 6` would
settle the first cases.

## (h) Rotation lines

Within a strong class, the `n` rotations of a closed tiling along its Hamiltonian cycle
are closed tilings of the same class (an orbit of size `n/d`, `d` = the order of the
rotational symmetry of the cycle inside `Aut`). These "rotation lines" are loops of the
metagraph; their orbits partition the closed tilings of every strong node. Question:
does the rotation structure (orbit sizes, their interaction with `sigma` and with the
complement) explain which closed symmetric tilings pair with open symmetric ones of the
same class (the blue self-lines)?
