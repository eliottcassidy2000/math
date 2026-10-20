---
id: THM-4729
title: "Closed and open tilings: the Hamiltonian-cycle bipartition of the metagraph lines and the forced arrangement of the merged metagraph. (1) A tiling is CLOSED when the tile (n,1) carries the arc 1 -> n, i.e. the base path closes to a Hamiltonian cycle, OPEN otherwise; the complement line flips that tile, so every line, blue or black, joins a closed tiling to an open one; there are 2^(m-1) closed tilings and 2^((m+f)/2-1) closed grid-symmetric ones. (2) Closure number: a class C has exactly n HC(C)/|Aut C| closed tilings (HC = number of Hamiltonian cycles; Aut acts freely on rooted Hamiltonian cycles), a merged node twice that; a class is strong iff its closure number is positive (Camion). (3) Forced arrangement: non-strong nodes are fully open, hence pairwise non-adjacent and loopless -- every line at a non-strong node ends at a strong node; fully closed nodes (every Hamiltonian path closes, H = n HC) are pairwise non-adjacent and loopless; every self-line lives at a strong node having both closed and open tilings; for n >= 6 every pure-blue node is non-strong, so all its blue lines end at strong mixed nodes. (4) The number of lines between the strong and the non-strong world equals the number of non-strong tilings NS(n) = sum_{k=1}^{n-1} a(k) 2^C(n-k-1,2) (first-cut decomposition), a(n) = 2^C(n-1,2) - NS(n) the number of strong tilings = sum over strong classes of H/|Aut|: a = 1,0,1,5,50,903,30773,2032504 and NS = 0,1,1,3,14,121,1995,64648 for n = 1..8 (neither in OEIS); strong-strong lines = 2^(m-1) - NS(n). (5) Convention independence: complementing the non-path arcs, reversing the path arcs, and complementing then mirroring induce the same line system on the MERGED metagraph (they differ on the unmerged class graph): the converse merge is exactly the quotient on which the line system is convention-free. (6) FINITE-EXACT n = 4..8: all laws, connectivity of the metagraph and of the blue graph on SC nodes, fully closed classes only C3 (n = 3) and the (H,|Aut|) = (15,3) class a -> C3 -> b -> a (n = 5), and the blue self-line classes n = 4..10 (one self-line per class, near-regular scores, 24 classes at n = 10 reproducing klein-S173)"
status: >
  AUDIT SCOPE UPDATE (2026-10-20): L4's all-order pure-blue conclusion is
  CONDITIONAL on repairing THM-4730's odd-order anti-automorphism proof gap.
  L1-L3, L5, F and K do not use that gap. Original pre-audit claims follow:
  PROVED (pre-audit): Laws L1-L5 and Proposition F (first-cut count), Proposition K (convention independence);
  the n = 4 metagraph is forced by the laws (worked example).
  FINITE-EXACT (closed_open_laws_n7.out, closed_open_laws_n8.out): every law checked on all 2^m tilings for
  n = 4..8 (3,528 nodes, 1,048,576 lines at n = 8): zero violations of the closure number, zero lines inside the
  non-strong world, zero loops at non-strong nodes, blue non-loop degree odd at every SC node, strong/non-strong
  line counts 3, 14, 121, 1995, 64648 = NS(n), metagraph and blue graph connected, fully closed nodes 0, 1, 0, 0, 0
  (n = 4..8) and none among the 178,133 strong classes at n = 9 (fully_closed_n9.out), pure-blue inventories as THM-2444/2454; blue_self_lines_classes.out: the blue self-line classes at
  n = 4, 6, 8, 10 with (H, |Aut|, HC, t, s, closure, scores).
  CITED: definitions.md (tiling model), LEM-003 (free action), THM-643 (T1-T6, parities), THM-644(a), THM-646
  (score-complement law on lines), THM-787 (the transitive node's single blue line), THM-793/THM-849 (self-line
  censuses), THM-796 (half-edge parities), THM-2444/2454/2477 and THM-4730 (pure-blue structure), klein-S173
  (24 blue self-loops at n = 10, HYP-4971/MISTAKE-124), Camion 1959, Redei 1934.
  NOT PROVED: that the metagraph (or the blue graph) is connected for every n; that no fully closed strong
  tournament of order >= 6 exists (H > n HC); any closed form for a(n); the determination problem (whether the
  known laws fix the metagraph for n >= 5) -- all in HYP-9345.
source: mac-mini-2026-10-20-metagraph (owner's prompt: comprehensive patterns of the merged metagraph, counts of nodes, tilings and coloured lines, new lines between tilings, forced arrangements, new terminology)
depends_on:
  - 01-canon/definitions.md
  - 01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-644-gridsym-fiber-law-and-anti-redei.md
related:
  - 01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md
  - 01-canon/theorems/THM-787-transitivity-flow-on-the-metagraph.md
  - 01-canon/theorems/THM-793-mode-b-line-tower.md
  - 01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md
  - 01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
scripts:
  - 04-computation/experiments/metagraph_20261020/closed_open_laws.py
  - 04-computation/experiments/metagraph_20261020/blue_self_lines_classes.py
  - 04-computation/experiments/metagraph_20261020/merged_metagraph.py
  - 04-computation/experiments/metagraph_20261020/fully_closed_n9.py
outputs:
  - 04-computation/experiments/metagraph_20261020/fully_closed_n9.out
  - 04-computation/experiments/metagraph_20261020/closed_open_laws_n7.out
  - 04-computation/experiments/metagraph_20261020/closed_open_laws_n8.out
  - 04-computation/experiments/metagraph_20261020/blue_self_lines_classes.out
  - 04-computation/experiments/metagraph_20261020/merged_metagraph_n7.out
  - 04-computation/experiments/metagraph_20261020/merged_metagraph_n8.out
---

# THM-4729 -- every line of the metagraph joins a closed tiling to an open one

**PROVED (pre-audit) + FINITE-EXACT** as itemized in the status.

## 0. Setting (01-canon/definitions.md, THM-643)

Base path `P: n -> n-1 -> ... -> 1`; tiles `(a,b)`, `a >= b+2`, `m = C(n-1,2)`; a tiling
`t in {0,1}^m` orients each tile (`0`: `a -> b`, `1`: `b -> a`) and gives the tournament
`T_t`. Grid involution `sigma` = the relabelling `i -> n+1-i`, with `f = floor((n-1)/2)`
fixed tiles; `t` is grid-symmetric (blue) iff `sigma t = t` iff `i -> n+1-i` is an
anti-automorphism of `T_t` (THM-644(a)); `T_{sigma t} = T_t^rev` up to relabelling.
Line: `{t, t^c}` with `t^c` the all-tiles flip; blue iff `t` is grid-symmetric (then so is
`t^c`). Node: isomorphism class of `T_t`, merged with the converse class when distinct
(SC = self-converse node, merged = converse pair). THM-643: the fibre of a class `C` has
`tc(C) = H(C)/|Aut C|` tilings (odd), a merged node `2 tc(C)`; blue tilings sit only at
SC nodes, an odd number per SC node; `sum_nodes tc = 2^m`, blue tilings `2^{(m+f)/2}`.

## 1. Closed and open tilings; the closure number

**Definition.** The tile `(n,1)` joins the two ends of the base path. A tiling is
**closed** when it orients `(n,1)` as `1 -> n` (bit `1`), so that `P` together with this
arc is a Hamiltonian cycle of `T_t`, and **open** when `n -> 1` (bit `0`). The
**closure number** of a class is the number of closed tilings in its fibre, `cl(C)`; of a
node, `c_u` (the sum over its one or two classes); `o_u = t_u - c_u` is its open count.

**L1 (one closed end per line).** `t -> t^c` flips every tile, in particular `(n,1)`.
So every line has exactly one closed and one open end, and there are `2^{m-1}` closed and
`2^{m-1}` open tilings. The mirror `sigma` fixes `(n,1)` (`(n,1) -> (n+1-1, n+1-n)`),
so complementation is a bijection between closed and open grid-symmetric tilings:
`2^{(m+f)/2-1}` of each; blue lines, too, join a closed to an open tiling.

**L2 (closure number).** `cl(C) = n * HC(C) / |Aut(C)|`, where `HC(C)` is the number of
Hamiltonian cycles of `C`. *Proof.* The tilings of `C` are the Hamiltonian paths of `C`
modulo `Aut(C)` (LEM-003: the action is free; the tiling is the labelled tournament
obtained by naming the path's vertices `n, n-1, ..., 1`). A tiling is closed iff its path
`(v_1, ..., v_n)` has `v_n -> v_1`, i.e. iff the path is a Hamiltonian cycle with a chosen
starting vertex; there are `n HC(C)` rooted Hamiltonian cycles, `Aut` acts freely on them,
and two rooted cycles give the same tiling iff an automorphism carries one to the other. `#`
For a merged node `c_u = 2 cl(C)` since `HC(C^rev) = HC(C)`. Camion's theorem (a
tournament is strong iff it has a Hamiltonian cycle, `n >= 3`) gives:
`c_u > 0 <=> u is strong`. Summing L2 over classes: `sum_C n HC(C)/|Aut C| = 2^{m-1}`.

## 2. The forced arrangement (L3-L5)

**L3 (the non-strong world is a loopless independent set).** A non-strong node has
`c_u = 0`: all its tilings are open. Every line at it therefore has its closed end at a
node with a closed tiling, i.e. at a strong node. Hence no line joins two non-strong
nodes, no self-line sits at a non-strong node, and the whole line system is a bipartite
pairing between the closed half and the open half of the tiling space. Dually, a
**fully closed** node (`o_u = 0`, i.e. every Hamiltonian path of its class closes,
`H = n HC`) is adjacent only to nodes carrying open tilings, two fully closed nodes are
never adjacent, and a fully closed node has no self-line. Every self-line (blue or black)
lives at a strong node with `c_u >= 1` and `o_u >= 1`.

**L4 (CONDITIONAL on the repaired odd-order part of THM-4730; points of symmetry
are leaves of the closed world for n >= 6).** The 2026-10-20 weighted-core audit
found an affine/nonaffine anti-automorphism gap in that dependency. The following
argument is retained conditionally. By THM-2454 and
THM-4730 a pure-blue class is a palindromic stack of `{K1, C3}` halves with centre in
`{none, K1, C3, T5}`; it is strong only for `C3` (`n = 3`) and `T5` (`n = 5`). So for
`n >= 6` every pure-blue node is non-strong, all its tilings are open, and all its blue
lines end at closed grid-symmetric tilings of strong SC nodes with both colours (mixed
nodes), one line per tiling (`t_u = 1` or `3`), an odd number of them non-loop
(THM-643 T4), hence all non-loop. THM-787's unique blue line of the transitive node is the
case `t_u = 1`: the all-`0` tiling is open, its complement the closed all-`1` tiling of
the strong tournament `R_n` (path plus every backward chord), which is mixed.

**L5 (closure degrees).** At a node `u` exactly `c_u` line-ends are closed and `o_u` are
open; a self-line uses one of each. Hence a node with `c_u = 0` or `o_u = 0` has no
self-line (L3), and blue self-lines at an SC node are at most
`min(sc_u, so_u)` where `sc_u` counts its closed grid-symmetric tilings and
`so_u = s_u - sc_u`. In the census every blue self-line class carries exactly one blue
self-line (Section 5).

## 3. Proposition F: the first-cut count of non-strong tilings

`T_t` is not strong iff some cut `k in [1, n-1]` has every tile between the blocks
`{n, ..., n-k+1}` and `{n-k, ..., 1}` oriented downward (bit `0`): a non-strong
tournament's strong components are linearly ordered and every Hamiltonian path, in
particular `P`, traverses them in order, so the initial component is an initial segment
of `P`; conversely such a cut prevents a return to the top block. Taking the FIRST cut
(the initial strong component, of size `k`) gives the bijective decomposition

```text
NS(n) = sum_{k=1}^{n-1} a(k) * 2^{C(n-k-1, 2)},      a(n) = 2^{C(n-1,2)} - NS(n),         (F)
```

where `a(k)` is the number of strong tilings on `k` vertices (the top block carries a
strong tiling of its sub-path, the cross tiles are forced, the bottom block is arbitrary;
`C(0,2) = 0`). By LEM-003, `a(n) = sum_{strong C} H(C)/|Aut(C)|` and
`NS(n) = sum_{non-strong C} H(C)/|Aut(C)|`:

```text
n     : 1  2  3  4   5    6      7        8
a(n)  : 1  0  1  5  50  903  30773  2032504     (strong tilings; not in OEIS)
NS(n) : 0  1  1  3  14  121   1995    64648     (non-strong tilings; not in OEIS)
```

By L3 every non-strong tiling is the open end of a line whose closed end is at a strong
node, so **the number of lines between the strong and the non-strong world is exactly
NS(n)**, and the number of strong-strong lines is `2^{m-1} - NS(n)`; the census values
`3, 14, 121, 1995, 64648` (n = 4..8) agree with (F). The non-strong world carries
`NS(n)/2^m -> 0` of the tilings (`3.1%` at n = 8) yet `12.6%` of the nodes at n = 8
(443 of 3,528): non-strong classes have small fibres.

## 4. Proposition K: the merge is the convention-free quotient

Three natural "flip everything but one thing" operations on tilings: (i) complement the
non-path arcs, `t -> t^c` (the explorer's line); (ii) reverse the path arcs only,
`T_t -> T_t with P reversed`; (iii) complement then mirror, `t -> sigma(t^c)`. Reversing
the path arcs of `T_t` gives `(T_{t^c})^rev` (reverse everything, then undo the non-path
reversal), whose class is the converse of `[T_{t^c}]` and whose tiling along the reversed
path is `sigma(t^c)`. So (ii) and (iii) join `[T_t]` to `[T_{t^c}]^rev` while (i) joins it
to `[T_{t^c}]`: on the unmerged class graph the three conventions give different line
systems, on the merged metagraph (converse classes identified) the same one. The owner's
merge is therefore forced by convention independence, not only by the `2 x odd` fibre
count of converse pairs. Under (i) the mirror `sigma` pairs black lines with equal merged
endpoints (THM-796 (4.7)); under (ii) the same pairing is realised by `t -> t^c`.

## 5. Census (FINITE-EXACT n = 4..8) and the worked examples

```text
n | nodes (strong / non-strong) | lines strong-strong / strong-nonstrong | fully closed | blue loops | black loops | pure-blue
4 |    3 (   1 /   2)           |      1 /     3                          | 0            | 1          |   0         | 1
5 |   10 (   6 /   4)           |     18 /    14                          | 1  (15,3)    | 0          |   4         | 3
6 |   34 (  22 /  12)           |    391 /   121                          | 0            | 2          |  24         | 2
7 |  272 ( 216 /  56)           |  14389 /  1995                          | 0            | 0          | 114         | 4
8 | 3528 (3085 / 443)           | 983928 / 64648                          | 0            | 4          | 618         | 3
```

All laws hold with zero violations at every `n`; the metagraph is connected and the blue
graph on SC nodes is connected for `n = 4..8`. Fully closed classes: `C3` (`n = 3`) and
at `n = 5` the class `Q5 = a -> C3 -> b -> a` with `(H, |Aut|, t, s) = (15, 3, 5, 1)`,
all five tilings closed; none at `n = 4, 6, 7, 8`, and none among the 178,133 strong
classes at `n = 9` (`fully_closed_n9.py`; `T5` has `HC = 2`: two closed and one open tiling).

**n = 4 is forced.** Nodes: `TT4` (SC, non-strong, `t = s = 1`), the pair
`{v -> C3, C3 -> v}` (merged, non-strong, `t = 2`), `S4` (strong, SC, `t = 5, s = 3,
c = 4`). L3 sends the three non-strong tilings to `S4`: the symmetric one as a blue line
`TT4 -- S4`, the two asymmetric ones as the mirror pair of black lines `pair -- S4`. `S4`'s
remaining two tilings are both symmetric (`3 - 1`) and one is closed, one open: a blue
self-line. Lines `2 blue + 2 black`, degrees `1, 2, 5`: the unique multigraph.

**n = 5 (table; closure column added to THM-643's data).** Nodes `0..9`; `c` = closed
tilings, `sc` = closed symmetric:

```text
node | type   | strong | H  | Aut | t  | s | c  | sc | loops (blue/black)
  0  | SC     | no     |  1 |  1  |  1 | 1 |  0 | 0  | 0/0   TT5
  1  | merged | no     |  3 |  3  |  2 | 0 |  0 | 0  | 0/0   (1,1,C3) pair
  2  | merged | no     |  5 |  1  | 10 | 0 |  0 | 0  | 0/0   (1,S4) pair
  3  | SC     | no     |  3 |  3  |  1 | 1 |  0 | 0  | 0/0   (1,C3,1)  pure-blue
  4  | SC     | yes    | 15 |  3  |  5 | 1 |  5 | 1  | 0/0   Q5, fully closed
  5  | SC     | yes    | 11 |  1  | 11 | 3 |  5 | 1  | 0/2
  6  | SC     | yes    |  9 |  1  |  9 | 1 |  5 | 1  | 0/0
  7  | SC     | yes    |  9 |  1  |  9 | 3 |  5 | 1  | 0/0
  8  | SC     | yes    | 13 |  1  | 13 | 3 | 10 | 2  | 0/2
  9  | SC     | yes    | 15 |  5  |  3 | 3 |  2 | 2  | 0/0   T5, pure-blue
blue : 0-7, 3-8, 4-9, 5-6, 5-7, 5-9, 7-8, 8-9
black: 1-6 x2, 2-4 x2, 2-5 x2, 2-6 x2, 2-7 x2, 2-8 x2, 4-7 x2, 5-5 x2, 5-8 x2, 6-7 x2, 6-8 x2, 8-8 x2
```

The 14 non-strong tilings (nodes 0-3) all attach to strong nodes (L3); the fully closed
`Q5` meets only open tilings (`T5`'s one open tiling, node 2, node 7); `T5` is adjacent to
the two closed symmetric tilings of `8` and `5` and to `Q5`. The laws fix all degrees and
the closed/open split but not the pairing: the blue system is one of several perfect
matchings between closed and open symmetric tilings compatible with the degrees -- the
determination problem of HYP-9345.

## 6. The blue self-line classes (n = 4, 6, 8, 10)

Blue self-lines exist only at even `n` (THM-023, THM-648; score argument: the complement
sends scores to `n-1-s` except `n - s(start)`, `n-2-s(end)`, and self-converseness forces
`s(start) = n/2`, `s(end) = n/2 - 1`). Counts `1, 2, 4, 24` (klein-S169/S173; the
`2^{n/2-2}` guess fails at `n = 10`, HYP-4971). New here: in every case **each class
carries exactly one blue self-line**, its two tilings one closed and one open; every
class is strong, near-regular (scores `(n/2-1)^{n/2} (n/2)^{n/2}`), and almost all
asymmetric: `|Aut| = 1` except one class each at `n = 6, 8, 10` with `|Aut| = 3, 3, 5`:

```text
n = 4 : (H,|Aut|,HC,t,s,c) = (5,1,1,5,3,4)                                     1 class  (S4)
n = 6 : (41,1,4,41,9,24), (45,3,5,15,9,10)                                      2 classes
n = 8 : (637,1,47,637,39,376), (657,1,52,657,37,416), (653,1,52,653,37,416), (621,3,45,207,39,120)
n = 10: 24 classes, H in [15201, 15565], HC in [897, 936], s in {179,...,189}, one with |Aut| = 5
        ((15565, 5, 936, 3113, 187, 1872)), the other 23 with |Aut| = 1.
```

## 7. Terminology fixed by this theorem

**closed / open tiling** (the `(n,1)` tile), **closure number** `cl(C) = n HC/|Aut|`,
**fully open node** (= non-strong), **fully closed node** (`H = n HC`), **mixed-closure
node** (both), **first cut** (Proposition F), **convention-free quotient** (Proposition K).
Owner's dictionary: point of symmetry = pure-blue node; for `n >= 6` a point of symmetry
is a fully open leaf of the closed world.

## 8. Files

`04-computation/experiments/metagraph_20261020/closed_open_laws.py` (`_n7.out`, `_n8.out`),
`blue_self_lines_classes.py` / `.out`, `merged_metagraph.py` (`_n7.out`, `_n8.out`:
per-node tables with `(H, |Aut|)`, symmetric counts, incidences, multiplicities),
`staircase_symmetries.py`, `half_flips.py` (negative results: the staircase transpose,
rotation and half flips carry no tournament law), `blue_self_lines.py`,
`blue_self_lines_sym.py` (counts `1, 0, 2, 0, 4, 0, 24` for `n = 4..10`).
