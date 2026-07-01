---
id: HYP-3808
title: THE BLUE/BLACK LINE STRUCTURE OF THE MERGED METAGRAPH IS A PARITY DECOMPOSITION -- the BLACK lines form an EVEN-degree (Eulerian) graph and the BLUE lines an all-ODD-degree graph, and tiling-count PARITY = the SC/NS type (SC merged node <=> ODD tiling count, NS <=> EVEN). All the owner's parity claims verified (n=4,5,6): pure-black nodes have EVEN black_other; mixed have EVEN black_other + ODD blue_other; pure-blue have ODD blue_other. Closed forms: total lines=2^(m-1), blue lines=2^((m+floor((n-1)/2))/2 - 1), m=C(n-1,2). Eligibility: BLUE lines only MIXED-MIXED / MIXED-PUREBLUE / MIXED-self (PURE_BLUE is always PENDANT, never blue-blue, never self); BLACK lines MIXED-MIXED / MIXED-PUREBLACK / MIXED-self at n<=5, PLUS PURE_BLACK-PURE_BLACK (sea) + PURE_BLACK-self at n>=6. CORRECTION to the owner's conjecture: self-loops are NOT only on mixed -- pure-black (NS) self-loops appear at n=6; the layered pure-black--mixed--pure-blue model is an n<=5 phenomenon.
status: CONFIRMED (exhaustive n=4,5,6 on the tournament-tiling-explorer merged metagraph). Verifies the owner's three parity claims; establishes the parity theorem (SC odd / NS even tiling count) and the even-graph(black)/odd-graph(blue) line decomposition; gives closed-form line counts; REFUTES the "self-loops only on mixed" conjecture at n=6 and locates the NS-NS sea + pure-black self-loop onset at n=6. A structural characterization of the pairing/assignment process, not yet a proof of all parts for general n.
source: mac-mini-2026-07-01-S83
related:
  - HYP-3799   # S82 even graphs -- the BLACK line subgraph IS an even-degree graph (this ties the two)
  - HYP-3798   # S81 arc-hypercube invariants (same tiling model)
results:
  - 04-computation/merged_metagraph_blue_black_lines_macmini_20260701.py
  - 05-knowledge/results/merged_metagraph_blue_black_lines_macmini_20260701.out
---

# HYP-3808 -- the merged metagraph's blue/black lines are a parity decomposition

Ground truth (tournament-tiling-explorer.html): a **line** = `{tiling t, flip(t)}` (`flip` = complement-
tiling, flip all `m=C(n-1,2)` tiles); **blue** iff `isGridSym(t)` (anti-diagonal reflection), **black**
else; `2^m` tilings pair into `2^{m-1}` lines. MERGE by transpose (SC = transpose-self; NS classes pair).
Each line connects `merged(class(t))` to `merged(class(flip(t)))`; a node's **tiling count** = its number
of tilings = its degree in the line-multigraph (`= blue_other + black_other + 2*(blue_self+black_self)`).

## The owner's three parity claims -- ALL VERIFIED (n=4,5,6)
- **PURE_BLACK** (merged NS) nodes: `black_other` is EVEN.
- **MIXED** (SC with both) nodes: `black_other` EVEN and `blue_other` ODD (odd total).
- **PURE_BLUE** (SC all grid-sym) nodes: `blue_other` ODD, only blue.

## The parity theorem (new, clean)
`tiling_count(node)` parity `= blue_other` parity (black_other and selfs are even). So:
**SC merged node `<=>` ODD tiling count; NS merged node `<=>` EVEN tiling count.**
Reason: an NS-merged node is a pair `{A, A^op}` with `|tilings(A)| = |tilings(A^op)|`, so count `= 2|A|`
(even); an SC node is a single self-complementary class, and it has an ODD number of grid-symmetric tilings
(hence odd `blue_other`, hence odd count). Corollary (blue-subgraph handshake): **the number of SC merged
nodes is EVEN** (verified: 2, 8, 12 for n=4,5,6).

## The even/odd decomposition of the lines (the deep frame)
- The **BLACK line subgraph** has EVEN degree at every node (`black_other` even everywhere) -- it is an
  **even-degree (Eulerian) graph**, decomposing into edge-disjoint cycles. *This IS an "even graph" in the
  cycle-space / A002854 sense (S82 notion I) -- the merged metagraph's black lines realize one.*
- The **BLUE line subgraph** has ODD degree at every SC node -- an **all-odd-degree graph** on the SC nodes
  (an odd/T-join-type object), forcing `#SC` even.
- So the line-multigraph splits as `EVEN(black) (+) ODD(blue)`, and `tiling_count = odd_blue(@SC) +
  even_black + 2*self`. The even/odd theme the project chases sits literally in the metagraph's own lines.

## Eligibility rules for the pairing process (refined + corrected)
The "process" = realize the prescribed tiling counts (the H-spectrum) by adding colored lines; each line
`+1` to two nodes (or `+2` to one, a self-loop). Precise eligibility (verified n<=6):
- **BLUE line** endpoints are both grid-sym tilings `=>` both SC `=>` allowed pairs are `MIXED-MIXED`,
  `MIXED-PURE_BLUE`, and `MIXED` self-loops. **NEVER `PURE_BLUE-PURE_BLUE`, never `PURE_BLUE` self-loop:
  every pure-blue node is a PENDANT hanging off the mixed core.**
- **BLACK line** endpoints are non-grid-sym `=>` both in `{NS, mixed}` `=>` `MIXED-MIXED`,
  `MIXED-PURE_BLACK`, `MIXED` self-loops; **and at `n>=6` also `PURE_BLACK-PURE_BLACK` (the NS-NS sea) and
  `PURE_BLACK` self-loops.**
- **CORRECTION to the owner's conjecture:** self-loops are NOT confined to mixed nodes. At `n=6`, pure-black
  (NS) nodes acquire black self-loops (24 of them) and a large NS-NS black "sea" (290 lines). The layered
  `pure-black -- mixed -- pure-blue` picture (and self-loops-only-on-mixed) is exact ONLY for `n<=5`; the
  sea and pure-black self-loops switch ON at `n=6` -- a genuine structural transition.

## Metrics (illuminating the process)
1. **Tiling-count parity** = the SC/NS 2-coloring (odd=SC, even=NS).
2. **`blue_other` = the "SC charge"** (odd for SC, 0 for NS); **`black_other` = even everywhere**.
3. **Closed forms**: `#lines = 2^{m-1}`; `#grid-sym tilings = 2^{(m+floor((n-1)/2))/2}`; `#blue lines =`
   half that `= 2, 8, 32, 256, 2048` (n=4..8); `#black lines = 2, 24, 480, 16128, 1046528`; blue fraction
   `-> 0`. Recursion `#lines(n)/#lines(n-1) = 2^{n-2}`.
4. **Pendant count** = #pure-blue leaves (`blue_other=1`).
5. **Interface degree** of a mixed node `= (blue_other, black_other) = (odd, even)`.
6. **NS-NS sea fraction** of black lines = `0` (n<=5), `290/480` (n=6) -- the transition order parameter.
7. **Self-loop census** by category: mixed (blue+black); pure-black (black, n>=6 only).
8. **Category counts**: `#PURE_BLUE, #PURE_BLACK, #MIXED` (n=6: 2, 22, 10); `#SC = #PURE_BLUE+#MIXED` even.

## Reframings
- The assignment is a **colored degree-realization / `f`-factor** on a category-typed host with parity
  constraints: BLACK = an Eulerian (even) factor on `{NS, mixed}`; BLUE = an odd factor on the SC nodes.
- Realizing the black lines = an **Eulerian decomposition** (cycles) of the even black graph; realizing blue
  = a **T-join / odd-degree** construction with pure-blue as leaves.

## Recursive structure + concrete next target
- Line counts recurse cleanly (`2^{n-2}`; grid-sym doubling). Mode-A (`n->n-1`) / Mode-B (`n->n-2`) tiling-
  count recursions should induce category recursions -- **open: how do PURE_BLUE/PURE_BLACK/MIXED counts and
  the sea-onset recurse?**
- **CONCRETE NEXT TARGET**: (a) PROVE the parity theorem (SC odd / NS even) via `|A|=|A^op|` + the odd-grid-
  sym-count lemma; (b) PROVE the eligibility rules (blue `=>` SC-SC, black `=>` non-blue) as theorems;
  (c) EXPLAIN the `n=6` onset of the NS-NS sea + pure-black self-loops (why `n>=6`?); (d) characterize the
  black Eulerian graph's cycle structure and whether it recovers the H-spectrum degrees.

## Honest scope
All parity claims, the parity theorem, the even/odd line decomposition, the eligibility rules, the closed
forms, and the `n=6` transition are VERIFIED exhaustively for `n=4,5,6`. The conjecture "self-loops only on
mixed" is REFUTED at `n=6`. General-`n` proofs of the parity theorem and eligibility (and the sea-onset
mechanism) are the open next targets.
