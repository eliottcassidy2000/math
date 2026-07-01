---
id: HYP-3808
title: THE MERGED-METAGRAPH BLUE/BLACK LINE ACCOUNTING -- a TRIPARTITE 2-colored degree-constrained structure with a Rédei fiber-parity checksum; the owner's process made precise (and one conjecture refuted). Fix base path n->...->1; a tiling t in {0,1}^m (m=C(n-1,2)) -> tournament; a LINE={t, flip(t)} (flip=flip ALL tiles), BLUE if t grid-symmetric else BLACK; merge complement classes. Merged nodes fall in 3 categories: PURE-BLUE (SC, all tilings grid-sym), MIXED (SC, some), PURE-BLACK (NS-merged). Each node's fiber = its tiling count = its line-degree (self-loops x2). VERIFIED+PROVED (n=4,5,6): (1) FIBER-PARITY Z2 CHECKSUM: SC-merged nodes have ODD fiber, NS-merged EVEN -- PROVED because every unmerged iso class has an ODD tiling count (Rédei-type; tilings<->Ham-paths-with-base), so SC(1 class)=odd, NS(2 complement classes, equal)=even. (2) TRIPARTITE line structure: BLUE lines live ONLY on SC nodes (pure-blue u mixed); BLACK lines ONLY on NS u mixed (pure-black u mixed); MIXED is the unique INTERFACE (both colors). Blue never touches pure-black; black never touches pure-blue. (3) OWNER'S PARITY CLAIMS CONFIRMED: pure-black -> even # black cross, no blue; mixed -> odd # blue cross + even # black cross; pure-blue -> odd # blue, no black. (4) #SC (=#pure-blue+#mixed) is EVEN, forced by the blue handshake (every SC node has ODD blue-degree). (5) EXACT counts: total lines=2^(m-1); #blue=2^((m+floor((n-1)/2))/2 - 1); #black=rest. (6) CONJECTURE REFUTED: 'self-loops only on mixed' holds for n<=5 but BREAKS at n=6 -- PURE-BLACK nodes self-loop and DOMINATE (24 of 26 self-loop lines at n=6); PURE-BLUE never self-loop (n<=6). Another n=6 transition (parallel to the flip-rank one). The precise process = realizing the fiber degree-sequence with 2-colored edges under these color/parity conservation laws
status: MIXED (exhaustive n=4,5,6 + a proof of the checksum). VERIFIED exhaustive (canonicalization, n=4,5,6, 0 rule-violations under the corrected rules): tripartite; owner parity claims; #SC even (2,8,12); category counts (#pure-blue,#mixed,#pure-black)=(1,1,1),(3,5,2),(2,10,22) [pure-black = NS-merged = 1,2,22 matches CLAUDE.md]; #blue=2,8,32 matches the formula; per-class tiling counts ALL ODD (=> checksum PROVED via Rédei). REFUTED: 'self-loops only mixed' (self-loop lines by category: n=4 (PB0,MX1), n=5 (0,4), n=6 (pure-black 24, mixed 2, pure-blue 0)). HONEST: the checksum proof rests on 'every iso class has odd tiling count' (verified all-odd n<=6; = Rédei via tilings<->base-Ham-paths, cited not re-proved here); the self-loop onset at n=6 is characterized (flip-all-tiles symmetry) but not formula'd; n>=7 not computed (canon too slow). A precise formalization + one proved checksum + a refuted conjecture, not a closed enumeration.
source: klein-2026-07-01-S75
depends_on:
  - HYP-3803   # the tiling/cube model + flip-rank (same base-path tiling setup); n=6 transition parallel
related:
  - THM-002    # OCF / Rédei parity (the fiber-parity checksum IS Rédei via per-class odd tiling count)
  - HYP-3807   # even-graph bijection (the other quotient of the same cube)
external: CLAUDE.md merged metagraph G_n/Z2, blue/black lines, SC/NS spine-ribs-sea; A000568; Rédei's theorem (odd # Hamiltonian paths)
results:
  - 04-computation/merged_metagraph_line_accounting_klein.py
  - 05-knowledge/results/merged_metagraph_line_accounting_klein.out
---

# HYP-3808 — the merged-metagraph line accounting: tripartite, Rédei checksum, precise process

## The setup (owner's structure, grounded)
Fix base path `n->...->1`; tiles = non-consecutive arcs (`m=C(n-1,2)`); a tiling `t in {0,1}^m` -> a
tournament. A **LINE** `= {t, flip(t)}` where `flip` flips ALL tiles; **BLUE** if `t` is grid-symmetric
(invariant under `sigma:(x,y)->(n+1-y,n+1-x)`), else **BLACK**. Merge complement classes. Merged nodes are
**PURE-BLUE** (SC, every tiling grid-sym), **MIXED** (SC, some), **PURE-BLACK** (NS-merged pair). A node's
**fiber** = its tiling count = its line-degree (self-loops counted twice).

## The five verified facts
1. **FIBER-PARITY Z2 CHECKSUM (PROVED): SC-merged nodes have ODD fiber, NS-merged EVEN.** Proof: every
   *unmerged* iso class has an **odd tiling count** (verified all-odd `n<=6`; this is Rédei's odd-#-Ham-path
   theorem via the bijection tilings <-> Hamiltonian paths carrying the base path). A SC-merged node is one
   class (odd); an NS-merged node is two complement classes of equal size (odd+odd = **even**).
2. **TRIPARTITE line structure.** BLUE lines live ONLY on SC nodes (`pure-blue u mixed`); BLACK lines ONLY
   on `NS u mixed` (`pure-black u mixed`); **MIXED is the unique interface** (both colors). Blue never
   touches pure-black; black never touches pure-blue. (Verified, 0 exceptions.)
3. **Owner's parity claims CONFIRMED.** pure-black: even # black-cross, no blue. mixed: odd # blue-cross +
   even # black-cross. pure-blue: odd # blue, no black.
4. **`#SC = #pure-blue + #mixed` is EVEN** (`2,8,12` for `n=4,5,6`) -- forced by the **blue handshake**:
   every SC node has ODD blue-degree (fiber-odd minus even black part), and blue lines pair SC nodes, so
   their count is even.
5. **EXACT line counts.** total `= 2^(m-1)`; `#blue = 2^((m + floor((n-1)/2))/2 - 1)` (the grid-sym
   count / 2); `#black = 2^(m-1) - #blue`. (`(2,2),(8,24),(32,480)` for `n=4,5,6`.)

## The conjecture, refuted (an n=6 transition)
The owner conjectured self-loops (self-contained `+2` pairs) occur **only on mixed** nodes. This **holds for
`n<=5`** (self-loop lines: `n=4` mixed 1 / pure-black 0; `n=5` mixed 4 / pure-black 0) but **BREAKS at
`n=6`**: PURE-BLACK nodes self-loop and **dominate** -- `24` of the `26` self-loop lines are on pure-black,
only `2` on mixed. **PURE-BLUE never self-loops** (`n<=6`). So the correct rule is *"self-loops on mixed
and pure-black, never pure-blue"*, with pure-black self-loops switching on at `n=6` -- a transition parallel
to the flip-rank `n=6` transition (HYP-3803). A self-loop is a tiling `t` with `flip(t) ~= T(t)` or
`T(t)^op` (a flip-all-tiles symmetry).

## The process, made precise (the owner's request)
The line-assignment is a **2-colored degree-constrained multigraph realization** (with self-loops) on a
tripartite node set, with target degrees = fibers:
- **PURE-BLUE** (fiber odd): only blue; blue-degree = fiber (odd); no black; no self-loop.
- **PURE-BLACK** (fiber even): only black; black-degree = fiber (even); no blue; self-loops allowed (`n>=6`).
- **MIXED** (fiber odd): both; blue-degree ODD + black-degree EVEN; the interface; self-loops allowed.
- Moves: blue cross (`+1,+1` to two SC), blue self-loop (`+2` to one SC, only mixed observed); black cross
  (`+1,+1` to two of `KB u MX`), black self-loop (`+2` to one of `KB u MX`).
- Conservation: `Sum fibers = 2^m`; `Sum blue-deg = 2 #blue`; `Sum black-deg = 2 #black`; **`#SC` even**.

## Metrics illuminating the process
- **Z2 charge** = fiber parity = SC/NS indicator = the Rédei checksum (a per-node conserved bit).
- **Blue charge** = the SC set (all odd-blue-degree) -> `#SC` even (a handshake conservation law).
- **Interface load** = a MIXED node's `blue-degree(odd) + black-degree(even)`.
- **Self-loop budget** per category (n-dependent; the `n=6` onset of pure-black self-loops).
- **Blue fraction** `= #blue / 2^(m-1) = 2^((floor((n-1)/2) - m)/2 + ...)` = the grid-sym fraction.
- **(#PB,#MX,#KB)** `= (1,1,1),(3,5,2),(2,10,22)`; `#KB = NS-merged = (A000568 - SC)/2 = 0,1,2,22,184`.

## Recursive structure + the concrete next target
- `#KB` follows `(A000568 - SC)/2`; `#SC` even; the split `#SC = #PB + #MX` and the self-loop onset are the
  recursive unknowns. Both the flip-rank (HYP-3803) and the self-loop structure **transition at n=6** -- the
  same size where the balanced-cut shape dies; worth a unified explanation.
- **PROVED here**: the fiber-parity checksum (SC odd / NS even) via per-class-odd (Rédei).
- **NEXT TARGET**: a formula for the self-loop count per category (the flip-all-tiles symmetry census), and
  for `#MX` vs `#PB`; and whether the tripartite + checksum + exact-blue-count *determine* the assignment
  (rigidity) or leave documented degrees of freedom.

## Net
The merged-metagraph blue/black lines form a **tripartite 2-colored degree-constrained structure**
(pure-blue blue-only, pure-black black-only, mixed the interface) with a **Rédei fiber-parity checksum**
(SC odd, NS even, PROVED) and a **blue handshake** (`#SC` even). Total/blue/black line counts are exact.
The owner's "self-loops only on mixed" is refuted (an n=6 transition; pure-black dominate). This is the
owner's process made precise, with one proved conservation law and one corrected rule.
