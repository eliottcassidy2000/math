---
id: THM-4953
title: "Lines of the order-n tiling metagraph between node classes (THM-643/4729 frame). (A) The strong components of each end of a line {t, t^c} are intervals of the base path and the two cut sets are disjoint, so k(T_t) + k(T_{t^c}) <= n + 1 with equality only on the transitive line; non-cut positions come in runs of >= 2, so no end has exactly n - 1 strong components. (B) The closed end is strong and the open end has cut set exactly S for N(S) = sum_{S' superset S} (-1)^|S' - S| N_sup(S') lines, N_sup(empty) = 2^(m-1), N_sup(S') = 2^(sum over blocks C(|B|-1, 2)); blue lines: the same inversion with 2^((m+f)/2 - 1) and 2^((e + f)/2) on the symmetrised cut set. Hence the share of lines whose open end is not strong is 2^(4-n) (1 + O(2^-n)), and lines between two non-strong nodes do not exist. (C) Both ends of a blue line are self-converse, so blue lines join only eternally-blue-valued nodes (THM-4930) and never touch an eternally-black value; the number of black lines between the self-converse nodes and the others has the parity of the number of non-self-converse nodes. (D) Value fibres: the transitive node's single line goes to the sole-exit tower (tribonacci count, THM-4865); for h < f(n) the fibre of h carries exactly sum_{H(C) = h} h/|Aut C| lines, so fibre 3 carries n - 2 lines (n >= 4) and fibre 5 carries 5(n - 3) (n >= 5). FINITE-EXACT census of lines by component, value-colour, first-creator and self-converse classes at n = 5..9."
status: >
  PROVED (mac-mini-2026-11-08-ternary, elementary): (A), (B) by the interval lemma for tournaments with a Hamiltonian
  path, the crossing description of cuts, and Moebius inversion on the Boolean lattice of cut sets; the formulas were
  checked against the exhaustive line census (lines_census.c, every closed tiling) at n = 5, 6, 7, 8, cut set by cut set, all lines and blue lines, and at n = 9 (134,217,728 lines in eight parts). (C) from the definitions and THM-4729 (both ends of a blue line
  are blue tilings), THM-4930 (self-converse value => eternally blue) and THM-643's per-node parities. (D) from
  THM-4729's closure (closed => strong => H >= f(n)) and the classification of the realisers of 3 and 5 (a 3-cycle,
  resp. the strong 4-tournament, with transitive padding; THM-1370 / THM-4744 for the strong floors), with the
  transitive partner identified by its arc rule. FINITE-EXACT: the node-level census with pynauty canonical forms at
  n = 5..8 (lines_nodes.py; THM-643's fibre sizes H/|Aut| and parities re-verified), the (H closed, H open) pair
  statistics (pairs_analysis.py). KNOWN/CITED: the tiling frame and the blue/black lines (THM-643, 644), closed and
  open tilings (THM-4729), SC = A000570 (THM-793), the tribonacci tower (THM-4865 (c), A000213), eternal colours
  (THM-4930), first creators (THM-4922). Concurrent: opus-2026-11-09-S87's blue-line H mod 4 law at order 11
  (PREREGISTRATION-2 of ternary_metagraph_20261109), independent of this. AUDIT: PASS WITH CORRECTIONS
  (05-knowledge/results/ternary_ribbon_metagraph_audit_20261108.md; minors here: THM-4729 (3)-(4) cited for the
  strong-strong count and the absence of lines between non-strong nodes, one transcription error in the census
  paragraph), applied on 2026-11-09; the audit re-derived (A), (B) and verified the formula by its own enumeration
  at n = 5, 6, 7 and the node census at n = 5, 6.
session: mac-mini-2026-11-08-ternary
source: 05-knowledge/results/ternary_ribbon_metagraph_lines_20261108.md
scripts:
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/lines_census.c (lines_n{5..9}_p*.txt)
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/lines_formula.py (formula_n{5..9}.txt)
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/lines_nodes.py (lines_nodes_n{4..8}.json)
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/lines_analysis.py (lines_analysis.txt)
  - 04-computation/experiments/ternary_ribbon_metagraph_20261108/pairs_analysis.py
related:
  - THM-643
  - THM-644
  - THM-4729
  - THM-4730
  - THM-785
  - THM-793
  - THM-4865
  - THM-4922
  - THM-4930
  - HYP-9520
  - HYP-9541
---

# THM-4953 - metagraph lines between node classes

Frame (01-canon/definitions.md, THM-643, THM-4729): order n, base path n -> n-1 -> ... -> 1, m = C(n-1, 2) tiles
(a, b) with a >= b + 2, bit 0: a -> b, bit 1: b -> a; nodes = isomorphism classes (merged with the converse in the
merged metagraph); lines {t, t^c}, blue iff t is grid-symmetric; closed iff tile (n, 1) carries 1 -> n; every line
has one closed end. Position i is the gap between path vertices i + 1 and i; a tile (a, b) crosses i iff b <= i < a.

**(A) Intervals, disjoint cuts, the missing component count.** The condensation of a tournament is a transitive
order, and a Hamiltonian path cannot return to an earlier component, so the strong components of T_t are intervals of
the base path; position i is a cut of T_t iff every crossing tile is 0 and a cut of T_{t^c} iff every crossing tile is
1. Every position is crossed (n >= 3), so the cut sets are disjoint and k(T_t) + k(T_{t^c}) <= n + 1; equality needs
every position cut on one side, and the tile (i + 2, i) crossing both i and i + 1 forces all tiles equal: the
transitive line {0^m, 1^m}. Every tile crosses >= 2 consecutive positions, so non-cut positions form runs of length
>= 2 and no end has exactly n - 1 components.

**(B) Exact counts by cut set.** The closed end has tile (n, 1) = 1, which crosses every position, so it is strong.
The open end's cut set contains S iff every tile crossing S is 1; the free tiles are those inside the blocks of S, so
N_sup(S) = 2^(e(S)), e(S) = sum over blocks B of C(|B| - 1, 2), for S nonempty, and 2^(m-1) for S empty. The number of
lines whose open end has cut set exactly S is N(S) = sum_{S' superset S} (-1)^|S' - S| N_sup(S'), and the number of
lines between the strong class and the class with k strong components is L_k = sum_{|S| = k - 1} N(S). That the
closed end is strong, that no line joins two non-strong nodes and that the strong-strong count is 2^(m-1) - NS(n)
(NS = 14, 121, 1,995, 64,648 at n = 5..8) are THM-4729 (3)-(4), KNOWN; the new content is the full cut-set
distribution N(S) and its blue refinement. Blue lines: a
grid-symmetric tiling is free on the orbit representatives of its free tiles; with sigma(i) = n - i and Sbar = S u
sigma(S), N_sup^blue(S) = 2^((m+f)/2 - 1) for S empty and 2^((e(Sbar) + f(Sbar))/2) otherwise (f(Sbar) = number of
sigma-fixed tiles, a + b = n + 1, inside the blocks), and the same inversion gives the blue counts (zero unless S is
symmetric). Values (all lines; blue in brackets), open end with k components:

| n | k = 1 | 2 | 3 | 4 | 5 | 6 | n |
|---|---|---|---|---|---|---|---|
| 5 | 18 (6) | 10 (0) | 3 (1) | 0 | 1 (1) | | |
| 6 | 391 (27) | 101 (1) | 15 (3) | 4 (0) | 0 | 1 (1) | |
| 7 | 14,389 (239) | 1,816 (0) | 153 (15) | 20 (0) | 5 (1) | 0 | 1 (1) |
| 8 | 983,928 (1,978) | 61,671 (5) | 2,739 (59) | 206 (2) | 25 (3) | 6 (0) | 1 (1) at k = 8 |
| 9 | 130,053,749 (32,249) | 4,067,314 (0) | 92,695 (501) | 3,672 (0) | 260 (16) | 30 (0) | 7 (1), 0, 1 (1) at k = 7, 8, 9 |

The leading terms of 2^(m-1) - N(empty) are the single cuts at positions 1 and n - 1 (an open end with a source or
a sink), 2^(C(n-2, 2)) = 2^(m-1) 2^(3-n) each, the next 2^(m-1) 2^(6-2n); so the share of lines whose open end is not
strong is 2^(4-n) (1 + O(2^-n)) (observed non-strong shares 0.4375, 0.2363, 0.1218, 0.0617, 0.0310 at n = 5..9
against 2^(4-n) = 0.5, 0.25, 0.125, 0.0625, 0.03125), and the blue k = 2 count vanishes at odd n (no symmetric single cut).

**(C) Colours.** A blue tiling is grid-symmetric, hence its tournament is self-converse and its value eternally blue
(THM-4930); both ends of a blue line are blue tilings (THM-4729). So every blue line joins two eternally-blue-valued
self-converse nodes and no blue line touches an eternally-black value. By THM-643 (blue odd and black even at
self-converse nodes, black odd elsewhere), the black lines between the self-converse class and its complement have
the parity of the number of non-self-converse nodes: 12, 156, 5,044, 53,072 lines against 4, 44, 368, 6,704 nodes at
n = 5..8.

**(D) Fibres.** The transitive node has one tiling (0^m) and one line, to 1^m: the tournament with i -> j iff
j = i - 1 or j >= i + 2, the sole-exit tower of THM-4865 with H = 1, 3, 5, 9, 17, 31, 57, 105, 193, 355, 653, 1201
(n = 2..13, A000213). For h < f(n) no h-tiling is closed (closed => strong => H >= f(n)), so the fibre of h carries
exactly sum_{H(C) = h, |C| = n} h/|Aut C| lines. The tournaments with H = 3 are TT_a => C3 => TT_b (|Aut| = 3,
n - 2 classes): fibre 3 carries n - 2 lines for n >= 4; those with H = 5 are TT_a => S4 => TT_b with S4 the strong
4-tournament (|Aut| = 1, n - 3 classes): fibre 5 carries 5(n - 3) lines for n >= 5 (census 3, 4, 5, 6 and 10, 15, 20,
25 at n = 5..8).

**Census by node class (FINITE-EXACT).** Iso-class nodes 12, 56, 456, 6,880 (self-converse 8, 12, 88, 176) at
n = 5..8. Self-converse: SC-SC blue 8, 32, 256, 2,048; SC-SC black 12, 10, 876, 772; SC-NS black 12, 156, 5,044,
53,072; NS-NS 0, 314, 10,208, 992,684. Value colours (nodes blue/black 12/0, 48/8, 370/86, 5,556/1,324): lines with
both values eternally blue 32, 316, 10,356, 704,608 (shares 1, 0.617, 0.632, 0.672), mixed 0, 176, 5,326, 308,080,
both eternally black 0, 20, 702, 35,888. First creators (THM-4922; nodes 6, 33, 328, 5,457): lines FC-FC 18, 371,
13,723, 932,047, FC-old 14, 131, 2,472, 104,099, old-old 0, 10, 189, 12,430. Multiplicities: the largest number of
lines between two distinct iso-class nodes is 2, 6, 8, 12, never above min(H/|Aut|) of the pair; self-lines 4, 8, 44,
206; same-value lines 7, 53, 470, 8,179; the closed end has the larger value in 0.75, 0.73, 0.72, 0.70 of the lines;
at n = 8 the correlation of the end values is 0.741. At n = 9 (134,217,728 lines, 32,768 blue; values only):
fibres 3 and 5 carry 7 and 30 lines, same-value lines 267,735 (0.0020), the closed end is larger in 0.686 of the
lines, the correlation is 0.800, the largest value is 3,357; among the 12,188,691 lines whose two values are both
below 1131 (decided colours), 9,921,395 join two eternally blue values, 115,704 two eternally black ones and
2,151,592 one of each (the other 91% of the lines carry a value >= 1131). Merged nodes (converse pairs identified): 10, 34, 272, 3,528 nodes; largest multiplicity between distinct
merged nodes 2, 14, 16, 26 at the value pairs (9, 9), (37, 43), (151, 153), (659, 653); largest blue multiplicity
1, 2, 2, 4; self-lines 4, 26, 114, 622 (blue 0, 2, 0, 4).
