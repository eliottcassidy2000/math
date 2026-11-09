---
id: HYP-9541
title: "Metagraph lines and ternary Chan-Pak (after THM-4952, THM-4953). (1) The share of lines of the order-n metagraph joining two eternally-black-valued nodes tends to 0 and the share joining two eternally-blue-valued nodes tends to 1 (the value-colour form of HYP-9520); blue lines are a vanishing share 2^((m+f)/2 - m) of all lines and all of them join eternally-blue values. (2) The largest number of lines between two distinct merged nodes grows without bound (2, 14, 16, 26 at n = 5..8; 2, 6, 8, 12 for iso-class nodes) but subexponentially, and is attained at pairs of large, nearly equal values; the largest blue multiplicity (1, 2, 2, 4) grows too. (3) The closed end carries the larger value in a share of lines that stays in (1/2, 1) (0.75, 0.73, 0.72, 0.70 at n = 5..8) and the correlation of the two end values converges to a limit in (0, 1). (4) Ternary: the strong base-3 spectrum is Spec_3 minus {1, 64}; every v = 1 (mod 3) other than 10 and 31 is H_3 of a span-3 (indeed span-2, except for the finitely many values whose strong realisers need width 3) ribbon word, which implies HYP-9530 (A). (5) For every span w the packing series' dimension at the Hamiltonian-path base x = 2 equals the nonempty-word rank of THM-4880's path series (13 for w = 2, 61 for w = 3) and lies strictly below the generic dimension (14, resp. >= 70): the path base is always a rank drop, by one for span 2 and by at least nine for span 3."
status: >
  OPEN (mac-mini-2026-11-08-ternary). (1) NUMERICAL from the n = 5..8 node census (black-black shares 0, 0.039,
  0.043, 0.034; blue-blue 1, 0.617, 0.632, 0.672; the eternally-black node share rises 0, 0.14, 0.19, 0.19, so the
  conjecture rests on HYP-9520's thinning, not on monotone data). (2) NUMERICAL (iso-class nodes 2, 6, 8, 12; merged nodes 2, 14, 16, 26 at adjacent
  large values (37, 43), (151, 153), (659, 653); blue 1, 2, 2, 4). (3) NUMERICAL (four orders). (4) The strong miss 64 is PROVED (THM-4952 (E)); the
  completeness clauses are FINITE-EXACT to 10^4 (span 2, all lengths) and PENDING-S3 (span 3); the "finitely many
  width-3 values" clause is a guess from base 2, where only 47 needed span 3 (THM-4880 (E)). (5) PROVED for w = 2
  (THM-4952 (B)); for w = 3 the exact ranks at x = 3 and x = 2 are 70 and 61 (the generic rank and the other drop
  points for w = 3 are not computed). Chan-Pak and ternary Chan-Pak
  remain OPEN (MISTAKE-559 applies).
session: mac-mini-2026-11-08-ternary
date: 2026-11-08
source: 05-knowledge/results/ternary_ribbon_metagraph_lines_20261108.md
related: [THM-4952, THM-4953, THM-4880, THM-4940, HYP-9530, HYP-9520, THM-4930, THM-643, THM-4729]
scripts: 04-computation/experiments/ternary_ribbon_metagraph_20261108/lines_analysis.py; pairs_analysis.py; ribbon_levels_mod.c; packing_intrep.py
---

# HYP-9541 - metagraph line shares and ternary ribbon completeness

**(1) Value colours.** THM-4953 (C) proves that blue lines join only eternally-blue values. The conjecture is about
black lines: as n grows, almost every line joins two nodes whose values are eternally blue, and the lines joining two
eternally-black values become negligible. The node census gives black-black shares 0, 0.039, 0.043, 0.034 and
blue-blue shares 1, 0.617, 0.632, 0.672 at n = 5..8, while the share of eternally-black-valued nodes is 0, 0.14, 0.19,
0.19: the conjecture is HYP-9520 (the eternally black values have density zero) transported to the metagraph, and a
proof needs the distribution of H over the tilings of order n (THM-643: the fibre of h has sum h/|Aut| tilings) to
concentrate on values where black has density zero.

**(2) Multiplicities.** The number of lines between two nodes is at most min(H/|Aut|) of the two (each tiling lies on
one line). The observed maximum between distinct iso-class nodes is 2, 6, 8, 12 and between distinct merged nodes 2, 14,
16, 26 (at (37, 43), (151, 153), (659, 653): large, nearly equal values); the conjecture is unbounded but
subexponential growth at such pairs (at n = 8 the iso-class maximum 12 is between the two classes of one converse pair
with H = 445, a self-line of the merged metagraph).

**(3) Direction of the complement.** The closed end (a strong tournament with a Hamiltonian cycle through the base
path) has the larger count in about seven lines out of ten; the complement flips every non-path arc, i.e. passes to
the converse and reverses one Hamiltonian path, and the conjecture is that this operation lowers H more often than it
raises it in a limiting proportion, with the two counts positively correlated in the limit (0.741 at n = 8).

**(4) Ternary.** 64 is the only value 1 (mod 3) attained by a tournament but by no strong tournament (proved for 64;
conjectured to be the only one). Ternary ribbon completeness: every v = 1 (mod 3) other than 10 and 31 is H_3 of a
span-3 ribbon word; by THM-4940 (1) (multiplicativity over order joins) and THM-4880 (C) (ribbon closure) this
reduces to the strong span-3 words and would prove HYP-9530 (A).

**(5) Rank drops.** For span 2 the packing series drops in dimension exactly at x = 0 and x = +-2 (THM-4952 (B)),
by one at x = 2. For span 3 the dimension is 70 at x = 3 and 61 at x = 2 (THM-4880's path rank 62 minus the empty
word). Conjecture: at every span the value at x = 2 is the nonempty-word path rank and is a strict drop, growing with
w; the extra relations at x = 2 are the OCF identities between packings and Hamiltonian paths that have no analogue
at other bases.
