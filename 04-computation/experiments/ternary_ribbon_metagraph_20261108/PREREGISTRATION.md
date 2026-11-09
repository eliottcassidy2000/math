# Pre-registration: base-3 ribbon transfer census and metagraph lines between node classes

mac-mini-2026-11-08-ternary. Written and pushed before any census program was run. Concurrent: opus-2026-11-09-S87
runs the base-3 order-10 strong census and the lane reduction (its PREREGISTRATION.md, Q1-Q3); nothing here
duplicates it. The two readings of "ternary Chan-Pak" in the repository are both addressed: (a) base 3, the values of
H_3(T) = I(Omega(T), 3) (THM-4940, HYP-9530 (A)); (b) the ternary return forest h -> 3h + {0, 2, 4} of
triple_reentry_interfaces_20261103 (readiness of near-join digits).

## Objects

- Ribbon words of span w (THM-4880): back arcs of span <= w on an ordered vertex set; letters = subsets of {1..w}.
- Odd-cycle packing polynomial I(Omega(T), x) = sum_k alpha_k(T) x^k; H_b = I(Omega, b); H_2 = H (OCF).
- Metagraph of order n (01-canon/definitions.md, THM-643/644/4729): base path n -> ... -> 1, m = C(n-1, 2) tiles,
  tilings t in {0,1}^m, nodes = isomorphism classes (merged with converse), lines {t, t^c} (all tiles flipped),
  blue iff t is grid-symmetric. Closed tiling: tile (n,1) carries 1 -> n.
- New node classifications: (K) number k of strong components (1..n); (V) eternal colour of the value H (THM-4930:
  eternally blue / eternally black); (F) first creator or not (THM-4922); (S) self-converse or not (THM-643).

## Predictions (to be scored, not edited)

- **P1 (finite-state).** I(Omega(T_u), x) is a finite-state series over ribbon words: a frontier automaton whose
  states are multisets of fragments (head age < w, tail age or free, length parity) with closed odd cycles weighted
  x; evaluating at x = 2 reproduces THM-4880's path count on every tested word (a consequence of OCF, so a test of
  both automata). Predicted reachable state counts: span 1 at most 10, span 2 at most 60, span 3 at most 600.
- **P2 (base-3 misses, span 2).** The values H_3 of span-2 ribbon words below 10^5 miss 10 and 31 and a finite set of
  further values, all below 2,000 (no miss in [2,000, 10^5]).
- **P3 (base-3 misses, span 3).** The span-3 misses below 10^4 are 10, 31 and at most three further values.
- **P4 (strong = all).** As in base 2, every attained value is attained by a strong word of the same length bound
  from length 4 on.
- **P5 (component formula, theorem to be proved before the census is read).** For a line {t, t^c} the cut sets of
  T_t and T_{t^c} are disjoint, so k(T_t) + k(T_{t^c}) <= n + 1 with equality only on the transitive line; the closed
  end is strong; the number of lines whose open end has cut set exactly S is N(S) = sum_{S' superset S}
  (-1)^{|S' - S|} 2^{e(S')}, e(S') = sum over the blocks B of S' of C(|B| - 1, 2); for blue lines the same with
  2^{(e(S') + f(S'))/2} and only sigma-symmetric S', f = number of sigma-fixed tiles inside blocks. The census at
  n = 7, 8, 9 will match these formulas exactly.
- **P6 (strong-strong share).** The share of lines with both ends strong exceeds 0.9 at n = 9 and increases with n.
- **P7 (value colours).** At every order 6 <= n <= 9 some line joins two nodes whose values are both eternally
  black; the share of lines with both ends eternally blue exceeds 0.5 at n = 9.
- **P8 (multiplicity).** The largest number of lines between two distinct nodes at order n is at most
  min(H/|Aut|) of the two nodes and is attained by a pair of self-converse nodes for n <= 8.
- **P9 (readiness, return forest).** For every odd h <= 99 attained at order <= 8, some tournament of order <= 8 with
  H = h has a vertex t with R_t = h + 1, and one with R_t = h + 2 (so 3h + 2 and 3h + 4 are near-join children).
  Expected failures: none for h >= 9; h in {1, 3, 5} may fail for R = h + 2.

A failure of P2/P3 says nothing about HYP-9530 (A): it measures the span-w family only.

Outcome: to be appended below without editing the predictions.

## Outcome (appended after the runs; the predictions above are unedited)

- **P1 CONFIRMED.** Reachable packing states 4, 20, 112 (predicted at most 10, 60, 600); 360 random words agree with
  the brute force and, at x = 2, with THM-4880's path automaton.
- **P2 CONFIRMED to 10^5.** The span-2 misses are exactly 10 and 31 (complete for all lengths to 10^4 by the fixed
  point; to 10^5 by length 15, state sets capped beyond): no further value at all (the predicted "finite set below
  2,000" is empty).
- **P3 CONFIRMED to 5,000** (not 10^4: 70 reduced coordinates, the run was capped at 5,000,000 states at length
  11): the span-3 misses are exactly 10 and 31, attained by length 10.
- **P4 REFUTED at one value.** Strong span-2 words attain every attained value except 64, which no strong tournament
  of any order attains (THM-4952 (E)); otherwise strong = all.
- **P5 CONFIRMED at n = 5, 6, 7, 8, 9** (formula = census, cut set by cut set, all and blue lines; at n = 9 over
  134,217,728 lines).
- **P6 CONFIRMED by theorem and by census to n = 8**: strong-strong shares 0.5625, 0.7637, 0.8782, 0.9383, and the
  exact formula gives 1 - 2^(4-n)(1 + O(2^-n)), above 0.9 from n = 8 on; n = 9 census: 0.9690.
- **P7 CONFIRMED at n = 6, 7, 8, half-decidable at n = 9.** Black-black lines 20, 702, 35,888, 115,704; blue-blue
  shares 0.617, 0.632, 0.672 at n = 6..8; at n = 9 the blue-blue share is 0.814 among the 9% of lines whose two
  values are below 1131 and undecidable overall (THM-4930's exact colours stop at 1131).
- **P8 half CONFIRMED, half REFUTED.** The multiplicity bound min(H/|Aut|) holds on every pair; the maximum is attained
  by a self-converse pair at n = 5 and 7 but by an SC-NS pair at n = 6 (41, 43) and by a converse pair (both
  non-self-converse, H = 445) at n = 8.
- **P9 REFUTED.** Near-join readiness fails for 47 of the 97 attainable odd h <= 199 at digit 2 and for 42 at digit 4
  (orders <= 8), including many h >= 9 (15, 25, 27, 35, 41, 43, ...); the expected failures at h <= 5 are a small part.
