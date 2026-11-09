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
