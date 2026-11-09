---
id: THM-4951
title: "Packing-number lines of the tiling metagraph. Let nu(T) be the maximum number of vertex-disjoint cycles of T and k = floor(n/3). On every line {t, t + 1}, nu(t) + nu(t + 1) >= k (Session26's consecutive-triples mechanism). (1) On a tight line (equality) both packing numbers are read off the triple word x_a = [tile (a, a + 2) reversed]: nu(t) = P(A) and nu(t + 1) = P(B), where P is the largest 3-separated subset. (2) The boundary is attained with every split: for each 0 <= j <= k there is a tight line with nu-pair {j, k - j} (prefix staircase), and for each j >= 1 with 6j + 4 <= n a BLUE tight line with nu-pair {2j, k - 2j} (double staircase). (3) At n = 3k + 2 the triple words of tight lines are exactly the 2^k block-monochromatic words on the aligned blocks {3b, 3b + 1, 3b + 2}, and each is realised (multi-staircase), with nu(t) = number of 1-blocks; the blue ones are the palindromic words, so their nu-pairs are the {j, k - j} with j even or k odd. Counts: tight colourings 2^(k-1), (k+1) 2^(k-2), (k^2 + 9k - 2) 2^(k-4) at n = 3k + 2, 3k + 1, 3k. (4) Exact census of all lines by nu-pair and colour for n <= 9, and of tight lines for n <= 12. (5) Blue lines by H mod 4 at odd n: (1,1) = (3,3) always (from THM-4922 (4)); the split is 1:2:1 at n = 5, 7, 9 but not at n = 11."
status: >
  PROVED: (1); (2); (3) for n = 3k + 2 (necessity, sufficiency, the blue pairs (3d)) and the count 2^(k-1);
  (5)(a), (b), (d).
  FINITE-EXACT:
  - the tight-colouring counts at n = 3k + 1 and 3k (DP for k <= 13; the closed forms are fitted, not proved);
  - (4) and (5)(c).
  Pre-registered and tested:
  - PREREGISTRATION-2: R1 CONFIRMED (and KNOWN, see below); R2, R3 REFUTED at n = 11.
  - PREREGISTRATION-3: T1-T5 (see the outcome there).
  KNOWN and credited:
  - The bound nu(t) + nu(t + 1) >= floor(n/3), with its consecutive-triples proof and the per-colour window
    certificates (P(A), P(B) here, called nu_c there), is Session26 (session26_random_bridges_20261022/session26-blue/
    RESULT.md section 1, audited). It is stated there for the feedback vertex number tau >= nu, but the proof exhibits
    disjoint 3-cycles. Session26 shows sharpness only on the transitive line. THM-554 has the consecutive-triple
    criterion.
  - The tiling frame, sigma, tau and blue/black lines: THM-643, THM-646, THM-4729.
  - Even n, all blue lines have H = 1 (mod 4): THM-791.
  - The odd-order mod-4 balance and its mechanism (H = 1 + 2F mod 4 with F the symmetric cycles; the involution
    reversing internal pair arcs; the parity of the pair path count via Redei on a quotient tournament): THM-4922 (4),
    proof supplied by its audit (first_creators_genealogy_20261107/audit/AUDIT.md section 4). This session
    re-derived the decomposition and the quotient-tournament parity before finding that proof. They are not claimed.
  NEW: the tight-line analysis (1)-(3) (equality case, staircase constructions, block law, counts); the nu-pair census
  (4); the line-level mod-4 data and the degeneracy count (5)(c), (d).
  NOT CLAIMED: the word law at n = 3k, 3k + 1 (HYP-9540 (A)); which non-tight pairs occur (HYP-9540 (B)).
session: opus-2026-11-09-S87
source: 05-knowledge/results/ternary_metagraph_20261109/README.md
related:
  - THM-643
  - THM-646
  - THM-4729
  - THM-791
  - THM-4922
  - THM-554
  - HYP-9540
---

# THM-4951 - packing-number lines of the tiling metagraph

## Setting

**Tilings and lines.**
- The base path is 0 -> 1 -> ... -> n-1. The tiles are the pairs (i, j) with j >= i + 2, and a tiling t marks the
  reversed tiles (j -> i).
- The line through t is {t, t + 1}, where t + 1 reverses every tile. sigma(i, j) = (n-1-j, n-1-i).
- A line is blue if t = sigma(t), and black otherwise. This is THM-643's frame with the vertex order reversed.

**Packing number and the triple word.**
- nu(T) is the maximum number of vertex-disjoint cycles. Every cycle contains a 3-cycle on its vertex set, so this is
  also the maximum number of vertex-disjoint 3-cycles.
- k = floor(n/3).
- The triple word of t is x_a = [tile (a, a + 2) reversed], for a = 0..n-3. A = {a : x_a = 1} and B is its
  complement; the triple word of t + 1 is the complement of that of t.
- The window {a, a + 1, a + 2} is a 3-cycle of t iff a is in A, and of t + 1 iff a is in B (THM-554).

**Packing of positions and tight colourings.**
- P(X) is the largest subset of X whose elements differ pairwise by >= 3: the maximum number of disjoint windows of
  X. P([0, n-3]) = k.
- A colouring of [0, n-3] is tight if P(A) + P(B) = k.

**Tight lines.**
- Session26's argument gives nu(t) >= P(A), nu(t + 1) >= P(B) and P(A) + P(B) >= k. Hence nu(t) + nu(t + 1) >= k.
- A line is tight if equality holds.

**Staircases.**
- For an interval I of vertices, the staircase on I is the tiling that reverses exactly the tiles inside I.
- On I it induces the tournament with i -> i + 1 and j -> i for j >= i + 2.
- Its cycles are exactly the subintervals of length >= 3: from its least vertex a cycle can only step +1, and it
  closes at the first backward arc. So its packing number is floor(|I|/3).

## Statement

**(1) Tight lines are read off their words.**
- On a tight line, nu(t) = P(A) and nu(t + 1) = P(B), and the triple word is a tight colouring.
- The transitive tiling 0 is the only tiling with nu = 0.
- Nodes with nu <= d form an independent set, loops included, when n >= 6d + 3.
- For n >= 6 no blue line has both ends with nu <= 1.

**(2) Every split of the boundary occurs.**
- **(a) Prefix staircase.** For n >= 3 and 0 <= j <= k - 1, let t_j be the staircase on [0, 3j + 1]. Then
  nu(t_j) = j and nu(t_j + 1) = k - j.
  - The line is black for j >= 1.
  - Together with the transitive line (j = 0), every pair {i, k - i} occurs on a tight line.
  - tau(t_j) = j and tau(t_j + 1) = k - j as well, so Session26's tau bound is attained with every split.
- **(b) Double staircase.** For j >= 1, a = 3j + 2, c >= 0 and n = 2a + c, let t be the union of the staircases on
  [0, a - 1] and [a + c, n - 1]. It is blue, with nu(t) = 2j and nu(t + 1) = k - 2j. In particular, blue tight lines
  other than the transitive line exist for every n >= 10. At n = 9 the census gives 38.

**(3) The block law at n = 3k + 2.**
- **(a) Necessity.** The tight colourings of [0, 3k - 1] are exactly the block-monochromatic words on the aligned
  blocks {3b, 3b + 1, 3b + 2}. So the triple word of every tight line at n = 3k + 2 is a block word beta in
  {0, 1}^k, and nu(t) = |beta|_1.
- **(b) Sufficiency (multi-staircase).** For each block word beta, let t_beta be the union of the staircases on the
  windows [3s, 3s + 3r + 1] of the maximal runs [s, s + r - 1] of 1-blocks. Then:
  - the triple word of t_beta is beta;
  - nu(t_beta) = |beta|_1 and nu(t_beta + 1) = |beta|_0, so the line is tight;
  - F = {3b + 2 : beta_b = 0} is a feedback vertex set of t_beta + 1;
  - the line is blue iff beta is a palindrome.
- **(c) Counts.** The tight colourings up to complement number 2^(k-1) at n = 3k + 2 (by (a)). The DP of
  `staircase_check.py` gives (k + 1) 2^(k-2) at n = 3k + 1 and (k^2 + 9k - 2) 2^(k-4) at n = 3k, for every k <= 13
  (closed forms fitted, not proved).
- **(d) Blue tight lines.** sigma reverses the triple word, so the word of a blue line is a palindrome.
  - At n = 3k + 2 the blue tight lines are therefore exactly those of the palindromic block words, each realised by a
    blue multi-staircase.
  - A palindrome of even length has an even number of 1s. So the blue tight nu-pairs at n = 3k + 2 are exactly the
    {j, k - j} with j even or k odd.
  - Example: at n = 14 there is no blue tight line with pair (1, 3).

**(4) Census (FINITE-EXACT).** All lines by unordered nu-pair, black + blue (`nu_lines.c`; totals 2^(m-1) lines and
2^(e-1) blue lines checked):

| n | (0,k) | (1,1) | (1,2) | (1,3) | (2,2) | (2,3) | (3,3) |
|---|---|---|---|---|---|---|---|
| 4 | 0 + 1 | 2 + 1 | | | | | |
| 5 | 0 + 1 | 24 + 7 | | | | | |
| 6 | 0 + 1 | 56 + 0 | 202 + 8 | | 222 + 23 | | |
| 7 | 0 + 1 | 30 + 0 | 1220 + 45 | | 14878 + 210 | | |
| 8 | 0 + 1 | 12 + 0 | 5024 + 17 | | 1041492 + 2030 | | |
| 9 | 0 + 1 | | 3634 + 38 | 15294 + 113 | 732080 + 237 | 11480854 + 2893 | 121953098 + 29486 |

Tight lines (`tight_lines.c`; it reproduces the census rows for n <= 9):

| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| tight lines | 1 | 1 | 1 | 57 | 31 | 13 | 3673 | 1802 | 541 | 563,609 |
| blue | 1 | 1 | 1 | 1 | 1 | 1 | 39 | 10 | 19 | 297 |
| distinct triple words (x_0 = 1) | 1 | 1 | 1 | 5 | 3 | 2 | 17 | 8 | 4 | 50 |

- For k >= 2 and n <= 11, every tight line other than the transitive one has nu-pair (1, k - 1).
- At n = 12 (k = 4) the pairs are (1, 3): 96,254, all black, and (2, 2): 467,058 black + 296 blue. So there are no
  blue (1, 3) tight lines at n = 12.
- For n <= 12 the set of triple words of tight lines equals the set of tight colourings. The set of words of blue
  tight lines equals the set of palindromic tight colourings: 1, 1, 1, 1, 1, 1, 3, 2, 2, 3 words for n = 3..12.
- At n = 13 and 14 every tight colouring is realised (first-hit search; 20 and 8 colourings).

**(5) Blue lines by H mod 4.**
- **(a) Even n.** Every blue line is (1, 1). This is THM-791.
- **(b) Odd n.** On blue tilings, H = 1 and H = 3 (mod 4) are equinumerous, and the blue lines satisfy
  #(1,1) = #(3,3). This is THM-4922 (4)'s mechanism, restricted to tilings containing the base path. The involution
  only reverses the mirror tiles (i, n-1-i), which are never path arcs. The line statement follows by counting the
  tilings of each class.
- **(c) Splits (FINITE-EXACT).** (1,1) : (1,3) : (3,3) is
  - 2 : 4 : 2 at n = 5;
  - 64 : 128 : 64 at n = 7;
  - 8192 : 16384 : 8192 at n = 9;
  - 4,014,080 : 8,749,056 : 4,014,080 at n = 11, where the (1,3)-share is 0.52148, not 1/2. This is
    PREREGISTRATION-2, R3 REFUTED.
- **(d) What decides the split.** Write n = 2c + 1, x for the c mirror tiles and y for the other blue coordinates, and
  f = (H - 1)/2 mod 2 = A(y) + x.B(y) (THM-4922 (4), step 3). Let
  - D = {y : B(y) = B(y + 1)}, and
  - D_e = {y in D : A(y) + A(y + 1) + 1 = e}.
  Then
      #(1,3) = 2^(N + c - 2) + 2^(c - 2) (|D_1| - |D_0|),   N = number of 2-element sigma-orbits of tiles.
  - At n = 5, 7, 9, D is empty.
  - At n = 11, |D| = 94,208 = 23 * 2^12, |D_1| = 69,632 = 17 * 2^12 and |D_0| = 24,576 = 6 * 2^12.

## Proof

**(1).**
- **Tight lines.** On a tight line, nu(t) >= P(A), nu(t + 1) >= P(B) and P(A) + P(B) >= k, with the sum of the left
  sides equal to k. So all three are equalities.
- **nu = 0.** nu(t) = 0 forces t transitive. The only transitive tournament containing the path 0 -> ... -> n-1 is the
  one with every tile forward, so t = 0.
- **Independence.** nu <= d at both ends gives 2d >= k, so n <= 6d + 2.
- **Blue lines.** On a blue line, tile (0, 2) and its sigma-image (n-3, n-1) carry the same bit. The end where that
  bit is 1 contains the 3-cycles on {0, 1, 2} and {n-3, n-2, n-1}, which are disjoint for n >= 6, so that end has
  nu >= 2.

**(2a).** Let a = 3j + 2, L = [0, a-1], R = [a, n-1], b = n - a >= 1 (since j <= k - 1), and u = t_j + 1.

- **The t side.** t_j is the order join of the staircase on L and a transitive tournament on R. So
  nu(t_j) = floor(a/3) = j, and tau(t_j) = j by deleting every third vertex of L.
- **Structure of u.** On L, u is transitive (path arcs and forward tiles). On R, it is the staircase. Every cross tile
  (l, r) is reversed in u, so r -> l. Hence a -> l for every l <= a - 2, and the only arc from L to R is the path arc
  a-1 -> a.
- **Upper bound.** A cycle meeting L must leave L, so it contains a. At most one cycle of a packing meets L. The
  others lie in R - {a} = [a + 1, n - 1], where u is a staircase with floor((b - 1)/3) disjoint cycles. If no cycle
  meets L, the bound is floor(b/3) <= 1 + floor((b - 1)/3).
- **Lower bound.** The 3-cycle a-2 -> a-1 -> a -> a-2, together with floor((b - 1)/3) consecutive 3-cycles in
  [a + 1, n - 1].
- So nu(u) = 1 + floor((n - 3j - 3)/3) = k - j. The vertex a, plus every third vertex of [a + 1, n - 1], is a
  feedback vertex set of that size, so tau(u) = k - j.
- **Colour.** sigma maps the tiles inside [0, a - 1] to the tiles inside [n - a, n - 1]. For 1 <= j <= k - 1 these
  sets differ, so the line is black.

**(2b).** Let L1 = [0, a - 1], M = [a, a + c - 1] and L2 = [a + c, n - 1].

- **The t side.** t is sigma-fixed because sigma swaps the tiles inside L1 and L2. It is the order join
  L1 => M => L2 of two staircases and a transitive middle, so nu(t) = 2j.
- **Structure of u.** In u = t + 1, L1 and L2 are transitive and M is a staircase. All cross arcs point from later to
  earlier vertices except the path arcs.
  - The only arc out of L1 to later vertices is a-1 -> a, so every cycle meeting L1 contains a.
  - The only arc into L2 from earlier vertices is a+c-1 -> a+c, so every cycle meeting L2 contains a + c - 1.
- **Case c >= 2.**
  - Upper bound: nu(u) <= 2 + floor((c - 2)/3), since the other cycles lie in M - {a, a + c - 1}.
  - Lower bound: the 3-cycles {a-2, a-1, a} and {a+c-1, a+c, a+c+1}, plus consecutive 3-cycles of M - {a, a + c - 1}.
  - Then nu(u) = 2 + floor((c - 2)/3) = k - 2j, because n = 6j + 4 + c.
- **Case c in {0, 1}.** Every cycle meets both L1 and L2 or contains the single middle vertex, so nu(u) = 1 = k - 2j.
- **Blue tight lines for n >= 10.** j = 1, a = 5, c = n - 10.

**(3a).**
- **Block words are tight.** If the colouring is block-monochromatic with j blocks coloured 1, a 3-separated subset
  takes at most one position per block. Picking 3b from each 1-block gives P(A) = j, and similarly P(B) = k - j.
- **Tight words are block words.** Suppose block 0 = {0, 1, 2} contains positions a in A and b in B. Then
  P(A) >= 1 + P(A ∩ [5, L - 1]) and P(B) >= 1 + P(B ∩ [5, L - 1]), because a, b <= 2. Also
  P(A ∩ [5, L - 1]) + P(B ∩ [5, L - 1]) >= P([5, 3k - 1]) = k - 1. So P(A) + P(B) >= k + 1, and the colouring is not
  tight.
- **Induction.** If block 0 lies in A, then P(A) = 1 + P(A ∩ [3, L - 1]) and P(B) = P(B ∩ [3, L - 1]). So the
  colouring is tight iff its restriction to [3, 3k - 1] is tight, and induction on k finishes the proof.
- **Packing numbers.** nu(t) = P(A) = |beta|_1 by (1).

**(3b).**
- **The t side.** Distinct runs are separated by at least one 0-block, so the windows are disjoint (gap 3l - 2 >= 1
  vertices). A tile (a, a + 2) lies inside a window iff a lies in the run, so the triple word is beta. t_beta is the
  order join of staircases of sizes 3r + 2 and single vertices, so nu(t_beta) = sum r = |beta|_1.
- **The u side, lower bound.** u = t_beta + 1 contains the consecutive 3-cycles {3b, 3b + 1, 3b + 2} of the 0-blocks,
  so nu(u) >= |beta|_0.
- **The u side, upper bound.**
  - 3b + 2 lies in the window of a run [s, s + r - 1] iff s <= b <= s + r - 1. So F avoids every window.
  - In u, an arc from a smaller to a larger vertex is a path arc or a tile inside a window. So between the maximal
    path-intervals of V - F every arc points backwards, and every cycle of u - F lies in one interval.
  - Each such interval is exactly the window of the run of 1-blocks between two consecutive 0-blocks, or before the
    first or after the last one. It may also be a 2-vertex interval with no tile. In u that window is transitive.
  - Hence u - F is acyclic, and nu(u) <= tau(u) <= |F| = |beta|_0.
- **Colour.** sigma reverses the triple word, so t_beta is sigma-fixed iff beta is a palindrome.

`staircase_check.py` verifies (2a) for n = 3..18 and all j, (2b) for j = 1, 2 and c = 0..6, and (3b) for all 62 block
words at n = 5, 8, 11, 14, 17, including the feedback-set claim and the colour.

**(5d).**
- If y is not in D, then f(t) + f(t + 1) = A(y) + A(y + 1) + 1 + x.(B(y) + B(y + 1)) is balanced over the 2^c values
  of x.
- If y is in D, it is constant, equal to e(y).
- Each line is counted twice over the pairs (y, x). This gives the formula.
- The n = 11 values follow from the exact joint counts of `mirror_check.c` (mirror_check_n11.out) and the
  (P3)-failure count 94,208.

## Remarks

- **Packing number versus feedback number.** nu <= tau always. Both staircase families are tight for tau as well.
- **What the tight boundary looks like.** Every tight line found so far is close to a union of staircases. At
  n = 3k + 2 the boundary is completely described by the 2^(k-1) block word pairs. At other n the words have more
  freedom (one short block, or local interleavings such as 1010 at n = 6).
- **A natural lens, not a new invariant of the classes.** nu is an isomorphism and converse invariant, so the line
  counts between nu-classes descend to the merged metagraph, exactly as Session26's tau did.
