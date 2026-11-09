---
id: THM-5020
title: "The word law for tight lines holds at every order. In the tiling line metagraph let w be a 0/1 word on the positions 0..n-3 with A = its 1-positions, B = its 0-positions, P = largest 3-separated subset, and R(w) the union of the staircases on the windows [p, q + 2] of its 1-runs [p, q]. (1) Junction lemma: if w has J junctions (an isolated letter between two runs of the other letter of length >= 2, i.e. a factor 11011 or 00100 centred on it) then P(A) + P(B) >= ceil((n - 2 + 2J)/3); so a tight colouring has at most one junction, and only at n = 0 (mod 3). (2) For every word, nu(R(w) + 1) = P(B). (3) If w avoids 11011, nu(R(w)) = P(A). (4) Hence at every n the triple words of tight lines (nu(t) + nu(t+1) = floor(n/3)) are exactly the tight colourings, realised by R(w) or by the line of R(complement w); the blue tight lines are exactly the palindromic ones. This proves HYP-9540 (A) and (C) at n = 3k, the open case. (5) Central-pair lemma: in a palindrome of even length, P(X) odd forces both central positions into X; so at even n no blue tight line has both packing numbers odd."
status: >
  PROVED: (1)-(5), for all n.
  FINITE-EXACT checks, independent of the proofs:
  - (1)-(3) on every word for n <= 15;
  - (4) on every tight colouring for n <= 24;
  - (5) on every palindrome of even length <= 24.
  The finite census of admissible words in (6) is also FINITE-EXACT.
  Pre-registered (PREREGISTRATION.md): Q1, Q2 REFUTED as registered. The run-window words realise every allowed pair
  of HYP-9540 (B) for 13 <= n <= 40 except black (k, k) at n = 3k, blue (k, k) at n = 0, 1 (mod 3), and blue
  both-odd non-tight pairs at even n. Q3 is not established.
  KNOWN and used:
  - THM-4951: the tiling frame, tight lines, (1) "tight lines are read off their words", the staircase lemma, and the
    run-window theorem for words whose runs all have length >= 2. This theorem extends that to every word without the
    factor 11011.
  - Session26: the line bound and its window certificates.
  - THM-554: consecutive triples.
  - THM-4731 (6), THM-337: the staircase T*_n.
  NEW: the junction lemma; the triangle classifications of R(w) and R(w) + 1 with the certificate argument (2) and
  the greedy feedback set (3); the word law and blue word law at n = 3k; the central-pair lemma.
  NOT CLAIMED: HYP-9540 (B) beyond the pairs realised in (6); HYP-9540 (D).
  AUDIT: PENDING.
session: opus-2026-11-09-S88
source: 05-knowledge/results/word_law_markov_20261109/README.md
related:
  - THM-4951
  - HYP-9540
  - THM-5021
  - THM-554
  - THM-4731
  - THM-5002
---

# THM-5020 - the word law for tight lines holds at every order

## Setting (as in THM-4951)

**Tilings, lines, packing number.**
- Path 0 -> 1 -> ... -> n-1. Tile (i, j), j >= i + 2, reversed iff j -> i. A line is {t, t + 1}, where t + 1 reverses
  every tile. It is blue iff t = sigma(t), with sigma(i, j) = (n-1-j, n-1-i).
- nu = maximum number of vertex-disjoint cycles = of vertex-disjoint 3-cycles. k = floor(n/3).

**Words and packings.**
- L = n - 2. A word w on the positions 0..L-1 has A = {a : w_a = 1} and B = {a : w_a = 0}. Its complement is w-bar.
- P(X) is the largest subset of X with pairwise differences >= 3. It is computed by the greedy algorithm.
- w is tight if P(A) + P(B) = k = ceil(L/3).
- The triple word of a tiling t is x_a = [tile (a, a + 2) reversed].
- On every line, nu(t) >= P(A), nu(t + 1) >= P(B), and P(A) + P(B) >= k (Session26).

**Runs, windows, junctions.**
- A run is a maximal interval of equal letters. The window of a 1-run [p, q] is the vertex interval [p, q + 2].
- R(w) is the tiling that reverses exactly the tiles lying inside some window. Windows of distinct 1-runs meet in at
  most one vertex: they share q + 2 exactly when a single 0 separates the runs. So every tile lies in at most one
  window.
- The triple word of R(w) is w: the tile (a, a + 2) lies in the window of [p, q] iff p <= a <= q.
- A junction is an isolated letter whose two neighbouring runs both have length >= 2. Equivalently it is the centre
  of a factor 11011 (a 0-junction) or 00100 (a 1-junction). J0, J1 and J = J0 + J1 count them.

## Statement

**(1) Junction lemma.** For every word, P(A) + P(B) >= ceil((L + 2J)/3).
- So a tight colouring has 2J <= 3k - L. This slack is 0, 1, 2 at n = 2, 1, 0 (mod 3).
- Hence J = 0 at n != 0 (mod 3), and J <= 1 at n = 0 (mod 3). No tight colouring has both a 0-junction and a
  1-junction.
- A palindromic tight colouring with a junction has it at the centre: sigma maps junctions to junctions, and there is
  only one.

**(2) Complement side.** For every word w, nu(R(w) + 1) = P(B).

**(3) Own side.** If w has no 0-junction (avoids 11011), then nu(R(w)) = P(A).

**(4) Word law and blue word law at every order.** For every n >= 3:
- a word is the triple word of a tight line iff it is a tight colouring;
- it is the triple word of a blue tight line iff it is a palindromic tight colouring.

A realisation:
- if w avoids 11011, the line of R(w);
- otherwise w avoids 00100 by (1), and the line of R(w-bar) has the endpoint R(w-bar) + 1 with triple word w.

The line is blue iff w is a palindrome. This proves HYP-9540 (A) and (C) at n = 3k, the case THM-4951 left open. At
n != 0 (mod 3) it reproves THM-4951 (3a)-(3f).

**(5) Central-pair lemma.** Let w be a palindrome of even length L and X = A or B. If P(X) is odd, then both central
positions L/2 - 1 and L/2 lie in X. Since they carry the same letter, P(A) and P(B) are not both odd. So:
- at even n, no blue tight line has both packing numbers odd (for example no blue (1, 3) tight line at n = 12, as
  THM-4951 (4) observed);
- no admissible palindromic word realises a both-odd pair at even n.

**(6) Lines with prescribed nu-pairs (consequence of (2), (3)).**
- Every word w without 11011 gives a line with nu-pair (P(A), P(B)). It is blue iff w is a palindrome. Words without
  00100 give such lines through their complements.
- An exact census of these admissible words for 13 <= n <= 40 (pairs_check.out) realises every pair allowed by
  HYP-9540 (B), except:
  - black (k, k) at n = 3k;
  - blue (k, k) at n = 0, 1 (mod 3);
  - blue non-tight pairs with i and j both odd at even n.
- The first of these is impossible for words: P(A) = k at n = 3k forces A to contain the unique maximum packing
  {0, 3, ..., 3k - 3}. The third is (5).
- Where complete censuses exist, such pairs are realised by other lines: black (3, 3) at n = 9 (THM-4951 (4)); blue
  (3, 3) and (4, 4) at n = 12 (S87 audit census of all blue lines). There is no complete black census at n >= 10.

## Proof

**(1).** Let z be a 0-junction, so w_(z-2) = w_(z-1) = w_(z+1) = w_(z+2) = 1, and let w' be w with position z
deleted.
- Take optimal packings S'_A, S'_B of w' and map them into w by p -> p (p < z), p -> p + 1 (p >= z). Distances only
  grow, so both images stay 3-separated.
- The B-positions of w within distance 2 of z are none, since z +- 1 and z +- 2 are 1-positions. So z can be added
  to the image of S'_B.
- Hence P(A) + P(B) >= P(A') + P(B') + 1, and symmetrically for 1-junctions.
- Deleting a junction letter merges its two neighbouring runs into a longer run of the other letter. Every other
  junction keeps two neighbouring runs of length >= 2, so it remains a junction.
- Deleting all J junction letters gives a word of length L - J, for which P + P >= ceil((L - J)/3). So
  P(A) + P(B) >= ceil((L - J)/3) + J = ceil((L + 2J)/3).

**Triangles of R(w) and of its complement.** Write t = R(w) and u = t + 1.
- Arcs of t: path arcs and tiles outside every window point forward; tiles inside a window point backward.
- Arcs of u: path arcs and tiles inside a window point forward; tiles outside every window point backward.
- A cycle contains a 3-cycle on its vertex set, so nu is the maximum number of disjoint 3-cycles.

**Triangles of t.** Let {a < b < c} be a 3-cycle of t.
- If it is oriented a -> b -> c -> a, then c -> a is a window tile, so b lies in the same window. Inside a window the
  only forward arcs are path arcs. So b = a + 1, c = a + 2 and a is in A: the window triangle {a, a + 1, a + 2}.
- If it is oriented a -> c -> b -> a, then (b, c) and (a, b) are tiles in windows W' and W''. If W' = W'', then
  (a, c) would lie in W'. So W' and W'' are different windows sharing b. Then b is the end of W'' = [p'', b] and the
  start of W' = [b, e']: a shared-vertex triangle {a, b, c} with a <= b - 2 and c >= b + 2.

**Triangles of u.** Let {a < b < c} be a 3-cycle of u.
- The orientation a -> c -> b -> a is impossible. a -> c forward with c >= a + 2 means (a, c) lies in a window, which
  then contains b, so b -> c.
- So the cycle is a -> b -> c -> a with (a, c) outside every window. Either b = a + 1 or (a, b) lies in a window;
  either c = b + 1 or (b, c) lies in a window. This gives four types:
  - (T1) {x, x + 1, x + 2} with x in B.
  - (T2) {p - 1, p, c}: p is the start of a 1-run, c lies in its window, c >= p + 2.
  - (T3) {a, q + 2, q + 3}: q is the end of a 1-run [p, q], a lies in [p, q], and q + 3 <= n - 1.
  - (T4) {a, b, c}: b is the vertex shared by two windows [p'', b] and [b, e'], with a < b < c.

**(2).**
- **Lower bound.** The T1 triangles at a maximum 3-separated subset of B are disjoint.
- **Certificates.** Give each triangle D of u a position x(D) in B:
  - T1: x;
  - T2: p - 1, a 0-position because p starts a 1-run;
  - T3: q + 1, a 0-position because q ends a 1-run (q + 1 <= L - 1);
  - T4: b - 1, the isolated 0 between the two 1-runs.

  In every case x(D) + 1 lies in D.
- **Claim.** If D and D' are disjoint, then |x(D) - x(D')| >= 3. Suppose x = x(D) < x' = x(D') <= x + 2.
  - If x' = x + 1: x + 1 is in B, so D is not T2 or T4, where x + 1 is the start of a 1-run. So D is T1 or T3, and
    both contain x + 2 = x' + 1, which lies in D'.
  - If x' = x + 2, and D is T1 or T3: D contains x + 1 and x + 2, and D' contains x + 3.
    - If D' is T1 or T2, it contains x + 2.
    - If D' is T3 or T4, its 1-run ends at x + 1 and starts there (x is in B). So the run is {x + 1}, and D' contains
      its member a' = x + 1, or a' in {x + 1, x + 2} for T4.
  - If x' = x + 2 and D is T2 or T4: the 1-run starting at x + 1 has length 1 (x + 2 is in B). So the window is
    [x + 1, x + 3], and D contains x + 3 or (for T4 with c = x + 2) x + 2.
    - In the first case D' contains x' + 1 = x + 3.
    - In the second, every type of D' contains x + 2 or the run member x + 1.

  So in every case D and D' meet.
- **Conclusion.** The certificates of disjoint triangles form a 3-separated subset of B. So nu(u) <= P(B).

**(3).**
- **Lower bound.** The window triangles at a maximum 3-separated subset of A.
- **Upper bound.** Let g_1 < g_2 < ... be the greedy packing of A, and F = {g_i + 2}.
  - Every a in A lies in some [g_i, g_i + 2], so g_i + 2 lies in the window triangle {a, a + 1, a + 2}.
  - Consider a shared-vertex triangle {a, b, c}, with the windows [p'', b] and [b, e']. Since w has no 0-junction, one
    of the two 1-runs at the isolated 0 (b - 1) has length 1.
    - If the run before is {b - 2}, then a = b - 2. The greedy element covering b - 2 is b - 4 or b - 2, since b - 3
      is in B. Its +2 is a or b.
    - If the run after is {b}, then c = b + 2. The greedy element covering b is b - 2 or b, since b - 1 is in B. Its
      +2 is b or c.
  - So every 3-cycle of t meets F. Hence t - F is acyclic, and nu(t) <= tau(t) <= |F| = P(A).

**(4).**
- If w avoids 11011, the line of R(w) has nu = (P(A), P(B)) by (2), (3), sum k, and triple word w.
- Otherwise J0 = 1 and J1 = 0 by (1), so w-bar avoids 11011. Then the line of R(w-bar) has packing numbers
  (P(B), P(A)), and its endpoint R(w-bar) + 1 has triple word w.
- **Colour.** Under v -> n - 1 - v the window of a run [p, q] goes to the window of the mirrored run
  [L - 1 - q, L - 1 - p]. So R(w) is sigma-fixed iff w is a palindrome; the same holds for w-bar.
- **Converse.** Every tight line's triple word is a tight colouring (THM-4951 (1)), and sigma reverses triple words.

**(5).** Let h = L/2, so X is mirror-symmetric about h - 1/2.
- A 3-separated subset of X has at most P(X ∩ [0, h - 1]) elements in each half. So
  P(X) <= 2 P(X ∩ [0, h - 1]).
- The elements <= h - 2 and their mirror images >= h + 1 are 3 apart. So P(X) >= 2 P(X ∩ [0, h - 2]).
- If P(X) is odd, these two bounds must differ. Then h - 1 lies in X, and by symmetry so does h.

**Checks.**
- `word_law_check.py` verifies (1)-(3) on every word with n <= 15, and (4) on every tight colouring with n <= 24
  (0 failures, word_law_check.out).
- `palindrome_parity_check.py` verifies (5) on every palindrome of even length <= 24.
- `pairs_check.py` performs the census for (6).

## Remarks

- **Why n = 3k was hard.** Only at n = 0 (mod 3) does a tight colouring have room (slack 2) for a junction. At a
  0-junction the two 1-windows share a vertex, and the shared-vertex triangles {a, b, c} with a <= b - 3 and
  c >= b + 3 escape the greedy feedback set; 1110111 at n = 9 is the smallest case. The cure is to put the
  staircases on the other colour's runs. The two kinds of junction cannot coexist, so one side is always safe.
- **Symmetry and asymmetry.** The construction commutes with sigma and with complementation. The one asymmetric
  choice, which colour carries the windows, is forced by the colour of the junction. A palindromic word can carry a
  junction only at its centre, so blue lines never need an asymmetric choice.
- **Packing versus covering.** In (2) a feedback set may not exist. For w = 010 the complement has three
  pairwise-intersecting 3-cycles with no common vertex (nu = 1 < tau = 2). That is why (2) is proved by
  certificates (a packing argument), not by a feedback set.
