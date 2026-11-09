---
id: THM-4951
title: "Packing-number lines of the tiling metagraph. Let nu(T) be the maximum number of vertex-disjoint cycles of T and k = floor(n/3). On every line {t, t + 1}, nu(t) + nu(t + 1) >= k (Session26's consecutive-triples mechanism). (1) On a tight line (equality) both packing numbers are read off the triple word x_a = [tile (a, a + 2) reversed]: nu(t) = P(A) and nu(t + 1) = P(B), where P is the largest 3-separated subset. (2) The boundary is attained with every split: for each 0 <= j <= k there is a tight line with nu-pair {j, k - j} (prefix staircase), and for each j >= 1 with 6j + 4 <= n a BLUE tight line with nu-pair {2j, k - 2j} (double staircase). (3) Run-window theorem: if every run of a colouring w of the positions has length >= 2, the union of the staircases on the windows [p, q + 2] of its 1-runs [p, q] has triple word w and packing numbers P(A), P(B), and it is blue iff w is a palindrome. At n = 3k + 2 the tight colourings are the 2^k block words on aligned blocks of 3, and at n = 3k + 1 the block words with one block of 2; so at every n != 0 (mod 3) the triple words of tight lines are exactly the tight colourings, and those of blue tight lines exactly the palindromic ones (at n = 3k + 2 the blue pairs are {j, k - j} with j even or k odd). Counts up to complement: 2^(k-1) at n = 3k + 2 and (k+1) 2^(k-2) at n = 3k + 1 (proved), (k^2 + 9k - 2) 2^(k-4) at n = 3k (fitted). (4) Exact census of all lines by nu-pair and colour for n <= 9, and of tight lines for n <= 12. (5) Blue lines by H mod 4 at odd n: (1,1) = (3,3) always (from THM-4922 (4)); the split is 1:2:1 at n = 5, 7, 9 but not at n = 11; the excess is governed by D = {y : B(y) = B(y + 1)}."
status: >
  PROVED:
  - (1): the first, second and fourth bullets;
  - (2);
  - (3a)-(3f), including the counts at n = 3k + 2 and n = 3k + 1. The n = 3k data quoted inside (3f) are
    FINITE-EXACT, not proved;
  - (5)(b)'s line corollary and (5d).
  KNOWN:
  - (1), third bullet: Session26's corollary, transferred to nu.
  - (5)(a): THM-791 (a).
  - (5)(b), the balance on blue tilings: the THM-4922 (4) audit proof, restricted to tilings containing the base path.
  FINITE-EXACT:
  - the count (k^2 + 9k - 2) 2^(k-4) at n = 3k (DP for k <= 13; the audit's brute force for k <= 7). The closed
    form was fitted here; it is PROVED by THM-5021 (2) (S88);
  - (4) and (5)(c);
  - the n = 3k run-window failure set (an isolated 0 between two 1-runs of length >= 2), for n <= 18 (follow-up
    audit).
  SAMPLED: (5)(e), the (1,3)-shares at n = 13, 15, 17 (PREREGISTRATION-4).
  Pre-registered and tested:
  - PREREGISTRATION-2: R1 CONFIRMED (and KNOWN, see below); R2, R3 REFUTED at n = 11.
  - PREREGISTRATION-3: T1-T5 CONFIRMED.
  - PREREGISTRATION-4: U1 REFUTED; U2-U4 CONFIRMED.
  KNOWN and credited:
  - The bound nu(t) + nu(t + 1) >= floor(n/3), with its consecutive-triples proof, the per-colour window certificates
    (P(A), P(B) here, nu_c there) and the independence corollary, is Session26
    (05-knowledge/results/session26_random_bridges_20261022/session26-blue/RESULT.md section 1, audited;
    session26-tournament/WIRING_REVIEW.md). It is stated there for the feedback vertex number tau >= nu, but the proof
    exhibits disjoint 3-cycles. Session26 shows sharpness only on the transitive line. THM-554 has the
    consecutive-triple criterion.
  - The staircase lemma (the cycles of a staircase are its intervals of length >= 3). The full staircase is the
    tiling-frame tournament T*_n, the transitive node's line partner: THM-4731 (6), with the structure lemma "the
    out-neighbours of i are i + 1 and all j <= i - 2"; THM-337, the "base-path staircase tournament". The concurrent
    THM-4952 (C) records that its odd-cycle packings are the disjoint odd intervals.
  - The tiling frame, sigma, the complement involution and blue/black lines: THM-643, THM-646, THM-4729.
  - Even n, all blue lines have H = 1 (mod 4): THM-791 (a). THM-791's "odd-n contrast" has the first record of the
    odd-order blue balance: at n = 5 the 16 blue tilings split 8/8 by H mod 4.
  - The odd-order mod-4 balance and its mechanism (H = 1 + 2F mod 4 with F the symmetric cycles; the involution
    reversing internal pair arcs; the parity (P2) of the pair path count via Redei on a quotient tournament):
    THM-4922 (4), proof supplied by its audit (first_creators_genealogy_20261107/audit/AUDIT.md section 4, steps 3-5).
    This session re-derived the decomposition and the quotient-tournament parity before finding that proof. They are
    not claimed (MISTAKE-742).
  NEW: the tight-line analysis (1)-(3) (equality case, staircase and run-window constructions, block laws, counts);
  the nu-pair census (4); the line-level mod-4 data and the degeneracy count (5)(c)-(e).
  NOT CLAIMED here: the word law at n = 3k (HYP-9540 (A)), now PROVED by THM-5020 (S88, junction lemma and the two
  sides of the run-window tiling), and which non-tight pairs occur (HYP-9540 (B); see THM-5020 (6), HYP-9590).
  AUDIT: independent audit PASS WITH CORRECTIONS (05-knowledge/results/ternary_metagraph_20261109/audit/AUDIT.md).
  - Every number was reproduced with independent code, including all lines for n <= 9, all blue lines for n <= 12
    and the tight lines for n <= 12.
  - The corrections are applied: credits m6, m7, m11; proof gaps m8 (2b), m9 (3b); notation m10; wording m14.
  - The n = 3k + 1 count was upgraded to PROVED, on the audit's suggestion m12.
  - (3e) and (3f) (the run-window theorem) were added after the audit. The follow-up audit (audit/FOLLOWUP.md) is
    PASS WITH CORRECTIONS: (3e) and (3f) are correct. They were re-verified for all 150,048 colourings with runs
    >= 2 at n <= 26 and by samples to n = 32. Its F1 is fixed: the n = 3k failure description was wrong.
  - MISTAKE-742.
session: opus-2026-11-09-S87
source: 05-knowledge/results/ternary_metagraph_20261109/README.md
related:
  - THM-5020
  - THM-5021
  - THM-643
  - THM-646
  - THM-4729
  - THM-4731
  - THM-337
  - THM-791
  - THM-4922
  - THM-554
  - THM-4952
  - HYP-9540
  - 05-knowledge/results/session26_random_bridges_20261022/session26-blue/RESULT.md
---

# THM-4951 - packing-number lines of the tiling metagraph

## Setting

**Tilings and lines.**
- The base path is 0 -> 1 -> ... -> n-1. The tiles are the pairs (i, j) with j >= i + 2, and a tiling t marks the
  reversed tiles (j -> i). There are m = C(n-1, 2) tiles.
- The line through t is {t, t + 1}, where t + 1 reverses every tile. sigma(i, j) = (n-1-j, n-1-i).
- A line is blue if t = sigma(t), and black otherwise. This is THM-643's frame with the vertex order reversed.
- At odd n = 2c + 1 the sigma-fixed tiles are the c mirror tiles (i, n-1-i). The remaining tiles form N two-element
  sigma-orbits. e = N + c is the number of sigma-orbits, so there are 2^e blue tilings and 2^(e-1) blue lines. The
  same count holds at even n with the corresponding fixed tiles.

**Packing number and the triple word.**
- nu(T) is the maximum number of vertex-disjoint cycles. Every cycle contains a 3-cycle on its vertex set, so this is
  also the maximum number of vertex-disjoint 3-cycles.
- k = floor(n/3), and L = n - 2 is the number of positions.
- The triple word of t is x_a = [tile (a, a + 2) reversed], for a = 0..L-1. A = {a : x_a = 1} and B is its
  complement; the triple word of t + 1 is the complement of that of t.
- The window {a, a + 1, a + 2} is a 3-cycle of t iff a is in A, and of t + 1 iff a is in B (THM-554).

**Packing of positions and tight colourings.**
- P(X) is the largest subset of X whose elements differ pairwise by >= 3: the maximum number of disjoint windows of
  X. P([0, L-1]) = k.
- A colouring of [0, L-1] is tight if P(A) + P(B) = k.

**Tight lines.**
- Session26's argument gives nu(t) >= P(A), nu(t + 1) >= P(B) and P(A) + P(B) >= k. Hence nu(t) + nu(t + 1) >= k.
- A line is tight if equality holds.

**Staircases (KNOWN; THM-4731 (6), THM-337).**
- For an interval I of vertices, the staircase on I is the tiling that reverses exactly the tiles inside I.
- On I it induces the tournament with i -> i + 1 and j -> i for j >= i + 2.
- Its cycles are exactly the subintervals of length >= 3: from its least vertex a cycle can only step +1, and it
  closes at the first backward arc. So its packing number is floor(|I|/3).
- In u = t + 1 the staircase on I becomes transitive on I.

## Statement

**(1) Tight lines are read off their words.**
- On a tight line, nu(t) = P(A) and nu(t + 1) = P(B), and the triple word is a tight colouring.
- The transitive tiling 0 is the only tiling with nu = 0.
- Nodes with nu <= d form an independent set, loops included, when n >= 6d + 3. This is Session26's corollary, with
  the same proof.
- For n >= 6 no blue line has both ends with nu <= 1.

**(2) Every split of the boundary occurs.**
- **(a) Prefix staircase.** For n >= 3 and 0 <= j <= k - 1, let t_j be the staircase on [0, 3j + 1]. Then
  nu(t_j) = j and nu(t_j + 1) = k - j.
  - j = 0 gives the transitive line; for j >= 1 the line is black.
  - So every pair {i, k - i} occurs on a tight line.
  - tau(t_j) = j and tau(t_j + 1) = k - j as well (tau the feedback vertex number), so Session26's tau bound is
    attained with every split.
- **(b) Double staircase.** For j >= 1, a = 3j + 2, c >= 0 and n = 2a + c, let t be the union of the staircases on
  [0, a - 1] and [a + c, n - 1]. It is blue, with nu(t) = 2j and nu(t + 1) = k - 2j.
  - So some blue tight line other than the transitive line exists for every n >= 10 (j = 1), and the pair
    {2j, k - 2j} occurs in blue for every n >= 6j + 4.
  - At n = 9 the census gives 38.

**(3) Block laws and the run-window theorem.**
- **(a) n = 3k + 2, necessity.** The tight colourings of [0, 3k - 1] are exactly the block-monochromatic words on the
  aligned blocks {3b, 3b + 1, 3b + 2}. So the triple word of every tight line at n = 3k + 2 is a block word beta in
  {0, 1}^k, and nu(t) = |beta|_1.
- **(b) n = 3k + 2, sufficiency (multi-staircase).** For each block word beta, let t_beta be the union of the
  staircases on the windows [3s, 3s + 3r + 1] of the maximal runs [s, s + r - 1] of 1-blocks. Then:
  - the triple word of t_beta is beta;
  - nu(t_beta) = |beta|_1 and nu(t_beta + 1) = |beta|_0, so the line is tight;
  - F = {3b + 2 : beta_b = 0} is a feedback vertex set of t_beta + 1;
  - the line is blue iff beta is a palindrome.
- **(c) Counts of tight colourings, up to complement.**
  - n = 3k + 2: 2^(k-1), by (a).
  - n = 3k + 1: (k + 1) 2^(k-2), by (e).
  - n = 3k: (k^2 + 9k - 2) 2^(k-4). This was fitted here: it agrees with the DP of `staircase_check.py` for k <= 13
    and with the audit's brute force for k <= 7. It is PROVED by THM-5021 (2).
- **(d) Blue tight lines at n = 3k + 2.** sigma reverses the triple word, so the word of a blue line is a palindrome.
  - The blue tight lines at n = 3k + 2 are therefore exactly those of the palindromic block words, each realised by a
    blue multi-staircase.
  - A palindrome of even length has an even number of 1s. So the blue tight nu-pairs at n = 3k + 2 are exactly the
    {j, k - j} with j even or k odd.
  - Example: at n = 14 there is no blue tight line with pair (1, 3).
- **(e) The word law at n = 3k + 1.**
  - The tight colourings of [0, 3k - 2] are exactly the block words for the decompositions into consecutive blocks
    of length 3, except one block of length 2.
  - Every run of such a word has length >= 2. So by (f) each is the triple word of a tight line, the run-window
    tiling, which is blue iff the word is a palindrome.
  - Up to complement there are exactly (k + 1) 2^(k-2) of them.
- **(f) Run-window theorem.** Let w be a colouring of [0, L - 1] in which every maximal run of 0s or of 1s has length
  >= 2. Let t_w be the union of the staircases on the windows [p, q + 2] of the 1-runs [p, q]. Then:
  - the triple word of t_w is w;
  - nu(t_w) = P(A) and nu(t_w + 1) = P(B), where P(X) is the sum of ceil(l/3) over the runs of X (l = run length);
  - so the line is tight iff w is a tight colouring, and t_w is blue iff w is a palindrome.

  At n = 3k + 2, t_w is the multi-staircase of (b). With (a) and (e) this gives HYP-9540 (A) and (C) at every
  n != 0 (mod 3).
  - At n = 3k the tight colourings may have runs of length 1. The run-window tiling still realises 5 of 5 at n = 6,
    16 of 17 at n = 9 and 47 of 50 at n = 12.
  - FINITE-EXACT for every n = 3k <= 18 (the follow-up audit, audit/runwindow_n3k_runs.out): it fails exactly when
    an isolated 0 lies between two 1-runs that both have length >= 2. That gives 0, 1, 3, 8, 20 failures at
    n = 6, 9, 12, 15, 18, and every failing word contains 1110111.
  - An isolated 0 next to a 1-run of length 1 is harmless. Examples: 1010 at n = 6; 1010000, 1111010, 1110101 and
    1010111 at n = 9.
  - All of them are realised by some tight line for n <= 12 (by (4)).

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

The audit also censused all blue lines for n <= 12 (audit_lines.out, audit_blue12.out). At n = 12 their nu-pairs are:

| pair | (0,4) | (1,4) | (2,2) | (2,3) | (2,4) | (3,3) | (3,4) | (4,4) |
|---|---|---|---|---|---|---|---|---|
| blue lines | 1 | 67 | 296 | 5,513 | 129,241 | 22,586 | 2,313,819 | 534,399,389 |

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
- At n = 13 and 14 every tight colouring is realised (first-hit search; 20 and 8 colourings). (e) and (f) now prove
  this for n = 13 and 14.

**(5) Blue lines by H mod 4.**
- **(a) Even n (KNOWN, THM-791 (a)).** Every blue line is (1, 1).
- **(b) Odd n (KNOWN mechanism; the line corollary is PROVED here).** On blue tilings, H = 1 and H = 3 (mod 4) are
  equinumerous. This is THM-4922 (4)'s involution restricted to tilings containing the base path: it only reverses
  the mirror tiles (i, n-1-i), which are never path arcs. Counting tilings of each class then gives
  #(1,1) = #(3,3) on blue lines.
- **(c) Splits (FINITE-EXACT).** (1,1) : (1,3) : (3,3) is
  - 2 : 4 : 2 at n = 5;
  - 64 : 128 : 64 at n = 7;
  - 8192 : 16384 : 8192 at n = 9;
  - 4,014,080 : 8,749,056 : 4,014,080 at n = 11, where the (1,3)-share is 0.52148, not 1/2. This is
    PREREGISTRATION-2, R3 REFUTED.
- **(d) What decides the split.** Write n = 2c + 1, x for the c mirror tiles and y for the N other blue coordinates,
  and f = (H - 1)/2 mod 2 = A(y) + x.B(y) (THM-4922 (4), step 3). Let
  - D = {y : B(y) = B(y + 1)}, and
  - D_e = {y in D : A(y) + A(y + 1) + 1 = e}.
  Then
      #(1,3) = 2^(N + c - 2) + 2^(c - 2) (|D_1| - |D_0|).
  - At n = 5, 7, 9, D is empty.
  - At n = 11, |D| = 94,208 = 23 * 2^12, |D_1| = 69,632 = 17 * 2^12 and |D_0| = 24,576 = 6 * 2^12.
- **(e) Sampled orders (PREREGISTRATION-4).** With s_n the (1,3)-share minus 1/2:

  | n | 13 | 15 | 17 |
  |---|---|---|---|
  | samples | 24,000,000 | 12,000,000 | 12,000,000 |
  | s_n | -0.000906 (se 0.000017) | +0.001689 (se 0.000022) | -0.000111 (se 0.000014) |
  | |D| / 2^N | 0.0274 | 0.0238 | 0.0092 |

  - `mirror_sample.c` reproduces n <= 11 exactly.
  - The audit's independent sampler (odd-cycle parities only) gives s_13 = -0.00109 (se 0.00011).
  - So the (1,3) excess is not one-signed: the conjecture "it persists at every odd n >= 11" is refuted at n = 13. Its
    replacement is HYP-9540 (D).

## Proof

**(1).**
- **Tight lines.** On a tight line, nu(t) >= P(A), nu(t + 1) >= P(B) and P(A) + P(B) >= k, with the sum of the left
  sides equal to k. So all three are equalities.
- **nu = 0.** nu(t) = 0 forces t transitive. The only transitive tournament containing the path 0 -> ... -> n-1 is the
  one with every tile forward, so t = 0.
- **Independence (Session26).** nu <= d at both ends gives 2d >= k, so n <= 6d + 2.
- **Blue lines.** On a blue line, tile (0, 2) and its sigma-image (n-3, n-1) carry the same bit. The end where that
  bit is 1 contains the 3-cycles on {0, 1, 2} and {n-3, n-2, n-1}, which are disjoint for n >= 6, so that end has
  nu >= 2.

**(2a).** Let a = 3j + 2, X = [0, a-1], R = [a, n-1], b = n - a >= 1 (since j <= k - 1), and u = t_j + 1.

- **The t side.** t_j is the order join of the staircase on X and a transitive tournament on R. So
  nu(t_j) = floor(a/3) = j, and tau(t_j) = j by deleting every third vertex of X.
- **Structure of u.** On X, u is transitive (path arcs and forward tiles). On R, it is the staircase. Every cross tile
  (l, r) is reversed in u, so r -> l. Hence a -> l for every l <= a - 2, and the only arc from X to R is the path arc
  a-1 -> a.
- **Upper bound.** A cycle meeting X must leave X, so it contains a. At most one cycle of a packing meets X. The
  others lie in R - {a} = [a + 1, n - 1], where u is a staircase with floor((b - 1)/3) disjoint cycles. If no cycle
  meets X, the bound is floor(b/3) <= 1 + floor((b - 1)/3).
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
- **Case c >= 2, upper bound.** Take a packing.
  - If one of its cycles meets L1 and another meets L2, the remaining cycles lie in M - {a, a + c - 1}. That gives
    at most 2 + floor((c - 2)/3) cycles.
  - If no cycle meets L1 (resp. L2), the remaining cycles lie in M - {a + c - 1} (resp. M - {a}) or in M. That gives
    at most 1 + floor((c - 1)/3) or floor(c/3) cycles, both <= 2 + floor((c - 2)/3).
  - A single cycle meeting both L1 and L2 contains both a and a + c - 1, and is covered by the same count.
- **Case c >= 2, lower bound.** The 3-cycles {a-2, a-1, a} and {a+c-1, a+c, a+c+1}, plus consecutive 3-cycles of
  M - {a, a + c - 1}.
- So nu(u) = 2 + floor((c - 2)/3) = k - 2j, because n = 6j + 4 + c.
- **Case c in {0, 1}.**
  - For c = 0, every cycle meets L1, since L2 is transitive, and hence contains a.
  - For c = 1, every cycle meeting L1 or L2 contains the middle vertex a, and no cycle lies in M = {a}.
  - So all cycles share a, and nu(u) = 1 = k - 2j.
- **Blue tight lines for n >= 10.** j = 1, a = 5, c = n - 10.

**(3a).** Here L = 3k.
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
- **Colour.**
  - sigma reverses the triple word, so a sigma-fixed t_beta has a palindromic beta.
  - Conversely, v -> n - 1 - v maps the window [3s, 3s + 3r + 1] of the run [s, s + r - 1] to the window of the
    mirrored run [k - s - r, k - s - 1]. So for a palindromic beta the set of tiles inside windows is
    sigma-invariant.

**(3e).**
- **Block words are tight.** A block word with j 1-blocks has P(A) <= j and P(B) <= k - j, because a block has
  diameter <= 2. Since P(A) + P(B) >= k, the colouring is tight.
- **Tight colourings are block words.** Induction on k, for colourings of [0, 3k - 2]. For k = 1 the statement is that
  {0, 1} is monochromatic. For k >= 2:
  - If positions 0 and 1 differ, then P(A) + P(B) >= 2 + P([4, 3k - 2]) = k + 1.
  - If 0, 1, 2 lie in A, then P(A) = 1 + P(A ∩ [3, ...]) and P(B) = P(B ∩ [3, ...]). The rest is a tight colouring of
    length 3(k - 1) - 1, and the induction applies.
  - If 0, 1 lie in A and 2 lies in B, then P(A) = 1 + P(A ∩ [3, ...]) and P(B) = P(B ∩ [2, ...]). The colouring of
    [2, 3k - 2], of length 3(k - 1), must be tight. By (a), shifted by 2, it is block-monochromatic on aligned blocks.
    With the block [0, 1] this is the decomposition (2, 3, ..., 3).
  - The case where 0, 1 lie in B is symmetric.
- **Count.** A block word determines its run lengths: one run has length = 2 (mod 3) and the rest = 0 (mod 3). So it
  determines its block word beta in {0, 1}^k and the run containing the short block. Conversely, each beta and each
  of its runs gives a word. The total is sum over beta of runs(beta) = 2^(k-1)(k + 1), or (k + 1) 2^(k-2) up to
  complement. The audit gave the same count by the recursion W(k) = 2^(k-1) + 2 W(k-1).

**(3f).**
- **Runs and P.** Distinct runs of the same colour are separated by >= 2 positions. A greedy choice in each run from
  its left end gives P(X) >= sum ceil(l/3), and subadditivity gives <=.
- **The t side.** The 1-windows are pairwise disjoint (p' >= q + 3), so t_w is the order join of staircases of sizes
  l + 2 and single vertices. Hence nu(t_w) = sum floor((l + 2)/3) = P(A).
- **The u side, lower bound.** For x in B the tile (x, x + 2) lies in no window, so it is reversed in u = t_w + 1, and
  disjoint B-windows {x, x + 1, x + 2} give disjoint 3-cycles. So nu(u) >= P(B).
- **The u side, upper bound.** For each 0-run [p0, q0], take the points min(p0 + 3i + 2, q0 + 1) for i >= 0 with
  p0 + 3i <= q0. Call their union F.
  - There are ceil(l/3) points per 0-run, and every B-window {x, x + 1, x + 2} with x in [p0, q0] contains one.
  - Each point is either free (outside all 1-windows) or q0 + 1, the first vertex of the next 1-window. So no point is
    interior to a 1-window.
  - Let I be an interval of consecutive vertices avoiding F, with |I| >= 3. It contains no B-window, so its positions
    [min I, max I - 2] are consecutive 1-positions. They lie in one 1-run, so I lies inside that run's window, where u
    is transitive.
  - Each window meets only one interval, and consecutive intervals are separated by F. So every arc of u between
    different intervals points from the later one to the earlier one.
  - Hence u - F is acyclic, and nu(u) <= tau(u) <= |F| = P(B).
- **Colour.** As in (3b): sigma maps the window of a run [p, q] to the window of its mirror run [L-1-q, L-1-p], so t_w
  is sigma-fixed iff w is a palindrome.

**Checks.**
- `staircase_check.py` verifies (2a) for n = 3..18 and all j, (2b) for j = 1, 2 and c = 0..6, and (3b) for all 62
  block words at n = 5, 8, 11, 14, 17, including the feedback-set claim and the colour.
- `run_window_theorem_check.py` verifies (3f) for all 284 colourings with all runs >= 2 at n = 5..13: the word, both
  packing numbers, |F|, and acyclicity of u - F. It verifies (3e) by brute force for k = 1..6.
- The audit re-verified (2a), (2b), (3b) and (3d) independently, including n <= 26 for (2b) and k <= 6 for (3b).

**(5d).**
- f(t + 1) = A(y + 1) + (x + 1).B(y + 1) = A(y + 1) + x.B(y + 1) + 1, using sum_i B_i(y + 1) = 1. That identity is
  (P2), step 5 of the THM-4922 (4) audit proof (|O| is odd).
- So if y is not in D, then f(t) + f(t + 1) = A(y) + A(y + 1) + 1 + x.(B(y) + B(y + 1)) is balanced over the 2^c
  values of x.
- If y is in D, it is constant, equal to e(y).
- Each line is counted twice over the pairs (y, x). This gives the formula.
- The n = 11 values follow from the exact joint counts of `mirror_check.c` (mirror_check_n11.out) and the
  (P3)-failure count 94,208. The audit reproduced them by odd-cycle parities (blue_parity.out).

## Remarks

- **Packing number versus feedback number.** nu <= tau always. The staircase and run-window families are tight for tau
  as well.
- **What the tight boundary looks like.** Every tight line found so far is close to a union of staircases.
  - At n != 0 (mod 3) the boundary words are exactly the block words, and the run-window tilings realise each.
  - At n = 3k the words allow runs of length 1, such as 1010 at n = 6. The words where an isolated 0 lies between two
    1-runs of length >= 2 (containing 1110111) need other constructions; this is FINITE-EXACT for n <= 18.
- **A natural lens, not a new invariant of the classes.** nu is an isomorphism and converse invariant, so the line
  counts between nu-classes descend to the merged metagraph, exactly as Session26's tau did.
