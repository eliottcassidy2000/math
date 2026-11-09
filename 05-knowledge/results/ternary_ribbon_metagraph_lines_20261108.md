# Ternary Chan-Pak on ribbon words (the odd-cycle packing polynomial is finite-state; base 3 misses exactly 10 and 31 in span 2 and 3; 64 has no strong realiser) and the metagraph lines between component, value-colour, first-creator and self-converse classes

**Session** mac-mini-2026-11-08-ternary (worktree `math-wt-ternary-metagraph-20261108`, a sibling of the repository,
removed at close-out). **Directive** (owner, 2026-11-08): "start another session proving ternary Chan-Pak and other
conjectures as you propose. think about the metagraph and black and blue lines and conjecture about the number of blue
or black lines that can be between new classifications of the nodes you can describe". **Concurrent:**
opus-2026-11-09-S87 received the same directive and runs the base-3 order-10 strong census and lane reduction
(THM-4950/4951, HYP-9540, PREREGISTRATION.md Q1-Q3) and the blue-line `H mod 4` law at order 11 (PREREGISTRATION-2:
the 1 : 2 : 1 split holds at n = 5, 7, 9 only, R1 proved by Redei); nothing below duplicates either. Both readings of
"ternary Chan-Pak" in the repository are treated: (a) base 3, the values of `H_3(T) = I(Omega(T), 3)` (THM-4940,
HYP-9530 (A)); (b) the ternary return forest `h -> 3h + {0, 2, 4}` (triple_reentry_interfaces_20261103).
**Pre-registration:** `04-computation/experiments/ternary_ribbon_metagraph_20261108/PREREGISTRATION.md` (P1-P9, pushed
at 360dfbac48 before any census ran); outcome appended there. **Canon:** THM-4952 (ribbon packing transfer and base-3
census), THM-4953 (metagraph lines), HYP-9541 (conjectures). Chan-Pak and ternary Chan-Pak remain OPEN.

**Status in one line.** PROVED: the odd-cycle packing polynomial of a ribbon word is a finite-state series (automaton
with 4 / 20 / 112 states for span 1 / 2 / 3), its span-2 minimal dimension is 14 for every base except 0 and +-2
(13 at +-2: the Hamiltonian-path base is the unique positive rank drop), the tower satisfies the x-bonacci recurrence
f(n) = f(n-1) + f(n-2) + (x-1) f(n-3) (exactly 2^n growth in base 3), every line of the order-n metagraph joins a
strong closed end to an open end whose cut set is counted by an explicit Moebius formula (so the open end is
non-strong for a share 2^(4-n)(1 + O(2^-n)) of the lines and never has exactly n - 1 components), blue lines join
only eternally-blue-valued nodes, the transitive node's line goes to the tribonacci tower, and the fibres 3 and 5
carry exactly n - 2 and 5(n - 3) lines. PROVED by complete finite enumeration: span-2 ribbon words attain every value
1 (mod 3) up to 10^4 except 10 and 31 (fixed point of the pruned state sets, so for all lengths), strong span-2 words
miss 10, 31, 64, and 64 is the base-3 count of no strong tournament (Busch + orders <= 8). FINITE-EXACT: the line
census by node class at n = 5..8 and by value at n = 9. NEGATIVE: near-join digits do not realise the return forest
(47 of 97 attainable h <= 199 lack R = h + 1 at order <= 8). OPEN: everything in HYP-9541.

## 1. Inheritance

THM-4880 (ribbon transfer for H; 27 profile states, rank 14; strong connectivity regular), THM-4881 (reachable set;
congruence coset; exact tripling), THM-4940 (H_b = I(Omega, b): multiplicative, 1 mod b, monotone; the 7- and 21-holes
in every base; base 3 exact to 5398 missing 10 and 31), HYP-9530 (ternary Chan-Pak), THM-4865 (c) (the tribonacci
tower), THM-4744 (Busch's floor f(n), CITED), THM-4873 (C2 digit heredity; no ternary parent covers everything),
HYP-9453 (near-join h + 2R), THM-643/644/4729/4730/785/793 (the tiling metagraph: tiles, lines, blue = grid-symmetric,
closed/open, per-node parities, SC = A000570), THM-4922/4930/HYP-9520 (first creators; eternal colours of values),
01-canon/definitions.md (tiles, tilings, base path). OCF / Grinberg-Stanley: I(Omega(T), 2) = H(T) (CITED). The
frontier dynamic programme for cycle packings is KNOWN as a method (as in THM-4880 (A)'s audit correction 4); the
explicit tournament automata and everything built on them are new.

## 2. Theorem A (the packing polynomial of a ribbon word is finite-state; PROVED)

Letters, words and the window are as in THM-4880. A packing is a set of vertex-disjoint odd cycles; I(Omega(T), x) =
sum over packings of x^(number of cycles). Restricted to the processed prefix {1..q}, a packing is a set of completed
cycles and a set of vertex-disjoint directed path fragments of unfinished cycles. The head of a fragment will be
entered from a later vertex by a back arc, so it lies in the window; the tail will be left by a forward arc, which
exists unless the later vertex has a back arc to it, so a tail outside the window is free. State: the multiset of
fragments (head age in [0, w), tail age or FREE, parity of the length). The new vertex is skipped, starts a fragment,
is appended after a tail (forward arc), prepended before a head (back arc with span in the letter), or both (two
fragments merge, or one fragment closes into a cycle, allowed iff the length is odd, weight x). Heads reaching age w
kill the state; tails reaching age w become free; the final state is empty. Hence
`I(Omega(T_u), x) = <empty| M_{a_n}(x) ... M_{a_1}(x) |empty>` with matrices over Z[x] whose entries are linear in x.
Reachable states: 4, 20, 112 for w = 1, 2, 3 (ribbon_packing.py). Verified on 360 random words of orders <= 9 (120
per span) against an independent set-partition brute force (Hamiltonian cycle counts of odd subsets), and at x = 2
against THM-4880's path automaton on the same words: the two automata agree, which is OCF checked on both.

**Minimal dimension (PROVED, exact).** For span 2 the Z-module of functionals u -> I(Omega(T_u v), x) has rank 14
over Q(x): a nonsingular 14 x 14 minor of the functional matrix has determinant `-x^13 (x - 2)(x + 2)`, and the rank
is 14 at more integer specialisations than the degree bound of a 15 x 15 minor (see the audit correction below). The rank drops exactly at x = 0 (to 1: the constant series) and at x = +-2 (to 13; one linear relation among
the fourteen generic basis functionals). **Audit correction (MISTAKE-740).** The session's first certificate of the
generic rank used eighteen specialisations and the bound "a nonzero 15 x 15 minor has degree <= 15", which is the
degree of a single letter matrix; the functional rows are products of up to four letter matrices, with entries of
degree <= 2 (at most one odd cycle, needing three vertices, closes in four steps), so the bound is 30 at depth 4 and
45 at depth 5. The audit found rank 14 at 65 integer points for the depth-4 matrix and 61 for the depth-5 matrix,
which proves that all 15 x 15 minors vanish identically and that the depth-4 span is letter-invariant;
rankdrop.py now tests 36 points against the computed degree bound. So the Hamiltonian-path base x = 2 is the unique positive base at which the
span-2 packing series has lower dimension, and 13 is THM-4880 (G)'s irreducible dimension (the path series' rank 14
counted the empty word, where the path automaton gives 0 and the packing automaton 1). Integer representations
(packing_intrep.py, echelon or LLL bases, verified on 400 random words each): span 2 at x = 2: rank 13, entries <= 2;
x = 3: rank 14, entries <= 5; x = 4: 14, <= 6; x = 5: 14, <= 21; x = 1, 6, 7: rank 14. Span 3 at x = 3: rank 70,
entries <= 4 after LLL (pk3_span3_intrep.txt); at x = 2: rank 61, which is THM-4880's path rank 62 without the empty
word: for span 3 the Hamiltonian-path base drops the dimension by 9, not by 1.

**Letters in base 3** (characteristic polynomials on the 14 reduced coordinates): `--`: x^13 (x - 1); `1-`:
x^11 (x - 2)(x^2 + x + 1); `-2`: x^3 (x - 1)(x + 1)(x^9 - x^8 - 6x^6 - 11x^4 - 5x^3 - 14x^2 - 5x - 6); `12`:
x (x - 1)(x + 1)(x^11 - x^10 - x^9 - 4x^8 - 6x^7 - 14x^6 - 23x^5 - 28x^4 - 24x^3 - 17x^2 - 9x - 2). Growth of the
constant-letter runs: 2.000000, 2.446327, 2.718009 (base 2: 1.839287, 2.288291, 2.562972). No two-letter block has a
characteristic polynomial x^13 (x - c) with c != 1 in base 3, whereas in base 2 six ordered blocks have c = 3
(THM-4881 (B)): the exact tripling is a base-2 phenomenon (an ordinal sum with C3 multiplies H_3 by 4, and no short
block does that exactly).

## 3. Theorem B (the tower; PROVED)

For the one-letter word `1-^n` (the sole-exit tower of THM-4865) the odd cycles are exactly the odd intervals of length
>= 3 of the vertex order, traversed i -> j -> j-1 -> ... -> i+1 -> i, and a packing is a set of disjoint odd intervals;
classifying by the last vertex (uncovered, or in an interval of odd length 2j+1) gives f(n) = f(n-1) + x sum_{j>=1}
f(n-2j-1), hence `f(n) = f(n-1) + f(n-2) + (x - 1) f(n-3)`, f(0) = f(1) = f(2) = 1 (the x-bonacci recurrence; verified
symbolically against the automaton for n <= 13): 1, 1, 1, x+1, 2x+1, 4x+1, x^2+6x+1, 3x^2+9x+1, 8x^2+12x+1, ... At
x = 2 this is tribonacci (THM-4865 (c), KNOWN). At x = 3 it is `1, 1, 1, 4, 7, 13, 28, 55, 109, 220, 439, 877, 1756,
3511, 7021, 14044` = OEIS A122552 (a(n) = a(n-1) + a(n-2) + 2a(n-3); KNOWN sequence, the tournament reading is new),
with characteristic polynomial (x - 2)(x^2 + x + 1): `H_3(1-^n) = (3 * 2^n + c_n)/7`, c_n = 4, 1, -5 for n = 0, 1, 2
(mod 3), for every n >= 0. In base 3 the tower grows exactly like 2^n.

## 4. Theorem C (base-3 ribbon census; PROVED by finite fixed point) and Theorem D (64)

ribbon_levels_mod.c is ribbon_levels64.c of THM-4880 with values indexed by (v - 1)/M (M = 3), run on the base-3
integer representation and the 17-state (span 2) or 121-state (span 3) condensation automaton. Pipeline check at
M = 2 with the base-2 packing representation: the misses below 5,000 are 7, 21, 47, as in THM-4880.

**Span 2, H_3 <= 10^4.** Distinct reduced states per length (deduplicated on the reduced state and the condensation
state), values per length and cumulative misses among v = 1 (mod 3):

| length n | distinct states | values at length <= n | solid interval | cumulative misses (all) | cumulative misses (strong) |
|---|---|---|---|---|---|
| 5 | 68 | 7 | [13, 22] | 3,327 | 3,327 |
| 6 | 196 | 19 | [43, 61] | 3,315 | 3,315 |
| 7 | 575 | 46 | [34, 61] | 3,288 | 3,288 |
| 8 | 1,778 | 121 | [34, 175] | 3,213 | 3,214 |
| 9 | 5,573 | 306 | [34, 373] | 3,028 | 3,029 |
| 10 | 17,710 | 763 | [34, 1426] | 2,571 | 2,572 |
| 11 | 56,013 | 1,864 | [34, 3172] | 1,470 | 1,471 |
| 12 | 172,412 | 3,286 | [34, 7267] | 48 | 49 |
| 13 | 432,274 | 3,332 | [34, 10000] | 2 | 3 |
| 14 | 780,552 | 3,332 | [34, 10000] | 2 | 3 |
| 15 | 1,050,904 | 3,332 | [34, 10000] | 2 | 3 |
| 16 | 1,171,582 | 3,332 | [34, 10000] | 2 | 3 |
| 17 | 1,196,745 | 3,332 | [34, 10000] | 2 | 3 |
| 18 | 1,197,728 | 3,332 | [34, 10000] | 2 | 3 |
| 19 | 1,197,559 | 3,332 | [34, 10000] | 2 | 3 |
| 20 | 1,197,559 | 3,332 | [34, 10000] | 2 | 3 |
| 21-24 | 1,197,559 | 3,332 | [34, 10000] | 2 | 3 |

The pruned state sets are a fixed point from length 19 (1,197,559 states, unchanged through length 24); since the
level n + 1 is the image of level n under the letters pruned by H_3 <= 10^4 (sound by THM-4940 (1), monotonicity
under vertex deletion), the state set of every length is this fixed point, and the census is complete for all
lengths: **every v = 1 (mod 3) with v <= 10,000 other than 10 and 31 is H_3 of a span-2 ribbon word; every such v
other than 10, 31, 64 is H_3 of a strong span-2 ribbon word.** To 10^5 (pk3_s2_1e5_levels.txt): states 5,276,750 at length 15 and 12,872,609 at length 16, then the cap of
14,000,000 from length 17; every v = 1 (mod 3) with v <= 100,000 except 10 and 31 is attained by length 15 and every
such v except 10, 31, 64 by a strong word; no fixed point is claimed at this bound. The audit verified the fixed point
at 10^4 as set equality (S_19 = S_20 = S_21 = S_22; the program records cardinalities only) and noted that 10 and
31 are excluded for every tournament by THM-4940 (3)-(4) regardless.

**Span 3, H_3 <= 5,000** (pk3_s3_5e3_levels.txt; rank 70, the 121-state automaton): 1, 2, 8, 64, 384, 2,144,
11,668, 66,424, 395,672, 2,187,270 states at lengths 1..10, then the cap of 5,000,000; every v = 1 (mod 3) with
v <= 5,000 except 10 and 31 is attained by length 10, and every such v except 10, 31, 64 by a strong word.

**Theorem D (64).** 64 = 1 + 3 * 21 is not the base-3 packing count of any strong tournament. Proof: for a strong T,
H(T) = I(Omega, 2) <= I(Omega, 3) = H_3(T) (termwise), so H_3(T) = 64 forces H(T) <= 64 < 75 = f(9) (Busch, THM-4744),
hence |T| <= 8; the complete enumeration of all tournaments of orders <= 8 (spec3_small.py, 6,880 classes at order 8,
packing polynomials by the set-partition brute force) gives the strong base-3 values, and 64 is not among them (the
only values 1 (mod 3) below 1069 attained at orders <= 8 but not by a strong tournament are 1 and 64). 64 is attained
as the ordinal sum C3 => T with H_3(T) = 16 (H_3(C3) = 4). In base 2 the strong spectrum at orders <= 8 omits only 1
(finite data; THM-4880 (C) is the ribbon-family statement). Across bases (audit correction 3): at orders <= 8 the
attained values without a strong realiser are {1} in bases 2, 4, 5 and {1, 175}, {1, 120} in bases 6, 7 (audit
census of all 12,357 tournaments of orders <= 8; 175 and 120 exceed f(9) = 75 and stay undecided), and for b >= 8
every attained value below 64 is 1 + ab with a <= 7 in the intercyclic spectrum A (THM-4940 (2)), each such a being
the odd-cycle count of a strong intercyclic tournament of order <= 5; hence 64 is the smallest value other than 1, in
any base, that is attained but has no strong realiser. For comparison the same enumeration reproduces THM-4940's gap tables: base 2 misses 7, 21 (then only values
above 611, the order limit); base 3 misses 10, 31 (then above 1069); base 4 misses 13, 57, 109 (then above 1669);
base 5 misses 16, 51, 136, 371 (then above 2286).

## 5. Negative result (the return forest; FINITE-EXACT, P9 REFUTED)

The ternary return forest `h -> 3h + 2j` (j = 0, 1, 2; roots 1, 23, 25, 63, 65, 67 after deleting 7 and 21,
triple_reentry_interfaces_20261103) would prove Chan-Pak if every attained h had realisers "ready" for the digits 2
and 4. The near-join of (T, t) (two new vertices forming a 3-cycle with t, both beaten by T - t) has H = h + 2R_t(T)
(HYP-9453, k = 1) and is strong when T is, so digit 2j is near-join-ready at h iff some (T, t) with H(T) = h has
R_t = h + j. readiness.py computes every (h, R_t) over all tournaments of orders <= 8 and all vertices: of the 97 odd
h <= 199 attained at orders <= 8, only 50 admit R_t = h + 1 and only 55 admit R_t = h + 2; h = 3, 15, 27, 35, 41, 43,
57, 67, 83, 101, 105, 111, 121, 131, 139, 141, 143, 149, 151, 153, 157, 159, 165, 189, 195, 199 admit neither (complete
lists in readiness_n8.txt). Larger orders can add response values, but the small orders already show that near-joins
are not a digit system; this is the mechanism-level form of THM-4873 (iii) (no ternary parent covers everything) and
of THM-4880 (H) (no bounded-residual digit block in the ribbon family). Where readiness exists it always comes from a
strong realiser (every "1" in the table has a strong witness).

## 6. Theorem E (lines of the metagraph between node classes; PROVED) and the census (FINITE-EXACT)

Setting (01-canon/definitions.md; THM-643, THM-4729). Order n, base path n -> n-1 -> ... -> 1, m = C(n-1, 2) tiles
(a, b) with a >= b + 2, tiling t in {0,1}^m (bit 0: a -> b, bit 1: b -> a), tournament T_t, nodes = isomorphism classes
(merged with the converse class in the merged metagraph), lines {t, t^c} (all tiles flipped), blue iff t is
grid-symmetric (i -> n+1-i an anti-automorphism), closed iff tile (n, 1) carries 1 -> n. Every line has exactly one
closed end (THM-4729). New node classifications: (K) the number k of strong components; (V) the eternal colour of the
value H (THM-4930: eternally blue if some self-converse tournament has that count, eternally black otherwise; exact
below 1131); (F) first creator (THM-4922: a class of the least order realising its value) or not; (S) self-converse
or not (THM-643).

**(E1) Intervals and disjoint cut sets.** The strong components of a tournament with a Hamiltonian path P are
intervals of P (the condensation is a transitive order and P cannot return to an earlier component). Position i
(between path vertices i+1 and i) is a cut of T_t iff no tile crossing it (b <= i < a) is 1, and a cut of T_{t^c} iff
every tile crossing it is 1. Each position is crossed by a tile (n >= 3), so the cut sets of the two ends are
disjoint: `k(T_t) + k(T_{t^c}) <= n + 1`, with equality iff every position is a cut of one end, which forces all tiles
equal (tile (i+2, i) crosses both i and i+1): the transitive line {0^m, 1^m} only. Every tile crosses at least two
consecutive positions, so the non-cut positions of either end form runs of length >= 2: no end has exactly n - 1
components (the zero in every k-table).

**(E2) The cut-set formula.** Let the closed end be t (tile (n,1) = 1, so t has no cut: THM-4729's closed end is
strong). The open end has cut set containing S iff every tile crossing S is 1; the free tiles are those inside the
blocks of S, so the number of closed tilings whose open end has cut set containing S is `N_sup(S) = 2^(m-1)` for
S empty and `2^(e(S))`, e(S) = sum over blocks B of C(|B| - 1, 2), otherwise (the tile (n,1) crosses every position
and is already fixed). Moebius inversion on the Boolean lattice of cut sets gives the exact count
`N(S) = sum_{S' superset S} (-1)^|S' - S| N_sup(S')` of lines whose open end has cut set exactly S, and the number of
lines between the strong class and the class with k components is `L_k = sum_{|S| = k-1} N(S)`. For blue lines the
same inversion applies to `N_sup^blue(S) = 2^((m+f)/2 - 1)` (S empty) and `2^((e(Sbar) + f(Sbar))/2)` otherwise,
Sbar = S u sigma(S), sigma(i) = n - i, f(Sbar) = number of sigma-fixed tiles (a + b = n + 1) inside the blocks of
Sbar (a grid-symmetric tiling is free on orbit representatives of the free tiles). lines_formula.py evaluates both
and lines_census.c enumerates every line; the k-tables and the full cut-set tables agree exactly at n = 5, 6, 7, 8 and
9 (134,217,728 lines in eight parts, merge_parts.py). That the closed end is strong, that no line joins two
non-strong nodes and that the strong-strong count is 2^(m-1) - NS(n) are THM-4729 (3)-(4) (KNOWN; NS = 14, 121,
1,995, 64,648 at n = 5..8 reproduced); the new content is the full cut-set distribution and its blue refinement.

| n | lines | blue | open end strong (blue) | k = 2 | k = 3 | k = 4 | k = 5 | k = 6 | k = n |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 32 | 8 | 18 (6) | 10 (0) | 3 (1) | 0 | - | - | 1 (1) |
| 6 | 512 | 32 | 391 (27) | 101 (1) | 15 (3) | 4 (0) | 0 | - | 1 (1) |
| 7 | 16,384 | 256 | 14,389 (239) | 1,816 (0) | 153 (15) | 20 (0) | 5 (1) | 0 | 1 (1) |
| 8 | 1,048,576 | 2,048 | 983,928 (1,978) | 61,671 (5) | 2,739 (59) | 206 (2) | 25 (3) | 6 (0) | 1 (1) |
| 9 | 134,217,728 | 32,768 | 130,053,749 (32,249) | 4,067,314 (0) | 92,695 (501) | 3,672 (0) | 260 (16) | 30 (0) | 7 (1) at k = 7; 0 at k = 8; 1 (1) |

**(E3) Asymptotics.** The two largest terms of 2^(m-1) - N(empty) are the single cuts at positions 1 and n - 1
(the open end has a source or a sink), each 2^(C(n-2, 2)) = 2^(m-1) * 2^(3-n); the next, positions 2 and n - 2, are
2^(m-1) * 2^(6-2n). Hence the share of lines whose open end is not strong is `2^(4-n) (1 + O(2^-n))`: 0.5625, 0.7637,
0.8782, 0.9383, 0.9690 strong-strong shares at n = 5..9 against 1 - 2^(4-n) = 0.5, 0.75, 0.875, 0.9375, 0.96875. For blue lines a
single symmetric cut exists only at even n (position n/2), so the blue k = 2 count vanishes at odd n and is
2^(C(n/2 - 1, 2)) minus corrections at even n (5 at n = 8).

**(E4) Value fibres.** The transitive node (H = 1, one tiling) has one line; its partner 1^m is the tournament
i -> j iff j = i - 1 or j >= i + 2, i.e. the sole-exit tower of THM-4865, with H = 1, 3, 5, 9, 17, 31, 57, 105, 193,
355, 653, 1201 for n = 2..13 (tribonacci, A000213; KNOWN count, new identification). For h < f(n) every h-tiling is an
open end, so the fibre of h carries exactly sum_{C: H(C) = h, |C| = n} h / |Aut C| lines: the realisers of 3 are
TT_a => C3 => TT_b (|Aut| = 3, n - 2 classes), so fibre 3 carries exactly n - 2 lines for n >= 4; the realisers of 5
are TT_a => S4 => TT_b with S4 the strong 4-tournament (|Aut| = 1, n - 3 classes), so fibre 5 carries exactly
5(n - 3) lines for n >= 5 (census: 3, 4, 5, 6 and 10, 15, 20, 25 at n = 5..8). Fibre 9: 8, 37, 57, 78; fibre 15:
1, 34, 94, 152.

**(E5) Colours.** A blue tiling is grid-symmetric, so its tournament is self-converse and its value is eternally
blue; both ends of a blue line are blue tilings (THM-4729). Hence no blue line touches an eternally-black-valued node,
and every blue line joins two eternally-blue-valued self-converse nodes. From THM-643's per-node parities (blue odd and
black even at self-converse nodes, black odd at the others), the number of black lines between the self-converse class
and its complement has the parity of the number of non-self-converse nodes (12, 156, 5,044, 53,072 at n = 5..8 against
4, 44, 368, 6,704 nodes).

**Census by class (FINITE-EXACT; lines_nodes.py with pynauty canonical forms; THM-643's fibre sizes H/|Aut| and all
its parities re-verified at n = 5..8).**

| n | nodes (SC) | S: SC-SC blue / black, SC-NS, NS-NS | V: nodes blue / black; lines blue-blue (blue) / blue-blue / blue-black / black-black | F: first-creator nodes (values); FC-FC, FC-old, old-old (blue in brackets) |
|---|---|---|---|---|
| 5 | 12 (8) | 8 / 12, 12, 0 | 12 / 0; 8 / 24 / 0 / 0 | 6 (4); 18 (6), 14 (2), 0 |
| 6 | 56 (12) | 32 / 10, 156, 314 | 48 / 8; 32 / 284 / 176 / 20 | 33 (12); 371 (27), 131 (5), 10 (0) |
| 7 | 456 (88) | 256 / 876, 5,044, 10,208 | 370 / 86; 256 / 10,100 / 5,326 / 702 | 328 (58); 13,723 (215), 2,472 (34), 189 (7) |
| 8 | 6,880 (176) | 2,048 / 772, 53,072, 992,684 | 5,556 / 1,324; 2,048 / 702,560 / 308,080 / 35,888 | 5,457 (243); 932,047 (1,807), 104,099 (201), 12,430 (40) |

Lines between two eternally-black-valued nodes: 0, 20, 702, 35,888 (shares 0, 0.039, 0.043, 0.034); between two
eternally-blue-valued nodes: 1.000, 0.617, 0.632, 0.672. Same-value lines (both ends in one fibre): 7, 53, 470,
8,179 (shares 0.219, 0.104, 0.029, 0.0078); the closed end has the larger value in 24, 373, 11,838, 737,773 lines
(0.75, 0.73, 0.72, 0.70) and the smaller in 1, 86, 4,076, 302,624. At n = 8 the mean values are 417.8 (closed) and
359.6 (open), the correlation of the two end values is 0.741, and the smallest value on any line is the transitive 1
(with 57). At n = 9 (values only, merge_parts.py): fibres 3 and 5 carry 7 and 30 lines, same-value lines 267,735
(0.0020), the closed end is larger in 0.686 of the lines, the correlation is 0.800, the largest value is 3,357; of the
12,188,691 lines whose two values are both below 1131 (decided colours), 9,921,395 join two eternally blue values,
115,704 two eternally black and 2,151,592 one of each, while 91% of the lines carry a value >= 1131 whose eternal
colour THM-4930 does not yet decide. Multiplicities (iso-class nodes): the largest number of lines between two distinct nodes is 2, 6, 8, 12 at
n = 5..8, always at most min(H/|Aut|) of the two nodes; self-lines (both ends in one class) 4, 8, 44, 206 (blue 0, 2,
0, 4); the nodes whose distinct neighbours exhaust their fibre number 7, 17, 52, 213. Merged nodes (converse pairs identified; 10, 34, 272, 3,528 nodes at n = 5..8, of which 8, 12, 88, 176
self-converse; a pair node's fibre is 2H/|Aut|, verified): the largest number of lines between two distinct merged
nodes is 2, 14, 16, 26, attained at the value pairs (9, 9), (37, 43), (151, 153), (659, 653); the largest number of
blue lines between two distinct nodes is 1, 2, 2, 4 (at (657, 661) and (657, 629) for n = 8); self-lines (both ends
in one merged node, including the lines between a class and its converse) 4, 26, 114, 622 (blue 0, 2, 0, 4); the
largest neighbourhoods have 7, 20, 94, 461 distinct neighbours against fibres of 13, 62, 218, 1,022.

## 7. Conjectures (HYP-9541)

(1) The share of lines joining two eternally-black-valued nodes tends to 0 and the share joining two
eternally-blue-valued nodes tends to 1 (the value-colour form of HYP-9520). (2) The largest multiplicity between two
distinct merged nodes grows without bound but subexponentially, and is attained at pairs of large, nearly equal
values. (3) The closed end has the larger value in a share of lines bounded away from 1/2 and from 1 (0.70-0.75
observed), and the correlation of the two end values converges. (4) Base 3: the strong base-3 spectrum is the base-3
spectrum minus {1, 64}; the span-3 ribbon family attains every v = 1 (mod 3) except 10 and 31 (ternary ribbon
completeness, which implies HYP-9530 (A)). (5) For every span w the packing series' dimension at x = 2 equals the nonempty-word rank of THM-4880's
path series (13 for w = 2, 61 for w = 3) and lies strictly below the generic dimension (14, resp. >= 70); the drop at
the path base grows with w.

## 8. Audit

Independent adversarial audit (agent; `ternary_ribbon_metagraph_audit_20261108.md`): **PASS WITH CORRECTIONS**,
one major and eight minors, all applied on 2026-11-09 (MISTAKE-740). Major: the certificate of the generic rank in
THM-4952 (B) used the degree bound of a single letter matrix (15) for the minors of the functional matrix, whose
entries are products of letters (degree 2 at depth 4, 3 at depth 5), and the committed script tested fewer points
than the text claimed; the audit proved the generic rank with 65 and 61 specialisations and the script now tests
more points than the computed bound. Minors: the fixed point at 10^4 was inferred from equal cardinalities (the
audit verified set equality, and THM-4940 (3)-(4) exclude 10 and 31 anyway); the "first composite without a strong
realiser in any base" claim needed the audit's bases-6-7 census and the block-lemma argument for b >= 8; FC-FC at
n = 5 is 18, not 24; THM-4729 (3)-(4) are cited for the strong-strong count and the absence of non-strong-non-strong
lines; the hypothesis and the note now agree on "subexponential" multiplicity growth; the x = 2 identification is
spelled out; the base-2 strong spectrum at orders <= 8 is cited as finite data; the capped span-3 run is described as
capped. The audit reproduced with its own code the automaton (200 words, own cycle enumeration and own path counts),
the tower, the base-3 letter polynomials and blocks, the order-<= 8 census in bases 2-7, the near-join formula and the
readiness counts, every structural claim of THM-4953 and the Moebius formula at n = 5, 6, 7, and the node census at
n = 5, 6.

## 9. Scripts and receipts (`04-computation/experiments/ternary_ribbon_metagraph_20261108/`)

PREREGISTRATION.md; ribbon_packing.py (automaton, brute force, self-test); packing_intrep.py (pk{x}_span{w}_intrep.txt
for x = 2..7, span 2, and x = 3, span 3); rankdrop.py; ribbon_levels_mod.c (pk3_s2_1e4_levels.txt, pk3_s2_1e5_levels.txt,
pk3_s3_5e3_levels.txt, pk2_s2_5e3_levels.txt); spec3_small.py (spec3_orders_le_8.txt); readiness.py (readiness_n8.txt);
lines_census.c (lines_n{5..9}_p*.txt); lines_formula.py (formula_n{5..9}.txt); lines_nodes.py (lines_nodes_n{4..8}.json);
lines_analysis.py (lines_analysis.txt); pairs_analysis.py.
