# Ribbon tournaments: bounded-span back arcs make the Hamiltonian-path count a finite-state series, and the Chan-Pak spectrum becomes the value set of a 14-dimensional integer automaton

**Session** mac-mini-2026-11-04-ribbon (worktree `math-wt-recursive-20261104`). **Directive** (owner, 2026-11-04):
"spend another session looking for recursive classifications of pieces of tournaments structure that can compose and
grow with size in order to come to the proof of the full {7,21} spectrum".
**Status in one line.** KNOWN as methods, new as explicit tournament automata: the ribbon transfer calculus (Theorem A)
and the regularity of strong connectivity (Theorem B). PROVED: the closure corollary (C), the end-twin lemma (D), the absence of 47 from span 2 and its presence in
span 3 (E), the run recurrences (F), the irreducibility of the span-2 representation (G) and two no-go statements
(H; the unipotency half re-proved by the audit, MISTAKE-719). FINITE-EXACT: the strong span-2 values below 10^5 miss exactly 7, 21, 47; the solid prefix of the value
set grows geometrically (section 7). HYPOTHESIS: HYP-9473 (span-3 completeness implies Chan-Pak; the ribbon tiling
law). The Chan-Pak conjecture itself remains OPEN; nothing here is a proof of completeness (MISTAKE-559 applies).
Canon: THM-4880 (this note is its source), HYP-9473. Scripts: `04-computation/experiments/ribbon_spectrum_20261104/`.

## 1. Inheritance and the point of the construction

The lane's compositional tools all have *unbounded span*: the one-exit law and the sole-exit tower (THM-4862 (a),
THM-4865, tribonacci), the near-join `h + 2^k B` (HYP-9453), the one-arc join `H(A => B) = h_A h_B - e_a s_b + R_a R'_b`
(SEAM.md, golden_hamilton RECURRENCE (G1), THM-4875 (B)), the one-apex joins to 2^32 (THM-4875 (D)), Moon's interval
chains (THM-4744), the C2 grammar (klein). Their states are unbounded integer triples `(h, e, R)` (THM-4865 (f)) and the
missing lemma is always a density statement about those triples (the B-sets of HYP-9453 (2), the ready pairs of
HYP-9454). The present session asks instead for pieces whose *types* form a finite set: tournaments on an ordered
vertex set in which every arc points forward except back arcs of bounded span. Then (i) the Hamiltonian-path count is
computed by a finite profile automaton (Theorem A), (ii) strong connectivity is a regular property (Theorem B), and
(iii) the whole family's value set is the multiplicative closure of the values of a finite-alphabet N-rational series
(Corollary C). The Chan-Pak statement for this family becomes a statement about one explicit 14-dimensional integer
representation (span 2) or 62-dimensional one (span 3). Reduction: THM-4094 (completeness <=> strong realisation of the
lanes p != 7 and 7p), THM-1370 (7, 21 omitted), THM-4876 (no single-realiser induction), HYP-9029 (the tiling law for
all strong tournaments; the ribbon tiling law of section 7 is its finite-state analogue).

## 2. Definitions

A **ribbon word of span w** is a sequence of letters `a_1, ..., a_n`, each `a_q` a subset of `{1, ..., w}`; its
tournament `R(a)` has the vertices `1, ..., n`, the back arcs `q -> q - s` for `s in a_q` with `s < q`, and the forward
arc `p -> q` for every other pair `p < q`. The *window* at time `q` is `{q - w + 1, ..., q}`. Letters are written as
span masks: for span 2, `--` = {}, `1-` = {1}, `-2` = {2}, `12` = {1, 2}. `H` = number of directed Hamiltonian paths.
The **ribbon width** of a tournament is the least `w` for which some vertex ordering has all back arcs of span `<= w`
(width 0 = transitive). Every tournament has finite width, so the ribbon families exhaust all tournaments.

## 3. Theorem A (profile transfer; PROVED)

Fix `w`. Process the vertices in order. A Hamiltonian path `P` of `R(a)` restricted to `{1, ..., q}` is a set of
vertex-disjoint directed paths (fragments), the maximal runs of `P` inside `{1, ..., q}`. The head of a fragment that
is not the first vertex of `P` is entered from a later vertex, hence by a back arc of span `<= w`: the head lies in the
window. The tail of a fragment that is not the last vertex of `P` is left by a forward arc to a later vertex; a tail
outside the window has every arc to the future forward, so only its presence matters. Hence the **profile state**

    list of fragments, each (head, tail) with head = age in [0, w) or "unfed" (it left the window without a
    predecessor; at most one such fragment, which must contain the first vertex of P) and tail = age or "free",

determines exactly which extensions by the next vertex are possible: a new singleton fragment; appending after a tail
`t` (free, or at age `j` with span `j + 1` not in the letter); prepending before a head at age `j` (span `j + 1` in
the letter); or both, joining two distinct fragments. Ageing moves heads of age `w - 1` to "unfed" and tails to
"free"; two unfed fragments kill the state. Writing `N_s(a_1 ... a_q)` for the number of partial configurations of
state `s`, `N(a_1 ... a_{q+1}) = M_{a_{q+1}} N(a_1 ... a_q)` with the nonnegative integer matrix `M_a` of transition
multiplicities, and `H(R(a)) = sum over the one-fragment states of N_s(a_1 ... a_n)`. The reachable profile states
number 7, 27, 119 for `w = 1, 2, 3` (ribbon.py; verified against the subset DP on 332 random words of orders `<= 9`
for each `w`, 996 in all; the four-letter transfer also reproduces the exact value sets of alphabets.py and
probe47.py). The one-letter sub-alphabet `{1}` is the sole-exit tower of THM-4865 (every vertex's only earlier
out-neighbour is its predecessor), so the tribonacci law is the `w = 1` run of this calculus. **Attribution (audit
correction 4).** This is the standard frontier-based / bounded-bandwidth transfer-matrix dynamic programme for
Hamiltonian paths (Knuth, TAOCP 4A, section 7.1.4, simpath; frontier-based search on graphs of bounded pathwidth),
KNOWN as a method; what is new is the tournament specialisation (tails outside the window are free, at most one unfed
head, the explicit 7 / 27 / 119 states) and everything built on it below.

**Minimal dimension and the integer quotient (PROVED, exact arithmetic).** The Hankel rank of the series (pairing of
the reachable vectors with the co-reachable functionals, minrep.py) is 4 for `w = 1` and 14 for `w = 2`. The
Z-module generated by the functionals `u -> H(u v)` over all suffixes `v` is a sublattice of `Z^27` of rank 14,
generated already by the suffixes of length `<= 3` (exact Q-ranks 1, 5, 11, 14, 14 for suffix lengths `<= 0, ..., 4`;
intrep.py's progress labels were one less than the true suffix length and have been renamed to generator depth),
invariant under every letter;
in its echelon basis the state of a prefix is `y(u) in Z^14` with `y(empty) = e_1`, `H(u) = y_2(u)`, and the letters
act by the integer matrices `A_a` of intrep.py with all entries in `[-5, 5]` (`span2_intrep.txt`; verified on 500
random words of length `<= 11`). The same construction for span 3 gives rank 62; the echelon basis has entries of
size `10^54`, but an LLL-reduced basis (lll_span3.py, sympy's `DomainMatrix.lll`) has entries `<= 6` and transition
matrices with entries in `[-4, 4]` (`span3_intrep.txt`). A **semantic basis** of the span-2 observables: the fourteen
functionals `u -> H(u v)` for `v` = (empty), `--`, `1-`, `-2`, `12`, `--.-2`, `1-.-2`, `-2.1-`, `-2.-2`, `-2.12`,
`12.-2`, `-2.-2.1-`, `-2.-2.-2`, `12.-2.1-` are a Q-basis; so the whole future of a prefix is determined by its
counts under these fourteen short extensions (algebra.py).

**Irreducibility (Theorem G, PROVED, exact).** The states of the nonempty words span a 13-dimensional space (the
empty word adds the fourteenth: `y_1 = y_2` on every nonempty word). On it the four letters generate the full matrix
algebra `M_13(Q)` (dimension 169), and so do the three nonempty letters alone (exact_algebra.py, python-flint, integer
arithmetic). Consequently there is no common invariant subspace and no lower-dimensional quotient of the span-2
series: every "numeration system" for the family must use all thirteen coordinates. Characteristic polynomials of the
letters on the 13-dimensional space: `--`: `x^12 (x - 1)`; `1-`: `x^10 (x^3 - x^2 - x - 1)`; `-2`:
`x^5 (x^8 - 2x^7 + x^6 - 4x^5 + 3x^4 - 6x^3 + 2x^2 - 4x + 1)`; `12`:
`x (x^12 - 2x^11 - 2x^9 - 2x^8 - 4x^7 - 4x^6 - 4x^5 + x^4 + 2x^3 + 3x^2 + 2x + 1)`.

## 4. Theorem B (strong connectivity is regular; PROVED)

Let `C_1 < ... < C_m` be the strong components of `R(a_1 ... a_q)` in condensation order. Add `v = q + 1` with past
out-neighbourhood `A = {q + 1 - s : s in a_{q+1}}` (inside the window); every other earlier vertex beats `v`.
(i) If `A` is empty, `v` is a new sink component after `C_m`. (ii) Otherwise let `i_0` be least with `C_{i_0}` meeting
`A`. The vertex `v` reaches exactly `C_{i_0} u ... u C_m` (from a vertex of `C_{i_0}` every later component is reachable;
nothing before `i_0` is reachable from `A u {v}`, since arcs into `C_i`, `i < i_0`, come only from components `C_j`
with `j < i`, none of which `A` meets). A component `C_j`, `j >= i_0`, reaches `v` iff some vertex of
`C_j u ... u C_m` lies outside `A`; these sets decrease with `j`, so the components merging with `v` form an initial
stretch `C_{i_0}, ..., C_{j_1}` of `[i_0, m]`, and `C_{j_1 + 1}, ..., C_m` (all inside `A`; for `w = 2` singletons, since a
strong component is never a pair; for general `w` components contained in `A`) remain as components after the merged one; components before `i_0` are untouched. (iii) A component
containing no window vertex and lying before the first window-containing component can never be merged again (future
back arcs land in the window), so the final word is not strong. (Abstractly the regularity is KNOWN: strong connectivity is MSO-definable and the tournament is interpreted from the
word, so the language is regular by Buechi-Elgot-Trakhtenbrot; the explicit automaton is the content.) Hence the
automaton whose state is the list of
components from the first window-containing one onward (per component: the set of window ages it contains and whether
it contains an older vertex, needed only for the test "inside A"), with a dead state for (iii), decides strong
connectivity: 4, 17, 121 states for `w = 1, 2, 3` (condensation.py; agreement with brute-force reachability on 254,
21,844 and 64,680 words). Strong components of a ribbon word need **not** be intervals: the letters `--.1-.-2.1-` give
the arcs `2 -> 1, 3 -> 1, 4 -> 3` and the condensation `{2} < {1, 3, 4}`, although every cut between consecutive
vertices is crossed by a back arc. (The session's pre-registration file first asserted the cut criterion; corrected
within the hour, before any result depended on it; MISTAKE-718.)

## 5. Corollary C (closure; PROVED) and Lemma D (end twins; PROVED)

**C.** `H` is multiplicative over strong components (THM-1862/THM-4094). An induced subtournament of a ribbon word of
span `w` on any vertex subset is a ribbon word of span `<= w` (spans only shrink), so every strong component of a ribbon
word is a strong ribbon word of span `<= w`; and the ordinal sum `u => v` of two ribbon words is the ribbon word `u . v'`
where the `k`-th letter of `v'` is the `k`-th letter of `v` intersected with `{1, ..., k - 1}` for `k <= w` (those bits
pointed at nothing inside `v` and would point into `u`), the rest unchanged; for span 2: `u . -- . (second letter of v
intersected with {1}) . v_3 ...`. Hence the value set of the span-`w` family is
exactly the multiplicative closure of its strong values, and its primes are exactly its strong primes.

**D.** Let `X` be a tournament with a vertex `b` whose only out-neighbour is `a`; add `c` with `c -> a`, every vertex of
`X` other than `a, b` beating `c`, and either `b -> c` or `c -> b`. The two tournaments have the same `H`.
*Proof.* With `b -> c`, `c` has the sole exit `a` and `H_1 = h - e_a + R_a` (THM-4862 (a)), `R_a = sum_{S subset X - a}
h(X[S]) s_a(X - S)`. With `c -> b`: paths ending at `c` number `h - e_a - e_b`; paths `... c a ...` number
`sum_{S subset X - a} (h(X[S]) - e_b(X[S])) s_a(X - S)`; paths `... c b ...` number
`sum_{S subset X - b} (h(X[S]) - e_a(X[S])) s_b(X - S)`. Since `b`'s only exit is `a`: `e_b(X[S]) = h(X[S - b]) -
e_a(X[S - b])` for `b in S`, and `s_b(X - S) = s_a(X - S - b)` for `a notin S`, `s_b({b}) = 1`. The two sums share
`sum_{S' subset X - a - b} (h(X[S']) - e_a(X[S'])) s_a(X - b - S')`, and the difference is
`H_1 - H_2 = e_b(X) - (h(X - b) - e_a(X - b)) = 0`. Checked on 300 random tournaments of orders 4-9 (probe47.py).
Consequence: after the letter `1-`, the final letters `-2` and `12` give the same count (the profile states differ,
and later letters separate them; twins_state.py reports state-equal or H-equivalent pairs only, so its receipt for this
pair is the absence of a line).

## 6. Propositions E and F

**E (47; PROVED).** A strong tournament with `H = 47` has order `<= 8` (Busch's floor `f(9) = 75`, THM-4744). All
strong span-2 words of orders `<= 9` were enumerated (probe47.py, four letters): none has `H = 47`; `47` is prime, so
by C no span-2 ribbon word has `H = 47`. The two strong classes of order 7 with `H = 47` and the four classes of order 8
with `H = 47` (not strong: ordinal sums of an order-7 realiser with one vertex) have ribbon width exactly 3 (width.py,
bandwidth DP). Every tournament with `H = 47`, of any order, has width exactly 3: by C its strong component with
count 47 is one of the two order-7 classes, width is monotone under induced subtournaments (an ordering with spans
`<= w` restricts), and ordinal sums with transitive pieces have width `<= 3` by the closure of C. Explicit span-3 words
realise the order-7 classes
(witness47.py, re-verified by the transfer and by brute force):
`---.1--.---.-2-.--3.12-.-23` and `---.---.---.-23.1--.-2-.--3` (span-3 masks: `-`, `1`, `2`, `3` in the three
positions). Width distribution of all tournaments (FINITE-EXACT, width.py):

| order | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| width 0 | 1 | 1 | 1 | 1 | 1 | 1 |
| width 1 | 1 | 3 | 6 | 11 | 20 | 36 |
| width 2 | - | - | 5 | 38 | 194 | 778 |
| width 3 | - | - | - | 6 | 238 | 4953 |
| width 4 | - | - | - | - | 3 | 1112 |
| strong, width 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| strong, width 2 | 0 | 0 | 5 | 28 | 123 | 418 |

(strong classes by width: order 6: 1/28/6; order 7: 1/123/226/3; order 8: 1/418/4483/1106 for widths 1/2/3/4.)

**F (run recurrences; PROVED).** For a constant letter `a` the Krylov space of the initial state under `M_a` is
finite-dimensional and its minimal polynomial annihilates every coordinate, hence `H(a^n)` (recurrences.py, exact):
`1-`: the tribonacci law (THM-4865 (c)); `-2`: a recurrence of order 16; `12`: of order 21. The minimal recurrences
(minrec.py; Berlekamp-Massey on 60 terms, verified to n = 119, which with the proved recurrences of orders 16 and 21
forces them for all n) are: `-2`: for `n >= 9`, `H(n) = 2H(n-1) - H(n-2) + 4H(n-3) - 3H(n-4) + 6H(n-5) - 2H(n-6) +
4H(n-7) - H(n-8)`, terms `1, 1, 3, 5, 11, 29, 65, 149, 343, 777, 1781, 4081, 9329, 21361, 48881, ...`, growth
2.288291 (the largest root of `x^8 - 2x^7 + x^6 - 4x^5 + 3x^4 - 6x^3 + 2x^2 - 4x + 1`); `12`: an order-12 recurrence
with characteristic polynomial `x^12 - 2x^11 - 2x^9 - 2x^8 - 4x^7 - 4x^6 - 4x^5 + x^4 + 2x^3 + 3x^2 + 2x + 1`, terms
`1, 1, 1, 5, 15, 41, 105, 261, 659, 1689, 4347, 11165, 28623, 73329, 187875, ...`, growth 2.562972. Neither sequence
is in the OEIS (2026-11-04). Span-3 runs (recurrences.py): `3`: 1, 1, 1, 5, 13, 29, 77, 217, 617, 1601, 4263, ...
(growth 2.68); `123`: 1, 1, 1, 1, 9, 41, 175, 629, 2121, 6889, ... (growth 3.32).

## 7. The value sets (FINITE-EXACT)

**Strong span-2 words.** The DFS enumerator ribbon_enum.c (27 profile coordinates, the 17-state condensation automaton,
pruning by `H(prefix) <= H(word)` - the insertion lemma - and by dead condensation states) over all strong words with
`H <= 5000`: 19,764,437 prefixes, 6,133,659 strong words, every odd number `<= 5000` attained except `7, 21, 47`
(test_s2_misses.txt; the file's length histogram 1:1 3:1 4:1 5:3 6:5 7:12 8:25 9:45 10:90 11:171 12:329 13:625 14:787
15:370 16:32 records the lengths at which the depth-first search first encountered each value, not minimal lengths - the
minimal-length counts are the differences of the level table below, e.g. 109 - 42 = 67 new values at length 8). The
audit's own undeduplicated DFS (21 levels, 1.64 billion nodes, 643 s) to 50,000 gives the misses 7, 21, 47 as well.

**Deduplicated level enumeration.** ribbon_levels.c runs breadth-first over the reduced states `(y(u), condensation
state)`, merging prefixes with equal reduced state (they have the same `H` under every continuation), for `H <= N`.
With `N = 50,000` (lv2_5e4_levels.txt; the levels are nested because the empty letter keeps `H`):

| length n | distinct states | values at length <= n | solid prefix [49, d_n] | cumulative odd misses <= 50,000 |
|---|---|---|---|---|
| 8 | 1,716 | 109 | [49, 165] | 24,892 |
| 9 | 5,349 | 249 | [49, 257] | 24,752 |
| 10 | 16,914 | 598 | [49, 639] | 24,403 |
| 11 | 53,284 | 1,397 | [49, 1833] | 23,604 |
| 12 | 167,709 | 3,284 | [49, 3033] | 21,717 |
| 13 | 525,630 | 7,731 | [49, 8609] | 17,270 |
| 14 | 1,646,567 | 18,022 | [49, 19097] | 6,979 |
| 15 | 5,046,447 | 24,989 | [49, 45457] | 12 |
| 16 | 13,566,206 | 24,997 | [49, 49999] | 4 |
| 17 | 27,797,197 | 24,997 | [49, 49999] | 4 |

The strong values coincide with the values at every level from `n = 4` on (every attained value `<= 50,000` is
attained by a strong word of the same length bound; the strong bitmap is pre-seeded with 1, 3, 5, so the strong column
is only meaningful from `n = 4`; the raw file's header names eleven columns for ten numbers, the unprinted one being
the per-level strong count). With the length cap 17 the run finishes: the odd numbers `<= 50,000` not attained by words of length `<= 17` are exactly `7, 21, 47` (the fourth entry of the raw file, 50001, is the bound itself), and the same holds for strong words; since 7 and 21 are never attained (THM-1370) and 47 never in span 2 (Proposition E), the strong span-2 value set below 50,000 is exactly the odd numbers minus {7, 21, 47}, every value realised at length `<= 17`. An earlier run of the same program to length 18 (lv2_5e4_levels_partial18.txt) gives the identical rows through 17 before it was stopped for memory (27.8 million states at length 17). The solid prefix `d_n` grows by factors 1.6-2.9 per vertex: this is the
**ribbon tiling law** (HYP-9473 (3)). Per-length value sets of the three-letter alphabet alone (no empty letter; not
nested), from alphabets.py/probe47.py for `n <= 12`: solid intervals [49,165], [49,257], [129,639], [249,1833],
[397,3033] for `n = 8..12`, each overlapping the next (pre-registered P3 for `n = 13..18`). The deduplicated level run restricted to that alphabet (lv2s_1e5_levels.txt; the
letters at vertices 1-2 effectively truncated as in alphabets.py; `H <= 100,000`) reproduces these five intervals
exactly and continues with `[725, 8609]`, `[1281, 19097]`, `[2715, 45457]`, `[4681, >= 99999]`, `[8743, >= 99999]`,
`[14463, >= 99999]` for `n = 13..18` (the upper ends from `n = 16` on are the bound): every consecutive pair overlaps,
`c_{n+1} <= d_n`. The same run shows that every odd number `<= 100,000` other than `7, 21, 47` is attained by a
**strong** word of this alphabet of length `<= 16` (the accepting-state column coincides with the all-words column
at every level; distinct states 684, 1,934, 5,517, 15,741, 44,858, 127,722, 363,493, 1,031,691, 2,596,270,
4,458,863, 4,769,028 at `n = 8..18`). The strict variant in which every vertex `q >= 2` carries a real back arc
(lv2s_strict_1e5_levels.txt) has the intervals `[49,151], [129,257], [249,639], [397,1833], [725,2817],
[1281,7037], [2715,18575], [4681,42369], [8743,89309], [14463, >= 99999]` for `n = 8..17`, again overlapping.

**Span 3.** 47 is a span-3 count (E). The deduplicated level enumeration of span 3 (62 reduced coordinates, the 121-state automaton, `H <= 5,000`; lv3_5e3_levels_partial.txt) gives 66,350 / 394,690 / 2,417,506 distinct states at lengths 8, 9, 10 and, already at length 10, every odd number `<= 5,000` except 7 and 21 (by words of length `<= 10`; the raw file lists 7, 21 and the bound 5001); since 7 and 21 are never attained, the strong span-3 value set below 5,000 is complete, and the solid prefix is `[23, 4999]`. Level 11 exceeded the 8-million-state memory cap and was stopped. The 119-coordinate DFS to 20,000 (s3_2e4) was still running at the close-out commit.

**Pre-registration.** PREDICTIONS.md (committed 32d1356382 before any C run): P1 (strong span-2 words attain every odd
`<= 10^6` except 7, 21, 47) - confirmed to 100,000 (three-letter level run) and to 50,000 with all four letters, and to 5,000
by the independent 27-coordinate DFS; not to 10^6 (memory); P2 (span 3: every odd `<= 10^5` except 7, 21) - confirmed to 5,000 by the level method, not to 10^5 (memory); P3 - its overlap clause `c_{n+1} <= d_n` is CONFIRMED for `n = 13..18`; its growth clause
(`c_n < 2c_{n-1}` and `d_n > 2d_{n-1}`) is REFUTED: `c_10 = 129 > 2 c_9 = 98` and `d_9 = 257 < 2 d_8 = 330` already in the
data the prediction was written from (a careless clause; MISTAKE-718), and `c_15 = 2715 > 2 c_14 = 2562` out of sample. The claim in the same
file that every three-letter word is strong was false (MISTAKE-718) and was corrected before the enumerations ran.

## 8. Junctions, and two no-go statements (Theorem H, PROVED by exact finite computation)

For a block `J` the bilinear form `(u, v) -> H(u J v)` has a **coupling rank** `rank(Q M_J P)` (`P` a basis of the
26-dimensional reachable space, `Q` of the 14 co-reachable functionals; junctions.py). Among the 340 junctions of
length `<= 4` the ranks are 1 (96 junctions: ordinal sums, e.g. `--.--`, `--.1-`), 2 (27: the near-join `h + 2R` at
the last vertex followed by a sink, e.g. `--.-2.--`, `--.12.--`), 3 (153: one-arc joins, e.g. `--`, `--.-2`, `1-.--`,
`-2.1-`), 4 (34: e.g. `1-`, `1-.-2`), 8 (9), 9 (6), 10 (11) and 12 (4: `12`, `12.12`, ...). So the lane's one-arc
calculus (rank 3) and the near-join (rank 2) are sub-structures of the ribbon transfer, while a single `12` letter
couples almost everything (rank 12 of 14).

**No bounded-residual digit block (H1).** For the 583 blocks `B` of length `<= 5` over the three nonempty letters or of
length `<= 4` containing an empty letter, and the eight multipliers `alpha in {2, ..., 9}`, the residual
`H(uB) - alpha H(u)` over the 8,040 distinct prefix vectors of lengths 3-9 (7,873 distinct reduced states) is constant
for exactly 22 of the 4,664 (block, alpha) pairs - the ordinal-sum blocks with `H(B) = alpha`, residual 0: two of length
3 with `alpha = 3` (`--.--.-2`, `--.1-.1-`), and of length 4 eight with `alpha = 3` and twelve with `alpha = 5` (`--`
followed by a C3, resp. by the strong 4-tournament) - and takes at least 227 distinct values for every other pair
(digitscan.py, digitscan2.py; the audit's check_G_H.py for the full distribution). In other words no short block acts
as "multiply by alpha and add a digit from a finite set": the residual of every genuine extension is of the order of `H`
itself (for the letter `1-` it is `sigma - e`, THM-4862 (a)).

**No identity-plus-nilpotent block (H2; proof by the audit, MISTAKE-719).** For every block `B` of length `<= 4`,
`rank(M_B - I)` on the 26-dimensional reachable space is `>= 25` (recurrences.py) - a datum that does not exclude a
single nilpotent Jordan block of size 26, which is what the session had claimed it did. The proof: on the
14-dimensional quotient the eigenvalue 1 of the block matrix `A_B` has algebraic multiplicity 0 for 289 and 1 for 51 of
the 340 blocks, never 2 or more (exact characteristic polynomials, python-flint, audit); the reachable space maps onto
`Q^14`, so an identity-plus-nilpotent action on it would force `charpoly(A_B) = (x - 1)^14`; and a non-constant
arithmetic progression `H(u B^k) = H(u) + kc` would need a Jordan block of size at least 2 at the eigenvalue 1. Hence no
block inserts an arithmetic progression of values. (Footnote, audit remark 13: of the 96 rank-1 junctions, 94 leave a
cut of the vertex order uncrossed; `--.12.-2.--` and `--.12.-2.1-` are ordinal sums along a non-interval cut, the
MISTAKE-718 phenomenon again.)

Together with Theorem G these say: a proof of completeness for the ribbon family cannot go through a quotient
numeration, a bounded-carry digit system, or arithmetic progressions - it has to be a statement about the whole
13-dimensional set of reachable states (section 9).

## 9. What a proof would need (the content of HYP-9473)

(1) Chan-Pak follows from: *every odd number other than 7 and 21 is the count of a strong span-3 ribbon word* (by
Corollary C the strong words suffice; 47 needs span 3). The span-2 family is conjectured to miss exactly 7, 21, 47.
(2) The ribbon tiling law: the solid prefix `[49, d_n]` of the span-2 values of length `<= n` has `d_n -> infinity`
(equivalently every odd number `>= 49` is attained). (3) Because the representation is irreducible (G) and no digit
block exists (H), the law is a *collective* property of the reachable set `Y_n = {y(u) : |u| = n} subset Z^13`: the
image of `Y_n` under the 14 extension functionals must cover every odd number in a window. The natural inductive
object is therefore a two-dimensional density statement: for each `h` in the solid prefix, the set
`{H(u a) : H(u) = h, a a letter}` is sparse (extdensity.py: 3-40 values per `h` at length 10, with gaps), so the
cover comes from the union over `h`, exactly as for the full spectrum (THM-4876 (5)): the finite-state structure
does not make the induction local, it makes the *obstruction* exact. (4) Open: an invariant of `Y_n` (a lattice box in
some two coordinates, or a sumset structure) that propagates under the three letters; the two-sided pairing
`H(uv) = <chi(v), y(u)>` suggests building `Y_n` as a sumset of prefix and junction contributions.

## 10. Audit

Independent adversarial audit (agent; report `ribbon_transfer_audit_20261104.md`): **PASS WITH CORRECTIONS**, one
major and twelve minors, all applied on 2026-11-04. Major: the session's (H2) "no identity-plus-nilpotent block" was
justified by `rank(M_B - I) >= 25`, which is compatible with a single nilpotent Jordan block; the audit proved the
statement with exact characteristic polynomials (eigenvalue 1 at most simple for all 340 blocks) and the proof above is
the audit's (MISTAKE-719). Minors: the theorem's file name and the index lines called the 14-dimensional representation
irreducible (it is reducible; the 13-dimensional nonempty-word subrepresentation is irreducible; file renamed); (A) and
(B) are KNOWN as methods (frontier DP; MSO regularity) and are so typed; (H1) restated as the computed fact (22
constant pairs, `>= 227` residuals otherwise; 583 blocks, not "340/342"); the ordinal-sum rule of C for general `w`;
"singletons or a pair" in Theorem B; the width table's strong entries at orders 3-4 (0, not 1); intrep.py's progress
labels; the level receipts' header and seed artefacts; the DFS length histogram is a first-encounter record; the
all-order width argument for 47 and the non-strongness of the order-8 realisers; twins_state.py's receipt by absence.
The audit reproduced with its own code: both self-tests, the transfer on 200 random words and the reduced
representation on 100, the merge-lemma automaton (0 mismatches on all 87,380 span-2 words of length `<= 8`), the
counterexample's condensation, Lemma D (0 failures in 300, control 273/273 differ), 47's absence from strong span-2
words of orders `<= 8`, both span-3 witnesses, the width rows, the run sequences and growth constants (including
2.68076 and 3.32260 for span 3), the exact charpolys, the algebra dimension 169, the junction ranks
96/27/153/34/9/6/11/4, the 5,000 census, an independent undeduplicated DFS to 50,000 (misses exactly 7, 21, 47), the
level table rows `n = 8..17`, and the OEIS absence of both run sequences.

## 11. Scripts and receipts (`04-computation/experiments/ribbon_spectrum_20261104/`)

ribbon.py (engine; brute force; profile DP; transfer matrices; self-test), values.py (first BFS, superseded),
runs.py, minrep.py (+ minrep_w3.log), alphabets.py, probe47.py, extdensity.py, defects.py, twins_state.py,
digitscan.py, digitscan2.py, width.py, witness47.py, export_matrices.py (span2/3_matrices.txt),
condensation.py (strong-connectivity automaton; self-test), export_automaton.py (span2/3_condensation.txt),
recurrences.py, minrec.py, intrep.py (span2_intrep.txt; progress labels renamed to generator depth after the audit),
lll_span3.py (span3_intrep.txt), algebra.py,
exact_algebra.py, junctions.py, semantic_basis.py, ribbon_enum.c (test_s2_misses.txt; s2_1e5_*, s3_2e4_*),
ribbon_levels.c / ribbon_levels64.c (lv2_5e4_levels.txt, lv3_2e4_levels.txt, lv2s_3e5_levels.txt),
PREDICTIONS.md (pre-registration and its correction).
