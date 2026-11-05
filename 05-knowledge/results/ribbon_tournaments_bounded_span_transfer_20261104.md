# Ribbon tournaments: bounded-span back arcs make the Hamiltonian-path count a finite-state series, and the Chan-Pak spectrum becomes the value set of a 14-dimensional integer automaton

**Session** mac-mini-2026-11-04-ribbon (worktree `math-wt-recursive-20261104`). **Directive** (owner, 2026-11-04):
"spend another session looking for recursive classifications of pieces of tournaments structure that can compose and
grow with size in order to come to the proof of the full {7,21} spectrum".
**Status in one line.** PROVED: the ribbon transfer calculus (Theorem A), the regularity of strong connectivity
(Theorem B), the closure corollary (C), the end-twin lemma (D), the absence of 47 from span 2 and its presence in
span 3 (E), the run recurrences (F), the irreducibility of the span-2 representation (G) and two no-go statements
(H). FINITE-EXACT: the strong span-2 values below 5 x 10^4 miss exactly 7, 21, 47; the solid prefix of the value
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
out-neighbour is its predecessor), so the tribonacci law is the `w = 1` run of this calculus.

**Minimal dimension and the integer quotient (PROVED, exact arithmetic).** The Hankel rank of the series (pairing of
the reachable vectors with the co-reachable functionals, minrep.py) is 4 for `w = 1` and 14 for `w = 2`. The
Z-module generated by the functionals `u -> H(u v)` over all suffixes `v` is a sublattice of `Z^27` of rank 14,
generated already by the suffixes of length `<= 3` (length `<= 2` gives rank 11), invariant under every letter;
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
stretch `C_{i_0}, ..., C_{j_1}` of `[i_0, m]`, and `C_{j_1 + 1}, ..., C_m` (all inside `A`, hence singletons or a pair
of window vertices) remain as components after the merged one; components before `i_0` are untouched. (iii) A component
containing no window vertex and lying before the first window-containing component can never be merged again (future
back arcs land in the window), so the final word is not strong. Hence the automaton whose state is the list of
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
word is a strong ribbon word of span `<= w`; and the ordinal sum `u => v` of two ribbon words is the ribbon word
`u . -- . (second letter of v intersected with {1}) . v_3 ...` (the first letter of `v` becomes empty and its second
letter loses the span-2 bit, which pointed at nothing inside `v`). Hence the value set of the span-`w` family is
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
and later letters separate them: twins_state.py).

## 6. Propositions E and F

**E (47; PROVED).** A strong tournament with `H = 47` has order `<= 8` (Busch's floor `f(9) = 75`, THM-4744). All
strong span-2 words of orders `<= 9` were enumerated (probe47.py, four letters): none has `H = 47`; `47` is prime, so
by C no span-2 ribbon word has `H = 47`. The two strong classes of order 7 with `H = 47` and the four classes of order 8
with `H = 47` have ribbon width exactly 3 (width.py, bandwidth DP), and explicit span-3 words realise them
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
| strong, width 2 | 1 (C3 has width 1) | 1 | 5 | 28 | 123 | 418 |

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
(test_s2_misses.txt; first-attainment lengths 1:1 3:1 4:1 5:3 6:5 7:12 8:25 9:45 10:90 11:171 12:329 13:625 14:787
15:370 16:32).

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

The strong values coincide with the values at every level (every attained value `<= 50,000` is attained by a strong
word of the same length bound). With the length cap 17 the run finishes: the odd numbers `<= 50,000` not attained by words of length `<= 17` are exactly `7, 21, 47` (the fourth entry of the raw file, 50001, is the bound itself), and the same holds for strong words; since 7 and 21 are never attained (THM-1370) and 47 never in span 2 (Proposition E), the strong span-2 value set below 50,000 is exactly the odd numbers minus {7, 21, 47}, every value realised at length `<= 17`. An earlier run of the same program to length 18 (lv2_5e4_levels_partial18.txt) gives the identical rows through 17 before it was stopped for memory (27.8 million states at length 17). The solid prefix `d_n` grows by factors 1.6-2.9 per vertex: this is the
**ribbon tiling law** (HYP-9473 (3)). Per-length value sets of the three-letter alphabet alone (no empty letter; not
nested), from alphabets.py/probe47.py for `n <= 12`: solid intervals [49,165], [49,257], [129,639], [249,1833],
[397,3033] for `n = 8..12`, each overlapping the next (pre-registered P3 for `n = 13..18`: the three-letter level run lv2s_3e5 was still queued at this commit; its outcome is recorded in the close-out commit).

**Span 3.** 47 is a span-3 count (E). The deduplicated level enumeration of span 3 to 20,000 (62 reduced coordinates, 121-state automaton; lv3_2e4) and the 119-coordinate DFS to 20,000 (s3_2e4) were still running at this commit; their outcomes are recorded in the close-out commit.

**Pre-registration.** PREDICTIONS.md (committed 32d1356382 before any C run): P1 (strong span-2 words attain every odd
`<= 10^6` except 7, 21, 47) - confirmed to 50,000 by the level method (and to 5,000 by the independent 27-coordinate DFS; the DFS to 10^5 was still running at this commit),
not to 10^6; P2 (span 3: every odd `<= 10^5` except 7, 21) - not yet tested beyond 47 at this commit (runs queued, see above); P3 - queued at this commit. The claim in the same
file that every three-letter word is strong was false (MISTAKE-718) and was corrected before the enumerations ran.

## 8. Junctions, and two no-go statements (Theorem H, PROVED by exact finite computation)

For a block `J` the bilinear form `(u, v) -> H(u J v)` has a **coupling rank** `rank(Q M_J P)` (`P` a basis of the
26-dimensional reachable space, `Q` of the 14 co-reachable functionals; junctions.py). Among the 340 junctions of
length `<= 4` the ranks are 1 (96 junctions: ordinal sums, e.g. `--.--`, `--.1-`), 2 (27: the near-join `h + 2R` at
the last vertex followed by a sink, e.g. `--.-2.--`, `--.12.--`), 3 (153: one-arc joins, e.g. `--`, `--.-2`, `1-.--`,
`-2.1-`), 4 (34: e.g. `1-`, `1-.-2`), 8 (9), 9 (6), 10 (11) and 12 (4: `12`, `12.12`, ...). So the lane's one-arc
calculus (rank 3) and the near-join (rank 2) are sub-structures of the ribbon transfer, while a single `12` letter
couples almost everything (rank 12 of 14).

**No bounded-residual digit block (H1).** For every block `B` of length `<= 5` over the three nonempty letters and every
block of length `<= 4` containing an empty letter, and every multiplier `alpha in {2, ..., 9}`, the residual
`H(uB) - alpha H(u)` over the 8,040 distinct prefix states of lengths 3-9 takes a bounded set of values only when it is
identically 0, which happens exactly for the ordinal-sum blocks (`alpha = 3`: blocks `--` followed by a C3, e.g.
`--.--.-2`, `--.1-.1-`; `alpha = 5`: `--` followed by the strong 4-tournament) (digitscan.py, digitscan2.py). In other
words no short block acts as "multiply by alpha and add a digit from a finite set": the residual of every genuine
extension is of the order of `H` itself (for the letter `1-` it is `sigma - e`, THM-4862 (a)).

**No identity-plus-nilpotent block (H2).** For every block of length `<= 4`, `rank(M_g - I)` on the reachable space is
`>= 25` of 26 (recurrences.py); no block inserts an arithmetic progression of values.

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

An independent adversarial audit (agent, launched 2026-11-04 before the pre-audit commit) was asked to recompute the self-tests, the integer quotient, the condensation automaton, the end-twin lemma, Busch's floor, the recurrences, the irreducibility, the no-go scans, the 5,000 census and the junction ranks; its report is `ribbon_transfer_audit_20261104.md` and its outcome, with every correction applied, is recorded here in the close-out commit.

## 11. Scripts and receipts (`04-computation/experiments/ribbon_spectrum_20261104/`)

ribbon.py (engine; brute force; profile DP; transfer matrices; self-test), values.py (first BFS, superseded),
runs.py, minrep.py (+ minrep_w3.log), alphabets.py, probe47.py, extdensity.py, defects.py, twins_state.py,
digitscan.py, digitscan2.py, width.py, witness47.py, export_matrices.py (span2/3_matrices.txt),
condensation.py (strong-connectivity automaton; self-test), export_automaton.py (span2/3_condensation.txt),
recurrences.py, minrec.py, intrep.py (span2_intrep.txt), lll_span3.py (span3_intrep.txt), algebra.py,
exact_algebra.py, junctions.py, semantic_basis.py, ribbon_enum.c (test_s2_misses.txt; s2_1e5_*, s3_2e4_*),
ribbon_levels.c / ribbon_levels64.c (lv2_5e4_levels.txt, lv3_2e4_levels.txt, lv2s_3e5_levels.txt),
PREDICTIONS.md (pre-registration and its correction).
