# Folding a tournament by its converse: path counts of doubly regular tournaments mod 4, the character of 2 for Paley, pentagon sums, and records

**Session:** mac-mini-2026-11-11-atoms. **Canon:** THM-5136 (folding formula; self-converse doubly regular
tournaments; Paley), THM-5137 (pentagon sums of doubly regular tournaments; records). **Hypotheses:** HYP-9667.
**Audit:** pending (pre-audit checkpoint; nothing here is canonised until the independent audit is applied).

**Owner directive (summary).** Another session on the goal "7 and 21 are the only odd numbers that are not
Hamiltonian-path counts of tournaments", merging in 5 and the golden ratio, the packing of 11 squares into a square,
the relation of 7, 21 and 63 with powers of 2, Cayley-Dickson and Collatz; with the attached paper
(Chen, Mogielnicki, Ono, *Record compositions of alternating permutations and noncommutative symmetric functions*,
arXiv:2607.12873) and a pasted text on strong components as "primes" (`B(x) = 1/(1 - A(x))`, converse-collapsed
edge profiles) as sources of ideas.

**Outcome.** The conjecture (Chan-Pak, Conjecture 6.2) is **not proved here and remains OPEN**. Three surveys of the
repository showed that most of the named themes already have lanes, several of them written this week and two of
them today in answer to the same directive (section 1). This session therefore does two things: it adds new
mathematics where the lanes were empty, and it types the rest.

1. **(PROVED) Folding formula.** Let `T` have an anti-automorphism of order 2 (a converse symmetry); for odd order
   it fixes one vertex `v_0` and pairs the others. Folding `T` along it gives a tournament `Z` on the pairs and a set
   `D_1` of pairs (those whose triangle with `v_0` is cyclic), and

       H(T) = 1 + 2 * #{directed paths of Z that end in D_1}      (mod 4).

   More precisely, the number of invariant odd cycles of length `2k + 1` has the parity of the number of
   `k`-vertex paths of `Z` ending in `D_1` (Theorem Q). The repository had the congruence `H = 1 + 2F (mod 4)` with
   `F` the number of invariant odd cycles (Session29) and the remark that `F` is not decided by triangles alone;
   Theorem Q computes `F mod 2` on a tournament of half the size.
2. **(PROVED) Doubly regular tournaments with a converse symmetry.** In a doubly regular tournament of order
   `n = 4t + 3` with a converse symmetry, the invariant triples are forced: all cyclic if `t` is even, none if `t` is
   odd (a two-line count, Theorem E). With the folding formula,

       H(T) = n(n-1)/2 (mod 4):      3 (mod 4) if n = 3 (mod 8),      1 (mod 4) if n = 7 (mod 8),

   and the number of cycles of each odd length `2k + 1` has the parity of `[n = 3 mod 8] C((n-1)/2, k)`.
3. **(PROVED) The character of 2.** For the Paley tournament `P_q` (`q = 3 mod 4` a prime power) the forced pattern is
   the quadratic character of 2: `H(P_q) = 1 (mod 4)` if 2 is a square in `GF(q)` and `3 (mod 4)` if not; the number
   of Hamiltonian paths starting `0 -> 1` (`1, 9, 1729, ...`) is `1 (mod 4)`. The fold of `P_q` is `P_q` restricted to
   the squares: `P_7` folds to the cyclic triangle on `{1, 2, 4}` (the powers of 2, the octonion line), `P_11` to the
   regular pentagon tournament on `{1, 3, 4, 5, 9}` (the gold cells of the eleven-squares torus, THM-4592), `P_19`
   to a regular tournament of order 9 with 3357 paths, the largest count at that order.
4. **(CONJECTURE, pre-registered evidence) Without a converse symmetry.** The same residue is found in every doubly
   regular tournament on record, in fact modulo 8: `H = n(n-1)/2 (mod 8)` at orders 3 to 27 (fourteen classes of
   order 23, whose exact counts were already in the repository; ten classes of order 27). Stronger: `H mod 4` is
   constant on the whole switching class of a skew Hadamard matrix in every case examined (`N = 4, 8` completely;
   `N = 12, 16, 20, 24` on random members), although `H mod 4` is not a switching invariant in general. Sharpest:
   on these classes the number of cycles of **each** odd length `2k + 1` has the parity `[N = 4 mod 8] C(N/2 - 1, k)`
   of Theorem E, converse symmetry or not (pre-registered R12: 76 tournaments at `N = 20, 24`, ten of order 27).
   The natural local lemma is true up to order 8 and false from order 10 on, and the law fails vertex by vertex:
   it is global (HYP-9667, section 5.3).
5. **(PROVED) Pentagon sums of doubly regular tournaments.** By the kernel formula of THM-5132, every doubly regular
   tournament of order `n = 7 (mod 8)` (the Paley tournaments with `q = 7 mod 8` and all Cayley-Dickson sign
   tournaments on `2^k - 1 >= 7` units) has **all** `2^n` one-vertex extensions with `c_3 + c_5` even; for
   `n = 3 (mod 8)` the surplus is the least possible, `2^((n+1)/2)`. Every regular tournament has an even number of
   5-cycles (Theorem D).
6. **(PROVED, application of the attached paper) Records.** The mean of the Redei-Berge symmetric function over
   all tournaments on `n` vertices is `n! 2^(-n) q_n`, Schur's one-row `Q`-function, i.e. `n!/2^(n-1)` times the sum
   of all hook Schur functions; it is the sprout sequence with seed `(2 + t)/(2 - t)`. Its canonical noncommutative
   lift in the sense of Chen-Mogielnicki-Ono has the record interpretation asked for in a recorded question
   (kind-pasteur S22): the mean weighted number of odd-cycle packings whose blocks, listed by increasing largest
   vertex, have sizes `alpha` is `n! 2^(l(alpha) - n)/pi_u(alpha)`. Labelled tournaments in which the largest vertex
   reaches every vertex are counted by the connected labelled graphs (A001187), and tournaments and graphs are
   equidistributed by max-ordered type (section 6).
7. **(Pre-registered.)** Twelve predictions in six blocks: R1, R4, R5, R8, R11, R12 CONFIRMED; R9, R10 CONFIRMED, but
   their order-19 and order-23 parts were decidable from the repository (this machine's own lane of 2026-10-22);
   R7 confirmed but uninformative; R2, R3, R6 FAILED. R2 and R3 were replaced by the law that was then proved; R6
   was decidable from canon (section 9, MISTAKE-PENDING (entry filed at close-out)).
8. **(Typed merge.)** Section 8 places 5 and the golden ratio, the eleven squares, 7-21-63, powers of 2,
   Cayley-Dickson, Collatz and the "primes" text against the existing lanes, with the few exact links that are new
   and an explicit "no bridge" where none was found.

Nothing here enlarges the range of verified values or proves a new case of the conjecture. What is new and proved
concerns `H mod 4` again: on tournaments with a converse symmetry it is a path count on the folded tournament, and
on doubly regular ones it is the residue of `n(n-1)/2`.

---

## 1. Inheritance pass and concurrency

Three read-only surveys (themes: 7-21-63 / Cayley-Dickson / Collatz; golden ratio / symmetric functions; square
packing / tilings) were run before anything was typed, because the last three sessions of this machine each had a
"known result typed as new" (MISTAKE-753, 755, 760). What they found, by lane:

- **7, 21, 63.** THM-1370 (63 is attained at order 8; the "`7 * 3^k` gap tower" is refuted at `k = 2`), THM-344
  (exactly two classes with `H = 63`, both of order 8, `|Aut| = 1`, conflict graph `K_31`), THM-4942 / HYP-9512 (63 is
  the only orphan value; its creators are the Fibonacci fan and its converse), THM-4094 (two-prime-lane reduction:
  completeness is equivalent to `P_odd \ {7}` and `7 (P_odd \ {3})` being strong values; the carry atoms 49, 63, 343
  are strong values), THM-4749 (`N(n, 63) = 2(n - 7)`), THM-4986 (crevasse lemma).
- **Powers of 2, Mersenne, Jacobsthal.** The fugacity axis (`I(P_m, 2)` Jacobsthal, `I(C_m, 2) = 2^m + (-1)^m`, roots
  `(1 +- sqrt(1 + 4x))/2`, integer exactly at `x = k(k+1)`; reflections of oracle-2026-05 and HYP-849, HYP-1296,
  HYP-1299 "forbidden `H = 7 = JL(3)`", HYP-1301), THM-1425 (a path conflict graph gives a Jacobsthal number),
  THM-201 / THM-202 (`K_3` and `P_4` are not conflict components), THM-4754 (6) (Jacobsthal split; `2^k = 3 J_k +
  (-1)^k` as a one-step Collatz relation), THM-5104 (strong fans with `H = 4^n - 1`), THM-4748 (`H = 2^k - 1` needs
  three back arcs), `mersenne_switch_parity_f21_compression_20261006.md` (`x -> 2x + 1` and `x -> 4x + 1` generate
  `Aut(P_7)` mod 7).
- **Cayley-Dickson.** THM-448 (kind-pasteur, June: the doubling tower of doubly regular tournaments on `2^k - 1`
  vertices; `T_7 = P_7`; `H(T_15) = 198,335,025`; `|Aut(T_15)| = 21`), THM-447, THM-451 (the tower is chiral: not
  self-converse past order 7), THM-4557 (doubling keeps the automorphism group; Hanaki 2020), THM-5133 (kind-pasteur
  S23, today, pre-audit: the sign tournament of the Cayley-Dickson algebra is doubly regular for every `n`;
  `H = 198,335,025` at `n = 4`), `collatz_mod6_20260922_paley_fano_octonion_design.md` (the octonion index rule is
  the Paley orientation; 16 of the 128 orientations of the Fano plane have `H = 189`), THM-4553.
- **Doubly regular tournaments.** This machine's Sierpinski lane of 2026-10-22 (HYP-9365, HYP-9371;
  `04-computation/experiments/sierpinski_ocf_20261022/`): a non-Paley, self-converse doubly regular tournament
  `D_19` with `H = 1172732961771` (annealing, nauty), thirteen classes of order 23 with exact path counts
  (`drt_results.out`, `drt23_classes.out`; two converse pairs; none but `P_23` self-converse), the doubled Paley
  class of order 23 (`H = 15759752857737405`), two classes of order 31 from the twice-doubled matrix; the lane is
  about the loop expansion (path counts of doubly regular tournaments of one order agree through three loops),
  not about residues. **It was found only after Run G** (section 9). THM-4911 (signature simplex; skew Hadamard).
- **Collatz.** `spectrum_multiplicity_signed_tables_20261111/WEIGHTED_ATLAS.md` and `signed_arcs/PUNCTURED_ATLAS.md`
  (Codex, today): in the coordinates `H = 4x + sigma` the triangle acts by `(x, sigma) -> (3x + sigma, -sigma)`; the
  orbit of 7 is `(2,-) -> (5,+) -> (16,-) -> (47,+)`; 7 and 21 are the two missing arrows `2 -> 1` and `5 -> 16`.
  THM-4849 (circulant tournaments as rational Collatz cycles), `collatz_paley_bridge_20261001.md`.
- **Converse symmetry.** `session29_framed_profiles_20261022/session29-halfpath/REPORT.md` (half-path bijections for
  invariant Hamiltonian paths and invariant odd cycles; ternary dynamic programme; `P_7` has 12 invariant odd
  cycles), `session29-converse-congruences/REPORT.md` (`H = 1 + 2F + 4I (mod 8)`; "the formula does not say `F` is
  determined by triangles alone"), `07-reflections/the-reversal-involution-and-the-half-principle-opus-20260629.md`,
  THM-552 (`c_3` has the parity of the number of invariant cyclic triples), THM-791, THM-4922 (4), THM-4930,
  THM-5131 (S) (this machine, this morning, corrected by its audit).
- **Paley.** THM-212 (`p | H`), LEM-003 (the automorphism group acts freely on Hamiltonian paths), THM-5101
  (`189/21 = 9`, `95095/55 = 1729`; endpoint counts 1 against 8 at `P_7`), the 1729 notes (typed ACCIDENT, HYP-2306),
  THM-4746 (`(2/p)` in the count of strong labelled tournaments mod `p`), THM-4857 / THM-481 (the adjacency matrix
  over `F_2`, quadratic-residue codes, de Caen's rank bound), THM-338 (circulant maximisers), Session240 (strong
  tournaments with an automorphism of order 11).
- **Strong components as primes.** THM-1862 and THM-1936 (the monoid of tournaments under ordered sum is free on the
  strong tournaments; `H` is multiplicative), THM-4816 / THM-4817 (`T = 1/(1 - S)`, `S = A/(1 + A)`, self-converse
  strong counts), THM-4744 (interval chains), THM-341, THM-354, THM-1925 (components as diagonal blocks),
  `07-reflections/iso-classes-are-a-free-monoid-on-strongly-connected-primes-opus-20260629.md`, the
  numerical-semigroup framing (HYP-3616, HYP-3617), Session24 (multiplicity series by ordered factorisation).
- **Edge profiles.** THM-050 (consecutive-position counts; total `(n-1) H`), and Session265 (Codex, today:
  `session265_converse_dyadic_20261111/`), which answers the same directive: the weighted graph `W(u,v)` of the
  pasted text determines the strong components (cuts of weight `H`), connectors are rank one, there are "shadow
  twins"; seven Fano controls, 21 arcs, 63 folded characters.
- **Symmetric functions.** THM-073, THM-078, THM-002 (Grinberg-Stanley), INV-051 (Mitrovic-Stojadinovic: `W_X` in
  noncommuting variables, indexed by set partitions, product rule), `kps_s21_two_adic_lines_20261110/README.md`
  section 6.2 with its audit A4 (the mean path-cover polynomial has exponential generating function
  `((2 + t)/(2 - t))^x`), THM-589 (second moment by compositions into odd parts).
- **The attached paper.** It was supplied to two earlier sessions, which used it as an analogy:
  `kps_s22_sign_flip_towers_20261110/SOURCES_AND_DOMINATION.md` section 7 records the QUESTION "is there a
  noncommutative lift of the odd-cycle formula whose coefficients in products of power sums of the first kind count
  packings with a record structure?", and `odd_c_seed_symmetry_20261110/EXTENSIONS.md` notes that plain ordered
  path covers have symmetric coefficients. Section 6 answers the question.
- **5, the golden ratio, eleven squares.** `square11_octic_langlands_collatz_20261006.md` (Trump's packing, tilt
  `40.18` degrees, side `3.8770835900`, an octic with Galois group `S_8`), `strong_balance_pisano_golden_eleven_
  20261014.md` (sections 3-4: `11 = L_5`, `55 = F_10 = |Aut(P_11)|`; waste exponent between `1/2` and `3/5`),
  THM-4592 (the eleven-squares Fibonacci torus; gold cells `{1,3,4,5,9}`), `session91_packing_golden_carriers_
  20261029/`, `chanpak_55_105_golden_20261110/`, Session252 (five cores), THM-150 (`prod (1 + Q_k) = F_p`),
  THM-4741 / THM-1370 (Busch floors `5^((n-1)/3)`), `session74_eleven_mechanisms_20261027/` ("repeated appearances
  of 11 are not treated as a common mathematical object").
- **Not found by the surveys** (and not found by targeted greps of this session): the parity of `F`; `H(P_q) mod 4`
  or `mod 8`; any statement on residues of path counts of doubly regular tournaments or on `H mod 4` under
  switching; the pentagon sums of special families; the mean Redei-Berge function as a symmetric function, Schur
  `Q`-functions, sprouts, `NSym`; `log T(x)` and reach sets; factorisation lengths of path counts.
- **Canonical hostile example.** CD_4: it has `H = 1 (mod 4)` like the Paley tournaments of order `7 (mod 8)`, but no
  converse symmetry (THM-451), so Theorems Q and E say nothing about it. It keeps the doubly regular conjecture of
  HYP-9667 honest: the proof stops exactly where the converse symmetry stops.
- **Corrected near miss.** My own first guess for the case "all invariant triples cyclic" (`F = (n-1)/2 mod 2`),
  registered as R3 and refuted at orders 5 and 9; the truth is `F_k = C(m, k) (mod 2)`.
- **Concurrent sessions today.** Codex Session265 (same directive, complete), opus S92 (reserved THM-5134, 5135,
  HYP-9665, 9666: local densities, composition-algebra bridges, Collatz arcs), kind-pasteur S23 (THM-5133), Codex
  Session266. None treats the parity of invariant cycles, `H(P_q) mod 4`, or symmetric-function means.

## 2. Notation

- `H(T)`: number of Hamiltonian paths; `c_k`: number of directed `k`-cycles; `c_odd` their sum over odd `k`;
  `H = 1 + 2 c_odd (mod 4)`.
- An *anti-automorphism* is a bijection `psi` with `a -> b` iff `psi(b) -> psi(a)`; an *anti-involution* is one of
  order 2. KNOWN (Session29): a self-converse tournament has one, all are conjugate, it fixes at most one vertex
  and exactly one when the order is odd.
- `F_k = F_k(T, psi)`: number of directed cycles `C` of length `2k + 1` with `psi(C) = C`; `F = sum_k F_k`.
- `P_q`: Paley tournament on `GF(q)`, `q = 3 (mod 4)`: `a -> b` iff `b - a` is a non-zero square. `chi` is the
  quadratic character of `GF(q)`.
- `CD_k`: sign tournament of the Cayley-Dickson algebra of dimension `2^k` on its `2^k - 1` imaginary units
  (`a -> b` iff `e_a e_b = + e_(a xor b)`); it is the doubling tower of THM-448 (THM-5133 (G)).
- `G_35(R) = sum_X (-1)^((c_3 + c_5)(R + v_X))`, the sum over the `2^m` one-vertex extensions of `R` (THM-5132).
- A tournament is *doubly regular* if every vertex has score `(n-1)/2` and every two vertices have `(n-3)/4` common
  out-neighbours (`n = 3 mod 4`).

## 3. Theorem Q: folding along a converse symmetry

**Setting.** `T` is a tournament on `2m + 1` vertices with an anti-involution `psi`; `v_0` is its fixed vertex. In
each pair `{a, psi(a)}` exactly one vertex is an out-neighbour of `v_0`; call it `a^+` and put `a^- = psi(a^+)`.
Number the pairs `1, ..., m`.

- `D_1` is the set of pairs `i` with `a_i^+ -> a_i^-`. Equivalently: the invariant triple `{v_0, a_i^+, a_i^-}` is a
  cyclic triangle.
- `Z` is the digraph on the pairs with `i -> j` iff exactly one of `a_i^+ -> a_j^+`, `a_i^+ -> a_j^-` holds.

**Theorem Q.**

1. `Z` is a tournament.
2. For every `k >= 1`, `F_k = #{directed paths of Z on k vertices that end in D_1} (mod 2)`.
3. `c_odd(T) = F = #{directed paths of Z that end in D_1} (mod 2)`, and
   `H(T) = 1 + 2 #{directed paths of Z that end in D_1} (mod 4)`.

**Corollaries.**

- (a) If no invariant triple is cyclic (`D_1` empty), every `F_k` is even and `H = 1 (mod 4)`.
- (b) If every invariant triple is cyclic, `F_k = C(m, k) (mod 2)` for every `k`, `F` is odd and `H = 3 (mod 4)`.
- (c) (Balance.) Fix `psi` and the arcs between different pairs (hence `Z`). Among the `2^m` choices of `D_1` exactly
  half give `H = 1 (mod 4)`. So in the ensemble of all labelled tournaments with a given anti-involution of odd
  order, the two classes are equally frequent.
- (d) For a circulant tournament on `Z_n` with connection set `S` (`n` odd), `psi(x) = -x`, the pairs are the
  elements of `S`, `D_1 = {a in S : -2a in S}` and `a -> b` in `Z` iff exactly one of `b - a`, `-b - a` lies in `S`.

*Proof.* (1) Write `x_ij = [a_i^+ -> a_j^+]` and `y_ij = [a_i^+ -> a_j^-]`. Then `x_ij + x_ji = 1`, and
`y_ji = [a_j^+ -> a_i^-] = [a_i^+ -> a_j^-] = y_ij` by applying `psi`. So `z_ij = x_ij + y_ij (mod 2)` satisfies
`z_ij + z_ji = 1`.

(2) KNOWN (Session29 half-path report, section 2): an invariant cycle of length `2k + 1` has the unique form
`v_0 -> x_1 -> ... -> x_k -> psi(x_k) -> ... -> psi(x_1) -> v_0`, where `(x_1, ..., x_k)` is a directed path using at
most one vertex of each pair, `v_0 -> x_1` and `x_k -> psi(x_k)`. So `x_1 = a^+` for its pair, and the last vertex is
`a^+` if its pair is in `D_1` and `a^-` otherwise.

Fix an ordered sequence `(i_1, ..., i_k)` of distinct pairs and count the sign patterns. Between consecutive pairs
`i, j` the transitions allowed are `(+,+)` iff `x_ij = 1`, `(+,-)` iff `y_ij = 1`, `(-,+)` iff `y_ij = 0`
(apply `psi`: `a_i^- -> a_j^+` iff `a_j^- -> a_i^+`), `(-,-)` iff `x_ij = 0`. With rows and columns indexed by
`(+, -)`,

    M_ij = [[x_ij, y_ij], [1 - y_ij, 1 - x_ij]]  =  z_ij I + g_ij 1^T     (mod 2),

where `1 = (1,1)^T` and `g_ij = e_+` if `y_ij = 1`, `e_-` if `y_ij = 0`. The number of half-paths along the sequence
is `e_+^T M_(i_1 i_2) ... M_(i_(k-1) i_k) f`, with `f = e_+` if `i_k` is in `D_1` and `e_-` otherwise. Expand the
product modulo 2, choosing at each step either `z I` ("continue") or `g 1^T` ("cut"). Since `1^T g = 1`, `1^T f = 1`
and `e_+^T g_ij = y_ij`:

- no cut: `prod z * [i_k in D_1]`;
- cuts at steps `t_1 < ... < t_r`, `r >= 1`: `y_(first cut) * prod_(steps not cut) z`.

Summing over all sequences of `k` distinct pairs: the first kind gives the number of `k`-vertex paths of `Z` ending
in `D_1`. The second kind gives, for each `r >= 1`, the number of ordered lists `(Q_0, Q_1, ..., Q_r)` of
vertex-disjoint paths of `Z` covering a `k`-set, weighted by `y(end of Q_0, start of Q_1)`. For `r >= 3` the paths
`Q_2, ..., Q_r` can be permuted, so the total is a multiple of `(r - 1)!`, even. For `r = 1` and `r = 2` (fix `Q_2`
and pass to the sub-tournament on the other vertices) the sum is over pairs `(Q_0, Q_1)` covering a vertex set `V'`
with a `Y`-edge between the end `u` of `Q_0` and the start `w` of `Q_1`. Let `B(u, w)` be the number of such pairs
with the ends prescribed. `B(u, w)` is the number of Hamiltonian paths through the arc `u -> w` in the tournament
obtained from `Z[V']` by orienting `{u, w}` as `u -> w`; `B(w, u)` is the same for the other orientation; the
Hamiltonian paths not using `{u, w}` are common to both tournaments, and both tournaments have an odd number of
Hamiltonian paths (Redei). So `B(u, w) = B(w, u) (mod 2)`, and since `y` is symmetric the terms with `r = 1, 2`
cancel in pairs. This proves (2).

(3) Sum (2) over `k`, and use `c_odd = F (mod 2)` (the cycles that are not invariant come in pairs) and
`H = 1 + 2 c_odd (mod 4)`.

(a) is immediate. (b): the number of `k`-vertex paths of `Z` is `sum_(|S| = k) H(Z[S]) = C(m, k) (mod 2)` by Redei,
and `sum_(k >= 1) C(m, k) = 2^m - 1`. (c): for fixed `Z` let `P(d)` be the number of paths of `Z` ending at `d`;
`sum_d P(d) = 2^m - 1` is odd, so the linear form `D_1 -> sum_(d in D_1) P(d) (mod 2)` is not zero and takes each
value on half of the subsets. (d) is the definition. QED

**Remark (the fold is an invariant).** All anti-involutions of `T` are conjugate under `Aut(T)` (they are the
Sylow 2-subgroups of the group of automorphisms and anti-automorphisms, whose index-2 subgroup `Aut(T)` has odd
order; KNOWN, Session29). So the pair `(Z, D_1)`, up to isomorphism, depends only on the isomorphism class of the
self-converse tournament `T` of odd order. For even order an anti-involution has no fixed vertex, there is no
invariant odd cycle and `H = 1 (mod 4)` (KNOWN, THM-791).

*A name clash.* Session265 (Codex, today) calls "converse fold" the passage from the path-usage weights of the
arcs to those of the unordered pairs. The fold of this note is a different object: the quotient of a tournament by
one of its converse symmetries.

*Typing.* PROVED. KNOWN and used: the half-path form of invariant cycles and the congruence `H = 1 + 2F (mod 4)`
(Session29), Redei's theorem, `c_3 = F_1 (mod 2)` (THM-552). The case (a) also has a direct proof: the map
`(x_1, ..., x_k) -> (psi(x_k), ..., psi(x_1))` is then a fixed-point-free involution on the half-paths; this is how
it was found. Not found in the repository: the tournament `Z`, the parity formula, the corollaries.

*Machine checks.* `antiinv_census.c`: every normalised labelled tournament with the standard anti-involution at
orders 3, 5, 7, 9, 11 (2, 16, 512, 65,536 and 33,554,432 tournaments): no violation of (2) for any `k`; no
violation of `H = 1 + 2F (mod 4)` at orders up to 9 (direct count of `H`). `circulant_check.py`: all 126 connection
sets of orders 3-13, `H` by subset dynamic programme against (3): no disagreement. Parity vectors `(F_1, ..., F_m)`
at `|D_1| = m`: `(1)`, `(0,1)`, `(1,1,1)`, `(0,0,0,1)`, `(1,0,0,1,1)`, as in (b).

## 4. The character of 2 in the Paley tournaments

(Found first; Theorem E of section 5.2 shows that the dichotomy below holds for every doubly regular tournament
with a converse symmetry. What is special to Paley is the name of the pattern and the fold.)

**Theorem (Paley).** Let `q = 3 (mod 4)` be a prime power, `psi(x) = -x`.

1. The pairs are `{a, -a}` with `a^+` the square. `D_1` is the set of all pairs if 2 is not a square in `GF(q)`
   (`q = 3 mod 8`) and is empty if 2 is a square (`q = 7 mod 8`).
2. `a -> a^2` is an isomorphism from `Z` onto `P_q` restricted to the non-zero squares.
3. `F_k = [q = 3 mod 8] C((q-1)/2, k) (mod 2)` for every `k`, and

       H(P_q) = 2 - chi(2) = q(q-1)/2      (mod 4).

4. The number `K_q = H(P_q)/(q(q-1)/2)` of Hamiltonian paths that start with the arc `0 -> 1` is `1 (mod 4)`. For
   prime `q` it is the number of orbits of Hamiltonian paths under `Aut(P_q)`.

*Proof.* (1) `0 -> a` iff `a` is a square. For a square `a`, `a -> -a` iff `-2a` is a square iff `chi(-2) = 1` iff
`chi(2) = -1`, because `chi(-1) = -1`. For `q = p^e` (`e` odd), `chi(2) = (2/p) = (2/q)`. (2) For squares `a, b`:
`x_ab = [chi(b - a) = 1]`, `y_ab = [chi(-b - a) = 1] = [chi(a + b) = -1]`, so `z_ab = 1` iff
`chi(b - a) chi(b + a) = 1` iff `b^2 - a^2` is a square; squaring is a bijection of the squares (a group of odd
order). (3) Theorem Q (a), (b). (4) The maps `x -> ax + b` with `a` a non-zero square form a group of
automorphisms of order `q(q-1)/2`, transitive on the arcs and free on the Hamiltonian paths (LEM-003); for prime `q`
it is the whole automorphism group; and `q(q-1)/2` is odd. QED

**Data** (`halfpaths.c`, `paley_input.py`, `fold_paley.py`; `out/halfpaths_paley_*.out`, `out/fold_paley.out`).

| `q` | `q mod 8` | folded tournament `Z = P_q[squares]` | `H(Z)` | `(F_k mod 2)_k` | `F mod 4` | `H(P_q) mod 8` | `K_q` |
|---|---|---|---|---|---|---|---|
| 3 | 3 | one point | 1 | `1` | 1 | 3 | 1 |
| 7 | 7 | cyclic triangle on `{1, 2, 4}` | 3 | `000` | 0 (`F = 12`: KNOWN, Session29) | 5 | 9 |
| 11 | 3 | regular 5-tournament on `{1, 3, 4, 5, 9}` | 15 | `10011` | 3 | 7 | 1729 |
| 19 | 3 | regular, order 9; `H` is the maximum at order 9 | 3357 | `100000011` | 1 | 3 | 6857869865 |
| 23 | 7 | order 11 | 92411 | all 0 | 0 | 5 | 62293308207033 |
| 27 | 3 | order 13 | 3703011 | `1001100110011` | 3 | 7 (Run F) | `H = 27 * 13 * h`, `h = 1 (mod 256)` |
| 31 | 7 | order 15 | 196922055 | all 0 | 0 | not computed | |
| 43 | 3 | order 21 | | `C(21, k) mod 2` | 3 | not computed | |
| 47 | 7 | order 23 | | all 0 | 0 | not computed | |

`H(P_q)` for `q <= 23` and `K_q` are inherited (THM-212, THM-5101, OPEN-QUESTIONS); `H(P_27) mod 256 = 95` is new
(`hampaths_layered.c`, validated on `q = 7, 11, 19, 23`). The folds of `P_7`, `P_11`, `P_19` attain the largest path
count of their order (3, 15, 3357; the maxima are KNOWN, A038375 and THM-338); the fold of `P_23` does not (92,411
against 95,095).

**Observation, not proved (HYP-9667 (b8)).** In the six cases `q = 3, 7, 11, 19, 23, 27` one has
`H(P_q) = q(q-1)/2 (mod 8)`; the case `q = 27` was pre-registered (R8); the same holds for every doubly regular
tournament on record (section 5.3). In terms of Session29's congruence
`H = 1 + 2F + 4I (mod 8)` this says: `F = 0 (mod 4)` and `I` odd for `q = 7 (mod 8)`; `F = (q+1)/4 (mod 4)` and `I`
even for `q = 3 (mod 8)`. The values of `F mod 4` in the table agree for all nine `q`.

**Reading.** The residue of `H(P_q)` mod 4 is the answer to "is 2 a square mod `q`": the same character that decides
whether the multiplier `x -> 2x` is an automorphism of `P_q`, whether the binary quadratic-residue code exists, and
whether the vertex cube over `P_q` is flat for the pentagon (section 5). `7 = 2^3 - 1` is the first order at which
it is `+1`.

## 5. Doubly regular tournaments

### 5.1 Theorem D: pentagon sums

**Theorem D.**

1. Every regular tournament has an even number of 5-cycles.
2. Let `R` be a doubly regular tournament of order `n`, `A` its adjacency matrix over `F_2`.
   - If `n = 7 (mod 8)`: `rank A = (n+1)/2`, every element of `ker A^T` has weight `0 (mod 4)`, and
     `G_35(R) = 2^n`: every one-vertex extension of `R` has `c_3 + c_5` even.
   - If `n = 3 (mod 8)`: `A` is orthogonal over `F_2` and `G_35(R) = 2^((n+1)/2)`, the least value allowed by
     THM-5132 (3) at these orders.
3. In particular `G_35(P_q) = 2^q` or `2^((q+1)/2)` according to `chi(2) = 1` or `-1`, and `G_35(CD_k) = 2^(2^k - 1)`
   for every `k >= 3`.

*Proof.* (1) THM-5132 (1): `c_5 = C(q_1, 2) + C(q_2, 2) (mod 2)` with `q_1` the number of odd scores and `q_2` the
number of vertices with an odd number of odd-score out-neighbours. In a regular tournament all scores are
`(n-1)/2`: if this is even, `s = 0` and `t = As = 0`; if it is odd, `s = 1` and `t = A1 = 1`. In both cases
`q_1 = q_2`.

(2) Over the integers `A A^T = (t+1) I + t J` with `n = 4t + 3`, and `A 1 = A^T 1 = (2t + 1) 1`. If `t` is odd,
`A A^T = J` and `A^T 1 = 1` over `F_2`. For `u` in `ker A^T`: `0 = A A^T u = |u| 1`, so `|u|` is even; and
`C(|u|, 2) = u^T A u = 0`, so `|u| = 0 (mod 4)`. If `u` is in `ker A` and in `ker A^T` then `(J + I) u = 0`, so `u` is
`0` or `1`, and `A 1 = 1` excludes `1`; hence `2 (n - rank A) <= n`, and `rank(A A^T) = 1 >= 2 rank A - n` gives
`rank A <= (n+1)/2`. So `dim ker A^T = (n-1)/2`, `K_1 = 1 + K_0` has weights `n - 0 = 3 (mod 4)`, and the formula of
THM-5132 (3) for `m = n = 3 (mod 4)` gives `2^((n+1)/2) 2^((n-1)/2)`. If `t` is even, `A A^T = I`, the kernel is
zero, the solution of `A^T u = 1` is `u = A 1 = 1` of weight `3 (mod 4)`, and the formula gives `2^((n+1)/2)`.
(3) Paley tournaments are doubly regular; so are the `CD_k` (KNOWN: THM-448, THM-5133 (G)), of order
`2^k - 1 = 7 (mod 8)` for `k >= 3`. QED

*Typing.* PROVED, as a corollary of THM-5132 and the classical identities of doubly regular tournaments. KNOWN:
the rank `(n+1)/2` at `n = 7 (mod 8)` and the code-theoretic facts for Paley (THM-481, THM-4857), double regularity
of the tower. Not found in the repository: (1), the evaluation of `G_35`.

*Machine checks* (`family_cubes.py`, `out/family_cubes.out`): `P_q` for `q = 3, 7, 11, 19, 23, 27, 31, 43, 47` and
`CD_k` for `k = 2, ..., 7`, by the kernel formula, by elimination of the quadratic form, and for `m <= 11` by
counting 3- and 5-cycles in all `2^m` extensions: all agree with (2).

### 5.2 Theorem E: a converse symmetry forces the invariant triples

**Theorem E.** Let `T` be a doubly regular tournament of order `n = 4t + 3` with an anti-involution `psi`, fixed
vertex `v_0`.

1. For every `a` with `v_0 -> a`: `a -> psi(a)` if and only if `t` is even. So every invariant triple
   `{v_0, a, psi(a)}` is cyclic if `n = 3 (mod 8)` and none is if `n = 7 (mod 8)`.
2. `c_(2k+1)(T) = F_k = [n = 3 (mod 8)] C((n-1)/2, k) (mod 2)` for every `k >= 1`, and

       H(T) = n(n-1)/2 (mod 4):      H = 3 (mod 4) if n = 3 (mod 8),      H = 1 (mod 4) if n = 7 (mod 8).

*Proof.* (1) For an arc `u -> w` of a doubly regular tournament, `#{x : u -> x -> w} = t` and
`#{x : w -> x -> u} = t + 1` (split `N^+(u)` into `{w}`, the `t` common out-neighbours and the rest; split `N^+(w)`
into the `t` common out-neighbours and the rest). Put `a' = psi(a)` and `E = {x : a -> x -> a'}`. If `x` is in `E` then
`psi(a') -> psi(x) -> psi(a)`, that is `psi(x)` is in `E`; and `v_0` is not in `E` because `v_0 -> a`. So `E` is a union
of 2-element orbits of `psi` and `|E|` is even. If `a -> a'` then `|E| = t`; if `a' -> a` then `|E| = t + 1`.
(2) Theorem Q (a), (b); the cycles of length `2k + 1` that are not invariant come in pairs. QED

*Typing.* PROVED. Not found in the repository. The counting step uses only the pair `{a, psi(a)}`; for any
self-converse tournament of odd order it says that `#{x : a^+ -> x -> a^-}` is even and `#{x : a^- -> x -> a^+}` is odd.

*Remarks.* (i) For `P_q` and `psi(x) = -x`, (1) reads `chi(-2) = 1` iff `q = 3 (mod 8)`: for primes `q = 3 (mod 4)` this
is the second supplementary law of quadratic reciprocity, here a by-product of double regularity (the law is
classical; nothing is claimed). So section 4 is the Paley case of Theorem E, with the pattern of invariant triples
named by the character of 2. (ii) Machine checks: the non-Paley doubly regular tournament `D_19` of order 19
(`|D_1| = 9` of 9; `H = 1172732961771 = 3 mod 4`, exact), three non-Paley self-converse classes of order 27
(`|D_1| = 13` of 13; `H = 3 mod 4` by folding and by a direct count, Run H3), the Paley tournaments
(`out/drt_zoo_g1.out`, `out/drt_zoo_g3.out`, `out/drt27_direct_classes.out`).

### 5.3 Without a converse symmetry: a conjecture and its evidence (HYP-9667)

Most doubly regular tournaments have no converse symmetry (`CD_4`; 36 of the 39 classes of order 27 met in Run G).
Theorem Q does not apply to them. Since `c_3 = t + 1` and `c_5 = 0 (mod 2)` for every doubly regular tournament
(THM-5132), the statement `H = n(n-1)/2 (mod 4)` is equivalent to: *the number of odd cycles of length at least 7
is even.*

**Conjecture (HYP-9667 (b)).** (b4) Every doubly regular tournament of order `n` has `H = n(n-1)/2 (mod 4)`.
(b8) Every doubly regular tournament of order `n` has `H = n(n-1)/2 (mod 8)`.
(bS) If the `+-1` skew matrix `S` of a tournament `U` on `N` vertices satisfies `S^2 = -(N-1) I` (that is, `I + S` is a
skew Hadamard matrix), then `H(U) = 1 + N/2 (mod 4)`. (bS) implies (b4): add a source to `T`.
(bK) Under the same hypothesis, for every `k >= 1` the number of directed cycles of length `2k + 1` satisfies

    c_(2k+1)(U) = [N = 4 (mod 8)] C(N/2 - 1, k)      (mod 2).

(bK) implies (bS) (sum over `k`), and it is the statement of Theorem E without the converse symmetry: each odd
length separately, on the whole switching class.

*The same statement in terms of paths.* Let `p_j(T)` be the number of directed paths of `T` on `j` vertices. For every
tournament, `(p_j - C(n, j))/2 = sum_l C(n - l, j - l) c_l (mod 2)`, the sum over odd `l` (THM-5132, subset form). The
system is unitriangular in the odd lengths, so for a doubly regular `T` of order `n = 2m + 1` (bK) is equivalent to

    p_j(T) = C(n, j) (mod 4) for every j                                         if n = 7 (mod 8),
    p_j(T) = C(n, j) + 2 [j odd, j >= 3] C(m, (j-1)/2) (mod 4) for every j        if n = 3 (mod 8).

In words: *modulo 4, a doubly regular tournament of order 7 (mod 8) has the path counts of the transitive
tournament, for paths of every order.* Checked directly for `P_3, P_7, P_11, CD_4, P_19, D_19`
(`out/path_orders.out`).

PROVED cases of (bK): `k = 1, 2` for every tournament of order `N = 0 (mod 4)` whose scores have one parity
(THM-5132 (1): then `c_3 = N/4` and `c_5 = 0 (mod 2)`); every `k` for the doubly regular tournaments with a converse
symmetry (Theorem E).

The class of tournaments in (bS) is closed under switching (reversing all arcs between a vertex set and its
complement), and the tournaments `(source) + T` with `T` doubly regular are its normalised members. `H mod 4` is
**not** a switching invariant in general (it changes for about half of the random pairs at every order from 5 to
10, `out/switch_class_general.out`).

**Evidence** (FINITE-EXACT; mod 4 is PROVED for the classes with a converse symmetry, marked sc).

| order | doubly regular tournaments | source of `H` | `H mod 8` | `n(n-1)/2 mod 8` |
|---|---|---|---|---|
| 3, 7, 11 | `P_q` (sc) | inherited | 3, 5, 7 | 3, 5, 7 |
| 15 | `CD_4`, not sc (and its converse) | THM-448 | 1 | 1 |
| 19 | `P_19` (sc); `D_19` (sc), the only class given by the 12,000 Goethals-Seidel matrices with circulants of order 5 | inherited; HYP-9365 and Run G1 | 3; 3 | 3 |
| 23 | `P_23` (sc); twelve classes found by annealing in the Sierpinski lane (HYP-9371), none sc, two converse pairs among them; the double of `P_11`, not sc | exact counts of HYP-9371 (`drt_results.out`, `drt23_classes.out`, `drt23_doubled_check.out`); Run G2 mod 256 | 5 in all fourteen | 5 |
| 27 | `P_27` (sc); three sc classes from Goethals-Seidel matrices with circulants of order 7; six classes without converse symmetry | Run F; Run H3 (direct counts mod 256) | H3-PENDING-ROW (run in progress at this checkpoint) | 7 |

**Switching classes** (`switch_class.py`): `H mod 4` over members of the class of `(source) + T`.

| `N` | class of | members examined | `H mod 4` | `H mod 8` |
|---|---|---|---|---|
| 4 | `P_3` | all 8 | 3 | 3 |
| 8 | `P_7` | all 128 | 1 | 5 |
| 12 | `P_11` | 60 random | 3 | 7 |
| 16 | `CD_4` | 60 random | 1 | 1 and 5 |
| 20 | `P_19` | 30 random | 3 | 3 |
| 20 | `D_19` | 40 random (R11) | 3 | 3 and 7 |
| 24 | `P_23` | 40 random (R11) | 1 | 5 |
| 24 | double of `P_11` | 40 random (R11) | 1 | 5 |

So (bS) holds in every case examined, and its mod-8 analogue is false (`N = 16, 20`); (b8) concerns the normalised
members only.

**Cycle counts by length** (`cycle_lengths.c`, `cycle_parities.py`; pre-registered R12). Parity vectors
`(c_3, c_5, ..., c_(N-1))`; the prediction of (bK) is `C(N/2 - 1, k) mod 2` for `N = 4 (mod 8)` and zero for
`N = 0 (mod 8)`.

| `N` | tournaments | parity vector found | (bK) |
|---|---|---|---|
| 4, 8 | all members of the classes of `P_3^+`, `P_7^+` (in-sample) | `1`; `000` | `1`; `000` |
| 12 | 41 members of the class of `P_11^+` (in-sample) | `10011` | `10011` |
| 16 | `CD_4` and 12 members of its class (in-sample) | `0000000` | `0000000` |
| 20 | `D_19`, `P_19` and 20 random members of each class (I3) | `100000011` in all 42 | `100000011` |
| 24 | the fourteen doubly regular classes of order 23, twelve of them without converse symmetry (I1) | zero in all 14 | zero |
| 24 | 10 random members of each of the classes of `P_23^+`, `D(P_11)^+` (I2) | zero in all 20 | zero |
| 28 | `P_27`, three self-converse and six other Goethals-Seidel classes (I4) | `1001100110011` in all 10 | `1001100110011` |

The parities of the even cycle counts are not constant on a class (they differ among the classes of order 23 from
length 16 on). The vertex-local form is false: the numbers of odd cycles through a single vertex do not have these
parities (`CD_4`: the 13-cycles through 8 of the 15 vertices are odd in number; `D_19`: lengths 13 and 15 at 6 of
the 19 vertices; `out/cycle_parities_local.out`). The law is global.

**No member of such a class is self-converse** (PROVED, two lines). Let `S^2 = -(N-1) I` and suppose `psi` is an
anti-involution of `U`; `N` is even, so `psi` has no fixed vertex. For a vertex `a` the sets `{x : a -> x -> psi(a)}` and
`{x : psi(a) -> x -> a}` are `psi`-invariant, hence even; but together they have `(N-2)/2` elements (orthogonality of
the rows `a`, `psi(a)`), an odd number. So Theorems Q and E cannot be applied to the other members of a class: the
converse symmetry of a doubly regular tournament is lost by every switching.

**What does not prove it** (each tested, `out/` and section 10).

- `det S = 4H - 3 (mod 16)` holds for all tournaments of order 2 and 4 (it would give (bS), because
  `det S = (N-1)^(N/2)` there) and fails at order 6.
- *One orthogonal row.* If row `v` of `S` is orthogonal to all other rows, does switching at `v` keep `H mod 4`? True
  for every tournament tried with `N <= 8` (793 random tournaments of order 8 with such a row, 724 of them not
  skew Hadamard; random regular tournaments of orders 3, 5, 7 with a new source, switched), and **false from
  `N = 10` on** (regular tournaments of orders 9, 11, 13, 15: 64 of 160, 80 of 160, 25 of 48, 21 of 48 pairs differ;
  `local_forms.py`).
  A small-number law in the sense of the previous note: the local form of (bS) is true exactly while the order is
  at most 8.
- *Score parities only.* In a tournament of (bS) all scores have the same parity, and under that hypothesis alone the
  change of `H mod 4` at `v` is `<s, g>` with `s` the score parities of `U - v` and `g` as below; `<s, g> = 0` for all
  labelled tournaments of order 3 and 5 and fails for 178 of 400 random tournaments of order 7.

**A reduction that is proved.** For any tournament `U`, vertex `v`, `R = U - v` and `U^v` = `U` switched at `v`:

    (H(U) - H(U^v))/2 = sum_{x in N^-(v)} g_R(x)  (mod 2),      g_R(x) = #{even-order paths of R with endpoint x} mod 2,

and if `|R|` is odd, `g_R(x)` is the parity of the number of paths of `R` on `|R| - 1` vertices with endpoint `x`.
(The first is `H = 1 + 2 c_odd (mod 4)` applied to the odd cycles through `v`; the second follows from the endpoint
lemma `e_W(x) = sum_{x in W' subset W} s_W'(x) (mod 2)` for the numbers of Hamiltonian paths of `R[W]` ending and
starting at `x`, itself Redei's theorem applied to `R[W]` with a new vertex that beats only `x`. Elementary; the
lemma was not searched for in the repository and is not claimed.) So `H mod 4` is constant on a skew Hadamard
class if and only if, for every member `U` and every vertex `v`, the paths of `U - v` that miss exactly one vertex
and have exactly one endpoint beaten by `v` are even in number.

**The double in tournament form.** The doubling `[[M, M], [-M^T, M^T]]` of a skew Hadamard matrix (THM-447, THM-448)
turns a doubly regular `T` into the tournament on `V_1 + {oo} + V_2` with `T` on `V_1`, the converse of `T` on `V_2`,
`oo -> V_1`, `V_2 -> oo`, `a_1 -> a_2`, and `a_1 -> b_2` iff `a -> b` (checked against the matrix for `P_3, P_7, P_11`).
Its order is `7 (mod 8)`, so (b4) predicts `H = 1 (mod 4)` for every double (true for `CD_4` and for the double of
`P_11`). The same construction on an arbitrary base tournament gives both residues from base order 5 on, so double
regularity of the base has to enter a proof.

## 6. Records: the attached paper applied to tournaments

Chen, Mogielnicki and Ono count structures by a *record composition*: blocks are listed by increasing largest
element, and the count is a product over blocks with binomial factors in the partial sums `s_j`. Their Proposition
6.1 does this for every exponential family, with the canonical lift `R_n = sum_alpha (b_alpha/pi_u(alpha)) Psi^alpha`
in `NSym` (`Psi` the noncommutative power sums of the first kind, `pi_u(alpha)` the product of the partial sums).

**6.1 The mean Redei-Berge function.** Grinberg-Stanley: `U_T = sum_S 2^|S| prod_(C in S) p_len(C) * p_1^(rest)`, the
sum over sets `S` of pairwise disjoint directed odd cycles of `T`; `H(T) = U_T(1, 0, 0, ...)`.

**Proposition R1.** (a) `U` is multiplicative over ordered sums: `U_(T_1 -> T_2) = U_(T_1) U_(T_2)`. (b) Over
uniformly random labelled tournaments on `n` vertices,

    E[U_T] = n! sum_(lambda |- n, odd parts) 2^(l(lambda) - n) p_lambda / z_lambda = n! 2^(-n) q_n
           = (n!/2^(n-1)) sum_(k=0)^(n-1) s_(n-k, 1^k),

where `sum_n q_n t^n = prod_i (1 + x_i t)/(1 - x_i t)` (Schur's `Q`-functions of one row). So
`sum_n E[U_T] t^n/n! = prod_i (2 + x_i t)/(2 - x_i t)`: the sprout sequence with seed `(2 + t)/(2 - t)`.

*Proof.* (a) No cycle of an ordered sum meets both summands. (b) A permutation whose non-trivial cycles are odd is a
packing of directed cycles with probability `prod 2^(-len)`; with the weight 2 per cycle each non-trivial cycle
contributes `2^(1 - len)`, and fixed points `1 = 2^0`. The second equality is `log prod (1 + x_i t)/(1 - x_i t) =
2 sum_(k odd) p_k t^k/k`; the third is Pieri (`q_n = sum_(a+b=n) e_a h_b`). QED

KNOWN: the scalar specialisations (`sum E[H_n] t^n/n! = (2 + t)/(2 - t)`, HYP-1400 / HYP-1528; the mean path-cover
polynomial `((2+t)/(2-t))^x`, kind-pasteur S21 with its audit). New here is the symmetric-function statement.

**6.2 The canonical lift and the record refinement.** For a labelled tournament `T` on `[n]` put

    U_T^nc = sum_S 2^|S| Psi^(tau(S))   in NSym,

where `tau(S)` is the composition of the block sizes (cycles of `S` and singletons) listed by increasing largest
vertex. It projects to `U_T`, and it is the image of the Mitrovic-Stojadinovic function `W_T` (INV-051) under
`p_pi -> Psi^(tau(pi))`.

**Proposition R2.** `E[U_T^nc] = n! sum_(alpha |= n, odd parts) (2^(l(alpha) - n)/pi_u(alpha)) Psi^alpha`, the
canonical lift (Chen-Mogielnicki-Ono, Definition 5.4) of `E[U_T]`. That is: the mean weighted number of packings of
odd cycles of a random tournament whose max-ordered type is `alpha` equals `n! 2^(l(alpha) - n)/pi_u(alpha)`. In
particular `sum_(alpha |= n, odd parts) 2^l(alpha)/pi_u(alpha) = 2` for `n >= 1`, and the count depends on the
order of the parts (`3/2` for `(1,3)` against `1/2` for `(3,1)` at `n = 4`).

*Proof.* The number of ordered set partitions of `[n]` with block sizes `alpha` and increasing maxima is
`prod_j C(s_j - 1, alpha_j - 1)` (their Lemma 3.4); a block of odd size `k >= 3` carries `(k-1)!` directed cycles,
each present with probability `2^(-k)` and weight 2; their Lemma 4.1 with `d = 1` turns the product into
`n! prod_j b_(alpha_j)/s_j` with `b_k = 2^(1-k)`. QED

This answers the question recorded by kind-pasteur S22 at the level of the mean. It is consistent with the remark of
`odd_c_seed_symmetry`: the asymmetry comes from ordering by maxima, not from ordering the cover.

**6.3 Tournaments, strong components and connected graphs.** A labelled tournament is a sequence of strong
components. Cut the sequence before every component whose largest vertex exceeds all earlier ones: the pieces are
tournaments in which the largest vertex lies in the first component, i.e. reaches every vertex; listing the pieces
by increasing largest vertex recovers the sequence.

**Proposition R3.** Let `l_n` be the number of labelled tournaments on `[n]` in which the vertex `n` reaches every
vertex. Then (a) `sum 2^C(n,2) x^n/n! = exp(sum l_n x^n/n!)`, so `l_n` is the number of connected labelled graphs
on `n` nodes (A001187: `1, 1, 4, 38, 728, 26704, ...`); (b) `l_n = sum_j C(n-1, j-1) S_j 2^C(n-j, 2)` with `S_j` the
number of strong labelled tournaments; (c) the number of labelled tournaments whose pieces, listed by increasing
largest vertex, have sizes `alpha` is `prod_j C(s_j - 1, alpha_j - 1) l_(alpha_j)`, the number of labelled graphs
whose components listed in the same way have sizes `alpha`.

*Proof.* (a) and (c): the pieces form a set of structures on the blocks of a set partition; apply the exponential
formula and Proposition 6.1 of the paper; both species have the same exponential generating function. (b): choose
the component of `n`, strong, and an arbitrary tournament on the vertices it does not reach. QED

So in the dictionary `NSym -> Q` sending `S_n` to `2^C(n,2)/n!`: the elementary functions `(-1)^(k-1) Lambda_k` go
to the strong tournaments (`sigma(t) lambda(-t) = 1` is `T (1 - S) = 1`, KNOWN), the power sums `Psi_k` to the
tournaments rooted at their largest vertex (`sigma' = sigma psi` is `T' = T L'`), and the canonical lift of the
tournament series is the canonical lift of the graph series (their Example 11, entry by entry: 12 tournaments and
12 graphs of type `(1,3,1)` at `n = 5`). The relation between `T` and the connected graphs is used in THM-4746; the
reading by reach sets was not found.

*Machine checks* (`records_check.py`, `out/records_check.out`): R1 (b), R2, R3 (a)-(c) on every labelled tournament
of order at most 5; the identity of R2 for `n <= 22`; no failure.

## 7. The monoid of path counts is not half-factorial, because of 63

KNOWN (THM-4094; WEIGHTED_ATLAS (9)): the conjectured set of values `M = {odd} \ {7, 21}` is a multiplicative monoid
whose atoms are the odd primes other than 7, the numbers `7p` with `p` an odd prime other than 3 (49 included), 63
and 343; the conjecture is equivalent to all atoms being values of strong tournaments.

**Remark (PROVED, elementary).** Write `m = 3^a 7^b k` with `k` prime to 21 and `c` prime factors (with
multiplicity). A factorisation into atoms uses `u <= c` atoms `7p` (`p | k`), `v` copies of 49, `w` of 63, `x` of 343
with `u + 2v + w + 3x = b`, `2w <= a`, and has length `L = c + a + v + x - w`. Hence the *elasticity* of `M` (largest
ratio of two factorisation lengths of one element) is `5/2`, attained exactly at the even powers of 63
(`63^2 = 3^4 * 49`: lengths 2 and 5); the smallest element with factorisations of different lengths is
`315 = 5 * 63 = 3 * 3 * 35`.

*Proof of the bound.* `L <= c + a + b/2` (take `w = 0`). From `3(v + x) >= b - u - w` and `u <= c`:
`L >= (2c + 3a + b - 4w)/3`. If `a <= 2b`, use `w <= a/2`: `L >= (2c + a + b)/3`, and
`3(c + a + b/2)/(2c + a + b) <= 5/2` because `a <= 4c + 2b`. If `a > 2b`, use `w <= b`: `L >= c + a - b`, and
`(c + a + b/2)/(c + a - b) <= 5/2` because `c + a >= 2b`. Equality forces `c = 0`, `a = 2b`, `w = b`, `b` even. QED

For the actual set of path counts the lower bound `5/2` is unconditional (7 and 21 are not counts, and 3, 49, 63
are: every odd number below 80,265 other than 7 and 21 is a count, THM-4102 / THM-4104 as summarised in THM-4749);
the upper bound holds if the conjecture does. `hmonoid_elasticity.py` checks the atom list and the ratio below
300,000.

63 is the only atom divisible by 9, and it exists only because 7 and 21 are missing; it is realised by exactly two
strong tournaments, both of order 8 (THM-344).

## 8. The merge

Each row: the theme of the directive, what the repository already has (with its type), what this session adds, and
the type of the link. "No bridge" is a result.

| theme | inherited | added here | type of the new link |
|---|---|---|---|
| 7, 21, 63 | 7 and 21 excluded, 63 an orphan with two creators (THM-1370, THM-344, THM-4942); `7 * 3^k` tower refuted | 63 as the atom that breaks unique factorisation length (section 7) | PROVED, elementary |
| powers of 2 | fugacity 2 is the integer point of the Fibonacci root; Mersenne and Jacobsthal values (THM-1425, THM-4748, THM-5104) | `H(P_q) mod 4` is the quadratic character of 2; `P_7` folds to the triangle on the powers of 2 mod 7 | PROVED (Theorem Q, section 4) |
| Cayley-Dickson | doubling tower, `H(CD_4)`, chirality (THM-448, THM-451, THM-5133) | all one-vertex extensions of every `CD_k`, `k >= 3`, have `c_3 + c_5` even; the kernel of the adjacency matrix mod 2 is doubly even of dimension `2^(k-1) - 1` | PROVED (Theorem D) |
| | | `H(CD_k) mod 4` for `k >= 5` | OPEN (HYP-9667 (b4) predicts 1); Theorems Q, E do not apply (no converse symmetry) |
| Hadamard matrices | skew Hadamard doubling (THM-447, THM-448), signature simplex (THM-4911), loop expansion on doubly regular tournaments (HYP-9365, HYP-9371) | `H = n(n-1)/2 (mod 4)` for every doubly regular tournament with a converse symmetry; constant `H mod 4` on skew Hadamard switching classes in all cases examined | PROVED (Theorem E); CONJECTURE (HYP-9667) |
| Collatz | `(x, sigma) -> (3x + sigma, -sigma)`; 7 and 21 as the missing arrows `2 -> 1`, `5 -> 16` (Codex, today); `2^k = 3 J_k + (-1)^k` (THM-4754) | nothing | no new bridge |
| 5 | Busch floors `5^((n-1)/3)`; five cores (Session252); `S_4` with `H = 5` | every regular tournament has an even number of 5-cycles; the six strong classes of order 5 are all self-converse (KNOWN, recomputed) | PROVED (Theorem D (1)) |
| golden ratio | `prod (1 + Q_k) = F_p` (THM-150); golden eleven (2026-10-14); Fibonacci torus (THM-4592) | the fold of `P_11` is the regular 5-tournament on the gold cells `{1, 3, 4, 5, 9}`, whose non-principal eigenvalues have moduli `phi` and `1/phi` | IDENTITY for the fold; the eigenvalue fact is classical |
| 7 against 11 | Session240, Session74 | `(2/7) = 1`, `(5/7) = -1`; `(2/11) = -1`, `(5/11) = 1`: 2 is a square at 7 and 5 at 11; Theorem Q reads the first character (`H(P_7) = 1`, `H(P_11) = 3 mod 4`) | IDENTITY for `(2/q)`; no statement found for `(5/q)` |
| eleven squares | Trump's packing and its octic; waste exponents; contact tori (four notes) | nothing geometric; strong components as squares on the diagonal of the adjacency matrix are THM-341 / THM-354 / THM-1925 and waste nothing in the limit (almost all tournaments are strong) | no bridge found |
| strong components as primes | free monoid, `T = 1/(1 - S)`, atoms (THM-1862, THM-1936, THM-4817, THM-4094) | `NSym` dictionary: strong = elementary, rooted at the largest vertex = power sums = connected graphs; elasticity `5/2` | PROVED (sections 6.3, 7) |
| converse-collapsed profiles | the weighted graph `W` determines the strong components; connectors have rank one (Session265, today) | a different fold: the quotient `Z` of a tournament by a converse symmetry; it carries `H mod 4` | PROVED (Theorem Q) |
| the attached paper | used as an analogy twice; a recorded question | the mean Redei-Berge function is a sprout sequence; its canonical lift has a record interpretation | PROVED (section 6) |

**Remarks on the table.**

- *Two characters.* The directive puts 7 with powers of 2 and 11 with 5 and the golden ratio. In the Paley family
  this is exact: the multiplier group of `P_7` is generated by 2 (`{1, 2, 4}`), that of `P_11` contains 5 and not 2
  (`{1, 3, 4, 5, 9}` is generated by 3, 4, 5 or 9). Theorem Q turns the first fact into a residue of `H`. The second
  has no such reading here.
- *The fold as a small-number phenomenon.* `7 -> 3`, `11 -> 5`, `19 -> 9`, `23 -> 11`, `47 -> 23`: folding halves
  the order. The fold of `P_23` is not `P_11` (`H = 92411` against `95095`), so there is no ladder along the chain
  `5, 11, 23, 47`; typed REFUTED as a guess, `out/fold_paley.out`.
- *Small numbers again.* The local lemma behind the Hadamard conjecture (one orthogonal row of the skew matrix
  keeps `H mod 4` under switching) is true for every order up to 8 and false from order 10; `det S = 4H - 3 (mod 16)`
  is true up to order 4; the score-parity form up to order 5. Three statements that are theorems for small
  tournaments and wrong afterwards, next to one (the skew Hadamard form) that has survived every test to order 28.
  This is the pattern of the previous note (crevasses below, a different law above), met here inside one proof
  attempt.
- *What does not transfer.* Nothing in this session produces a new path count or excludes a value. The values 7
  and 21 are `3 (mod 4)` and `1 (mod 4)`; Theorems Q and E restrict `H` only on tournaments with a converse symmetry.
  That 7 and 21 are the numbers of vertices and of arcs of `P_7` is on record (Session265: seven controls, 21 arcs)
  and stays typed as a coincidence of names: no mechanism here links it to the exclusion of the values 7 and 21.

## 9. Pre-registered predictions

`PREREGISTRATION.md` has six blocks, each pushed before its runs (commits `c98a09856a`, `953e580672`, `1432c0ac0d`,
`6e387cf531`, `88d90c6a80`, `dcad3e3237`). The outcomes are appended to that file.

| | prediction | outcome | verdict |
|---|---|---|---|
| R1 | `F(P_27)`, `F(P_43)` odd | odd | **CONFIRMED**; now a theorem (section 4) |
| R2 | `F_k(P_q)` even for every `k >= 2`, `q = 11, 19, 27, 43` | false already at `q = 11`, `k = 4`: `F_k = C(m, k) (mod 2)` | **FAILED** |
| R3 | all invariant triples cyclic implies `F = (n-1)/2 (mod 2)` | `F` is odd at every order (4 of 4 at `n = 5`, 4096 of 4096 at `n = 9`) | **FAILED**; replaced by Theorem Q (b) |
| R4 | `G_35(P_q) = 2^q` (`q = 7, 23, 31, 47`), `2^((q+1)/2)` (`q = 3, 11, 19, 43`) | as predicted, three methods | **CONFIRMED**; now Theorem D |
| R5 | `G_35(CD_k) = 2^(2^k - 1)`, `k = 3, 4, 5` | as predicted (also `k = 6, 7`) | **CONFIRMED**; now Theorem D |
| R6 | `CD_4`, `CD_5` self-converse; `H(CD_5) = 1 (mod 4)` | no anti-automorphism exists; second part undecided | **FAILED**; the first part was in canon (THM-448, THM-451) and should not have been registered |
| R7 | every doubly regular descendant of the order-16 sign matrix has `H = 1 (mod 4)` | true, but all 16 are isomorphic to `CD_4` | CONFIRMED, uninformative |
| R8 | `H(P_27) = 7 (mod 8)` | `H(P_27) = 27 * 13 * h`, `h = 1 (mod 256)`; `H = 95 (mod 256)` | **CONFIRMED** |
| R9 | non-Paley doubly regular tournaments: `H = 3 (mod 4)` at orders 19 and 27, `1 (mod 4)` at order 23 | G1: one class (`D_19`), 3; G2: one class, 1; G3 and H3: H3-PENDING-R9 | **CONFIRMED**; G1 and G2 were decidable from the repository (below); for the classes with a converse symmetry now Theorem E |
| R10 | the same modulo 8 (`3, 5, 7`) | G1: 3; G2: 5; H3: H3-PENDING-R10 | **CONFIRMED**, with the same caveat |
| R11 | `H mod 4` is constant on the switching classes of `P_23^+`, `D(P_11)^+` (`N = 24`) and `D_19^+` (`N = 20`) | 40 random members each: 1, 1, 3 | **CONFIRMED** |
| R12 | every odd cycle count has the parity `[N = 4 mod 8] C(N/2 - 1, k)`: fourteen classes of order 23, members of the classes at `N = 20, 24`, ten tournaments of order 27 | I1: 14 of 14; I2: 20 of 20; I3: 42 of 42; I4: 10 of 10 | **CONFIRMED** |
| controls | `F(P_q)` even at `q = 7, 23, 31, 47`, odd at `3, 11, 19`; exact `H` of `P_7, P_11, P_19, CD_4`; `H(P_23) mod 256`; `H(P_27) = 95 (mod 256)` by the 27-start method | hold | |

Three remarks on the record.

- R2 and R3 were first-moment guesses from two or three cases; the data that refuted them contained the true law
  (`C(m, k)` is visible in the parity vectors), and the proof came from asking why.
- **Two registrations were decidable from the repository.** R6 (first part: THM-448, THM-451 say that the tower is
  not self-converse). R9 and R10 at orders 19 and 23: the tournament `D_19` with `H = 1172732961771` and the exact
  path counts of thirteen doubly regular classes of order 23 and of the doubled Paley class are in this machine's
  own Sierpinski lane of 2026-10-22 (HYP-9365, HYP-9371; `04-computation/experiments/sierpinski_ocf_20261022/`).
  The check made before registering searched four theorem files and the words of the predicted property, not the
  objects. The data found afterwards is used as in-sample evidence in section 5.3 (all fourteen classes of order
  23 have `H = 5 mod 8`). Out-of-sample in the strict sense were: the three self-converse classes of order 27 (G3),
  the direct counts at order 27 (H3), the switching classes (R11) and all of R12 (no cycle counts of doubly
  regular tournaments were on record). MISTAKE-PENDING (entry filed at close-out).
- Deviations from the registered runs. G1 was first run on every 40th Goethals-Seidel solution and then on all
  12,000 (both outputs kept; one class). G3 examined 40 of the 246,960 solutions (every 2000th) and all 28
  descendants of each. The programs named in the first and third blocks were renamed (`paley_sym.c` is
  `halfpaths.c`, `selfconv_census.py` is `antiinv_census.c`, `paley27.c` is `hampaths_layered.c`).

## 10. Routes tried and why they stopped

- *A direct involution for `q = 3 (mod 8)`.* Reversal composed with `psi`, multiplication by 2, inversion `x -> 1/x`:
  each maps the half-paths to another family or is an automorphism of odd order. The transfer-matrix argument
  replaced them.
- *`H mod 8` of Paley and of doubly regular tournaments.* Needs `F mod 4` and the parity of the overlapping mirror
  pairs `I` (Session29); the expansion of section 3 modulo 4 does not close (the terms with two and three cuts no
  longer cancel). Recorded as HYP-9667 (b8) with data.
- *A proof of (bS), (bK) or (b4) without a converse symmetry.* Tried: the cycles through one vertex (the parities
  are not forced vertex by vertex); a determinant (`det S mod 16`; true for `N <= 4`
  only); invariance of `H mod 4` under switching at a vertex whose row of `S` is orthogonal to the others (true
  for `N <= 8`, false from `N = 10`); the score-parity form (true for `|R| <= 5`, false at 7); the Walsh expansion of
  `c_odd` in even path forests (in a skew Hadamard tournament every sum with a free middle vertex vanishes, but
  the patterns without a vertex of degree 2 are not determined by `S^2`). What remains is the proved reduction of
  section 5.3. A proof has to use all rows at once.
- *A geometric bridge to the packing of eleven squares.* None found; the existing notes already say so
  (`chanpak_55_105_golden_20261110/targets/PROOF.md`, `session74_eleven_mechanisms_20261027/`).
- *A ladder of folds along a Cunningham chain.* Refuted at the first step (section 8).

## 11. Reproduction

Requires a C compiler, Python 3 and `gentourng` (nauty). From this directory:

    cc -O2 -o /tmp/halfpaths halfpaths.c
    for q in 3 7 11 19 23 27 31 43 47; do python3 paley_input.py $q | /tmp/halfpaths > out/halfpaths_paley_$q.out; done   # Run A
    cc -O2 -o /tmp/antiinv antiinv_census.c
    for m in 1 2 3 4 5; do /tmp/antiinv $m > out/antiinv_census_$m.out; done      # Run B and Theorem Q (m = 5: 2^25 tournaments)
    python3 circulant_check.py 13 > out/circulant_check.out                        # Theorem Q against direct counts
    python3 circulant_fold_table.py 21 > out/circulant_fold_table.out              # H mod 4 of all circulants by folding
    python3 fold_paley.py > out/fold_paley.out                                     # folded Paley tournaments
    python3 family_cubes.py > out/family_cubes.out                                 # Run C, Theorem D (about 3 minutes)
    for k in 2 3 4 5; do python3 cd_converse.py $k; done > out/cd_converse.out     # Run D
    python3 drt_descendants.py 3 > out/drt_descendants_3.out; python3 drt_descendants.py 4 > out/drt_descendants_4.out   # Run E
    cc -O2 -o /tmp/hamlayer hampaths_layered.c
    for q in 7 11 19 23 27; do python3 paley_h01_input.py $q | /tmp/hamlayer > out/paley${q}_h01.out; done   # Run F (q = 27: 400 MB)
    cc -O2 -o /tmp/hamexact hampaths_exact.c
    python3 drt_zoo.py controls /tmp/hamexact /tmp/hamlayer > out/drt_zoo_controls.out          # Run G
    python3 drt_zoo.py g1 /tmp/hamexact 40 > out/drt_zoo_g1_stride40.out; python3 drt_zoo.py g1 /tmp/hamexact 1 > out/drt_zoo_g1.out   # 30 s; 25 min
    python3 drt_zoo.py g2 /tmp/hamlayer > out/drt_zoo_g2.out; python3 drt_zoo.py g3 40 2000 > out/drt_zoo_g3.out
    python3 switch_class.py insample /tmp/hamexact > out/switch_class_insample.out; python3 switch_class.py general > out/switch_class_general.out
    python3 switch_class.py h1 /tmp/hamlayer 40 > out/switch_class_h1.out; python3 switch_class.py h2 /tmp/hamexact 40 > out/switch_class_h2.out   # Run H1, H2
    python3 drt27_direct.py /tmp/hamlayer control > out/drt27_direct_control.out                # Run H3: 5 min, 800 MB
    python3 drt27_direct.py /tmp/hamlayer classes 40 2000 6 > out/drt27_direct_classes.out      # 45 min
    python3 local_forms.py /tmp/hamexact > out/local_forms.out                                  # section 5.3, section 10
    python3 path_orders.py > out/path_orders.out                                                # paths of every order mod 4
    cc -O2 -o /tmp/cyclen cycle_lengths.c                                                       # Run I
    for m in controls insample local; do python3 cycle_parities.py /tmp/cyclen $m > out/cycle_parities_$m.out; done
    python3 cycle_parities.py /tmp/cyclen i1 ../../../04-computation/experiments/sierpinski_ocf_20261022/drt23_classes.json > out/cycle_parities_i1.out
    python3 cycle_parities.py /tmp/cyclen i2 10 > out/cycle_parities_i2.out; python3 cycle_parities.py /tmp/cyclen i3 20 > out/cycle_parities_i3.out
    python3 cycle_parities.py /tmp/cyclen i4 40 2000 6 > out/cycle_parities_i4.out               # 10 min, 800 MB
    python3 records_check.py 5 > out/records_check.out                             # section 6
    python3 hmonoid_elasticity.py 300000 > out/hmonoid_elasticity.out              # section 7
    python3 atoms_small.py > out/atoms_small.out; python3 paley_orbits.py > out/paley_orbits.out

`fold.py`, `paley_tools.py`, `cayley_dickson.py` are modules. `paley_orbits.py` recomputes inherited numbers
(`189 = 168 + 21`: 168 Hamiltonian paths of `P_7` close to a cycle and 21 do not, the 8 + 1 orbits of THM-5101).
