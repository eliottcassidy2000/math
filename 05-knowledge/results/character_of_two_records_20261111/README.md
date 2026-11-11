# Folding a tournament by its converse: the character of 2 in the path count of Paley tournaments, pentagon sums of doubly regular tournaments, and records

**Session:** mac-mini-2026-11-11-atoms. **Canon:** THM-5136 (folding formula, Paley), THM-5137 (pentagon sums of doubly
regular tournaments; records). **Hypotheses:** HYP-9667.
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
mathematics in two places where the lanes were empty, and it types the rest.

1. **(PROVED) Folding formula.** Let `T` have an anti-automorphism of order 2 (a converse symmetry); for odd order
   it fixes one vertex `v_0` and pairs the others. Folding `T` along it gives a tournament `Z` on the pairs and a set
   `D_1` of pairs (those whose triangle with `v_0` is cyclic), and

       H(T) = 1 + 2 * #{directed paths of Z that end in D_1}      (mod 4).

   More precisely, the number of invariant odd cycles of length `2k + 1` has the parity of the number of
   `k`-vertex paths of `Z` ending in `D_1` (Theorem Q). The repository had the congruence `H = 1 + 2F (mod 4)` with
   `F` the number of invariant odd cycles (Session29) and the remark that `F` is not decided by triangles alone;
   Theorem Q computes `F mod 2` on a tournament of half the size.
2. **(PROVED) The character of 2.** For the Paley tournament `P_q` (`q = 3 mod 4` a prime power), `D_1` is everything
   or nothing according to whether 2 is a non-square or a square in `GF(q)`, and `Z` is `P_q` restricted to the
   squares. Hence

       H(P_q) = 1 (mod 4) if q = 7 (mod 8),        H(P_q) = 3 (mod 4) if q = 3 (mod 8),

   that is `H(P_q) = q(q-1)/2 (mod 4)`; for prime `q` the number of orbits of Hamiltonian paths under the
   automorphism group (`1, 9, 1729, ...`) is `1 (mod 4)`. `P_7` folds to the cyclic triangle on `{1, 2, 4}` (the
   powers of 2, the octonion line), `P_11` to the regular pentagon tournament on `{1, 3, 4, 5, 9}` (the gold cells of
   the eleven-squares torus, THM-4592), `P_19` to a regular tournament of order 9 with 3357 paths, the largest count
   at that order.
3. **(PROVED) Pentagon sums of doubly regular tournaments.** By the kernel formula of THM-5132, every doubly regular
   tournament of order `n = 7 (mod 8)` (the Paley tournaments with `q = 7 mod 8` and all Cayley-Dickson sign
   tournaments on `2^k - 1 >= 7` units) has **all** `2^n` one-vertex extensions with `c_3 + c_5` even; for
   `n = 3 (mod 8)` the surplus is the least possible, `2^((n+1)/2)`. Every regular tournament has an even number of
   5-cycles (Theorem D).
4. **(PROVED, application of the attached paper) Records.** The mean of the Redei-Berge symmetric function over
   all tournaments on `n` vertices is `n! 2^(-n) q_n`, Schur's one-row `Q`-function, i.e. `n!/2^(n-1)` times the sum
   of all hook Schur functions; it is the sprout sequence with seed `(2 + t)/(2 - t)`. Its canonical noncommutative
   lift in the sense of Chen-Mogielnicki-Ono has the record interpretation asked for in a recorded question
   (kind-pasteur S22): the mean weighted number of odd-cycle packings whose blocks, listed by increasing largest
   vertex, have sizes `alpha` is `n! 2^(l(alpha) - n)/pi_u(alpha)`. Labelled tournaments in which the largest vertex
   reaches every vertex are counted by the connected labelled graphs (A001187), and tournaments and graphs are
   equidistributed by max-ordered type (section 6).
5. **(Pre-registered.)** Eight predictions: R1, R4, R5, R8 CONFIRMED; R7 confirmed but uninformative; R2, R3, R6
   FAILED. R2 and R3 were replaced by the law that was then proved (Theorem Q); R6 was decidable from canon and
   should not have been registered (section 9).
6. **(Typed merge.)** Section 8 places 5 and the golden ratio, the eleven squares, 7-21-63, powers of 2,
   Cayley-Dickson, Collatz and the "primes" text against the existing lanes, with the few exact links that are new
   and an explicit "no bridge" where none was found.

Nothing here enlarges the range of verified values or proves a new case of the conjecture. What is new and proved
concerns `H mod 4` again: on tournaments with a converse symmetry it is now a path count on the folded tournament.

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
  or `mod 8`; the pentagon sums of special families; the mean Redei-Berge function as a symmetric function, Schur
  `Q`-functions, sprouts, `NSym`; `log T(x)` and reach sets; factorisation lengths of path counts.
- **Canonical hostile example.** CD_4: it has `H = 1 (mod 4)` like the Paley tournaments of order `7 (mod 8)`, but no
  converse symmetry (THM-451), so Theorem Q says nothing about it. It keeps the doubly regular question of HYP-9667
  honest.
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

**Observation, not proved (HYP-9667 (a)).** In the six cases `q = 3, 7, 11, 19, 23, 27` one has
`H(P_q) = q(q-1)/2 (mod 8)`; the case `q = 27` was pre-registered (R8). In terms of Session29's congruence
`H = 1 + 2F + 4I (mod 8)` this says: `F = 0 (mod 4)` and `I` odd for `q = 7 (mod 8)`; `F = (q+1)/4 (mod 4)` and `I`
even for `q = 3 (mod 8)`. The values of `F mod 4` in the table agree for all nine `q`.

**Reading.** The residue of `H(P_q)` mod 4 is the answer to "is 2 a square mod `q`": the same character that decides
whether the multiplier `x -> 2x` is an automorphism of `P_q`, whether the binary quadratic-residue code exists, and
whether the vertex cube over `P_q` is flat for the pentagon (section 5). `7 = 2^3 - 1` is the first order at which
it is `+1`.

## 5. Theorem D: pentagon sums of doubly regular tournaments

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

**Question (HYP-9667 (b)).** Does every doubly regular tournament of order `n` have `H = 1 (mod 4)` for
`n = 7 (mod 8)` and `H = 3 (mod 4)` for `n = 3 (mod 8)`? True for the Paley tournaments (section 4) and for `CD_4`
(`H = 198,335,025`), which has no converse symmetry. The 16 doubly regular tournaments obtained from the order-16
Cayley-Dickson sign matrix are all isomorphic to `CD_4` (`drt_descendants.py`), so this is two families, not a
census.

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
| | | `H(CD_k) mod 4` for `k >= 5` | OPEN (HYP-9667 (b)); Theorem Q does not apply (no converse symmetry) |
| Collatz | `(x, sigma) -> (3x + sigma, -sigma)`; 7 and 21 as the missing arrows `2 -> 1`, `5 -> 16` (Codex, today); `2^k = 3 J_k + (-1)^k` (THM-4754) | nothing | no new bridge |
| 5 | Busch floors `5^((n-1)/3)`; five cores (Session252); `S_4` with `H = 5` | every regular tournament has an even number of 5-cycles; the six strong classes of order 5 are all self-converse (KNOWN, recomputed) | PROVED (Theorem D (1)) |
| golden ratio | `prod (1 + Q_k) = F_p` (THM-150); golden eleven (2026-10-14); Fibonacci torus (THM-4592) | the fold of `P_11` is the regular 5-tournament on the gold cells `{1, 3, 4, 5, 9}`, whose non-principal eigenvalues have moduli `phi` and `1/phi` | IDENTITY for the fold; the eigenvalue fact is classical |
| 7 against 11 | Session240, Session74 | `(2/7) = 1`, `(5/7) = -1`; `(2/11) = -1`, `(5/11) = 1`: 2 is a square at 7 and 5 at 11; Theorem Q reads the first character (`H(P_7) = 1`, `H(P_11) = 3 mod 4`) | IDENTITY for `(2/q)`; no statement found for `(5/q)` |
| eleven squares | Trump's packing and its octic; waste exponents; contact tori (four notes) | nothing geometric; strong components as squares on the diagonal of the adjacency matrix are THM-341 / THM-354 / THM-1925 and waste nothing in the limit (almost all tournaments are strong) | no bridge found |
| strong components as primes | free monoid, `T = 1/(1 - S)`, atoms (THM-1862, THM-1936, THM-4817, THM-4094) | `NSym` dictionary: strong = elementary, rooted at the largest vertex = power sums = connected graphs; elasticity `5/2` | PROVED (sections 6.3, 7) |
| converse-collapsed profiles | the weighted graph `W` determines the strong components; connectors have rank one (Session265, today) | a different fold: the quotient `Z` of a tournament by a converse symmetry; it carries `H mod 4` | PROVED (Theorem Q) |
| the attached paper | used as an analogy twice; a recorded question | the mean Redei-Berge function is a sprout sequence; its canonical lift has a record interpretation | PROVED (section 6) |

**Three remarks on the table.**

- *Two characters.* The directive puts 7 with powers of 2 and 11 with 5 and the golden ratio. In the Paley family
  this is exact: the multiplier group of `P_7` is generated by 2 (`{1, 2, 4}`), that of `P_11` contains 5 and not 2
  (`{1, 3, 4, 5, 9}` is generated by 3, 4, 5 or 9). Theorem Q turns the first fact into a residue of `H`. The second
  has no such reading here.
- *The fold as a small-number phenomenon.* `7 -> 3`, `11 -> 5`, `19 -> 9`, `23 -> 11`, `47 -> 23`: folding halves
  the order. The fold of `P_23` is not `P_11` (`H = 92411` against `95095`), so there is no ladder along the chain
  `5, 11, 23, 47`; typed REFUTED as a guess, `out/fold_paley.out`.
- *What does not transfer.* Nothing in this session produces a new path count or excludes a value. The values 7
  and 21 are `3 (mod 4)` and `1 (mod 4)`; Theorem Q restricts `H` only on tournaments that have a converse symmetry.

## 9. Pre-registered predictions

`PREREGISTRATION.md` has three blocks, each pushed before its runs (commits `c98a09856a`, `953e580672`,
`1432c0ac0d`).

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
| controls | `F(P_q)` even at `q = 7, 23, 31, 47`, odd at `3, 11, 19` | hold | |

The two failures R2 and R3 were first-moment guesses from two or three cases; the data that refuted them contained
the true law (`C(m, k)` is visible in the parity vectors), and the proof came from asking why.

## 10. Routes tried and why they stopped

- *A direct involution for `q = 3 (mod 8)`.* Reversal composed with `psi`, multiplication by 2, inversion `x -> 1/x`:
  each maps the half-paths to another family or is an automorphism of odd order. The transfer-matrix argument
  replaced them.
- *`H mod 8` of Paley tournaments.* Needs `F mod 4` and the parity of the overlapping mirror pairs `I`; the
  expansion of section 3 modulo 4 does not close (the terms with two and three cuts no longer cancel). Recorded as
  HYP-9667 (a) with data.
- *A census of doubly regular tournaments.* Only the Paley and the Cayley-Dickson families were available at
  orders where `H` can be counted; the descendants of one Hadamard matrix are one class.
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
    python3 records_check.py 5 > out/records_check.out                             # section 6
    python3 hmonoid_elasticity.py 300000 > out/hmonoid_elasticity.out              # section 7
    python3 atoms_small.py > out/atoms_small.out; python3 paley_orbits.py > out/paley_orbits.out

`fold.py`, `paley_tools.py`, `cayley_dickson.py` are modules. `paley_orbits.py` recomputes inherited numbers
(`189 = 168 + 21`: 168 Hamiltonian paths of `P_7` close to a cycle and 21 do not, the 8 + 1 orbits of THM-5101).
