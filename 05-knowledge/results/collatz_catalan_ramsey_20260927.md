# Catalan numbers, the Ramsey rows and the Collatz spine: the critical spine blocks are Dyck paths, Catalan parity lives at the Mersenne tower orders, the Paley rows of both Ramsey tables, an exact law for the Seidel doubling of graphs, and the orbit as a permutation tournament

**Session:** opus, `collatz-posets-zeta5-20260927` (seventh note), 2026-09-27.
**Owner's directive:** "merge in the fact that 1, 2, 5, 14, 42 counts binary
trees, and how those numbers connect to our previous work; work upcoming
directions as they emerge; go back and look at our work on the Ramsey
numbers, particularly R(5,5), and see how extensions can be made and ideas
combined."
**Inherits (cited):** THM-4495 (no-descent counts `W_k` of order `2^(h* k)
k^(-3/2)`; the ladder decomposition `W = 1/(1 - P)` into positive min-ending
words), the parallel spine note (tight spine blocks with Beatty lengths),
the fourth note (the family dichotomy at `q <= 3` versus `q >= 5`),
THM-447/448 (the skew-Sylvester doubling `D(T)` and the Mersenne tower),
THM-455 (transitive subtournaments of the tower: `3, 5, 7, 11` at orders
`7, 15, 31, 63`; the tournament Ramsey numbers `R(2..6) = 2, 4, 8, 14, 28`,
`34 <= R(7) <= 47`; Momihara–Suda's `trans(QR_q)`), THM-483 (the zigzag law
`trans(D(T)) = z(T)`, refuting the `+2` sandwich), THM-464 (`R_b(1) = R(3, b)`,
the classical row `R(3,5) = 14`), HYP-9162 (the Sierpinski tournament tower,
orders `2^k - 1`), S10 (Gilbreath's sea kernel and Pascal mod 2). Classical
facts from memory: Erdős–Moser 1964, Reid–Parker 1970, Sanchez-Flores 1994/98,
Kalbfleisch's `R(6,6) >= 102` from the Paley graph `P_101`, Exoo 1989 and
Angeltveit–McKay 2024 for `R(5,5)` (bounds verified today on the Wikipedia
table), Kummer's theorem, Erdős's 1979 conjecture on the ternary digits of
`2^k`.

**Status: PROVED elementary (Propositions 1–3, Theorem 4, Proposition 5) +
FINITE-EXACT (ladder counts to length 15; Catalan parity to `n = 4000` and
the ternary condition to `k = 200`; Paley clique numbers to `q = 101` and
Paley tournament trans to `q = 31`; the doubling law exhaustively for `n <=
6` and on random graphs for `n = 7, 8`; split numbers of `P_5..P_37`) +
CITED + NUMEROLOGY typed + DIRECTION. The repo had no `R(5,5)` work; its
Ramsey work is the tournament thread, which this note extends to graphs.
Collatz OPEN.** Script
`04-computation/experiments/collatz_catalan_ramsey_20260927.py`, output
beside it. Independent audit OWED (as for the session's other notes).

## 0. The answer in one paragraph

The Catalan numbers are the spine-block counts of the critical member of the
Collatz family: at multiplier `q = 4` (drift zero, slope one) the positive
min-ending words of THM-4495's ladder decomposition are a rise followed by a
Dyck path, so there are `C_n` of them at length `2n + 1` (`1, 2, 5, 14, 42,
132, 429` at lengths `3, 5, ..., 15`), the positive words are the central
binomials, and `W = 1/(1 - P)` holds verbatim (Proposition 1). At that
critical point the fourth note's certification density still tends to one,
but only like `1 - 1/sqrt(2 pi k)`; and the `k^(-3/2)` in THM-4495's `W_k ~
2^(h* k) k^(-3/2)` is the `n^(-3/2)` of `C_n ~ 4^n n^(-3/2)/sqrt(pi)`, the
universal first-passage exponent, with the base `4` replaced by `2^(h*)` because
the Collatz slope is irrational and the drift negative. Catalan parity lives
at the tower's orders: `C_n` is odd iff `n = 2^k - 1` (the Mersenne orders of
the Sierpinski tournament tower), and `C_n` is prime to `6` iff moreover `2^k -
1` has no ternary digit `2`, which below `k = 200` happens only for `k = 1, 2,
5, 8` (Proposition 2; the shifted form of Erdős's ternary-digit conjecture on
`2^k`). The repo's Ramsey work is the tournament thread (THM-455/483), not
`R(5,5)`; the classical facts are `43 <= R(5,5) <= 46` (Exoo 1989;
Angeltveit–McKay 2024). The Paley rows of both tables are computed here:
Paley graphs give `R(3,3) > 5`, `R(4,4) > 17`, `R(5,5) > 37`, `R(6,6) > 101`
(`P_41` already has clique number 5, so Exoo's 42-vertex graph is not
Paley), Paley tournaments give `R(4) > 7`, `R(5) > 11`, `R(6) > 23`, `R(8) >
31`. The extension: the Seidel doubling `D(G)` of a graph (the analogue of
THM-447's skew-Sylvester doubling) obeys an exact law, `omega(D(G))` is the
largest induced complete split subgraph of `G` and `alpha(D(G)) =
max(alpha(G) + 1, sep(G))` with `sep` the largest induced clique-plus-
independent-set with no edges between (Theorem 4, the graph analogue of the
zigzag law), verified exhaustively to `n = 6`; the doubled Paley graphs give
`R(5,5) > 26` and `R(6,6) > 58`, far from extremal, because Paley graphs
contain large split subgraphs; the increment of `max(omega, alpha)` under
doubling is `0, 1` or `2` on every graph with at most 6 vertices and reaches
`3` at `n = 9`, so the graph sandwich fails as the tournament one did. The
Catalan values in the Ramsey tables (`R(3,5) = R(5) = 14 = C_4`, Exoo's `42 =
C_5`) are coincidences; the exact Catalan–Ramsey link is Erdős–Szekeres
(`123`-avoiding permutations are counted by `C_n`), and the orbit of `27`,
read as a permutation tournament, has largest transitive subtournament `18`
(Proposition 5).

## 1. Catalan numbers are the critical spine blocks (PROVED + FINITE-EXACT)

In the family `T_q(n) = n/2, (qn+1)/2` the multiplier after `k` steps with
`o` odd letters is `q^o/2^k`, and the no-descent (positive) condition
`q^(o_j) > 2^j` for all `j` is a lattice-path condition with steps `log_2 q
- 1` and `-1`. At `q = 4` the steps are `+1, -1`.

**Proposition 1.** At `q = 4`: (a) the positive words of length `k` (all
prefix sums positive) number `C(k-1, floor((k-1)/2))`; (b) the positive
min-ending words (THM-4495's ladder blocks: positive, ending at the minimum
of the prefix sums, which is `1`) exist only at odd lengths `2n + 1` and
number `C_n = (1/(n+1)) C(2n, n)`, a rise followed by a Dyck path of length
`2n`; (c) `W = 1/(1 - P)` holds; (d) the fraction of words of length `k` not
certified by coefficient descent is `C(k-1, floor((k-1)/2))/2^k ~ 1/sqrt(2
pi k)`.

*Proof.* (a) is the ballot theorem; (b) a positive word ends at its minimum
`1` iff the path after its first step is a Dyck path returning to level 1;
(c) is THM-4495's decomposition, whose two factors are here the Catalan and
central-binomial generating functions (`1/(1 - t C(t^2)) = sum C(k-1,
floor((k-1)/2)) t^k`); (d) Stirling. Checked to length 15 (script part 1). ∎

Reading. `1, 2, 5, 14, 42` are the numbers of tight spine blocks of lengths
`3, 5, 7, 9, 11` for the critical map, where the spine blocks of the parallel
note (positive min-ending words) are exactly the Dyck excursions. At the
Collatz slope `log_2 3` the blocks have the Beatty lengths of the parallel
note's Theorem 3 and are counted by THM-4495's `B_n`, and the law `W_k ~
2^(h* k) k^(-3/2)` keeps the Catalan exponent `-3/2` while the base drops
from `4` (all `4^n` words of length `2n` at slope one... i.e. `2` per step) to
`2^(h*) = 1.93` per step, because the Collatz walk has drift `-0.415` and only
an exponentially small fraction of words stay positive. The fourth note's
dichotomy (certification density `-> 1` iff `q <= 3`) has `q = 4` as its
boundary, where the density still tends to one but polynomially: the
Catalan regime is the critical case between Terras (`q = 3`, exponential) and
the positive-density undecided sets (`q >= 5`).

## 2. Catalan parity and the tower orders (PROVED + FINITE-EXACT)

**Proposition 2.** (a) `C_n` is odd iff `n = 2^k - 1`. (b) `C_n` is coprime
to `6` iff `n = 2^k - 1` and `2^k - 1` has no digit `2` in base 3 (i.e. is a
sum of distinct powers of 3). Below `n = 4000` this holds for `n = 0, 1, 3, 31,
255`, and below `k = 200` the admissible `k` are `1, 2, 5, 8`.

*Proof.* Kummer: `v_2(C(2n, n))` is the number of carries in `n + n` in base
2, i.e. `s_2(n)`, so `v_2(C_n) = s_2(n) - v_2(n+1) = s_2(n+1) - 1`, zero iff
`n + 1 = 2^k`. Likewise `v_3(C(2n,n))` is the number of base-3 carries of
`n + n`, which is the number of ternary digits `2` of `n`, and `v_3(n+1) = 0`
for `n + 1 = 2^k`. Checked directly to `n = 4000` and `k = 200` (script part
2). ∎

Reading. The odd Catalan numbers sit exactly at the orders `3, 7, 15, 31,
63, ...` of the Sierpinski tournament tower (HYP-9162; THM-447/448), the
same Pascal-mod-2 mechanism as Gilbreath's sea kernel (S10): `C_n = C(2n,n) -
C(2n, n+1)` and Pascal's triangle mod 2 is Sierpinski's. The coprime-to-6
condition is the shifted form of Erdős's conjecture that `2^k` has a ternary
digit `2` for every `k > 8`: here it is `2^k - 1`, with the known small
solutions `k = 1, 2, 5, 8` (`1, 3, 31 = 27 + 3 + 1, 255 = 243 + 9 + 3`) and no
others below 200. Whether the list is complete is an Erdős-type question
this note does not touch; the repo's Collatz thread already records that
Erdős's problem needs the top ternary digits of `2^k` and that the Collatz
endgame reads low digits instead (synthesis, wave 3).

## 3. The Ramsey rows (CITED + FINITE-EXACT; coincidences typed)

**What the repo has.** No `R(5,5)` work. The Ramsey thread is the tournament
Ramsey function `R(t)` (least `n` forcing a transitive subtournament `TT_t`):
`R(2..6) = 2, 4, 8, 14, 28`, `34 <= R(7) <= 47` (Erdős–Moser, Reid–Parker,
Sanchez-Flores, Neiman–Mackey–Heule; THM-455), the tower values `trans(T_7,
T_15, T_31, T_63) = 3, 5, 7, 11` with the extremality window closing between
31 and 63, and the exact zigzag law `trans(D(T)) = z(T)` (THM-483), which
refuted the `+2` sandwich by the alternating family `A_l`; THM-464 places
the classical row `R(3, b)` (`3, 6, 9, 14, ...`) as the height-1 row of the
subgrid games.

**The classical facts.** `43 <= R(5,5) <= 46` (Exoo 1989; Angeltveit–McKay
2024; verified today), with McKay–Radziszowski–Exoo's 656 `(5,5,42)`-graphs,
none extendable, behind the conjecture `R(5,5) = 43`; `R(6,6) >= 102`.

**The Paley rows (script part 3).** Clique = independence number of the
Paley graph `P_q`, `q = 1 mod 4`:

| `q` | 5 | 13 | 17 | 29 | 37 | 41 | 53 | 61 | 73 | 89 | 97 | 101 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `omega(P_q)` | 2 | 3 | 3 | 4 | 4 | 5 | 5 | 5 | 5 | 5 | 6 | 5 |

so `R(3,3) > 5`, `R(4,4) > 17` (both sharp), `R(5,5) > 37` (Exoo's 42-vertex
witness is not Paley: `P_41` already contains `K_5`), `R(6,6) > 101`
(Kalbfleisch's bound; `P_97` has `omega = 6`). Largest transitive
subtournament of the Paley tournament `P_q`, `q = 3 mod 4`: `3, 4, 5, 5, 7`
for `q = 7, 11, 19, 23, 31` (matching Momihara–Suda and THM-455), so `R(4) >
7` (sharp), `R(5) > 11`, `R(6) > 23`, `R(8) > 31`; the extremal objects for
`R(5) = 14` and `R(6) = 28` are the sporadic circulant `ST_13` and the
quadratic-residue tournament on `GF(27)`, not prime Paley tournaments.

**Coincidences (NUMEROLOGY).** `R(3,5) = 14 = C_4 = R(5)` and Exoo's `42 =
C_5` vertices. The exact Catalan–Ramsey link is elsewhere: Erdős–Szekeres
(a sequence of length `(r-1)(s-1) + 1` has an increasing `r`-subsequence or a
decreasing `s`-subsequence) is the Ramsey theorem for permutation
tournaments, and the permutations of `[n]` with no increasing 3-subsequence
are counted by `C_n` (Knuth): the Catalan numbers count the extremal
witnesses of the smallest Erdős–Szekeres instance, not Ramsey numbers.

## 4. The Seidel doubling of graphs and its exact law (PROVED + FINITE-EXACT)

THM-447 doubles a tournament through its skew Seidel matrix `M`, `M' = [[M,
M + I], [M - I, -M]]`. The graph analogue: for a graph `G` on `n` vertices
with Seidel matrix `S` (`S_ij = -1` if adjacent, `+1` if not, `0` on the
diagonal), let `D(G)` be the graph on `2n` vertices with Seidel matrix

    S' = [[ S,      S + I ],
          [ S + I,   -S   ]],

i.e. the plain copy carries `G`, the prime copy carries the complement of
`G`, cross pairs `(u, v')` with `u != v` are adjacent iff `uv` is an edge of
`G`, and twins `(u, u')` are never adjacent.

**Theorem 4.** For every graph `G`:
(a) `omega(D(G)) = cs(G)`, the largest `|K| + |I|` over a clique `K` and an
independent set `I` of `G`, disjoint, with every vertex of `K` adjacent to
every vertex of `I` (the largest induced complete split subgraph `K_a ∨
bar K_b`);
(b) `alpha(D(G)) = max(alpha(G) + 1, sep(G))`, where `sep(G)` is the largest
`|I| + |K|` over an independent set `I` and a clique `K`, disjoint, with no
edge between them (the largest induced `bar K_a + K_b`).
Consequently `alpha(D(G)) >= alpha(G) + 1`, `omega(D(G)) >= max(omega(G),
alpha(G))`, and `cs(G) = sep(G)` for self-complementary `G`.

*Proof.* (a) A clique of `D(G)` consists of plains, which form a clique `K`
of `G`, and primes, which form a clique of the complement, i.e. an
independent set `I` of `G`; twins are never adjacent, so `K` and `I` are
disjoint; the cross rule makes every plain adjacent to every prime iff
`K × I ⊆ E(G)`. (b) An independent set of `D(G)` consists of an independent
`I` on the plain copy and a clique `K` on the prime copy; a cross pair `(u,
w')` with `u != w` is non-adjacent iff `uw` is not an edge, and twins are
always non-adjacent. If `u in I ∩ K` and `|K| >= 2`, then `u` and any other
`w in K` are adjacent in `G` (clique) yet `(u, w')` must be non-adjacent:
impossible; so either `I ∩ K = ∅` (giving `sep(G)`) or `K = {u} ⊆ I` (giving
`alpha(G) + 1`). Verified exhaustively for all labeled graphs on `3..6`
vertices (`33864` graphs) and on 300 random graphs on 7 and 120 on 8 vertices
(script part 6). ∎

This is the graph face of THM-483's zigzag law (there a transitive chain of
`D(T)` is a shuffle of an ascending chain of `T` and a descending one; here a
clique of `D(G)` is a clique of `G` joined to an independent set). The
`max(omega, alpha)` increment under doubling is `0, 1` or `2` on every graph
with at most 6 vertices (exhaustive: `{0: 3016, 1: 27282, 2: 2470}` at `n =
6`) and reaches `3` on random graphs at `n = 9`, so the graph analogue of
the `+2` sandwich fails as the tournament one did, one size later. Ramsey
use: `D(P_5), D(P_13), D(P_17), D(P_29), D(P_37)` have `max(omega, alpha) =
3, 4, 5, 5, 7` on `10, 26, 34, 58, 74` vertices, giving `R(4,4) > 10`, `R(5,5) >
26`, `R(6,6) > 34`, `R(6,6) > 58`, `R(8,8) > 74`: far below the records,
because the Seidel double of `G` contains every complete split subgraph of
`G` as a clique and Paley graphs have `cs(P_q) = 3, 4, 5, 5, 7` for `q = 5,
13, 17, 29, 37`. The tower idea (THM-455) leaves extremality between orders
31 and 63; the graph doubling never reaches it. What the two exact laws
share is the mechanism: a doubling that carries a copy and a reversed
(complemented) copy cannot suppress the mixed substructures, and the mixed
substructures are exactly the objects the Ramsey question counts.

## 5. The orbit as a permutation tournament (PROVED, elementary; FINITE-EXACT)

The value-time poset of an orbit (parallel note: `i` below `j` iff `i < j`
and `x_i > x_j`) is a permutation poset; the associated permutation
tournament (`i -> j` iff `i < j` and `x_i < x_j`, else `j -> i`) has
transitive subtournaments equal to monotone subsequences of the value
sequence, so its largest transitive subtournament is `max(LIS, LDS)`, and
Erdős–Szekeres guarantees a monotone subsequence of length `ceil(sqrt(L))`
in `L` values. For the orbits of `27, 703, 6171, 77031` (`42, 63, 97, 130` odd
values): `(LIS, LDS) = (18, 9), (18, 17), (17, 21), (18, 32)`, largest
transitive subtournaments `18, 18, 21, 32`, against the guaranteed `7, 8, 10,
12`. A divergent orbit's permutation tournament would have unbounded
transitive subtournaments (its records form an increasing subsequence); the
tournament Ramsey numbers add nothing to Erdős–Szekeres here.

## 6. Directions (DIRECTION; none pursued)

* **D18. Iterated doubling of graphs.** Theorem 4 gives `omega(D(G))` in
  terms of `G`; iterating, `omega(D^k(G))` is a `k`-fold nested split
  structure. Compute `D^k(P_5)` (`10, 20, 40, 80` vertices) and compare with
  the tower's `3, 5, 7, 11`; a closed form for `cs(D(G))` in terms of `cs(G)`,
  `sep(G)` and `omega, alpha` would be the graph analogue of THM-455's
  tower reading.
* **D19. Catalan numbers prime to 6.** Are `k = 1, 2, 5, 8` the only `k` with
  `2^k - 1` a sum of distinct powers of 3? This is the `-1` shift of Erdős's
  ternary conjecture; the repo's Erdős face (wave 3) may have the tools to
  say whether the two questions are equivalent in difficulty.
* **D20. The critical map as a model.** At `q = 4` every statement of the
  Collatz thread (spine blocks, ladder decomposition, certification
  density, shadow prices) has an exact Catalan/ballot form; the critical map
  is not integer-valued, but its combinatorics is the reference point
  against which the Collatz slope's irrationality (Beatty lengths, `h* < 1`,
  drift) is measured. A short dictionary "Collatz statement / its `q = 4`
  Catalan form" would make the thread's asymptotics legible at a glance.

## 7. Reproduction

    cd 04-computation/experiments
    python3 collatz_catalan_ramsey_20260927.py > collatz_catalan_ramsey_20260927.out

Standard library only; about eight minutes (the exhaustive `n = 6` doubling
census and the Paley clique searches).
