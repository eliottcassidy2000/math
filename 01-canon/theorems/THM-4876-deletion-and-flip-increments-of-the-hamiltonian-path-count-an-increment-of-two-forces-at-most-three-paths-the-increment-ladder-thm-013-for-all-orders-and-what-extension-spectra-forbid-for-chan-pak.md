---
id: THM-4876
title: "Deletion and flip increments of the Hamiltonian-path count H. (1) For a tournament T and a vertex v, H(T) - H(T - v) = 2 sum over the odd directed cycles C through v of H(T - V(C)) (Claim A, KNOWN from the odd-cycle formula), and for a tournament X with a new vertex v of cut (A, B), H(X + v) = sum_S e_B(S) s_A(V - S) + s_A(V) + e_B(V) (path pairs). (2) H(T) - H(T - v) = 2 iff v lies on exactly one odd cycle, a triangle (v, a, b), and T - {v, a, b} is transitive; then T - v has at most one odd cycle, H(T - v) in {1, 3} and H(T) in {3, 5}: no deletion lowers a count of at least 5 by exactly 2. (3) H(T) - H(T - v) = 4 iff v lies on exactly two odd cycles, both with transitive complements (and H(T - v) in {1, 5} through order 9). (4) The arc-flip increment is the unsigned subset convolution H(T') - H(T) = U(L_j, R_i) - U(L_i, R_j) whose signed companion is THM-017's identity, and THM-013's Delta H = Delta I(Omega, 2) holds for all n as a corollary of the odd-cycle formula H = I(Omega, 2). (5) FINITE-EXACT: the one-vertex extension spectra E(X) of all tournaments of orders 3-9 have no value h + 2 for h >= 5, smallest positive increments k_min(h) = 2, 3, 4, 5, 5, 5, 6, 7, 9, 9, 8, 10, 9, 11, 12, 11 for h = 5, 9, 11, ..., 39 (about h/3.5), and windows of at most 24 consecutive odd values, never covering [2h + 1, 3h]: no induction that extends a fixed realiser by one vertex can prove Chan-Pak. (6) The one-centre family TT_m + u has H = 1 + 2 sum_{a in A, b notin A, a < b} c(b - a), c(1) = 1, c(d) = 2^(d-2), hence contains 2^k + 1 and (2^i - 1)(2^l - 1) + 2 for all i, l >= 1 and omits 7, 13, 21, 29, ..."
status: >
  PROVED (mac-mini-2026-11-04-chanpak, elementary): statements 2, 3 (the structural halves), 4, 6, and the path-pair
  formula in 1. KNOWN/CITED: Claim A and the odd-cycle formula H(T) = I(Omega(T), 2) (Grinberg-Stanley, arXiv:2307.05569;
  Lass; the repository's THM-070 derives Claim A from it), Redei's theorem (H odd), Moon (every vertex of a strong
  tournament lies on a 3-cycle), THM-017 (signed subset convolution, PROVED for all n in March 2026), THM-013 (the
  same identity proved for n <= 7 in March 2026, now a corollary for all n), THM-012b (insertion decomposition),
  THM-4873 (end-count lemma, opus-S77). FINITE-EXACT: statement 5 and the value lists in 3 and 6
  (chanpak_spine_20261104/ext_spectrum_true.py, ext_windows.py, increment_ladder.py, increment2_order9.py,
  spine_numeration.py; formulas brute-force checked: path pairs on 468 (X, cut) pairs of orders <= 5, Claim A at v on
  every cut with increment <= 8 of orders <= 8, the session124 spine transfer on 69,174 triples). Not independently
  audited yet. The Chan-Pak conjecture itself (every odd number other than 7 and 21 is a Hamiltonian-path count) remains
  OPEN; this theorem closes one route (statement 5) and sharpens the local structure of the spectrum.
session: mac-mini-2026-11-04-chanpak
source: 05-knowledge/results/chanpak_increments_thm013_thm017_merge_20261104.md
scripts:
  - 04-computation/experiments/chanpak_spine_20261104/ext_spectrum_true.py, ext_windows.py, increment_ladder.py (+ _n8.out), increment2_order9.py (+ .out), spine_numeration.py, two_centre_family.py
related:
  - 01-canon/theorems/THM-013-arc-flip-delta-I.md (upgraded by statement 4), THM-017-even-odd-split-proof.md, THM-016-claim-b-alternating-sum.md, THM-012b-insertion-decomposition.md, THM-070-claim-a-from-ocf.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md, THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md, THM-4862-*, THM-4865-*, THM-4873-*
  - 05-knowledge/hypotheses/HYP-9469-*, HYP-9453-*, HYP-9420-*, HYP-9029-*
  - Chan-Pak, Computational complexity of counting coincidences, Theoret. Comput. Sci. 1015 (2024) 114776, arXiv:2308.10214, Conjecture 6.2; Grinberg-Stanley, The Redei-Berge symmetric function of a directed graph, arXiv:2307.05569
---

# THM-4876 — increments of the Hamiltonian-path count

## Setting

`H(T)` is the number of directed Hamiltonian paths of a tournament `T`; it is odd (Rédei). For a vertex `v`, `T - v`
is the deletion. For a tournament `X` and a new vertex `v`, the *cut* `(A, B)` is the partition of `V(X)` into the
out-neighbours `A` (`v -> A`) and in-neighbours `B` (`B -> v`) of `v`; `X + v` is the extension. For `U <= V(X)` and
`u in U`, `e_U(u)` and `s_U(u)` count the Hamiltonian paths of `X[U]` ending, resp. starting, at `u`;
`e_B(S) = sum_{u in S cap B} e_S(u)`, `s_A(S) = sum_{u in S cap A} s_S(u)`. The *odd-cycle formula* (Lass;
Grinberg-Stanley) is `H(T) = I(Omega(T), 2) = sum over families F of pairwise vertex-disjoint odd directed cycles of
2^|F|`. The one-vertex *extension spectrum* is `E(X) = {H(X + v) : all 2^|X| cuts}`; `min E(X) = H(X)` (a dominating
`v`), and `(H(X + v) - H(X))/2` is the *increment* of the cut.

## Statements

1. **Increment at a vertex (KNOWN) and the path-pair formula (PROVED).** For every tournament `T` and vertex `v`,

       H(T) - H(T - v) = 2 sum_{C odd cycle, v in C} H(T - V(C)),            (Claim A)

   and for every `X` and cut `(A, B)`,

       H(X + v) = sum_{S nonempty proper} e_B(S) s_A(V - S) + s_A(V) + e_B(V).

2. **Increment two (PROVED).** `H(T) - H(T - v) = 2` iff `v` lies on exactly one odd cycle, that cycle is a triangle
   `(v, a, b)`, and `T - {v, a, b}` is transitive. In that case `T - v` contains at most one odd cycle, so
   `H(T - v) in {1, 3}` and `H(T) in {3, 5}`. Consequently a vertex deletion never lowers a count `H(T) >= 7` by
   exactly `2`, and `h + 2 notin E(X)` whenever `H(X) >= 5`. (Lemma used: a vertex that lies on an odd cycle of
   length `>= 5` lies on a shorter odd cycle, hence on a triangle; a vertex on exactly one odd cycle is on a triangle.)
3. **Increment four (PROVED structure; FINITE-EXACT value).** `H(T) - H(T - v) = 4` iff `v` lies on exactly two odd
   cycles, each with transitive complement `T - V(C)`. Through order `9` this forces `H(T - v) in {1, 5}` (and
   increment `6` forces `H(T - v) in {3, 5, 9}`, increment `8` forces `H(T - v) in {1, 5, 9, 11}`).
4. **Flip increment (PROVED) and THM-013 for all orders.** Let `i -> j` be an arc of `T`, `T'` the tournament with that
   arc reversed, `W = V - {i, j}`, and for `U >= W`, `L_x(U)`, `R_x(U)` the numbers of Hamiltonian paths of `T[U]` ending,
   resp. starting, at `x`. Then

       H(T') - H(T) = U(L_j, R_i) - U(L_i, R_j),   U(L_x, R_y) = sum_{S <= W} L_x(S cup {x}) R_y({y} cup (W - S)),

   the two unsigned subset convolutions counting the Hamiltonian paths of `T'` in which `j` immediately precedes `i`
   and those of `T` in which `i` immediately precedes `j`; the signed convolutions `B(L_i, R_j)`, `B(L_j, R_i)`
   (alternating in `|S|`) coincide (THM-017). Since `H = I(Omega, 2)` for every tournament, `H(T) - H(T') =
   I(Omega(T), 2) - I(Omega(T'), 2) = sum_k 2^k Delta(alpha_k)` for every `n`: THM-013's "general formula", proved
   there for `n <= 7` and verified to `n = 9`, holds for all `n`.
5. **Extension spectra (FINITE-EXACT, orders 3-9).** Over all tournaments `X` of order `n <= 8` (all cuts), the
   smallest positive increment `k_min(h) = min {(H(X + v) - h)/2 > 0 : H(X) = h}` is

   | h | 1 | 3 | 5 | 9 | 11 | 13 | 15 | 17 | 19 | 23 | 25 | 27 | 29 | 31 | 33 | 35 | 37 | 39 |
   |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
   | k_min | 1 | 1 | 2 | 3 | 4 | 5 | 5 | 5 | 6 | 7 | 9 | 9 | 8 | 10 | 9 | 11 | 12 | 11 |

   (the same table at orders 7 and 8; order 9 extends it, `increment2_order9.out`), i.e. `k_min(h)/h` stays in
   `[0.27, 0.40]`. Every increment-`2` cut has `h in {1, 3}` (statement 2), every increment-`4` cut `h in {1, 5}`.
   The longest run of consecutive odd values inside `E(X)` has length at most `24` at order `8` and is never an
   interval containing `[2h + 1, 3h]`; among the `6880` tournaments of order `8` none has `E(X) >= odd[h, 2h - 1]`.
   Hence the two natural single-realiser inductions for Chan-Pak (reach every odd number in `[h, 2h - 1]`, or a
   window `[ch, c'h]`, from one realiser of `h`) are impossible.
6. **One-centre family (PROVED).** Let `X = TT_m + u` with `u -> A`, `[m] - A -> u`. Every odd cycle passes through `u`,
   so `H(X) = 1 + 2 N(A)` with `N(A) = sum_{a in A, b notin A, a < b} c(b - a)`, `c(1) = 1`, `c(d) = 2^(d-2)`. In
   particular `A = {1}` gives `2^(m-1) + 1` (the reversed-extreme-arc family of THM-4749) and `A = {1, ..., i}` with
   `m = i + l` gives `(2^i - 1)(2^l - 1) + 2`; `N(A) = 3` is impossible (so `7` is not a one-centre count, in line
   with THM-1370), and the family realises about half of the odd numbers (`50,006` of the `100,000` below `200,000`;
   missing `7, 13, 21, 29, 35, 41, 43, 49, ...`). A second centre whose cut is monotone in the transitive order is
   again one-centre (its deletion of `u` is transitive), so it adds nothing.

## Proofs

**1.** Claim A: in `I(Omega(T), 2)` split the families `F` by whether they contain a cycle through `v`; those that do
not are the families of `T - v`, and those that do contain exactly one such cycle `C` (two odd cycles through `v`
intersect), the rest being a family of `T - V(C)`; so `I(Omega(T), 2) = I(Omega(T - v), 2) + 2 sum_C I(Omega(T - V(C)), 2)`,
and the odd-cycle formula converts each term to a path count. Path pairs: deleting `v` from a Hamiltonian path of
`X + v` leaves an ordered pair `(P_1, P_2)` of vertex-disjoint paths of `X` covering `V(X)`, `P_1` possibly empty and
otherwise ending in `B`, `P_2` possibly empty and otherwise starting in `A`; conversely every such pair gives one
Hamiltonian path of `X + v`. Summing over `S = V(P_1)` gives the formula.

**2.** *Lemma.* If `v` lies on an odd cycle `C = v -> a_1 -> ... -> a_k -> v` with `k >= 4` (`k` even), then `v` lies on
a shorter odd cycle: if `v -> a_2` fails then `a_2 -> v` and `(v, a_1, a_2)` is a triangle; if `a_(k-1) -> v` fails then
`v -> a_(k-1)` and `(v, a_(k-1), a_k)` is a triangle; otherwise `v -> a_2 -> ... -> a_(k-1) -> v` is a cycle of length
`k - 1`, odd. So a vertex on exactly one odd cycle is on a triangle. Now let `H(T) - H(T - v) = 2`. By Claim A the sum
`sum_C H(T - V(C))` equals `1` with positive odd terms, so `v` lies on exactly one odd cycle `C`, which is a triangle
`(v, a, b)` by the lemma, and `T - V(C)` is transitive. Write `A`, `B` for the out- and in-neighbourhoods of `v`, so
`a in A`, `b in B`, and let `R = V - {v, a, b}` (transitive), `A' = A - {a}`, `B' = B - {b}`. A 2-vertex path `x -> y`
with `x in A`, `y in B` closes a triangle `(v, x, y)`, which must be `C`: hence the only arc from `A` to `B` is `a -> b`.
Every path from `A` to `B` therefore passes through the consecutive arc `a -> b`, and the cycle through `v` it closes is
odd iff the path has an even number of vertices, i.e. iff its `A`-part (a path inside `A` ending at `a`, with `j`
vertices) and `B`-part (inside `B` from `b`, `j'` vertices) satisfy `j + j'` even. Uniqueness of `C` forbids `(j, j') =
(2, 2)`: not both an in-neighbour `x in A'` of `a` and an out-neighbour `y in B'` of `b` exist. By converse symmetry
(`H` is invariant under reversing all arcs, which swaps the roles of `A` and `B`) assume `a -> A'`. Then every `A`-part
is `(a)` alone (`j = 1`), so uniqueness forbids every odd `j' >= 3`: an out-neighbour `y` of `b` in `B'` has no
out-neighbour in `B'`, i.e. is the sink of the transitive `B'`; so `b` beats at most one vertex `s` of `B'`, and
`B' - s -> b`. Now `T - v = (A' cup {a}) cup (B' cup {b})` with all arcs from the `B`-side to the `A`-side except
`a -> b`; the `A`-side is transitive (`a -> A'`, `A'` transitive), the `B`-side has no cycle (a cycle there would contain `b`, since `B'` is
transitive, leave `b` to its only out-neighbour `s` in `B'`, and return to `b` from a vertex of `B'` reachable from `s`,
but the sink `s` has no out-neighbour in `B'`), and a cycle meeting both sides uses `a -> b`,
returns to the `A`-side by some `y -> x`, and must then reach `a` inside the `A`-side, where `a` has no in-neighbour:
so `x = a`, the cycle is `a -> b -> ... -> y -> a` with the `B`-part a path from `b`, i.e. `(b)` (impossible, `b -> a`
fails) or `(b, s)`: the only possible cycle of `T - v` is the triangle `(a, b, s)`. Hence `T - v` has at most one odd
cycle, `H(T - v) = 1 + 2 [b -> s] in {1, 3}`, and `H(T) = H(T - v) + 2 in {3, 5}`. Conversely a triangle through `v`
with transitive complement and no other odd cycle through `v` gives increment `2` by Claim A (realised by `C_3`, `v`
any vertex: `3 - 1`; and by the strong 4-tournament, `v` its vertex of out-degree 1 or in-degree 1: `5 - 3`).

**3.** Claim A with `sum_C H(T - V(C)) = 2` and odd terms: exactly two cycles, both terms `1`. The value lists are
`increment_ladder.py` (orders of `T` up to `9`).

**4.** The Hamiltonian paths of `T` in which `i` and `j` are not consecutive are exactly those of `T'` with the same
property (they use no arc between `i` and `j`). A Hamiltonian path of `T` with `i` immediately followed by `j`
decomposes uniquely as a Hamiltonian path of `T[S cup i]` ending at `i` followed by one of `T[j cup (W - S)]` starting
at `j`, for a unique `S <= W`; the two pieces avoid the pair `{i, j}`, so the same count in `T'` is `U(L_i, R_j)`
computed in `T` as well. Likewise the paths of `T'` with `j` immediately followed by `i` number `U(L_j, R_i)`. Subtract.
The signed statement is THM-017; `Delta H = Delta I(Omega, 2)` is the odd-cycle formula applied to `T` and `T'`.

**5.** Computation (`ext_spectrum_true.py`: the path-pair formula vectorised over all cuts, brute-force checked;
`ext_windows.py`; `increment_ladder.py`, which also enumerates the odd cycles through `v` and checks Claim A on every
cut of increment `<= 8`; `increment2_order9.py` for order `9`).

**6.** `T - u = TT_m` is acyclic, so every odd cycle of `X` passes through `u` and any two intersect: the independent
families of `Omega(X)` are the empty family and the singletons, `H(X) = 1 + 2 #(odd cycles)`. An odd cycle through `u`
is `u -> a -> (increasing path) -> b -> u` with `a in A`, `b notin A`, `a < b`, and the inner path is `a`, a subset `S` of
the `b - a - 1` vertices strictly between, then `b`; the cycle has `|S| + 3` arcs, odd iff `|S|` is even, which gives
`2^(b-a-2)` choices for `b - a >= 2` and `1` for `b - a = 1`. For `A = {1, ..., i}`, `m = i + l`: with
`sum_{d=1}^D c(d) = 2^(D-1)` one gets `N = 2^(l-1) + (2^(i-1) - 1)(2^l - 1)` and `1 + 2N = (2^i - 1)(2^l - 1) + 2`.
`N(A) = 3` would be three adjacent `10`-pairs (then two of them at positions `p < p'` give the extra pair `(p, p' + 1)`)
or a pair at distance `3` plus an adjacent pair (then the letters between create a further pair); both exceed `3`.
The coverage counts are `spine_numeration.py` (the same family as the session124 spine with the transitive parent,
whose transfer `(F, q, t)` is re-derived and brute-force checked there). ∎

## Scope and boundary

* Statement 2 is a theorem about *all* tournaments; statement 3's "`H(T - v) in {1, 5}`" and the whole of statement 5
  are census facts (orders of `T` up to `9`, resp. `10` for `k_min` at order `9` when that run finishes), typed
  FINITE-EXACT. The linear growth of `k_min(h)` is HYP-9469.
* Statement 5 does not touch the repository's actual coverage machinery (the digital bank to `10^8`, the sole-exit
  tower, the near-join), which extends *different* realisers at each step; it only shows that the extension spectrum
  of one realiser is too sparse near `h` to carry an induction.
* The one-centre count in 6 is the session124 "transitive parent" spine family in closed form; its value set is not
  characterised here (HYP-9469 (3)).

## Reading

Claim A says that adding a vertex raises `H` by twice a sum of path counts of complements of odd cycles through it.
The smallest conceivable raise, `2`, needs a single triangle with a transitive complement, which traps the whole
tournament at five paths or fewer; the census shows the smallest raise above a count `h` is about `h/3.5`, not `2`.
So the Hamiltonian-path spectrum is not "locally connected" at any realiser: Chan-Pak completeness, if true, is a
global statement about the union of many families, as the repository's banks already treat it. The March-2026
theorems THM-013 and THM-017 are the unsigned and signed faces of one arc-flip decomposition, and the odd-cycle formula
of Grinberg-Stanley closes the gap THM-013 left open.
