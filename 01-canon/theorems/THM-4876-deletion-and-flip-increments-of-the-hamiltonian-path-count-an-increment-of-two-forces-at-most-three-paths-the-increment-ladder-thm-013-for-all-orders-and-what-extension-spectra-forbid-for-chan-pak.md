---
id: THM-4876
title: "Deletion and flip increments of the Hamiltonian-path count H. (1) For a tournament T and a vertex v, H(T) - H(T - v) = 2 sum over the odd directed cycles C through v of H(T - V(C)) (Claim A, KNOWN from the odd-cycle formula), and for a tournament X with a new vertex v of cut (A, B), H(X + v) = sum_S e_B(S) s_A(V - S) + s_A(V) + e_B(V) (path pairs). (2) H(T) - H(T - v) = 2 iff v lies on exactly one odd cycle, a triangle (v, a, b), and T - {v, a, b} is transitive; then T - v has at most one odd cycle, H(T - v) in {1, 3} and H(T) in {3, 5}: no deletion lowers a count of at least 5 by exactly 2. (3) H(T) - H(T - v) = 4 iff v lies on exactly two odd cycles, both with transitive complements (and H(T - v) in {1, 5} through order 10). (4) (KNOWN by repository since 2026-03-05: PROP-001 and THM-013/014/017, MISTAKE-008; restated here) The arc-flip increment is the unsigned subset convolution H(T') - H(T) = U(L_j, R_i) - U(L_i, R_j) whose signed companion is THM-017's identity, and THM-013's Delta H = Delta I(Omega, 2) holds for all n as a corollary of the odd-cycle formula H = I(Omega, 2). (5) FINITE-EXACT: the one-vertex extension spectra E(X) of all tournaments of orders 3-9 have no value h + 2 for h >= 5, smallest positive increments k_min(h) = 2, 3, 4, 5, 5, 5, 6, 7, 9, 9, 8, 10, 9, 11, 12, 11 for h = 5, 9, 11, ..., 39 (between h/3.7 and h/2.5; audit C: k_min(h) >= h/4 for every h attained at order 10), and windows of at most 24 consecutive odd values which, for h notin {1, 5}, never cover odd[2h + 1, 3h] at orders <= 8: the induction through [h, 2h - 1] from one realiser is impossible at every order by (2), and the window induction is unavailable through order 8. (6) The one-centre family TT_m + u has H = 1 + 2 sum_{a in A, b notin A, a < b} c(b - a), c(1) = 1, c(d) = 2^(d-2), hence contains 2^k + 1 and (2^i - 1)(2^l - 1) + 2 for all i, l >= 1 and omits 7, 13, 21, 29, ... (7) A vertex that is not a sink ends at most (h - 1)/2 of the h Hamiltonian paths (dually, a non-source starts at most (h - 1)/2), with equality at the penultimate vertex of TT_(k+2) with its extreme arc reversed (k >= 1, h = 2^k + 1); so for h >= 5 no vertex ends h - 1 or h - 2 paths, and the sole-exit extension at a non-sink raises the count by at least (h + 1)/2 (through THM-4862 (a)'s law 2h - e_a + sigma), i.e. k >= (h + 1)/4 in the normalisation of (5)"
status: >
  PROVED (mac-mini-2026-11-04-chanpak, elementary): statements 2, 3 (the structural halves), 6, 7, and the
  path-pair formula in 1 (statement 7 was added after audit C was launched and is covered by audit D). KNOWN by
  repository (audit C, correction C1): statement 4 - the unsigned flip decomposition is THM-013 "Significance",
  THM-014 and THM-017/MISTAKE-008 (2026-03-05), and Delta H = Delta I(Omega, 2) for all n is PROP-001 (with
  THM-002/CONJ-001, kind-pasteur S12, 2026-03-05); it is restated here only to place the two faces side by side.
  KNOWN/CITED: Claim A and the odd-cycle formula H(T) = I(Omega(T), 2) (Grinberg-Stanley, arXiv:2307.05569,
  Theorem 1.39, restated as Irving-Omar's Corollary 20 per THM-002; the repository's THM-070 derives Claim A from it),
  Redei's theorem (H odd), Moon (every vertex of a strong tournament lies on a 3-cycle), THM-017 (signed subset
  convolution, PROVED for all n in March 2026), THM-013 (the same identity, PROVED for all n as a corollary of OCF
  since PROP-001), THM-012b (insertion decomposition), THM-4873 (end-count lemma, opus-S77), THM-4862 (a) (one-exit
  extension law, used in 7). FINITE-EXACT: statement 5 and the value lists in 3 and 6
  (chanpak_spine_20261104/ext_spectrum_true.py, ext_windows.py, increment_ladder.py, increment2_order9.py,
  spine_numeration.py; formulas brute-force checked: path pairs on 468 (X, cut) pairs of orders <= 5, Claim A at v on
  every cut with increment <= 8 of orders <= 8, the session124 spine transfer on 69,174 triples). Audited: audit C (independent, adversarial,
  05-knowledge/results/chanpak_increments_audit_C_20261104.md, its own C census through order 10) PASS WITH
  CORRECTIONS, all fourteen applied on 2026-11-04 (MISTAKE-715): statement 4 retyped KNOWN, two sentences of
  statement 5 repaired (the h in {1, 5} and h in {1, 3} exceptions), the Lass attribution dropped, the N(A) = 3
  argument completed, the strong 4-tournament realiser corrected; audit C also verified statements 2 and 3
  exhaustively through order 10 and the k_min table against its order-10 census. Audit D (statement 7, independent, 05-knowledge/results/chanpak_endcount_audit_D_20261104.md, own C code): PASS WITH
  CORRECTIONS, seven minor items applied (D1-D7: the h = 3 wording, the empty-prefix convention, k >= 1, committed
  receipts for the order-9 check, the scripts field, the dual's wording, the sharper 2 e_a + e_O <= h); audit D verified
  the bound, its dual, 2 e_a + e_O <= h, the sole-exit law and the floor (3h + 1)/2 on every class of orders <= 9. The Chan-Pak conjecture itself (every odd number other than 7 and 21 is a Hamiltonian-path count) remains
  OPEN; this theorem closes one route (statement 5) and sharpens the local structure of the spectrum.
session: mac-mini-2026-11-04-chanpak
source: 05-knowledge/results/chanpak_increments_thm013_thm017_merge_20261104.md
scripts:
  - 04-computation/experiments/chanpak_spine_20261104/ext_spectrum_true.py, ext_windows.py, increment_ladder.py (+ _n8.out), increment2_order9.py (+ .out), spine_numeration.py, two_centre_family.py, endcount_ratio.py (+ _n8.out), verify_stmt7_consequences.py (+ .out), ordered_paths_lemma.py (+ .out)
related:
  - 01-canon/theorems/THM-013-arc-flip-delta-I.md (cross-referenced by statement 4; the all-orders identity is PROP-001-arc-flip-identity.md), THM-017-even-odd-split-proof.md, THM-016-claim-b-alternating-sum.md, THM-012b-insertion-decomposition.md, THM-070-claim-a-from-ocf.md, THM-4862 (one-exit law), THM-4873 (end-count lemma)
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
`e_B(S) = sum_{u in S cap B} e_S(u)`, `s_A(S) = sum_{u in S cap A} s_S(u)`. The *odd-cycle formula* (Grinberg-Stanley,
arXiv:2307.05569, Theorem 1.39) is `H(T) = I(Omega(T), 2) = sum over families F of pairwise vertex-disjoint odd directed cycles of
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
   cycles, each with transitive complement `T - V(C)`. Through order `10` (order `9` by the session, `10` by audit C's
   census) this forces `H(T - v) in {1, 5}` (and
   increment `6` forces `H(T - v) in {3, 5, 9}`, increment `8` forces `H(T - v) in {1, 5, 9, 11}`).
4. **Flip increment and THM-013 for all orders (KNOWN by repository: THM-013/014/017 and MISTAKE-008 for the unsigned
   decomposition, PROP-001 with THM-002/CONJ-001 for the all-orders identity, 2026-03-05; restated, audit C C1).** Let `i -> j` be an arc of `T`, `T'` the tournament with that
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

   (the same table at orders 7, 8 and 9 of `X`; the order-9 run, `increment2_order9.out`, all 191,536 tournaments
   with all 512 cuts, adds `k_min = 14, 15, 14, 13, 13, 15, 14, 18, 18, 17` at `h = 41, 43, 45, 47, 49, 51, 53, 55, 57,
   59`), i.e. `k_min(h)/h` stays in `[0.26, 0.40]`; audit C's independent order-10 deletion census (`audit_C/del10_merged.out`)
   reproduces the table and gives `k_min(h) >= h/4` for every `h` attained at order `10` (`h <= 3357`), least ratio
   `33/129`, the minima at `h = 2^k + 1` being the sole-exit cuts at the penultimate vertex of `R_(k+2)` with digit `1`
   (so `k_min(2^k + 1) <= 2^(k-2) + 1` there, orders `4-10`). Every increment-`2` cut has `h in {1, 3}` (statement 2;
   confirmed through order `9` of `X`, and by audit C through order `10`), every increment-`4` cut `h in {1, 5}`
   (audit C: through order `10`).
   The longest run of consecutive odd values inside `E(X)` has length at most `24` at order `8`, and for
   `h notin {1, 5}` no tournament of order `<= 8` has `odd[2h + 1, 3h] <= E(X)` (direct containment test,
   `audit_C/auditC_window_direct.out`; the realisers of `h = 5`, the strong 4-tournament with `E = {5, 9, 11, 13, 15}`
   and its paddings, do contain `{11, 13, 15}`, and `h = 1` does trivially); among the `6880` tournaments of order `8`
   none with `h >= 5` has `E(X) >= odd[h, 2h - 1]` (the seven with `h in {1, 3}` do), which is statement 2's
   `h + 2 notin E(X)` at every order. Hence of the two natural single-realiser inductions for Chan-Pak, reaching
   every odd number in `[h, 2h - 1]` from one realiser of `h` is impossible at every order (PROVED, statement 2), and
   reaching the window `[2h + 1, 3h]` (so any `[ch, c'h]` with `c <= 2 < 3 <= c'`) from one realiser is unavailable
   through order `8` except at `h in {1, 5}` (FINITE-EXACT).
6. **One-centre family (PROVED).** Let `X = TT_m + u` with `u -> A`, `[m] - A -> u`. Every odd cycle passes through `u`,
   so `H(X) = 1 + 2 N(A)` with `N(A) = sum_{a in A, b notin A, a < b} c(b - a)`, `c(1) = 1`, `c(d) = 2^(d-2)`. In
   particular `A = {1}` gives `2^(m-1) + 1` (the reversed-extreme-arc family of THM-4749) and `A = {1, ..., i}` with
   `m = i + l` gives `(2^i - 1)(2^l - 1) + 2`; `N(A) = 3` is impossible (so `7` is not a one-centre count, in line
   with THM-1370), and the family realises about half of the odd numbers (`50,006` of the `100,000` below `200,000`;
   missing `7, 13, 21, 29, 35, 41, 43, 49, ...`). A second centre whose cut is monotone in the transitive order is
   again one-centre (its deletion of `u` is transitive), so it adds nothing.
7. **End-count bound (PROVED; added after audit C).** Let `a` be a vertex of a tournament `X` with `h = H(X)` and let
   `e_a` be the number of Hamiltonian paths ending at `a`. If `a` is not a sink, then `e_a <= (h - 1)/2`; dually, if `a`
   is not a source, the number `s_a` of paths starting at `a` is at most `(h - 1)/2`. Equality holds at the penultimate
   vertex `k + 1` of `TT_(k+2)` with the arc `(1, k + 2)` reversed (`k >= 1`, `h = 2^k + 1`, `e_(k+1) = 2^(k-1)`; its
   single out-neighbour is the former sink `k + 2`). More precisely `2 e_a + e_O <= h`, where `e_O` counts the paths
   ending in `N+(a)` (audit D, D7). Consequences: for `h >= 5` no vertex ends `h - 1` or `h - 2` paths (the case `h = 3`
   of THM-4873 (i), `e_a = h - 2 = 1`, occurs only for `TT_j => C3` with `a` in the `C3`), and the
   sole-exit extension `v -> a` (`a` a non-sink; `v -> a` and every other vertex `-> v`) has `H(X + v) = 2h - e_a + sigma`
   with `sigma = [P A_a](X) - h >= 0` (THM-4862 (a); `sigma >= 0` because splitting a Hamiltonian path of `X` at `a` is an
   injection into the pairs counted by `[P A_a]`, with `h(X[empty]) = 1` - the convention of THM-4862's source packet -
   so that the paths starting at `a` land in the `S = empty` term), hence `H(X + v) >= (3h + 1)/2`: the count rises by at least
   `(h + 1)/2`, i.e. `k >= (h + 1)/4` in the normalisation of statement 5 (half the raw difference): the sole-exit half
   of HYP-9469 (1); the general cut reduces to the ordered-path lemma HYP-9469 (6), of which this is the case `|A| = 1`.
   (Checked exhaustively on all tournaments of orders `<= 9` by audits C and D with their own code,
   `audit_C/del9.out` and `audit_D/classes_n9.out`: no violation of the bound, of its dual or of `2 e_a + e_O <= h`; the
   maximum of `e_a/h` over non-sinks at order `9` is `64/129`, and `128/257` at order `10` by audit C; the session's own
   receipt, `endcount_ratio_n8.out`, stops at order `8`. Audit D's census observations, FINITE-EXACT through order `9`:
   the floor `(3h + 1)/2` is attained only at `h = 3` (`TT_j => C3`); for `h >= 5` the least sole-exit count is
   `(3h + 3)/2`, on `R_(k+2)` where `sigma = 1`; and the equality cases of the bound with `h >= 3` are exactly
   `TT_j => R_(k+2)` at the vertex `k + 1`, one class per `(j, k)`.)

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
any vertex: `3 - 1`; and by the strong 4-tournament, `v` either of its two vertices that lie on only one of its two triangles - one of
out-degree `1`, one of in-degree `1`; the other two vertices lie on both triangles and give `5 - 1 = 4`: `5 - 3`).

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
`N(A) = 3` is impossible: let `a` be the first `1` and `b` the last `0` with `a < b` (if there is none, `N = 0`); every
contributing pair lies inside `[a, b]`. If `b - a >= 4` the pair `(a, b)` alone contributes `c(b - a) >= 4`. If
`b - a = 3`, `(a, b)` contributes `2` and each filling of the two middle letters adds at least two further pairs
(`1000`: `(a, a+1), (a, a+2)`; `1100`: `(a, a+2), (a+1, a+2), (a+1, b)`; `1010`: `(a, a+1), (a+2, b)`; `1110`:
`(a+1, b), (a+2, b)`), so `N >= 4`. If `b - a <= 2` the window is `10`, `100` or `110`, with `N = 1, 2, 2`. So
`N in {0, 1, 2} cup [4, infinity)` (attained values below `12`: `0, 1, 2, 4, 5, 7, 8, 9, 11`; audit C C7).
The coverage counts are `spine_numeration.py` (the same family as the session124 spine with the transitive parent,
whose transfer `(F, q, t)` is re-derived and brute-force checked there). ∎

**7.** Let `O` be the out-neighbourhood of `a` (nonempty) and `E_a` the set of Hamiltonian paths ending at `a`. For
`P = (x_1, ..., x_(n-1), a) in E_a` let `i` be the least index with `x_i in O` (it exists since `O != empty` and `a` is last).
Then `x_(i-1) -> a` if `i >= 2` (as `x_(i-1) notin O`), so `phi(P) = (x_1, ..., x_(i-1), a, x_i, ..., x_(n-1))` is a Hamiltonian
path; it does not end at `a`, and its last vertex `x_(n-1)` beats `a`. `phi` is injective, since `P` is recovered by
deleting `a` and appending it. Hence `e_a <= |F|`, where `F` is the set of Hamiltonian paths not ending at `a` whose last
vertex beats `a`; `|F| <= (h - e_a) - e_O` with `e_O` the number of paths ending in `O`. If `e_a = 0` there is nothing to
prove. If `e_a >= 1`, `a` lies in the last strong component `X_s` of `X`, which has at least three vertices (`a` is not
a sink), so `a` has an out-neighbour `z` in `X_s`; every vertex of a strong tournament is the end of a Hamiltonian path
of it (break a Hamiltonian cycle, Camion), and prefixing Hamiltonian paths of the earlier components gives a Hamiltonian
path of `X` ending at `z`: `e_O >= 1`. So `2 e_a <= h - 1`. (Rédei's parity alone gives `2 e_a <= h - 1` from `2 e_a <= h - e_O <= h`; the paragraph
proves the sharper `2 e_a + e_O <= h`, tight on the family below, where `e_O = 1`.) The dual statement is the converse
tournament. Equality: in
`TT_(k+2)` with the arc `(1, k + 2)` reversed, every Hamiltonian path other than the transitive one uses the only
backward arc `k + 2 -> 1`, so it reads (`S` increasing), `k + 2`, `1`, (complement of `S` increasing) for a subset
`S` of `{2, ..., k + 1}`, and every `S` works: `h = 2^k + 1`. The vertex `k + 1` (single out-neighbour `k + 2`) is last
exactly when `k + 1 notin S`, i.e. for `2^(k-1) = (h - 1)/2` paths. (Brute force, orders `3..8`: `e_(m-1) = 2^(m-3)` with
`h = 2^(m-2) + 1`, and `m - 1` is the only non-sink vertex attaining the bound for `m >= 4`.) ∎

## Scope and boundary

* Statement 2 is a theorem about *all* tournaments; statement 3's "`H(T - v) in {1, 5}`" and the whole of statement 5
  are census facts (orders of `T` up to `9`, resp. `10` for `k_min` at order `9` when that run finishes), typed
  FINITE-EXACT. The linear growth of `k_min(h)` is HYP-9469.
* Statement 5 does not touch the repository's actual coverage machinery (the digital bank to `10^8`, the sole-exit
  tower, the near-join), which extends *different* realisers at each step; it only shows that the extension spectrum
  of one realiser is too sparse near `h` to carry an induction.
* The one-centre count in 6 is the session124 "transitive parent" spine family in closed form; its value set is not
  characterised here (HYP-9469 (3)).
* Statement 7 bounds the sole-exit increment only; the linear bound for arbitrary cuts (HYP-9469 (1)) is open.

## Reading

Claim A says that adding a vertex raises `H` by twice a sum of path counts of complements of odd cycles through it.
The smallest conceivable raise, `2`, needs a single triangle with a transitive complement, which traps the whole
tournament at five paths or fewer; the census shows the smallest increment `k_min(h)` of statement 5 lies between `h/3.7` and `h/2.5` for
`5 <= h <= 59` (audit C: `k_min(h) >= h/4` for every `h` attained at order `10`), far above the `1` of a single triangle.
So the Hamiltonian-path spectrum is not "locally connected" at any realiser: Chan-Pak completeness, if true, is a
global statement about the union of many families, as the repository's banks already treat it. The March-2026
theorems THM-013 and THM-017 are the unsigned and signed faces of one arc-flip decomposition, and the odd-cycle formula
of Grinberg-Stanley closes the gap THM-013 left open.
