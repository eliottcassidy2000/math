---
id: THM-4812
title: "Presentation fibres 7 and 21 are absent at every order: no tournament has tc = H/|Aut| in {7, 21} (the tc version of THM-1370; HYP-9352 (2) PROVED). (a) Order 12 by exhaustive invariant scans over the order-12 type lemma (THM-4790 had n <= 11); (b) orbit product lemma: tc(T) >= prod over the Aut(T)-orbits O of tc(T[O]), each T[O] vertex-transitive; (c) vertex-transitive lemma: every vertex-transitive tournament other than K1, C3, R5 and P7 (tc = 1, 1, 3, 9) has tc >= 23; (d) so a counterexample would contain a strong S with tc(S) = 21, 13 <= |S| <= 20 or |S| = 22, all of whose Aut-orbits induce K1, C3, R5 or P7, and Aut(S) a subgroup of index 1 or 3 of the product of the orbit groups; an exhaustive search of the 83 resulting (configuration, subgroup) cases finds none; (e) the prime untouchables have relative density zero among the primes; (f) every strong presentation fibre is finite: a strong T with tc(T) = v has order at most 3 ln(v^(1 + M(v)/6)/0.5848)/ln(5/3), M(v) the largest order of a vertex-transitive tournament with tc <= v, so O(log^2 v); (g) the strong fibres for odd v <= 21 have sizes 2, 1, 2, 0, 5, 4, 1, 6, 6, 4, 0 (v = 1, 3, ..., 21), all of order at most 8."
status: >
  PROVED (all orders), given the cited theorems and the finite certificates below: tc(T) is never 7 or 21.
  PROVED: (b) from THM-002 (OCF: H(T) >= prod H(T[B_i]) over any vertex partition) and the faithful restriction of Aut(T) to
  its orbits; (c) for prime orders (a transitive group of prime degree contains a p-cycle, so the tournament is a circulant),
  for order 15 (an odd-order transitive group of degree 15 is solvable, hence imprimitive, and either contains a regular Z_15
  or contains one of two explicit groups G_A = A_3 x| Z_5, G_B = A_5 x| Z_3), and for orders m >= 19 (greedy packing of
  disjoint directed triangles in the OCF sum against Dixon's bound |Aut| <= 3^((m-1)/2)); (d) the reduction from (b), (c), THM-4790 (a),
  the strong minimum f(n) (Busch 2006, Corollary 1, as in THM-4744) and the normality of index-3 subgroups of odd-order groups; (e) from Montgomery-Vaughan's exceptional-set theorem for Goldbach and the sieve bound for primes p = 2q + 1.
  FINITE-EXACT: (a) by exhaustive scans of all tournaments invariant under an element of each type of the order-12 type
  lemma (proved below): (11) none strong; (7) 178,560 strong; (5)(5) 9,024; (5)(3) 178,560; (5)(3)(3) 2,112;
  (3)(3)(3)(3) 4,126,848; (3)(3)(3) 13,643,136; 3,491,072 automorphism-group computations on the candidates H = 7a or 21a,
  no hit; the same program reproduces THM-4790's counts 1,088 (n = 7, (3)), 3,760 (n = 9, (3)^3), 23,520 (n = 10, (3)^3),
  480 (n = 10, (5)(5)) and 392,000 (n = 11, (3)^3). (c) at orders 3-17: Cayley tournaments on Z_3, Z_5, Z_7, Z_9, Z_3^2,
  Z_11, Z_13, Z_15, Z_17 (minimum tc 1, 3, 9, 39, 39, 1729, 94949, 162165, 792878221; tc <= 21 only for C3, R5, P7); the
  complete census of the 15 regular 9-tournaments (vertex-transitive ones: tc 363, 373, 39); every G_A- and G_B-invariant
  tournament (tc 162165 and 515241); the greedy bound for every odd m in [19, 401], with an analytic tail beyond.
  CITED: THM-4790 (reduction to strong components, n <= 11), THM-1370 (H never 7 or 21), THM-002 (OCF), THM-4741 (tc = 1 iff
  K1/C3 stack), Dixon 1967 (|Aut T| <= 3^((n-1)/2)), Moon 1972 / Busch 2006 / THM-4744 (strong minimum f); Feit-Thompson 1963 and Galois's theorem on solvable primitive groups (or the list of primitive groups of degree 15),
  Montgomery-Vaughan 1975 (Acta Arith. 27) and the Brun-Selberg upper-bound sieve for (e). The {2, 5} prime-untouchable completeness remains conditional on the
  distinct-Goldbach route of THM-4805.
  FINITE-EXACT (d): window_configs.py + tc_window_scan.cpp: 24 orbit configurations survive the counting, all with tc = 21;
  index-1 and index-3 subgroups with transitive projections give 83 cases; all 185,642 strong invariant tournaments with the
  prescribed orbit structure were scanned (H by the block-run substitution DP for the 73 module cases, cross-checked against
  Held-Karp on all 93,646 module tournaments of order <= 16, and by Held-Karp for the 10 non-module cases), no tc in {7, 21}.
  PROVED (f) from (b), (c)'s greedy bound, THM-4741 and Busch's f(n).
  FINITE-EXACT (g): the strong census of orders 3-10 (tc_strong_fibres.cpp) and the generalised window search
  (fibre_configs.py, fibre_scan.cpp: 567 configuration-subgroup cases, 811,496 strong invariant tournaments of orders >= 11,
  no tc value <= 21) with the order bounds of the same counting.
  OPEN: whether every odd v not in {7, 21} is a tc value (HYP-9352 (1)-type spectrum question); an audit of the window search.
session: mac-mini-2026-10-27-potentials
source: 05-knowledge/results/dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md
depends_on:
  - 01-canon/theorems/THM-4790-presentation-fibres-7-and-21-absent-through-eleven-vertices.md
  - 01-canon/theorems/THM-002-ocf.md
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
related:
  - 01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md
  - 05-knowledge/hypotheses/HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md
scripts:
  - 04-computation/experiments/tc_fibres_n12_20261027/tc_invariant_scan.cpp (+ tc_invariant_scan_n12.out)
  - 04-computation/experiments/tc_fibres_n12_20261027/cayley_tc.cpp (+ cayley_tc.out)
  - 04-computation/experiments/tc_fibres_n12_20261027/regular_vertex_transitive.py (+ regular9_vertex_transitive.out)
  - 04-computation/experiments/tc_fibres_n12_20261027/vt15_noncayley.py (+ vt15_noncayley.out)
  - 04-computation/experiments/tc_fibres_n12_20261027/type_lemma_checks.py (+ .out)
  - 04-computation/experiments/tc_fibres_n12_20261027/window_configs.py, tc_window_scan.cpp (+ window_scan.out)
  - 04-computation/experiments/tc_fibres_n12_20261027/tc_strong_fibres.cpp, fibre_configs.py, fibre_scan.cpp, strong_fibre_table.py (+ fibre_scan_v5_21.out, strong_fibre_table.out)
  - 04-computation/experiments/tc_fibres_n12_20261027/tc_spectrum_sweep.cpp (+ tc_spectrum_summary.out)
---

# THM-4812 — the {7, 21} presentation gaps hold at every order

`tc(T) = H(T)/|Aut(T)|` is the number of marked-path presentations of `T` (THM-4741); it is an odd integer (LEM-003).
THM-1370 excludes `H in {7, 21}` at every order; THM-4790 excludes `tc in {7, 21}` through order 11 and reduces any
counterexample to a strong one. The two exclusions had different scopes (all orders for H, at most 11 vertices for tc);
this theorem gives the tc exclusion the same scope as THM-1370.

## (a) Order twelve

**Type lemma (n = 12).** An odd-order `G <= S_12` with `|G| >= 19` contains an element of type `(11)`, `(7)`, `(5)(5)`,
`(5)(3)`, `(5)(3)(3)`, `(3)(3)(3)` or `(3)(3)(3)(3)` (fixed points omitted).

*Proof.* Primes 11 and 7 give 11- and 7-cycles (Cauchy). If `25 | |G|` a Sylow 5-subgroup is `Z_5^2` (no 25-cycle on 12
points); two non-proportional single 5-cycles in it commute, so their supports are disjoint (a 5-cycle commuting with a
5-cycle on an overlapping support preserves that support and acts there as a power of it), and the product has type
`(5)(5)`. If `G` is a 3-group, `|G| >= 27`. Its orbits have sizes 1, 3 or 9. An orbit of size 9 carries a derangement of
3-power order (Jordan); a suitable power of order 3 acts on that orbit as `(3)^3`, giving type `(3)^3` or `(3)^4`. With
only 3-orbits `O_1..O_r`, `r <= 4`, `G` embeds in `F_3^r` with dimension at least 3, and every subspace of dimension at
least 3 of `F_3^3` or `F_3^4` contains a vector of weight at least 3 (direct check over the 40 hyperplanes of `F_3^4`):
type `(3)^3` or `(3)^4`. Finally let `|G| = 5 * 3^j` (`j >= 2`). The number of Sylow 5-subgroups divides `3^j` and is
`1 mod 5`, so it is 1 or 81. If 81, `G` has a 3-subgroup of order at least 81 and the 3-group case applies. If the Sylow
5-subgroup `<a>` is normal and `a` is a single 5-cycle, a Sylow 3-subgroup acts on `<a>` through `Aut(Z_5) = Z_4`, hence
trivially; its elements centralise `a`, so they fix the support of `a` pointwise (the centraliser of a 5-cycle in its
support is `<a>`), and an element of order 3 times `a` has type `(5)(3)` or `(5)(3)(3)`. If `a` has type `(5)(5)` we are
done. `QED`

For a strong `S` on 12 vertices with `tc(S) in {7, 21}`, `|Aut(S)| >= f(12)/21 = 375/21`, so `|Aut(S)| >= 19`. The scans of
the status line (all invariant tournaments of each type, Gray-code over pair orbits, exact Held-Karp `H`, nauty `|Aut|`
whenever `H = 7a` or `21a` with `a` odd, `a >= 19`, `ord(sigma) | a`) find no hit. With THM-4790 (a), no tournament on at
most 12 vertices has `tc in {7, 21}`.

## (b) The orbit product lemma

For any partition `V = B_1 u ... u B_r`, the OCF (THM-002) gives `H(T) >= prod H(T[B_i])`: collections of disjoint odd
cycles each inside one block are independent in the conflict graph, and their total weight is the product of the block
OCF sums; every other collection adds a positive weight. Apply this to the orbits of `G = Aut(T)`. `G` preserves each orbit
`O` and acts on `T[O]` by automorphisms; since `G` acts faithfully on `V`, it embeds in the product of its restrictions,
so `|G| <= prod |Aut(T[O])|`. Dividing,

    tc(T) >= prod_(O in V/Aut(T)) tc(T[O]),                                    (b)

and each `T[O]` is vertex-transitive (hence regular, of odd order).

## (c) The vertex-transitive lemma

**Lemma.** A vertex-transitive tournament `T` of order `m` has `tc(T) >= 23` unless `T` is `K1, C3, R5` or `P7`
(`tc = 1, 1, 3, 9`).

*Proof.* `Aut(T)` is transitive of odd order.
- `m` prime: a transitive group of prime degree contains a `p`-cycle, so `T` is a circulant. The Cayley census gives the
  minimum `tc` at `m = 7, 11, 13, 17` as `9` (Paley; the other circulant class has `tc = 25`), `1729`, `94949`,
  `792878221`.
- `m = 9`: the complete census of the 15 regular 9-tournaments has exactly three vertex-transitive classes, with `tc =
  363, 373, 39` (the last is `C3[C3]`, `|Aut| = 81`).
- `m = 15`: an odd-order group is solvable (Feit-Thompson) and a solvable primitive group has prime-power degree, so the
  group is imprimitive with blocks of size 3 or 5. With five blocks of size 3 the action on blocks is `Z_5` (the only
  odd-order transitive group of degree 5), the kernel `K` acts on each block as `Z_3`, and transporting labels by a Sylow
  5-subgroup `P` identifies `K` with a nonzero `F_3[Z_5]`-submodule of `F_3^5` projecting onto every coordinate: the
  diagonal, the augmentation module `A_3`, or everything (`x^4+x^3+x^2+x+1` is irreducible over `F_3`). The diagonal and
  the full module give a regular `Z_15`; `K = A_3` gives `G_A`. Three blocks of size 5 give likewise the diagonal, all, or
  the augmentation module `A_5` of `F_5^3` (`x^2+x+1` irreducible over `F_5`), i.e. `G_B`. So `T` is a Cayley tournament
  of `Z_15` (minimum `tc` 162165) or is `G_A`- or `G_B`-invariant; the 8 + 8 invariant tournaments of those groups are
  `R5[C3]` (`tc` 162165) and `C3[R5]` (`tc` 515241).
- `m >= 19`: in the OCF sum keep only collections of disjoint directed triangles. A vertex-transitive tournament has
  `c3 = m(m^2-1)/24` triangles, `(m^2-1)/8` through each vertex, and each triangle meets at most `D = 3((m^2-1)/8 - 1)`
  others, so `i_k >= prod_(i<k) (c3 - i(D+1))/k!` and `H >= sum_k 2^k i_k`. With Dixon's bound
  `|Aut| <= 3^((m-1)/2)` this exceeds `23 |Aut|` for every odd `m` in `[19, 401]` (exact integers), and for `m > 401` the
  single term `k = floor(m/20)` suffices: there `c3 - k(D+1) >= 0.0229 m^3` and `2^k i_k >= (0.0458 m^3)^k/k!`, whose
  logarithm exceeds `log 23 + ((m-1)/2) log 3` with a margin increasing in `m` (33 at `m = 401`). `QED`

## (d) The window, and its closure

Let `T` have `tc(T) in {7, 21}`; by THM-4790 (a) a strong component `S` has `tc(S) in {7, 21}`. By (b) and (c) every
orbit of `Aut(S)` induces `K1`, `C3`, `R5` or `P7`, and the product of their `tc` values (`1, 1, 3, 9`) is at most
`tc(S)`. Hence `G = Aut(S)` is a subgroup with transitive projections of `P = Z_3^c x Z_5^r x (Z_7 x| Z_3)^p` (the product
of the orbit groups `Aut(C3)`, `Aut(R5)`, `Aut(P7)`), with `3c + 5r + 7p <= |S|`, `(r, p)` in `{(0,0), (1,0)}` for
`tc = 7` and in `{(0,0), (1,0), (2,0), (0,1)}` for `tc = 21`, and `|G| >= f(|S|)/tc(S)` because `H(S) >= f(|S|)`
(`f(3k) = 3 5^(k-1)`, `f(3k+1) = 5^k`, `f(3k+2) = 9 5^(k-1)`; Busch). The inequality `f(n) <= tc |P|` fails at `n = 21`,
for every `n >= 23`, and for `tc = 7` at every `n >= 13` (`type_lemma_checks.py`); with (a) only `tc(S) = 21` and
`13 <= |S| <= 20` or `|S| = 22` remain. There the index `[P : G] <= 21|P|/f(|S|)` is at most 6, and being odd it is 1 or
3. An index-3 subgroup of an odd-order group is normal (3 is the least prime divisor), hence the kernel of a surjection
`P -> Z_3`, which factors through `P^ab = Z_3^(p+c) x Z_5^r`. This leaves 24 orbit configurations and 83 subgroups. For
index 1 every orbit is a module (the rotation of one orbit alone is an automorphism), so `S` is a substitution of the
blocks into a quotient; for index 3 the between-orbit orbitals may split. Fix the standard labelled copies of the orbit
tournaments (any `S` can be relabelled so) and enumerate every orientation of every pair of mutually inverse
between-orbit orbitals of `G`: the 185,642 strong tournaments obtained have no `tc in {7, 21}` (status line). Hence no
tournament has `tc` equal to 7 or 21.

## (e) Prime untouchables

If `p >= 7` is prime and `p - 1 = q + r` with distinct primes `q, r` (`q = 2` allowed), then `s(qr) = 1 + q + r = p`. So a
prime untouchable `p > 5` has either `p - 1` outside the Goldbach set, or `p - 1 = 2q` as its only representation. The
first set has `O(x^(1-delta))` elements up to `x` (Montgomery-Vaughan); the second lies in `{p = 2q + 1 : q prime}`, which
has `O(x/log^2 x)` elements (upper-bound sieve). Hence the prime untouchables up to `x` number `O(x/log^2 x) = o(pi(x))`:
almost all primes are touchable, unconditionally. That 2 and 5 are the only ones remains the conditional statement of
THM-4805.

## (f) Strong presentation fibres are finite

Let `T` be strong with `tc(T) = v`, with `Aut(T)`-orbit pieces `P_i = T[O_i]` (vertex-transitive) and `t_i = tc(P_i)`.
By (b), `prod t_i <= v`. A piece with `t_i = 1` is `K1` or `C3` (THM-4741: `tc = 1` exactly for `K1`/`C3`-stacks, and a
vertex-transitive stack has one block). The other pieces have `t_i >= 3`, so there are at most `log_3 v` of them, each of
order at most `M(v)`, the largest order of a vertex-transitive tournament with `tc <= v`; `M(v)` is finite because
`tc >= greedy(m)/3^((m-1)/2)` for vertex-transitive tournaments of order `m >= 19` (proof of (c)), and the greedy bound
grows faster than any exponential. Dixon's bound on those pieces and `|Aut(C3)| = 3` give
`|Aut(T)| <= 3^(n/3) v^(M(v)/6)`, while `H(T) >= f(n) >= 5^(-1/3) 5^(n/3)`. So `v = H/|Aut| >= 0.5848 (5/3)^(n/3) v^(-M(v)/6)`,
i.e. `n <= 3 ln(v^(1 + M(v)/6)/0.5848)/ln(5/3)`. Since `M(v) = O(log v/log log v)`, the bound is `O(log^2 v)`. (For
`v = 21`, `M = 7` and the bound is 41; the sharper accounting of (d) gives 22.) The full fibre of `v` (all tournaments)
is infinite, because stacking `K1` and `C3` components does not change `tc`; finiteness is a statement about strong
classes.

## (g) The strong fibres below 22

For odd `v <= 21` the counting of (d), with the vertex-transitive pieces `K1, C3, R5, P7` (the only ones with
`tc <= 21`) and the product constraint `prod tc(T[O]) <= v`, bounds the order of a strong `S` with `tc(S) = v` (by 12
for `v = 5`, 16 for `v = 9`, 19 for `v = 11, 13, 15`, 22 for `v = 17, 19, 21`), and the admissible subgroups again have
index 1 or 3 (index 5 would need a factor 5 of `|P|` in a configuration where it does not occur). The generalised search
over orders `>= 11` (567 cases, 811,496 strong invariant tournaments) finds no strong tournament with `tc <= 21`, and
the strong census of orders `<= 10` gives the fibres:

| `v` | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| strong classes with `tc = v` | 2 | 1 | 2 | 0 | 5 | 4 | 1 | 6 | 6 | 4 | 0 |
| orders | 1, 3 | 5 | 4, 5 | - | 5-7 | 5-7 | 5 | 6-7 | 6-8 | 6-7 | - |

(`v = 1`: `K1`, `C3`; `v = 3`: `R5`; per-order counts in `strong_fibre_table.out`.) Every strong presentation fibre
below 22 is therefore finite, explicit, and supported on orders at most 8.
