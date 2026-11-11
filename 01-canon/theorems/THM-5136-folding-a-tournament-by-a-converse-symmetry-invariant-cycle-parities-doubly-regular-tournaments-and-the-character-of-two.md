---
id: THM-5136
title: "Folding a tournament by a converse symmetry. Let T have an anti-involution psi with fixed vertex v_0 (odd order 2m+1). The fold is the tournament Z on the m pairs {a, psi(a)} (i -> j iff exactly one of a_i^+ -> a_j^+, a_i^+ -> a_j^- holds, a^+ the out-neighbour of v_0 in its pair) with the set D_1 of pairs whose triangle with v_0 is cyclic. PROVED: (Q) the number of psi-invariant cycles of length 2k+1 has the parity of the number of k-vertex paths of Z ending in D_1; hence H(T) = 1 + 2 #{paths of Z ending in D_1} (mod 4). (E) In a doubly regular tournament of order n = 4t+3 with a converse symmetry every invariant triple is cyclic if t is even and none is if t is odd; hence H = n(n-1)/2 (mod 4) and c_(2k+1) = [n = 3 mod 8] C((n-1)/2, k) (mod 2). (P) Paley tournaments: H(P_q) = 1 (mod 4) iff 2 is a square in GF(q); the fold of P_q is P_q on the squares."
status: >
  PROVED: (Q), (E), (P) and the corollaries (source note, sections 3, 4, 5.2). (Q): the half-paths of an invariant
  cycle are counted by transfer matrices M_ij = [[x, y], [1 - y, 1 - x]], which are z I + g 1^T modulo 2; the terms
  with cuts cancel by Redei's theorem (the number of Hamiltonian paths through the arc between two vertices has the
  same parity for both orientations) and, from three cuts on, by permuting segments. (E): the set of vertices x with
  a -> x -> psi(a) is psi-invariant and misses v_0, so it is even, and has t or t + 1 elements.
  Machine checks: (Q) on every normalised labelled tournament with an anti-involution of order 3, 5, 7, 9, 11 (2, 16,
  512, 65,536, 33,554,432 tournaments; every k), H = 1 + 2F (mod 4) by direct count for orders <= 9; all 126 circulant
  tournaments of orders 3-13 against a direct count of H; (E) on the non-Paley doubly regular tournament D_19
  (exact H) and on three non-Paley self-converse classes of order 27 (direct counts mod 256); (P) for
  q = 3, 7, 11, 19, 23, 27, 31, 43, 47 (invariant cycles of every length modulo 256) and against the known H for q <= 27.
  KNOWN, credited: the half-path form of invariant cycles and H = 1 + 2F (mod 4), H = 1 + 2F + 4I (mod 8), the
  remark that F is not decided by triangles alone, conjugacy of anti-involutions (Session29 half-path and
  converse-congruences reports, 2026-10-22); c_3 = F_1 (mod 2) (THM-552); H = 1 (mod 4) for self-converse tournaments
  of even order (THM-791); H = 1 + 2 c_odd (mod 4) (THM-466, THM-478); Redei's theorem; p | H(P_p), free action of
  Aut on Hamiltonian paths, the values H(P_q) for q <= 23 and the orbit numbers 1, 9, 1729, ... (THM-212, LEM-003,
  THM-5101, THM-4747); D_19 and its path count (HYP-9365, this machine, 2026-10-22); maxima of H at orders <= 9
  (THM-338).
  NEW here (the audit searched the repository as well): the fold (Z, D_1), statement (Q) and its corollaries, (E),
  (P) including the residue of H(P_q) modulo 4 for every prime power q = 3 (mod 4), H(P_27) modulo 256.
  NOT CLAIMED: novelty in the literature; anything about tournaments without a converse symmetry (HYP-9667); the
  Chan-Pak conjecture, which remains OPEN.
  Pre-registered: R1 (F(P_27), F(P_43) odd) CONFIRMED and now a corollary; R2, R3 FAILED and were replaced by (Q);
  R8 CONFIRMED (H(P_27) = 7 mod 8). MISTAKE-PENDING (entry filed at close-out) (two registrations were decidable from the repository).
  AUDIT: AUDIT-PENDING (pre-audit checkpoint)
session: mac-mini-2026-11-11-atoms
source: 05-knowledge/results/character_of_two_records_20261111/README.md
related:
  - THM-5137
  - HYP-9667
  - THM-5131
  - THM-5132
  - THM-552 (c3-parity-dichotomy-self-converse)
  - THM-791
  - THM-4930
  - THM-212 (paley-divisibility)
  - THM-5101
  - THM-448 (drt-mersenne-doubling-tower)
  - THM-451
  - HYP-9365
  - HYP-9371
---

# THM-5136 - folding a tournament by a converse symmetry

Source: [the note](../../05-knowledge/results/character_of_two_records_20261111/README.md), sections 3, 4, 5.2.
Independent audit: AUDIT-PENDING (pre-audit checkpoint)

## Setting

`T` is a tournament on `2m + 1` vertices with an *anti-involution* `psi` (`a -> b` iff `psi(b) -> psi(a)`,
`psi^2 = 1`); `v_0` is its fixed vertex. In each pair `{a, psi(a)}` exactly one vertex is an out-neighbour of `v_0`;
call it `a^+`, and `a^- = psi(a^+)`. Number the pairs `1, ..., m`.

- `D_1` = the set of pairs `i` with `a_i^+ -> a_i^-`, i.e. the invariant triple `{v_0, a_i^+, a_i^-}` is cyclic.
- `Z` = the digraph on the pairs with `i -> j` iff exactly one of `a_i^+ -> a_j^+`, `a_i^+ -> a_j^-` holds.
- `F_k` = the number of directed cycles `C` of length `2k + 1` with `psi(C) = C`; `F = sum_k F_k`.

All anti-involutions of `T` are conjugate under `Aut(T)` (KNOWN), so `(Z, D_1)` is an invariant of a self-converse
tournament of odd order.

## Statement (Q): the folding formula

1. `Z` is a tournament.
2. `F_k = #{directed paths of Z on k vertices that end in D_1} (mod 2)` for every `k >= 1`.
3. `c_odd(T) = F = #{directed paths of Z that end in D_1} (mod 2)` and
   `H(T) = 1 + 2 #{directed paths of Z that end in D_1} (mod 4)`.

Corollaries. (a) `D_1` empty: every `F_k` is even, `H = 1 (mod 4)`. (b) `D_1` everything: `F_k = C(m, k) (mod 2)`,
`F` is odd, `H = 3 (mod 4)`. (c) For fixed `psi` and fixed arcs between different pairs, exactly half of the `2^m`
choices of `D_1` give `H = 1 (mod 4)`. (d) Circulant tournament on `Z_n` with connection set `S`, `psi(x) = -x`: the
pairs are the elements of `S`, `D_1 = {a in S : -2a in S}`, `a -> b` in `Z` iff exactly one of `b - a`, `-b - a` is
in `S`.

## Proof of (Q)

Write `x_ij = [a_i^+ -> a_j^+]`, `y_ij = [a_i^+ -> a_j^-]`. Then `x_ij + x_ji = 1` and `y_ij = y_ji` (apply `psi`), so
`z_ij = x_ij + y_ij (mod 2)` is a tournament.

An invariant cycle of length `2k + 1` has the unique form
`v_0 -> x_1 -> ... -> x_k -> psi(x_k) -> ... -> psi(x_1) -> v_0` with `(x_1, ..., x_k)` a directed path meeting each
pair at most once (KNOWN). So `x_1` is the `+` vertex of its pair, and `x_k` is the `+` vertex iff its pair is in `D_1`.
For a sequence `(i_1, ..., i_k)` of distinct pairs the number of such paths is

    e_+^T M_(i_1 i_2) ... M_(i_(k-1) i_k) f,      M_ij = [[x_ij, y_ij], [1 - y_ij, 1 - x_ij]],

rows and columns indexed by the signs `+, -`, `f = e_+` if `i_k` is in `D_1` and `e_-` otherwise. Modulo 2,
`M_ij = z_ij I + g_ij 1^T` with `g_ij = e_+` if `y_ij = 1` and `e_-` if `y_ij = 0`. Expanding the product: the term
without `g` is `prod z * [i_k in D_1]`; a term with cuts at `r >= 1` steps is `y(first cut) * prod_(other steps) z`.
Summed over all sequences, the first kind is the number of `k`-vertex paths of `Z` ending in `D_1`. The second kind
counts lists `(Q_0, ..., Q_r)` of disjoint paths of `Z` covering a `k`-set, with weight `y(end Q_0, start Q_1)`. For
`r >= 3` permute `Q_2, ..., Q_r`: a multiple of `(r - 1)!`. For `r = 1, 2` (fix `Q_2`): let `B(u, w)` be the number of
pairs `(Q_0, Q_1)` covering a vertex set `V'` with `Q_0` ending at `u` and `Q_1` starting at `w`. `B(u, w)` is the number
of Hamiltonian paths through the arc `u -> w` of the tournament obtained from `Z[V']` by orienting `{u, w}` as
`u -> w`; the Hamiltonian paths avoiding `{u, w}` are the same for both orientations and both totals are odd (Redei);
so `B(u, w) = B(w, u) (mod 2)`, and as `y` is symmetric these terms cancel in pairs. This proves 2. Statement 3
follows from `c_odd = F (mod 2)` and `H = 1 + 2 c_odd (mod 4)`. (b): the number of `k`-vertex paths of `Z` is
`sum_(|S| = k) H(Z[S]) = C(m, k) (mod 2)`. (c): the total number of paths of `Z` is `2^m - 1`, odd, so the linear form
`D_1 -> #{paths ending in D_1} mod 2` is not zero.

## Statement (E): doubly regular tournaments with a converse symmetry

Let `T` be doubly regular of order `n = 4t + 3` (scores `2t + 1`, any two vertices have `t` common out-neighbours)
with an anti-involution `psi`.

1. For every `a` with `v_0 -> a`: `a -> psi(a)` iff `t` is even. So `D_1` is everything if `n = 3 (mod 8)` and empty if
   `n = 7 (mod 8)`.
2. `c_(2k+1)(T) = F_k = [n = 3 (mod 8)] C((n-1)/2, k) (mod 2)` for every `k`, and `H(T) = n(n-1)/2 (mod 4)`.

Proof. For an arc `u -> w` of a doubly regular tournament `#{x : u -> x -> w} = t` and `#{x : w -> x -> u} = t + 1`.
The set `E = {x : a -> x -> psi(a)}` is `psi`-invariant and does not contain `v_0` (since `v_0 -> a`), so `|E|` is even;
`|E| = t` if `a -> psi(a)` and `t + 1` otherwise. Then (Q) (a), (b), and the cycles that are not invariant come in
pairs.

## Statement (P): Paley tournaments

Let `q = 3 (mod 4)` be a prime power, `P_q` the Paley tournament on `GF(q)` (`a -> b` iff `b - a` is a non-zero
square), `psi(x) = -x`, `chi` the quadratic character.

1. `a^+` is the square of its pair; `D_1` is everything if `chi(2) = -1` (`q = 3 mod 8`) and empty if `chi(2) = 1`
   (`q = 7 mod 8`) (directly: `a -> -a` iff `chi(-2) = 1`).
2. `a -> a^2` is an isomorphism from `Z` onto `P_q` restricted to the non-zero squares (`z_ab = 1` iff `b^2 - a^2` is a
   square).
3. `H(P_q) = 2 - chi(2) = q(q-1)/2 (mod 4)`; the number of Hamiltonian paths starting with the arc `0 -> 1` is
   `1 (mod 4)` (for prime `q`: the number of orbits under `Aut(P_q)`).

| `q` | `q mod 8` | fold `Z` | `H(Z)` | `(F_k mod 2)_k` | `F mod 4` | `H(P_q) mod 8` |
|---|---|---|---|---|---|---|
| 3 | 3 | a point | 1 | `1` | 1 | 3 |
| 7 | 7 | cyclic triangle on `{1, 2, 4}` | 3 | `000` | 0 | 5 |
| 11 | 3 | regular 5-tournament on `{1, 3, 4, 5, 9}` | 15 | `10011` | 3 | 7 |
| 19 | 3 | regular, order 9 | 3357 | `100000011` | 1 | 3 |
| 23 | 7 | order 11 | 92411 | all 0 | 0 | 5 |
| 27 | 3 | order 13 | 3703011 | `1001100110011` | 3 | 7 |
| 31 | 7 | order 15 | 196922055 | all 0 | 0 | |
| 43 | 3 | order 21 | | `C(21, k) mod 2` | 3 | |
| 47 | 7 | order 23 | | all 0 | 0 | |

`H(P_27) = 27 * 13 * h` with `h = 1 (mod 256)`, so `H(P_27) = 95 (mod 256)`.

## What is not true

- `F = (n-1)/2 (mod 2)` when every invariant triple is cyclic (pre-registered R3): `F` is odd at every order.
- `F_k(P_q)` even for `k >= 2` when `q = 3 (mod 8)` (pre-registered R2): `F_k = C((q-1)/2, k) (mod 2)`.
- The fold of `P_23` is not `P_11` (`H = 92411` against `95095`): no ladder of folds along `5, 11, 23, 47`.

## Falsifier

A tournament with an anti-involution violating (Q) 2 for some `k`; a doubly regular tournament with a converse
symmetry and `H != n(n-1)/2 (mod 4)`.
