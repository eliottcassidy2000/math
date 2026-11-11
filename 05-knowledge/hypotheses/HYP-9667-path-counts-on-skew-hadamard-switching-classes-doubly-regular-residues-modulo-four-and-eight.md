---
id: HYP-9667
title: "Path counts and cycle counts on skew Hadamard switching classes. (bK) If the +-1 skew matrix S of a tournament U on N vertices satisfies S^2 = -(N-1) I (I + S is a skew Hadamard matrix), then for every k the number of directed cycles of length 2k+1 has the parity of [N = 4 mod 8] C(N/2 - 1, k). (bS) Under the same hypothesis H(U) = 1 + N/2 (mod 4); equivalently H mod 4 is constant on the switching class of every skew Hadamard matrix. (b4) Every doubly regular tournament of order n has H = n(n-1)/2 (mod 4), i.e. its number of odd cycles of length at least 7 is even (PROVED when it has a converse symmetry: THM-5136 (E)). (b8) Every doubly regular tournament of order n has H = n(n-1)/2 (mod 8). (c) For the Paley tournaments the number F of invariant odd cycles is 0 (mod 4) if q = 7 (mod 8) and (q+1)/4 (mod 4) if q = 3 (mod 8)."
status: >
  OPEN (CONJECTURE with FINITE-EXACT evidence).
  PROVED special case of (b4): doubly regular tournaments with a converse symmetry, in particular all Paley
  tournaments (THM-5136 (E), (P)).
  (b4), (b8) FINITE-EXACT on every doubly regular tournament whose path count is on record: orders 3, 7, 11 (Paley);
  15 (CD_4 and its converse, THM-448); 19 (P_19 and D_19); 23 (P_23, twelve annealed classes and the doubled Paley
  class: exact counts of HYP-9371, all = 5 mod 8); 27 (P_27 and nine Goethals-Seidel classes, six of them without
  converse symmetry: direct counts modulo 256, H3-PENDING-ROW (run in progress at this checkpoint)).
  (bS) FINITE-EXACT: all members of the classes at N = 4, 8; random members at N = 12, 16, 20 (two classes), 24 (two
  classes). Its analogue modulo 8 is FALSE (N = 16, 20).
  (bK) PROVED for k = 1, 2 (every tournament of order N = 0 mod 4 whose scores have one parity, THM-5132) and for
  every k on doubly regular tournaments with a converse symmetry (THM-5136 (E)). FINITE-EXACT for every odd length:
  all members at N = 4, 8; 41 and 13 members at N = 12, 16; pre-registered R12: 42 tournaments at N = 20, the
  fourteen doubly regular classes of order 23 (twelve without converse symmetry) and 20 members at N = 24, ten
  tournaments of order 27 (six without converse symmetry), 86 of 86. The vertex-local form (cycles through one vertex) is FALSE.
  (c) FINITE-EXACT for q = 3, 7, 11, 19, 23, 27, 31, 43, 47.
  Pre-registered: R9, R10 (non-Paley doubly regular tournaments of orders 19, 23, 27), R11 (switching classes at
  N = 20, 24), R12 (cycle counts by length) CONFIRMED; the order-19 and order-23 parts of R9, R10 were decidable from the repository (MISTAKE-PENDING (entry filed at close-out)).
  REFUTED local forms (each a theorem for small orders only): det S = 4H - 3 (mod 16) (true for N <= 4); "one row of
  S orthogonal to the others keeps H mod 4 under switching at that vertex" (true for N <= 8 on every tournament
  tried, false from N = 10); "equal score parities keep it" (true for N <= 6, false at N = 8).
  PROVED reduction: (H(U) - H(U^v))/2 = sum over x in N^-(v) of g_R(x) (mod 2), R = U - v, g_R(x) = parity of the
  number of even-order paths of R with endpoint x; for |R| odd, of the paths on |R| - 1 vertices with endpoint x.
  NOT CLAIMED: that these are all doubly regular tournaments of the orders listed (the classification of skew
  Hadamard matrices of orders 16 to 28 is in the literature and was not checked here).
  AUDIT: AUDIT-PENDING (pre-audit checkpoint)
session: mac-mini-2026-11-11-atoms
source: 05-knowledge/results/character_of_two_records_20261111/README.md
related:
  - THM-5136
  - THM-5137
  - THM-5132
  - THM-448 (drt-mersenne-doubling-tower)
  - THM-451
  - THM-4557
  - THM-4911
  - HYP-9365
  - HYP-9371
  - HYP-9663
---

# HYP-9667 - path counts on skew Hadamard switching classes

Source: [the note](../../05-knowledge/results/character_of_two_records_20261111/README.md), sections 4, 5.2, 5.3, 10.

## Statements

A tournament `U` on `N` vertices has the skew matrix `S` (`S_ab = 1` if `a -> b`, `-1` if `b -> a`, `0` on the
diagonal). *Switching* `U` at a vertex set reverses the arcs between the set and its complement; it conjugates `S`
by a diagonal sign matrix. The tournaments with `S^2 = -(N-1) I` form classes closed under switching; a tournament
`T` on `n = N - 1` vertices is doubly regular iff `(source) + T` is in such a class.

- **(bK)** `S^2 = -(N-1) I` implies `c_(2k+1)(U) = [N = 4 (mod 8)] C(N/2 - 1, k) (mod 2)` for every `k >= 1`: all odd
  cycle counts are even if `N = 0 (mod 8)`, and the parity vector is that of the binomial coefficients of
  `m = N/2 - 1` if `N = 4 (mod 8)`. This is THM-5136 (E) 2 without the converse symmetry and on the whole switching
  class. (bK) implies (bS), because `H = 1 + 2 c_odd (mod 4)` and `sum_k C(m, k) = 2^m - 1` is odd.
  Equivalent form for a doubly regular `T` of order `n = 2m + 1` (by the subset identity of THM-5132), with `p_j` the
  number of directed paths on `j` vertices: `p_j = C(n, j) (mod 4)` for every `j` if `n = 7 (mod 8)` (the path counts of
  the transitive tournament modulo 4), and `p_j = C(n, j) + 2 [j odd, j >= 3] C(m, (j-1)/2) (mod 4)` if `n = 3 (mod 8)`.
- **(bS)** `S^2 = -(N-1) I` implies `H(U) = 1 + N/2 (mod 4)`.
- **(b4)** Every doubly regular tournament of order `n` has `H = n(n-1)/2 (mod 4)`: `3` if `n = 3 (mod 8)`, `1` if
  `n = 7 (mod 8)`. Since `c_3 = (n+1)/4` and `c_5 = 0 (mod 2)` for every doubly regular tournament (THM-5132), this is
  equivalent to: the number of odd cycles of length at least 7 is even. (bS) implies (b4).
- **(b8)** Every doubly regular tournament of order `n` has `H = n(n-1)/2 (mod 8)`.
- **(c)** For `P_q`: `F = 0 (mod 4)` if `q = 7 (mod 8)`, `F = (q+1)/4 (mod 4)` if `q = 3 (mod 8)` (with
  `H = 1 + 2F + 4I (mod 8)` of Session29, (b8) for Paley is then a statement about the parity of `I`).

`H mod 4` is not a switching invariant in general: it changes for about half of the random pairs at every order
from 5 to 10.

## What is proved

- (b4) for every doubly regular tournament with a converse symmetry (THM-5136 (E)): the invariant triples are all
  cyclic or all transitive according to `n mod 8`, and the folding formula does the rest.
- The reduction of the status line. So (bS) is equivalent to: in a skew Hadamard tournament, for every vertex `v`,
  the paths of `U - v` that miss exactly one vertex and have exactly one endpoint beaten by `v` are even in number.

## Evidence

| order | doubly regular tournaments (sc = with a converse symmetry) | `H mod 8` | `n(n-1)/2 mod 8` |
|---|---|---|---|
| 3, 7, 11 | `P_q` (sc) | 3, 5, 7 | 3, 5, 7 |
| 15 | `CD_4` (not sc) and its converse | 1 | 1 |
| 19 | `P_19` (sc), `D_19` (sc) | 3, 3 | 3 |
| 23 | `P_23` (sc); twelve annealed classes (not sc; two converse pairs); the double of `P_11` (not sc) | 5 (all fourteen) | 5 |
| 27 | `P_27` (sc); three sc and six non-sc Goethals-Seidel classes | H3-PENDING-ROW (run in progress at this checkpoint) | 7 |

| `N` | switching class of `(source) +` | members | `H mod 4` | `H mod 8` |
|---|---|---|---|---|
| 4 | `P_3` | all 8 | 3 | 3 |
| 8 | `P_7` | all 128 | 1 | 5 |
| 12 | `P_11` | 60 random | 3 | 7 |
| 16 | `CD_4` | 60 random | 1 | 1, 5 |
| 20 | `P_19` | 30 random | 3 | 3 |
| 20 | `D_19` | 40 random | 3 | 3, 7 |
| 24 | `P_23` | 40 random | 1 | 5 |
| 24 | double of `P_11` | 40 random | 1 | 5 |

Cycle counts by length (parity vectors `(c_3, c_5, ...)`; `cycle_lengths.c`, `cycle_parities.py`):

| `N` | tournaments | found | (bK) |
|---|---|---|---|
| 4, 8 | all members of the classes of `P_3^+`, `P_7^+` | `1`; `000` | `1`; `000` |
| 12 | 41 members of the class of `P_11^+` | `10011` | `10011` |
| 16 | `CD_4` and 12 members of its class | `0000000` | `0000000` |
| 20 | `D_19`, `P_19`, 20 random members of each class | `100000011` (42 of 42) | `100000011` |
| 24 | fourteen doubly regular classes of order 23 | zero (14 of 14) | zero |
| 24 | 10 random members of each of the classes of `P_23^+`, `D(P_11)^+` | zero (20 of 20) | zero |
| 28 | `P_27`, three self-converse and six other classes of order 27 | `1001100110011` (10 of 10) | `1001100110011` |

The even cycle counts have no such law (their parities differ among the classes of order 23). Nor is (b8) a
term-by-term law: among the fourteen classes of order 23 the odd cycle counts modulo 4 differ at lengths 19 and 23
(`c_odd = 0, 2, 4 (mod 8)` occur) and the parity of the number of disjoint pairs of odd cycles differs, while
`H = 5 (mod 8)` for all (`out/cycle_parities_i1mod4.out` of the source packet). The numbers of odd
cycles through a single vertex do not have the parities of (bK) (`CD_4`, `D_19`): the statement is global.

Exact path counts at order 23: `04-computation/experiments/sierpinski_ocf_20261022/drt_results.out`,
`drt23_classes.out`, `drt23_doubled_check.out` (HYP-9371). Everything else: `out/` of the source packet.

## Local forms that fail (small-number laws)

| statement | true for | first failure |
|---|---|---|
| `det S = 4H - 3 (mod 16)` (would give (bS)) | all tournaments of order 2, 4 | order 6 |
| row `v` of `S` orthogonal to all other rows implies `H(U^v) = H(U) (mod 4)` | 793 random tournaments of order 8 with such a row (724 not skew Hadamard); `(source) +` regular tournament of order 3, 5, 7, switched | `N = 10` (64 of 160 pairs differ), `N = 12, 14, 16` |
| all vertices but `v` have scores of one parity implies the same | `N = 4, 6` (all labelled tournaments `U - v`) | `N = 8` |

So a proof of (bS) must use the orthogonality of all rows at once.

## Falsifiers

- A doubly regular tournament with `H != n(n-1)/2 (mod 4)` refutes (b4) and (bS); with `H != n(n-1)/2 (mod 8)`, (b8).
  The smallest untested objects: the remaining classes of order 23 (the literature counts 37) and of order 27;
  `CD_5` of order 31 (predicted `H = 1 mod 8`; a direct count needs about 10 GB).
- A member of a skew Hadamard switching class with the other residue modulo 4 refutes (bS); one with a single odd
  length of the wrong parity refutes (bK).
- A prime power `q = 3 (mod 4)` violating (c).

## Routes

1. A converse symmetry up to switching (a signed permutation carrying `S` to `-S`). When its permutation part fixes
   no vertex, no normalised member need be self-converse; whether the folding argument of THM-5136 extends to this
   case was not examined.
2. The Walsh expansion of `c_odd` in even path forests: in a skew Hadamard tournament every sum with a free vertex
   of degree 2 vanishes; what remains are the patterns without such a vertex, and the claim is that their total
   contribution is even.
3. The doubling `[[M, M], [-M^T, M^T]]` (THM-447, THM-448) in tournament form: the double of `T` lives on
   `V_1 + {oo} + V_2` with `T` on `V_1`, the converse of `T` on `V_2`, `oo -> V_1`, `V_2 -> oo`, `a_1 -> a_2`, and
   `a_1 -> b_2` iff `a -> b` (checked against the matrix for `P_3, P_7, P_11`). (b4) for doubles is the case
   `n = 7 (mod 8)`. For an arbitrary base tournament the residue of the double is not determined (both residues
   occur from base order 5 on), so double regularity of the base has to enter.
