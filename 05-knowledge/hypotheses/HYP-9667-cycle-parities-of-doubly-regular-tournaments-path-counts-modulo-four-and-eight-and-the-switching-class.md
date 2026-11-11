---
id: HYP-9667
title: "Cycle counts of doubly regular tournaments modulo 2, their path counts modulo 4 and 8, and what happens on the switching class. (bK) For every doubly regular tournament of order n and every k, c_(2k+1) = [n = 3 mod 8] C((n-1)/2, k) (mod 2) (PROVED for k <= 3, and for every k when the tournament has a converse symmetry: THM-5136 (E)); equivalently, for n = 7 (mod 8) its numbers of paths of every order are those of the transitive tournament modulo 4. (b4) H = n(n-1)/2 (mod 4). (b8) H = n(n-1)/2 (mod 8). REFUTED by the independent audit: the same statements for every member of the switching class of a skew Hadamard matrix (S^2 = -(N-1) I); they are true for every such tournament of order at most 16 and false at order 24. (QF) OBSERVATION: on every switching class of tournaments, H mod 4 and the parity of each cycle count are quadratic functions of the switching set."
status: >
  OPEN (CONJECTURE with FINITE-EXACT evidence) for (bK), (b4), (b8) on doubly regular tournaments; REFUTED for
  the switching-class forms (bS), (bK-class); OBSERVATION for (QF) and (bT).
  PROVED cases of (bK): k = 1, 2 (every tournament of order N = 0 mod 4 whose scores have one parity, THM-5132);
  k = 3 for every doubly regular tournament (c_7 = n(n+1)(n-1)(n-3)(n-5)(n^2-6n+1)/896 is divisible by 8; the count
  is due to Savchenko 2017 and is re-derived in the source note, section 5.4); every k for doubly regular
  tournaments with a converse symmetry, in particular all Paley tournaments (THM-5136 (E), (P)). So (b4) is
  equivalent to: the number of odd cycles of length at least 9 of a doubly regular tournament is even.
  FINITE-EXACT for (bK), (b4), (b8): every doubly regular tournament of order at most 23 - by the classification
  1, 1, 1, 2, 2, 37 tournaments at orders 3, ..., 23 (Spence; McKay's catalogue), and the audit's nauty census shows
  that the recorded classes of order 23, their converses and the descendants of their matrices are all 37. At
  order 27 (722 tournaments in the literature): (bK), (b4) for P_27 and all 39 Goethals-Seidel classes met, 36 of
  them without converse symmetry; (b8) for ten of them (direct counts).
  REFUTED (audit/AUDIT.md of the source packet; re-checked with the author's programs): (bS) "S^2 = -(N-1) I
  implies H = 1 + N/2 (mod 4)" and the class form of (bK). Counterexample: (source) + D23_0 switched at one vertex
  has H = 3 (mod 4) and an odd number of 23-cycles; 18 of its 23 single-vertex switchings violate the class form
  of (bK). Both had been pre-registered and CONFIRMED (R11, R12) on classes that turned out to be special.
  WHERE THE CLASS FORMS HOLD (FINITE-EXACT): for every tournament with S^2 = -(N-1) I and N <= 16 (all members of
  all classes; audit); in every test at N = 20 (all single and double switchings and about 200 random members of
  each of the two classes); at N = 24 on the classes of P_23 and of the double of P_11 (all 300 single and double switchings of each; random
  members) and on none of the twelve annealed classes. No criterion is known: "a single class of descendants"
  holds for the good classes and also for D23_8, whose class is not constant.
  (QF) FINITE-EXACT: all labelled tournaments of order 4, 5, 6 (degree exactly 2 at orders 5, 6; identically zero
  at order 4); NUMERICAL: random tournaments of orders 7-20 and members of skew Hadamard classes (third
  differences); the per-length form at orders 5-8. Not proved. For N <= 5 it follows from the Ising form of H on a
  switching class (THM-4524 (A3) (d)).
  Pre-registered: R9, R10 (non-Paley doubly regular tournaments of orders 19, 23, 27), R11 (switching classes of
  P_23^+, D(P_11)^+, D_19^+), R12 (cycle counts by length) CONFIRMED as registered; the order-19 and order-23 parts
  of R9, R10 were decidable from the repository, and the general conjectures recorded from R11, R12 are the ones
  the audit refuted (MISTAKE-763).
  REFUTED local forms (each true for small orders only): det S = 4H - 3 (mod 16) (N <= 4; the order-4 case is the
  recorded |Pf| = H mod 4); "one row of S orthogonal to the others keeps H mod 4 under switching at that vertex"
  (N <= 8, exhaustive; false at N = 10 for 1152 of 3840 pairs); "equal score parities keep it" (N <= 6); the
  vertex-local form of (bK) (false for CD_4, D_19).
  PROVED: the reduction (H(U) - H(U^v))/2 = sum over x in N^-(v) of g_R(x) (mod 2); no tournament with
  S^2 = -(N-1) I, N = 0 (mod 4), is self-converse.
  RELATION: for Paley tournaments (b8) is HYP-9273 (B) (opus S33), which is open; (b4) for Paley is THM-5136 (P).
  NOT CLAIMED: that the tournaments listed are all doubly regular tournaments of their orders.
  AUDIT: the independent audit (05-knowledge/results/character_of_two_records_20261111/audit/AUDIT.md) REFUTED the
  first form of this hypothesis - the switching-class statements (bS) and (bK-class), counterexamples at N = 24
  (its M1) - and confirmed the doubly regular form on all 37 tournaments of order 23 (its own census) and in
  every recomputation; this file was rewritten accordingly and the minor points 1, 2, 3, 6, 7, 8, 14 applied.
  MISTAKE-763.
session: mac-mini-2026-11-11-atoms
source: 05-knowledge/results/character_of_two_records_20261111/README.md
related:
  - THM-5136
  - THM-5137
  - THM-5132
  - THM-5127
  - THM-4524
  - THM-448 (drt-mersenne-doubling-tower)
  - THM-451
  - THM-4557
  - THM-4911
  - HYP-9365
  - HYP-9371
  - HYP-9663
---

# HYP-9667 - cycle parities of doubly regular tournaments; switching classes

Source: [the note](../../05-knowledge/results/character_of_two_records_20261111/README.md), sections 4, 5.2-5.4, 10.

## Statements

`T` is a doubly regular tournament of order `n = 4t + 3 = 2m + 1`; `c_l` is its number of directed `l`-cycles, `p_j`
its number of directed paths on `j` vertices, `H = p_n`.

- **(bK)** `c_(2k+1)(T) = [n = 3 (mod 8)] C(m, k) (mod 2)` for every `k >= 1`: all odd cycle counts are even if
  `n = 7 (mod 8)`, and the parity vector is that of the binomial coefficients of `m` if `n = 3 (mod 8)`. This is
  THM-5136 (E) 2 without the converse symmetry.
  Equivalent form (subset identity of THM-5132): `p_j = C(n, j) (mod 4)` for every `j` if `n = 7 (mod 8)` (the path
  counts of the transitive tournament modulo 4), and `p_j = C(n, j) + 2 [j odd, j >= 3] C(m, (j-1)/2) (mod 4)` if
  `n = 3 (mod 8)`.
- **(b4)** `H(T) = n(n-1)/2 (mod 4)`. It follows from (bK) (`H = 1 + 2 c_odd (mod 4)`, `sum_k C(m, k) = 2^m - 1`).
  Since `c_3 = (n+1)/4`, `c_5 = 0`, `c_7 = 0 (mod 2)` are proved, (b4) says that the number of odd cycles of length at
  least 9 is even.
- **(b8)** `H(T) = n(n-1)/2 (mod 8)`. For Paley tournaments this is HYP-9273 (B). (The orbit form `H/|Aut| = 1 (mod 8)`
  of that hypothesis fails for the doubly regular `T_15`; the form here holds for it.) Not a term-by-term law: among the classes of order 23 the odd cycle counts
  modulo 4 differ at the lengths 19 and 23 and the parity of the number of disjoint pairs of odd cycles differs,
  while `H = 5 (mod 8)` for all.
- **(c)** For `P_q`: the number `F` of invariant odd cycles is `0 (mod 4)` if `q = 7 (mod 8)` and `(q+1)/4 (mod 4)` if
  `q = 3 (mod 8)` (FINITE-EXACT for `q <= 47`).

## What is refuted: the switching class

Switching a tournament at a vertex set reverses the arcs between the set and its complement. The tournaments `U`
on `N` vertices with `S^2 = -(N-1) I` form classes closed under switching, and `(source) + T` with `T` doubly regular
are their normalised members (the descendants of the skew Hadamard matrix `I + S`).

- **(bS), REFUTED:** "`S^2 = -(N-1) I`, `N = 0 (mod 4)`, implies `H(U) = 1 + N/2 (mod 4)`".
- **(bK-class), REFUTED:** "`S^2 = -(N-1) I` implies `c_(2k+1)(U) = [N = 4 (mod 8)] C(N/2 - 1, k) (mod 2)`".

Counterexample (independent audit; `audit/counterexample_N24_outmasks.txt` of the source packet): for the annealed
doubly regular tournament `D23_0` of HYP-9371, `(source) + D23_0` switched at one vertex has `H = 3 (mod 4)`
against `1`, and its odd cycle parities are `00000000001` (lengths 3, ..., 23) against zero. Of the 23 single-vertex
switchings, 18 change the parity vector and 12 change `H mod 4`. The same happens on the other annealed classes
(45 of 60 random members violate the class form of (bK), 30 of 60 violate (bS)); only the lengths 19, 21, 23 ever
have the wrong parity. The 24 descendants of each of these matrices nevertheless satisfy (bK) and (b8).

Both statements were pre-registered by the author (R11, R12) and confirmed on the classes of `P_23^+`, `D(P_11)^+`,
`D_19^+`, `P_19^+`: the test objects were the most symmetric classes available, and the twelve generic classes of
order 24 in the repository were used only through their normalised members.

- **Where the class forms hold.** For every skew Hadamard tournament of order at most 16 (exhaustive; by the
  classification these are the classes of `P_3`, `P_7`, `P_11` and of `CD_4` with its converse); in every test at
  `N = 20` (`P_19`, `D_19`); at `N = 24` on the classes of `P_23` and of the double of `P_11`, and on none of the twelve
  annealed classes. No criterion is known ("one class of descendants" fails: `D23_8`), and no mechanism.

## (QF): the switching class carries a quadratic form (observation)

For a tournament `U` and a vertex set `X` let `U_X` be `U` switched at `X`.

- **(QF)** `X -> (H(U_X) - H(U))/2 mod 2` is a polynomial of degree at most 2 over `F_2`. Equivalently, for any three
  vertices the eight tournaments `U_S`, `S` a subset of the triple, satisfy `sum_S H(U_S) = 0 (mod 4)`.
- **(QF-length)** The same for `X -> c_l(U_X) mod 2`, for every length `l`.

For a single vertex of a tournament this is the quadratic form of a vertex cube (THM-5127); the switching cube is a
different cube (one coordinate per vertex, each move a complementation in a vertex cube). If (QF) holds, a class is
constant as soon as all single and double switchings of one member are; and (b4) says that the rows of a skew
Hadamard matrix are zeros of the quadratic form of its switching class.

## Evidence

Doubly regular tournaments (sc = with a converse symmetry, where (b4) and (bK) are PROVED):

| order | tournaments | `H mod 8` | odd cycle parities `(c_3, c_5, ...)` | (bK) |
|---|---|---|---|---|
| 3, 7, 11 | `P_q` (sc) | 3, 5, 7 | `1`; `000`; `10011` | the same |
| 15 | `CD_4` (not sc) | 1 | `0000000` | `0000000` |
| 19 | `P_19` (sc), `D_19` (sc) | 3, 3 | `100000011` | `100000011` |
| 23 | all 37: `P_23` (sc) and 36 without converse symmetry (the recorded classes, their converses, and the other descendants of their matrices) | 5 (all) | zero (all) | zero |
| 27 | `P_27` (sc); all 39 Goethals-Seidel classes met (3 sc, 36 not) | 7 (the ten counted) | `1001100110011` (all 40) | `1001100110011` |

Switching classes (members of the class of `(source) + T`):

| `N` | `T` | members | odd cycle parities / `H mod 4` |
|---|---|---|---|
| 4, 8, 12, 16 | `P_3`, `P_7`, `P_11`, `CD_4` | all (8, 128, 2048, 32,768; audit) | constant |
| 20 | `P_19`, `D_19` | all single and double switchings; about 200 random members of each | constant |
| 24 | `P_23`, double of `P_11` | all 300 single and double switchings of each; 40 random members each | constant |
| 24 | the annealed classes `D23_0, ..., D23_11` | single-vertex switchings; random members | **not constant** (audit) |

`H mod 8` is not constant even on the classes of `CD_4` and `D_19`.

Exact path counts at order 23: `04-computation/experiments/sierpinski_ocf_20261022/drt_results.out`,
`drt23_classes.out`, `drt23_doubled_check.out` (HYP-9371). Everything else: `out/` and `audit/` of the source packet.

## Local forms that fail (small-number laws)

| statement | true for | first failure |
|---|---|---|
| `det S = 4H - 3 (mod 16)` | all tournaments of order 2, 4 | order 6 |
| row `v` of `S` orthogonal to all other rows implies `H(U^v) = H(U) (mod 4)` | every tournament with such a row for `N <= 8` (exhaustive) | `N = 10` (1152 of 3840 pairs) |
| all vertices but `v` have scores of one parity implies the same | `N = 4, 6` | `N = 8` |
| the class form (bS) itself | every skew Hadamard tournament for `N <= 16`; every test at `N = 20` | `N = 24` |
| the parities of (bK) for the cycles through one vertex | Paley tournaments | `CD_4` (order 15) |

## Falsifiers

- A doubly regular tournament with one odd cycle count of the wrong parity refutes (bK); with
  `H != n(n-1)/2 (mod 4)`, (b4); `(mod 8)`, (b8). The smallest untested objects: the remaining classes of order 27;
  `CD_5` of order 31 (predicted `H = 1 mod 8`).
- A tournament and a triple of vertices with `sum_S H(U_S) != 0 (mod 4)` refutes (QF).

## Routes

0. The length 9. In the literature the number of `m`-cycles of a doubly regular tournament depends only on the
   order for `m <= 8` (Savchenko 2016, 2017, arXiv:2403.07629) and conjecturally for `m = 9`; here `c_3, ..., c_9`
   coincide for all 37 tournaments of order 23 (`c_9 = 97,774,380`; audit) and the two of order 19 and differ from
   length 10 on. A formula for `c_9(n)` (a paper on the 9-cycles of doubly regular tournaments is announced by
   Savchenko) would settle (bK) for `k = 4`: the prediction is `c_9` odd exactly for `n = 11 (mod 16)`.
1. A proof of (QF) (three vertices; eight switchings), then the algebra of the quadratic form on a skew Hadamard
   class: (b4) is the statement that the 2N vectors `+- row` are singular.
2. The Walsh expansion of `c_odd` in even path forests (THM-4524 (A3) on the switching cube): in a doubly regular
   tournament every sum with a free vertex of degree 2 is determined; what remains are the patterns without such
   a vertex, and the claim is that their total contribution is even.
3. The doubling `[[M, M], [-M^T, M^T]]` (THM-447, THM-448) in tournament form: the double of `T` lives on
   `V_1 + {oo} + V_2` with `T` on `V_1`, the converse of `T` on `V_2`, `oo -> V_1`, `V_2 -> oo`, `a_1 -> a_2`, and
   `a_1 -> b_2` iff `a -> b`. (b4) for doubles is the case `n = 7 (mod 8)`. For an arbitrary base tournament the
   residue of the double is not determined (both residues occur from base order 5 on).
