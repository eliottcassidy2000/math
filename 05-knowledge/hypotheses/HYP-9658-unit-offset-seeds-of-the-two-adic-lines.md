---
id: HYP-9658
title: "Unit-offset seeds of the 2-adic lines of the Hamiltonian-path spectrum: (A) every odd c has a rooted tournament with offset c0 = 1 and multiplier c 2^j - equivalently a tournament B and a vertex set O with W(O) = H(B) - 1 and U_B(O) = c 2^j - so that every line {1 + c 2^s} is eventually in the spectrum by one explicit ray; (B) the 33 odd c <= 4097 without a seed of order <= 10 (least 245) acquire one by order 12; (C) every odd number of binary weight three other than 7 and 21 is a path count; (D) the number of unit-offset roots grows geometrically with the order (2, 5, 7, 16, 31, 68, 207, 1333 class-vertex pairs for orders 2..9)."
status: >
  OPEN (conjectures (A)-(D)). FINITE-EXACT support: every rooted class of order <= 9 (seeds.c) and every vertex set of
  size 2-4 over every class of order 9 (unitsets.c, i.e. the roots of order 10 with out-degree 2-4): 5,529 odd c have
  a seed, among them all 65 odd c <= 129 and 2,016 of the 2,049 odd c <= 4097. PRE-REGISTERED tests, text pushed
  before each run (commits 2c27353a14 and 050845a98d on origin/main): order nine P2, P3 CONFIRMED, P1, P4 FAILED;
  order ten Q1-Q4 CONFIRMED; gates G1, G2 hold. (A) is neither implied by nor implies the Chan-Pak conjecture: a
  seed gives a ray of a line, not the finitely many points below it, and Chan-Pak gives every point but no ray. The
  earlier form of (B) (the eleven odd c <= 129 without a seed of order <= 9) is CONFIRMED by the order-ten run and is
  now part of THM-5126 (B3). The Chan-Pak conjecture remains OPEN.
session: kind-pasteur-2026-11-10-S21
source: 05-knowledge/results/kps_s21_two_adic_lines_20261110/README.md
scripts:
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/scripts/seeds.c, unitsets.c, unitpairs.py, finalize.py
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/unit_offset_roots_orders_2_to_9.txt, line_seeds.txt, below_rays.txt
related:
  - 01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-5122-shared-core-exchange-intervals-and-negative-inertia-capacity.md (the matrix W)
  - 01-canon/theorems/THM-5085-*, THM-5109-* (fixed-root tails, dyadic cones), THM-4888-*
  - 05-knowledge/hypotheses/HYP-9656-response-bounded-self-feeding-closures-are-finite.md
---

# HYP-9658 - unit-offset seeds of the 2-adic lines

Notation of THM-5126: for a rooted tournament `(L, v)` with `B = L - v` and `O = N+(v)`,
`H(L o 0^t) = c0 + M 2^t` for `t >= 1`, with `c0 = H(B) - W(O)` (`W(O)` the sum of the two-path-cover couplings
`W_ef` over the pairs of `O`) and `M = U_B(O)` (the sum of the responses). A *seed* of the line
`L_c = {1 + c 2^s : s >= 1}` is a unit-offset root (`W(O) = H(B) - 1`) with `U_B(O) = c 2^j`, or a doubling vertex
with `H - 1 = c 2^s0`.

## Conjectures

- **(A) Seeds exist.** Every odd `c >= 1` has a seed: there are a tournament `B` and a vertex set `O` with
  `W(O) = H(B) - 1` and `U_B(O) = c 2^j`. Equivalently, for every odd `c` the line `L_c` is in the spectrum from
  some point on, by one explicit ray.
- **(B) The thirty-three.** Each of the 33 odd `c <= 4097` without a seed of order `<= 10` (see README section 2
  for the list) has a seed of order at most twelve.
- **(C) Weight three.** Every number `1 + 2^a + 2^b` (`a > b >= 1`) other than `7` and `21` is a Hamiltonian-path
  count. (Rows `a - b <= 7` are THM-5126 (B4). This is a consequence of Chan-Pak.)
- **(D) Growth.** The number `u(n)` of unit-offset roots of order `n` (class-vertex pairs) satisfies
  `u(n + 1) >= 3 u(n)` for `n >= 7`; observed `2, 5, 7, 16, 31, 68, 207, 1333` for `n = 2..9`.

## What a proof of (A) would and would not give

A seed for `c` with multiplier `c 2^j` puts `1 + c 2^s` into the spectrum for `s >= j + 1`. The `j` points below
are a separate finite question. So (A) is a statement about the *2-adic ends* of the spectrum: every line has an
end inside `S`. Chan-Pak in addition needs every point below the first ray, uniformly in `c`. The least `s0(c)`
found is at most `8` for `c <= 129` (attained at `c = 19`); no bound in `c` is conjectured.

The condition `W(O) = H(B) - 1` is a single equation between a quadratic form of the indicator of `O` (the cut
matrix `W` of THM-5122, all of whose off-diagonal entries are even and positive) and the path count of the base.
For a transitive base it forces `|O| <= 1` (fans, `c = 1`). The roots found have out-degree at most `2` up to
order seven; out-degree `3` first occurs at order eight, and no unit set of size four exists over a base of order
nine. For a root with `H(B) = E` (every Hamiltonian path of `B` ends at an in-neighbour of `v`) the multiplier is
`M = H - 1` and the ray is the doubling tower `1 + (H - 1) 2^t`: 116 of the 336 unit-offset roots of order `<= 8`
are of this kind (197 of the 1,333 at order nine); 8 of the 336 are trivial (`M = 0`).

## Pre-registered tests and outcomes (texts fixed before each run in PREREGISTRATION.md)

Order nine (direct, every root of every class):
- P1 (the number of odd `c` with a tail seed of order `<= 9` lies in `[250, 420]`): **FAILED**, 636.
- P2 (at least 6 of `9, 13, 19, 27, 33, 35, 41, 53, 59, 63` acquire a seed): **CONFIRMED**, 8 (not `19, 33`).
- P3 (`c = 13`): **CONFIRMED**, `(H, d, E, U) = (183, 61, 35, 269)`.
- P4 (`c = 33`): **FAILED** at order nine.
- G1 (no ray through 7 or 21): holds.

Order ten (through the cut matrix, root out-degree 2-4; out-degree `>= 5` not searched):
- Q1 (at least 6 of the eleven odd `c <= 129` without a seed acquire one): **CONFIRMED**, all eleven.
- Q2 (`c = 19`): **CONFIRMED**, `s0 = 8`. Q3 (`c = 33`): **CONFIRMED**, `s0 = 6`.
- Q4 (at least 450 of the 513 odd `c <= 1025` have a seed): **CONFIRMED**, 510.
- G2: holds.

## Next tests (not run)

- Order eleven through the cut matrix (9,733,056 bases of order ten): (B); the count of unit sets of size four.
- A structural source of seeds: a grammar for the pairs and triples `O` with `W(O) = H(B) - 1`.
