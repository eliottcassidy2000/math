---
id: HYP-9658
title: "Unit-offset seeds of the 2-adic lines of the Hamiltonian-path spectrum: (A) every odd c has a rooted tournament with offset c0 = 1 and multiplier c 2^j, so that every line {1 + c 2^s} is eventually in the spectrum by one explicit ray; (B) the eleven odd c <= 129 without a seed of order <= 9 (19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115) acquire one by order 11; (C) every odd number of binary weight three other than 7 and 21 is a path count; (D) the number of unit-offset roots of order n grows geometrically (2, 5, 7, 16, 31, 68, 207, 1333 for n = 2..9)."
status: >
  OPEN (conjectures (A)-(D)). FINITE-EXACT support: all rooted classes of order <= 9 (seeds.c on the class lists of
  classgen.c); 636 odd c have a seed of order <= 9, among them 54 of the 65 odd c <= 129. PRE-REGISTERED tests P1-P4
  and gate G1 for the order-nine run (pre-registration commit 2c27353a14 on origin/main before the run): P2 and P3
  CONFIRMED, P1 and P4 FAILED, G1 holds (section "Outcomes"). (A) is neither implied by nor implies the Chan-Pak
  conjecture: a seed gives a ray of a line, not the finitely many points below it, and Chan-Pak gives every point
  but no ray. The Chan-Pak conjecture remains OPEN.
session: kind-pasteur-2026-11-10-S21
source: 05-knowledge/results/kps_s21_two_adic_lines_20261110/README.md
scripts:
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/scripts/seeds.c, finalize.py
  - 05-knowledge/results/kps_s21_two_adic_lines_20261110/unit_offset_roots_orders_2_to_9.txt, line_seeds.txt
related:
  - 01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-5085-*, THM-5109-* (fixed-root tails, dyadic cones), THM-4888-*
  - 05-knowledge/hypotheses/HYP-9656-response-bounded-self-feeding-closures-are-finite.md
---

# HYP-9658 - unit-offset seeds of the 2-adic lines

Notation of THM-5126: for a rooted tournament `(L, v)`, `c0 = H(L - v) - delta'` and `M = H - E + delta'`, and
`H(L o 0^t) = c0 + M 2^t` for `t >= 1`. A *seed* of the line `L_c = {1 + c 2^s : s >= 1}` is a unit-offset root
(`c0 = 1`) with `M = c 2^j`, or a doubling vertex with `H - 1 = c 2^s0`.

## Conjectures

- **(A) Seeds exist.** Every odd `c >= 1` has a seed. Equivalently: for every odd `c` the line `L_c` is in the
  spectrum from some point on, by one explicit ray.
- **(B) The eleven.** Each of `c = 19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115` has a unit-offset root of order
  at most eleven.
- **(C) Weight three.** Every number `1 + 2^a + 2^b` (`a > b >= 1`) other than `7` and `21` is a Hamiltonian-path
  count. (Rows `a - b in {1, 2, 3, 4, 6, 7}` are THM-5126 (B4). This is a consequence of Chan-Pak.)
- **(D) Growth.** The number `u(n)` of unit-offset roots of order `n` (rooted isomorphism classes counted as
  class-vertex pairs) satisfies `u(n + 1) >= 3 u(n)` for `n >= 7`; observed `2, 5, 7, 16, 31, 68, 207, 1333`.

## What a proof of (A) would and would not give

A seed for `c` with multiplier `c 2^j` puts `1 + c 2^s` into the spectrum for `s >= j + 1`. The `j` points below
are a separate finite question. So (A) is a statement about the *2-adic ends* of the spectrum: every line has an
end inside `S`. Chan-Pak in addition needs the first point above every missing prefix, uniformly in `c`. The
least `s0(c)` found is `<= 6` for `c <= 129`; no bound in `c` is conjectured.

A unit-offset root has `delta' = H(L - v) - 1` crossed covers. For one-apex roots this forces a fan (`c = 1`). For
a root with `H(L - v) = E` (every Hamiltonian path of `L - v` ends at an in-neighbour of `v`) it gives
`M = H - 1`: the ray is the doubling tower `1 + (H - 1) 2^t`. 116 of the 336 unit-offset roots of order `<= 8`
have `d = E` (197 of the 1,333 at order nine); 8 of the 336 are trivial (`M = 0`).

## Pre-registered tests (order nine; text fixed before the run in PREREGISTRATION.md)

- P1: the number of odd `c` with a tail seed of order `<= 9` lies in `[250, 420]`.
- P2: at least 6 of the ten odd `c <= 63` without a tail seed at order `<= 8`
  (`9, 13, 19, 27, 33, 35, 41, 53, 59, 63`) acquire one at order nine.
- P3: `c = 13` acquires a tail seed at order nine.
- P4: `c = 33` acquires a tail seed at order nine.
- G1 (gate, forced by THM-1370): no ray contains 7 or 21.

## Outcomes

- P1 **FAILED**: 636 odd `c` (the count grew by a factor 7.1 from order eight, not 3.5).
- P2 **CONFIRMED**: 8 of 10 (`9, 13, 27, 35, 41, 53, 59, 63`); `19` and `33` remain.
- P3 **CONFIRMED** (`c = 13`: order nine, `(H, d, E, U) = (183, 61, 35, 269)`, ray `1 + 13 * 2^s`, `s >= 5`).
- P4 **FAILED**: no unit-offset root of order `<= 9` has `M = 33 * 2^j`.
- G1 holds: the least `s0` is `2` for `c = 3` and `3` for `c = 5`.

Unit-offset roots at order nine: 1,333 among 1,723,824 class-vertex pairs. Odd `c <= 129` without a seed after the
run: the eleven of (B).

## Next tests (not run)

- Order ten (9,733,056 classes): (B) for at least six of the eleven; `u(10) >= 4000`.
- A structural source of seeds: describe the unit-offset roots with `d = E` (doubling towers) as a grammar.
