---
id: HYP-9458
title: "Cut arcs of the unique Hamiltonian cycle: the number T(n, k) of uniquely Hamiltonian tournaments of order n whose unique Hamiltonian cycle has exactly k cut arcs (arcs whose deletion destroys strong connectivity) satisfies T(n, 2) = 2^(n-5), T(n, 3) = n 2^(n-7), T(n, 4) = (n^2 + n - 8) 2^(n-10) and T(n, n-1-j) = C(n-5+j, j) for j = 0, 1, 2; pre-registered at orders 14, 15, 16 together with an independent (Shelor-free) generation check of Garey's counts"
status: OPEN (proof) - P1-P6 CONFIRMED at orders 14, 15, 16 (opus-2026-11-04-S74; run after the pre-registration commit 525d6c2d2e)
session: opus-2026-11-04-S74
date: 2026-11-04
related: [THM-4867, THM-4866, THM-4862, THM-4846, HYP-9457]
note: 05-knowledge/results/uniquely_hamiltonian_one_arc_cut_shelor_20261104.md
script: 04-computation/experiments/uh_one_arc_cut_20261104.c
---

# HYP-9458 - cut arcs of uniquely Hamiltonian tournaments

## Setting

U(n) is the set of isomorphism classes of tournaments of order n with exactly one directed Hamiltonian cycle C. Its size
is Garey's F(2n - 6).

An arc e of C is a cut arc if T - e is not strong. Equivalently, e is the only backward arc of a vertex cut.

**Classical fact.** Shelor 1977 gives every T in U(n) a vertex of in-degree 1, so cut arcs exist.

**Definition.** T(n, k) is the number of classes in U(n) with exactly k cut arcs on C.

## Exploratory evidence (orders 5-13, before this commit)

Program `uh_one_arc_cut_20261104.c`. It generates U(n) independently of Shelor's grammar, by Moon deletion and
cyclic-interval extension. Its counts are 1, 1, 3, 8, 21, 55, 144, 377, 987, 2584, 6765 at orders 3-13, which are
Garey's numbers.

Rows of T(n, k), for k = 2, 3, ...:

| n | T(n, 2), T(n, 3), ... |
|---:|---|
| 5 | 1, 1, 1 |
| 6 | 2, 3, 2, 1 |
| 7 | 4, 7, 6, 3, 1 |
| 8 | 8, 16, 16, 10, 4, 1 |
| 9 | 16, 36, 41, 30, 15, 5, 1 |
| 10 | 32, 80, 102, 85, 50, 21, 6, 1 |
| 11 | 64, 176, 248, 231, 155, 77, 28, 7, 1 |
| 12 | 128, 384, 592, 608, 456, 259, 112, 36, 8, 1 |
| 13 | 256, 832, 1392, 1560, 1289, 819, 406, 156, 45, 9, 1 |

## Fitted closed forms

- T(n, 2) = 2^(n-5), for n >= 5;
- T(n, 3) = n 2^(n-7), for n >= 6;
- T(n, 4) = (n^2 + n - 8) 2^(n-10), for n >= 7;
- T(n, n-1) = 1, T(n, n-2) = n - 4, and T(n, n-3) = C(n-3, 2).

## Pre-registered predictions

Confirmatory run: `uh_one_arc_cut 16`.

- **(P1) Sanity: the generator.** |U(14)| = 17,711, |U(15)| = 46,368 and |U(16)| = 121,393 (Garey). This checks an
  independent generator; it is not a prediction.
- **(P2) Sanity: classical structure.** Every class at orders 14-16 has a vertex of in-degree 1 and a vertex of
  out-degree 1, by Shelor and converse symmetry. Hence every class has a one-arc cut.
- **(P3) Two cut arcs.** T(14, 2) = 512, T(15, 2) = 1,024, T(16, 2) = 2,048.
- **(P4) Three cut arcs.** T(14, 3) = 1,792, T(15, 3) = 3,840, T(16, 3) = 8,192.
- **(P5) Four cut arcs.** T(14, 4) = 3,232, T(15, 4) = 7,424, T(16, 4) = 16,896.
- **(P6) Top diagonals.**

  | n | T(n, n-1) | T(n, n-2) | T(n, n-3) |
  |---:|---:|---:|---:|
  | 14 | 1 | 10 | 55 |
  | 15 | 1 | 11 | 66 |
  | 16 | 1 | 12 | 78 |

Any failure is reported with the observed value.

## Results (2026-11-04, after the pre-registration commit 525d6c2d2e)

The run was `uh_one_arc_cut 16` (4 min). Its output is `04-computation/experiments/uh_one_arc_cut_20261104.out`.

- **(P1)** The class counts are 17,711, 46,368 and 121,393, which are Garey's numbers.
- **(P2)** Every class has a vertex of in-degree 1 and a vertex of out-degree 1.
- **(P3)-(P6) CONFIRMED:**

| n | T(n, 2) | T(n, 3) | T(n, 4) | top diagonals |
|---:|---:|---:|---:|---|
| 14 | 512 | 1,792 | 3,232 | 1, 10, 55 |
| 15 | 1,024 | 3,840 | 7,424 | 1, 11, 66 |
| 16 | 2,048 | 8,192 | 16,896 | 1, 12, 78 |

**Open.** A proof of the closed forms, and the bivariate generating function of T(n, k). The natural route is a transfer
matrix on Shelor's five-letter port grammar (operation_closed_synthesis section II).
