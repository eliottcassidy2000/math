---
id: THM-4817
title: "The order-12 ledger and the two after it: P(12) = 142,980,435,791 prime tournaments (A280772(12)), 1,170,949 self-converse primes, c(12) = 76,601,960,166 intransitivity classes, P_lab(12) = 68,487,766,269,146,864,640 labelled primes, l(12) = 35,524,962,710,023,558,019; the order-13 and order-14 ledgers from the order-13 types (P(13) = 46,453,867,277,461, P(14) = 27,687,212,017,612,715); every strong and labelled-strong count agrees with OEIS A051337 / A054946, and the order-14 Burnside defect 3,574,801,896/7 is pre-registered"
status: >
  FINITE-EXACT (mac-mini-2026-10-27-dipole; not independently audited at filing). The ledger at order n uses the prime
  cycle-index types of orders < n (THM-4816: klein's census through 10 reproduced, the tower at 11-13) and A(n), SC(n)
  computed by Burnside's lemma (equal to OEIS A000568 and A002785 through n = 15). Agreeing cross-checks: the strong
  counts S(n) equal OEIS A051337(n) and the labelled strong counts equal A054946(n) at n = 12, 13, 14 (consistency of the inputs only: S = A/(1 + A) is fixed
  by A000568 and S_lab = 1 - 1/T by 2^C(n,2), so these two agreements do not test the prime types); at n = 12 and 13
  the labelled-count identity of THM-4816 A5 holds against the full non-rigid census. The order-14 non-rigid sum was
  pre-registered; its 7- and 5-adic parts were then checked against the order-14 primes with an automorphism of order
  7 or 5 (98 and 4,200 classes; both agree). The 3-part is untested. Independently audited: PASS WITH CORRECTIONS (18
  edits applied; audit_thm4816_4817/REPORT.md).
session: mac-mini-2026-10-27-dipole
source: 05-knowledge/results/dipole_principle_and_prime_cycle_index_tower_20261027.md (section 7)
scripts:
  - 04-computation/experiments/prime_cycle_index_tower_20261027/tower.py (tower_n14.out)
  - 04-computation/experiments/prime_cycle_index_tower_20261027/burnside_counts.py
related:
  - 01-canon/theorems/THM-4816-prime-cycle-index-tower-order-eleven-twelve-thirteen-types-from-non-rigid-primes-and-the-ledger.md
  - 01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md
  - 01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md
  - 05-knowledge/hypotheses/HYP-9413-minimum-intransitivity-chain-family-fibonacci-primes-and-the-intransitivity-count-asymptotics.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
  - 05-knowledge/hypotheses/HYP-9423-dipole-endpoint-conjectures-two-connected-descent-and-the-prime-tower-predictions.md
---

# THM-4817 — The order-12, 13 and 14 ledgers

The identities are those of THM-4803 and THM-4814 (B), evaluated from the types of THM-4816:

- `A = SEQ(S)`, `S = X + sum_Q Z_{Aut Q}[A]` (all and strong classes);
- `SC(n) = SC_nonstrong(n) + SC_strong-nonprime(n) + SCP(n)` (self-converse classes; Burnside over `Aut±(Q)` with
  self-converse contents on the odd anti-automorphism cycles);
- `c = MSET(c_s)`, `c_s = X + sum_{Q up to duality} Z_{Aut±(Q)}[c]` (intransitivity classes, THM-4803);
- `T_lab = 1/(1 - S_lab)` with `T_lab = sum_{n >= 0} 2^C(n,2) x^n/n!`, labelled strong non-primes `n! sum_k P_lab(k)/k! [x^n] (T_lab - 1)^k`,
  and `l_s = x + (1/2) sum_k P_lab(k) l^k/k!`, `1 + l = exp(l_s)` (labelled intransitivity structures).

## The ledgers (FINITE-EXACT)

| quantity | n = 11 | n = 12 | n = 13 | n = 14 |
|---|---|---|---|---|
| all classes `A(n)` (A000568) | 903,753,248 | 154,108,311,168 | 48,542,114,686,912 | 28,401,423,719,122,304 |
| non-strong | 19,288,658 | 1,798,161,433 | 307,332,423,619 | 96,931,930,964,248 |
| strong non-prime | 87,588,161 | 9,329,713,944 | 1,780,914,985,832 | 617,279,770,545,341 |
| **prime `P(n)`** (A280772) | 796,876,429 | **142,980,435,791** | **46,453,867,277,461** | **27,687,212,017,612,715** |
| strong `S(n)` (= A051337) | 884,464,590 | 152,310,149,735 | 48,234,782,263,293 | 28,304,491,788,158,056 |
| self-converse `SC(n)` (A002785) | 279,968 | 1,492,288 | 95,458,560 | 872,687,552 |
| self-converse non-strong | 2,768 | 8,839 | 280,111 | 1,492,876 |
| self-converse strong non-prime | 17,893 | 312,500 | 2,182,256 | 103,568,601 |
| **self-converse primes `SCP(n)`** | 259,307 | **1,170,949** | **92,996,193** | **767,626,075** |
| strong intransitivity `c_s(n)` | 442,336,617 | 76,154,750,809 | 24,117,379,959,402 | 14,152,240,453,790,335 |
| **intransitivity classes `c(n)`** | 447,114,461 | **76,601,960,166** | **24,193,986,713,786** | **14,176,434,888,156,918** |
| labelled strong `s_lab(n)` (= A054946) | 35,262,092,417,856,512 | 72,926,863,133,112,198,144 | 300,318,571,786,159,783,496,704 | 2,467,430,973,323,656,141,183,549,440 |
| **labelled primes `P_lab(n)`** | 31,808,433,417,949,440 | **68,487,766,269,146,864,640** | **289,269,151,048,145,181,957,120** | **2,413,723,787,266,849,272,549,734,400** |
| **labelled intransitivity `l(n)`** | 16,901,802,400,214,078 | **35,524,962,710,023,558,019** | **147,811,543,365,879,760,556,872** | **1,222,222,358,932,702,742,749,843,027** |

The order-11 column is THM-4814 (B) (klein), reproduced here from the tower's own types; the bold entries of the other
columns are new. `P(n) - SCP(n)` is even in every column, as converse pairs require.

**Sequences.** Prime (simple) tournaments, OEIS A280772 (which stops at order 8): 1, 0, 3, 15, 197, 4008, 136,655,
7,906,161, 796,876,429, 142,980,435,791, 46,453,867,277,461, 27,687,212,017,612,715 (orders 3-14; orders 9-10 from
THM-4803 / THM-4814, 11 from THM-4814 (B)). Self-converse prime tournaments: 1, 0, 3, 3, 55, 72, 2,235, 5,525, 259,307,
1,170,949, 92,996,193, 767,626,075; labelled primes 2, 0, 264, 10,320, 973,680, 161,280,000, ... ; intransitivity
classes `c(n)`, `c_s(n)`, `l(n)` as above. None of these four is in the OEIS (search 2026-10-27).

**Ratios.** `P(n)/A(n)` = 0.882, 0.928, 0.957, 0.975 at n = 11-14 (THM-4803: tends to 1); `P_lab(n)/2^C(n,2)` = 0.883,
0.928, 0.957, 0.975; `c(n)/A(n)` = 0.4947, 0.4971, 0.4984, 0.4991 (THM-4803 and HYP-9413: tends to 1/2).

**Cross-check.** klein's THM-4818 (B) computed the order-12 column independently (the same identities, its own
order-11 types): `P(12)`, `SCP(12)`, `c_s(12)` and `c(12)` agree exactly. Its `P_lab(12)` and `l(12)` differ by
1,024 and 512 from the values above. The values here are exact, as confirmed by the labelled identity against the
order-12 types; THM-4818's `P_lab(12)` is the nearest IEEE double to the exact value, and its `l(12)`, computed exactly from that
`P_lab(12)`, inherits half of the error (MISTAKE-683).

## The Burnside defects

`n! P(n) - P_lab(n) = n! · sum over non-rigid prime classes of (1 - 1/|Aut|)`:

| n | sum of (1 - 1/\|Aut\|) | from the census |
|---|---|---|
| 11 | 1,337,728/165 = 8,107.44 | 12,128 · 2/3 + 23 · 4/5 + 3 · 10/11 + 54/55 (agrees) |
| 12 | 8,571,902/45 = 190,486.71 | 285,601 · 2/3 + 72 · 4/5 + 32 · 8/9 (agrees) |
| 13 | 4,389,979,016/585 = 7,504,237.63 | 11,255,339 · 2/3 + 432 · 4/5 + 368 · 8/9 + 5 · 12/13 + 38/39 (agrees) |
| 14 | **3,574,801,896/7 = 510,685,985.14** | prediction |

The order-14 value is pre-registered (HYP-9423 (6)). Its reduced denominator is 7, so some order-14 prime has an
automorphism group of order divisible by 7 (an element of type `7^2`), and, if every such group has order exactly 7,
their number is congruent to 6 modulo 7. Audit (2026-10-27, `04-computation/experiments/prime_cycle_index_tower_20261027/audit_thm4816_4817/`; the same as the later `out/nr14_p7.txt`,
`out/nr14_p5.txt`): the 8,192 invariant tournaments of type `7^2` give 98 prime classes, 93 with `|Aut| = 7` and 5 with
`|Aut| = 21`, so the hypothesis of the last clause fails; their 7-adic part `93 · 6/7 + 5 · 20/21 = 1774/21 ≡ 1/7`
(mod `Z_(7)`) agrees with `3,574,801,896/7`. The 524,288 invariant tournaments of type `1^4 5^2` give 4,200 classes, all
with `|Aut| = 5`, and `4,200 · 4/5` is an integer, as the denominator requires. The rest, `1,532,047,622/3`, comes from
the groups of 3-power order, so the order-14 non-rigid count is at most `98 + 4,200 + (3/2) · 1,532,047,622/3 = 766,028,109`. With nearly all non-rigid primes of `|Aut| = 3`, the order-14 non-rigid
count is about `(3/2) · 510,685,985 ≈ 7.7 · 10^8`.
