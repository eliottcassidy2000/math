---
id: HYP-9467
title: "Pre-registered test of THM-4874 at orders 10-15. For every sole-exit-free tournament T with hc(T) = k in {3, 4, 5}, let k1 = #{v : hc(T - v) = 1} and k2 = #{v : hc(T - v) = 2}. THM-4874 predicts the joint counts of (k1, k2, |Aut|) exactly from the block of T's terminal core: hc = 3 gives (3, 0, 3) for F(2n - 10) classes (C3 block) and (3, 1, 1) for F(2n - 11) (W4 block); hc = 4 gives (2, 1, 1) and (2, 2, 1) for F(2n - 11) each; hc = 5 gives (2, 0, 1) for F(2n - 12) ((W4, all)), and (3, 0, 1), (2, 1, 1), (2, 2, 1) for F(2n - 13) each (the three 5-blocks, the first being (AT_5, its two sole exits)); nothing else occurs"
status: CONFIRMED (pre-registered in 671f7fb87e; all predictions at orders 10-15 held exactly, no other combination occurred); recorded in THM-4874
session: kind-pasteur-2026-11-04-S10
date: 2026-11-04
related: [THM-4874, HYP-9464, THM-4869, THM-4872, THM-4868]
script: 04-computation/experiments/fewcycle_deletion_census_20261104.c, 04-computation/experiments/fewcycle_block_deletions_20261104.py
---

# HYP-9467 - the uniquely Hamiltonian deletion law from the core block: order-15 test

## The claim being tested (THM-4874, reserved; proof to be filled this session)

Let T have no vertex of out-degree 1 and exactly k Hamiltonian cycles, where k = 3 and n >= 6, or k = 4, 5 and
n >= 7. THM-4869 (C2) writes T = Q[z <- (D, u)]: a uniquely Hamiltonian rooted prefix (or one vertex l), and a core
D = chain w_1 => ... => w_p => block (W, A) plus u, with only the arc u -> l leaving D. Then:

| v | hc(T - v) |
|---|---|
| v in Q - z, v != l | k * (paths of Q - z - v from l), a multiple of k |
| v = u | 0 |
| v = l | 0, unless T - D = {l}: then k (p >= 1) or h0(W, A) (p = 0) |
| v = w_i (chain) | k |
| v in W | delta(v) = sum over x in A - v of end_x(W - v) |

Only W contributes values 1 and 2, apart from the minimal case T - D = {l}, p = 0, which needs n = |W| + 2 <= 7.
The block data (`fewcycle_block_deletions_20261104.out`, all blocks of type <= 5, re-deriving THM-4869 B5) are:

| type | block (W, A) | delta over W | k1 | k2 | h0 | classes (THM-4869 C4 cores) |
|---|---|---|---|---|---|---|
| 3 | (C3, all) | 1, 1, 1 | 3 | 0 | - (p >= 1) | F(2n - 10); Aut = Z_3 (THM-4872 (B)) |
| 3 | (W4, {c, d}) | a 1, b 2, c 1, d 1 | 3 | 1 | 2 | F(2n - 11) |
| 4 | (W4, {a, c, d}) | a 1, b 3, c 2, d 1 | 2 | 1 | 2 | F(2n - 11) |
| 4 | (W4, {b, c, d}) | a 1, b 2, c 2, d 1 | 2 | 2 | 1 | F(2n - 11) |
| 5 | (W4, all) | a 1, b 3, c 3, d 1 | 2 | 0 | - (p >= 1) | F(2n - 12) |
| 5 | (AT_5, its sole exits {4, 5}) | 1 at 1, 4, 5; 3 at 2, 3 | 3 | 0 | 4 | F(2n - 13) |
| 5 | (W_5, A_1), W_5 with scores 1, 2, 2, 2, 3 | 1, 1, 2, 3, 3 | 2 | 1 | 4 | F(2n - 13) |
| 5 | (W_5, A_2), same W_5 | 1, 1, 2, 3, 2 | 2 | 2 | 4 | F(2n - 13) |

## Pre-registered predictions (orders 10-15)

`fewcycle_deletion_census_20261104.c` is the HYP-9464 auditor's independent census `aud_fewcycle.c` (every strong
class with hc <= 5, by vertex extension), changed only by NM = 15 and the joint statistic. Run with `15`. For each
order n = 10, ..., 15 its `JOINT` lines will show exactly these counts and no other (k1, k2, |Aut|) combination:

| n | hc 3: (3,0,3) | hc 3: (3,1,1) | hc 4: (2,1,1) | hc 4: (2,2,1) | hc 5: (2,0,1) | hc 5: (3,0,1) = (2,1,1) = (2,2,1) |
|---|---|---|---|---|---|---|
| 10 | 55 | 34 | 34 | 34 | 21 | 13 each |
| 11 | 144 | 89 | 89 | 89 | 55 | 34 each |
| 12 | 377 | 233 | 233 | 233 | 144 | 89 each |
| 13 | 987 | 610 | 610 | 610 | 377 | 233 each |
| 14 | 2,584 | 1,597 | 1,597 | 1,597 | 987 | 610 each |
| 15 | 6,765 | 4,181 | 4,181 | 4,181 | 2,584 | 1,597 each |

So u_3(15) = 10,946, u_4(15) = 8,362 and u_5(15) = 7,375, and the deletion law holds at order 15:
- 3 deletions for hc = 3;
- 2 deletions for hc = 4;
- for hc = 5, 3 deletions for F(17) = 1,597 classes and 2 for L(18) = 5,778.

**What is blind.**
- Order 15 is new for everything.
- At orders 10-14 the k1 counts and |Aut| were known from the HYP-9464 audit. The k2 counts and the joint
  correlation are new there.

**Exploratory record.** A smoke test to order 9 ran before this commit (`fewcycle_deletion_census_20261104_explore.out`)
and matches the table's formulas at orders 7-9 (n = 6 is the stated minimal exception). Every outcome will be
reported, with the class if a prediction fails.

## Results

**CONFIRMED: every prediction holds, at every order 10-15, with no other (k1, k2, |Aut|) combination.**
- The run was `fewcycle_deletion_census_20261104 15`, made after the pre-registration commit 671f7fb87e. It took
  16 min wall time, 826 s of it at order 15.
- The record is `fewcycle_deletion_census_20261104.out`. Its `JOINT` lines were compared mechanically with the table
  above (18 rows, all MATCH).

Order 15 has 1,090,032 strong classes with hc <= 5:
- u_3, u_4, u_5 = 10,946, 8,362, 7,375;
- t_1 = 46,368 = F(24), Garey's uniquely Hamiltonian count.

Its joint counts:
- hc 3: (3,0,3) 6,765 and (3,1,1) 4,181;
- hc 4: (2,1,1) 4,181 and (2,2,1) 4,181;
- hc 5: (2,0,1) 2,584, and (3,0,1), (2,1,1), (2,2,1) 1,597 each.

So the deletion law holds at order 15: hc 3 gives 3 deletions, hc 4 gives 2, and hc 5 gives 3 for 1,597 classes and 2
for 5,778 = L(18).
