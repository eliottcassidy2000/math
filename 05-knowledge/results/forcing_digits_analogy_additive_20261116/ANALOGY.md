# 13 and 41 against 7 and 21: the analogy made as strong as the arithmetic allows

Session klein-2026-11-16-force. Builds on wave 50 of chanpak_waves_20261113 (mac-mini-2026-11-16-tower), which is
the source of the owner's observations and is cited, not repeated; new content is the charge and digit census of
FORCING.md read on these numbers. Receipts `digits_4..9.out`, `strongB_3..9.out`.

## 1. What wave 50 proved (cited)

In count coordinates `H = 2x + 1` the shortcut Collatz inverse has the even child `2H − 1` for every `H` and the odd
child `(2H − 1)/3` for the forks `H ≡ 5 (mod 6)`; the five states read modulo 18, and the first sterile fork
(`H ≡ 11 (mod 18)`, odd child a multiple of 3 in `x`) is `H = 11`. **The holes are exactly the two children of the
first sterile fork: 21 = 2·11 − 1 and 7 = 21/3**, no other count has a hole as a child, and the next sterile forks 29,
47, 65, … have children (57, 19), (93, 31), (129, 43), … that are all counts. The even children of the holes are
13 = 2·7 − 1 and 41 = 2·21 − 1, the ready seeds where the rows `L_3` and `L_5` restart; among forks up to 45 those with
no charge-0 vertex in one-parent growth (`n ≤ 10`) are 11, 29, 35 — THM-5272's list of doubling failures.

## 2. The rows and the tree are the same object (PROVED, from COLLATZ_ROWS and THM-5259)

`L_3 = {7, 13, 25, 49, …}` is the sterile ray through the hole 7 (`x = 3·2^k`, never a fork); `L_5 = {11, 21, 41,
81, …}` is the ray of the first sterile fork 11 (`x = 5·2^k`), on which 21 is the doubling child and 41 the
grandchild — a fork again (`x = 20 ≡ 2 (mod 3)`) whose odd child is `27` (`(2·41 − 1)/3`), THM-5272's arithmetic
number. So: the hole 7 heads a ray that never forks; the hole 21 sits one step below a fork that does; and "the rows
through the holes start up again" at 13 and 41 is the literal statement that the next rungs of `L_3` and `L_5` are
counts.

## 3. What the charge census adds (FINITE-EXACT; THM-5260 (B) for the all-order statements)

With the padding laws of FORCING.md §2 every count below 125 with one nontrivial strong component has its full
all-order charge and digit spectrum determined by its strong realisers of order ≤ 9:

| count | digit 0 | digit 1 (any order) | digit 2 | rooted doubling child | role |
|---|---|---|---|---|---|
| 7 | — | — | — | — | hole |
| 21 | — | — | — | — | hole |
| 11 | 0 | 1 | 1 | 0 | first sterile fork; doubling fails (THM-5272) |
| 13 | 0 | 1 | 1 | 1 | even child of 7; restart of L_3; ready seed |
| 41 | 1 | 0 | 0 | 1 | even child of 21; restart of L_5; ready seed |
| 29 | 1 | 1 | 1 | 0 | sterile fork; doubling fails |
| 35 | 1 | 0 | 0 | 0 | hidden generator 5·7; doubling fails |
| 47 | 0 | 1 | 1 | 0 | sterile fork |
| 65 | 1 | 1 | 1 | 1 | sterile fork |
| 83 | 1 | 0 | 0 | 1 | sterile fork |
| 101 | 1 | 0 | 0 | 0 | sterile fork |
| 119 | 1 | 1 | 1 | 1 | sterile fork |
| 19 | 1 | 1 | 1 | 0 | odd child of 29 |
| 57 | 1 | 0 | 0 | 1 | even child of 29 |
| 27 | 1 | 0 | 0 | 1 | odd child of 41 |

So the restarts 13 and 41 differ: 13 has digits 1 and 2 and a rooted doubling child but no digit-0 vertex (it is one of the four counts 11, 13, 47, 67 with no `B_t = h`); 41 has digit 0 and a rooted doubling child but neither digit 1 nor digit 2 — `3·41 + 2 = 125` and `127` are not tower-forced from 41 at any order. The first sterile fork 11 has every digit but no rooted doubling child; the sterile forks 29, 47, 101 likewise lack the rooted doubling child, while 65, 83, 119 have one.

## 4. The analogy, stated as strongly as the data allow

- **Same shape, different mechanism.** 7, 21 are the two children of the sterile fork 11 and are holes by THM-1370
  (counts multiply over strong components and 7, 21 fall between consecutive strong minima). 13, 41 are the doubling
  children of the holes; they are counts, they are singleton-strong (13: THM-4884 (2)) or near-singleton (41: one class
  at order 6), and they carry ready roots. The analogy "13, 41 : 7, 21" is the shift by one rung along `L_3` and
  `L_5`: the hole is where the rung is missing, the restart is where the count-doubling `2H − 1` becomes physical again.
- **Sterility is about odd children, doubling failure about even children.** The sterile forks `11, 29, 47, 65, 83,
  …` (`H ≡ 11 mod 18`) are defined by their odd child being a multiple of 3 in `x`; THM-5272's doubling failures
  `11, 29, 35` are about the even child `2H − 1` not containing a realiser of `H`. They coincide at 11 and 29 and part
  at 35 (not a sterile fork; `35 = 5·7` is a hidden generator, THM-5258 (B)) and at 47 (a sterile fork where doubling
  succeeds). The rooted-append census (charge 0) gives a third list that contains both: counts below 125 with no rooted doubling child at any order are 11, 19, 29, 31, 35, 39, 43, 47, 53, 55, 61, 63, 67, 71, 73, 77, 85, 87, 101, 113 (THM-5260 (C) with the order-9 strong census), of which only 11, 29, 47, 101 are sterile forks and only 11, 29, 35 are containment failures. The strongest true statement is wave 50's: the holes are the two children of the first sterile fork, and the restarts are the holes' doubling children; the three lists (sterile forks, containment failures, rooted-append failures) agree at 11 and 29 and nowhere else below 125.
- **P111 of wave 50 in the rooted-append sense.** P111 asks that for every prime `p ≡ 1 (mod 4)` other than 13 and 41 the count `(p + 1)/2` have a realisation where the doubling `D` is physical. In the rooted-append sense (a charge-0 root at some order, which by THM-5260 (C) is decided by the strong realisers of orders ≤ 9 for `(p + 1)/2 < 125`) this **fails at p = 37, 61, 109, 173** (`(p + 1)/2 = 19, 31, 55, 87` have no rooted doubling child at any order); P111 can only hold with a weaker notion of "physical doubling" (any one-vertex extension of a realiser, THM-5272's containment), which the port spectra of THM-5272 (E) say does succeed at 19 and 31.
