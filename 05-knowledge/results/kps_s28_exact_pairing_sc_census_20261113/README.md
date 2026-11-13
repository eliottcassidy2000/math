# Exact pairings of merged classes, with self-converse structure, through order 9

kind-pasteur-2026-11-13-S28 · [THM-5153](../../../01-canon/theorems/THM-5153-exact-slot-pairings-of-merged-tournament-classes-and-self-converse-structure-to-order-nine.md) · **Chan–Pak remains open.**

The owner asked to find exact pairings, including self-converse classes, toward a complete spectrum and the
A000568 recurrence. Sessions 314 and 315 hold the merged-metagraph recurrences and generating functions, the
pointed Burnside recurrence for A000568, the self-converse gaps below 75, the square-sum ports, and the
Monster/Ogg filters. This note supplies a complementary census that neither has.

An **exact pairing** is two merged classes (tournaments up to converse) with the same full slot polynomial Q.
Since Q(T^op) = Q(T), a converse pair is never counted as a pairing, and every pairing is invisible to every uniform
exterior (Session 310 NEUTRAL_TAILS).

## Census (FINITE-EXACT, every class of orders 3-9)

| n | classes (A000568) | self-converse (A002785) | merged | Q-classes | Q-classes with >= 2 merged classes | of which SC-SC | of which SC-NSC | strong-twin Q-classes | with an SC member | SC classes alone in their Q-class |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| 4 | 4 | 2 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 2 |
| 5 | 12 | 8 | 10 | 8 | 2 | 1 | 1 | 1 | 1 | 5 |
| 6 | 56 | 12 | 34 | 24 | 9 | 3 | 2 | 5 | 3 | 5 |
| 7 | 456 | 88 | 272 | 156 | 69 | 15 | 28 | 49 | 27 | 37 |
| 8 | 6,880 | 176 | 3,528 | 1,728 | 827 | 34 | 69 | 707 | 73 | 53 |
| 9 | 191,536 | 2,752 | 97,144 | 50,511 | 21,217 | 158 | 1,175 | 19,528 | 1,161 | 1,348 |

**How to read the table.**

- "Merged" equals (A000568 + A002785)/2, a built-in check.
- A Q-class is the set of merged classes sharing Q. "SC-SC" means the Q-class contains at least two self-converse
  classes. "SC-NSC" means it contains a self-converse class and a converse pair.

**Findings.**

- The full slot polynomial does not determine self-converseness: SC-NSC pairings exist from order 5.
- The Q-class counts 2, 3, 8, 24, 156, 1,728, 50,511 were recorded by the S26 audit, and 24 at order 6 by Session
  261 MODULE_REBASE. The merged and self-converse refinement is new.

## Two cross-checks with today's packets

- **Session 315's self-converse gaps.** The odd values below 75 attained by no self-converse class of order <= 9 are
  exactly 7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63. This is Session 315's global list, here from an
  independent class generator (first_creators classgen) for the orders where it is finite.
- **S26's exact pairs.** Among the seven strong members of strong/decomposable exact pairs through order 9
  (THM-5151 (4)), exactly one is self-converse: the order-8 class with H = 45 (classgen code 167804930, rigid),
  whose partner is C3 => Y5'. The other six form three converse pairs: H = 15, H = 95 (Session 310's B and its
  converse) and the other two H = 45 members.

## Replay

Run `python scpairs.py` in a scratch directory holding `cls_n.txt` (classgen output: code, H, c3, strong, SC, aut,
converse code) and `prof_n.txt` (fprof profiles) for n = 3..9. The kps_s26 packet's `scripts/README.md` explains how
to build both.
