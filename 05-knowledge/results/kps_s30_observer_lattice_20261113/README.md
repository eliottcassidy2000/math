# Constraining the merged metagraph from several observers at once

kind-pasteur-2026-11-13-S30 · [THM-5155](../../../01-canon/theorems/THM-5155-score-sequence-and-arc-usage-multiset-pin-merged-classes-except-symmetric-core-flips.md) · [HYP-9681](../../hypotheses/HYP-9681-score-and-arc-usage-coincidences-are-symmetric-score-module-flips.md) · **Chan–Pak remains open.**

The owner asked for complete understanding of the merged metagraph by constraining it from all angles. Sessions
314, 315, 319 and 321 hold its recurrences, the Fano wiring, BSD parity and selective returns, and Session 322 holds
the prescribed-target inverse.

This packet measures how much of a merged class (a tournament up to converse) each standard observer recovers,
alone and in pairs, over every class of orders 3-9.

## The observer lattice (FINITE-EXACT)

| n | merged classes | H | scores (symmetrised) | slot polynomial Q | arc multiset | (H, C, \|Aut\|, c3) | scores + Q | scores + arcs |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 10 | 7 | 7 | 8 | 10 | 8 | 9 | 10 |
| 6 | 34 | 19 | 14 | 24 | 34 | 27 | 27 | 34 |
| 7 | 272 | 77 | 33 | 156 | 271 | 162 | 182 | 272 |
| 8 | 3,528 | 320 | 76 | 1,728 | 3,526 | 1,633 | 2,093 | 3,528 |
| 9 | 97,144 | 1,520 | 184 | 50,511 | 97,121 | 29,742 | 63,686 | 97,138 |

How to read the table:
- C is the number of Hamiltonian cycles, and c3 the number of directed triangles.
- "Symmetrised" scores fold s with n − 1 − s. The exact converse-invariant observer, the pair {scores(T),
  scores(T^op)}, gives the same last column.
- The arc multiset is THM-5154's per-arc Hamiltonian-path count multiset, zeros included.

## The two observers are almost complete together (THM-5155)

Each alone fails. The arc multiset cannot tell run slides and run-flanked core flips apart (THM-5154). The score
sequence cannot separate most classes.

Together they separate every merged class through order 8. At order 9 exactly six pairs remain. Each is
TT2 => S => TT1 against TT2 => S^op => TT1, where S is a rigid, non-self-converse strong 6-tournament whose score
sequence reads the same backwards:
- (1,2,2,3,3,4) four times;
- (2,2,2,3,3,3) once;
- (1,1,2,3,4,4) once.

This is forced. Slides shift the core's scores by the number of singletons after it, so the score sequence sees
them. A flip maps the core's internal scores s to |S| − 1 − s, which is invisible exactly when that multiset is
self-reversed.

HYP-9681 conjectures that this is the only failure mode at every order: reversing a module whose own score
sequence is self-reversed and whose seam weights are preserved. Session 320's strong order-24 pair, which reverses a
regular vertex-transitive module, fits the pattern.

## Replay

In a directory with classgen's cls_n.txt, fprof's prof_n.txt and kps_s29's arcms_n.txt, run
`python scripts/observers.py` and `python scripts/score_arcs.py`.
