# The tournament "magic measure": H-spread within a fixed score class (kps-S22, HYP-2707)

Operationalizing the Clifford/magic lens (HYP-2707) on the central invariant H = OCF.

DEFINITION. For a score sequence s, let spread(s) = max H − min H over all labeled tournaments with score s.
A score class with **spread 0** is a STABILIZER class: H is score-determined there, hence computable from the
degree-<=2 (Clifford / c3-level, THM-555) data. A class with **spread > 0** is MAGIC: H carries cycle-space
(degree>=3) content the scores cannot see.

EXACT (all labeled tournaments, n=4,5,6):
```text
n  #score-classes  #stabilizer(spread 0)  #magic  max H-spread  argmax score
4        4                4                  0          0         (all)
5        9                8                  1          4         (1,2,2,2,3)
6       22               13                  9         14         (1,2,2,3,3,4)
```

FINDINGS.
- **Magic onsets at n=5** (n=4 is ALL stabilizer — H fully score-determined = fully Clifford), and the magic
  fraction grows (0 -> 1/9 -> 9/22). This matches the c5/alpha_2 (degree-4) onset and the homology onsets
  (beta_3@n=6, the repo's "n=5/6 thresholds").
- **The H-maximizer (regular/Paley) is NOT the max-magic.** At n=5 the regular score (2,2,2,2,2) has spread 0
  (H=15 always — a STABILIZER class!), while the max-magic is the UNBALANCED (1,2,2,2,3) (spread 4). So
  "magic" (H-unpredictability from scores) is orthogonal to "H-large": the maximally-cyclic Paley tournament
  lives in a low-magic (nearly H-determined) score class.

CONSEQUENCE for H-computation. Stabilizer score classes (spread 0) admit a closed/score-level H value (Clifford-
cheap); magic classes need the full cycle-space. A poly-time pre-classifier "is this score class stabilizer?"
would route the cheap cases. The "magic measure" spread(s) is a candidate tournament magic monotone.

OPEN. Characterize the stabilizer (spread-0) score classes (which scores force H?); the growth rate of the
magic fraction; whether spread(s) has a closed form for the extremal unbalanced classes. -> HYP-2707, THM-027,
THM-555.
