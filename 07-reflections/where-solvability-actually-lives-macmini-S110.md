# Where solvability actually lives — the ladder, the shelves, and the fiber monodromy

*mac-mini-2026-07-15-S110. Companion canon: THM-868 (geometric ladder), THM-869 (Morse
shelves); probes in ryser_fiber_shelves_macmini_S110. The owner asked how the +8 walk
relates to the Cayley–Dickson tower and the geometry that kills the quintic. The honest
answer has three layers, and the discipline is in keeping them apart.*

## 1. The walk is the solvable part — literally graded

THM-866's climb is a composition series in miniature: one conserved step (+8), abelian
all the way, terminating at the total order. THM-868 now shows the same shape on the
integer side: the run-filtration is a GEOMETRIC series in one substituted variable
(v = x²/(1−x)² for rows, u = x³/((1−x)²(1+x)) for diagonals), and 2^n and Fibonacci are
its telescoped limits. A tower is "solvable" in this sense when its layers are powers of
one thing — geometric tails, conserved steps, radicals. Moser's missing 32nd region is
the first coefficient of a geometric tail; the deficit thereafter is not an anomaly but
a *lawful shadow* (GF x^{2j+1}/((1−x)^{2j+1}(1−2x))). Where the tower is geometric,
nothing interesting can hide: every level is populated, every deficit is predictable.

## 2. The shelves are NOT the quintic — a seduction corrected

The descent obstruction (kps cont.22's false peaks) looked like Abel–Ruffini appearing at
n = 5. The n = 7 full-move census (16,800 stuck) and the n = 6 zero kill that reading:
stuckness is **odd-n floor-parity geometry** — regular floors exist only at odd n, and
the obstruction is the upset-saturated {0,±2} configuration (THM-869), with shelf a = 2
first opening at n = 9, not at any Galois-theoretic threshold. Lesson repeated from
MISTAKE-137, one level up: a threshold at n = 5 is only evidence about 5 if the n = 6
and n = 7 cases are computed in the SAME invariant language. Unit-step stuckness and
full-move stuckness diverge exactly where the seduction lived.

## 3. What actually changes at n = 5 — the fiber monodromy

The precise anchor: score fibers are connected under directed-3-cycle reversals (Ryser —
machine-confirmed on all 291 vector fibers at n = 5), so "loop in the fiber" is a real
monodromy notion. At n = 4 every reversal preserves α₁ mod 2 (0 flipping edges out of
32): digit 1 of H is fiber-constant, H is score-determined — the abelian shadow (scores)
sees everything. At n = 5 there are exactly 240 parity-flipping reversal edges, all
inside the relabelings of the (1,2,2,2,3) class: **digit 1 becomes monodromy-carried at
n = 5** — the first n where the fiber's move-group can smuggle an invariant past the
scores, and the same n where S_n first fails to be solvable. Local-law hunt: three
candidate one-triangle statistics refuted; one necessary condition survives (a flipping
reversal always has an ODD number of 5-cycles sharing exactly one arc with the reversed
triangle — the (1,0) cell of the census is empty). This is Arnold's shape of argument in
tournament clothing: obstructions that live on loops and refuse to be functions of local
data. What we do NOT have (and should not claim): a homomorphism from fiber loops onto
A₅, or any causal link from S₅'s simplicity to the 240. Logged as the sharpest open
question of the thread: is the flipping-edge set of the Ryser complex a coboundary for
some richer local functional, or genuinely a monodromy character? Compute the cycle
space of the Ryser graph on the (1,2,2,2,3) fiber and test whether parity is constant on
commutator loops — the literal Arnold test, and it is finite.

## 4. The tower taxonomy (the abstraction the owner asked to keep)

Three tower species now sit side by side in the repo:

| species | layers | tail behaviour | examples |
|---|---|---|---|
| geometric (solvable) | powers of one substitution | exact product GFs, lawful shadows | THM-868 ladders; binary/Zeckendorf numeration; THM-866's +8 chain |
| property-decay | each level loses one structure | decay point, then lawless | Cayley–Dickson R→C→H→O→S (order, comm., assoc., division); OCF 2-adic tower (digit 0 Rédei-constant, digit 1 dies at n = 7 / score-wise at n = 5); THM-865's locker family |
| periodicity | layers repeat mod k | the period is the law | Bott/Clifford mod 8; our mod-8 level lattice; the odd-support split's (1+i)^N period-8 law (A038505) — the one honest "8" |

The Cayley–Dickson connection is genus, not mechanism: CD and the OCF digit tower are
both property-decay towers (kps cont.22's dictionary), and the CD nodes n = 2^k + 1 do
sit on our level lattice (counts 2, 6, 31, 205 — but so does every n). The quintic
connection is also genus: Galois solvability is the statement that a number's tower is
geometric-species; Abel–Ruffini/Arnold is the statement that S₅'s tower is not, and the
tournament world crosses the same line at n = 5 in the fiber-monodromy sense of §3.
The mod-8 of the axis and the mod-8 of Bott are, on present evidence, homonyms — except
for the one computed bridge: the support filtration's odd/even split oscillates by
(1+i)^N, the same Gaussian-integer periodicity that drives Clifford's period 8. That is
a real shared mechanism (fourth roots of unity in a binomial transform), and exactly as
much as the evidence supports.

## 5. Beyond (the program, one line each)

- Arnold test on the Ryser complex (finite, n = 5): parity on commutator loops.
- The shelf census sequence (40, 16800, …; a = 2 counts at n = 9) — closed form.
- The u-substitution grammar: which other repo towers (Farey-14 depth, OCF layers) admit
  a substitution variable making them geometric? Where none exists, THAT is the theorem.
- OEIS: seven new sequences from one filtration (THM-868 table) — submit with b-files.
- The A362193 coincidence (d = 3 Moser rows = Grassmannian pattern-avoiders): find the
  bijection; Grassmannian permutations have exactly one descent — suspiciously "rank 2".
