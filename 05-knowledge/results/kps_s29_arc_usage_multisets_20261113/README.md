# Arcs replaced by their path counts: what the multiset remembers

kind-pasteur-2026-11-13-S29 · [THM-5154](../../../01-canon/theorems/THM-5154-the-arc-usage-multiset-is-complete-up-to-run-slides-core-flips-and-converse-to-order-nine.md) · [HYP-9680](../../hypotheses/HYP-9680-the-arc-usage-multiset-determines-strong-tournaments-up-to-converse-and-the-slot-polynomial.md) · **Chan–Pak remains open.**

The owner asked what remains, and which arcs are left behind, when each arc of a tournament class is replaced by
the number of Hamiltonian paths through it. The resulting multisets at order 4 are:

| class | H | arc-usage multiset | arcs left behind |
|---|---:|---|---:|
| TT4 | 1 | {1, 1, 1} | 3 |
| C3 => v and v => C3 (one class up to converse) | 3 | {2, 2, 2, 1, 1, 1} | 0 |
| the strong 4-tournament | 5 | {4, 3, 3, 3, 1, 1} | 0 |

Session 319 works the same prompt in parallel. It has the labelled shadow, the exact strong-component product of
the weights, dark (nontraceable) arcs with the Busch–Jacobson–Reid classification credited, and the reusable
two-vertex gate. Session 318 has the reservoir results the owner mentions. This packet studies the **unlabelled**
multiset, the owner's "set", which forgets which arc carries which count.

## What the multiset keeps

**Left behind.** The arcs left behind are those between non-consecutive strong components, plus the dark arcs
inside strong components. Dark arcs first occur at order 5. Session 319 has both facts.

**What the multiset depends on (THM-5154 (A), proved).** By the component product, the multiset depends only on:
- the strong cores, in order;
- which singleton runs between them are empty;
- the total number of singletons;
- for each core with nonempty runs on both sides, only the union of its start and end weights.

So three moves leave both the multiset and the slot polynomial Q unchanged:
- (i) sliding singletons between nonempty runs, keeping which runs are empty;
- (ii) replacing a run-flanked core by its converse;
- (iii) global converse.

**Nothing else is forgotten through order 9 (THM-5154 (B), FINITE-EXACT, every class).**

| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|
| merged classes | 2 | 3 | 10 | 34 | 272 | 3,528 | 97,144 |
| distinct multisets | 2 | 3 | 10 | 34 | 271 | 3,526 | 97,121 |

- Every one of the 1 + 2 + 22 coincidences is a move.
- The first is TT1 => C3 => TT3 = TT2 => C3 => TT2 at order 7, a slide onto a self-converse class.
- The order-9 list contains 11 single-core slides, one two-core slide (TT2 => C3 => TT1 => C3 against
  TT1 => C3 => TT2 => C3) and 10 core flips TT2 => S => TT1 against TT2 => S^op => TT1.

**Consequences, through order 9.**
- The multiset determines the full slot polynomial of every tournament.
- It determines every strong tournament up to converse.
- It is far finer than Q: of the 21,217 order-9 slot-polynomial classes holding at least two merged classes,
  21,199 are fully separated (THM-5153).

**Beyond order 9 (pre-registered P1, CONFIRMED).** In the uniquely Hamiltonian strong slice, the multisets of
orders 10, 11, 12 are pairwise distinct up to converse: 195, 517 and 1,309 merged classes, of which 13, 47 and 34
are self-converse.

**What fails at higher order (credited).** Session 319's labelled-shadow collisions at orders 22 and 42 flip a
vertex-transitive core whose start and end weights are flat. Those flips preserve the seam multiset even next to an
empty run. So the general move set needs seam-preserving flips, which is HYP-9680 (C).

## The 3-and-2 thread the owner asked to merge (synthesis, credited)

- **Collatz cycles.** On the positive integers the only known cycle is 1 -> 4 -> 2. On the negative integers the
  three known cycles are those of −1, −5 and −17. Count the odd steps k and the halvings L in each cycle:
  - the positive cycle has (k, L) = (1, 2);
  - the negative cycles have (k, L) = (1, 1), (2, 3) and (7, 11).

  Each cycle element is a quotient by 2^L − 3^k, which here is 1, −1, −1 and −139. The ratios L/k = 2, 1, 3/2,
  11/7 bracket log_2 3 = 1.58496...:
  - 2 lies above log_2 3, so the cycle is positive;
  - the others lie below, so the cycles are negative;
  - 11/7 is the semiconvergent between 3/2 and 8/5.

  So every known cycle is a near-coincidence 2^L ≈ 3^k with a tiny denominator. The repository's THM-4591 (integer
  Terras cycle periods {1, 2, 3, 11}) and THM-4849 (circulants as rational Terras cycles with denominator
  2^q − 3^k) hold this. The same small values |2^x − 3^y| govern the Pillai chain law of THM-5109 (Ellison's
  bound).
- **Reservoirs.** Session 318 shows two strong reservoirs realizing 9·2^(k−2) + 9 at every k >= 4, and reservoirs
  5, 3, 1 realizing 21·2^k + 15, +17, +19. These are dyadic rays whose intercepts are small 3-smooth multiples:
  9 = 3^2, and 21 = 3·7 contains the forbidden 7.

  In the arc-usage language, a reservoir change moves weight between seam arcs and inner run arcs. That is exactly
  the data the multiset keeps only up to moves (i) and (ii).

## Replay

In a scratch directory with classgen's cls_n.txt (kps_s26 scripts/README), run `scripts/arcms.c cls_n.txt n`, then
`python scripts/arcanal.py` (census) and `python scripts/collision_moves.py` (moves). For P1, run
`scripts/arcms_rows < c1_n.txt` on the C=1 lists, then `python scripts/c1arccheck.py`.
