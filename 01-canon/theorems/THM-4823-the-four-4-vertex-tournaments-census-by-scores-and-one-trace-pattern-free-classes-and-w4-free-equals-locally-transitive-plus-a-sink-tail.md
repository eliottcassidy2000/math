---
id: THM-4823
title: "The four 4-vertex tournaments inside a tournament: their counts are the score sums sum C(s,3), sum C(n-1-s,3), C(n,4) and one trace tr(A^4)/4, and the trace is not a function of the scores from order 5 on (counting strong 4-subtournaments costs one matrix multiplication); the pattern-free classes are K1/C3 stacks (S4-free, Narayana), orders <= 7 only (TT4-free, P7 last), and the W4-free class, which is exactly the locally transitive tournaments followed by a transitive tail of successive sinks, so W(n) = W(n-1) + A000016(n) - 1 = 1, 1, 2, 3, 6, 11, 20, 35, 64, 115 (not in OEIS)"
status: >
  PROVED (mac-mini-2026-10-27-slab; independent audit pending at filing): (A) the census identities (the trace
  identity is THM-118's tr(A^4) = 4 c_4, kind-pasteur); (B) the S4-free and TT4-free descriptions (classical inputs:
  Moon's pancyclicity, Erdős–Moser and Reid–Parker via THM-455/THM-4526); (C) the W4-free structure theorem and the
  recursion. The count of locally transitive classes is A000016 (Brouwer 1980, as attributed in THM-468). FINITE-EXACT:
  the census checks at orders 4-8 (pattern4.py) and the W4-free and locally transitive counts at orders 3-10 (w4free.c).
  No literature novelty claim for (C); OEIS has no entry for W(n) (searched 2026-10-27).
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/regular_slab_prescribed_black_core_eisenstein_catalyst_20261027.md
scripts:
  - 04-computation/experiments/regular_transitive_slab_20261027/pattern4.py (+ pattern4.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/w4free.c (+ w4free.out)
related:
  - 01-canon/theorems/THM-118-trace-cycle-identity.md
  - 01-canon/theorems/THM-584-complement-is-antipodal-map-level-parity-spectrum.md
  - 01-canon/theorems/THM-468-tournament-determinant-floor.md
  - 01-canon/theorems/THM-402-round-tournaments-are-2-dichromatic.md
  - 01-canon/theorems/THM-455-trans-subtournament-doubling-and-tower.md
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 01-canon/theorems/THM-1780-H-leaves-the-spectral-ladder-at-n6.md
---

# THM-4823 — the four 4-vertex tournaments: scores plus one trace, and the pattern-free classes

There are four tournaments on four vertices:

| type | scores | description |
|---|---|---|
| **TT4** | 0, 1, 2, 3 | transitive |
| **W4** | 1, 1, 1, 3 | a source over a 3-cycle |
| **L4** | 0, 2, 2, 2 | a 3-cycle over a sink |
| **S4** | 1, 1, 2, 2 | strong; exactly one Hamilton cycle |

THM-584 reads them as the Klein four-group of two bits, "source present" and "sink present". TT4 has both, W4 only a source, L4 only a sink, S4 neither.

## (A) The census is scores plus one trace (PROVED)

Let T have n vertices and scores s_v. Put P = Σ C(s_v, 3) and Q = Σ C(n − 1 − s_v, 3).
- **Sources and sinks.** A 4-set has at most one source, and v is the source of exactly C(s_v, 3) of them. So P counts the 4-sets of type TT4 or W4, and likewise Q counts those of type TT4 or L4.
- **The trace.** S4 is the only 4-tournament with a directed 4-cycle, and it has exactly one. So #S4 = c₄ = tr(A⁴)/4 (THM-118).

Hence

    #S4 = tr(A⁴)/4,   #TT4 = P + Q + #S4 − C(n,4),   #W4 = C(n,4) − Q − #S4,   #L4 = C(n,4) − P − #S4.

**The trace is genuinely needed.** #S4 is not a function of the score sequence from order 5 on. The number of score sequences on which #S4 takes more than one value is:

| n | 5 | 6 | 7 | 8 |
|---|---|---|---|---|
| sequences with non-constant #S4 | 2 of 9 | 11 of 22 | 41 of 59 | 136 of 167 |

So the 4-vertex census costs exactly one matrix product, O(n^ω), beyond the scores. This is one step up from the level of degree statistics. It matches the fact that H is a spectral invariant at orders 4 and 5 and stops being one at order 6 (THM-1780). `pattern4.py` checks every identity on every class of orders 4 to 8.

## (B) S4-free and TT4-free (PROVED from classical inputs)

- **S4-free.** A strong tournament on k ≥ 4 vertices has a 4-cycle (Moon: strong tournaments are vertex-pancyclic). So T is S4-free iff every strong component has at most 3 vertices. Equivalently, T is a transitive stack of K1 and C3 blocks.
  - Count: the compositions of n into parts 1 and 3, which is Narayana's cows A000930 (as in THM-4741 and THM-2454).
  - Values: 3, 4, 6, 9, 13 at n = 4, …, 8.
- **TT4-free.** Every tournament on 8 vertices contains TT4 (Erdős–Moser), and P₇ is the unique TT4-free class on 7 vertices (THM-455, THM-4526).
  - Counts: 3, 3, 1, 1, 0 at n = 4, …, 8.

## (C) W4-free = locally transitive, followed by a sink tail (PROVED)

T is **W4-free** iff every out-neighbourhood N⁺(x) is transitive. By the converse, L4-free means every in-neighbourhood is transitive, and the two classes have equal counts. **Locally transitive** means both; these are the round, vortex-free tournaments of THM-402 and THM-468.

**Theorem C.** A W4-free tournament without a sink is locally transitive. Consequently, T is W4-free iff T = L ⇒ TT_j: a locally transitive tournament L, with no sink unless |L| ≤ 1, followed by j ≥ 0 vertices each of which is a sink of everything before it. The number of W4-free classes is

    W(n) = W(n−1) + LT(n) − 1 = 1 + Σ_{j=2..n} (A000016(j) − 1),

and LT(n) − 1 counts the locally transitive classes without a sink; TT_n is the only one with a sink.

*Proof.*

**Setup.** Let T be W4-free with no sink, and suppose some N⁻(x) contains a 3-cycle C = (a → b → c → a).

**Step 1: C beats every out-neighbour of x.** Let y ∈ N⁺(x).
- **Some vertex of C beats y.** Otherwise C ⊆ N⁺(y), a 3-cycle inside an out-neighbourhood.
- **Then all three do.** By the cyclic symmetry of C, say a → y. Then N⁺(a) ⊇ {b, x, y} with b → x → y. Transitivity of N⁺(a) forces b → y. The same argument in N⁺(b) ⊇ {c, x, y} gives c → y.

**Step 2: C beats the whole out-closure.** Step 1 shows C ⊆ N⁻(y) for every y ∈ N⁺(x), so it applies again with y in place of x. By induction C beats every vertex of the out-closure R of x, the set of vertices reachable from x along arcs. Hence R ∩ C = ∅.

**Step 3: contradiction.** R is closed under out-neighbours, so a sink of T[R] would be a sink of T. Hence T[R] has no sink. But R ⊆ N⁺(a), and N⁺(a) is transitive, so T[R] is transitive and therefore has a sink. Contradiction.

**The structure statement.** A W4-free tournament with a sink z has T − z W4-free. Conversely, adding a sink keeps out-neighbourhoods transitive. Repeated removal of sinks ends at a W4-free tournament without a sink, which is locally transitive by the first part, or at a single vertex. Conversely, L ⇒ TT_j is W4-free whenever L is locally transitive: an out-neighbourhood N⁺_L(x) ∪ TT_j is transitive because N⁺_L(x) dominates TT_j.

**The count.** A locally transitive tournament with a sink z has N⁻(z) = V − z transitive, so it is TT_n. ∎

| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| W4-free classes, W(n) | 2 | 3 | 6 | 11 | 20 | 35 | 64 | 115 |
| locally transitive classes, LT(n) = A000016(n) | 2 | 2 | 4 | 6 | 10 | 16 | 30 | 52 |

Both rows come from the census (`w4free.c` over gentourng), and both agree with the formula. W(1) = W(2) = 1.

## Reading (typed)

- **IDENTITY.** The 4-vertex census sits exactly at the level of degree sums plus one trace, and three of the four pattern-free classes are counted by explicit recurrences or formulas.
- **NO TRANSFER.** A subexponential counter for forbidden local substructures over the C(n,2) arc bits, such as the external ETH preprint's 2^{o(n²)}, has nothing to add for 4-vertex patterns. TT4-free is empty from order 8, S4-free is a linear recurrence, and W4-free and L4-free are polynomial-time from A000016.
- **What would change that.** Richer pattern classes, with non-recursive counts, begin with 5-vertex patterns. That is where such a counter could first matter.
