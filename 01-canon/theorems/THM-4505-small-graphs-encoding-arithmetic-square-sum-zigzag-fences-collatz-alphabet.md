---
id: THM-4505
title: "Small graphs that encode arithmetic: the square-sum problem first becomes possible at 15, with a unique path from 8 to 9 forced by a degree law. This is the k = 4 case of the Zigzag Threshold Theorem: for targets j^2 + k(4-k), the first chain appears at n = 4k-1, uniquely, from 2k to 2k+1. Friedman's fences satisfy an exact Euler field count and an isoperimetric bound, and complete grids exist iff 2n+1 is composite (Sundaram's sieve). The powers-of-2-or-3 sum graph is never Hamiltonian as a cycle, and its admissible windows are cut out by 3^a and the gaps |2^p - 3^a|"
status: >
  PROVED (elementary) + INDEPENDENTLY AUDITED (Lemma 1, Theorems 1-5,
  Propositions F1-F5, Theorem F2, Theorems C1-C2); FINITE-EXACT (tables and
  the zigzag converse over -300 <= c <= 100); CITED (OEIS A071983,
  A071984, A090461 (Gerbicz), A304120; Friedman's Maximum Fence Area table,
  whose values are best known, not proven optima); ANALOGY and NUMEROLOGY
  typed. Collatz is not addressed.
  Q_n is the square-sum graph on 1..n (x ~ y iff x+y is a square).
  (1) Degree law: deg_n(x) = floor(sqrt(x+n)) - floor(sqrt x) - [x = 2s^2].
  The vertices of degree <= 1, for every n, are given by an explicit
  finite table (x <= 18). Q_n has minimum degree >= 2 iff n >= 31.
  (2) Q_n has no Hamiltonian path for 2 <= n <= 14. Q_15 has exactly one,
  8,1,15,10,6,3,13,12,4,5,11,14,2,7,9. Its ends are forced: 8 (since
  8 + 8 = 16, its only partner is 1) and 9 (only partner 7) are the only
  leaves. Q_15 has exactly 15 edges, and the one unused edge is {1,3}.
  (3) Leaf Lemma (any target set): if at most two vertices have degree
  <= 1, each is a target, a half-target, or the successor of a target
  leaf.
  (4) Zigzag Threshold Theorem. Let k >= 2 and targets
  S_c = {j^2 + c}, c = k(4-k). The sum graph on 1..n has no Hamiltonian
  path for 2 <= n <= 4k-2, and exactly one at n = 4k-1: the zigzag Z_k
  from 2k to 2k+1, which uses the targets 2k+1, 4k, 6k+1. The squares are
  the case k = 4 (c = 0): n = 15, ends 8 and 9. FINITE-EXACT converse: for
  -300 <= c <= 100, a unique first chain whose ends are a half-target and
  the next integer occurs iff c = k(4-k).
  (5) Typing of the owner's "8 at one end, 9 at the other". For squares
  the zigzag needs 2k+1, 4k, 6k+1 all square. That is the simultaneous
  Pell system t^2 - 2s^2 = 1, u^2 - 6s^2 = 1, whose solution s = 2 gives
  3^2 - 2*2^2 = 1 and 5^2 - 6*2^2 = 1. This is the LAW. Reading 9 - 8 as
  the Catalan/Gersonides pair 3^2 - 2^3 (the free Collatz cycle
  -5,-7,-10 of THM-4484) is NUMEROLOGY for square sums, since the leaf
  property of 8 uses 2*8 = 4^2, never 8 = 2^3. The pair (8,9) is both
  Pell and Catalan only because t^2 - 2^m = 1 forces (3,3). The Catalan
  pair is a genuine mechanism where 8 and 9 are both targets (perfect-power
  sums; powers-of-2-or-3 sums): consecutive targets give zigzag chains.
  (6) Fences. n unit segments; two fences meet in at most one point, an
  end of at least one; every end lies on another fence; every field has
  area <= 1. With T-points (a fence passes through) and V0-points (only
  ends meet):
  - #fields = n + (#components) - V0, so A(n) <= n - 2 (F1).
  - A(n) <= ((sqrt(1 + 4n/sqrt(pi)) - 1)/2)^2 < n/sqrt(pi) (F2).
  - lambda = lim A(n)/n = sup A(n)/n exists and lies in
    [1/2, 1/sqrt(pi)] (F3).
  - A complete i x j grid uses 2ij + i + j fences, so it exists with n
    fences iff 2n+1 is an odd composite (Sundaram's sieve, F4). The records
    equal the best grid exactly at n = 4, 7, 12 (2n+1 = 9, 15, 25).
  - A(2k + ceil(2 sqrt k)) >= k (F5).
  - A(3) = sqrt(3)/4 is optimal.
  No record was improved.
  (7) The Collatz-alphabet sum graph C_n (x ~ y iff x+y is a power of 2
  or of 3). The largest power of 2 that is <= n has degree <= 1, so C_n
  is never Hamiltonian as a cycle (C1). Window Theorem (C2): the n with no
  isolated vertex and at most two leaves are exactly [2,3] and one window
  W_a for each a >= 2. With B = 3^(a-1), W_a is
  [max(3^a - 2^p, 2^(p+1) - 3^(a-1)), 3^a - 1] if (B, 3B) contains two
  powers of 2 (2^p < 2^(p+1)), and [2*3^(a-1), 3^a] if it contains one.
  Which case occurs follows the Beatty word of log_2 3. The left ends
  are the gaps |2^q - 3^b|, the same numbers as Collatz cycle denominators
  (ANALOGY: no map). Hamiltonian paths are tabulated for n <= 2200.
  OPEN:
  - whether lambda = 1/2 (wave 19 fence lane);
  - [CLOSED, see UPDATE below: s = 2 is the unique positive solution of
    the simultaneous Pell system (CITED, Anglin 1996)];
  - which n in each window W_a admit a Hamiltonian path of C_n;
  - the zigzag converse for all c.
source: collatz-procgen-20260922 session, smallgraph lane (2026-09-26), answering the owner's prompts on Friedman's fence problem, small graphs that encode arithmetic, and the square-sum problem at 15 with 8 at one end and 9 at the other; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (Gersonides/Catalan; used only for the typing in (5))
related:
  - 05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md
  - 05-knowledge/results/crossroads_family_20260926_squares.md (independent recovery of the unique n = 15 path: Q_15 is a five-cycle with two pendant paths at adjacent cycle vertices, and deleting their common cycle edge {1,3} leaves the path; first four-cycle and first two-edge switch at 46)
  - 05-knowledge/results/procgen_kuratowski_20260925_tait_kempe_triples.md (Q_n is planar iff n <= 24; the Kuratowski/Euler counting that F1 turns into an identity)
note: 05-knowledge/results/procgen_smallgraph_20260926_small_graphs_arithmetic.md
scripts: 04-computation/experiments/procgen_smallgraph_20260926_{lib,t1,t2,t3,fence6,run}.py and procgen_smallgraph_20260926_ham.c
script_audit: 04-computation/experiments/procgen_smallgraph_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_smallgraph_20260926.out
output_sha256: 7ef1341de08d565b5609d3dda75a21a2e3f7ee7f476191e344a7247b00350357
output_audit: 05-knowledge/results/procgen_smallgraph_20260926_orchestrator_check.out
output_audit_sha256: 0443df9f7fd33e3cd7218aad39291028a4465736c0bcd959bc81d9bbfa1ee65f
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of:
  - Lemma 1 and Theorems 1, 2, 4, 5 (including the degree bookkeeping at
    n = 4k-1: leaves {2k, 2k+1}, deg 1 = deg 3 = 3, all others 2, exactly
    n edges, so the omitted edge is {1,3} and the path is unique);
  - Propositions F1-F5 and Theorem F2;
  - Theorems C1 and C2.
  Independent code (procgen_smallgraph_20260926_orchestrator_check.py,
  written from the note's statements without reading the lane's scripts)
  confirms:
  - the leaf table and "min degree >= 2 iff n >= 31" for n <= 300;
  - no square-sum path for n <= 14, and a unique one at 15 from 8 to 9;
  - the Zigzag Threshold Theorem for k = 2..9, by exhaustive Hamiltonian
    path counting;
  - Sundaram's grid criterion for n <= 2000. A first version of this check
    failed at n = 316 because of the orchestrator's own range bug (grids
    with j < 100 only). Fixed, it passes.
  The lane's full pipeline was re-run (wall 750 s, peak RSS 330 MB,
  119 checks). Its output is identical to the committed .out up to timing
  and memory lines.
  Scope notes:
  - Friedman's table values are best-known lower bounds; nothing here
    proves any of them optimal except n = 3.
  - The zigzag converse is finite (-300 <= c <= 100).
  - The window left ends coincide with Collatz cycle denominators as
    numbers only.
---

# THM-4505 — small graphs that encode arithmetic: square sums, zigzags, fences, and the Collatz alphabet

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_smallgraph_20260926_small_graphs_arithmetic](../../05-knowledge/results/procgen_smallgraph_20260926_small_graphs_arithmetic.md).

**UPDATE 2026-09-26 (orchestrator; CITED).** The simultaneous Pell system `t^2 - 2s^2 = 1`, `u^2 - 6s^2 = 1` has the unique positive solution `s = 2`.
- **Source.** W. S. Anglin, "Simultaneous Pell equations", Math. Comp. 65 (1996) 355–359, proved that `x^2 - a z^2 = 1`, `y^2 - b z^2 = 1` (with `a != b`) has at most one positive solution whenever `max{a,b} <= 200`. This is quoted from M. A. Bennett, "On the number of solutions of simultaneous Pell equations" (J. reine angew. Math. 1998). Bennett's text was read 2026-09-26 from the author's page; Anglin's table itself was not read.
- **Consequence for square sums.** With `(a,b) = (2,6)` and the solution `(3,5,2)`, the three targets `2k+1, 4k, 6k+1` are all squares only for `k = 4`, i.e. `9, 16, 25`. So the zigzag of the square-sum problem is unique, and the owner's "8 at one end, 9 at the other" happens for squares exactly once.
- **Equivalent classical problem** (our deduction). Writing `s = 2m`, the system says that `m^2` is triangular (`8m^2 + 1 = t^2`) and generalized pentagonal (`24m^2 + 1 = u^2`). The pair `(a,b) = (8,24)` is also in Anglin's range, so `m = 1` is the only solution: 1 is the only square that is both triangular and generalized pentagonal. Wikipedia's "Pentagonal number" article (read 2026-09-26) says that no formal proof of the pentagonal square triangular case "has yet appeared in print", and credits unpublished work applying Anglin's paper. The deduction here rests on Bennett's statement of Anglin's theorem.

## 1. Why the square-sum problem first becomes possible at 15, from 8 to 9

In the square-sum graph, a vertex's degree is the number of squares in `(x, x+n]`, minus one if `2x` is itself a square. Every vertex of degree at most 1, for every `n`, lies in a small explicit table.

For `n ≤ 14` there is always an isolated vertex or three leaves, so no arrangement exists.

At `n = 15`:
- `8` (its double `16` is a square, so its only partner is `1`) and `9` (its only partner is `7`) are the only leaves, so they must be the ends.
- The graph has exactly one edge more than a path needs. That edge is `{1,3}`, and removing it leaves a single path.

The square case is the member `k = 4` of a law. For targets `j^2 + k(4-k)`, the first possible `n` is always `4k - 1`, with exactly one chain, a zigzag from `2k` to `2k+1`. For squares, the zigzag needs the three squares `9, 16, 25`. That is the first solution of the simultaneous Pell system `t^2 - 2s^2 = 1`, `u^2 - 6s^2 = 1`, so `3^2 - 2·2^2 = 1` is the law behind "8 at one end and 9 at the other".

Reading `9 - 8` as `3^2 - 2^3`, the identity behind the free Collatz cycle `-5, -7, -10`, is a coincidence of the pair `(8, 9)` with no mechanism connecting the problems. It becomes a real mechanism only in sum problems where 8 and 9 are both targets.

## 2. Fences as plane graphs

A fence configuration is a plane graph whose junctions are its vertices and whose bounded faces are the fields.
- **Field count.** Euler's formula counts fields exactly: `#fields = n + c - V0`.
- **Area bound.** Isoperimetry bounds the area: `A(n) < n/sqrt(pi)`.
- **Limit.** `A(n)/n` converges to a constant `λ` in `[1/2, 1/sqrt(pi)]`.
- **Grids.** Complete grids exist iff `2n+1` is composite: Sundaram's sieve reads fence counts as factorizations. Grids are record-optimal only at `n = 4, 7, 12`.
- **Beyond grids.** Irregular configurations win. The comparison with free versus sporadic cycles is an ANALOGY.

## 3. The Collatz alphabet

`C_n` joins `x, y` when `x + y` is a power of 2 or of 3. It is a small graph built from the two Collatz primes.
- It never closes into a Hamiltonian cycle.
- Its windows of admissible `n` are cut out by the powers of 3 and the gaps `|2^p - 3^a|`, in the order of the Beatty word of `log_2 3`.
