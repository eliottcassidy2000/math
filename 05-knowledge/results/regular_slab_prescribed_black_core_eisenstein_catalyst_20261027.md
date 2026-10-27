# The transitive slab, the prescribed black core, the Eisenstein catalyst, and the four 4-vertex tournaments

2026-10-27, mac-mini-2026-10-27-slab. **Status:** THM-4819, THM-4820, THM-4823 (PROVED parts plus FINITE-EXACT
censuses) and HYP-9426 (OPEN), plus a ledger addendum and a typed reflection. All of it has been independently audited
(three audits; reports under `04-computation/experiments/regular_transitive_slab_20261027/audit_*/`). Every proof was
confirmed and every number reproduced; the corrections are applied, and MISTAKE-687 records them.

## Owner's directive

The core text of the directive, as given:
- "merge in prior work on 11 and 15 and think about the fact every regular tournament of order 11 is exactly 15 moves
  from the transitive class";
- the ETH-false counting paragraph from klein's reflection, with its two next steps, "add the branch-count coordinate
  to the classification ledger, and place every atlas row at the strongest invariant level its partial results reach";
- "Exponential Time Hypothesis itself quantifies over algorithms and sits at roughly Σ⁰₄, so think of how its recent
  disproof relates to extensive prior attempts at connecting the tournament on 4 vertices with other concepts. and
  connect recent progress on galactic algorithm speedups matrix multiplication etc with our repo concepts and find
  extensions through synthesis. think how the first few handfuls of tournaments [...] and the micro/macrocosm duality
  relate to the below source of inspiration", followed by a pasted text on the arithmetic hierarchy at levels 4–10,
  the analytical hierarchy and set theory;
- "consider how [every] prescribed vertex of a regular tournament on 2m+1 ≥ 7 vertices can be isolated from a strong
  core in exactly m black moves, with both native source/sink roles handled. The Eisenstein work also gives a sharp
  catalyst threshold: the minimum merge cost for charges (11+k,k,k) is 41 exactly when k ≥ 36."

The pasted inspiration text was treated as data.

## What is new

1. **The slab (THM-4819 (A), PROVED).** L + R ≥ (n²−1)/8 for every tournament of odd order.
   - Equality holds iff every sorted score lies between its index and the middle.
   - Every TT–regular geodesic stays in this slab.
   - The slab classes on geodesics are exactly those reachable from TT by R-descending moves.

   The censuses at orders 5, 7 and 9 put every slab class with R ≤ m − 1 on a geodesic (HYP-9426 (1)). Samples of slab
   classes built from regular tournaments, rather than from TT, are all exactly L from TT: 1,000 of 1,000 at order 11 and
   300 of 300 at order 13 (`slab_sample_search.c`).
   - "Order 11 is exactly 15" is the realized score-transport distance 1 + 2 + 3 + 4 + 5.
   - Moves only rearrange score units monotonically.
   - The hard part is realizing each descending pair by a Hamilton path while avoiding sources and sinks.
2. **Exact distance (n²−1)/8 (THM-4819 (C), FINITE-EXACT).** Every word below was found by search and replayed.
   - **Order 11:** every regular class (THM-4815 B8').
   - **Order 13:** all 64 circulants and 1,000 random regular tournaments.
   - **Order 15:** all 128 circulants and 200 random regular tournaments, plus R5[C3], C3[R5] and D(P₇).
   - **Order 19:** P₁₉, at exactly 45 moves.

   Prior work merged: Session74's P₁₁ certificate, THM-4810 (C) (lower bound and the order-15 forest-frame threshold),
   HYP-9420 (4), THM-4812 and THM-4557 (the special order-15 classes).
3. **The prescribed black core (THM-4819 (D), PROVED).**
   - **Lemma B.** For a regular T of order 2m+1 ≥ 7, a vertex v and *any* Hamilton cycle of T − v, some rotation of the
     cycle closed at v is black. If all m rotations were blue, the anti-automorphisms σ_t would force
     N⁺(c_t) = σ_t(N⁻(v)), which says that U contains no two elements 2 apart and that exactly one of
     2t − t′ + 1 and 2t′ − t + 1 lies in U. A counting argument mod 4 refutes this for every m ≥ 3; at m = 2 it is
     consistent, matching the pure-blue T5.
   - **Consequence.** With Session75's audited construction from its second move on, and the converse for the other
     role, every prescribed vertex and both native roles are reached in exactly m black moves. This resolves the
     strengthening HYP-9422 (1) left OPEN. An independent search confirms it on every (class, vertex, role) at orders 7,
     9 and 11: 42, 270 and 26,906 instances.
   - **Even balanced orders (THM-4819 (E)).** At n = 2m the dipole law forces each vertex's role: score m − 1 versus
     m, combined with the parity of m.
     - **The even lemma.** Lemma B′ halves the start set mod 2m − 1, which forces a step-2 arithmetic progression. The
       pair (0, 2) then breaks (b).
     - **The theorem.** With Session75 §4, every vertex reaches its forced role in exactly m − 1 black moves, for every
       n ≥ 10.
     - **Evidence and exceptions.** Censuses confirm n = 6, 8, 10 (30, 680 and 133,330 instances). n = 4 fails at S4,
       whose only Hamilton path ending at a score-2 vertex is blue.
4. **The Eisenstein catalyst (THM-4820; independently audited, every number PASS).** At filing the owner's statement
   had no source in the repository; the search covered every text file and the dyadic, paid-ROOT, renewal and necklace
   branches. Session75's `session75_final_mathematics_20261027/EISENSTEIN_INVOICES.md`, committed minutes later, derives
   the same result independently under the same reading, and is evidently the owner's source. Under that reading the
   numbers hold exactly:
   - **The reading.** Unit charges 1, ω, ω²; free binary merges; each merge costs β(x,y)₊, Session75's dyadic pair
     interaction.
   - **PROVED.**
     - cost = ⌊Q(root)/3⌋ + Σ slack, so g(x) ≥ ⌊Q/3⌋.
     - Neutral triples are free catalysts.
     - 11 has no tight split, because 11 is inert and (11²+2)/3 = 41 is not of the form b² + 3t². So g(11) = 41.
   - **FINITE-EXACT.**
     - A cost-41 tree on (47,36,36), checked independently.
     - f(46,35,35) = 42, by exact dynamic programming over all count vectors.

   So **f(11+k,k,k) = 41 exactly when k ≥ 36**. Further results:
   - **The lattice table.** g is ⌊Q/3⌋ or ⌊Q/3⌋ + 1 on all norms ≤ 2000, with twenty exceptional norms.
   - **Thresholds.** k*(n) for n ≤ 14 is 0, 1, 2, 4, 7, 10, 14, 19, 24, 31, 36, 44, 52, 61.
   - **The larger window (n ≤ 86).**
     - **Exact values.** By the audit's window lemma (a cost-C tree has every node at norm ≤ 3C + 1), these values are
       exact.
     - **Predictions.** Every (C) prediction held: 83 and 85 are exceptional, and 53, 55, 71 are not.
     - **New exceptions.** 51 is a seed of norm ≡ 0, with no tight split at all. 59 is an exception by propagation, and
       it refutes the first-drafted converse.
   - **Beyond (audit).** g(101) = ⌊101²/3⌋ + 2, which refutes HYP-9426 (3) as drafted. The excess is 0 or 1 below norm
     10,201.
5. **The four 4-vertex tournaments (THM-4823).**
   - **The census.** It is Σ C(s,3), Σ C(n−1−s,3), C(n,4) and one trace, tr(A⁴)/4 (THM-118). The trace is not a
     function of the scores from order 5 on, so the census costs exactly one matrix product.
   - **S4-free and TT4-free.** The S4-free classes are K1/C3 stacks (Narayana). TT4-free classes exist only at
     orders ≤ 7.
   - **W4-free (PROVED).** A W4-free tournament without a sink is locally transitive, so
     W(n) = W(n−1) + A000016(n) − 1 = 1, 1, 2, 3, 6, 11, 20, 35, 64, 115. This sequence is not in OEIS.
6. **Ledger addendum.** PROBLEM_CLASSIFICATION §5 adds, per klein's reflection §6:
   - the branch-count coordinate, the three invariant levels and census horizons, with every row placed;
   - the ETH level remark (corrected by the audit: ETH proper is Σ⁰₃, the 2^{o(n)} form is Π⁰₄, the preprint's claim is Σ⁰₂);
   - one new tournament row, whose strongest invariant is the I2 slab.
7. **Reflection.** [eth-sigma-four-the-four-vertex-tournament-galactic-algorithms-and-the-micro-macrocosm-macmini-20261027](../../07-reflections/eth-sigma-four-the-four-vertex-tournament-galactic-algorithms-and-the-micro-macrocosm-macmini-20261027.md).

## Reproduction

All code is in `04-computation/experiments/regular_transitive_slab_20261027/`; binaries are gitignored. Build each C
program with `cc -O3 -march=native -I/opt/homebrew/opt/nauty/include file.c /opt/homebrew/opt/nauty/lib/libnauty.a`
(or `libnauty1.a` for the C++ files).

- **Slab.** `python3 slab_census.py n` for n = 5, 7 (graph files of the dipole packet); `slab_closure` (order 9);
  `forward_closure n K`.
- **Exact distances.** `regular_geodesic_search n mode count seed budget cands`. Then
  `paley_geodesic_search n {2|3} 1 seed budget cands` (mode 2 is Paley; mode 3 takes `ROWS` in hex) and
  `rerun_undecided.sh`.
- **Black core.**
  - `gentourng -q -d<m>D<m> n | ./prescribed_isolation n 400000`;
  - `python3 start_set_conditions.py 13`;
  - `python3 gap_parity.py 2 40`;
  - `python3 black_endpoint_check.py 7`.
- **Catalyst** (in `catalyst/`).
  - `./merge_debt 11 60`; `./merge_debt 11 36 tree`;
  - `python3 verify_tree.py tree_47_36_36.txt 47,36,36`;
  - `./lattice_debt 2500 50 dump`; `./lattice_debt 7500 86`;
  - `python3 tight_splits.py lattice_W2500.txt`; `./thresholds.sh`.
- **Four vertices.** `python3 pattern4.py 8`; `gentourng -q n | ./w4free n`.

## Not done

- **Larger thresholds.** k*(15) = 70 needed a box with K = 90; n ≥ 16 needs larger boxes.
- **Larger Paley tournaments.** P₂₃ and P₃₁ were not attempted: the subset-DP memory at n = 23 needs a per-start
  variant.
- **Five vertices.** No five-vertex census was attempted (the reflection, item 3).

## Audit outcome (2026-10-27)

- **THM-4819** (`audit_slab/REPORT_THM4819.md`). Every proof was PROVED as written. The audit added:
  - the strong form of Lemma B, checked on every Hamilton cycle: all regular classes at orders 7–11 and every
    circulant at orders 13 and 15, about 1.97·10⁹ instances;
  - the proof's construction executed and replayed;
  - 127 archived geodesic words.

  Corrections:
  - R5[C3] and C3[R5] are Z₁₅ circulants, not non-Cayley;
  - one side clause of a remark was valid only for the weak form, so m = 4 joins its list;
  - (E) now holds from n = 8, after an exhaustive Lemma B′ check;
  - wording fixes.
- **THM-4820** (`audit_catalyst/REPORT_THM4820.md`). Every number passed, including f(46,35,35) = 42, certified over
  5·10⁸ pairs. Corrections:
  - **Exact, not window-dependent.** The windowed values are exact, by the window lemma (proved).
  - **51** has no tight split at all.
  - **Overclaims removed.** The inertness reading is restricted to the (C) class, and the two-thirds asymptotics is
    marked conjecture.
  - **HYP-9426 (3)** is refuted by g(101) = ⌊101²/3⌋ + 2.
  - **Concurrent derivation.** Session75's EISENSTEIN_INVOICES derives 41 at k ≥ 36 independently, and is credited.
- **THM-4823, the ledger and the reflection** (`audit_four/REPORT_THM4823_LEDGER.md`).
  - THM-4823 passed, with wording fixes and added formulas: the labelled count W_lab and W(n) = A133670(n) − n.
  - Side finding: THM-402's verification count of 544 is the W4-free count; a correction note is added there.
  - In the ledger §5:
    - the ETH level is corrected to Σ⁰₃ (a note is added upstream in klein's reflection);
    - the horizons are taken from canon (Collatz 2⁷¹, cycles above 1.375·10¹¹ odd elements, JC(2) degrees 101 and 108,
      H ≥ disc through order 7);
    - the definitions and the Reading are corrected;
    - missing rows are added.
  - In the reflection:
    - the rigidity mechanism is trace duality, not Lipschitz;
    - "rank is not Lipschitz" is deleted;
    - the Klein-four sentence is fixed;
    - EM is typed as an ANALOGY, and its citations were web-checked.

