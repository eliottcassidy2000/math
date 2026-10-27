---
id: HYP-9426
title: "After THM-4819 and THM-4820: (1) every slab class with R <= m-1 lies on a geodesic between the transitive and the regular class (which implies HYP-9423 (7)); (2) the prescribed black core at even balanced orders (RESOLVED, THM-4819 (E)); (3) Eisenstein merge debt at most floor(Q/3) + 1 everywhere (REFUTED by g(101)), replaced by the question whether the excess is bounded; (4) the seeds of the exceptional set, the points without a tight split (the first-drafted 'only if' is refuted by 59); (5) catalyst thresholds k*(n) ~ c n^2; (6) W4-free tournaments without a sink are locally transitive (PROVED, THM-4823); (7) at most two blue rotations per Hamilton cycle"
status: >
  OPEN (mac-mini-2026-10-27-slab, 2026-10-27; items (3)-(5) independently audited with THM-4820, corrections applied).
  (1) FINITE-EXACT at orders 5, 7, 9; sampled at 11 (1,000 of 1,000) and 13 (300 of 300). (2) RESOLVED the same day:
  PROVED for n >= 8 (THM-4819 (E); n = 8 by the audit's exhaustive Lemma B' check), FINITE-EXACT at 6, FALSE at 4 (S4). (3) as first drafted REFUTED by the
  audit: g(101) = floor(101^2/3) + 2; the excess is 0 or 1 on every norm below 10,201. (4) THM-4820 (C) PROVED for the
  points without a tight split; the first-drafted converse REFUTED by 59; the propagation form FINITE-EXACT for n <= 86
  (exact values, THM-4820 (D)). (5) OPEN; PROVED only k*(n) >= (n^2 - n - g(n))/3; values for n <= 15. (6) PROVED
  (THM-4823 (C)). (7) FINITE-EXACT on every Hamilton cycle at orders 7-11 (all regular classes) and on every circulant at 13 and 15
  (about 1.97 * 10^9 instances, THM-4819 audit); the weaker form, that not all m rotations are blue, is THM-4819 Lemma B.
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/regular_slab_prescribed_black_core_eisenstein_catalyst_20261027.md
related:
  - 01-canon/theorems/THM-4819-regular-tournaments-under-native-moves-the-transitive-slab-exact-distance-and-the-prescribed-black-core.md
  - 01-canon/theorems/THM-4820-eisenstein-merge-debt-floor-q-over-three-free-neutral-catalysts-and-the-eleven-forty-one-threshold.md
  - 05-knowledge/hypotheses/HYP-9423-dipole-endpoint-conjectures-two-connected-descent-and-the-prime-tower-predictions.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
---

# HYP-9426 — slab exactness, the even prescribed core, Eisenstein exceptions, catalyst thresholds, W4-free tournaments

1. **Slab exactness.** Let n = 2m + 1. Every tournament with L + R = (n² − 1)/8 and R ≤ m − 1 is exactly L moves from the transitive class.
   - **What it implies.** Every such class lies on a geodesic to the regular set. This includes every regular class, so it implies HYP-9423 (7).
   - **Evidence.** THM-4819 (B), at orders 5, 7 and 9 (complete censuses).
   - **Sampled at orders 11 and 13** (`slab_sample_search.c`). The samples are slab classes that are *not* built from TT: a random regular tournament followed by k random L-descending moves, so R = k. Each was searched for a chain of exactly L L-descending moves to TT, and every chain found was replayed.
     - Order 11: all 1,000 samples with R = 1, …, 4 are exactly L from TT.
     - Order 13: all 300 samples with R = 1, …, 5 are, one of them after raising the node budget to 300,000.
   - **Why the bound R ≤ m − 1.** A source or sink forces R ≥ m and blocks every L-descending move (THM-4819 (A4)). At order 7, strong λ = 1 failures already occur at R = m, so the bound cannot be relaxed to R ≤ m.
2. **Even balanced orders, prescribed vertex — RESOLVED the same day (THM-4819 (E)).** Let n = 2m ≥ 6 and let T have all scores in {m − 1, m}. For every vertex v there is a black word of length m − 1 into C(v, ρ), where ρ is the role the dipole law forces.
   - **Which role is forced.** After m − 1 moves:
     - a vertex of score m − 1 ends as a source if m is even and as a sink if m is odd;
     - a vertex of score m ends in the opposite role.

     Only that role is reachable in m − 1 moves.
   - **What is known.** Session75 §4 proves the unprescribed version.
   - **Evidence** (`prescribed_isolation.c` in even mode).
     - n = 6: all 5 classes, 30 instances.
     - n = 8: all 85 classes, 680 instances.
     - n = 10: all 13,333 classes, 133,330 instances.
   - **Order four fails.** At n = 4 only 2 of S4's 4 vertices succeed. The only Hamilton path of S4 = "100111" ending at its score-2 vertex 2 is 3 → 0 → 1 → 2, and it is blue. This is the even analogue of the T5 exception.
   - **Resolution.** The even black endpoint lemma (THM-4819 Lemma B′) halves the start set mod 2m − 1, which forces an arithmetic progression of step 2. The pair (0, 2) then breaks (b) for every m ≥ 7, and m = 5, 6 fall to exhaustive search.
   - **Orders 6 and 8.** The audit's exhaustive check of every Hamilton cycle settles n = 8 (m = 4). At n = 6 the lemma is false (3 all-blue cycles), and the census alone covers it.
   - **Order 4** is the exception.
3. **Eisenstein debt: REFUTED as first drafted; what remains open.** The draft claimed g(x) ≤ ⌊Q(x)/3⌋ + 1 for every
   x. The THM-4820 audit refutes this: g(101) = 3402 = ⌊101²/3⌋ + 2.
   - **Lower bound.** Every candidate first split with β ∈ {−4, …, 2} fails arithmetically.
   - **Upper bound.** An explicit tree of cost 3402.
   - **Below 101².** The excess is 0 or 1 on every norm below 10,201.
   - **More excess-2 points.** The same arithmetic gives an excess of at least 2 for 78 odd n ≤ 3000, beginning 101, 115,
     197, 235, 257, 289.

   Open, in its place:
   - **(3a)** Is the excess g(x) − ⌊Q(x)/3⌋ bounded? Each further unit of excess needs a wider band of norm equations to
     fail at once.
   - **(3b)** Does E = {x : excess ≥ 1} have density zero among all lattice points? It is infinite: Session75's
     EISENSTEIN_INVOICES proves 11^(1+820k) ∈ E for every k.
4. **The structure of E (corrected).** The draft conjectured that, for odd n with 3 ∤ n, n ∈ E iff n has only inert
   prime factors and (n² + 2)/3 is not of the form b² + 3t².
   - **"If"** is THM-4820 (C).
   - **"Only if" is REFUTED by n = 59.** The integer 59 is inert and (59² + 2)/3 = 3³·43, so 59 has a tight β = −1 split.
     Yet 59 is exceptional, and the value is exact (THM-4820 (D) window lemma): every tight split of 59 has an
     exceptional part.

   **Corrected statement.** Two facts are not conjectural:
   - **Propagation is a definition.** x ∈ E iff every tight split of x has a part in E. This is Session75's lifting
     criterion, and THM-4820 (D)'s window lemma keeps the recursion inside norms ≤ Q(x) + 1.
   - **The rational seeds are known.** The seeds are the points with no tight split at all. Among odd rational n with
     3 ∤ n they are exactly THM-4820 (C)'s class, because (C)'s two split criteria are equivalences (PROVED).

   What is open:
   - **(4a)** The seeds in general. 51 = 3·17 is a seed with norm ≡ 0 (mod 3).
   - **(4b)** A closed criterion for the exceptions that arise by propagation, such as 59 and the fifteen split-norm
     exceptions up to norm 2000.

   - **Odd exceptions to 86 (exact).** 11, 23, 25, 29, 41, 47, 51, 59, 83, 85.
   - **No tight split at all.** 11, 23, 25, 29, 41, 47, 83, 85 (by (C)), and also 51. Its norm 2601 is ≡ 0 (mod 3), and
     none of its four tight norm equations has a nontrivial solution.
   - **By propagation.** 59.
   - **Predictions that held.**
     - Every (C) prediction (11, 23, 25, 29, 41, 47, 83, 85) is exceptional.
     - 5, 17, 53, 55 and 71 are not.
     - The even predictions 92, 94 and 102 hold, and so do 88 and 100 (audit).
   - **A new odd exception.** 99 (audit).
   - **Even exceptions to 86.** 22, 44, 46, 50, 58, 82, each 2^j times an odd exception.
   - **Further predictions from (C).** 101, 107, 113, 115 are exceptional; the audit confirms 101 and 115, with excess 2.
5. **Catalyst thresholds.** g(n) = n²/3 + o(n²), and k*(n)/n² tends to a constant c.
   - **What is proved.** Only k*(n) ≥ (n² − n − g(n))/3 (THM-4820 (E)). The lower constant 2/9 would follow from the
     first claim, which is open, so it is not proved.
   - **Data.**
     - k*(n) for n = 1, …, 15: 0, 1, 2, 4, 7, 10, 14, 19, 24, 31, 36, 44, 52, 61, 70.
     - k*(15)/15² = 0.311.
   - **Test.** k*(n) for n ≤ 20 with K = 6n.
6. **W4-free tournaments — PROVED the same day (THM-4823 (C)).** Every tournament without a sink in which every
   out-neighbourhood is transitive is locally transitive. So W(n) = W(n−1) + LT(n) − 1, with LT agreeing with A000016
   at n = 3, …, 10.
7. **At most two blue rotations.** Let T be regular of order at least 7, let v be a vertex and let C be a Hamilton cycle of T − v. Then at most two of the m rotations of C closed at v are blue.
   - **What is proved.** Fewer than m (THM-4819 Lemma B).
   - **Evidence.** The THM-4819 audit (`audit_slab/lemmaB_all_cycles.c`) checked every Hamilton cycle of T − v for every vertex: all regular classes at orders 7, 9, 11, and every circulant at orders 13 and 15, about 1.97·10⁹ instances. The maximum is 2 at every order, including 13 (4,992 circulant cycles). The producer's random sample saw at most 1 at order 13.
