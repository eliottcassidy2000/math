# Edge multiset dimension of hypercubes: edim_m(Q_6) = 15, and how edim_m(Q_d) grows

**Lane:** procgen_edim, 2026-10-01 (worktree collatz-procgen-20260922, machine mac-mini).
**Source problem:** J. Allikvere, *The edge multiset dimension of hypercubes*, arXiv:2608.09983v1
(2026-08-05), Section 10, Open Problems 1, 2 and 4.
**Runner:** `04-computation/experiments/procgen_edim_20261001_run.py`. Its stdout is in
`05-knowledge/results/procgen_edim_20261001.out`: 127 checks, ending with ALL CHECKS PASSED.
**Data:** `05-knowledge/results/procgen_edim_20261001_q6_resolving15_orbits.txt` lists the 229 extremal orbits.
**Promoted:** [THM-4525](../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md) (FINITE-EXACT + INDEPENDENTLY AUDITED, orchestrator 2026-10-01: a third exhaustive search with its own reduction confirms no resolving set of size <= 14 and the same 229 orbits at size 15; audit `procgen_edim_20261001_orchestrator_check.out`); the growth conjecture is [HYP-9169](../hypotheses/HYP-9169-edge-multiset-dimension-of-hypercubes-grows-like-exp-cube-root.md).

## Status header

| # | Claim | Label |
|---|-------|-------|
| 1 | **edim_m(Q_6) = 15.** This settles Open Problem 1. | FINITE-EXACT (two independent exhaustive searches, cross-checked) |
| 2 | Q_6 has exactly **229** Aut(Q_6)-orbits of edge-multiset resolving 15-sets. Every one has trivial stabilizer, so there are exactly 229 * 46080 = **10,552,320** resolving 15-sets. | FINITE-EXACT |
| 3 | For k = 14 the minimum defect (192 - #distinct histograms) is **1**. Exactly one orbit attains it. Exactly 16 orbits have defect 2. | FINITE-EXACT |
| 4 | Lemmas L1-L4: trivial stabilizer; antipodal reversal with the odd-size parity corollary; Walsh alternating-sum identity; edim_m(Q_6) >= 7 by counting. | PROVED |
| 5 | L5 (entropy bound): edim_m(Q_d) >= exp((c - o(1)) d^(1/3)) with c = (ln 2 / sqrt 2)^(2/3) = 0.6215. **So edim_m(Q_d) grows faster than any polynomial in d.** | PROVED |
| 6 | Upper bounds from explicit sets: edim_m(Q_7) <= **19**, Q_8 <= 26, Q_9 <= 38, Q_10 <= 48, Q_11 <= 65, Q_12 <= 76. The paper had 63, 115, 246, 492 for d = 7..10. | VERIFIED (exact check of explicit sets) |
| 7 | The random-flow forest lemma holds for every inclusion probability q. With sparse q the union bound gives edim_m(Q_d) <= M_d with ln M_d / d^(1/3) = 2.76-2.80 for d = 11..32 (for example Q_11 <= 511, Q_16 <= 1056, Q_32 <= 6638). Density 1/2 gives about 2^(d-1). | VERIFIED (double precision, 1e-6 relative safety margin per factor, no interval arithmetic) |
| 8 | No landmark set built from the tournament structure of the n = 5 tiling cube resolves: unions of isomorphism classes, of H(T)-levels, or of score classes. | FINITE-EXACT |
| 9 | ln edim_m(Q_d) = Theta(d^(1/3)). | OPEN (conjecture; lower half PROVED as L5; upper half supported by item 7 numerics and a heuristic) |

## 0. Definitions and cited facts (CITED: arXiv:2608.09983v1)

Let G be a graph, S a nonempty subset of V(G), and e = uv an edge. Put d(e,s) = min(d(u,s), d(v,s)).
The edge multiset representation of e is the multiset {d(e,s) : s in S}. Equivalently it is the
histogram H_e(r) = #{s in S : d(e,s) = r} for 0 <= r <= d-1 (paper Sec. 2). S is edge-multiset
resolving if no two edges have the same representation. edim_m(G) is the least such |S|, or infinity
if there is none.

Results cited from the paper:
* Lemma 3 (projection lemma): for e = {u, u + e_i}, d(e,w) = d_H(u', w'), where ' deletes coordinate i.
* Prop. 2: |{w : d(e,w) = r}| = 2 C(d-1, r), and complement symmetry.
* Thm. 6: edim_m(Q_d) is infinite exactly when 2 <= d <= 5.
* Prop. 9 / Table 1: certificates of sizes 15, 63, 115, 246, 492 for d = 6..10.
* Prop. 10: 6 <= edim_m(Q_6) <= 15.
* Lemma 11: the random-flow forest lemma.
* Sec. 10 remark: 96 annealing restarts at size 14 found nothing. This was heuristic only, and is now explained by Theorem 1.

The runner (S1) re-verifies the Q_6, Q_7 and Q_8 certificates of Table 1 with two independent checkers.

## 1. Theorem 1: edim_m(Q_6) = 15

**Upper bound.** The paper's Table 1 set (mask 0x02283022a042a00a) is resolving (VERIFIED). So are
representatives of 228 further orbits (Section 2).

**Lower bound.** No edge-multiset resolving set of size k <= 14 exists. For k <= 6 this also follows
from Prop. 10 or L4; the computation covers every 1 <= k <= 14 anyway.

### 1.0 Why naive class-splitting pruning is invalid here
Multiset representations do not refine as landmarks are added. Take S = {17, 8, 32} and s = 7. The
edges {0,2} and {1,3} have histograms (0,2,1,0,0,0) and (0,1,2,0,0,0) under S, but both have
(0,2,2,0,0,0) under S + s. Classes of equal representations can therefore merge, so a
"split each class into at most 6 parts" bound cannot be used to prune. Instead we use symmetry
normal forms and enumerate completely.

### 1.1 Reduction (PROVED)
Write V(Q_6) = V(Q_5) x {0,1}, with coordinate 5 as the last coordinate. Let
beta_i(S) = #{s : s_i = 0} - #{s : s_i = 1}. Edge-multiset resolvability is Aut(Q_6)-invariant,
because d(sigma e, sigma s) = d(e, s).

* **Normal form A.** Every S is Aut(Q_6)-equivalent to A x {0} u B x {1}, where:
  * A is the canonical (minimum-mask) representative of its Aut(Q_5)-orbit;
  * |A| = a >= b = |B| and a - b = max_i |beta_i(S)|;
  * |beta_i| <= a - b for i = 0..4.

  Proof: transpose a coordinate of maximal |beta_i| with coordinate 5. Flip coordinate 5 if beta_5 < 0.
  Then apply the h in Aut(Q_5), acting on coordinates 0..4, that makes A canonical. This h only permutes
  and negates beta_0..beta_4.
* **Normal form B.** The same construction with a coordinate of *minimal* |beta_i|. Then
  |beta_i| >= a - b for i = 0..4.

### 1.2 Computation (FINITE-EXACT)
* **Orbit representatives.** `procgen_edim_20261001_q5orbits.c` lists the Aut(Q_5)-orbit
  representatives of a-subsets of Q_5 for a = 0..15. The counts are 1, 1, 5, 10, 47, 131, 472, 1326,
  3779, 9013, 19963, 38073, 65664, 98804, 133576, 158658. A numpy check certifies the lists complete
  and irredundant: the representatives are pairwise distinct, each is the minimum of its orbit under
  all 3840 automorphisms, and each count equals the Burnside count.
* **Method A** (`procgen_edim_20261001_searchA.c`). For each k, each a >= k/2, each representative A
  and each b-subset B that passes the A-filter, it tests S.
  * Histogram key: sum_s 16^d(e,s) over levels 0..4. This is injective for k <= 15, because level 5 is
    determined by the others.
  * Duplicates are detected with a 1024-slot probing table, exiting at the first collision.
* **Method B** (`procgen_edim_20261001_searchB.c`). It uses the independent normal form B and
  independent code:
  * distances via the projection lemma;
  * edges ordered by direction;
  * DFS in decreasing vertex order;
  * a direct-mapped 2^20 stamp array instead of hashing.
* **Leaf counts.** For every (k, a) and both methods, an independent numpy dynamic program over the
  imbalance vector recomputes the number of leaves. The counts agree exactly.

| k | leaves A | leaves B | resolving leaves (A / B) |
|---|---------:|---------:|:---:|
| 1-6 | 15,358 | 60,205 | 0 / 0 |
| 7 | 70,347 | 234,240 | 0 / 0 |
| 8 | 375,809 | 1,808,073 | 0 / 0 |
| 9 | 1,953,411 | 4,767,589 | 0 / 0 |
| 10 | 9,222,226 | 29,955,834 | 0 / 0 |
| 11 | 39,735,306 | 96,667,005 | 0 / 0 |
| 12 | 161,182,358 | 491,865,822 | 0 / 0 |
| 13 | 602,815,147 | 1,234,172,511 | 0 / 0 |
| 14 | 2,075,614,912 | 5,359,115,934 | 0 / 0 |
| 15 | 6,609,712,486 | 13,146,176,599 | 314 / 585 |

Running time on one Apple M2 core: about 50 ns per leaf.

### 1.3 Cross-checks (all re-run by the runner)
1. **Positive control at k = 15.** A gives 314 resolving leaves and B gives 585. Canonicalized under
   Aut(Q_6) (46080 elements, `procgen_edim_20261001_canon6.c`), both reduce to the **same 229 orbits**.
   The paper's set is one of them. All 229 representatives were re-verified by a pure-Python checker
   that shares no code with the C searches.
2. **Near-miss control at k = 14.** With defect output switched on, A and B find the same 17 orbits of
   defect <= 2: 1 orbit of defect 1 and 16 of defect 2. Python re-verifies the minimum defect.
3. **WLOG test.** 100 random 14- and 15-sets were mapped to both normal forms. Each landed in the
   searched domain and in the same Aut(Q_6)-orbit as the original set.
4. **Leaf counts.** Every (k, a) count matches the numpy DP count, for both methods.

## 2. Structure of the extremal sets
* **The 229 orbits.** They are listed in `procgen_edim_20261001_q6_resolving15_orbits.txt`, as the
  canonical minimum 64-bit mask with bit v = vertex v. All stabilizers are trivial, as L1 forces.
  Complementation (paper Prop. 2) gives the 229 orbits of resolving 49-sets.
* **Imbalance.** The distribution of max_i |beta_i(S)| over the 229 orbits is:
  5 (16 orbits), 7 (159), 9 (51), 11 (3). So every resolving 15-set has a coordinate with imbalance
  at least 5; method A finds nothing with a - b <= 3. The paper's set has beta = (-11, 5, -3, 1, 1, 1).

  Heuristic reason (degree-1 analogue of L3): an edge e = (i, x) has total landmark distance
  sum_s d(e,s) = (5k - sum_{j != i} (-1)^{x_j} beta_j) / 2. Large, generic imbalances spread these
  totals and separate edges.
* **The k = 14 near miss.** The unique defect-1 orbit has canonical mask 0x000001810690226d. It fails
  on exactly one antipodal edge pair, {24,26} and {37,39}. Their common histogram (1,1,5,5,1,1) is a
  palindrome, which is exactly what L2 requires of a lone collision.

## 3. Lemmas (PROVED)
**L1 (trivial stabilizer).** Let S be edge-multiset resolving in Q_d with d >= 2, and let
sigma in Aut(Q_d) fix S setwise. Then sigma = id.
* Proof: H^S_{sigma e} = H^{sigma^-1 S}_e = H^S_e. A nontrivial automorphism moves some edge: an
  automorphism fixing every edge as a set fixes every vertex of degree >= 2.
* Consequences: every resolving set has an orbit of size exactly 2^d d!. No set invariant under a
  nontrivial subgroup of Aut(Q_d) can resolve; this rules out unions of subgroup orbits, symmetric
  codes, and similar constructions.

**L2 (antipodal reversal).** For the antipodal edge e-bar, d(e-bar, s) = d - 1 - d(e, s). Hence
H_{e-bar}(r) = H_e(d-1-r).
* Collisions come in pairs {e,f} <-> {e-bar, f-bar}. A collision paired with itself is {e, e-bar},
  and needs H_e to be a palindrome.
* If d is even and |S| is odd, no histogram is a palindrome, so the number of colliding pairs is even.

**L3 (alternating sum).** For e = {u, u + e_i},
sum_r (-1)^r H_e(r) = chi_U(u) * S^(U), where U = (1,...,1) + e_i and S^(U) = sum_{s in S} (-1)^{U.s}.
* Within one direction the alternating sum takes only the two values +S^(U) and -S^(U).
* More generally, the histogram of e is equivalent to the vector of edge-averages of the Walsh layers
  g_w = sum_{|U|=w} chi_U S^(U), for w = 0..d-1.

**L4 (counting).** This improves Prop. 10 from 6 to 7.
* At most 6m edges contain a landmark, i.e. have H(0) >= 1.
* At most 6m edges contain the antipode of a landmark, i.e. have H(5) >= 1 (by L2).
* The remaining >= 192 - 12m edges have histograms on levels 1..4, which allow at most C(m+3, 3) values.
* Since 192 - 12m > C(m+3, 3) for m <= 6, edim_m(Q_6) >= 7.
* For general d the condition is d 2^(d-1) - 2dm <= C(m+d-3, d-3).

**L5 (entropy bound).** If S is resolving with |S| = m, then

log2(d 2^(d-1)) <= sum_{r=0}^{d-2} g(m C(d-1,r) / 2^(d-1)),   where g(mu) = (mu+1) log2(mu+1) - mu log2(mu).

*Proof.* Take e uniform among the edges.
1. The map e -> H_e is injective, so the entropy H(H_e) equals log2 |E|.
2. H_e is determined by its levels 0..d-2, because the levels sum to m.
3. By subadditivity, H(H_e) <= sum over r <= d-2 of H(H_e(r)).
4. E[H_e(r)] = m C(d-1,r) / 2^(d-1), because each vertex is at distance r from exactly d C(d-1,r) edges.
5. Among nonnegative-integer laws with a given mean, the geometric law maximizes entropy, with value g(mean).

*Asymptotics.* Put lambda = ln m.
* At most 1 + sqrt(2 d lambda) levels have mu_r >= 1 (Hoeffding), and each contributes
  <= lambda / ln 2 + 2 bits.
* Levels within sqrt(2 d lambda) of the middle that have mu_r < 1 contribute <= 2.45 bits each,
  O(sqrt(d lambda)) in total.
* All other levels have mu_r <= m^-3 and contribute o(1).
* So the right-hand side is <= sqrt(2 d lambda) * lambda / ln 2 * (1 + o(1)) + O(sqrt(d lambda)).
* This is < d once lambda <= (c - delta) d^(1/3), where c = (ln 2 / sqrt 2)^(2/3) = 0.6215.

Therefore **edim_m(Q_d) >= exp((0.6215 - o(1)) d^(1/3))**, which is superpolynomial in d.

Exact values of the bound (runner S2):

| d | 6 | 7 | 8 | 10 | 12 | 16 | 20 | 32 | 64 | 128 | 256 | 1024 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bound | 4 | 5 | 5 | 6 | 8 | 11 | 15 | 28 | 81 | 275 | 1159 | 50116 |

L5 overtakes L4 from d = 15 on.

## 4. Tournament-structured candidates (FINITE-EXACT, negative)
Q_6 is the n = 5 tiling cube: m = C(4,2) = 6 tiles, and by THM-474 tilings are switching classes.
* **Isomorphism classes.** The 64 tilings fall into 12 tournament isomorphism classes, of sizes
  1, 1, 1, 1, 3, 5, 5, 5, 9, 9, 11, 13. None of the 4095 unions of classes is resolving. The best
  union has 4 colliding pairs, at size 30.
* **H-levels and score classes.** No union of H(T)-levels (7 levels) or of score classes (9 classes)
  is resolving either. For H-levels the reason is structural:
  * the grid transpose G of THM-022 is a nontrivial coordinate permutation of Q_6;
  * H(G T) = H(T), because G T is isomorphic to T^op;
  * so every union of H-levels has nontrivial stabilizer, and L1 excludes it.

Resolving sets are asymmetric objects, and the tournament structure gives no shortcut.

## 5. Q_7 and the growth problem (Open Problems 2 and 4)

| d | edges | lower bound max(L4, L5) | explicit set (VERIFIED) | sparse union bound (VERIFIED) | paper |
|---|------:|:---:|---:|---:|---:|
| 6 | 192 | 7 | **15 = exact** | - | 15 |
| 7 | 448 | 8 | **19** | - | 63 |
| 8 | 1024 | 8 | 26 | - | 115 |
| 9 | 2304 | 8 | 38 | - | 246 |
| 10 | 5120 | 8 | 48 | - | 492 |
| 11 | 11264 | 8 | 65 | 511 | existence |
| 12 | 24576 | 9 | 76 | 576 | existence |
| 13-18 | | 9-13 | - | 667, 781, 910, 1056, 1221, 1408 | existence |
| 20, 24, 28, 32 | | 15, 19, 23, 28 | - | 1837, 2969, 4537, 6638 | existence |

### Open Problem 4: density 1/2 versus sparse
Density-1/2 sets are far from optimal.
* **Explicit sets.** Descending simulated annealing (`procgen_edim_20261001_sad.c`) found sets of
  density 15% (d = 7) down to 1.9% (d = 12). Every set is re-verified exactly; the sets are in the
  runner (S10).
* **Annealing effort.** At d = 7, 2 of 20 long runs reached 19; the best attempt at 18 ended with
  3 colliding pairs. At d = 8..12 the runs were short, so those numbers are only upper bounds.
* **Rigorous sparse bound.** The paper's forest lemma (Lemma 11) holds verbatim for any inclusion
  probability q: the cells are still independent binomials, and forest flows are still forced. We
  computed exact cell formulas for both edge-pair orbit families: parallel pairs with Hamming
  parameter h = 1..d-1, and crossing pairs with h = 0..d-2. Cells and orbit sizes were validated
  against brute force over all edge pairs for d = 4, 5, 6.
* **Size bound.** For given q, let U(q) be the union bound. Then edim_m(Q_d) <= M as soon as
  U(q) < 1 and P(Bin(2^d, q) > M) < 1 - U(q).
* **Check against the paper.** At q = 1/2 our code reproduces the paper's threshold:
  U(10) = 1.307 > 1 > U(11) = 0.156.

### Open Problem 2: growth
* **Lower side.** L5 shows superpolynomial growth.
* **Upper side.** The optimized sparse union bound gives ln M_d / d^(1/3) between 2.76 and 2.80 for
  every d = 11..32 tested. Over this range the curve is strikingly flat.
* **Heuristic.** For random S of size m, a typical edge pair's collision probability is a product of
  about sqrt(2 d ln m) central atoms, each of order m^(-1/2). The union bound over about d^2 4^d pairs
  then succeeds once (ln m)^(3/2) >> sqrt d.
* **Conjecture (OPEN).** ln edim_m(Q_d) = Theta(d^(1/3)). L5 is the lower half. The upper half would
  follow from an asymptotic version of the sparse union bound, the analogue of the paper's Secs. 7-9
  with q -> 0.

## 6. OPEN
* edim_m(Q_7) lies between 8 and 19; edim_m(Q_d) is unknown exactly for every d >= 7. Exhaustive
  search in the style of Sec. 1 is out of reach for Q_7 near k = 19.
* The upper half of the conjecture ln edim_m(Q_d) = Theta(d^(1/3)): a uniform proof that sparse random
  sets work, together with an interval-arithmetic certification of item 7.
* Open Problem 3 (a uniform analytic estimate starting at d = 11) was not addressed.

## 7. Reproduction
Run `python3 -u 04-computation/experiments/procgen_edim_20261001_run.py`.
* Time and memory: the certification run took 2073 s (about 35 min) on one Apple M2 core, with 127 checks and RSS below 150 MB.
* What it does: compiles `procgen_edim_20261001_{q5orbits,searchA,searchB,canon6}.c`, imports
  `procgen_edim_20261001_lib.py`, and prints only to stdout, ending with ALL CHECKS PASSED.
* `--quick` runs a 100-second smoke test with k <= 12. It is not a certification run.
* `procgen_edim_20261001_sad.c` is the annealer; the runner does not need it.
* Its output is `05-knowledge/results/procgen_edim_20261001.out`.
