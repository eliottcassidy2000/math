# Edge multiset dimension of Q_7: no resolving set of size <= 12, so 13 <= edim_m(Q_7) <= 19

**Provenance:** Claude thread session (project thread, branch `claude/project-thread-96889h`), 2026-10-01/02;
temporary identity, no `.machine-id`. Part of the edim2 lane; the growth results are in
[`edge_multiset_dimension_growth_20261002.md`](edge_multiset_dimension_growth_20261002.md), section 7 of which
summarizes this note. Previous bounds: 8 <= edim_m(Q_7) <= 19
([THM-4525](../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md); Allikvere arXiv:2608.09983
had 63 as the upper bound).
**Status:** FINITE-EXACT (computer-assisted exhaustive search, validated as described in Sections 3 and 6).
A blind independent audit (2026-10-02) found no false claim and no gap; it re-ran parts of the search with its own
code, and its wording and provenance fixes are applied here. **Promoted:** nothing (no THM/HYP IDs reserved).
**Labels used:** PROVED (Lemmas 1-3), FINITE-EXACT, EMPIRICAL (annealing).
**Code and evidence:** [`04-computation/edge_multiset_dimension_q7_20261002/`](../../04-computation/edge_multiset_dimension_q7_20261002/)
(C search, Python driver and checkers, run scripts, run records and logs; paths below are relative to it).
The growth runner `04-computation/experiments/edge_multiset_dimension_growth_20261002_run.py` (S10) re-checks the
two 19-sets and the orbit counts from scratch and audits the deposited records; with `--q7` it rebuilds the search
and re-runs the fast subset (Section 7.3).
**Parallel work.** The mac-mini lane note
[`procgen_edim2_20261001_growth_and_uniform_bounds.md`](procgen_edim2_20261001_growth_and_uniform_bounds.md) §6
landed on `origin/main` (12d7b60c, 00:49:53 UTC) while this note was being finished, and was later promoted as
[THM-4534](../../01-canon/theorems/THM-4534-edge-multiset-dimension-grows-like-exp-theta-cube-root-d.md), which
records its Q_7 bound as FINITE-EXACT (lane only). Independently, it proves 11 <= edim_m(Q_7) by an exhaustive
search in a maximal-imbalance normal form (k <= 10), and leaves one case of k = 11 unfinished (a = 11, b = 0 in its
normal form). The search here covers every set of size k <= 12, so it settles that case as well. On Q_6 the
parallel search ran k = 15 with a = 13, 14, 15 of its normal form and recovered the three THM-4525 orbits of
maximal imbalance 11; the search here reproduces all of THM-4525 on Q_6 (nothing for k <= 14, and all 229 orbits at
k = 15; Section 3.4). Its Burnside counts (Q_6, a = 1..8; Q_5, a = 13..15) agree with the ones here where they
overlap.

**Status of k = 13.** The k = 13 search (1,507,949,500,245 leaves) was started at 2026-10-02 00:09:40 UTC and was
still running when this note was committed. Nothing about k = 13 is claimed here; Section 4.4 says how to read
its outcome.

## 1. Summary

**Main result (FINITE-EXACT, computer-assisted).**
* Q_7 has no edge-multiset resolving set of any size k = 1, ..., 12. Hence **edim_m(Q_7) ≥ 13**.
* With the known resolving 19-set: **13 ≤ edim_m(Q_7) ≤ 19**. The previous bounds were 8 ≤ edim_m(Q_7) ≤ 19
  (11 ≤ edim_m(Q_7) in the parallel procgen_edim2 note).

**Method.**
* Every k-set is reduced to the minimum-imbalance normal form (Lemma 1, proved in Section 2.2): one
  Aut(Q_6)-orbit representative A in the half x_6 = 0, and an arbitrary b-set B in the half x_6 = 1.
* For every (k, a), every representative A and every admissible B is tested exactly.
* No pruning on partial collisions is used, because resolvability is not monotone.

**Validation (all passed).**
* Q_6 reproduction:
  * the same code, run on Q_6 for k = 1..15, finds nothing for k ≤ 14;
  * for k = 15 it finds 585 resolving leaves, which canonicalize under all 46080 automorphisms to exactly the
    229 orbits of the deposited list, all with trivial stabilizer;
  * the alternative hot-list engine (engine 1) reproduces k = 14 and 15 exactly. It replaces only the hot-list
    scan; the search, the keys and the full test are shared with engine 0.
* Leaf counts: for every (k, a), with d = 6 (k ≤ 15) and d = 7 (k ≤ 12), the C leaf count
  equals an independent numpy dynamic-programming count.
* Python checker:
  * it confirms the known Q_7 19-set and the paper's Q_6 15-set;
  * it agrees with the C checker on 1,636 sets.
  * No Q_7 search found a set, so there was nothing else to re-verify.
* Further checks:
  * the orbit-representative files equal the Burnside counts and pass an independent brute-force canonicity
    check;
  * verification mode recomputed 154,155,615 Q_7 leaves from scratch with 0 discrepancies;
  * 1,250 random sets confirm the normal-form construction.

**Timing (2 cores, shared machine).**
* Q_7, k ≤ 12: 4,376.5 CPU s in 37.4 min wall. Of this, k = 12 took 4,029 CPU s in 33.7 min wall.
* Speed: 11.7-14.7 ns per leaf; 3.04e11 leaves in total.
* k = 13 has 1.508e12 leaves (DP). At the observed rate of about 100 of the 19,735 a = 7 representatives
  per minute (14-15 ns per leaf), completion is expected around 03:10-03:40 UTC.
* k = 14 has 1.429e13 leaves, about 55 CPU hours. It was not attempted.

**Annealing.** Only the calibration runs were done:
* Q_6, k = 15: 55 resolving sets found in 60 s, in 44 distinct orbits. All 44 lie in the deposited list of
  229 orbits.
* Q_7, k = 19: one resolving 19-set found in 120 s. It is **not** equivalent to the known 19-set.
* The production runs for sizes 18 and 17 were **not run**; the CPU was kept for k = 13. The script is
  ready (Section 5).

## 2. Method

### 2.1 Notation and invariance
* Q_d has vertex set {0,1}^d, written as d-bit integers. Its edges are {u, u + 2^i} with bit i of u equal to 0.
* For an edge e and a vertex s, d(e,s) = min(d_H(u,s), d_H(u',s)), where u, u' are the endpoints of e.
* The histogram of e with respect to S is H_e(r) = #{s in S : d(e,s) = r}, for r = 0..d-1.
* S is (edge-multiset) resolving if the d 2^(d-1) histograms are pairwise distinct.
* Aut(Q_d) = { x -> pi(x xor t) }, where pi permutes coordinates. It has order 2^d d!.
* Every sigma in Aut(Q_d) maps edges to edges bijectively and preserves Hamming distance, so
  d(sigma e, sigma s) = d(e, s) and H^{sigma S}_{sigma e} = H^S_e. Hence **S is resolving iff sigma S is resolving.**
* Imbalances: beta_i(S) = #{s in S : s_i = 0} - #{s in S : s_i = 1}.
  * A translation by t replaces beta_i by (-1)^{t_i} beta_i.
  * A coordinate permutation permutes the beta_i.
  * So every automorphism permutes the multiset {|beta_i|}, and every automorphism acting only on coordinates
    0..d-2 maps beta_{d-1} to itself.
* Note beta_i(S) ≡ |S| (mod 2).

### 2.2 Lemma 1 (normal form; minimum-imbalance form "B")
**Statement.** Write a vertex of Q_d as (x, x_{d-1}) with x in Q_{d-1}. Let R_a be a set containing one
representative of every Aut(Q_{d-1})-orbit of a-subsets of Q_{d-1}. Every k-set S ⊆ V(Q_d) is
Aut(Q_d)-equivalent to a set T = A x {0} ∪ B x {1} with all of the following:
1. a = |A| ≥ b = |B| and a + b = k;
2. a - b = min_i |beta_i(T)|, so |beta_i(T)| ≥ a - b for i = 0..d-2;
3. A ∈ R_a.

**Proof.**
1. Let m = min_i |beta_i(S)|, attained at the coordinate j.
2. Let S_1 = tau S, where tau is the coordinate transposition (j, d-1) (the identity if j = d-1).
   Then beta_{d-1}(S_1) = beta_j(S), and the multiset {|beta_i(S_1)|} equals {|beta_i(S)|}.
3. If beta_{d-1}(S_1) < 0, let S_2 = S_1 xor 2^{d-1}, i.e. flip the last coordinate. This negates
   beta_{d-1} and leaves every other beta_i unchanged. Otherwise let S_2 = S_1.
4. Now beta_{d-1}(S_2) = m ≥ 0. Write S_2 = A' x {0} ∪ B' x {1}. Then |A'| - |B'| = beta_{d-1}(S_2) = m.
   So a := |A'| ≥ b := |B'|, a - b = m, and |beta_i(S_2)| ≥ m for every i.
5. Choose h ∈ Aut(Q_{d-1}) with h(A') ∈ R_a. Let ĥ(x, x_{d-1}) = (h(x), x_{d-1}); this is an element of Aut(Q_d).
6. Then ĥ S_2 = h(A') x {0} ∪ h(B') x {1}. Since ĥ acts only on coordinates 0..d-2, it permutes and negates
   beta_0..beta_{d-2} and fixes beta_{d-1}. So T := ĥ S_2 satisfies conditions 1-3.
∎

**Corollary.** If no pair (A, B) with A ∈ R_a, |B| = b and condition 2 is resolving, for any
a = ceil(k/2)..k, then Q_d has no resolving k-set. Since resolvability is not monotone in S, the search
tests every such pair; nothing is pruned because of a partial collision.

**Parity remark.** delta := a - b ≡ k (mod 2). When delta ≤ 1, condition 2 is vacuous:
* for delta = 0 there is nothing to check;
* for delta = 1, k is odd, so every beta_i is odd and |beta_i| ≥ 1 automatically.

### 2.3 Lemma 2 (weight lemma: some (k, a) need no representatives)
Let t = a - 2b > 0.
* Condition 2 together with |beta_i(B)| ≤ b forces |beta_i(A)| ≥ t for every i ≤ d-2.
* Then the minority value occurs at most (a - t)/2 = b times in each coordinate.
* Translate A so that the majority value is 0 in each coordinate. The total weight of A is then at most (d-1) b.
* Any a distinct vertices of Q_{d-1} have total weight at least W_min(a), obtained by taking the lightest
  vertices first.
* So **W_min(a) > (d-1) b excludes the pair (k, a).**

For d = 7: W_min(11) = 14, W_min(12) = 16, W_min(13) = 18 and W_min(14) = 20. Hence:
* for k ≤ 13, every a ≥ 11 is excluded;
* for k = 14, every a ≥ 12 is excluded.

So for k ≤ 13 only the representative files with a ≤ 10 are needed. The file for a = 11 was also run
whenever it exists; it gives 0 leaves, which agrees with the lemma.

### 2.4 Orbit representatives: orderly generation (orbreps.c)
**Canonical form.** can(T) is the lexicographically smallest sorted tuple in the Aut(Q_n)-orbit of T.

**Lemma 3 (heredity).** If T = (t_1 < ... < t_a) is canonical, then T' = (t_1, ..., t_{a-1}) is canonical.

*Proof.*
1. Suppose some g ∈ Aut(Q_n) has sorted g(T') = W <_lex T'. Let j be the first position with
   w_j ≠ t_j; then w_j < t_j and w_i = t_i for i < j.
2. sorted g(T) is the merge of W with the single element g(t_a). For i ≤ j, its i-th entry is
   ≤ w_i ≤ t_i, and its j-th entry is ≤ w_j < t_j.
3. So sorted g(T) differs from T at some first position i ≤ j and is smaller there. Hence g(T) <_lex T,
   which contradicts the canonicity of T. ∎

**Consequence.** Every canonical (a+1)-tuple equals a canonical a-tuple plus one element above its maximum.
* Level a+1 is therefore generated by testing all extensions T + x, with x > max T, of the level-a list.
* By induction, every orbit appears exactly once.

**Canonicity test.** Let M be the bit mask of T.
* Lex order of equal-size sets: X <_lex Y iff the lowest element of X Δ Y lies in X.
* If 0 ∉ T, then T is not canonical, because translating any element to 0 gives a smaller set.
* If 0 ∈ T, only images containing 0 can be smaller. These are exactly the maps x -> pi(x xor s) with s ∈ T.
  So the test loops over s ∈ T and over all n! coordinate permutations pi. The permutations are enumerated
  by the Steinhaus-Johnson-Trotter sequence of adjacent transpositions, each applied to the 2^n-bit mask
  as one delta swap; a startup self-check confirms that all n! permutations are visited.
* Exact shortcut. Let t_2 be the second-smallest element of T, and let w_s be the distance from s to its
  nearest other element of T. The second-smallest element of any image under translation s is at least
  2^{w_s} - 1, and that value is attained.
  * If 2^{w_s} - 1 > t_2, translation s cannot give a smaller image, so s is skipped.
  * If 2^{w_s} - 1 < t_2, T is not canonical.

**Checks.**
1. A self-test compares `is_canon` with a brute-force canonical form over all group elements on random sets:
   3000 sets for n = 4 and n = 5, and 1500 for n = 6, with 0 disagreements.
2. The level counts equal the Burnside counts computed by `py/burnside.py`, which uses the cycle index of the
   explicit group and shares no code with the generator. This holds for n = 5, a = 0..32, and n = 6, a = 0..11.
3. `src/repcheck.c` re-checks every listed representative independently, from explicit vertex tables. It
   uses all 2^n n! images, or in "fast" mode the images that contain 0, by the argument above. It also checks
   distinctness (see Section 3.1).
4. A list of pairwise distinct canonical sets is pairwise inequivalent. If in addition its length equals the
   Burnside count, it contains every orbit.

### 2.5 Search (edimsearch.c)
**Enumeration.** For each a and each representative A, the program enumerates every b-subset B of Q_{d-1}
by a DFS in increasing vertex order.
* A subtree is skipped only when some coordinate i satisfies beta_i + r < delta and beta_i - r > -delta,
  where r elements remain. Then every completion has |beta_i| < delta, so the subtree contains no leaf
  of the domain.
* For the last element, the admissible vertices form the exact mask of vertices that give |beta_i| ≥ delta
  for every i.
* The number of leaves visited is therefore exactly the size of the domain. This is checked against the
  independent DP count for every (k, a).

**Keys.**
* key(e) = sum over s of C(d(e,s)), where C(t) = 32^t for t ≤ d-2 and C(d-1) = 0.
* This stores levels 0..d-2 in 5-bit fields. Level d-1 is implied, because every histogram sums to k.
* For k ≤ 31, key(e) = key(f) iff H_e = H_f.
* Keys are updated incrementally along the DFS, as one vectorised row addition per DFS node.

**Leaf test (exact).** At a node where only the last element v remains, the candidates form a 64-bit mask.
* The pair (e, f) collides after adding v iff C(d(f,v)) - C(d(e,v)) = de, where de = key(e) - key(f) at the node.
* (x, y) -> C(y) - C(x) is injective on x ≠ y, because base-32 digits do not carry. So the killed set is:
  * {v : d(e,v) = d(f,v)} if de = 0;
  * the single class {v : d(e,v) = x, d(f,v) = y} if de = C(y) - C(x);
  * empty otherwise.
* These sets come from precomputed masks DM[e][x]. A perfect hash decodes de branch-free.

**Engines.** Each engine evaluates a hot list of edge pairs (capacity 1024).
* Engine 0 ("bulk"): the list is scanned in branch-free blocks of 8 until all candidates are dead. The pair
  that completes the kill is moved to the front.
* Engine 1 ("per-leaf"): for each candidate separately, the list is scanned with the direct comparison
  key(e) + C(d(e,v)) == key(f) + C(d(f,v)), with move-to-front.

**Fallback and reporting.**
* Every candidate the hot list does not kill gets a full test: all edge keys go into a hash table.
* The first colliding pair found is inserted into the hot list. The program asserts that this pair kills v
  under the kill logic.
* A candidate without a collision is resolving. Before it is printed, it is re-verified from scratch, with keys
  recomputed from the definition (minimum over both endpoints, no projection lemma).

So every leaf is either refuted by an exactly computed collision or tested completely.

### 2.6 Driver, checkpointing, validation scripts
* `py/driver.py D K` splits each representative file into chunks of about 120 s of CPU.
* It runs at most 2 chunks at once, each as `nice -n 10 src/edimsearch D K a repfile first last -e 0 -H 1024`.
* After each completed chunk it appends one JSON record to `runs/<run>/chunks.jsonl` and calls fsync.
* On a restart, only the parts of each file that are not yet covered are scheduled.
* At the end it checks that the chunks of every a tile [0, nrep) exactly, then writes
  `summary_k<K>.json` and logs `COMPLETE ...`.
* `py/validate.py D RUN k...` repeats the tiling check and compares the C leaf count of every (k, a) with
  `py/leafcount_dp.py`. That script counts the normal-form domain by a DP over imbalance vectors, without
  enumerating B and without sharing code with the C search.
* For a missing representative file, `validate.py` re-derives the weight-lemma exclusion.
* It re-verifies every RESOLVING set with the pure-Python checker `py/pycheck.py`, canonicalizes it with
  `src/canon` (all 2^d d! automorphisms), and writes `validation_k<k>.json`.

## 3. Validation

### 3.1 Orbit representatives
* **Burnside.** `src/orbreps` level counts equal the Burnside counts from `py/burnside.py`
  (`data/burnside_5.txt`, `data/burnside_6.txt`):
  * Q_5: every level a = 0..32;
  * Q_6, a = 0..11: 1, 1, 6, 16, 103, 497, 3253, 19735, 120843, 681474, 3561696, 16938566.
* **Self-test of `is_canon`** against a brute-force canonical form: 3000 random sets for n = 4 and n = 5, and
  1500 for n = 6, with 0 disagreements (`logs/iscanon_selftest.out`; the blind audit re-ran all three).
* **`src/repcheck`** (`logs/repcheck.out`, 7.3 min) checked every listed set for: correct size, lex-minimal in
  its orbit, and no duplicates.
  * Group used: the full group (3840 or 46080 explicit elements) for Q_5, a = 0..32, and Q_6, a = 0..8. For
    Q_6, a = 9 and 10 it used the "fast" mode: the a·720 images that contain 0, which by the elementary
    argument in Section 2.4 are the only ones that can be lex-smaller.
  * Result: noncanonical = wrongsize = duplicates = 0 for all 44 files.
  * Negative control (`logs/repcheck_negative.out`): a file with {0,1,3}, {0,1,2} twice and {1,2,4}. Both
    modes correctly report 2 non-canonical sets and 1 duplicate.
* **Conclusion.** Each file used is pairwise inequivalent (distinct canonical sets) and complete (count =
  Burnside). The Q_6 file for a = 11 was not repchecked. It is not needed for k ≤ 13 (Lemma 2), and its
  search gave 0 leaves anyway.

### 3.2 Checkers
* **`py/pycheck.py selftest`** (`logs/pycheck_selftest.out`). It computes histograms as tuples straight from
  the definition, with no projection lemma and no hashing.
  * The known Q_7 19-set {5,57,54,104,35,109,115,49,6,39,55,102,15,32,85,97,21,113,41} is resolving (defect 0).
  * The paper's Q_6 15-set is resolving.
  * Removing any single landmark from the 19-set gives defect ≥ 14:
    [33,27,54,45,31,30,26,30,28,45,33,39,32,30,16,17,34,14,24].
* **`py/xcheck.py 400 7`** (`logs/xcheck.out`) compares the C check mode (`edimsearch d 0 -x`) with `pycheck`:
  * The C check mode recomputes keys from scratch and asserts that its hash test and a sort test agree.
  * Sets: 1,636 in total, for d = 4, 5, 6, 7. They include random sets, the two resolving sets, and their
    34 one-landmark-deleted near misses.
  * The defect (#edges − #distinct histograms) agrees on every set: `XCHECK ALL AGREE`.

### 3.3 Verification mode (Q_7, production engine 0)
**What `-v` does at every last-level node.**
* It recomputes the keys of the current set from scratch, using the minimum over both endpoints, and
  compares them with the incremental keys.
* It computes the true collision status of every candidate leaf from scratch keys.
* It compares the kill mask of every hot pair with a direct per-candidate recomputation.
* It aborts if any kill hits a leaf that has no collision.
* It aborts if any full test disagrees with the truth.

All samples finished without FATAL.

| k | a | b | delta | reps (indices) | leaves verified | sec | log |
|---|---|---|---|---|---:|---:|---|
| 12 | 6 | 6 | 0 | 1 (1777) | 74,974,368 | 632.9 | logs/verify_d7_big1.out |
| 13 | 7 | 6 | 1 | 1 (12345) | 74,974,368 | 678.5 | logs/verify_d7_big2.out |
| 12 | 7 | 5 | 2 | 2 (9000-9001) | 2,641,978 | 25.8 | logs/verify_d7_small.out |
| 12 | 8 | 4 | 4 | 20 (70000-70019) | 8,820 | 0.1 | " |
| 13 | 8 | 5 | 3 | 10 (60000-60009) | 1,550,378 | 19.3 | " |
| 13 | 9 | 4 | 5 | 200 (300000-300199) | 5,252 | 0.1 | " |
| 13 | 10 | 3 | 7 | 3000 (0-2999; 34 feasible) | 451 | 0.0 | " |
| | | | | **total** | **154,155,615** | | |

An earlier sample also passed for both engines: d = 7, k = 11, a = 6, two representatives, 15,249,024 leaves
each. It has no log file, and nothing rests on it.

### 3.4 Reproduction of Q_6 (edim_m(Q_6) = 15, THM-4525)
Run with engine 0: `runs/d6_e0`, `logs/d6_e0.out`, `logs/validate_d6_e0.out`.

| k | leaves (C search) | leaves (DP) | equal | resolving leaves | orbits | CPU s | wall s | ns/leaf | validation |
|---|---:|---:|:-:|---:|---:|---:|---:|---:|:-:|
| 1 | 1 | 1 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 2 | 32 | 32 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 3 | 160 | 160 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 4 | 2,506 | 2,506 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 5 | 4,971 | 4,971 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 6 | 52,535 | 52,535 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 7 | 234,240 | 234,240 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 8 | 1,808,073 | 1,808,073 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 9 | 4,767,589 | 4,767,589 | yes | 0 | 0 | 0.1 | 0 | - | ok |
| 10 | 29,955,834 | 29,955,834 | yes | 0 | 0 | 0.6 | 1 | - | ok |
| 11 | 96,667,005 | 96,667,005 | yes | 0 | 0 | 2.1 | 2 | 21.72 | ok |
| 12 | 491,865,822 | 491,865,822 | yes | 0 | 0 | 13.0 | 10 | 26.43 | ok |
| 13 | 1,234,172,511 | 1,234,172,511 | yes | 0 | 0 | 31.0 | 28 | 25.12 | ok |
| 14 | 5,359,115,934 | 5,359,115,934 | yes | 0 | 0 | 158.1 | 93 | 29.50 | ok |
| 15 | 13,146,176,599 | 13,146,176,599 | yes | 585 | 229 | 397.8 | 214 | 30.26 | ok |

Per-(k, a) detail for k = 14 and 15 (every (k, a) for k = 1..15 is in `runs/d6_e0/validation_k*.json`):

| k | a | b | delta | reps | chunks | leaves (C) | leaves (DP) | resolving | CPU s |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| 14 | 7 | 7 | 0 | 1,326 | 5 | 4,463,125,056 | 4,463,125,056 | 0 | 117.2 |
| 14 | 8 | 6 | 2 | 3,779 | 2 | 885,632,884 | 885,632,884 | 0 | 38.8 |
| 14 | 9 | 5 | 4 | 9,013 | 2 | 10,349,317 | 10,349,317 | 0 | 2.1 |
| 14 | 10 | 4 | 6 | 19,963 | 2 | 8,677 | 8,677 | 0 | 0.0 |
| 14 | 11..14 | 3..0 | 8..14 | 38,073 / 65,664 / 98,804 / 133,576 | 2 each | 0 | 0 | 0 | 0.0 |
| 15 | 8 | 7 | 1 | 3,779 | 8 | 12,719,569,824 | 12,719,569,824 | 572 | 352.5 |
| 15 | 9 | 6 | 3 | 9,013 | 3 | 424,545,625 | 424,545,625 | 13 | 44.3 |
| 15 | 10 | 5 | 5 | 19,963 | 2 | 2,060,825 | 2,060,825 | 0 | 0.9 |
| 15 | 11 | 4 | 7 | 38,073 | 2 | 325 | 325 | 0 | 0.0 |
| 15 | 12..15 | 3..0 | 9..15 | 65,664 / 98,804 / 133,576 / 158,658 | 2 each | 0 | 0 | 0 | 0.0 |

**Findings.**
* No resolving set exists for k ≤ 14.
* For k = 15 there are 585 resolving leaves (572 with a = 8, 13 with a = 9).
  * All 585 were re-verified by `pycheck`.
  * Under all 46080 automorphisms of Q_6 they canonicalize to exactly **229 orbits**, every one with trivial
    stabilizer.
  * The recomputed list has the same 229 masks as the deposited
    [`procgen_edim_20261001_q6_resolving15_orbits.txt`](procgen_edim_20261001_q6_resolving15_orbits.txt)
    (which adds a 4-line comment header), so it is not stored again.
* Engine 1 (`runs/d6_e1`, `logs/validate_d6_e1.out`) re-ran k = 14 and 15. It gives the same per-(k, a) leaf
  counts and the same 585 resolving leaves, with an identical 229-orbit list.

| k | leaves (C search) | leaves (DP) | equal | resolving leaves | orbits | CPU s | wall s | ns/leaf | validation |
|---|---:|---:|:-:|---:|---:|---:|---:|---:|:-:|
| 14 (engine 1) | 5,359,115,934 | 5,359,115,934 | yes | 0 | 0 | 183.1 | 97 | 34.17 | ok |
| 15 (engine 1) | 13,146,176,599 | 13,146,176,599 | yes | 585 | 229 | 445.0 | 229 | 33.85 | ok |

`leafcount_dp.py` was later changed to read large files in chunks, to keep RAM low for the 16.9M-set Q_6
a = 11 file. After the change, all Q_6 validations were re-run (`logs/validate_d6_e*_rerun.out`), and the
validation JSON files are byte-identical to the earlier ones.

### 3.5 Normal-form reduction, empirically
**What `py/nf_test.py` does** (`logs/nf_test.out`), for random k-sets S:
* It applies the construction of the proof of Lemma 1. The lex-min of A' is computed by brute force over
  all of Aut(Q_{d-1}).
* It checks that A is in the representative file (binary search), that a ≥ b, that a − b = min |beta|, and
  that |beta_i| ≥ a − b.
* It checks with `src/canon` that the normal form lies in the same Aut(Q_d)-orbit as S.

**Results** (0 failures in 1,250 sets):

| d | k | sets | domain failures | same orbit | a-distribution |
|---|---|---:|---:|---|---|
| 6 | 14 | 400 | 0 | all | {7: 324, 8: 76} |
| 6 | 15 | 400 | 0 | all | {8: 395, 9: 5} |
| 7 | 11 | 150 | 0 | all | {6: 148, 7: 2} |
| 7 | 12 | 150 | 0 | all | {6: 129, 7: 21} |
| 7 | 13 | 150 | 0 | all | {7: 150} |

## 4. Q_7 results

### 4.1 Per k (run `runs/d7_e0`, engine 0, hot-list capacity 1024, 2 workers; `logs/d7_e0_k1-12.out`, `logs/validate_d7_e0.out`)

| k | leaves (C search) | leaves (DP) | equal | resolving leaves | orbits | CPU s | wall s | ns/leaf | validation |
|---|---:|---:|:-:|---:|---:|---:|---:|---:|:-:|
| 1 | 1 | 1 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 2 | 64 | 64 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 3 | 384 | 384 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 4 | 12,154 | 12,154 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 5 | 32,283 | 32,283 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 6 | 686,390 | 686,390 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 7 | 4,300,936 | 4,300,936 | yes | 0 | 0 | 0.0 | 0 | - | ok |
| 8 | 68,340,374 | 68,340,374 | yes | 0 | 0 | 0.6 | 1 | - | ok |
| 9 | 317,832,905 | 317,832,905 | yes | 0 | 0 | 2.7 | 2 | 8.50 | ok |
| 10 | 4,141,968,201 | 4,141,968,201 | yes | 0 | 0 | 50.1 | 42 | 12.10 | ok |
| 11 | 25,086,535,792 | 25,086,535,792 | yes | 0 | 0 | 293.9 | 177 | 11.72 | ok |
| 12 | 273,963,417,196 | 273,963,417,196 | yes | 0 | 0 | 4029.2 | 2022 | 14.71 | ok |
| **1-12** | **303,583,126,680** | | yes | **0** | | **4376.5** | **2244** | | ok |

**Notes.**
* Wall time is the time from the `start` line to the `COMPLETE` line in `logs/d7_e0_k1-12.out`.
* The whole run k = 1..12 took 23:10:44 to 23:48:08 UTC: 37.4 min wall and 4,376.5 CPU s.
* No resolving leaf was found for any k ≤ 12, so no set needed re-verification.

### 4.2 Per (k, a)
A row marked "weight lemma" has no representative file. Lemma 2 excludes it, and `validate.py` re-derived the
exclusion. The Q_6 file for a = 11 exists and was searched for k = 11 and 12: it gives 0 leaves, in
agreement with the DP and with Lemma 2.

| k | a | b | delta | reps | chunks | leaves (C) | leaves (DP) | resolving | CPU s |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 0 | 0.0 |
| 2 | 1 | 1 | 0 | 1 | 1 | 64 | 64 | 0 | 0.0 |
| 2 | 2 | 0 | 2 | 6 | 2 | 0 | 0 | 0 | 0.0 |
| 3 | 2 | 1 | 1 | 6 | 3 | 384 | 384 | 0 | 0.0 |
| 3 | 3 | 0 | 3 | 16 | 2 | 0 | 0 | 0 | 0.0 |
| 4 | 2 | 2 | 0 | 6 | 3 | 12,096 | 12,096 | 0 | 0.0 |
| 4 | 3 | 1 | 2 | 16 | 2 | 58 | 58 | 0 | 0.0 |
| 4 | 4 | 0 | 4 | 103 | 2 | 0 | 0 | 0 | 0.0 |
| 5 | 3 | 2 | 1 | 16 | 3 | 32,256 | 32,256 | 0 | 0.0 |
| 5 | 4 | 1 | 3 | 103 | 3 | 27 | 27 | 0 | 0.0 |
| 5 | 5 | 0 | 5 | 497 | 2 | 0 | 0 | 0 | 0.0 |
| 6 | 3 | 3 | 0 | 16 | 3 | 666,624 | 666,624 | 0 | 0.0 |
| 6 | 4 | 2 | 2 | 103 | 2 | 19,755 | 19,755 | 0 | 0.0 |
| 6 | 5 | 1 | 4 | 497 | 3 | 11 | 11 | 0 | 0.0 |
| 6 | 6 | 0 | 6 | 3,253 | 3 | 0 | 0 | 0 | 0.0 |
| 7 | 4 | 3 | 1 | 103 | 3 | 4,291,392 | 4,291,392 | 0 | 0.0 |
| 7 | 5 | 2 | 3 | 497 | 2 | 9,540 | 9,540 | 0 | 0.0 |
| 7 | 6 | 1 | 5 | 3,253 | 2 | 4 | 4 | 0 | 0.0 |
| 7 | 7 | 0 | 7 | 19,735 | 2 | 0 | 0 | 0 | 0.0 |
| 8 | 4 | 4 | 0 | 103 | 3 | 65,443,728 | 65,443,728 | 0 | 0.5 |
| 8 | 5 | 3 | 2 | 497 | 2 | 2,892,245 | 2,892,245 | 0 | 0.1 |
| 8 | 6 | 2 | 4 | 3,253 | 2 | 4,400 | 4,400 | 0 | 0.0 |
| 8 | 7 | 1 | 6 | 19,735 | 2 | 1 | 1 | 0 | 0.0 |
| 8 | 8 | 0 | 8 | 120,843 | 2 | 0 | 0 | 0 | 0.0 |
| 9 | 5 | 4 | 1 | 497 | 3 | 315,781,872 | 315,781,872 | 0 | 2.6 |
| 9 | 6 | 3 | 3 | 3,253 | 2 | 2,049,726 | 2,049,726 | 0 | 0.1 |
| 9 | 7 | 2 | 5 | 19,735 | 2 | 1,307 | 1,307 | 0 | 0.0 |
| 9 | 8 | 1 | 7 | 120,843 | 2 | 0 | 0 | 0 | 0.0 |
| 9 | 9 | 0 | 9 | 681,474 | 2 | 0 | 0 | 0 | 0.0 |
| 10 | 5 | 5 | 0 | 497 | 3 | 3,789,382,464 | 3,789,382,464 | 0 | 41.7 |
| 10 | 6 | 4 | 2 | 3,253 | 2 | 351,552,834 | 351,552,834 | 0 | 8.1 |
| 10 | 7 | 3 | 4 | 19,735 | 2 | 1,032,755 | 1,032,755 | 0 | 0.2 |
| 10 | 8 | 2 | 6 | 120,843 | 2 | 148 | 148 | 0 | 0.0 |
| 10 | 9 | 1 | 8 | 681,474 | 2 | 0 | 0 | 0 | 0.0 |
| 10 | 10 | 0 | 10 | 3,561,696 | 2 | 0 | 0 | 0 | 0.1 |
| 11 | 6 | 5 | 1 | 3,253 | 5 | 24,802,537,536 | 24,802,537,536 | 0 | 273.4 |
| 11 | 7 | 4 | 3 | 19,735 | 2 | 283,657,681 | 283,657,681 | 0 | 19.5 |
| 11 | 8 | 3 | 5 | 120,843 | 2 | 340,568 | 340,568 | 0 | 0.1 |
| 11 | 9 | 2 | 7 | 681,474 | 2 | 7 | 7 | 0 | 0.0 |
| 11 | 10 | 1 | 9 | 3,561,696 | 2 | 0 | 0 | 0 | 0.1 |
| 11 | 11 | 0 | 11 | 16,938,566 | 2 | 0 | 0 | 0 | 0.8 |
| 12 | 6 | 6 | 0 | 3,253 | 31 | 243,891,619,104 | 243,891,619,104 | 0 | 3256.7 |
| 12 | 7 | 5 | 2 | 19,735 | 9 | 29,900,496,510 | 29,900,496,510 | 0 | 749.5 |
| 12 | 8 | 4 | 4 | 120,843 | 2 | 171,234,804 | 171,234,804 | 0 | 22.1 |
| 12 | 9 | 3 | 6 | 681,474 | 2 | 66,778 | 66,778 | 0 | 0.1 |
| 12 | 10 | 2 | 8 | 3,561,696 | 2 | 0 | 0 | 0 | 0.1 |
| 12 | 11 | 1 | 10 | 16,938,566 | 2 | 0 | 0 | 0 | 0.8 |
| 12 | 12 | 0 | 12 | - | - | 0 (weight lemma: W_min(12)=16 > 0) | 0 | 0 | 0 |

### 4.3 Performance, and the cost of k = 13 and 14
**Engine 0 speed.**
* Engine 0 runs at 11.7-14.7 ns per leaf on Q_7 (about 10.7 leaves per last-level node, about 28 hot pairs
  scanned per node).
* The machine was shared: load average about 3 on 4 cores from other jobs. Clean single-chunk
  benchmarks give 11.7-13.3 ns.
* Engine 1 is about 1.1x slower on Q_6 (34 vs 30 ns per leaf). Its speed on Q_7 was not measured on a deposited
  run (the planned engine-1 re-run of Q_7 was cancelled; Section 8).

**Block-size benchmark.** Block sizes 1, 2, 3, 4 and 8 were benchmarked for the engine-0 scan, as was swapping
instead of move-to-front. The differences were within ±5% noise, so the validated binary was kept unchanged.

**Domain sizes from the DP** (`logs/dp_d7_k13.out`, `logs/dp_d7_k14.out`):

| k | a=7 | a=8 | a=9 | a=10 | a=11 | total | est. CPU (13-15 ns/leaf) | est. wall, 2 cores |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 | 1,479,619,152,480 | 28,261,148,980 | 69,191,803 | 6,982 | 0 | **1,507,949,500,245** | 5.4-6.3 h | 2.7-3.2 h |
| 14 | 12,259,701,549,120 | 2,016,262,554,454 | 18,959,624,460 | 17,799,808 | 293 | **14,294,941,528,135** | 52-60 h | 26-30 h |

For k = 14, a ≥ 12 is excluded by Lemma 2. Only k = 13 was started.

### 4.4 k = 13 (running at commit time; nothing claimed)
* `runs/post12.sh` started k = 13 at 2026-10-02 00:09:40 UTC, only after every check of Sections 3 and 4
  passed (`logs/post12.out`), as `runs/run_d7.sh 13 13`: 2 workers on the same binary, into `runs/d7_e0`.
* The a = 7 file has 19,735 representatives and carries 98% of the work. Progress lines look like
  `2026-10-02 00:11:42 done a=7 [2,108) leaves=7947283008 found=0 sec=121.0 ns/leaf=15.23 ...`.
* The records committed here stop at k = 12; the k = 13 chunk records are not included.

**How to read the outcome.**
* Success is the line `COMPLETE D=7 K=13 leaves=1507949500245 found=0 ...`; the leaves value must be exactly the
  DP total 1,507,949,500,245 (`logs/dp_d7_k13.out`). Then `python3 py/validate.py 7 runs/d7_e0 13` (about 2 min)
  must print `VALIDATE D=7 k=13 ok=True leaves C=1507949500245 DP=1507949500245 resolving_leaves=0 ...`.
  Only then does edim_m(Q_7) >= 14 follow.
* A resolving 13-set would be printed as `FOUND RESOLVING d=7 k=13 set=v1,...,v13` (vertices 0..127, bit i of v
  = coordinate i), already re-verified from scratch by the C program. `python3 py/pycheck.py 7 <<< "v1,...,v13"`
  must then print `RESOLVING k=13 defect=0`, and edim_m(Q_7) = 13 would follow.
* A failure shows as `ERROR chunk failed ...` or `ERROR inconsistent chunk summary ...`, then `FAILED k=13`.

**Re-running k = 13 from scratch.** `runs/run_d7.sh 13 13` after the build and the representative files of
Section 7.3 (about 5.4-6.3 CPU hours, 2.7-3.2 h on 2 cores). The driver is resumable: completed chunks are
skipped and a torn last JSON line is redone. Never run two drivers on the same run directory; overlapping chunks
fail the tiling check.

## 5. Annealing (EMPIRICAL)

**`src/anneal.c`.**
* Keys are exact, as in the search. The cost is the number of colliding edge pairs, sum over classes of
  C(m, 2).
* A move swaps a landmark with a non-landmark. With probability 0.5 the move is "focused": the new vertex
  separates a recorded colliding pair.
* Acceptance is Metropolis, with geometric cooling from T0 = 2 to T1 = 0.05 over 2M moves per restart,
  followed by a random restart.
* A cost-0 set is re-verified from scratch, by sorting full histograms, before it is printed.
* `py/anneal_summary.py` re-verifies every reported set with `pycheck` and canonicalizes it with `src/canon`.

**Calibration** (`logs/anneal_*_calib.out`, `logs/anneal_summary_calib.out`):

| graph | k | time | restarts | resolving sets found | Python-verified | distinct orbits |
|---|---|---:|---:|---:|---|---:|
| Q_6 | 15 | 60 s | 92 | 55 | yes | 44 (all among the 229 deposited orbits) |
| Q_7 | 19 | 120 s | 29 | 1 | yes | 1 |

**Calibration findings.**
* The new Q_7 19-set is {2,4,21,22,32,38,42,44,47,54,70,79,84,90,110,114,116,120,126}.
* Its canonical form is `00000000000101802802410483841c48`, with stabilizer order 1. The known 19-set has
  canonical form `000000000001018a20020ae900142005`, so the two sets are inequivalent.
* The other Q_7 k = 19 restarts ended with 3-5 colliding pairs. So the annealer finds Q_7 resolving 19-sets,
  but only in about 1 restart out of 30.

**Sizes 18 and 17: not run.**
* `runs/anneal.sh 900` is the planned run: 2 × 15 min for each size, one with 2M and one with 8M moves per
  restart, about 1 core-hour in total.
* It was cancelled to keep the CPU for k = 13. So there is **no** annealing
  result for 18 or 17, and no best defect to report.
* To run it later, after k = 13:
  1. `runs/anneal.sh 900`
  2. `python3 py/anneal_summary.py logs/anneal_q7_k1[78]_*.out`

Whatever the outcome, annealing proves nothing about non-existence.

## 6. What is proved

**Theorem (FINITE-EXACT, computer-assisted).** Q_7 has no edge-multiset resolving set of size k for any
k = 1, 2, ..., 12. Hence **edim_m(Q_7) ≥ 13**. The known resolving 19-set, re-verified here by the pure-Python
checker, gives the upper bound: **13 ≤ edim_m(Q_7) ≤ 19**. A second, inequivalent resolving 19-set was found
by annealing (Section 5).

Every size k ≤ 12 is excluded separately, because resolvability is not monotone.

**Exactly what was checked for each k = 1..12** (all logged in `runs/d7_e0/validation_k<k>.json` and
`logs/validate_d7_e0.out`):
1. **Coverage.** For every a = ⌈k/2⌉..k with a representative file (all a ≤ 11), the completed chunks tile
   [0, nrep) exactly. The single remaining pair (k, a) = (12, 12) is excluded by Lemma 2; trivially, 12
   distinct vertices of Q_6 cannot agree in all 6 coordinates.
2. **Leaf count.** The C leaf count equals the independent DP count, for every (k, a) and in total
   (303,583,126,680 leaves for k ≤ 12).
3. **Result.** 0 resolving leaves in every chunk. Every leaf was either killed by an exactly computed collision
   or fully tested.
4. **Representatives.** The files used (Q_6, a ≤ 11) have the Burnside counts. The files for a ≤ 10 passed
   `repcheck`, with the full group for a ≤ 8 and the exact images-containing-0 mode for a = 9 and 10 (a = 11
   contributes 0 leaves and is unnecessary by Lemma 2).
5. **Same binary.** The same binary passed:
   * the Q_6 end-to-end reproduction (k ≤ 14: nothing; k = 15: exactly the 229 known orbits; both engines);
   * the Q_7 verification-mode samples (154,155,615 leaves, 0 discrepancies, including the k = 12 and k = 13
     configurations with b = 6);
   * the checker cross-check against pure Python (1,636 sets).
6. **Normal form.** Lemma 1 (proof in Section 2.2) reduces every k-set to the searched domain. This was
   confirmed empirically on 1,250 random sets.

**What the proof rests on.** Lemma 1 and Lemma 2 (hand proofs, Section 2); the completeness of the representative
lists R_a (their lengths equal the Burnside counts, and `repcheck` shows the listed sets pairwise inequivalent, item
4); the correctness of the C search code as supported by the checks above; and the usual trust in compiler and
hardware. The sha256 values of the production binaries, the deposited sources and the representative files are in
`logs/provenance_sha256.txt`.

**Not proved.**
* Anything for k ≥ 13. The k = 13 run is in progress and must not be counted until its COMPLETE line shows
  `leaves=1507949500245 found=0` and `validate.py` reports ok=True (Section 4.4).
* The annealing gives no non-existence information.

## 7. Files and reproduction

### 7.1 Source files needed (sizes in bytes)
| file | bytes | role |
|---|---:|---|
| `src/edimsearch.c` | 18,567 | the search: normal form B, engines 0/1, `-v` verification mode, `-x` check mode |
| `src/orbreps.c` | 7,406 | Aut(Q_n)-orbit representatives by orderly generation; `orbreps n amax outdir [selftest_count]` |
| `src/canon.c` | 2,319 | canonical form under the full Aut(Q_n), plus stabilizer order (used by validate, nf_test, anneal_summary) |
| `src/repcheck.c` | 3,738 | independent brute-force check of representative files |
| `src/anneal.c` | 5,634 | simulated annealing |
| `py/driver.py` | 10,035 | resumable, checkpointed 2-worker driver |
| `py/validate.py` | 4,608 | tiling, C leaves = DP, Python re-verification, orbits |
| `py/leafcount_dp.py` | 4,904 | independent DP leaf counts (numpy) |
| `py/pycheck.py` | 2,377 | pure-Python checker; `selftest` covers the Q_7 19-set and the Q_6 15-set |
| `py/burnside.py` | 1,995 | Burnside orbit counts |
| `py/xcheck.py` | 1,990 | C check mode vs pure Python |
| `py/nf_test.py` | 3,263 | empirical normal-form test |
| `py/make_tables.py` | 2,801 | the tables of this report, from the JSON files |
| `py/anneal_summary.py` | 2,157 | annealing summaries and re-verification |
| `runs/run_d7.sh`, `run_d6.sh`, `checks.sh`, `verify_d7.sh`, `post12.sh`, `chain2.sh`, `anneal.sh` | 401 / 514 / 1,323 / 1,264 / 2,389 / 598 / 819 | the run scripts used. Their `cd` line was changed to the package directory; `checks.sh` also reads its negative control from `data/` instead of `tmp/`, and one comment line of `anneal.sh` changed |

* All C sources compile with `gcc -O3 -march=native -Wall [-Wextra]`; `anneal` also needs `-lm`.
* Dependencies: python3 with numpy only.
* Superseded development versions (hot-list policy experiments, block-size benchmark) are not included.

### 7.2 Data: regenerated, not stored; evidence that is stored
* `data/q5/reps_n5_a*.bin`: 9.5 MB. Regenerate with `./src/orbreps 5 32 data/q5` (6 s).
* `data/q6/reps_n6_a0..10.bin`: 34 MB. Regenerate with `./src/orbreps 6 10 data/q6` (43 s).
* `data/q6/reps_n6_a11.bin`: 135.5 MB, `./src/orbreps 6 11 data/q6` (4 min). It is **not needed** for k ≤ 13.
  Without the file, the driver records the a = 11 exclusion by Lemma 2 (b ≤ 2 gives W_min(11) = 14 > 6b), and
  `validate.py` accepts it.
* `data/burnside_5.txt` and `data/burnside_6.txt` (stored as a reference): regenerate with
  `python3 py/burnside.py 5` and `python3 py/burnside.py 6`.
* `data/neg_n5_a3.bin`: the negative control of `runs/checks.sh` ({0,1,3}, {0,1,2} twice, {1,2,4}).
* Stored evidence:
  * `runs/d7_e0/`: `chunks.jsonl` (k <= 12), `summary_k1..12.json`, `validation_k1..12.json`;
  * `runs/d6_e0/` (k = 1..15) and `runs/d6_e1/` (k = 14, 15): `chunks.jsonl`, summaries, validations;
  * `logs/`: the k <= 12 driver log `d7_e0_k1-12.out`, and the validation, verification, check, DP, `orbreps`
    and annealing-calibration logs cited above. The Q_6 progress logs and the drivers' `log.txt` copies are not
    stored; the chunk records carry the same data.
* To re-validate the stored records, regenerate the representative files (including a = 11) and run
  `python3 py/validate.py 7 runs/d7_e0 $(seq 1 12)`. It rewrites the `validation_k*.json` files, so compare them
  with the stored ones.

### 7.3 Commands
Run from `04-computation/edge_multiset_dimension_q7_20261002/`.
```bash
# build (seconds) and representative files (about 50 s)
gcc -O3 -march=native -Wall -Wextra -o src/edimsearch src/edimsearch.c
gcc -O3 -march=native -Wall -o src/orbreps src/orbreps.c
gcc -O3 -march=native -Wall -o src/canon src/canon.c
gcc -O3 -march=native -Wall -Wextra -o src/repcheck src/repcheck.c
gcc -O3 -march=native -Wall -o src/anneal src/anneal.c -lm
mkdir -p data/q5 data/q6 logs runs
./src/orbreps 5 32 data/q5 > logs/orbreps_q5.log          # 6 s
./src/orbreps 6 10 data/q6 > logs/orbreps_q6.log          # 43 s
python3 py/burnside.py 5; python3 py/burnside.py 6        # compare with the 'reps=' counts in the orbreps logs
```

**Full Q_7 reproduction of k ≤ 12.** About 73 CPU min, about 38 min wall on 2 cores, of which k = 12 is
about 34 min. Use a fresh run directory: `runs/d7_e0` holds the stored records, and the driver skips chunks it
finds there.
```bash
for k in $(seq 1 12); do python3 py/driver.py 7 $k --workers 2 --target 120 --engine 0 --run runs/repro_d7; done
python3 py/validate.py 7 runs/repro_d7 $(seq 1 12)        # ~2 min; expect 12 lines "ok=True" with the totals of 4.1
```
Without `data/q6/reps_n6_a11.bin` the driver records the a = 11 pairs as excluded by Lemma 2 instead of
searching them (they have 0 leaves either way).

**Fast subset** (what the growth runner's `--q7` option runs). About 3-5 min wall on 2 cores in total,
including the builds and `orbreps` above:
```bash
for k in $(seq 1 10); do python3 py/driver.py 7 $k --workers 2 --target 60 --engine 0 --run runs/fast_d7; done   # ~55 CPU s
python3 py/validate.py 7 runs/fast_d7 $(seq 1 10)          # DP cross-check: ok=True; k=10 leaves C=4141968201 DP=4141968201
for k in $(seq 1 13); do python3 py/driver.py 6 $k --workers 2 --target 60 --engine 0 --run runs/fast_d6; done   # ~50 CPU s
python3 py/validate.py 6 runs/fast_d6 $(seq 1 13)          # ok=True; k=13 leaves C=DP=1234172511, 0 resolving
python3 py/pycheck.py selftest                             # 19-set and Q6 15-set resolving
python3 py/xcheck.py 100 7                                 # "XCHECK ALL AGREE"
```

**Strongest end-to-end check.** About 10 more CPU min, about 5 min wall. It adds Q_6 k = 14 and 15:
```bash
for k in 14 15; do python3 py/driver.py 6 $k --workers 2 --target 60 --engine 0 --run runs/fast_d6; done
python3 py/validate.py 6 runs/fast_d6 14 15                # k=15: resolving_leaves=585 orbits=229 stab=[1]
```

**Further optional checks.**
* `runs/checks.sh`: pycheck, xcheck, repcheck, nf_test; about 9 min on 1 core.
* `runs/verify_d7.sh big1|big2|small`: 11, 11 and 1 min.
* `./src/orbreps 4 16 data 3000`, `./src/orbreps 5 32 data 3000` and `./src/orbreps 6 6 data 1500` (self-test mode
  writes no files): the `is_canon` self-tests of Sections 2.4 and 3.1, logged in `logs/iscanon_selftest.out`.

The fast subset was executed through the growth runner's `--q7` option, which builds the sources and runs
these commands in a temporary directory (161 s with 2 of 4 shared cores): Q_6, k <= 13, 1,859,531,279 leaves,
and Q_7, k <= 10, 4,533,173,692 leaves; C = DP for every (k, a), nothing found. The end-to-end block (Q_6,
k = 14, 15) was not re-executed under a fresh directory; it is the same command that produced `runs/d6_e0`, so
the expected numbers are those of Section 3.4.

**Notes on `driver.py`.**
* It refuses to start if another driver holds `runs/<run>/driver_k<K>.pid`.
* It never reuses records of a different (D, K).
* Re-running a finished K only rewrites its summary.

## 8. Points of uncertainty
1. **Annealing is incomplete.** Only the calibration was run before this note was first committed; the 18/17
   annealing was held back to keep the CPU for k = 13. A 90-minute run at size 18 on one spare core started at
   01:31 UTC on 2026-10-02; its outcome will be recorded here.
2. **One production run per Q_7 size.** Each k was run once, with engine 0.
   * Independent support: the DP leaf counts (enumeration), 154M leaves in verification mode (0.05% of the
     3.04e11 leaves), and the end-to-end Q_6 reproduction with both engines.
   * The planned engine-1 re-run of Q_7 for k ≤ 11 was cancelled to keep CPU for k = 13.
3. **Shared hash test.** The "full test" hash function (open addressing with exact key comparison) is shared
   by both engines, the search and the verification mode's ground truth.
   * It is cross-checked independently by the sort-based test in `-x` mode (agreement is asserted on every
     set), by the Python checker (1,636 sets), and by the Q_6 reproduction.
   * A fault that shows up only on Q_7 key patterns outside the verified samples cannot be excluded
     completely by testing.
4. **Hand proofs.** The normal form (Lemma 1) and the weight lemma (Lemma 2) are proved by hand. For k ≤ 12,
   Lemma 2 is used only for (12, 12), which is trivial. Lemma 1 was also tested empirically on 1,250 sets.
5. **Fast-mode repcheck.** It was used for Q_6 a = 9 and 10 and rests on the elementary argument that only
   images containing 0 can be lex-smaller. The full-group mode was used for every other file.
6. **Timing.** The k = 13 completion estimate (03:10-03:40 UTC) depends on the load from other jobs.
   Observed speeds were 12-15 ns per leaf.
7. **Earlier verify sample.** The earlier k = 11 verify-mode sample (15.2M leaves per engine)
   has no log file.
8. **Usual caveats.** Trust in compiler, hardware and memory. No ECC information is available.
