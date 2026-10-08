# Mersenne-line barriers: the t = 23 fan escapes, the unpaid progression is reached, and the stopping-time invariant confines deletion descent

Session opus-2026-10-07-S21.
* **Owner directive:** take codex-grounding's status as the brief: t = 23 reached; giant children ungrounded; "the next decisive target is a terminating, source-preserving rule for one of those actual children". Search past repo work for connections while testing hypotheses. Three PDFs were attached.
* **Scripts.** All are in `04-computation/experiments/`, each with a `.out`, and each prints ALL CHECKS PASSED.

| script | sections | runtime |
|---|---|---|
| `mersenne_line_barriers_20261007.py` | (A) chain vs orbits, (A2) no-meeting lemma, (B) reset rule, (C) t = 23, (D) fan, (E) progression, (F) tails, (M) stopping-time invariant | 25 s, 12 cores |
| `mersenne_fan_escape_20261007.py` | (P) fast parity vectors, (G) escape, (H) residue check, (I) contiguous edge | 2.5 min, 2 GB |
| `mersenne_fan_block_20261007.py` | (J) full block D ≤ 4000, (K) its clusters and the level-2 merge, (N) odd-step counts | 2.5 min |
| `mersenne_line_coalescence_20261007.py` | (L) local hardness, (Q)/(Q2) deep progression exits, (S) coalescence profiles, (H2) exact horizon clusters | 6 min |
| `mersenne_fan_level2_20261007.py` | level 2 to `2^29` | 21 min |

* **Canon:** THM-4605, HYP-9243.
* **Audit:** an independent adversarial audit followed checkpoint `023a49bcf`; its corrections are recorded in MISTAKE-588.

**Status.**
* PROVED: certificate soundness, the no-meeting lemma and the fast parity vectors (§1); the stopping-time invariant and its corollaries (§4.1); density-one exits for each D, a corollary of THM-4581 (3) (§4.4).
* KNOWN: the reset rule (THM-4601 (i), THM-4556); the chain itself (THM-4581 (6d), used here with the source as reference); exponent 1/2 of the no-exit tail for finite lag sets (THM-4593).
* FINITE-EXACT: every certificate and negative bound.
* NUMERICAL: tails and coalescence profiles.
* HEURISTIC: the mechanism of diffusive coalescence (HYP-9243); the O(√E) cluster span.
* ANALOGY: the PDFs and the two-walls reading.
* Collatz OPEN.

## 0. Results

1. **A certificate engine.**
   * A Mersenne exit `M_E ⇝ M_(E−D)` is the absorption of one pair chain: THM-4581's Mersenne-switch chain (6d), run with the source's own orbit as reference.
   * Its coins are the Terras parities of `x = 2·3^(E−1) − 1`, which need only `3^(E−1) mod 2^W`.
   * It reproduces codex-grounding's t = 23 result exactly: deletions `{3, 4, 7, …, 28}`, all at Terras time 6553.
   * Coins come from divide-and-conquer parity vectors (python-flint); the `2^25` coins of the fan bottom take about 12 s.
   * The recursion composes the affine maps of the two halves of a parity block. It is the recursion Elsenhans used to verify random numbers of up to 10^10 decimal digits ([arXiv:2502.16743](https://arxiv.org/abs/2502.16743)). The affine form and the parity-vector correspondence are due to Terras (1976), Everett (1977) and Lagarias (1985).
2. **The unpaid progression `t ≡ 13847 mod 27648`.**
   * 1999 of its first 2000 members have certified exits with fixed source.
   * The remaining member `t = 3525143` has no exit with D ≤ 256 at any Terras time ≤ `2^26 − 8`.
3. **The fan escapes.**
   * The source and its 24 children form a closed fan: each child's exits are exactly the fan members below it.
   * The bottom `M_99708993677` has no exit with D ≤ 256 before Terras time 19,000,765.
   * At that time the chains of 1916 deletions are absorbed together: D = 1–1910 and D = 1929–1932, 1935, 1936. The absorbed group is not an interval.
   * So `M_99708993677 ⇝ M_(99708993677 − D)` holds for those 1916 values. The deepest is `M_99708991741`.
   * **1949 exponents now share one certified future:**
     * the parent `M_99708993713`;
     * the seven exponents between it and the source;
     * the source `M_99708993705`;
     * the 24 children;
     * the 1916 exponents of the absorbed group.
4. **Level 2 lies deeper than `5.4·10⁸`.**
   * The other 2084 deletions with D ≤ 4000 form two interleaved clusters, with D ranges 1911–2902 and 2895–4000. They merge at step 4,474,989.
   * The merged cluster is not absorbed at any Terras time ≤ `2^29 − 9`.
5. **The stopping-time invariant (PROVED, §4.1).**
   * A certificate with merge value ≥ 3 preserves two quantities: `σ_T(x_K)`, the Terras time from `x_K` to 1, and the odd-step count of `M_K`. Conversely, equal invariants force a certificate.
   * So, among exponents whose orbits reach 1, deletion clusters are exactly the level sets of `(σ_T(M_K) − K, odd count of M_K)`. This is HYP-9242's partner criterion.
   * HYP-9242's corrected status (mac-mini, after its audit C, pushed concurrently) states the same equivalence for every K whose orbit reaches 1, with the merge at a value ≥ 8. Proposition 4.1 (a) also covers orbits that never reach 1, which Corollary 4.2 needs.
   * `σ_T(x_E) > (E − 1)·log₂3`. Hence **no route of certificates with merge values ≥ 3 from any `E > 61821` reaches any `K ≤ 12800`.**
   * For the supplied children, every deletion route ends at an exponent whose remaining orbit is still longer than `1.58·10^11` steps.
   * In words: each deleted bit removes exactly one rising step from `σ_T(M_K)` and leaves `σ_T(x_K)` unchanged, while an orbit costs at least `1 + log₂3 = 2.585` Terras steps per bit of the exponent, about 8.64 typically. The typical deficit is 7.64 steps per bit, the residual-horizon constant `6.95212·ln 3`.
6. **The answer to the brief.**
   * No source-preserving *deletion* rule terminates for the supplied children.
   * A terminating route must carry total time offset `σ(M_E) − σ(M_K)`. Deletions carry exactly D each: the free offset of the rising phase `3^j 2^(E−j) − 1`.
   * The remaining roughly `7.6·10^11` Terras steps (a typical value) must come from certificates whose offset is paid by computing the orbit, or by a closed form not yet known.
   * Constructively, **one orbit grounds the fan.**
     * `σ_T` is constant on the 1949-exponent class. So computing the orbit of any one member to 1 grounds all 1949.
     * That needs a starting number of `1.58·10^11` bits (`4.76·10^10` decimal digits). A typical orbit of that size has about `7.6·10^11` Terras steps; the proved lower bound is `1.58·10^11`.
     * This is 4.8 times the largest number Elsenhans verified (10^10 digits, about 14 h on one core of an i7-12700).
     * By his empirical `O(N^(1+ε))` scaling it is of the order of three days on one core. Memory for several copies of a 20 GB integer is our estimate; the paper reports no memory figures.
7. **Diffusive coalescence (HYP-9243, reworked around cluster size).**
   * On six fresh odd E near 10^11, the source's cluster below it has upper-median size N(T) = 28, 138 and 332 at T = `2^12`, `2^17` and `2^20`. Median `N(T)/√T` lies in 0.23–0.52.
   * At the exact horizon (K ≤ 12800) the per-cluster mean size is consistent with `√K`: 2/(mean size)·`√K` = 0.90–1.13 for K ≥ 800.
   * This is HYP-9242's orphan fraction, which equals 2/(mean cluster size) per odd K (§4.2). HYP-9242 now fits its exponent as 0.59, with 95% CI [0.42, 0.75].
   * The size-biased median below K grows a little faster in this range (0.25 → 0.44·`√T`).
   * Claimed: the exponent 1/2 at both scales. Not claimed: equal constants.
   * THM-4593 already pins the exponent 1/2 for the probability that a source has no exit at all, for finite lag sets. HYP-9243 concerns the growth of the cluster itself.
8. **The reset rule (KNOWN).**
   * For even E, `M_E ⇝ M_(E−1)` at Terras time 3: with `w = 3^(E−2) ≡ 1 mod 8`, both points reach `(9w − 1)/4`.
   * This is THM-4601 (i), in its lag-one form. In OEIS [A075485](https://oeis.org/A075485) and A193688 (A. Herrera, 2023) it appears as a(2n) = a(2n − 1) + 1.
   * Every hard exponent is odd.

## 1. The certificate

**Setting.**
* `T` is the Terras map and `M_E = 2^E − 1`, with `T^j(M_E) = 3^j 2^(E−j) − 1` for `j ≤ E`.
* So the source sits at `x = 2·3^(E−1) − 1` after `E − 1` steps, and the child `M_(E−D)` at `y = 2·3^(E−D−1) − 1`, with `y + 1 = 3^(−D)(x + 1)`.
* Write `y_n = 3^(k_n) x_n + e_n` and run THM-4581's pair chain with the source as reference:
  * coin = parity of `x_n`;
  * integer state `(k, a)` with `a = 3^max(0,−k) e`, so `3^(−k) y_n = x_n + a` for `k < 0`;
  * start `(−D, 1 − 3^D)`;
  * the four updates are the table of THM-4605.
* THM-4581 (6d) runs the same chain with the child as reference, from `(D, 3^D − 1)`.

**Theorem 1.1 (soundness; PROVED).** If the chain is at `(0, 0)` at time n, then `T^n(x) = T^n(y)`. So `M_E` and `M_(E−D)` have a common future, and `M_E` reaches 1 iff `M_(E−D)` does.

*Proof.*
* The relation `y_i = 3^(k_i) x_i + e_i` holds exactly at every step, by the four cases of the table (as in THM-4581). The audit re-derived all four updates for both signs of k.
* `(k, e) = (0, 0)` is equality, and `(0, 0)` is a fixed point of the update.
* The coins for `i < n` depend only on `x mod 2^(i+1)`.
* Relation absorption does not claim trivial-cycle value coincidences (S19's lesson, MISTAKE-580). ∎

**Lemma 1.2 (no meeting; PROVED).** Let `N ≥ 0`, `D ≥ 1` and `E ≥ (1 + log₃4)·N + D + 2`.
* If the chain is not absorbed at any time ≤ N, then `T^m(x) ≠ T^(m')(y)` for all `m, m' ≤ N`.
* If it is first absorbed at `n ≤ N`, every meeting with `m, m' ≤ N` has `m = m' ≥ n`.

*Proof.*
* Write `T^m(x) = (3^A x + C)/2^m` with `0 ≤ C < 3^A 2^m`, and `T^(m')(y) = (3^B y + C')/2^(m')` likewise, where `A, B ≤ N` are odd-step counts.
* Since `3^D (y + 1) = x + 1`, a meeting gives `x·(2^(m') 3^(A+D) − 2^m 3^B) = R` with `R = 2^m 3^B (1 − 3^D) + 2^m 3^D C' − 2^(m') 3^D C`. So `|R| < 3^(N+D+1) 4^N`.
* We have `x > 3^(E−1) ≥ 3^(N+D+1) 4^N`, so the bracket vanishes. Hence `m = m'` and `B = A + D`, i.e. the level `k_m = B − A − D` is 0.
* Then `e_m = y_m − x_m = 0`, so the chain is at `(0, 0)` at time m. ∎

**Consequences.**
* The hypothesis holds with a wide margin in every negative, with N the last checked time:
  * the fan bottom (E ≈ 10^11): N = `2^29 − 9` for level 2 (D ≤ 4000), 19,000,764 for D ≤ 256 before the escape, and 16375 for the fan closure (D ≤ 128);
  * `t = 3525143` (E ≈ 1.5·10^16): N = `2^26 − 8`, D ≤ 256.
* So each negative below means more than "no chain absorption": the first N + 1 points of the two orbits share no value.
* (A2) checks the lemma on 800 cases with E in [1500, 1520), D ≤ 40 and N = 600: 556 unabsorbed chains with disjoint orbit segments, and 244 absorbed chains meeting only at equal times, first at the absorption time.

**Convention.** "Within W" means at every Terras time ≤ W − 9, since chains run on W − 8 coins. Where the final state is also checked (`deep_exit` in the coalescence script), it means ≤ W − 8.

**Fast parity vectors (PROVED).**
* For x known mod `2^n`, `par(x, n)` returns `(V, a, c)`: V is the parity vector, a its weight, and `T^n(x') = (3^a x' + c)/2^n` for all `x' ≡ x mod 2^n`.
* The base case appends one step; an odd step updates `c ↦ 3c + 2^i`.
* The two halves combine by composing their affine maps: `x_mid = (3^(a1) x + c1) >> n1`, `a = a1 + a2`, `c = 3^(a2) c1 + 2^(n1) c2`.
* The cost is `O(M(n) log n)`.
* The recursion is Elsenhans' (see §5 for prior art).
* (P) compares it with naive iteration, including `n = 2^17` steps on the fan bottom.

**Validation.**
* 2400 large-E cases agree with direct orbits (202 merges).
* All 306 small-E disagreements are cycle coincidences at value ≤ 2.
* A numba prototype was used for speed only. Its int64 overflow guard abandons chains rather than claiming them, and no certificate uses it.

## 2. The supplied residual and its fan (FINITE-EXACT)

**The source.** `E = 99708993705` (t = 23 in `E(t) = 924745897 + 2^32 t`) has exits exactly `D ∈ {3, 4, 7, …, 28}`, all at Terras time 6553 (W = 7200, D ≤ 40). These are codex-grounding's 24 children.

**The original parent** `M_99708993713` (W = 8000, D ≤ 64) has these exits:
* D = 1–4 at Terras time 19;
* D = 7 and 8 at 35, where D = 8 is the source (codex-grounding's 8-bit receipt);
* D = 5, 6, 11, 12 and 15–36 at 6553. Here D = 11, 12, 15–36 are the 24 children, and D = 36 is the fan bottom, matching codex-grounding's composed 36-bit deletion.

**The fan is closed.**
* Within W = 16384 and D ≤ 128, each child's exits are exactly the fan members below it, and the bottom `E − 28` has none.
* The internal merge times are:
  * 3, for twelve exits: the (even, odd) reset pairs;
  * then 9, 25, 57, 76, 111, 123, 131, 276, 362, 512 and 644 (counts in the `.out`);
  * everything merges with the source at 6553.

**The bottom is deep** ((G)).
* No exit with D ≤ 256 occurs before Terras time 19,000,765.
* (G) returns the first absorption among the 256 chains, run with state collapse. Identical states have identical futures, so this first absorption is the global minimum.
* The orbit is ordinary (odd density 0.50001 over `2^25` steps). The barrier is the relation level, not a rising phase.

**The full block** ((J), (K): all 4000 chains with D ≤ 4000, with collapse).
* Merges of the marked chains:
  * D = 4000 into 2895 at step 773,555;
  * 2902 into 1911 at 1,338,070;
  * 1910 into 1 at 2,741,933, the last merge of the D = 1 group;
  * 2895 into 1911 at 4,474,989.
* At step 2,741,940:
  * the D = 1 group has its final 1916 members, `{1, …, 1910} ∪ {1929, 1930, 1931, 1932, 1935, 1936}`;
  * the other chains form two clusters, with D ranges 1911–2902 (982 members) and 2895–4000 (1102). The ranges interleave.
* Levels of the D = 1 group:
  * −1131 at 2,741,940;
  * −2335 at 4,474,989;
  * −2608 at `2^24`;
  * 0 at 19,000,765.
  * Its collapsed chain (D ≤ 256, from step 2,741,933 on) reaches minimum level −3867 ((G)).

**The escape.**
* The D = 1 group is absorbed at **Terras time 19,000,765** with all 1916 members.
* So `M_99708993677 ⇝ M_(99708993677 − D)` for each of them; the deepest is `M_99708991741` (D = 1936).
* The contiguous extent is 1910: (I) shows that D = 1911 stays separate.

**Independent check** ((H), (N)).
* `T^n(x) ≡ T^n(y_D) mod 2^(2^20)` at n = 19,000,765 for D = 1 and D = 256, computed from the affine maps of the fast parity vectors.
* The odd-step counts are 9500295, 9500296 and 9500551 for x, y_1 and y_256. The differences are 1 and 256: level 0.
* The residues differ at n − 1.

**Level 2** ((K) and the level-2 script).
* The remaining 2084 chains form one cluster from step 4,474,989 on. Its level is −6133 at 19,000,765.
* It is not absorbed at any Terras time ≤ `2^29 − 9`. Over that run the level of the D = 1911 chain, counted from step 0, stays in [−31536, −1910].
* The maximum, −1910, is one above the start of the D = 1911 chain and is attained before 19,000,765.
* That the cluster never climbs back later is unremarkable (HEURISTIC). A symmetric walk from −6133, at about one flip per two Terras steps, stays below −1910 for the remaining `5.2·10^8` steps with probability about 0.2 (reflection principle).
* Chains with D > 4000 were not run.

**What this does for the owner's obligation.**
* The 24 child ports are no longer trapped. They have certified strictly smaller descendants down to `M_99708991741`: a source-preserving child rule, with rank the exponent.
* They are not grounded, and by §4.1 no deletion route can ground them.
* Grounding needs either one orbit computation (cost in §0.6) or a non-deletion child whose time offset is not paid by the rising phase.

## 3. The unpaid progression `t ≡ 13847 mod 27648` (FINITE-EXACT)

* The source `M_E(t)`, with `E(t) = 924745897 + 2^32 t`, stays fixed. Each certificate deletes D bits to reach a smaller Mersenne number.
* Examples (largest D among the cheapest exits):
  * t = 13847: D = 6 at Terras time 273;
  * t = 41495: D = 36 at 957;
  * t = 69143: D = 4 at 73;
  * t = 96791: D = 12 at 129.
* **First 2000 members.**
  * 1970 have certificates within W ≤ 16384 (D ≤ 64). By first successful W: 256 (1258), 1024 (489), 4096 (169), 16384 (54). The median cheapest merge time is 157, the maximum 16151.
  * 29 more are certified by the deep search (D ≤ 256, W ≤ `2^22`). Their merge times run from 16615 to 2,520,494, and absorbed groups contain up to 256 deletions (the cap).
  * The remaining member `t = 3525143` (`E = 15140374823469225`) has no exit with D ≤ 256 at any Terras time ≤ `2^26 − 8`.
* Each certificate covers only a 2-adic cell of mass about `2^(−(cost−34))` (codex-grounding's anchored phase lemma). These are pointwise exits, the requested "explicit member", not new density.

## 4. The stopping-time invariant, horizon clusters and diffusive coalescence

### 4.1 The invariant (PROVED)

For K ≥ 2:
* `s(K) = σ_T(x_K)` is the Terras time from `x_K = 2·3^(K−1) − 1` to 1 (∞ if the orbit never reaches 1);
* when finite, `o(K)` is the number of odd steps of `M_K`'s orbit to 1;
* thus `σ_T(M_K) = K − 1 + s(K)`.

**Proposition 4.1.** Let `1 ≤ D ≤ E − 2`.
* (a) If the (E, D) chain is absorbed at time n with merge value `v = T^n(x_E) ≥ 3`, then `s(E) = s(E − D)`, and `o(E) = o(E − D)` when finite.
* (b) Conversely, if `s(E) = s(E − D) < ∞` and `o(E) = o(E − D)`, the chain is absorbed by time `s(E) − 3`, with merge value ≥ 8. HYP-9242's corrected status states the same.
* (c) An absorption at a value ≤ 2 occurs at value 2, with one orbit already at 1 and the other one step from 1. It happens iff `s(E) ≠ s(E − D)` and `2o − s` agree. No pair with K ≤ 12800 qualifies ((M)).

*Proof.*
* The level k counts odd steps of y minus odd steps of x, minus D. So after both orbits reach 1 at a common time, the level equals `o(E − D) − o(E)`.
* (a) From time n on the two orbits coincide. Neither reached 1 earlier, since its value at time n would then lie in {1, 2}. So both reach 1 at the same time, or neither does. The level is 0 at absorption and stays 0.
* (b) Both orbits are at 8 at time `s − 3`. The only Terras preimage of 1 is 2, the only preimage of 2 other than 1 is 4, and the only preimage of 4 is 8.
  * The level at time s is `o(E − D) − o(E) = 0`. The last three steps are even, so the level is also 0 at time `s − 3`.
  * There `e = y − x = 0`.
  * In (a), "merge value ≥ 3" is the same as "≥ 4": `x_K ≡ 2 mod 3`, and the Terras preimages of a multiple of 3 are multiples of 3, so no orbit of an `x_K` visits a multiple of 3.
* (c) If a first absorption at time n had value 1, both orbits were at 2 at time `n − 1`. Equal values with equal parities leave the level unchanged, so the chain was already at (0, 0). The same holds at value 2 with both predecessors at 4. So the predecessors are 1 and 4.
  * On the cycle each turn adds two steps and one odd step. So once both orbits are on the cycle in phase, the level is `[(2o − s)(E − D) − (2o − s)(E)]/2`, which vanishes iff `2o − s` agree. ∎

**Corollary 4.2 (deletion routes; PROVED).**
* Along any route of certificates with merge values ≥ 3, `(s(K), o(K))` is constant, so a route from E can end at K only if `s(K) = s(E)`.
* For a giant E this excludes every K whose orbit is known, unless some certificate on the route merges at the trivial cycle. Such a certificate is itself a complete orbit computation.
* `s(E) ≥ log₂ x_E > (E − 1)·log₂3`, because one Terras step at most halves.
* For K ≤ 12800 the maximum of s(K) is 97982 (mac-mini's exact table, cross-checked for K ≤ 2000 in (M)).
* Hence no route from any E > 61821 reaches any K ≤ 12800. For K ≤ 2000 the threshold is E > 9828.
* For the supplied exponents (E ≈ 10^11), every deletion route ends at an exponent with `s(K) = s(E) > 1.58·10^11`.
* (M) confirms (a) and (b) on all 919 absorptions with 4 ≤ E ≤ 160. None of them merges at the cycle.

**Offset accounting (PROVED; elementary).**
* Any certified merge `T^(t_s)(s) = T^(t_c)(c)` at a value ≥ 3 gives `σ(c) = σ(s) − (t_s − t_c)`. Offsets add along routes.
* In the Mersenne clock a deletion carries offset exactly D. Its chain part has offset 0, since it compares equal times n. So a deletion route from E to K carries offset `E − K`: one Terras step per deleted bit.
* An orbit costs at least `1 + log₂3 = 2.585` steps per bit of the exponent: `E − 1` rising odd steps plus at least `(E − 1)·log₂3` to come down from x.
* Typically it costs about 8.64 steps per bit:
  * the residual horizon is the classical `6.95212·ln 3 = 7.6377` (Lagarias; Kontorovich–Lagarias);
  * (M) measures 7.6438 over 1000 ≤ K ≤ 12800;
  * Ren's `2^100000 − 1` has 863,323 Terras steps and 481,603 odd steps, recomputed in (M): `σ_T/K = 8.633`.
* The deficit must be paid by computing orbit steps, or by a closed form other than the rising phase.

**Heuristic span (HEURISTIC).**
* `s(K)` fluctuates around `7.64K` by about `10.6·√K`: this is the standard deviation of `(s(K) − cK)/√K` over 1000 ≤ K ≤ 12800, from (M).
* Within a cluster s is constant, against a trend of 7.64 per unit K. That confines a cluster to a window of a few `√K` exponents.
* Measured spans average `2.6·√K` and reach `6.1·√K` (83 clusters with bottoms ≥ 1000).
* For E ≈ 10^11 this suggests that deletion routes from the fan never leave a window of order 10^6 exponents.

### 4.2 Exact horizon clusters (FINITE-EXACT; (H2), mac-mini's table)

* For K ≤ 12800 every orbit reaches 1. By Proposition 4.1 the deletion clusters are the level sets of `(σ_T(M_K) − K, o(K))`.
* There are 136 clusters among K = 2..12800.
* Every cluster bottom with K ≥ 4 is odd. This is PROVED by the reset rule: an even K ≥ 4 merges with K − 1 at time 3, at a value ≥ 20. (H2) confirms it.
* Consecutive members K, K + 1 of one cluster have standard step counts differing by exactly 1. In OEIS A075485 this is the identity a(2n) = a(2n − 1) + 1 together with runs of +1. For example, K = 33–42 form one cluster.

| window | odd K | median `T = σ_T(x_K)` | median `S/√T` | median (cluster size below K)/`√T` | orphan fraction · `√K` | `2/(mean size)` · `√K` |
|---|---|---|---|---|---|---|
| [400, 800) | 200 | 4652 | 0.078 | 0.245 | 1.546 | 1.487 |
| [800, 1600) | 400 | 8942 | 0.098 | 0.358 | 1.093 | 1.024 |
| [1600, 3200) | 800 | 17447 | 0.100 | 0.395 | 1.011 | 1.053 |
| [3200, 6400) | 1600 | 36671 | 0.219 | 0.341 | 1.135 | 1.134 |
| [6400, 12800] | 3200 | 73529 | 0.271 | 0.443 | 0.892 | 0.898 |

* The orphan fraction per odd K (HYP-9242) is 2/(mean cluster size), up to window edges. This is an identity, not an independent agreement. The checkpoint's comparison of a predicted band with HYP-9242's constant was a bookkeeping artifact (MISTAKE-588).
* The contiguous extent S is the wrong statistic: its median grows faster than `√T` at the horizon, with `S/√T` rising from 0.08 to 0.27. The cluster size below K stays at 0.25–0.44·`√T`.

### 4.3 Coalescence at E ≈ 10^11 (NUMERICAL; (S), (F), (L))

**Profiles.**
* Setting: six fresh odd E (seed 77), deletions D ≤ 2048.
* The table gives the contiguous absorbed extent S(T) and the cluster size N(T), the number of D ≤ 2048 absorbed by T.

| E | S, N at T = 4096 | at `2^17` | at `2^20` |
|---|---|---|---|
| 47791912765 | 0, 2 | 0, 28 | 0, 28 |
| 37169560031 | 36, 42 | 78, 84 | 78, 84 |
| 36802715565 | 2, 2 | 312, 332 | 312, 332 |
| 49148895517 | 8, 20 | 50, 72 | 304, 358 |
| 42690426595 | 10, 28 | 212, 214 | 212, 214 |
| 83639239217 | 68, 78 | 136, 138 | 1796, 1798 |
| upper median | 10, 28 | 136, 138 | 304, 332 |

* Median `N(T)/√T` ranges over 0.23–0.52 between `2^12` and `2^20`, and N grows in heavy-tailed jumps.
* 1–4 unabsorbed groups remain among the 2048 deletions at `2^20`.
* The t = 23 fan bottom is one heavy-tailed jump of this kind: N(T) = 0 until T = 19,000,765, then N = 1916 (D ≤ 4000) and S = 1910, against `√T = 4359` (ratio 0.44). A single selected source is not evidence for the scaling.

**Exit-depth tails** (600 odd E in [10^10, 10^11]; even E always exit at time 3).
* `P(no exit with D ≤ 64 by W)` = 0.540, 0.383, 0.258, 0.160, 0.092, 0.047 at W = 128, 256, …, 4096.
* The local exponents are 0.50, 0.57, 0.69, 0.80 and 0.97, about 0.70 overall. The last two rest on 96, 55 and 28 survivors of 600, with standard errors of about 0.13 and 0.19, so they do not show a steepening.
* Up to the cap this is THM-4593's statistic: even-D exits pair with D − 1 exits (THM-4556), and THM-4593 treats the Mersenne lag set D ≤ 61.
  * There the local slopes are 0.70, 0.69, 0.64, 0.59, 0.51 and 0.48 on [4·10², 10⁶]. The 0.69 is a coalescence transient.
  * The exponent is 1/2: PROVED from below by witness classes, and from above at sketch level.
  * Our numbers are consistent with that transient.
* mac-mini's depth spectrum for D ∈ {1, 3} on the Mersenne line (`runcompress_orphans_cayley_20261007.md`) shows the single-chain tail directly: `P(depth > s)·√s ≈ 10` for s in [10³, 10⁴].
* Locally below the t = 23 source, 581 of 20000 exponents (2.9%) have no exit with D ≤ 64 within 4096, all odd (reset rule).

### 4.4 The sources are 2-adically Haar

* The coins `β_0, …, β_(n−1)` depend on E only through E mod `2^(n−3)`, the order of 3 mod `2^(n−1)`. This is the 2-adic periodicity of THM-4556 (iv).
* For odd `E = 2b + 1`, `x_E = 2·9^b − 1`, and `b ↦ 9^b` is a measure isomorphism of `Z_2` onto `1 + 8Z_2` (the device of THM-4581 (6d)'s proof). So `x_E` is Haar on `1 + 16Z_2`, and on `5 + 16Z_2` for even E.
* **Corollary (PROVED, from THM-4581 (3)).** For every D ≥ 1, the exponents E with a certified exit `M_E ⇝ M_(E−D)` have natural density 1.
  * The (E, D) chain is a THM-4581 relation (`k_0 = −D`, `e_0 = (1 − 3^D)/3^D ∈ Z[1/3]`), absorbed almost surely on that positive-measure set.
  * Absorption by time n depends only on E mod `2^(n−3)`.
  * (6d) states this for D = 1.
* So depth statistics over exponents are Haar statistics in the limit. The 600 exponents of (F) are a pseudo-random sample, not an exact Monte Carlo: equidistribution of E < `2^37` at depths beyond about `log₂E` is not proved (THM-4601 (vi), HYP-9242).
* Depth is unbounded on the line: by THM-4601 (ii), no chain with D ≤ 8000 of an odd E is absorbed at any Terras time ≤ `3 + v₂(E − 1)` (≤ 5 for the fan bottom).
* A uniform depth bound would therefore be a 2-adic statement about E, not a real-place problem. The checkpoint's comparison with Mahler's 3/2 problem is withdrawn.

### 4.5 HYP-9243

* The hypothesis is now stated for the cluster size, not the contiguous extent. It asserts exponent 1/2 for two statistics:
  * the size-biased cluster `N_E(T)` below a source, from `T = 2^12` at E ≈ 10^11 up to the horizon;
  * the per-cluster mean size at the horizon, where the clusters are the exact level sets of §4.2. This one is HYP-9242's orphan fraction; its fitted exponent, 0.59 with 95% CI [0.42, 0.75], includes 1/2.
* THM-4593 pins the exponent of the no-exit probability; HYP-9243 concerns the growth of the cluster.
* Its mechanism (a symmetric level walk between clusters) is an assumption for integers: it is PROVED only in the Haar model (THM-4581 (1b)).
* The agreement with HYP-9242 is claimed only in the exponent −1/2.
* The former "consequence" is now Corollary 4.2 and is no longer part of the hypothesis.

## 5. Connections

**THM-4581** (mac-mini, Haar coalescence).
* Our certificates are its (6d) absorption events, with the source as reference.
* By §4.4 their law over exponents is its Haar law.
* Its rate (4), `T^(−1/2)` at sketch level, is the single-chain tail.
* Its one-big-jump heuristic (7), `P(no merge by T) ≈ (|k_0| + E[J])·√(4/(πT))`, gives meeting times of order `D²` for a start at level −D.

**HYP-9242** (orphans).
* Its partner criterion is the invariant of Proposition 4.1.
* Its corrected status (pushed concurrently, after its audit C) proves the equivalence with an equal-time merge for every K whose orbit reaches 1, at a value ≥ 8. That is Proposition 4.1 (b).
* It now fits the Mersenne orphan exponent as α̂ = 0.59 (95% CI [0.42, 0.75]).
* The orphan fraction is 2/(mean horizon cluster size) (§4.2). So HYP-9243's per-cluster exponent 1/2 is one of the readings HYP-9242 leaves open.

**THM-4593** (mac-mini, partner coalescence).
* It treats the no-exit probability for finite lag sets, including Mersenne lags D ≤ 61: the exponent is 1/2 and 0.69 is a coalescence transient. With all translation lags the decay is exponential.
* Its "one coin per step" (1): all relations to a common reference move in the same direction when they move. So the relative level of two clusters changes only when exactly one of them flips.
* Our (F) is its statistic up to the cap.

**THM-4601** (two anchors). (i) is the reset rule. By the +1 barrier (ii), no chain with D ≤ 8000 is absorbed at Terras times ≤ `3 + v₂(K − 1)` (≤ 5 for the fan bottom). The deep barriers are statistical, not that one.

**THM-4556.** Level sets of σ are unions of reset pairs (ii); the 2-adic periodicity of §4.4 is (iv).

**S18** (parity pairs `{D, D+1}`). They appear as the (even, odd) 3-step pairs inside every fan.

**S19/S20** (debt walk; two readers).
* The relation level between clusters is a debt walk. But S19's L (variance 4 per odd step, `L ≈ 2k`) is normalized differently from the chain level k (variance about 1/2 per Terras step).
* S20's long-time covariance limit concerns one coupled pair, not clusters.

**codex-grounding / codex-complement.**
* Their 0.512 paid density and anchored phase lemma describe 2-adic cells; our pointwise exits are single points of those cells.
* The fixed-head finiteness theorem (Beukers–Schlickewei) concerns first-hit ROOT words of `M_E` with a fixed middle head, and gives finitely many exponents per head. By codex-grounding's own scope statement it says nothing about ranked chains of child substitutions, which is what our certificates are.
* Proposition 4.1 covers that remaining class when the substitutions are deletions. A decreasing exponent rank does not give termination, because `s(K)` never decreases along the chain. The obstruction concerns time offsets, not head complexity.

**The two walls** (ANALOGY; a reading of `collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md` §0.1).
* codex-grounding's fixed-head finiteness sits at the INTEGRAL wall.
* The offset deficit of §4.1 sits at the DRIFT wall: the 7.64 steps per bit are the drift constant.
* We do not claim that variable heads escape the first wall.

**Prior art.**
* **Fast parity vectors.**
  * A.-S. Elsenhans, *Numerical verification of the Collatz conjecture for billion digit random numbers* ([arXiv:2502.16743](https://arxiv.org/abs/2502.16743), 2025; ANTS XVII). It uses the same recursive halving on the low bits, and its run time is empirically `O(N^(1+ε))`.
  * The affine form and the parity-vector bijection: R. Terras, Acta Arith. 30 (1976); C. Everett, Adv. Math. 25 (1977); J. Lagarias, Amer. Math. Monthly 92 (1985).
  * Divide and conquer on the low bits is also the device of D. Bernstein and B.-Y. Yang's constant-time gcd ([2019](https://eprint.iacr.org/2019/266)).
* **Merging trajectories.**
  * L. Garner, Discrete Math. 55 (1985) 57–64, on consecutive integers with equal heights.
  * M. Elia and A. Tucker, Integers 15 (2015) ([arXiv:1511.09141](https://arxiv.org/abs/1511.09141)), on counterexamples to Garner's criterion.
  * M. LaDue, *Clusters of integers with equal total stopping times in the 3x + 1 problem* ([arXiv:1709.02979](https://arxiv.org/abs/1709.02979), 2017): clusters of equal total stopping time, built from merging pairs of consecutive integers.
  * The classical pair 8n + 4, 8n + 5, which merges after three steps.
  * The repo's THM-4581, THM-4556 and THM-4593.
* **Stopping-time constants.** The total stopping time is about `6.95212 ln n` in Terras steps (Lagarias 1985; A. Kontorovich and J. Lagarias' stochastic models). Hence `σ_T(M_K) ≈ (1 + 6.95212 ln 3)K = 8.6377K`.
* **Trajectories of `2^n − 1`.**
  * OEIS [A075485](https://oeis.org/A075485) (T. D. Noe's table), A390816, and A193688 (A. Herrera, 2023: a(2n) = 1 + a(2n − 1)).
  * W. Ren et al. (2018) verified `2^100000 − 1` (863,323 halvings = Terras steps, 481,603 odd steps; recomputed in (M)).
* **Verification bound.** D. Barina, J. Supercomputing (2025) ([link](https://link.springer.com/article/10.1007/s11227-025-07337-0)): all n < `2^71`.
* **Coalescing random walks.** R. Arratia, thesis (1979) and [Ann. Probab. 9 (1981)](https://projecteuclid.org/euclid.aop/1176994264); M. Bramson and D. Griffeath (1980), where the one-dimensional density decays like `(πt)^(−1/2)`.

**The attached PDFs** (ANALOGY; their claims are recorded as claims).
* *Polynomial removal fails for ordered binary matrices* (OpenAI, Sep 25): a fixed 66×66 pattern whose `ε`-far matrices have only a `2^(−h)` fraction of copies, built in h levels.
  * Here the parameters left uncovered at depth W have measure about `W^(−α)`, with α ≈ 0.5–0.7 ((F)). Each certificate covers a cell of mass `2^(−depth)`.
  * So covering all but ε needs depth about `ε^(−1/α)`, hence about `2^(ε^(−1/α))` certificates: no polynomial removal for Mersenne certificate banks.
  * Their h-level construction mirrors the cluster levels.
* *Pointwise convergence of fourfold ergodic averages for mixing transformations* (OpenAI, Oct 4).
  * Mixing alone gives almost-everywhere convergence; here, THM-4581 gives almost-sure coalescence for Haar parameters.
  * The obligations concern specific integers, a null set. That is exactly the gap between an almost-everywhere theorem and the actual children.
* *Projective Hodge lines and ordinary Iitaka subadditivity* (OpenAI) proves the inequality `κ(X) ≥ κ(F) + κ(Z)` for fibrations. We found no inequality of that shape here, and record the paper as read, without a transfer.

## 6. Status and next steps

| item | status |
|---|---|
| certificate soundness; no-meeting lemma; fast parity vectors | PROVED |
| stopping-time invariant (Prop. 4.1); deletion routes (Cor. 4.2); offset accounting | PROVED |
| density-one exits for each D | PROVED (corollary of THM-4581 (3), §4.4) |
| reset rule; the chain; exponent 1/2 of the no-exit tail for finite lag sets | KNOWN (THM-4601 (i), THM-4581 (6d), THM-4593) |
| t = 23 deletions; fan closure; escape at 19,000,765 with 1916 deletions (deepest `M_99708991741`); 1949-exponent class; residue check; level 2 negative to `2^29 − 9` (D ≤ 4000) | FINITE-EXACT |
| progression: 1999/2000 first members certified; t = 3525143 negative to `2^26 − 8` | FINITE-EXACT |
| horizon clusters K ≤ 12800 (level sets) | FINITE-EXACT (from mac-mini's exact table) |
| exit-depth tails; coalescence profiles; N(T) ≈ (0.23–0.52)·√T | NUMERICAL |
| diffusive coalescence `N(T) ≍ √T` up to the horizon; O(√E) cluster span | HEURISTIC (HYP-9243) |
| PDFs; two-walls reading | ANALOGY |
| grounding of the 24 children; universal Collatz | OPEN |

**Next steps.**
1. **The one-orbit ground.** Compute `σ_T` for one member of the 1949-exponent class. By our estimate it needs more than 60 GB of memory and a few days on one core, and it grounds all 1949 members.
2. **Child types with large offsets.**
   * Candidates: non-deletion children whose time offsets are not paid by the rising phase (codex-complement's ternary entries and auxiliary bridges).
   * A child type is useful exactly when its offset `t_s − t_c` is large relative to its computation cost.
3. **Horizon clusters beyond K = 12800.** Exact orbits up to K ~ 10^5 would test HYP-9243's exponent.
