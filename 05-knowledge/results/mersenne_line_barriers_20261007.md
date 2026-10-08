# Mersenne-line barriers: the t = 23 fan escapes, the unpaid progression is reached, and the line coalesces diffusively

Session opus-2026-10-07-S21.
* **Owner directive:** take codex-grounding's status (t = 23 reached; giant children ungrounded; "the next decisive target is a terminating, source-preserving rule for one of those actual children") as the brief. Search past repo work for connections while testing hypotheses. Three PDFs were attached.
* **Scripts:**
  * `04-computation/experiments/mersenne_line_barriers_20261007.py` (+ `.out`; about 3 min);
  * `mersenne_fan_escape_20261007.py` (+ `.out`; about 6 min, 2 GB).
  * Both print ALL CHECKS PASSED.
* **Canon:** THM-4605 and HYP-9243.

**Status.**
* PROVED: certificate soundness; correctness of the fast parity vectors. The reset rule on the Mersenne line is KNOWN (THM-4601 (i) / THM-4556).
* FINITE-EXACT: every certificate and negative bound below.
* NUMERICAL: tails and coalescence profiles.
* HEURISTIC: the diffusive-coalescence explanation of the orphan law and its consequence for grounding.
* ANALOGY: the PDFs.
* Collatz OPEN.

## 0. Results

1. **An independent certificate engine, orders of magnitude cheaper.**
   * A Mersenne exit `M_E ⇝ M_(E−D)` is the absorption of one pair chain (THM-4581's table, with the source's own orbit as the reference). Its coins are the Terras parities of `x = 2·3^(E−1) − 1`, which need only `3^(E−1) mod 2^W`.
   * The engine reproduces codex-grounding's t = 23 result exactly: deletions `{3, 4, 7, …, 28}`, all absorbed at Terras time 6553.
   * Divide-and-conquer parity vectors (python-flint) compute the `2^25` coins of the fan bottom in about 12 s.
2. **The unpaid progression `t ≡ 13847 mod 27648` is reached member by member.**
   * 1999 of its first 2000 members have certified exits with fixed source.
     * 1970 have exits within 16384 steps; median cheapest merge time 157.
     * 29 needed depth up to 2.5·10⁶.
   * One member, `t = 3525143`, has no exit with `D ≤ 256` within `2^26`.
3. **The supplied giant children were trapped in a closed fan, and the fan escapes at Terras time 19,000,765.**
   * The source and its 24 children form a merge clique: every exit of every child lands on another fan member.
   * The fan bottom `M_99708993677` has no exit with D ≤ 256 before Terras time 19,000,765. Below it, all deletion chains collapse into one chain at step 2,741,933.
   * That chain's relation level makes a long negative excursion (minimum level −3867, after the collapse at step 2,741,933) and returns: it is absorbed at Terras time 19,000,765.
   * Certified: `M_99708993677 ⇝ M_(99708993677 − D)` for all D ≤ 256, and for D = 1910 (the block edge; D = 1911 is outside this event).
   * Hence the supplied source `M_99708993705`, all 24 children and the original parent `M_99708993713` have a certified common future with `M_99708991767`.
   * Independent check: the two orbits agree mod `2^(2^20)` at that time, from affine maps rather than the chain.
4. **The next barrier is deeper than 1.3·10⁸.** The cluster below `D = 1911` (tested `D = 1911, 1912, 2500, 4000`) collapses into one chain whose level reaches −12297, and it is not absorbed within `2^27` steps.
5. **Why the giant children stay ungrounded (HYP-9243).**
   * The Mersenne line coalesces **diffusively**. On six fresh odd exponents near 10^11, the contiguous block of smaller exponents absorbed by time T has median `S(T) ≈ (0.3–0.4)·√T` at `T = 2^17–2^20`, growing in heavy-tailed jumps. By `T = 2^20` only 1–4 clusters cover 2048 consecutive deletions.
   * Truncate at the orbit horizon `σ_T(M_K) ≈ 8.64 K` (mac-mini's exact data). The line is then left in about `√K` clusters with no deletion route between them.
   * That predicts an orphan fraction of `(0.76–1.1)/√K`. mac-mini's HYP-9242, measured exactly to K = 12800, gives fraction·√K ≈ 0.9–1.1. Two independent computations at scales 10^4 and 10^11 agree, within the poorly determined constant.
   * **Consequence (HEURISTIC).** Deletion certificates alone cannot ground a generic giant Mersenne number. Its horizon cluster has size about √K and almost surely contains no small exponent. Barrier depths along any descent reach the orbit horizon. Pointwise grounding of the supplied children by this route would cost about as much as their orbits.
6. **The reset rule (KNOWN).** For even E, `M_E ⇝ M_(E−1)` at Terras time 3: with `w = 3^(E−2) ≡ 1 mod 8`, both points reach `(9w − 1)/4`. This is THM-4601 (i), the lag-one form. Every hard exponent observed is odd.

## 1. The certificate

* `T` is the Terras map, `M_E = 2^E − 1`, and `T^j(M_E) = 3^j 2^(E−j) − 1` for `j ≤ E`.
* So the source sits at `x = 2·3^(E−1) − 1` after `E − 1` steps, and the child `M_(E−D)` at `y = 2·3^(E−D−1) − 1`, with `y + 1 = 3^(−D)(x + 1)`.
* Write `y = 3^k x + e` and run THM-4581's pair chain with the source as reference:
  * coin = parity of `x_n`;
  * integer state `(k, a)` with `a = 3^max(0,−k) e`, starting at `(−D, 1 − 3^D)`;
  * the four integer updates are listed in the script header.

**Theorem (soundness; PROVED).** If the chain reaches `(0, 0)` at time n, then `T^n(x) = T^n(y)`, so `M_E` and `M_(E−D)` have a common future: `M_E` reaches 1 iff `M_(E−D)` does.

*Proof.*
* The relation `T^i(y) = 3^(k_i) T^i(x) + e_i` holds exactly at every step: the four cases of the table, as in THM-4581.
* `(k_n, e_n) = (0, 0)` is equality.
* The coins for `i < n` depend only on `x mod 2^(i+1)`.
* Relation absorption excludes trivial-cycle value coincidences, which it does not claim (S19's lesson, MISTAKE-580). ∎

**Validation.**
* 2400 large-E cases agree with direct orbits (202 merges).
* All 306 small-E disagreements are cycle coincidences at value ≤ 2.
* A numba prototype was explored for speed only. Its int64 overflow guard abandons chains rather than claiming them, and it is not used for any certificate.

**Fast parity vectors.**
* `par(x, n)` returns the parity vector V, its weight a, and c with `T^n(x') = (3^a x' + c)/2^n` for all `x' ≡ x mod 2^n`.
* The two halves combine by `x_mid = (3^(a1) x + c1) >> n1`, `a = a1 + a2`, `c = 3^(a2) c1 + 2^(n1) c2`.
* Cost is `O(M(n) log n)`. It agrees with naive iteration, including on the fan bottom (script (P)).

## 2. The supplied residual and its fan (FINITE-EXACT)

* **The source.** `E = 99708993705` (t = 23 in `E(t) = 924745897 + 2^32 t`) has exits exactly `D ∈ {3, 4, 7, …, 28}`, all at Terras time 6553 (W = 7200, D ≤ 40). These are codex-grounding's 24 children.
* **The original parent** `M_99708993713` exits to the fan bottom with D = 36 at 6553, matching codex-grounding's composed 36-bit deletion. It also exits to D = 1–4 at time 19.
* **The fan is closed.**
  * With W = 16384 and D ≤ 128, every exit of every child lands on another fan member. Child `E − D` has exactly `25 − (index)` exits, all inside the fan.
  * The bottom `E − 28` has none.
  * The internal cascade:
    * each (even, odd) pair merges at time 3 (the reset rule);
    * pairs merge with neighbouring pairs later;
    * everything merges with the source at 6553.
* **The bottom is deep.** No exit with `D ≤ 256` occurs before Terras time 19,000,765: (G) returns the first absorption among all 256 chains.
* **Below it, the chains collapse.**
  * All 256 chains are one from step 2,741,933 on, the step at which the whole lower block coalesced.
  * The orbit itself is ordinary (odd density 0.50001 over `2^25` steps).
  * The barrier is the relation level, not a rising phase.
* **The escape** (`mersenne_fan_escape_20261007.py` (G)). The collapsed chain's level makes a long negative excursion (minimum level −3867, after the collapse at step 2,741,933) and returns. It is absorbed at **Terras time 19,000,765**. The absorbed group is all of `D = 1, …, 256`.
* **Independent check (H).**
  * `T^n(x) ≡ T^n(y_D) mod 2^(2^20)` at n = 19,000,765 for D = 1 and D = 256, computed from the affine maps of the fast parity vectors.
  * The odd-step counts differ by exactly D (9500296 − 9500295 and 9500551 − 9500295).
  * The residues differ at n − 1.
* **Block edge (I).** The chains for D = 1500 and 1910 collapse into the absorbed chain before 19,000,765, at step 2,741,933, where the whole lower block coalesced. D = 1911 stays separate. The deepest certified jump is therefore `M_99708993677 ⇝ M_99708991767`.
* **Level 2.** The next cluster (D = 1911, 1912, 2500, 4000 collapse together) is not absorbed by `2^27 = 134,217,728` steps; its level reaches −12297.

**What this does for the owner's obligation.**
* The 24 child ports are no longer trapped. They have certified strictly smaller descendants down to `M_99708991767`: a further source-preserving child rule, with rank the exponent.
* They are still not grounded. Grounding needs every lower barrier crossed, and §4 explains why that route cannot finish for a generic giant exponent.

## 3. The unpaid progression `t ≡ 13847 mod 27648` (FINITE-EXACT)

* Source `M_E(t)` with `E(t) = 924745897 + 2^32 t` stays fixed. Each certificate deletes D bits to a smaller Mersenne.
* Examples (largest D among the cheapest exits): t = 13847 (D = 6 at Terras time 273), 41495 (D = 36 at 957), 69143 (D = 4 at 73), 96791 (D = 12 at 129).
* **First 2000 members.**
  * 1970 have certificates within W ≤ 16384 (D ≤ 64), by first successful W: 256 (1258), 1024 (489), 4096 (169), 16384 (54). Median cheapest merge time is 157.
  * 29 more are certified by the deep search (D ≤ 256, W ≤ 2^22). Their merge times run from 16615 to 2,520,494, and absorbed groups contain up to 256 consecutive deletions.
  * The remaining member `t = 3525143` (`E = 15140374823469225`) has no exit with D ≤ 256 within `2^26`. It is a deep barrier like the fan bottom.
* Each certificate covers only a 2-adic cell of mass about `2^(−(cost−34))` (codex-grounding's anchored phase lemma). So these are pointwise exits: the requested "explicit member", not new density.

## 4. Diffusive coalescence and the orphan law (HYP-9243)

**Random giant exponents** (NUMERICAL; 600 odd E in [10^10, 10^11]; even E always exit at time 3):
* `P(no exit with D ≤ 64 by W)` = 0.540, 0.383, 0.258, 0.160, 0.092, 0.047 at W = 128 … 4096. That is about `W^(−0.9)`, steeper than one chain's `W^(−1/2)` because several deletions compete before they collapse.
* Locally below the t = 23 source, 2.9% of exponents are hard at 4096, all odd (reset rule).

**Coalescence profile** (`mersenne_line_coalescence_20261007.py` (S); NUMERICAL). On six fresh odd E near 10^11, chains for D = 1…2048 collapse and are absorbed. The contiguous absorbed extent below the source:

| T | S(T) per exponent | upper median of six |
|---|---|---|
| 4096 | 68, 10, 8, 2, 0, 36 | 10 |
| 2^17 | 136, 212, 50, 312, 0, 78 | 136 |
| 2^20 | 1796, 212, 304, 312, 0, 78 | 304 |

* The local slope of the median is 0.4–0.6, and the median `S(T)/√T` is 0.38 at `2^17` and 0.30 at `2^20` (individual values 0–1.8).
* By `2^20`, 1–4 live clusters remain among 2048 deletions.
* The t = 23 fan fits the same picture: block extent 1910 at T = 2.74·10⁶, against √T ≈ 1655.

**Mechanism.**
* Between two clusters the relation level k is a recurrent walk, the debt walk of S19/S20 (THM-4581), with no drift.
* Exponents D apart start at level −D, so they coalesce on time scale about D², heavy-tailed. Clusters at time T then have size about √T.

**The orphan law follows (HEURISTIC).**
* Orbits of `M_K` end after `σ_T ≈ 8.64 K` Terras steps (mean over 1000 ≤ K ≤ 12800, from mac-mini's exact `mersenne_sigma_12800.txt`).
* At that horizon the line is cut into clusters of size about `c√(8.64 K)`. Each cluster bottom is an orphan: no deletion partner merges before 1.
* So the orphan fraction is about `1/(c√(8.64K))`, i.e. `(0.76–1.1)/√K` for c between 0.30 and 0.45.
* HYP-9242 measures 0.89–1.14/√K exactly for K ≤ 12800, inside that band. The constant is not fitted, and c is poorly determined by six profiles.

**Consequence (HEURISTIC).**
* A deletion descent from a giant M_E must cross every barrier below it. Barrier depths grow with the cluster scale up to the horizon.
* The final cluster of a generic giant exponent has size about √E and contains no verified exponent.
* So deletion certificates alone cannot ground generic giant Mersennes. The t = 23 children face a level-2 barrier beyond 1.3·10⁸, and further levels beyond that.
* This is a second stopping reason, next to codex-grounding's fixed-head finiteness:
  * fixed heads give finitely many exponents (Beukers–Schlickewei);
  * variable heads exist pointwise, but their depth is unbounded and heavy-tailed, reaching the orbit horizon.
* Universality needs either non-deletion children whose structure survives recursion, or a uniform statement about the 2-adic digits of `3^E` controlling every barrier. The second is a Mahler-3/2-type wall (`00-navigation/MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md`).

## 5. Connections

* **THM-4581** (mac-mini, Haar coalescence). Our certificates are its absorption events, on actual integers. The almost-sure statement is the measure-one shadow of the barrier hierarchy, and its `T^(−1/2)` rate is the single-chain tail behind HYP-9243.
* **HYP-9242** (orphans). Its law is the horizon limit of diffusive coalescence (§4). Its heuristic (debt martingale; time ∝ log n) is the same mechanism, with the coalescence profile measured directly here.
* **THM-4601** (two anchors).
  * (i) is the reset rule (§0.6).
  * The +1 barrier delays odd Mersenne exponents by only `j = 3 + v_2(K−1)` steps (5 for the fan bottom). Deep barriers are statistical, not that one.
* **S18** (parity pairs `{D, D+1}`) appears as the (even, odd) 3-step pairs inside every fan.
* **S19/S20** (debt walk; two readers).
  * The relation level between clusters is the debt walk L. Its diffusive spread is the "variance 4 per step" of S19, and THM-4565's sibling merges are the fan's internal cascade.
  * S20's long-time limit `Cov → 2·1{d = −1}` (given THM-4581) is the statement that clusters eventually merge.
* **codex-grounding / codex-complement.**
  * Their 0.512 paid density and anchored phase lemma describe 2-adic cells. Our pointwise exits are single points of those cells.
  * The fixed-head finiteness theorem is consistent with the observed heads. Every certificate is a different head, of length up to 1.9·10⁷.
* **Mahler 3/2.** A uniform bound on barrier depths would be a uniform regularity statement for the low binary digits of `3^E` along all integers E. This is the same kind of statement as the Z-number and Mahler problems.
* **The two walls** (`collatz_coalescence_20261004_proof_strategies_and_bold_predictions.md`, §0.1). The two stopping reasons for the Mersenne-line program are the two walls of that atlas.
  * codex-grounding's fixed-head finiteness is the INTEGRAL wall: Beukers–Schlickewei S-units allow finitely many exponents per fixed head, i.e. bounded complexity.
  * HYP-9243 is the DRIFT wall: certificates exist pointwise almost surely (THM-4581), but barrier depths reach the orbit horizon, and the about √K orphan clusters left at the horizon are exactly the pointwise residue.
  * As the atlas predicts, the program sees each wall separately: variable heads escape the first, and the second is not escaped.

**The attached PDFs** (ANALOGY; claims as claims):
* *Polynomial removal fails for ordered binary matrices* (OpenAI, Sep 25): a fixed 66×66 pattern with `ε`-far matrices having only `2^(−h)`-fraction copies, built in h levels.
  * Here, parameters left uncovered at depth W have measure about C/W, polynomial in W. Each certificate covers a cell of measure `2^(−W)`.
  * So covering all but ε of the parameter space needs about `2^(C/ε)` certificates: no polynomial removal for Mersenne certificate banks. Their h-level construction mirrors our coalescence levels.
* *Pointwise convergence of fourfold ergodic averages for mixing transformations* (OpenAI, Oct 4).
  * Mixing alone gives almost-everywhere convergence. Here, THM-4581 gives almost-sure coalescence for Haar parameters.
  * The obligations concern specific integers, a null set. That is exactly the gap between an almost-everywhere theorem and the actual children.
* *Projective Hodge lines and ordinary Iitaka subadditivity* (OpenAI): `κ(X) ≥ κ(F) + κ(Z)` for fibrations.
  * The grounding cost of a source splits over the cluster hierarchy: within-cluster cost (fibre) plus between-cluster barriers (base).
  * The "ordinary" (non-log) setting corresponds to deletion children only.

## 6. Status and next steps

| item | status |
|---|---|
| certificate soundness; fast parity vectors | PROVED |
| reset rule for even E | KNOWN (THM-4601 (i)) |
| t = 23 deletions; fan closure; fan-bottom negatives to 2^20; escape at 19,000,765 (all D ≤ 256, and D = 1910); residue check; level-2 negative to 2^27 | FINITE-EXACT |
| progression: 1999/2000 first members certified; t = 3525143 negative to 2^26 | FINITE-EXACT |
| exit-depth tails; coalescence profiles; S(T) ≈ (0.2–0.5)√T | NUMERICAL |
| diffusive coalescence ⇒ orphan law ⇒ deletion descent cannot ground generic giant Mersennes | HEURISTIC (HYP-9243) |
| PDFs | ANALOGY |
| grounding of the 24 children; universal Collatz | OPEN |

**Next steps.**
1. Level 2 for the fan: depth `2^29`–`2^31` needs a streaming parity-vector computation (memory).
2. A child type that survives recursion and can enter orphan clusters (codex-complement's ternary entries / auxiliary bridges), tested against the cluster structure.
3. Measure `S(T)` to `2^25` on more exponents, and fit the constant against HYP-9242's.
