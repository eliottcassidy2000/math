# Lane `core`: pair chains, merge certificates and the −1 obstruction (golden session, 2026-10-07)

## 0. Bottom line

* **Collatz stays OPEN.** Nothing below proves it, and every claim carries a type.
* **A. Nonadjacent lags buy constants, not exponents, unless there are exponentially many of them.**
  * **Finite families.** For the partners y+1, …, y+R, q_R(T) ≥ c_R T^(−1/2).
    * PROVED for every R ≤ 32, and for S19's Mersenne lag sets D ≤ 7 and D ≤ 61, from explicit coalescence-first witness cylinders, each checked exactly with the pair chains.
    * For a general finite family it is a CONJECTURE: a witness has always been found when searched for.
    * With THM-4581(4) this gives α_R = 1/2, with the upper bound at sketch level.
    * So **"α_R > 1/2" is false.**
  * **The obstruction is partner coalescence.** With positive probability the partners merge with each other first. From then on they are one orbit, and the adjacent ±1 k-walk is back.
  * **Measured constants** (NUMERICAL, 10^5–2·10^5 paths to T = 10^6): q_R ≈ (c/R)·T^(−1/2), with c ≈ 11–14 for R ≤ 16 and 8–10 for R = 24–64.
    * R lags buy roughly a factor R (ANALOGY: CrocSwap's one factor of d), but no new exponent.
  * **S19's any-lag α ≈ 0.69 is reproduced**: 0.70 and 0.69 on [4·10², 4·10³] and [10³, 10⁴] for D ≤ 61.
    * It is a coalescence transient: the same family relaxes to 0.48 ± 0.03 on [10^5, 10^6].
  * **All translation lags.** With every translation lag allowed, q_∞(T) ≤ (T+1)·2^(−0.0346T) (PROVED: pigeonhole on (a, x_T mod 3^a)).
* **B. Merge certificates.**
  * Translation joins are classical (Roosendaal; Barina 2020). Depth-1 predecessor merging is in Angeltveit 2026.
  * **New: the maximal class-decided sieve**, which uses every backward-tree element below n.
    * Exact counts go to K = 34. It leaves 0.493% of residues mod 2^34, against 0.884% for descent.
    * At K = 30 it leaves 0.54× descent, 0.66× Roosendaal's sieve, and 0.82× the best classical combination.
  * It subsumes every join (FINITE-EXACT, K ≤ 30), because never-descending classes never merge with each other (FINITE-EXACT, s ≤ 30).
  * The rate C·K^(−3/2)·2^(−0.0500K) is unchanged; only C drops from about 5.6 to 3.1 (NUMERICAL).
* **C. Minimal counterexample.**
  * A minimal counterexample n₀ (> 2^71) lies in 6,915,181 of the 2^30 classes (0.644%). Every certificate threshold is ≤ 2^6.8, so this is PROVED given Barina's bound.
  * New explicit exclusions start at 2^10, e.g. n₀ ≢ 539 mod 1024 via m = (27n+7)/32.
  * The pair chains add nothing else.
* **D. Proof architecture.**
  * Bounded-memory Lyapunov weights and pathwise chain weights are impossible (PROVED).
  * The sieve leaves −1 uncertified at every depth (PROVED), and −5 and −17 to depth ≥ 88 (FINITE-EXACT, exhaustive search at each depth).
    * In both negative cycles these are the only uncertified points: the grand-orbit minima.
  * So a proof must use the sign, i.e. the 2-adic digit tail …000 versus …111.

## 1. What the sources say

* **THM-4581 and THM-4569** (worktree).
  * The chain (k, e) is driven by one fair bit per Terras step and is absorbed almost surely.
  * c T^(−1/2) ≤ P(no merge) ≤ C T^(−1/2)(log T)^2, with the upper bound at sketch level.
  * The integer forms are density-one statements.
* **Roosendaal**, ericr.nl/wondrous/techpage.html (read with curl).
  * The convergence sieve checks only numbers that "can not be shown to either converge or join the path of a lower number".
  * 2^16 leaves 1720 of 65536 residues. At 2^32, more than 99.2% are skipped. Residues 2, 4, 5, 8 mod 9 are skipped as well.
* **Barina**, J. Supercomput. 2020, §4: the same sieve, footnote 4 "Enhancement proposed by Eric Roosendaal". Barina 2025 verifies to 2^71.
* **Angeltveit**, arXiv:2602.10466 (2026).
  * Uses descent, mod-9 preimages, and a depth-1 "path-merging" sieve: if m = T^k(n) ≡ 2 mod 3 and (2m−1)/3 < n, skip n.
  * The deeper version "is less effective and we did not implement it".
  * Mod 27 adds no cases; mod 81 adds one.
* **CrocSwap**: nonadjacent swaps cost O(d) against O(d²) for adjacent moves.

## 2. Derivations and computations

### A. Nonadjacent lags

* **A1 (PROVED). One coin per step.**
  * Every chain is driven by the same bit β; at each flip, k ↦ k + 1 − 2β.
  * Once two partners merge, their chains are identical forever.
* **A2 (PROVED given a witness).**
  * Suppose that on a cylinder E = {y ≡ ρ mod 2^t₀} all partners have merged with each other by time t₀ and y with none.
  * Then for T ≥ t₀, q(T) ≥ 2^(−t₀)·P_{s₀}(the single chain y vs. cluster is not absorbed by T) ≥ c T^(−1/2). This uses the SRW-skeleton lower bound of THM-4581(4).
  * Witness cylinders exist for every R ≤ 32 (`witnesses_full.txt`; all 31 re-verified exactly with the pair chains) and for the Mersenne lag sets D ≤ 7 and D ≤ 61 (`A/mers_D*_witness.txt`).
    * D ≤ 7: 13% of paths are witnesses by T = 4000.
    * D ≤ 61: 0.26% of paths; the first witness is at time 2376.
  * Explicit case R = 2, y ≡ 21 mod 32:
    * at time 5, y+1 and y+2 both equal 27q+20, while y is at 3q+2;
    * the state is (k, e) = (2, 2), so q₂(T) ≥ 2^(−5)·2√(2/(πT)) ≈ 0.05 T^(−1/2).
* **A3 (PROVED). All lags.**
  * Let w be y's word (length T, weight a) and c(w) = 2^T T^T(y) − 3^a y.
  * The word map is a bijection of Z₂. So y merges with some y+r (r ≥ 1) by time T iff some w′ of weight a has c(w′) ≡ c(w) mod 3^a and c(w′) < c(w). The partner is r = (c(w) − c(w′))/3^a < 2^T. Coincidences with unequal odd counts are Haar-null.
  * Hence q_∞(T) = N_T/2^T, where N_T = #{(a, c mod 3^a)} ≤ Σ_a min(C(T,a), 3^a) ≤ (T+1)·2^(H(x*)T).
  * Here x* = 0.60909 solves H(x) = x·log₂3, so the rate is 0.03462 bits per step.
  * Exact values: q_∞(20) = 0.235 and q_∞(24) = 0.192, against q₁(20) = 0.598 (`allags.c`).
* **A4 (PROVED). The −1 shadow.**
  * During a run of L odd steps of y, no chain is absorbed and every flip lowers k.
  * Afterwards k^(r) = O_L(r−1) − L exactly. Partners r ≥ 2 sit in the band {−L/2, −L/2 − 1}, already coalesced in k; partner 1 sits at −L.
  * No positive-measure set of permanent joint failure exists, since each chain is absorbed almost surely.
  * Near −1, the relations fail together on measure 2^(−L), which only adds to the constant. Coalescence (A2) is what pins the exponent at 1/2.
* **A5 (NUMERICAL).** Translations, with √T·q averaged over T ∈ [10^5, 10^6]:

| R | 1 | 2 | 3 | 4 | 6 | 8 | 12 | 16 | 24 | 32 | 48 | 64 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| c_R | 11.20 | 6.64 | 4.60 | 3.31 | 2.31 | 1.47 | 1.04 | 0.82 | 0.43 | 0.28 | 0.22 | 0.13 |
| R·c_R | 11.2 | 13.3 | 13.8 | 13.3 | 13.9 | 11.7 | 12.5 | 13.1 | 10.2 | 9.1 | 10.5 | 8.1 |
| α[10²,10³] | 0.24 | 0.30 | 0.35 | 0.38 | 0.44 | 0.49 | 0.55 | 0.59 | 0.67 | 0.72 | 0.78 | 0.82 |
| α[10⁵,10⁶] | .49±.01 | .50±.01 | .52±.02 | .51±.02 | .48±.02 | .49±.03 | .44±.05 | .43±.05 | .44±.07 | .51±.09 | .53±.11 | .48±.14 |

* **Mersenne (S19) families** (Terras time).
  * Lag 1 alone: √T·q₁ → 16.7. This validates against S19 and audit A.
  * D ≤ 61: local slopes 0.70, 0.69, 0.64, 0.59, 0.51, 0.48 ± 0.03 across [4·10², 10⁶]. The plateau is √T·q ≈ 2.6.
  * D ≤ 241: slopes 0.71 down to 0.59 ± 0.06 on [10^5, 10^6], with √T·q still falling (1.8 at 10^6).
  * D ≤ 961 matches D ≤ 241 up to T = 10^5.
* **Heuristic picture.**
  * Translation partners start "stacked" at k = 0 and gain a factor of about R. This fits Arratia-type coalescence: y must stay apart for the partners' coalescence time, of order R².
  * Mersenne lags start spread out at k = −D, so the nearest lag shields the others.
  * Pre-asymptotic slopes above 1/2 last until the partners have coalesced.

### B. Exact sieves

* **Setting and certificates.**
  * A class is a parity word of length K, optionally with n mod 3^J. On the class, x_s = (3^(a_s) n + c_s)/2^s.
  * Descent: 3^(a_s) < 2^s.
  * Join: the chain from (0, −d) is absorbed. This covers all d, since d < 2^(s−a_s).
  * Branch: m = (2^i x_s − c(u))/3^b for a backward word u with b ≤ a_s + J odd steps and ratio 3^(a_s−b)·2^(i−s) < 1.
* **Maximality (PROVED).** Every certificate given by one affine m(n) on the class is of this form. Ratio exactly 1 means a translation, i.e. a join.
* **Validation.**
  * Descent reproduces A076227 to K = 36.
  * Joins give 1720 at 2^16 and 99.211% skipped at 2^32, matching Roosendaal.
  * Branches equal an independent big-integer brute force for every K ≤ 11.

| K | descent | +joins (Roosendaal) | +depth-1 (Angeltveit-type) | depth-1+joins | **maximal** |
|---|---|---|---|---|---|
| 16 | 2114 | 1720 | 1856 | 1562 | **1363** |
| 24 | 286581 | 234156 | 244392 | 206402 | **172868** |
| 30 | 12771274 | 10446423 | 9910223 | 8385079 | **6915181** |
| 32 | 41347483 | 33880411 | 34512882 | — | **23797887** |
| 34 | 151917636 | — | — | — | **84720794** |

* **B2 (FINITE-EXACT). Joins are subsumed.**
  * There are 0 join-only classes for every K ≤ 30.
  * Mechanism: if n merges with n−d at time s and n−d descends at j < s, then T^j(n−d) is a branch certificate (ratio 3^(a_j)/2^j < 1, b = a_s − a_j).
  * So subsumption reduces to the **survivor lemma**: no two never-descending words of equal length and weight collide mod 3^a.
    * Checked: 0 collisions for every s ≤ 30 (12.8M words at s = 30).
    * Beyond s = 30 it is a CONJECTURE.
* **B3 (NUMERICAL). Rates.**
  * With ρ = 2^(H(log₃2)−1) = 0.965907, C(K) = f·K^(3/2)/ρ^K stays flat:
    * descent 5.4–5.8; joins 4.3–4.6; depth-1 plus joins 3.5–3.6; maximal 2.9–3.2.
  * The ratio maximal/descent lies in 0.54–0.60 for K = 24–34.
  * HEURISTIC reason: the expected number of certifying backward words at walk depth h is ≤ C_λ e^(−λh) for 1 < λ < 2, so certificates fire only near the ends of the survivor excursion.
* **B4 (FINITE-EXACT). 3-adic refinement.**
  * J = 1 and J = 2 give exactly the classical factors 2/3 and 5/9. The mixed extras number 0, 2 and 139 at K = 20, 24, 28.
  * J = 3 adds nothing. J = 4 removes 2.3% more than the classical "10 mod 81" case.
  * At K = 28 with n mod 9: 0.409% remain, against 0.597% (Roosendaal ×5/9), 0.579% (depth-1 ×5/9) and 0.490% (both).

### C. Minimal counterexample

* **C1 (PROVED, given Barina 2025, KNOWN).**
  * Every certificate at depth ≤ 30 holds for n > 2^6.76 (`sieve_bounds.c`), and m ≥ 1 once n > 2^30.
  * Hence a least non-convergent n₀ has n₀ mod 2^30 in a listed set of 6,915,181 residues, and n₀ mod 9 ∈ {0, 1, 3, 6, 7}.
* **C2. Exclusions beyond the sources checked.**
  * The first ones are 539 and 615 mod 2^10, via depth-5 words: m = (27n+7)/32 and (27n+3)/32.
  * Example: 455 → … → 1154 = T^10(539). There are 17 such classes at 2^12.
* **C3 (PROVED, nothing new).**
  * For a counterexample, every chain (n₀, n₀−d) with d < n₀ is never absorbed.
  * Its k drifts to −∞ with slope ≤ 1/2 − log₃2 = −0.131, and e_t tends to the values of the 1–2 cycle.
  * This restates the classical odd-density bound.

### D. Proof architecture

* **D1. The conjecture that remains.**
  * Every positive integer n > 2^71 has, at some depth K, a class-decided certificate with threshold < n.
  * This is equivalent to Collatz. Equivalently U_∞ ∩ Z_{>0} = ∅, where U_∞ is the closed, Haar-null set of 2-adic integers whose class is never certified.
* **D2. Why measure methods stop.**
  * Haar measure is blind to the sign: both Z_{>0} and Z_{<0} are dense and null.
  * U_∞ ∋ −1 (PROVED: on the class −1 mod 2^K every backward word with b ≤ a has ratio 3^(j−b)·2^(i−j) ≥ 1).
  * The classes of −5 and −17 stay uncertified to depths 93 and 88; the search is exhaustive at each of those depths.
    * Every other point of their cycles is certified by depth 9 (`path_cert.py`). In both cycles the sieve isolates exactly the grand-orbit minimum.
* **D3. Bridges ruled out (PROVED).**
  * (i) Weights log n + g(n mod M), for any M and g. Take n ≡ −1 mod 2^(L+P) on a period-P orbit of r ↦ (3r+1)/2 mod M_odd: the g-terms telescope while the log terms add P·log(3/2).
  * (ii) Pathwise weights V(k, e) that decrease on every transition.
    * The positive pair (2, 1) gives the cycle (0, 1) → (−1, 1/3) → (0, 1) under bits (10)^∞.
    * The states (h, 3^h − 1), the Mersenne relations at −1, are fixed under all-odd bits.
  * (iii) Any criterion stated only for rational or eventually periodic 2-adic points; −1, −5 and −17 satisfy it too.
* **D4. The bridge must use information beyond the digit horizon.**
  * Test (FINITE-EXACT up to the search budget; 37 stopping-time records, 5 of them rechecked unchanged with a 50× budget): the certificate depth K(n) is 0.77–1.00 × σ(n), and up to 12.4 × log₂ n (n = 27: K = σ = 59).
  * So certificates are not confined to the log₂ n free bits, which is exactly where density methods stop (cf. Tao 2019).
  * Most promising direction: bound the certificate depth of a positive integer by archimedean data on its 0-tail. No such bound is known, and any bound must fail for −5 and −17.

## 3. Suggested canon filings

1. **THM. Partner coalescence pins the diffusive exponent.**
   * A1–A4: PROVED.
   * α_R = 1/2: sketch level, given the THM-4581(4) upper bound.
   * A5: NUMERICAL.
   * Evidence: `witnesses_full.txt`, `A/`, `allags.c`, `qsurv.c`.
2. **THM. The maximal class-decided Collatz sieve.**
   * Definition and maximality: PROVED. Counts to K = 34 and thresholds to K = 30: FINITE-EXACT.
   * Constraints on n₀: PROVED given Barina 2025. Join subsumption: FINITE-EXACT. Rate: NUMERICAL.
   * Prior art: Roosendaal; Barina 2020; Angeltveit 2026 (depth 1 only).
3. **HYP. Survivor lemma.** FINITE-EXACT for s ≤ 30; it implies joins ⊂ branches.
4. **PROP. The −1 and sign obstructions D2–D3.** PROVED, with FINITE-EXACT depths for −5 and −17.
5. **Update to HYP-9217 (2).**
   * Every finite lag set has α = 1/2: lower bound PROVED via witnesses (D ≤ 7, D ≤ 61); upper bound at sketch level.
   * The measured 0.69 is the coalescence transient (reproduced).
   * Unbounded D: OPEN. The local slope is 0.59 ± 0.06 on [10^5, 10^6] and still falling. A heuristic coalescing-front argument predicts 1/2.

## 4. Files (`scratchpad/golden/core/`)

* `sieve.c` (outputs `out_*.txt`), `sieve_bounds.c` (`out_bounds_K30.txt`).
* `verify_branch.py`, `classify_small.py`, `survcollide.c`.
* `allags.c`, `witness_search{,2}.py` → `witnesses_full.txt`.
* `qsurv.c`, `qwit.c` → `A/`. Checks: `trace_harness.c`, `check_qsurv_chain.py`. Analysis: `slopes.py`, `fit_rates.py`.
* `path_cert.py`, `cert_depth.py` (`.out`, `_bigbudget.out`).
* Sources: `barina2020.txt`, `roosendaal_techpage.txt`.
