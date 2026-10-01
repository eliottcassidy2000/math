# LRC(14), Collatz and the primes: "local versus global" made precise

**Status.** SYNTHESIS with PROVED (elementary) lemmas, FINITE-EXACT censuses, CITED literature (arXiv, 2017–2026) and EMPIRICAL statistics. No HYP or THM file was created and no canonical status is changed. **Collatz is OPEN. LRC(14) is OPEN in the repo**, but see the flag in §6.1: an arXiv preprint of 2026-09-02 claims a computer-assisted proof of exactly the repo's LRC(14) (and LRC(15)). That claim is CITED here and was not verified.

* **Labels.** Every claim carries one of PROVED, CITED, FINITE-EXACT, VERIFIED, EMPIRICAL, CONDITIONAL, OPEN, REFUTED. Every bridge is typed REAL, ANALOGY or NUMEROLOGY, and the test that decided the type is named.
* **Session and lane.** Session `collatz-procgen-20260922`, lane "localglobal", 2026-09-30.
* **Brief.** This lane answers the owner's thesis of 2026-09-30:
  * LRC and Collatz are "finely tuned to be aligned with the major subtle patterns in the primes";
  * LRC is *local* and Collatz is *global*, so "each lacks information about the other's type".
* **Scripts** (`04-computation/experiments/`):
  * `procgen_localglobal_20260930_lrc.py`: exact `M(V)`, lonely counts at denominator `q`, coverage;
  * `procgen_localglobal_20260930_gates.py`: exact gate counts by meet-in-the-middle and residue DP, primitive-word counts, the census, a float DP for large `p`;
  * `procgen_localglobal_20260930_duality.py`: the Poisson identity, the single-relation lemma, local badness modulo `ell`, the tight lines;
  * `procgen_localglobal_20260930_run.py`: the runner. It re-checks every finite claim and ends `ALL CHECKS PASSED`.
* **Output:** [`procgen_localglobal_20260930.out`](procgen_localglobal_20260930.out). Section tags such as A1.5 or B5 below point into it.
* **Conventions.**
  * LRC: a speed set `V` has `k` distinct positive integer speeds, and `delta = 1/(k+1)`; the repo's LRC(14) is `k = 13`.
  * `M(V) = max_t min_v ||t v||`. The slack is `sigma = (k+1)M - 1`, and `sigma >= 0` is LRC for `V`.
  * `V` *covers* `q` when `q` divides some speed. `V` is *divisor-complete* (DC, the repo's "covering" class) when it covers every `q` in `{2..k+1}`.
  * `lld(V)` is the least denominator of a lonely time. `u(V)` is the least integer dividing no speed.
  * Collatz: a word `w` of length `p` has `a` ones at `s_0 < ... < s_{a-1}`. Its carry is `c_w = sum_i 3^{a-1-i} 2^{s_i}`. The gate is `D = 2^p - 3^a`.
  * The rational periodic point of `x/2, (3x+1)/2` with word `w` is `x_w = c_w/D`, so a cycle of shape `(p,a)` exists iff `|D|` divides some `c_w` (THM-4471 §4, gates lane).

---

## 0. Answers in brief

### 0.1 Is LRC aligned with the primes?

**Yes in its first layer, which is the trivial one. No in its hard layer, which is additive.**

* **(PROVED, classical.)** If `V` misses some `q <= k+1`, then `t = 1/q` is a lonely time. So the class of sets that are hard to *prove* is exactly the divisor-complete class. In that weak sense "hard = covering the small primes" is true, and trivial.
* **(REFUTED as an extremal statement; FINITE-EXACT.)** The hardest instances do not cover:
  * The tight sets (`M = 1/(k+1)`) for `k <= 7` with speeds `<= 21` are exactly `{1..k}`, `{1,3,4,7}`, `{1,3,4,5,9}`, `{1,2,3,4,5,7,12}` and `{1,4,5,6,7,11,13}` (the sets of Wills and of Goddyn–Wong). All of them **fail** to cover `k+1`, and all are lonely at the first uncovered denominator `k+1`.
  * In the same exhaustive range, divisor-complete primitive sets have `M >= 2/(2k+1) > 1/(k+1)`.
  * Over 21 random cells (A1.5), the depth of prime coverage is only weakly correlated with the slack: Spearman lies between −0.26 and +0.05.
* **(EMPIRICAL / FINITE-EXACT.) What the primes do decide is where the certificate lives.**
  * `lld(V) = u(V)`, the first integer dividing no speed, for 87–98.5% of random rows (A1.5) and for every tight row.
  * Denominators that fail without being covered average 0.015–0.25 per random row.
  * The genuine LRC difficulty starts where divisibility stops. Inside the covering class the extra failing denominators in `(k+1, 2k+1)` are blocked *additively*, by an AP block of the speeds acting through Dirichlet pigeonhole.
  * Modulo a prime (`k <= 6`, small primes), every bad residue vector carries a repetition, a doubling or a Schur-type relation, of `l1`-norm `<= 4` and almost always `<= 3`. The test is informative only at `k = 3`, where most vectors carry no such relation (base rates, §3.2).

### 0.2 Is Collatz aligned with the primes?

**No at the gate primes. Yes at the archimedean place.**

* **Census (FINITE-EXACT).** All 780 gates `(p,a)` with `p <= 40` were counted exactly, with primitive words:
  * the local densities at the primes of `D` equal `1/q`: the ratio is 1.006 ± 0.115 wherever `L/q >= 5`;
  * empty zero classes match a necklace-Poisson model;
  * CRT gluing is independent in the bulk (1.08 ± 0.46).
* **Convergent gates up to `p = 1054` (FINITE-EXACT, float DP):** `|q nu_q - 1| <= 7·10^-7` at every small gate prime.
* **Every structured exception is explained** (§2.4):
  * cyclotomic: repeated words are automatically divisible by the cofactors `D(p,a)/D(p/d,a/d)`;
  * multiplicative order: a new **zero-carry identity** gives `N_5 = 0` for every shape `a = p-2` with `5 | D` (PROVED);
  * two isolated empty zero classes at *sub-gate* primes. These are EMPIRICAL, and an open micro-question.
* **No Hasse-type obstruction exists.** In 113 gates every prime power of `D` is locally solvable and the gate is not. In every one of them the local-product prediction is at most 0.353 primitive necklaces, so each is a failure of *size*, not of arithmetic.
* **The obstruction is archimedean.**
  * On the dyadic side the perigee bound alone removes 380 of the 497 dyadic gates with `p <= 40` (exactly those with `a < p/2`).
  * On the 3-adic side (`q = 3`) there is *no* perigee obstruction at all, and cycles exist there (−1, −5, −17).
* **What does matter** is the size `|2^p − 3^a|`: the continued fraction of `log_2 3` and Baker–Rhin. The prime factorisation of `D` does not.

### 0.3 The duality, typed (§3)

1. **Local-information asymmetry (PROVED).**
   * *LRC.* A counterexample is constrained at **every** prime `ell`: `V mod ell` must be "bad", meaning no lonely `m/ell`. The bad density is positive (the tight lines are always bad) and, at large `ell`, equals the tight-line density `~ c_k ell^{1-k}` (B5).
   * *Collatz.* A cycle word is constrained **only** at `ell | D` and at infinity.
     * For every `ell` not dividing `6D`, `x_w` is `ell`-integral.
     * At 2 the word *is* the local datum.
     * At 3, `x_w` is a unit: `c_w = 2^{s_{a-1}} (mod 3)`.
     * By CRT, "solvable at every `q^e || D` for the same `w`" is the same as global solvability.
   * This is the precise content of "LRC local, Collatz global". **The owner's asymmetry is CONFIRMED.**
   * The 2025–26 computer proofs of LRC for 8–15 runners are exactly sieves that are local at every prime, plus an archimedean bound (CITED, §1.5).
2. **Exact dictionary (PROVED; REAL, but it transfers nothing).**
   * A cycle of shape `(p,a)` exists ⟺ the carry set `W(p,a) = {c_w}` covers `|D|` ⟺ LRC's trivial witness `1/|D|` fails for `W(p,a)`.
   * The trivial lemma has force only when `|D| <= C(p,a)+1`. This *coverage regime* is finite. For `p <= 5000` it is exactly (FINITE-EXACT):

     `(4,2), (5,3), (8,5), (11,7), (16,10), (19,12), (27,17), (38,24), (46,29), (84,53)`.
   * The sporadic cycle {−17, ...} at `(11,7)` lives in this regime.
   * For every larger gate, Collatz asks whether a modulus `q > k+1` is covered. LRC theory never needs that information.
3. **Collatz objects belong to LRC's trivial class (PROVED).**
   * Carries are prime to 3, so the carry set is lonely at `t = 1/3`.
   * Elements of an integer cycle are prime to 3.
   * The runners `{2^s 3^r}` are prime to 5, and **`kappa(2,3) := sup_t inf_{s,r} ||2^s 3^r t|| = 1/5` exactly**. The extremisers are only `±1/5, ±2/5`.
   * More generally `kappa(S) = 1/P(S)`, where `P(S)` is the least prime not in `S` and `S` contains every prime below it.
   * LRC-hard sets must cover every `q <= k+1`, and a `{2,3}`-structured set never does.
4. **No transference theorem (PROVED for the decisive tests).** Neither "LRC for `{2^j 3^k}` modulo `D` ⟹ no cycle at `D`" nor its converse holds.
   * *Single-relation lemma:* one relation `xi` of the relation lattice `Lambda_q(V)` excludes every lonely `m/q` by itself iff `||xi||_1 = 1`, that is, coverage.
   * A Collatz cycle is a relation of norm `a >= 2` among the runners, so LRC tests cannot see it.
   * Conversely, the only "times" that certify the absence of a relation of norm `a` are huddles, with all runners near 0. The archimedean point and its `<2,3>`-orbit, which is the gates lane's structured spectrum, are the only ones known.
5. **Information budget (PROVED heuristic arithmetic).** At `a/p = log_3 2`:
   * a gate offers `log_2|D| = p − O(log p)` bits of local conditions (the `O(log p)` is Baker–Rhin);
   * these face `log_2 C(p,a) = 0.94996 p` bits of candidates;
   * the deficit is **0.050 bits per step**.

   LRC has an unbounded supply of primes, each worth about `log ell` bits.

### 0.4 Joint target (§5)

**The tight-line principle, TLP(k).**
* **Statement.** For all large primes `ell`, a residue speed vector `u in (F_ell^*)^k` has no lonely time `m/ell` iff it is, up to order and signs, a multiple of the reduction of a primitive set `T` with `M(T) <= 1/(k+1)`.
* **Evidence (FINITE-EXACT).**
  * `k = 3`: every prime `47 <= ell <= 199`. The bad set is exactly the line of `{1,2,3}`, of size `24(ell−1)`.
  * `k = 4`: every prime `137 <= ell <= 199`, two lines; the last exception found is `ell = 131`.
  * `k = 5`: `ell = 101` and `151`.
* **What it gives LRC.** One prime beyond the archimedean bound decides LRC(k). It also predicts the mechanism: a local obstruction forces short relations of rank `k−1`.
* **What it gives Collatz, as a contrast.** At gate primes nothing stabilises, because the local solution sets have density `1/q` whatever the cycles do (A2.1).
  * Cycles are visible only 2-adically, as *all* rational cycles (THM-4508 Lemma C).
  * Integrality is decided only at infinity.
  * So the Collatz analogue of TLP is the archimedean selection: side-aware perigee counting and Baker.

---

## 1. LRC side: what the primes decide

### 1.1 The covering lemma (PROVED, classical) and its converse layer

* **Lemma.** If no speed is divisible by `q <= k+1`, then `||v/q|| >= 1/q >= 1/(k+1)` for every `v`, so `t = 1/q` is lonely. This is THM-366/523 in the repo.
* Hence the class that needs proof is the divisor-complete class, which is how the LRC(14) finish map reduces the problem.
* Localised at a single denominator: `V` has a lonely time `m/q` with `q` prime and `q <= k+1` **iff** `q` divides no speed. Writing `Bad_q(k)` for the set of bad residue vectors:

  `Bad_q(k) = {u : some u_i = 0}`  for prime `q <= k+1`.

### 1.2 The hardest instances are not covering (A1.2, A1.4; FINITE-EXACT)

| k | primitive sets (max speed `<= B`) | tight sets (`sigma = 0`) | DC? | least covering `M` |
|---|---|---|---|---|
| 3 | B = 16 | `{1,2,3}` | no | 2/7 = 2/(2k+1) |
| 4 | B = 18 | `{1,2,3,4}`, `{1,3,4,7}` | no | 2/9 |
| 5 | B = 20 | `{1..5}`, `{1,3,4,5,9}` | no | 2/11 |
| 6 | B = 21 | `{1..6}` | no | 2/13 |
| 7 | B = 21 | `{1..7}`, `{1,2,3,4,5,7,12}`, `{1,4,5,6,7,11,13}` | no | 2/15 |

* **Tight sets.** Every tight set is lonely at `j/(k+1)` and nowhere with a smaller denominator. It covers every `q <= k` and misses `k+1`. Its `lld` equals its least uncovered integer, `k+1`, and no denominator fails without being covered.
  * The list agrees with the sporadic tight sets of Wills, `{1,3,4,7}`, `{1,3,4,5,9}` and `{1,4,5,6,7,11,13}`, and with the Goddyn–Wong multiples. See CITED: Zhang, arXiv:2608.13599, which classifies all one-entry modifications of `[n−1]`.
* **Consecutive rows** `{w..w+12}` with `w <= 400` (`k = 13`, A1.4):
  * the hardest row is `w = 1`, the AP: `sigma = 0`, not DC;
  * 371 of the 400 rows are DC, and among them the smallest slack is `sigma = 3/4`, at `w = 2`.
* **Random rows.** The medians of `sigma` for DC and non-DC rows are close (A1.5, last column).

**Verdict.** Prime coverage *raises* `M` above the tight value rather than lowering it. "Hard instances are exactly those whose speeds cover the small primes" is **REFUTED** in the extremal sense, and true only in the proof-theoretic sense of §1.1.

### 1.3 The least lonely denominator is the least uncovered integer (A1.5; EMPIRICAL, exact per row)

Each cell has 200 random primitive rows, seed 20260930. The table is copied verbatim from the `.out`.

| k | Vmax | DC % | P(lld = u) | P(llp = u_prime) | joint failures / row | ρ(depth, σ) | ρ(depth, lld) | median σ DC / non-DC |
|---|---|---|---|---|---|---|---|---|
| 8 | 18 | 20 | 0.870 | 0.735 | 0.250 | -0.19 | +0.56 | 0.636/0.862 |
| 8 | 60 | 13 | 0.925 | 0.845 | 0.105 | -0.15 | +0.59 | 1.204/1.323 |
| 8 | 1000 | 22 | 0.955 | 0.880 | 0.050 | +0.05 | +0.58 | 2.100/2.043 |
| 9 | 20 | 15 | 0.900 | 0.730 | 0.140 | -0.18 | +0.69 | 0.775/0.905 |
| 9 | 60 | 14 | 0.930 | 0.785 | 0.095 | -0.05 | +0.56 | 1.429/1.349 |
| 9 | 1000 | 19 | 0.945 | 0.845 | 0.070 | -0.11 | +0.64 | 2.127/2.129 |
| 10 | 22 | 15 | 0.910 | 0.820 | 0.205 | -0.25 | +0.35 | 0.804/0.958 |
| 10 | 60 | 14 | 0.955 | 0.940 | 0.050 | -0.10 | +0.55 | 1.444/1.391 |
| 10 | 1000 | 17 | 0.935 | 0.930 | 0.080 | +0.02 | +0.58 | 2.289/2.306 |
| 11 | 24 | 22 | 0.895 | 0.795 | 0.170 | -0.12 | +0.52 | 0.846/0.935 |
| 11 | 60 | 8 | 0.985 | 0.870 | 0.015 | -0.07 | +0.47 | 1.542/1.400 |
| 11 | 1000 | 14 | 0.955 | 0.900 | 0.060 | -0.14 | +0.53 | 2.414/2.385 |
| 12 | 26 | 6 | 0.970 | 0.920 | 0.045 | -0.25 | +0.40 | 0.857/1.000 |
| 12 | 60 | 9 | 0.970 | 0.960 | 0.035 | -0.19 | +0.43 | 1.432/1.404 |
| 12 | 1000 | 15 | 0.955 | 0.955 | 0.060 | +0.03 | +0.58 | 2.436/2.457 |
| 13 | 28 | 8 | 0.940 | 0.870 | 0.085 | -0.26 | +0.36 | 0.867/1.000 |
| 13 | 60 | 9 | 0.965 | 0.950 | 0.045 | -0.17 | +0.46 | 1.685/1.438 |
| 13 | 1000 | 12 | 0.955 | 0.955 | 0.045 | +0.01 | +0.48 | 2.510/2.552 |
| 14 | 30 | 9 | 0.940 | 0.855 | 0.105 | -0.10 | +0.28 | 1.010/1.069 |
| 14 | 60 | 8 | 0.965 | 0.905 | 0.045 | -0.02 | +0.44 | 1.500/1.500 |
| 14 | 1000 | 11 | 0.970 | 0.960 | 0.035 | +0.02 | +0.32 | 2.574/2.650 |

All 21 cells satisfy:
* `P(lld = u) >= 0.8`, and `>= 0.9` when `Vmax >= 60`;
* `|ρ(depth, σ)| <= 0.4`, while `ρ(depth, lld) > 0.2`.

Hence:
* **Prime coverage governs where a lonely time first appears. It does not govern how lonely the set is.**
* "Joint failures" are denominators `q < lld` that are covered by no speed and still fail. In random rows they are rare.

### 1.4 Deep wells and the LRC(14) rows (A1.3; FINITE-EXACT)

* **Deep wells** `{1..k−1, k(k+1)}`, for every `k = 4..14`:
  * `M = (k+1)/Phi6(k+1)` with a *unique* maximiser `(k+1)/Phi6(k+1)`;
  * but `lld = 2k+1`, because `t = 2/(2k+1)` is lonely (PROVED: `2v mod (2k+1)` avoids `0, ±1` for `v <= k−1`, and `2k(k+1) = k (mod 2k+1)`).
* **Correction to an OPEN-QUESTIONS line.**
  * The `{1..12,182}` paragraph (HYP-4047, opus S52) says the deep well is "lonely ONLY at 14/183 … the small-q census (q≤45) fails (best 1/15)".
  * Under the standard threshold 1/14 this is **false**. `t = 2/27` gives `min_v ||2v/27|| = 2/27 > 1/14`. Near `14/183` the lonely set is exactly `[(13+1/14)/182, (13+13/14)/182] ∪ [(14+1/14)/182, 13/168]`, about `[0.071821, 0.076531] ∪ [0.077315, 0.077381]`; a `10^-6` grid gives the span `[0.071822, 0.077380]` (A1.3).
  * What is true: the *maximiser* of `M` is unique.
  * See §6.2.
* **The three LRC(14) rows** of the bridges lane:

| row | M | DC | lld | least uncovered | llp | joint failures (q < lld, uncovered) |
|---|---|---|---|---|---|---|
| `{1..12,182}` | 14/183 | yes | 27 | 15 | 41 | 15–25 (11 values) |
| `{1..12,5460}` | 420/5461 | yes | 27 | 16 | 41 | 16,17,18,19,22,23,24,25 |
| `26·{1..12} ∪ {339}` | 1/13 (attained at 132 points of `(0,1/2]`) | yes | 27 | 15 | 41 | 15,17,19,21,23,25 |

* **Hardest primitive covering sets** (`k = 4..7`): `lld = 2k+1`, and 3–5 joint failures each.
* **The additive mechanism.**
  * An AP `{1..k−1}` contained in `V` blocks every denominator `q` in `(k+1, 2k)`. Indeed `q/(k+1) in (1,2]` forces `|v m mod q| >= 2`, that is `±m^{-1} notin [1, k−1]`, which needs `q >= 2k`. PROVED.
  * This is Dirichlet blocking, not divisibility.

### 1.5 Where the primes do come back in LRC (CITED + FINITE-EXACT)

* **`k+1` prime against composite** (A1.1, FINITE-EXACT).
  * The AP `{1..k}` has a lonely time with a *prime* denominator iff `k+1` is prime: `k = 10, 12` yes; `k = 7, 8, 9, 11, 13, 14` no.
  * So the repo's LRC(14), with `14 = 2·7`, cannot be certified by prime denominators.
  * The same split appears in the literature. Sungkawichai–Trakulthongchai (arXiv:2604.23906) use a polynomial argument when `k+1` and `p > k^2+k` are odd primes. Allikvere (arXiv:2609.02604) notes that this argument does not apply directly because 14 and 15 are composite; he finishes 14 by a direct search and 15 through `15 = 3·5`.
* **Computer proofs = local-everywhere sieves (CITED).**
  * Rosenfeld: 8 runners (arXiv:2509.14111) and 9 runners (arXiv:2512.01912).
  * Trakulthongchai: 9 and 10 runners (arXiv:2511.22427).
  * Sungkawichai–Trakulthongchai: 11, 12 and 13 runners (arXiv:2604.23906).
  * Allikvere: 14 and 15 runners (arXiv:2609.02604).
  * All of them combine a bound on a minimal primitive counterexample with exhaustive modular computations at many primes:
    * the bound comes from Tao (arXiv:1701.02048, speeds `n^{O(n^2)}`), from Malikiosis–Santos–Schymura (arXiv:2411.06903, Forum Math. Sigma 13 (2025) e164, speeds `<= binom(n+1,2)^{n−1}`), and from Allikvere's product bound;
    * the modular step is a covering search followed by lifting.
* **Tight-instance growth is a Jacobsthal phenomenon (CITED).**
  * Zhang (arXiv:2608.13599) shows that tight one-entry modifications satisfy `max V <= 0.60 n log n + 52 n`.
  * The constant 1/2 is sharp along `n = p# + 2` (primorials), and the proof uses the Jacobsthal function, which measures how the primes dividing `r` can cover runs of integers.
  * This is the one place where *extremal* LRC structure really is "aligned with the primes".
* **Transference to linear forms (CITED).**
  * Beck–Everett (arXiv:2609.06259): every counterexample or tight instance satisfies `m·V = 0` with `0 < ||m||_1 <= min(2k+3, (k+1)/(k−1) flt(k))`. Here `flt(k)` is Khinchine's flatness constant.
  * This is the Dirichlet/geometry-of-numbers side of the bridges lane's verdict.

---

## 2. Collatz side: local data at the gate primes

### 2.1 Classes (A2.1; q = 3; all 780 gates with p ≤ 40; exact counts of primitive words)

Each gate `(p,a)` is classified as follows. `L` is the number of primitive necklaces. A class marked "!" contains a local zero whose necklace-Poisson expectation `L/q^e` is at least 3.

| class | meaning | dyadic (`D > 0`) | 3-adic (`D < 0`) |
|---|---|---|---|
| GLOBAL | a primitive integral cycle | 1: (2,1) {1,2} | 2: (3,2) {−5,−7,−10}; (11,7) {−17,…} |
| REPEAT | only repetitions of these | 19 | 14 |
| ARCH | perigee: no word can be integral (`2^p > 4^a`) | 380 | 0 |
| PRIME-GATE | `|D|` a prime power: local = global, no solution | 14 | 32 |
| LOCAL | a proper prime-power divisor already has no solution | 19 | 180, plus 6 marked "!" |
| HASSE | every `q^e || D` solvable, the gate not | 64 | 49 |

* **ARCH.** The class is PROVED. Least point `<= 1/(2^{p/a} − 3) < 1` (Belaga). On the 3-adic side `c_min = 3^a − 2^a >= |D|`, so every word is eligible and there is no perigee obstruction for `q = 3`.
* **HASSE.** The largest local-product prediction is **0.353 primitive necklaces**, at `(19,12)` with `D = −7153 = −23·311`. Every Hasse-type failure is therefore a failure of size.

### 2.2 Local statistics (A2.1)

* **Prime types.** Each prime `q || D` (eligible gates, no primitive cycle) is one of:
  * *sub-gate*: it divides `D(p/d, a/d)` for some `d > 1`, `d | gcd(p,a)`;
  * *cofactor*: it divides `D(p,a)/D(p/d,a/d)`;
  * *primitive*: neither.

| prime type | pairs | empty zero class: observed / necklace-Poisson | `N_q / (Cp/q)` where `L/q >= 5`: mean ± sd |
|---|---|---|---|
| primitive | 542 | 277 / 263.1 | 1.006 ± 0.115 |
| cofactor | 279 | 79 / 81.7 | 0.998 ± 0.196 |
| sub-gate | 262 | 48 / 43.2 | 0.984 ± 0.201 |

* **CRT gluing.** The ratio of `N_m` to the CRT prediction `Cp·prod nu` is:
  * in the bulk (`0.5 <= a/p <= 0.8`, eligible, prediction `>= 20`): mean 1.081, sd 0.455 (n = 188);
  * over all gates: mean 1.020, sd 0.452 (n = 417).
* **Repetitions included.** If all words are counted, repetitions included, the z-scores against `C/q` reach **76.5**. The cause is the cyclotomic effect of §2.4(i), for example at `(40,22)`.

### 2.3 Verdict: arithmetic or archimedean?

* **Arithmetic (driven by `ord_q 2` and `ord_q 3`).** The local densities at the primes of `D` are `1/q` up to Poisson noise, after removing the explained structures of §2.4. They carry no global information:
  * the necklace-Poisson model predicts the empty zero classes to within one standard deviation for every prime type.
* **Archimedean.**
  * It decides most of the dyadic side outright: the perigee bound removes 380 of the 497 dyadic gates.
  * Elsewhere, size decides: no gate has local-product prediction above 0.353.
  * The uniform prediction `sum L/|D|` over the 116 eligible dyadic gates without a cycle is 2.442 primitive necklaces, and 0.561 after the local data are used. Over the 281 such 3-adic gates the figures are 0.929 and 0.513. The observed number is 0 in both.
* **The prime pattern that matters is `|2^p − 3^a|` itself.** Its size is governed by the convergents of `log_2 3` and by Baker–Rhin lower bounds. The coverage regime of §3.3 consists of convergents, semiconvergents and their doubles.

### 2.4 Structured exceptions, all explained or isolated

**(i) Cyclotomic repetitions (PROVED; A2.2, all 148 repetitions with `p <= 12` checked).**
* For a `d`-fold repetition,

  `c_{u^d} = c_u · (2^p − 3^a)/(2^{p/d} − 3^{a/d})`.

  So every prime of the cofactor divides `c_{u^d}` *for free*.
* Example: `2^40 − 3^22 = (2^20 − 3^11)(2^20 + 3^11)`, and every square word `u u` is `0 mod 2^20 + 3^11`.
* These local solutions never glue. Gluing them would need `D(p/d, a/d) | c_u`, which is a smaller gate's global problem.
* They are removed by counting primitive words only (Moebius over `d`). After that the z-anomaly disappears: the primitive ratio is 1.006 ± 0.115.

**(ii) Zero-carry identity (PROVED; A2.2, 8166 words checked).** Let the zeros of `w` be at `z_1 < ... < z_{p−a}`. Then

`c_w = −(2^p − 3^a) + sum_{j=1}^{p−a} 3^{a+j−1−z_j} 2^{z_j}`.

* *Proof.* Sum `(2/3)^s 3^{#zeros before s}` segment by segment, and use `3·(2/3) − 1 = 1`.
* Modulo any `q | D`, words with few zeros are therefore short sums of `{2,3}`-units.
* *Corollary (PROVED).* For `a = p − 2` and `5 | D` (which happens exactly for odd `p`), we have `2/3 = −1 (mod 5)` and hence

  `c_w = 3^a((−1)^{z_1} + 3(−1)^{z_2}) != 0 (mod 5)`.

  So `N_5 = 0`, which was checked for all 29 odd `p` in `[5, 61]`.
* This is an order-driven arithmetic obstruction (`ord_5(2/3) = 2`). It lives in the far 3-adic corner (`a/p -> 1`), where the archimedean count is already small. It accounts for four of the six "!" gates.

**(iii) Two isolated empty zero classes at sub-gate primes (EMPIRICAL; OPEN micro-question).**
* At `(36,32)` with `q = 263 | D(9,8) = −23·263`, the expected count is `L/q = 6.2`.
* At `(40,35)` with `q = 1931 = |D(8,7)|`, the expected count is `L/q = 8.5`.
* In both cases the residue histogram of `c_w mod q` has exactly one empty class, and it is 0.
* Under necklace-Poisson the probabilities are about `2·10^-3` and `2·10^-4`. Across the 262 sub-gate pairs the look-elsewhere expectation is about 0.05 for the second.
* For gates with `gcd(p,a) = 1` and few zeros, empty zero classes match the model (`z = 3..6`: 4/3.3, 4/2.5, 1/1.7, 0/0.4; scratch run). So the anomaly, if real, belongs to sub-gate primes, where the shape is a multiple of a smaller gate whose modulus is `q`.

### 2.5 Convergent and semiconvergent gates up to p = 1054 (A2.3)

| (p,a) | side | amp | log2 C − log2 |D| | small gate primes (ord_q 2, ord_q 3) | `|q nu_q − 1|` |
|---|---|---|---|---|---|
| (46,29) | dyadic | 40.5 | +0.01 | 39409 (9852, 4926) | — |
| (65,41) | dyadic | 87.7 | −0.08 | 19 (18,18), 29 (28,28), 17021 (17020, 3404) | `<= 7.2·10^-7` |
| (84,53) | 3-adic | 479 | +1.22 | 11 (10,5), 467 (466,233) | `<= 5.7·10^-11` |
| (149,94) | dyadic | 107 | −4.59 | 7^3 (3,6), 30809 (15404, 30808); fully factored | `2·10^-16` |
| (233,147) | dyadic | 138 | −8.76 | 5 (4,4) | 0 |
| (317,200) | dyadic | 193 | −12.70 | 2593 (81,648), 6791 (679,3395) | `<= 10^-16` |
| (401,253) | dyadic | 322 | −16.33 | 53 (52,52) | `2·10^-16` |
| (485,306) | dyadic | 979 | −19.07 | 929 (464,928); 123-digit composite cofactor | `10^-16` |
| (569,359) | 3-adic | 939 | −23.45 | 5, 582551; 143-digit composite cofactor | 0 |
| (1054,665) | 3-adic | 2.3·10^4 | −43.56 | 37 (36,18); 312-digit composite cofactor | `2·10^-16` |

* The factorisations are partial: trial division to `10^6` plus `sympy`'s rho and `p−1`.
* `log2 C − log2|D|` changes sign after `(84,53)` and then falls by about `0.05p`. That decline is the information-budget deficit of §3.6.

### 2.6 Control: 5x+1 (A2.5, p ≤ 22)

* The census finds the sporadic positive cycles at `(5,2)` (`D = 7`, five words) and `(7,3)` (`D = 3`, 14 words, two cycles), and the free cycle {−1,−2} at `(2,1)`.
* So the same machinery detects global solutions when they exist.

---

## 3. The duality, made precise (B)

### 3.1 Local-information asymmetry (PROVED)

**(a) LRC.**
* Let `V` be a counterexample to LRC(k). Then for **every** modulus `ell`, `V mod ell` is in `Bad_ell(k)`, the set of vectors with no lonely `m/ell`.
* For prime `ell <= k+1`, `Bad_ell(k) = {some u_i = 0}`.
* For prime `ell >= k+2` the bad density is a nontrivial codimension condition:
  * `ell·beta_ell` lies in `[0.2, 5.3]` for `k <= 6` and small `ell` (B3);
  * at large `ell` the density falls to `~c_k ell^{1−k}` (B5).
* So every prime carries information. The 2025–26 computer proofs are exactly sieves over many primes, closed by an archimedean (size) bound (§1.5).

**(b) Collatz.**
* For a word `w` of shape `(p,a)`, `x_w = c_w/D` is `ell`-integral for every prime `ell` not dividing `6D`.
* At `ell = 3`, `c_w = 2^{s_{a−1}} (mod 3)` and `D = 2^p (mod 3)`, so `x_w` is a 3-adic unit.
* At `ell = 2` the word is determined by `x mod 2^p` (Everett–Terras), so 2-adically every word occurs.
* The cycle condition is therefore local only at `ell | D` and at infinity.
* By CRT, the conditions at all `q^e || D` for one `w` together are the global condition. The only Hasse-type question is *uniformity in `w`* (gluing), and §2.2 answers it: in the bulk the events are CRT-independent (gluing ratio 1.08 ± 0.46).

**Type.** The owner's "local against global" is **CONFIRMED** as a theorem-level statement.

### 3.2 The relation lattice: Poisson, single relation, flatness, short relations (B1–B3)

* **Poisson identity (PROVED; checked exactly in B1).**

  `N_lone(V,q) = q · sum_{xi in Lambda_q(V)} prod_i ghat(xi_i)`,

  where `Lambda_q(V) = {xi : xi·V = 0 (mod q)}` and `g` is the indicator of `||x/q|| >= 1/(k+1)`.
  * LRC is therefore a *weighted sum over all relations*. The main term comes from `xi = 0`, and the corrections from short relations.
  * The Collatz count `N = |D|^{-1} sum_h S(h)` counts *one* family of relations: the 0/1 staircase vectors `epsilon_w` with `epsilon_w·R = 0 (mod D)`, where `R = (3^{a−1−i} 2^s)`.
* **Single-relation lemma (PROVED; B2, all vectors with entries in [−3,3], k = 3..8).**
  * The dual points `y = (m v_i/q mod 1)` satisfy `<xi, y> in Z`. On the lonely box `[delta, 1−delta]^k` the form `<xi, y>` ranges over an interval of length `(1−2 delta)||xi||_1`.
  * For `||xi||_1 = 1` that interval contains no integer.
  * For `||xi||_1 >= 2` and `k >= 3` its length is `>= 1`, so it does contain one.
  * Hence **coverage (`q | v_i`) is the only relation that blocks loneliness on its own**, and a Collatz cycle, a relation of norm `a`, never does.
* **Flatness transference (PROVED here from Khinchine's flatness theorem; CITED Beck–Everett for the real-time version).**
  * If `V mod ell` is bad, the lattice `Lambda_ell(V)^*` misses the box.
  * So the box has lattice width `<= flt(k)` in some direction `m in Lambda_ell(V)`.
  * Therefore `||m||_1 <= flt(k)/(1 − 2/(k+1)) = (k+1)/(k−1) · flt(k)`, the same constant as Beck–Everett's second bound.
* **Short relations (FINITE-EXACT; B3).** Every bad class for `k = 3..6`, at every prime up to 61/41/29/23, carries a relation of `l1`-norm `<= 4`.
  * The least norm is 2 or 3 for every class except 5 classes at `k = 4`, `ell = 41`, which need norm 4.
  * The relations are repetitions `u_i = ±u_j`, doublings `2u_i = ±u_j`, or Schur triples `u_i ± u_j ± u_l = 0`.
  * The base rate is the honest control (B3, "base rate").
    * For `k = 3` and `ell >= 37` only 33–52% of *all* vectors carry a norm-3 relation, yet all bad vectors do. The test is informative.
    * For `k >= 4` in the tested ranges every vector carries a relation of norm `<= 4`, so there the test is vacuous.
  * Local LRC badness is therefore additive (Schur/E3), which matches route [B] of the LRC(14) finish map.

### 3.3 The exact dictionary and the coverage regime (PROVED + FINITE-EXACT; A2.1, A2.4)

* **Dictionary.** A cycle of shape `(p,a)` exists ⟺ `|D|` divides some carry ⟺ the carry set `W(p,a)`, an LRC speed set with `C(p,a)` speeds, covers `|D|` ⟺ the trivial witness `t = 1/|D|` is *not* lonely for `W(p,a)`.
* **When it has force.** If `|D| <= C+1` and there is no cycle, the trivial lemma makes `t = 1/|D|` an LRC-lonely time for `W(p,a)`.
* **The coverage regime** `1 < |D| <= C(p,a)+1`:
  * for `p <= 5000` it is exactly `(4,2), (5,3), (8,5), (11,7), (16,10), (19,12), (27,17), (38,24), (46,29), (84,53)` (FINITE-EXACT);
  * beyond `p = 5000` it is empty, given an effective irrationality measure of `log_2 3` of moderate size (Rhin 1987; exact constants UNCITED-RECOLLECTION). This part is CONDITIONAL.
  * Of the ten, only `(11,7)` carries a cycle (the sporadic −17 cycle), and `(4,2)` carries the repeated {1,2}.
* **Type.** REAL as an exact reformulation; NO TRANSFER. The only LRC result it activates is the trivial lemma. For all other gates, the Collatz question ("is `|D| > k+1` covered?") is one about which LRC theory has nothing to say.

### 3.4 Collatz objects are LRC-trivial (PROVED; A1.6, B4)

| object | coverage | LRC status |
|---|---|---|
| carries `c_w` | never divisible by 3 (`c_w = 2^{s_{a−1}} mod 3`) | lonely at `t = 1/3`, level 1/3, any size |
| integer cycle elements | never divisible by 3 (a multiple of 3 has no odd-step preimage) | lonely at `t = 1/3` |
| runners `{2^s 3^r}` | never divisible by 5 | `kappa(2,3) = 1/5` exactly, extremisers `±1/5, ±2/5` only |

**Proof of `kappa(2,3) = 1/5`.**
* `{1,2,3,4}` is contained in `<2,3>` and is the tight AP for `k = 4`, so `kappa <= 1/5`, with equality only at `j/5`.
* `5` divides no `2^s 3^r`, so `kappa >= 1/5`.
* In general, `kappa(S) = 1/P(S)` when `S` contains every prime below `P(S)`.

**Comparison with BLMV.**
* Bourgain–Lindenstrauss–Michel–Venkatesh, *Some effective results for ×a×b*, Ergodic Theory Dynam. Systems 29 (2009) 1705–1722 (journal data UNCITED-RECOLLECTION).
* Their effective version of Furstenberg's theorem makes the orbit `{2^s 3^r m/q}` dense at an explicit rate, iterated-logarithmic in `q` as recalled (UNCITED-RECOLLECTION). So lonely times at a fixed level occur only at bounded denominators.
* The exact multiplicative LRC constant is trivial: 1/5.

**Type.** "Each lacks information about the other's type" is **CONFIRMED** in the strong form. The `{2,3}`-structured objects all sit in LRC's non-covering class.

### 3.5 Is there a transference theorem? (typed answer)

**"LRC-type loneliness of the multiplicative runner set modulo `D`" ⟹ "no cycle at `D`": NO (REFUTED as an implication).**
* Loneliness of the runners is generic. For the full semigroup it is decided by the prime 5, for every gate.
* Cycle gates and non-cycle gates do not differ in it. Decisive test: `kappa(2,3) = 1/5` does not depend on `D`.
* The single-relation lemma shows that a cycle relation of norm `a >= 2` cannot obstruct loneliness by itself.

**"No cycle at `D`" ⟹ loneliness: NO.** Non-coverage of `|D| > k+1` has no LRC consequence (§3.3).

**What does transfer is the dual statement.** A non-cycle certificate has to be a dual vector `(n, m)` for which

`sum_i (n_{i,s_i} + m 3^{a−1−i} 2^{s_i}/D)`

avoids `Z` for every word. A sufficient condition is a *huddle* whose signed sums never vanish: every runner phase within `1/(2a)` of 0. That is the Bohr-set, or simultaneous-approximation, side, which is the opposite of loneliness.
* The only certificates of this kind known are `m = 1` (the archimedean point) and its transports.
  * `m = 1` works when every necklace has a rotation with `0 < x_w < 1`: that is the perigee bound.
  * Its transports by `<2,3>`, through the rotation law, are exactly the gates lane's large spectrum, `h = u 2^{−j} 3^k` with `|u| <= 3`.
* So the structured spectrum of `S(h)` is the archimedean place seen through the rotation law. LRC has no analogue of it.

**Type.**
* REAL: the shared object (the relation lattice and its dual cloud), the Poisson identity, the dictionary of §3.3.
* ANALOGY: "LRC ↔ Collatz" as problems. LRC is point-in-box (Dirichlet). Collatz is vector-not-in-lattice (Liouville/Baker).
* This agrees with the bridges lane (`procgen_bridges_20260923`, B1/B5) and now has the explicit duality behind it.

### 3.6 Information budget (PROVED arithmetic; B4)

At the critical ratio `a/p = log_3 2 = 0.630930`:
* `log_2 C(p,a)/p -> H = 0.94996` bits;
* `log_2|D|/p -> 1`, with the loss `O(log p / p)` controlled by Baker–Rhin;
* so the CRT sieve at a gate has budget `log_2|D|` against `log_2 C`, a deficit of **0.05004 bits per step**.

Sample values: `(84,53)` has 76.3 against 75.1, `(485,306)` has 456.0 against 475.1, and `(1054,665)` has 996.0 against 1039.5.

So the local sieve at a gate is "exactly critical up to 5%". An LRC sieve gains about `log ell` bits per extra prime and is never critical.

---

## 4. Joint table: LRC(14) tools and their Collatz counterparts (C)

| tool | LRC(14) side | Collatz side | decisive test | label |
|---|---|---|---|---|
| covering / divisor-completeness | non-covering ⟹ witness `1/q` (THM-366/523) | cycle ⟺ the carry set covers `|D|` | §3.3 dictionary; coverage regime finite (10 gates, `p <= 5000`) | REAL dictionary, no transfer |
| covering-min / deep well | `M = n/Phi6(n)`, unique maximiser (THM-724/726, `k = 13`) | best approximations of `log_2 3` (convergent gates) | both are "most resonant" extremals. For `k = 4..7` the primitive covering minimum is `2/(2k+1) <` deep well (A1.2), so the cyclotomic form is not generic; nothing maps `Phi6(n)` to a gate | ANALOGY (cyclotomic values on both sides); a shared prime is NUMEROLOGY: every prime `>= 5` divides some `2^p − 3^a` |
| excised Bonferroni | localised defect at the pack clocks (bridges lane) | the drift's AM–GM defect is spread out (dimension 0.95) | bridges lane §3: defect localised (LRC) against spread (Collatz) | ANALOGY |
| density floors (THM-661) | covering-moment floor `mu_{1/7}(E) >= bar_k` | eligible-necklace density (side-aware perigee counting, gates lane) | both are size or measure bounds that close a case. Collatz's is heuristic (`sum L_elig/|D|`) | ANALOGY |
| Schur deficit / product sets (route [B], THM-730) | Schur deficit ⟹ `L > 0` | carries are staircase sums over the product set `{3^r 2^s}` | B3: local LRC badness is carried by Schur relations; Collatz runner sets are full of norm-3 relations (`2ρ = ρ'`) yet are lonely at 1/5 | REAL for LRC (local face of [B]); ANALOGY across |
| Weyl sums / soft Weyl bound (route [A]) | `Q_s = o(r^2)` on arc midpoints | `S(h)` off the structured set; square-root barrier (gates lane) | both need power saving off a structured set. On the Collatz side, the square-root barrier is a proved limit of second-moment methods (gates lane), and the structured set is the archimedean orbit | ANALOGY |
| min–max density game (THM-4447 capacity rule) | small clocks `{2,3,4}` are pigeonhole clocks | `{2,3}` are the Collatz primes | bridges lane S596: the clocks are `{2,3,4}` for every `n`, so they are pigeonhole clocks, not powers of 3 | NUMEROLOGY |
| Pillai / Gersonides (THM-4484, THM-4510) | — | free cycles ⟺ `|D| = 1` (`2−1`, `3−2`, `4−3`, `9−8`) | through the §3.3 dictionary a free gate (`|D| = 1`) is the carry set covering the modulus 1, which is automatic; Gersonides lists when it happens. LRC has nothing beyond the trivial lemma here | REAL through the dictionary, but trivial |
| Hankel / product formula (bridges B1) | Dirichlet/Siegel pigeonhole; transference (THM-2052, THM-4009) | Theorems R, D, Y; theta H1 | bridges lane AD1 | ANALOGY (LRC has no Liouville step) |
| sign strategies (THM-4508) | — | parity-graph cycles = 2-adic periodic points | Lemma C is a *2-adic stabilisation*: the closed walks at level `2^k` are exactly the rational cycles | ANALOGY to §5's tight lines (2-adic for Collatz, at large primes for LRC) |
| computer proofs modulo primes (new, CITED) | 8–15 runners: sieve at many primes plus an archimedean bound | cycle exclusion: archimedean (Baker, perigee) plus verified range, no prime sieve | §3.1: auxiliary primes carry no Collatz information | the asymmetry is REAL (PROVED) |

---

## 5. Joint target (D): local stabilisation, the tight-line principle

**TLP(k).**
* **Statement.** There is `ell_1(k)` such that for every prime `ell >= ell_1(k)`, a vector `u in (F_ell^*)^k` has no lonely time `m/ell` (for every `m` some `||m u_i/ell|| < 1/(k+1)`) **iff** `u = λ·(T mod ell)` up to order and the signs of coordinates, for some `λ in F_ell^*` and some primitive `k`-set `T` with `M(T) <= 1/(k+1)`.
* **Counting form `TLP*(k)`, assuming LRC(k).** `|Bad_ell(k)| = #(tight k-sets) · k! · 2^{k−1} (ell−1)`.

**Evidence (B5 and the scratch runs; FINITE-EXACT).**

| k | stabilised primes (bad set = tight lines exactly) | exceptions found | tight lines |
|---|---|---|---|
| 3 | every prime in [47, 199], and 31 | 29, 37, 41, 43 | `{1,2,3}`: `|Bad| = 24(ell−1)` |
| 4 | every prime in [137, 199]; also 59, 79, 83, 97, 107–127 | 61–73, 89, 101, 103, 131 | `{1,2,3,4}`, `{1,3,4,7}` |
| 5 | 101, 151 | 53–71, 89, 127 | `{1..5}`, `{1,3,4,5,9}` |

* The exceptions are reductions of sets with *short lonely intervals*, such as `{1,3,4,5}` (`M = 2/9`), `{1,x,x+1,x+2}` and `{1,2,3,5,8}`. Their intervals fall between the points of the grid `(1/ell)Z`.
* They thin out as `ell` grows. One necessary ingredient is in the literature: by Giri–Kravitz (arXiv:2304.01462) the accumulation points of the `k`-speed spectrum are the `(k−1)`-speed values. Given LRC(k−1), those are `>= 1/k > 1/(k+1)`, so `1/(k+1)` is isolated from above.
* That is not sufficient. A set with large speeds has a short lonely interval even when `M` is far above `1/(k+1)`. TLP asserts that, modulo a large prime, such sets never stay bad.
* Sampling at `ell = 211` (`k = 4`) and `ell = 101, 151` (`k = 5`) found only tight-line vectors.

**Status.** OPEN.
* PROVED: the "if" direction for tight `T` and primes `ell > (k+1) max T`. A lonely time `t` of a tight set has `||t v|| = 1/(k+1)` for some speed `v`, so its denominator divides `(k+1)v`; hence no lonely time of `T` has denominator `ell`, and every `λ·(T mod ell)` is bad. For the AP the lonely times are exactly `j/(k+1)`, so its line is bad at every prime `ell > k+1`.
* PROVED: the reduction below.
* FINITE-EXACT: the table.

**What a proof gives LRC (PROVED reduction).**
* Let `S_k` bound the speeds of a minimal primitive counterexample, for example `binom(k+1,2)^{k−1}` (Malikiosis–Santos–Schymura).
* Let `ell > max(ell_1(k), 2 S_k max_T max T)` be a prime with `TLP*(k)`.
* Then LRC(k) holds. A counterexample `V` would be bad mod `ell`, hence `V = λT` mod `ell`. Then `v_i t_j = ±v_j t_i (mod ell)` with both sides below `ell` in size, so `V ∝ T` over `Z`, which contradicts `M(T) = 1/(k+1)`.
* So LRC(k) is decided by *one* sufficiently large prime: the strongest form of "LRC is local".
* It also predicts the mechanism: local badness at large `ell` forces a rank-`(k−1)` lattice of short relations. That joins Beck–Everett's single relation (rank 1) to the repo's Schur/E3 route.

**What it gives Collatz (the mirror; PROVED + FINITE-EXACT).** Collatz has no analogue at primes. At the primes of `D` the local solution sets have density `1/q` whether or not a cycle exists (A2.1). Cycle and non-cycle gates look alike locally. The Collatz stabilisation happens in two other places:
* 2-adically: THM-4508 Lemma C says that at level `2^k` the closed walks are exactly the rational cycles `c_w/D`, *all* of them;
* archimedeanly: integrality is selected at infinity, by side-aware perigee counting (gates lane) and by Baker.

The transferable lesson for Collatz is directional. The "tight lines" of Collatz are the rational cycles, and the open problem is purely the archimedean selection among them. Gate-prime methods are provably blind: the census, and §3.1(b).

**Why this target, and not others.**
* It is strictly smaller than LRC(k): it needs the archimedean bound to conclude.
* Every instance of it is finite and checkable.
* It turns the owner's slogan into a dichotomy that can be falsified: *global extremal points are visible locally* (LRC, TLP) against *global points are invisible locally* (Collatz, census).
* Runner-up targets:
  * a mod-`ell` Beck–Everett with constant 4: supported only at `k = 3` (B3 base rates);
  * power-saving equidistribution at gate primes, uniform in `q <= C^{1/2}`: it would close off local methods for Collatz, but it says nothing about LRC.

---

## 6. Flags for the owner and for audit

1. **LRC(14) may be settled in the literature (CITED, not verified).**
   * J. Allikvere, *Fourteen and fifteen lonely runners*, arXiv:2609.02604, posted 2026-09-02. The abstract announces a proof of the conjecture for fourteen and fifteen runners.
   * The method combines a stronger bound on the speed product of a primitive counterexample with modular covering searches at many primes, binary lifting, and a direct search for 14. The code and certificates are archived publicly, according to the abstract.
   * The convention matches the repo: `k` speeds with threshold `1/(k+1)` are "`k+1` runners", as in Sungkawichai–Trakulthongchai for `k in {10,11,12}`. So "fourteen runners" is the repo's LRC(14), with 13 speeds and threshold 1/14.
   * Also new since the repo's map:
     * Sungkawichai–Trakulthongchai, 11–13 runners, arXiv:2604.23906;
     * Beck–Everett, relation theorem, arXiv:2609.06259;
     * Zhang, tight one-entry classification, arXiv:2608.13599;
     * Cordella, odd denominators in the six-speed spectrum, arXiv:2609.03444.
   * **Recommendation:** a dedicated audit lane should read 2609.02604, and its archived certificates if feasible, before any further LRC(14) work.
2. **OPEN-QUESTIONS, HYP-4047 line, about `{1..12,182}`.**
   * The text says "lonely ONLY at 14/183 … best 1/15". This is false at threshold 1/14: `t = 2/27` is lonely, and so is an interval of width about 0.0047.
   * True statement: `14/183` is the *unique maximiser* of `min_v ||t v||`. This lane did not edit the file.
3. **Covering-minimum pattern at small k (FINITE-EXACT, k <= 7, speeds <= 21).**
   * The primitive covering minimum there is `2/(2k+1)`, which is *below* the deep-well value `(k+1)/Phi6(k+1)` for `k = 4..7`.
   * This does not contradict THM-724/726, which are specific to `k = 13`. A 150 s local search at `k = 13` found nothing below `14/183`; its best was `1/13`.
   * It does mean that "deep well = covering minimum" is not a uniform-in-`k` phenomenon, and the `k = 13` proof must use something special about 13 and 14.
4. **Sub-gate zero-class anomalies** at `(36,32)` mod 263 and `(40,35)` mod 1931 (§2.4(iii)) are the only unexplained local events in the census.

## 7. Reproduction

* Run `nice python3 -u 04-computation/experiments/procgen_localglobal_20260930_run.py`.
* It takes about 6 min (345 s, 29198 checks). Memory stays far below the 500 MB cap; the heaviest step, an exact `p = 40` gate count, peaked near 130 MB when measured alone.
* Its stdout is `05-knowledge/results/procgen_localglobal_20260930.out`, which ends `ALL CHECKS PASSED`.
* Exploratory scripts live in `scratch/procgen_localglobal/` (untracked):
  * `explore_*`, `analyze_census*`, `subgate.py`, `fewzeros.py`;
  * `tightline*.py`: the stabilisation runs to `ell = 199`;
  * `srlp_sample.py`, `covmin13.py`.
* arXiv metadata was fetched through the export API with a generic User-Agent; no personal data was sent.
