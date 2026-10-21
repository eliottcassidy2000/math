---
id: HYP-9356
title: "Busch's minimum is attained exactly by the optimal interval chains: for every m >= 3 the strong m-tournaments with f(m) = min 3^a 5^b (2a + 3b = m - 1) Hamiltonian paths are exactly the interval chains IC(p_1, ..., p_k) - the transitive order 0 -> 1 -> ... -> m-1 with one back arc from the last to the first vertex of each block [l_i, l_i + p_i], consecutive blocks sharing an end vertex - whose parts lie in {2, 3, 4} with prod (2^(p_i - 1) + 1) = f(m). Hence E(3q) = q, E(3q+1) = 1, E(3q+2) = q(q+3)/2 classes; the self-converse ones are the palindromic compositions; all are rigid, with scores (1,1,2,3,...,m-2,m-2) and m - 2 cyclic triples"
status: >
  OPEN for m >= 13. Pre-registered predictions 1-3 PASSED (strong_ears.py 375: E(12) = 4, the four chains, none
  self-converse; no strong 12-tournament below 375; strong spectrum below 375 misses only 1, 7, 21).
  REDUCED (THM-4745, PROVED) to a one-way lemma: every exact gluing of two extremal chains of compatible orders has
  all cross arcs in one direction. That lemma is FINITE-EXACT for small chain pairs (exact_gluings.py).
  Through order 12 the classification is verified directly. FINITE-EXACT for m = 3..11: the census for m <= 9 (h_multiplicities.py) and exact ears for m = 10, 11
  (strong_ears.py 225, complete by Moon's strong-deletion theorem and the insertion lemma of THM-1370). The extremal
  classes are exactly the optimal chains, distinct optimal compositions give distinct classes, and self-converse
  is equivalent to palindromic. PROVED for all m (THM-4744, to be filed): every interval chain has exactly
  prod (2^(p_i - 1) + 1) Hamiltonian paths, so every optimal chain attains f(m). The count formula was read off
  m <= 9 and then confirmed at m = 10, 11 in the same session; those two orders were computed but not committed
  as predictions beforehand. The predictions below ARE committed before computation.
source: opus-2026-10-21-S50, 05-knowledge/results/h_spectrum_multiplicities_20261021.md (to be written)
related:
  - 05-knowledge/results/spectrum_basins_20261020e/spectrum.md (Busch's minimum, CITED there; the class series)
  - 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md (concurrent: fixed-H multiplicities; uses no equality characterization)
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (Moon's monotonicity, the insertion lemma)
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears.py (+ strong_ears_225.out)
---

# HYP-9356 — the extremal strong tournaments of Busch's theorem

## Statement

* Busch (EJC 13 (2006) N3, CITED) proves that a strong tournament on `m ≥ 3` vertices has at least `f(m)` Hamiltonian paths:
  * `f(3q) = 3·5^(q−1)`;
  * `f(3q+1) = 5^q`;
  * `f(3q+2) = 9·5^(q−1)`.

  He shows that Moon's construction attains this bound. He does not characterize or count the tournaments attaining it.
* **Interval chain.** For a composition `m − 1 = p_1 + … + p_k` with every `p_i ≥ 2`, the chain `IC(p_1, …, p_k)` has vertices `0, …, m−1`.
  * Put `l_1 = 0` and `l_(i+1) = l_i + p_i`.
  * Every arc goes forward except the back arcs `l_(i+1) → l_i`, one per block.
  * Block `i` alone is the nearly transitive tournament on `p_i + 1` vertices, which has `2^(p_i − 1) + 1` paths.
* **Conjecture.**
  * The strong `m`-tournaments with exactly `f(m)` Hamiltonian paths are exactly the chains `IC(p)` with parts in `{2, 3, 4}` and `∏ (2^(p_i − 1) + 1) = f(m)`. These are the optimal compositions:
    * `m = 3q`: one 2 and `q − 1` 3s;
    * `m = 3q + 1`: `q` 3s;
    * `m = 3q + 2`: two 2s and `q − 1` 3s, or one 4 and `q − 1` 3s.
  * Different compositions give non-isomorphic chains.
  * So the number of extremal classes is `E(3q) = q`, `E(3q+1) = 1` and `E(3q+2) = q(q+3)/2`.
  * The converse reverses the composition, so the self-converse extremal classes are the palindromes.
  * Every extremal class of order `m ≥ 4` has trivial automorphism group, score sequence `(1, 1, 2, 3, …, m−3, m−2, m−2)` and exactly `m − 2` cyclic triples. This is the minimum for a strong tournament (Moon).

## Evidence

| `m` | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|
| `f(m)` | 3 | 5 | 9 | 15 | 25 | 45 | 75 | 125 | 225 |
| `E(m)` computed | 1 | 1 | 2 | 2 | 1 | 5 | 3 | 1 | 9 |
| self-converse | 1 | 1 | 2 | 0 | 1 | 1 | 1 | 1 | 3 |

* Each computed extremal class was identified with an optimal chain by canonical forms.
* `m ≤ 9` comes from the exhaustive class census.
* `m = 10, 11` comes from exact ears.
  * By Moon, every strong `m`-tournament has a vertex `v` with `T − v` strong.
  * By the insertion lemma, `H(T − v) ≤ H(T)`.
  * So the ears of the strong `(m − 1)`-classes with `H ≤ f(m)` contain every strong `m`-class with `H = f(m)`.

## Why it should hold

Busch's Lemma 1 gives the decomposition behind the conjecture.
* A strong tournament of order at least 5 is nearly transitive, or it is the union of two strong subtournaments `A`, `B` with `|A ∩ B| = 1`.
* His proof gives `H(T) ≥ H(A)·H(B)`. So at equality:
  * both `A` and `B` are extremal;
  * the product `f(|A|) f(|B|)` equals `f(m)`;
  * the gluing is exact.
* The conjecture says that two optimal chains glue exactly only end to start. A nearly transitive tournament is extremal only for `m ≤ 5`, where it is the one-block chain.

## Pre-registered predictions

These were committed before the order-12 computation (`strong_ears.py 375`).
1. `E(12) = 4`: the compositions `(2,3,3,3)`, `(3,2,3,3)`, `(3,3,2,3)`, `(3,3,3,2)`. None of them is self-converse, since they form two converse pairs.
2. The four order-12 extremal classes are those chains, and they are rigid.
3. As a consequence of exact ears at the 375 cap: no strong 12-tournament has `H < 375`, and the strong spectrum below 375 omits exactly 1, 7 and 21.
4. For later orders: `E(13) = 1` (self-converse), `E(14) = 14` (2 self-converse), `E(15) = 5` (1 self-converse), `E(16) = 1`, `E(17) = 20`.

## What would refute it

Any of the following would refute it:
* an extremal strong tournament that is not an interval chain;
* two optimal compositions with isomorphic chains;
* `E(12) ≠ 4`.

## Results of the pre-registered predictions (computed after `fad6a40581`)

1. **Passed.** `E(12) = 4`: the strong 12-classes with `H = 375` are exactly the chains `(2,3,3,3)`, `(3,2,3,3)`, `(3,3,2,3)`, `(3,3,3,2)`. None is self-converse (`strong_ears_375.out`).
2. **Passed.** They are rigid.
3. **Passed.** No strong 12-tournament has `H < 375`, and the odd `d < 375` with `P_d = 0` are exactly 1, 7, 21.
4. **Untested.** The predictions `E(13) = 1`, `E(14) = 14`, `E(15) = 5`, `E(16) = 1`, `E(17) = 20` stand.

## Reduction (THM-4745)

* By Busch's Lemma 1, an extremal tournament of order `m ≥ 6` glues two extremal strong tournaments `A`, `B` exactly at one vertex, with `f(|A|)f(|B|) = f(m)`.
* THM-4745 (2) proves that a one-way exact gluing of two chains is the end-to-start concatenation.
* So the hypothesis for all orders follows from the **one-way lemma**: an exact gluing of two extremal chains has all its cross arcs in one direction.
* The lemma is exhaustive for chain pairs of at most 6 vertices with `(|A|−1)(|B|−1) ≤ 16`, and for `C_3`, `R_4` with the optimal 7-chain `IC(3,3)`. This covers every induction step up to order 10.
* `gluing_failure_modes.py` shows that a proof must use both the restriction condition and the merge condition.
