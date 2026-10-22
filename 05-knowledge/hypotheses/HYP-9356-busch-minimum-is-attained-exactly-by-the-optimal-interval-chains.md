---
id: HYP-9356
title: "Busch's minimum is attained exactly by the optimal interval chains: for every m >= 3 the strong m-tournaments with f(m) = min 3^a 5^b (2a + 3b = m - 1) Hamiltonian paths are exactly the interval chains IC(p_1, ..., p_k) - the transitive order 0 -> 1 -> ... -> m-1 with one back arc from the last to the first vertex of each block [l_i, l_i + p_i], consecutive blocks sharing an end vertex (Moon's chains S_n) - whose parts lie in {2, 3, 4} with prod (2^(p_i - 1) + 1) = f(m). Hence E(3q) = q, E(3q+1) = 1, E(3q+2) = q(q+3)/2 classes; the self-converse ones are the palindromic compositions; all are rigid, with scores (1,1,2,3,...,m-2,m-2) and m - 2 cyclic triples"
status: >
  KNOWN, not open. The characterization is J. W. Moon and L. L. M. Yang, "Strong tournaments with the fewest Hamiltonian
  paths", Utilitas Math. 72 (2007) 193-198 (CITED via El Sahili and Ghazo Hanna, arXiv:2101.00713, J. Graph Theory 2023:
  Moon and Yang constructed the "special chains" linking nearly transitive tournaments, proved they attain the minimum,
  and proved they are the only tournaments attaining it). The primary paper was not accessed.
  The chains, their product formula prod (1 + 2^(p_i - 1)) and the minimizing patterns are J. W. Moon, "The minimum
  number of spanning paths in a strong tournament", Publ. Math. Debrecen 19 (1972) 101-104, p. 102, as read by the S50
  audit. THM-1370 already recorded that Moon gave the chain construction.
  This file was filed as an open hypothesis by S50 before the literature was checked (MISTAKE-652).
  What S50 adds:
  - an independent finite confirmation through order 13: E(12) = 4 and E(13) = 1 were committed as predictions in
    fad6a40581 and passed (strong_ears.py 375, strong_ears_625.py);
  - an independent written proof of Moon's product formula and of the recognition and rigidity of chains (THM-4744);
  - an independent, computer-assisted reproof of the uniqueness theorem (THM-4745 (5), (6); pre-audit at the time of
    writing).
  The count formula follows from the characterization together with THM-4744 (2); whether Moon and Yang state it is not
  known here.
source: opus-2026-10-21-S50, 05-knowledge/results/h_spectrum_multiplicities_20261021.md (section 2)
related:
  - 01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md
  - 01-canon/theorems/THM-4745-exact-gluings-of-interval-chains-and-the-reduction-of-busch-equality-to-a-one-way-lemma.md
  - 05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md (S15: E(7) = 1 with scores (1,1,2,3,4,5,5))
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (Moon's chain construction mentioned)
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears.py (+ strong_ears_225.out, strong_ears_375.out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears_625.py (+ .out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/exact_gluings.py (+ .out)
---

# HYP-9356 — the extremal strong tournaments of Busch's theorem (known: Moon–Yang 2007)

## Statement (a known theorem)

* Busch (*EJC* 13 (2006) N3) proves that a strong `m`-tournament has at least `f(m)` Hamiltonian paths:
  * `f(3q) = 3·5^(q−1)`;
  * `f(3q+1) = 5^q`;
  * `f(3q+2) = 9·5^(q−1)`.

  Moon (1972) had shown that his chains `S_n` attain this value.
* **The chains.** For `m − 1 = p_1 + … + p_k` with every `p_i ≥ 2`, the chain `IC(p)` is the transitive order `0 → 1 → … → m−1` with one back arc per block. Each block is nearly transitive and has `2^(p_i − 1) + 1` paths.
* **Moon–Yang (2007), as reported.** The strong `m`-tournaments with exactly `f(m)` paths are exactly the optimal chains. (The secondary source does not define Moon–Yang's "special chains"; identifying them with the optimal `IC(p)` is an inference, consistent with Moon 1972 and with THM-4745 (6).)
  * `m = 3q`: one 2 and `q − 1` 3s;
  * `m = 3q + 1`: `q` 3s;
  * `m = 3q + 2`: two 2s and `q − 1` 3s, or one 4 and `q − 1` 3s.
* **Count and symmetry.** Different compositions give non-isomorphic chains (THM-4744 (2)). So the number of extremal classes is `E(3q) = q`, `E(3q+1) = 1`, `E(3q+2) = q(q+3)/2`, and the self-converse ones are the palindromes.

## How S50 met it

* The statement was found from the census at orders up to 9, before the literature was checked, and filed here as an open hypothesis.
* Predictions were then committed (`fad6a40581`) and tested. They are kept as an independent confirmation.

| prediction | result |
|---|---|
| `E(12) = 4`: the chains `(2,3,3,3)`, `(3,2,3,3)`, `(3,3,2,3)`, `(3,3,3,2)`, none self-converse | passed (`strong_ears_375.out`) |
| the order-12 extremal classes are rigid | passed |
| no strong 12-tournament below 375; strong spectrum below 375 misses only 1, 7, 21 | passed |
| `E(13) = 1`, the chain `(3,3,3,3)`, self-converse | passed (`strong_ears_625.out`) |
| `E(14) = 14`, `E(15) = 5`, `E(16) = 1`, `E(17) = 20` | not computed; implied by Moon–Yang and THM-4744 (2) |

* **Caveat.** Git cannot certify the order of the order-12 run. `strong_ears.py` already supported the 375 cap in `fad6a40581`, although its output was first committed afterwards. `strong_ears_625.py` did not exist at `fad6a40581`.
* The independent reproof of the uniqueness theorem is in THM-4745 (5), (6).
