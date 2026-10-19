# Rational quadratic cycles: denominator-first certificates

Research date and primary-source check: **2026-10-19**.

## Main conclusions

- Every rational parameter `c` with reduced denominator at most **1,000,000**, with no bound on its numerator, has no rational cycle of `x²+c` longer than three. This is a finite exact certificate using the proved Eliahou–Fares denominator and period theorems. Their cited theorem already covers the infinite family with at most two denominator primes; the remaining 101 denominators comprise 96 three-prime and five four-prime cases. The complete census yields 62 oriented cycles up to cyclic rotation, all of period three.
- Given the prime factorization of an arbitrary point denominator `b`, all cycles of length at least three can be enumerated across every parameter numerator in `F(ω(b))·poly(log b)` bit operations. The universal two-edge method has at most `4(2N−1)^ω(b)` signed seeds, `N=2^ω(b)+2`, and an explicit escape guard controls integer sizes. Factoring `b` is outside this complexity assertion.
- A simple-resonance valuation lemma excludes noncentral periodic points in a specified p-adic disk. At `c=−29/16`, every shell in the three-cycle tube modulo `43^k` has exactly 21 cycles of the explicitly proved period. The old period-six cycles modulo `43²` grow to period 258 modulo `43³`; different, closer period-six cycles remain at that larger modulus.
- Standard critical-cycle certificates give infinite parameter families: `c∈Z_5`, `c≡1 mod5` allows only period three over `Q_5`; `c∈Z_3`, `c≡2 mod3` allows only period two over `Q_3`.

The full Flynn–Poonen–Schaefer conjecture remains **open**. There is no claimed world priority, unbounded-denominator exclusion, or inference of rationality from local compatibility. The real parabolic three-cycle at `c=−7/4` and the rational three-cycle at `c=−29/16` are explicitly distinguished.

## Read first

1. [Proof, literature boundary, and exact theorem scope](session7-rational-cycles/denominator_signature_theorem.md).
2. [Independent mathematical and implementation audit](session7-cycle-review/INDEPENDENT_REVIEW.md).
3. [Complete primary 101-denominator output](session7-rational-cycles/denominators_through_1000.json).
4. [Independent comparisons](session7-cycle-review/comparison_receipt.json).

The two sibling directories preserve the original script imports and comparison paths. All computations use Python's standard library and exact integer/rational arithmetic. The independent audit imports no author implementation. Runtimes in independent outputs are descriptive and are not mathematical data.

## Quick copied-artifact replay

Run from this result directory. These commands compare all recorded exact cycles, rerun the efficient all-101-denominator solver, and rerun the finite local controls:

```sh
python3 session7-cycle-review/compare_receipts.py
python3 session7-rational-cycles/rational_cycle_pair_probe.py --max-b 1000 --compare session7-rational-cycles/denominators_through_1000.json --output /tmp/rational_pair_replay.json
python3 session7-rational-cycles/rational_cycle_pair_probe.py --factorization 2:2 3:20 5:20 --output /tmp/rational_large_exponent_replay.json
python3 session7-rational-cycles/critical_cycle_controls.py
python3 session7-rational-cycles/simple_resonance_controls.py
python3 session7-cycle-review/check_general_resonance.py
python3 session7-cycle-review/check_43_lift.py
```

For a full independent regeneration, including all five four-prime cases:

```sh
python3 session7-cycle-review/check_words.py
python3 session7-cycle-review/independent_search.py
python3 -O session7-cycle-review/independent_search.py --output session7-cycle-review/independent_results_optimized.json
python3 session7-cycle-review/word_projected_search.py
python3 -O session7-cycle-review/word_projected_search.py --output session7-cycle-review/projected_results_optimized.json
python3 session7-cycle-review/universal_pair_search.py
python3 session7-cycle-review/compare_receipts.py
```

The complete primary three-edge regeneration, slower than the two-edge version, is:

```sh
python3 session7-rational-cycles/rational_cycle_signature_probe.py --max-b 1000 --max-omega 4 --output /tmp/rational_primary_replay.json
```

The supplementary author-written `rational_cycle_signature_independent.py` checks actual primitive-word tables through `b=400`; its `b≤1000` output explicitly skips the five four-prime searches while checking their three-cycle parametrizations. It is not used as the complete independent certificate. The separate audit's projected-word and universal-pair runs both check all 101 denominators.
