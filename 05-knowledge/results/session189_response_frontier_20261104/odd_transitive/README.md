# Odd transitive mixed-quartic tomography

This is a separately reviewable new packet. It proves an exact binomial mixed-quartic identity for every odd transitive module and a root-family extension under the earlier ordered-witness/equal-log-concave-profile hypotheses. Its explicit even-order residual is not claimed nonzero in general.

- `PROOF.md`: general two-channel identity; exact transitive kernel; odd cancellation; even boundary; scoped fixed-ambient original-module root-family theorem.
- `prepare.py`, `controls.json`, `graphs.txt`:43 physical module/exterior cases and89 named source, local-converse and quotient-only-converse graphs. The three witness cases include a fully covered original-module family in a strong fixed ambient.
- `bank.cpp`, `bank.part0.json` through `bank.part4.json`:133,382 complete named response coefficients and their supporting actual small-support query counts, through13 old names. Every named cubic, quartic, and zero is retained. The parts are consecutive graph-record lists and reassemble losslessly by concatenating those lists.
- `verify.py`, `inner_coefficients.json`, `binomial_controls.json`, `physical_receipt.json`: independent literal inner two-/three-cover enumeration, finite-difference kernel checks, exact exterior support algebra, contracted capacities and28 independent joint three-cover coordinate checks.
- `query_controls.py`, `count_queries.cpp`, `physical_queries.json`, `query_input.txt`, `query_output.txt`:672 actual grouped finite-difference query counts, with full physical masks, including prescribed nonzero exterior masks in the strong ambient. Source, target quotient-only, and local-converse values remain separate.
- `receipt.json`, `receipt.optimized.json`, `query_receipt.json`, `query_receipt.optimized.json`: identical normal/optimized results.

The strict grouped contrasts are176 in the TT5/TT7 three-block cases with TT3 witness children, and240 in the TT5 strong-ambient case with TT2 witness children. Genuine third-excursion weights are36 or12, so cancellation is not an artifact of an excursion-free control.

Reproduce with existing standard tools:

    python prepare.py
    g++ -O3 -std=c++17 bank.cpp -o /tmp/s189_odd_bank
    /tmp/s189_odd_bank regenerated.json < graphs.txt
    python verify.py split regenerated.json
    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
    cmp receipt.json normal.json
    python query_controls.py prepare
    g++ -O3 -std=c++17 count_queries.cpp -o /tmp/s189_odd_queries
    /tmp/s189_odd_queries < query_input.txt > query_output.txt
    python query_controls.py > queries.normal.json
    python -O query_controls.py > queries.optimized.json
    cmp queries.normal.json queries.optimized.json
    cmp query_receipt.json queries.normal.json

The verifier's split step regenerates the five canonical complete-data parts and their receipt hashes. No packages are installed. The two-port graph has at most15 vertices in these bank controls; every Hamilton count and finite-difference coefficient fits signed64-bit arithmetic.

The all-order proofs do not depend on extrapolating the finite checks. The root-family theorem relies on the separately reviewed Session189 distance reconstruction and the inherited compressed-anchor and exterior-skew theorems. Independent review of this new packet remains required before promotion.
