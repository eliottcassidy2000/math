# Exact receiver-only route thresholds

The original nineteen-slot triple-equality receiver has sharp simultaneous event budgets **16 slots / 34 clocks**, both for ordinary connectivity and strict resource non-reacquisition, uniformly over all exact grades and charges with exterior records frozen.

The central grades distinguish the notions: grade 9 or 10 requires only 14/30 ordinarily but 16/34 strictly. A retained grade-9 pair has individual optimum 8/18 ordinarily versus 16/34 strictly.

Read `THEOREM.md` for scope, proof and sharp witnesses. `certificate.json` retains every primitive trade and exact grade statistics. `fibres.tsv` retains all 23444 fibre summaries. The checker recomputes all 131072 states and exact ordinary/strict optima for all 474896 endpoint pairs, with a frozen independently matched streaming digest.

Reproduce:

    python -B reproduce.py
    python -B verify_manifest.py --relocate

The packet is manifest-closed, standard-library-only and includes the frozen original native carrier and full native/dyadic authentication. It does not claim unrestricted-helper or all-size wiring thresholds.
