# Perfect-core twin-prime aliquot edges

19 October 2026. Status: **PROVED**, independently audited, with bounded exact controls.

Read [the proof](PERFECT_CORE_TWIN_EDGES.md) and [independent review](review/INDEPENDENT_REVIEW.md). For a fixed positive integer core, perfectness is exactly the guarded +2 translation case. Twin pairs with lower prime at least 5 correspond bijectively to actual aliquot edges 6p -> 6(p+2). Each individual fixed-6 prime lane exits its admissible class within two steps. The source and target prime/coprimality guards remain explicit. This is an elementary specialization of the inherited guarded identity; it proves no prime-infinitude or general aliquot-convergence claim.

## Reproduction

Use Python 3 and its standard library, from this package directory:

    python3 -B verify_perfect_core_edges.py
    python3 -B -O verify_perfect_core_edges.py
    python3 -B review/independent_controls.py
    python3 -B -O review/independent_controls.py

Compare the first two outputs with `verification.out` and `verification_optimized.out`, and the last two with their corresponding review receipts. All checks use explicit failures and remain effective under optimization. The author checks 901 guarded edges. The independent controls check shift rigidity for 1,000 positive cores, four specified perfect cores, and 2,260 fixed-6 prime starts through lower prime 20,000: 1,919 exit after one step and 341 after two. These finite regressions support the proofs without extending their scope.

`SOURCE_PROVENANCE.json` records the reviewed source hash, the two scope-only wording clarifications, the approved final source hash, and published artifact hashes. Publication adds relative dependency links, this reproduction guide, and review-path clarification; mathematical tests are unchanged.
