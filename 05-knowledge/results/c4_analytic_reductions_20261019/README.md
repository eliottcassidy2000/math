# Frozen analytic c=4 coverage, 2026-10-19

This packet proves sufficient reductions for seven-thin mixed-(1,2), c=4 cores. It does not prove the complete c=4 exclusion.

## Exact scopes

1. The one-defect exchange lemma bounds p(u<v) by 2/3 when the directed predecessor/successor defect is at most one, and by 1/2 for zero defects. The two-sided defect/slack criterion gives an actual ambient balanced pair on 49,592 of 79,121 degree-valid cores.
2. Under the separate maximum-index-span strict-majority hypothesis, every common vertex lies between the roots. A one-sided defect/slack contradiction covers 54,235 cores. The union is 56,630 cores, with 22,491 residual cores.

The combined compressed manifest tags each chosen certificate as AMBIENT_BALANCED_PAIR or MAXIMUM_SPAN_MAJORITY_CONTRADICTION. The first certificate is preferred where both kinds apply. Every core record includes its canonical index and predecessor masks. The canonical degree-valid catalogue hash is recorded, making union with residual cone coverage exact.

## Reproduction (Python standard library only)

From this extracted directory:

```sh
python -O verify_saturated_near_twin_lemma_20261019.py
python -O verify_one_defect_coverage_20261019.py
python -O build_combined_c4_analytic_cover_20261019.py
```

The first command checks the specialized exchange maps and regenerates its 31,568-core coverage. The second checks the general injections on all 5,231 naturally labeled posets through six vertices: 3,863,727 extension/pair controls, 37,169 zero-defect ordered pairs, 12,873 extra-successor pairs, and 12,873 extra-predecessor pairs. It independently regenerates the 49,592-core global criterion. The last command verifies those exact provenance hashes and builds the separately scoped combined coverage without repeating the injection audit.

All checks use explicit exceptions and remain active under Python optimization. Receipts and compressed manifests are deterministic and were replayed under Python 3.12.14. The proof is the written injection argument; the finite audits are controls and the c=4 counts are exhaustive finite morphology enumerations. No conjecture counterexample or historical novelty claim is made.

`manifest.json` records SHA-256 and size for every included file except itself. Originals remain unchanged outside this frozen directory. The separate cone scan can reuse this analytic manifest but its status is not included in this packet.
