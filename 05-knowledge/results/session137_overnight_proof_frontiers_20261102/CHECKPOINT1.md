# Checkpoint1 verification

The main proof was reviewed independently through determinant cycle expansion, matrix positivity and exact vertex-cycle identities. The global audit pins the final proof SHA256f7e79e9aa935533f2a0444ce2bf6b220f5ff1cb29fea16adc1fa19a7cc71ba6f. A further fresh referee and a separate7,412-class cycle census corroborate it; those additional artifacts will accompany the final session packet.

## Reproduce root checks

    python3 -B discriminant/check_global.py > root.normal.json
    python3 -B -O discriminant/check_global.py > root.optimized.json
    cmp root.normal.json discriminant/global.normal.json
    cmp root.optimized.json discriminant/global.optimized.json

The helper probe.py supplies exact integer Bareiss determinants. Python checks raise exceptions and stay active under optimization. Four explicit hostile controls retain the empty-domain exception, distinguished-cycle sign, OCF factor2, and weighted one-vertex paths.

## Independent checks

See independent_global/GLOBAL_DETERMINANT_AUDIT.md for exact C++ and Python commands. C++ Gaussian elimination uses a prime and a proved Hadamard reconstruction bound; a UBSan build repeats the ordinary census. The rational weighted route is separate from the producer and includes its complete equality locus. Proper-cycle deficits are checked separately.

All proof and receipt hashes are checked before publication, and published blob bytes are compared against the intended files. Repository documentation and route checks are separate from mathematical validation. No CI-success claim follows from local checks or absent GitHub statuses.

The six updated current surfaces are THM1950, THM3729, THM1860, CURRENT-FRONTIER, PROBLEM-LEDGER, and ACTIVE-GUARDRAILS. THM1860 is a status repair only: H>=1+2c3 was already an immediate inherited OCF corollary, not a new theorem. Historical reductions remain labelled as historical; their valid product laws are retained. Unrelated repository content is preserved.

A fresh isolated packet completed14 principal executions: normal/optimized Python routes and ordinary/UBSan C++ routes, including the certificate implementation-under-test audit. Every output matched its frozen receipt. [Portable replay record](REPLAY_CHECKPOINT1.json).
