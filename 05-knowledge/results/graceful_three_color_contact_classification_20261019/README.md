# Three-color contact classification

Packet 101. Exact finite, computer-assisted lemma: among three rank-at-most-two contact supports on four components, the only nonempty disconnected width-two fibers are cyclic two-choice supports on one spanning tree. Two independent standard-library checks cover all 15,625 ordered support triples, identifying exactly 96 disconnected cases.

Replay from this directory:

    python classify_rank2.py > /tmp/rank2-original.json
    cmp classification.json /tmp/rank2-original.json
    python verify_rank2_independent.py > /tmp/rank2-independent.json
    cmp independent_rank2_receipt.json /tmp/rank2-independent.json
    sha256sum -c SHA256SUMS

All runtime dependencies are co-located. The [rank-three criterion](../graceful_three_bad_rank_contacts_20261019/) handles the remaining higher-rank cases; it is contextual, not a dependency of the finite check. Arithmetic realizability is a separate question, resolved by a subsequent m=12 counterexample.
