# Partition meet projectors: signed versus physical

THEOREM.md gives self-contained proofs, explicit Pi_3 projectors, the all-finite-semilattice probability idempotence classification, and the exact CRT/overlap and symmetry boundaries. Frozen inputs are copied and hash-pinned under sources/. No source packet or remote was changed.

Reproduce:

    python3 verify.py > certificate.replay.json
    cmp certificate.json certificate.replay.json
    python3 -O verify.py > certificate.optimized.replay.json
    cmp certificate.json certificate.optimized.replay.json
    python3 -O verify.py --self-test-failure

The last command intentionally exits nonzero with RuntimeError. All earlier commands must succeed. Certificates use exact integer and rational arithmetic. Finite checks supplement, rather than replace, the all-size proof.

Core result: every nonnegative real idempotent under independent meet convolution on a finite meet-semilattice is zero or Dirac. Arbitrary signed idempotents instead correspond to 0/1 upper-zeta vectors. Independent resampling and correlated repetition are different operations; neither signed projectors nor naive normalization legitimizes overlapping physical alternatives.
