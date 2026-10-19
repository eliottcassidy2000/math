# Deletion–Frobenius exactness

A proof-level extension of the earlier characteristic-two tournament padding complex.

- All integer tournament statistics have a canonical unimodular basis of source-free induced-subtournament counts. In that basis every subset-deletion operator is diagonal, giving its complete integral cokernel and Smith torsion.
- Every prime p gives a p-contractible deletion complex, with explicit preimages for every kernel ker(d^k)=im(d^(p-k)).
- Every Hasse digit D_(p^e) has its own exact p-contraction with source-adjunction h^(p^e).
- Lucas carries give an exact finite transducer for nonzero deletion compositions; higher Hasse digits recover information erased by ordinary iteration.
- Mod-p spanning embeddings see exactly the hosts whose automorphism order is prime to p. The directed triangle is invisible mod3; higher Hasse deletion escapes this restricted observable space.

Proof: [PADDING_FROBENIUS_PROOF.md](PADDING_FROBENIUS_PROOF.md).

Reproduce: `python3 verify_padding_frobenius.py` and `python3 -O verify_padding_frobenius.py`. The two supplied receipts match on 46,839 exact checks. Independent review is separate.

Scope: elementary identities and exact reconstruction, not a Rule-30 prize, Gilbreath proof, good-tree classification, or a claimed novelty result for general N-complexes.

## Independent review and exact replay

The [independent mathematical review](review/PADDING_PROOF_REVIEW.md) passes the final proof, including all grading, bottom-degree, saturation, and Smith-form boundaries. Its separate class generator and matrix checker are included. The session's coordinating reader also checked the integral projection, cokernel, rank, and string arguments independently before approving this freeze.

From the repository root, these commands reproduce and compare every receipt:

```bash
cd 05-knowledge/results/tournament_hasse_frobenius_20261019
python3 -B verify_padding_frobenius.py > /tmp/padding-frobenius-replay.json
python3 -B -O verify_padding_frobenius.py > /tmp/padding-frobenius-replay-optimized.json
cmp receipt.json /tmp/padding-frobenius-replay.json
cmp receipt-optimized.json /tmp/padding-frobenius-replay-optimized.json
python3 -B review/verify_padding_independent.py > /tmp/padding-frobenius-review.json
python3 -B -O review/verify_padding_independent.py > /tmp/padding-frobenius-review-optimized.json
cmp review/padding_independent_receipt.json /tmp/padding-frobenius-review.json
cmp review/padding_independent_receipt_optimized.json /tmp/padding-frobenius-review-optimized.json
```

The producer uses only the Python standard library. The independent checker uses SymPy; the frozen run used SymPy 1.14.0. SHA-256 pins are in `SHA256SUMS` and the publication manifest.
