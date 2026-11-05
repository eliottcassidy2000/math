# Independent audit: complete dyadic fiber and same-order C1 replacement

Verdict: PASS. A provenance path reversal in DEPENDENCIES.md was corrected
before final release; all mathematical/code/receipt bytes stayed unchanged.
Final manifest SHA256: dbd103b4943af67b1db7e30d2a0d708b6c883edaabfc41facedab0d9e78721df.

- PROOF.md SHA256: c627cc79107440d50133c63593c16bcf5e9a70cb48f389a5b9f6f0f833e713a5
- verify.py SHA256: d9c35f84bfb10673cb5f2d558cda2d4d89d1106465533ae552c1b6078af1796d
- Both receipts SHA256: 160fd79c40fe95772f0ac36395713846fce79fb426fde4458d3d82882c4e145d

The full-fiber inverse is valid. Every weak start-one source has a canonical
strong/singleton core and a singleton SCC tail. A strong core's ratio U/D lies
in [1,2), strictly below2. If the full ratio is2^k, dividing by its tail factor
2^r forces r=k and a balanced core. This proves both directions on actual
graphs, not only on profiles or selected words.

The canonical balanced-source inverse then proves the entire stated ratio
fiber W=J(O^(k+1)(T)). Its core T is actual, balanced or singleton, and has
Hamilton count e. The singleton boundary is handled correctly.

For the residual replacement, b=2^(k+1) and h=(b+1)e give
(2b-1)h-(b-2)e=(2b^2+1)e. The alternative sink-tail construction has exactly
that count. Each word adds2k+4 vertices, so order is preserved. Both graphs
are strong C1. The old graph has no balanced root by the reviewed residual
dichotomy, while the new graph does, proving nonisomorphism. This is a scalar
and order-preserving replacement, not a false all-edge identity.

Isolated normal/optimized runs reproduce both receipts. All dependency bytes
match the previously reviewed packets. Independent controls in
../dyadic-independent.py and ../DYADIC_INDEPENDENT.json reproduce:
- all32,548 raw source profiles under U2000 and exactly3,126 dyadic profiles;
- all416 canonical tail recoveries through order5, with every reconstructed edge;
- all144 balanced-ratio recoveries through order5, also checking every edge
  under the full canonical root/sink/tail map, rather than only its profile;
- all15 literal graph-pair controls through order15 using a separate forward
  subset recurrence, including exclusion of every balanced root on the old
  side and the required balanced root on the replacement side.

No mathematical correction was needed. The wrapper correctly notes that these
scalar values already had sink-tail witnesses; new content is the complete
physical entry condition and same-order residual rebase. The result does not
settle arbitrary residual inputs, the balanced cofinite lemma, or odd coverage.
