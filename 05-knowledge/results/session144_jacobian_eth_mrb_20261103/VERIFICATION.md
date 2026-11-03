# Session144 verification and trust boundaries

## Frozen checkpoints

- Checkpoint1 commit: `b687cd28cac8eaa3ac9c0b96b959c1ccd2e89aa5`. Its frozen [Jacobian manifest](jacobian/MANIFEST.json) and all science files are unchanged by checkpoint2.
- Checkpoint2: [manifest](MANIFEST.json), [replay program](verify_checkpoint2.py), and [synthesis](SYNTHESIS.md). The mutable packet overview and navigation are outside the manifests. Manifest files themselves are byte-pinned by Git publication metadata rather than self-hashed.
- The original [MRB manifest](mrb/MANIFEST.json) is preserved exactly, as are all its 26 listed files. The additional independent analysis and rational-precision receipt are separate.

From this packet directory:

```
python3 -B jacobian/verify_frozen_packet.py
python3 -B -O jacobian/verify_frozen_packet.py
python3 -B verify_checkpoint2.py
python3 -B -O verify_checkpoint2.py
```

The second program validates all frozen hashes before copying inputs into a temporary directory. It replays the certificates there, so no frozen output is overwritten. Normal and optimized integration receipts are identical after the explicitly defined normalization of rational-precision timing/mode metadata. Python, SymPy and mpmath are needed; the exact rational precision evaluator itself uses only the standard library.

## Executed scope

- Jacobian checkpoint1: 37 producer checks, 46 independent checks, seven rejected corruptions, and independent enumeration/reconstruction of all 1,799 targets over F5/F7/F11.
- Jacobian checkpoint2: seven direct padded quotient-determinant cases, 35 equivariant-gauge checks, 30 coordinate-fiber reconstruction checks, and 13 quartic boundary checks with four rejected corruptions. All universal conclusions rely on the written proofs.
- ETH: every one of 1,434 CNFs in the stated finite encoding universe, an exact eight-vertex correction/permanent experiment, independent permutation evaluation at five points, exact interpolation and undersized-grid hostiles. The large imported algorithm and proof-assistant dependency graph were not rebuilt.
- MRB: 2,046 binary periodic masks, 40 unit masks, six certified nonzero perturbation-shift intervals, five hostile mutations, and eight hostile precision inputs. Fully rational p=1,8,16,32 certificates match independent reruns after removing only elapsed time and interpreter mode. No p=64 completion is claimed.
- Source provenance: all eight ETH input files independently match the pinned upstream raw bytes and both SHA256/Git blob hashes. The initial local extra-newline materialization and canonical correction are explicitly recorded in [SOURCES](eth/SOURCES.json).

Historical inherited MRB programs are preserved for provenance. Their Python assertions are not used as optimized-mode acceptance gates. Numerical illustrations and scalar quadrature checks are labelled separately from exact rational or directed-interval certificates.

## Documentation and publication scope

Only owned files are overlaid onto the exact parent tree; no repository index is used. API payloads are individually below 85,000 bytes. The final returned tree must match the independently built Merkle tree, and publication checks the expected parent before moving main.

Scoped documentation checks cover all owned file hashes, conflict markers, relative links and unchanged startup source pins. The full repository checker was not run because it reads the restricted central mistakes ledger. That ledger was neither read nor edited. The mandatory startup surface remains unchanged at 179,983 of its 180,000-byte budget. No empty reserved claim was introduced.
