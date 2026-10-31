# Session 125: rooted forests and even-Fibonacci numeration

**PROVED explicit bijections and determinant identities; FINITE-EXACT
implementation controls.** This packet gives a real object-level bridge:

`legal CG words of length L <-> tilings of 2L+1 cells <-> rooted forests of P_(L+1)`.

It preserves **nonzero digit count = domino count = forest edge count**, not
digit sum. Leading zeros retain the same integer but add rooted singleton
components, so fixed length must be retained. Direction is also part of the
state: the forest decoder reads LSD-first, unlike an MSB append-digit state.

Start with [PROOFS.md](PROOFS.md), then the independent
[proof review](audit/PROOF_REVIEW.md). Classical Chung-Graham and matrix-forest
results are credited there; no historical novelty claim is made.
The companion [numeration packet](../session125_origin32_numeration_20261031/numeration/PROOFS.md)
contains the arithmetic language proofs. Runtime checks here are standalone
and use Python's standard library only; no external source data is required.

Right-appending a zero appends a rooted singleton and realizes the exact
Beatty shift `n -> ceil(((3+sqrt(5))/2)*n)`. The one-edge canonical forests
have ranks `A_k` or `2A_k`; these ranks are not total path-forest counts.

The named-port Green matrix gives the exact next rooted-forest count.
For a leaf birth on P3, an endpoint gives 21 while the middle gives 20,
despite the same initial graph/count and arrival degree. This boundary state
is not the connectivity/endpoint state needed for Hamiltonian continuation.
**Neither CG digit normalization nor forest counts prove Hamiltonicity.**
The positive G32 origin and Gerbicz threshold are separate arithmetic and
graph-theoretic statements in the companion origin packet.

## Reproduction

From this packet directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python test_forest_bridge.py verification.normal.json
PYTHONDONTWRITEBYTECODE=1 python -O test_forest_bridge.py verification.optimized.json
cmp verification.normal.json verification.optimized.json
PYTHONDONTWRITEBYTECODE=1 python audit/independent_audit.py audit/independent.normal.json
PYTHONDONTWRITEBYTECODE=1 python -O audit/independent_audit.py audit/independent.optimized.json
cmp audit/independent.normal.json audit/independent.optimized.json
python verify_packet.py
python -O verify_packet.py
```

The producer performs **1,245,858** checked conditions per run; the independent
auditor performs **129,144**, importing no producer implementation. Both
include the word-201 reading-direction regression and right-singleton/Beatty
operation controls. Full JSON receipts are
retained in normal and optimized forms, byte-identical within each pair.

`verify_packet.py` verifies the exact publication whitelist, all file hashes,
source bindings in both receipt families, and normal/optimized equality.
`PUBLICATION_FILES.json` is the whitelist. `SHA256SUMS` covers every published
file except itself; `PUBLICATION_MANIFEST.json` records all content hashes
except its own and that checksum file to avoid a circular hash definition.
Generated `__pycache__` files are not part of the packet; the commands above
prevent them. Relocation replay and intentional corruption rejection were
checked before publication; the parent session retains that release receipt.
