# Reproduce and verify Session71

Run from this packet directory:

```sh
python3 -I -B verify_packet.py --replay
python3 -I -B -O verify_packet.py --replay
python3 -I -B verify_provenance.py
python3 -I -B test_packet_corruption.py .
```

Dependencies used: Python 3.12.14; SymPy 1.14.0 for the independent norm/rank
and golden-area controls; Matplotlib 3.10.8 only for regenerating the original
illustration. The remaining mathematical checkers use the standard library.
The drawing configures a writable temporary cache; its floating-point display
is not mathematical evidence.

## Frozen bytes and provenance

`FILE_MANIFEST.json` freezes every payload file's size, SHA-256 and mode.
`NATIVE_COPY_MANIFEST.json` records 107 byte-preserving native copies. All 39
first-checkpoint native files remain unchanged. Producer and auditor manifests
are copied intact; their declared scope is not silently enlarged.

`PROVENANCE_RULES.json` actively checks native hash lists/manifests, 25 inert
repository-source snapshots, exact approved source pins and 125 reviewed byte
bindings. Git origin paths were independently compared against the pinned
commits before publication. The sole full external source snapshot is the
unchanged licensed v3 HTML identified in `SOURCE_LICENSES.md`; other paper
PDFs/HTML, renderings and OCR remain outside the packet.

The source-version ledger distinguishes requested v3, the v4 correction notice
and the ITCS2026 publication. The external minimum-four theorem remains CITED
for polygonal pieces without flips; no full independent case-proof certification
is claimed.

## Portable replays

The envelope runs eleven declared jobs in a temporary copy, each normally and
under `python -O`, and compares stdout byte-for-byte with stored expectations:

1. Arithmetic/Gaussian/Farey controls (824,939 explicit checks).
2. Independent norm/labelled-rank audit (38,815 checks).
3. Independent golden-template area controls (five checks).
4. Active provenance verifier (1,247 checks in this final packet).
5. Dudeney abstract partial-affine chains, charge arithmetic and transcribed lists.
6. Native Penrose/involution/placement checker.
7. Penrose lane verifier, including its own normal/optimized and hostile gates.
8. Independent matching checker (6,186,654 checks; 60,398 involution pairs).
9. Independent geometric checker (468,013 checks; 1,364 placements).
10. Matching-audit byte-integrity verifier and four actual copied-byte corruptions.
11. Drawing regeneration; final SVG/PNG bytes must remain unchanged.

Seven additional semantic negative controls are required to fail with their
intended messages in **both** Python modes: the producer's injected guard fault,
invalid involution, erased overlap, product-only faithfulness, endpoint-color
loss, wrong rotation determinant, and staged-band overlap. These intentionally
false assertions are distinct from actual file-corruption tests.

The portable verifier rejects modified, added or deleted frozen inputs after
replay. It never executes source snapshots. The separate package mutation test
checks fourteen baseline/corruption cases in both modes (28 checks), including
size/hash/mode changes, missing/extra files, symlinks, path escape, duplicate
inventory, wrong source pin/schema, and a broken reviewed-byte binding.

## Publication verification

The durable session journal freezes each local commit's bundle, patch, tar,
Git blobs, modes, tree and parent. GitHub publication uses an expected-head
lease and `force=false`; concurrent main commits are preserved. Remote read-back
checks whole-tree equality, every changed path and all packet bytes/modes,
required ancestry, clean local/main alignment, and exact-commit CI queries.
Empty CI is not a pass. The source HTML and deterministic Matplotlib SVG each have a precise
`.gitattributes` whitespace/newline exemption to preserve reviewed original
bytes; all other changed paths are subject to the strict whitespace gate.
