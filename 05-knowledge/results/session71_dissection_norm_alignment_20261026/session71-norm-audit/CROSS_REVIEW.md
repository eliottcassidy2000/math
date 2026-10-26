# Exact-byte producer cross-review

Producer: `/workspace/scratch/07d5532f2fce/session71-norms`.

Accepted final files:

| File | SHA-256 |
|---|---|
| `REPORT.md` | `a1a83aba7162a076a169320d79a28d14f41dd4cbe09171c3d223f36ab2ed5efd` |
| `check.py` | `f7e1a568f90073a4b60ab5aecd01981941eacc4c8dde3415e0ac36a6756c6c11` |
| `SOURCE_PINS.json` | `f264554358e8d639b226f34015dfc000f2b381455f6f99c0525c3b77258d16b6` |

Independent replay commands:

```sh
python3 -B ../session71-norms/check.py > producer-normal.json
python3 -O -B ../session71-norms/check.py > producer-optimized.json
cmp producer-normal.json producer-optimized.json
```

Both PASS, **824,939 explicit checks**, byte-identical; receipt SHA-256 `acef4b41d34572c303f078a080e3d105a769307def567e2358f3acab776a3af0`.

Mathematical verdict: **accepted in its stated scopes**. The rational norm criterion, Gaussian identities, raw-versus-normalized Farey hostiles, real coordinate-field obstruction, and similarity counterexample are correct. Historical44/5 candidates are explicitly non-identified, not falsely unified. The paper's count arithmetic is typed as a reconstruction of lists, not an independent completeness proof.

The recovered **PS-D2** proof was read directly in producer source snapshot10. For ordinary polygonal regions (closures of their interiors), its supporting-line argument is valid. Every exterior square-side segment contributed by a complete piece edge lies wholly along that side; two such edges from distinct regular pieces cannot overlap in positive length, because their local interiors would occupy the same inward side. Thus each target side is a sum of complete piece-edge lengths. Each square side of the normalized Dudeney target requires an edge length outside Q(√3). Four distinct such exposed edges do **not** imply four pieces. This inherited field-of-lengths mechanism is more directly relevant than an unsupported Gaussian-minimum isomorphism.

Two review requests were fixed before acceptance:

1. Replace “unrestricted four-piece optimum” with the actual no-flip polygonal category.
2. Scope the boundary sum to regular polygonal regions, excluding dangling zero-area appendages. The final draft also distinguishes edge lengths from vertex coordinates and four exposed edges from four pieces.

The producer's minimal **real Cartesian** coordinate degree4 is consistent with this audit's minimal **complex vertex-coordinate** degree8. The latter requires i in the field; the former does not. Both concern target realizability, not the coordinates of any common dissection's cuts or a piece-count lower bound.

No producer file was edited by this audit. The independent checker and additional rank theorem are in this lane's `check_norm_rank.py` and `REPORT.md`.
