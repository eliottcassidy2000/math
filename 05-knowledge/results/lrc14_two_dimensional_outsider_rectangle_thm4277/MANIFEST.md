# THM-4277 promotion manifest

## Status and exact scope

This is the final finite-exact promotion packet for the reserved theorem
`THM-4277-uniform-two-dimensional-outsider-rectangle-common-deck-closure.md`.
It proves the displayed inequality for all 2,550 pairs
`450<=q<=499, 600<=r<=650` and all `binom(30,9)=14,307,150` labelled bodies
in the fixed thirty-label pool. It does not prove a neighbouring pair,
monotonicity, physical entry, or LRC(14).

The proposed canonical roots are

```text
SRC=04-computation/lrc14_two_dimensional_outsider_rectangle_thm4277
OUT=05-knowledge/results/lrc14_two_dimensional_outsider_rectangle_thm4277
```

The `repo-ready/` subtree in this packet materializes those exact destinations.

## Load-bearing and audit artifacts

All hashes are SHA-256 of raw LF bytes. C++ semantic transcripts have only the
nondeterministic line beginning `SECONDS ` removed. No mathematical line is
normalized or rewritten.

| canonical-relative path under `SRC` or `OUT` | SHA-256 |
|---|---|
| `SRC/lrc_rectangle_common.hpp` | `b0764f8e66fb7909755eb8ff3f56806d4cb776ec501ff5b168865892dcd4de4b` |
| `SRC/primary_endpoint.cpp` | `371da2c67b9079c254ece3f6e37508c5c250db4ca7649ec4a0d231a1d761a163` |
| `OUT/primary_endpoint_O3.out` | `33992e5f8e0a92a9a04c6974e432c0bf3e47cdde47d95c122824fc10942574d0` |
| `OUT/primary_endpoint_O3_NDEBUG.out` | `33992e5f8e0a92a9a04c6974e432c0bf3e47cdde47d95c122824fc10942574d0` |
| `OUT/primary_endpoint_UBSAN.out` | `33992e5f8e0a92a9a04c6974e432c0bf3e47cdde47d95c122824fc10942574d0` |
| `OUT/primary_UBSAN.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `SRC/independent_literal.cpp` | `edab6ee879681fc2cbbf249067ddb8ed169efbd9d8734fab9e4b384dd2e7f42c` |
| `OUT/independent_literal_O3.out` | `eb6e3c5f8f96c6eb47f0940c9e9ccdafc5ccc3349d1d634faa49350b2ed27f30` |
| `OUT/independent_literal_O3_NDEBUG.out` | `eb6e3c5f8f96c6eb47f0940c9e9ccdafc5ccc3349d1d634faa49350b2ed27f30` |
| `OUT/independent_literal_UBSAN.out` | `eb6e3c5f8f96c6eb47f0940c9e9ccdafc5ccc3349d1d634faa49350b2ed27f30` |
| `OUT/independent_UBSAN.stderr` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `SRC/verify_transcripts.py` | `d8a6e31c4f653ada95e43a068736e0d7d773ce57a839f03dcd5ffe41e7ccd70d` |
| `OUT/verify_transcripts.out` | `9b23fbddb2225f9afb0579357ecfb1d612e816a8cd1d78652cfb8c036a8d86ba` |
| `SRC/u256_comparator_audit.py` | `6343d641803c0f4b9c05bf518997755f4404d63d490ff041551184549ad5d98c` |
| `OUT/u256_comparator_audit.out` | `d7985b78d3b8f7ec9769c6f3e8ec7760c0025d776b7b490889fc99a8f704f48b` |
| `SRC/postprocess_current.py` | `3a7c8aa6c79d26d391125725b470c1b455353ef1a994b023c80d29959561fe95` |
| `OUT/postprocess_current.out` | `aed1a6329115fdb8da5ad511e221653c814b2517c202aa5e0dbb2083819b09da` |
| `OUT/REPRODUCTION.md` | `6ab742e9278205d5aa84ec2eecaaa716adcd4b35fc11f849af517ab0f7dcaedf` |

`SHA256SUMS` also includes `OUT/MANIFEST.md`; the manifest deliberately does
not include its own hash.

## Frozen mathematical ledgers

```text
rectangle pairs                         2,550
primitive ratios                        2,500
candidate budget                       16,384
activation cells                   41,779,200
nonnegative / equality       36,657,425 / 0
activation-sign byte FNV     9d3e995e23a7695a
common deck count / FNV      5,257 / 60f329212844f8ac
all bodies / failures        14,307,150 / 0
body checks / max prefix     508,103,822 / 4,129
half-budget failures                     7
weakest pair                         462,626
weakest repair mask                0x6042229
weakest mass          16477591782853/259521949879920
weakest gap                    7663493/259521949879920
post-THM-4276 residual                174,741
rectangle-new edges                    2,419
final residual                        172,322
final top layer              (256,670),(384,670)
```

The activity-independent covering-design argument additionally proves that
any disjoint-covering deck has at least 52 distinct eight-subset repair masks;
it does not claim that 52 is attainable.

## Independent routes and controls

- The primary engine fully sorts the repair universe and uses primitive
  endpoint-prefix integration with `(u,v,g)` retained.
- The literal engine selects candidates with a bounded heap, builds direct
  `lcm(D,14q,14r)` walls, intersects literal interval lists, and recursively
  enumerates bodies. The fixed-pool geometry and basic mask/FNV utilities are
  intentionally shared; pair integration, candidate selection, and body
  enumeration are detached.
- Both engines use the repaired exact unsigned-256 normalized-gap comparator.
  A separate Python-bigint path checks 81 boundary products and 100,000
  deterministic random comparisons.
- O3, O3+`NDEBUG`, and O1+UBSan semantic transcripts are byte-identical within
  each engine. Both UBSan stderr files are empty.
- Normal and `-O` Python runs are byte-identical for the transcript verifier,
  comparator oracle, and proof-graph postprocessor.
- `-Wall -Wextra -Wshadow -fsyntax-only` emits no diagnostics for either C++
  source. Runtime correctness gates use exceptions and survive `NDEBUG`; no
  C/C++ `assert` is load-bearing.

Replay environment:

```text
Apple clang 17.0.0 (clang-1700.0.13.5)
Target arm64-apple-darwin24.3.0
Python 3.10.0
Darwin 24.3.0 arm64
```

Raw wall times were recorded on a shared, concurrently loaded machine and are
not benchmarks:

| route | O3 | O3 + NDEBUG | UBSan |
|---|---:|---:|---:|
| primary endpoint | 1727.049644 s | 2642.793511 s | 3672.050935 s |
| independent literal | 1483.123519 s | 2439.110956 s | 2970.845444 s |

## Correction lineage and excluded files

An exploratory signed-128 cross-product selector could overflow and produced
the invalid primary diagnostic weakest cell `(499,640)`. It did not feed
activation, deck intersection, or body closure. The final sources replace it
with exact unsigned-256 products and both engines now agree on `(462,626)`.

The canonical packet excludes all compiled binaries, raw timed transcripts,
warning scratch files, optimized duplicate Python transcripts, `__pycache__`,
and the local `obsolete_exploratory_overflow_unsafe/` quarantine. The separate
post-freeze endpoint-670 bridge scouts are exploratory and are not THM-4277
dependencies or artifacts.
