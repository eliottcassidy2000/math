# Corrected, independently replayed endpoint-four finite head

2026-10-25. **FINITE-EXACT:** this packet revalidates the complete closed
endpoint-four marked-source universes at orders 8, 9 and 10. Both the
ordinary and optimized runs passed every source. It does not
repeat the legacy all-layer census or use that wider scope as a dependency.

## Inheritance and the concrete decoder repair

The inherited mechanism is THM-4732 section 3, at
`01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md`.
Its original k=4 receipts are
`04-computation/experiments/blue_line_black_detours_20261020_n8.out`,
`_n9.out`, and `_n10.out`: respectively 640, 10240 and 143360 sources,
all repaired in three black moves. The source iterator includes every
closed reflection-fixed standard-path presentation in that endpoint layer.
It does not filter by clones, strongness, mixedness, isomorphism class, or
target type. In particular, its pure-blue test counts failures rather than
silently removing them from the source denominator.

The inspected original source was
`04-computation/experiments/blue_line_black_detours_20261020.cpp`, with
SHA256 values before the correction:

```
raw: a7ce0404c5e4146f9f5dd0fbd94945c6ee1f1bdbc98d02445fe85d4ea875e2fd
LF:  97c9b16d3a73aea2146ae13c8667e65198f5683c003f2546a0951ef6627638c4
```

Original line 61 calculated `nxt` before checking whether a degree-one
endpoint was reached prematurely. On a disconnected path-plus-cycle
candidate, that calculation can read uninitialized `nb[cur][1]`.
The guard then rejects the candidate, but C++ still makes the earlier read
undefined behavior. This is a concrete decoder defect, not evidence that
any recorded detour or count was false.

The [corrected copy](blue_line_black_detours_corrected.cpp) moves only that
guard before the neighbor read. It additionally inserts explicit decoder
controls and emits every successful literal certificate. Its LF SHA256 is

```
4499981f6b76688d9477e40bb1aed3312a4bd226495cbab729d84d31212ff232
```

[legacy_source.patch](legacy_source.patch) records the complete difference
from the original LF-normalized source. The original file was not edited
by this lane; its hashes and Git history preserve the pre-correction
version if the coordinator applies the same guard repair there.

After these runs the coordinator applied precisely the same guard-only
repair to the canonical source. An independent byte comparison against
`git show 2a9a9c342e5940f0e50b97fba87a7cfb7a24de7b:<original-path>`
confirmed that this is its only LF-normalized change. The repaired
canonical LF SHA256 is
`304d0be86bf48bb91361f7930e68c0083e37459cc8487f5f557b47905c682b2e`.
The local corrected copy remains separately instrumented and hash-pinned.

The minimal disconnected candidate surviving the degree and edge-count
prefilters has five vertices and edges
`{01,23,34,24}`: an isolated directed edge and a triangle. Starting at 0
in the transitive tournament, the decoder reaches endpoint 1 at step 2.
It must reject before reading the nonexistent second neighbor. The
[minimal control receipt](decoder_control.json) executes this case,
alongside accepted forward and reversed directed paths. Each of the three
main censuses also runs the corresponding edge-plus-cycle control at its
actual order. Five is minimal because the path component needs at least
two vertices and a simple cycle needs at least three.

## Independent coverage and literal checks

[check.py](check.py) uses no imported tournament engine. It streams each
certificate from the corrected executable and checks:

1. The source rows form a tournament, P=(0,...,n-1) is directed and
   closed, P's reversal is an anti-automorphism, and endpoint scores
   are exactly (n-5,4).
2. Each of the three words is a permutation, directed at its actual
   successive source, and fails a full anti-reversal test there.
3. Every move is recomputed pair by pair from its definition; the final
   complete row tuple equals the independently computed K_P(T).
4. The three unordered path-edge sets XOR to E(P).
5. No marked source appears twice. The number of unique validated
   sources equals an independently counted complete universe.

For that last test, write the actual incoming row of vertex 0 as
`{n-1} union I`. Exactly three members of `{2,...,n-2}` form I.
The full reversal identity determines all pairs meeting the last vertex
from this first row. The remaining freedom consists precisely of the
reflection orbits of nonpath pairs entirely in the middle vertices
`{1,...,n-2}`. The checker explicitly constructs these pair-orbits:

| n | Port sets | Independent middle bits | Complete marked sources |
|---|---:|---:|---:|
| 8 | 10 | 6 | 640 |
| 9 | 20 | 9 | 10240 |
| 10 | 35 | 12 | 143360 |

Thus the total is `binomial(n-3,3) * 2^floor((n-3)^2/4)`.
Every admissible port set and middle bit choice defines one source, and
different choices give different literal source rows. Validity, uniqueness
and this exact cardinality prove coverage independently of the producer's
iterator and Hamilton-path search. No completeness assumption about the
search algorithm is needed after a certificate for every source has passed.

For any marked word P on arbitrary vertex names, relabel in its path
order, use the standard-word certificate, and transport all three words
back. This preserves all arcs, ports, clone mergers, and the literal
target; it proves the all-marked-word statement rather than a statement
only about unmarked isomorphism classes.

The canonical hostile is the order-seven source
`(6,4,8,19,103,79,15)`, with exact target
`(122,124,121,112,32,64,16)`. The independent checker enumerates the
target's nine Hamiltonian words and proves all are blue. It also rejects
a proposed repair that uses the source's blue base word. The corrected
near miss is the early-endpoint read above. The least-used relevant
sidecar is the complete labelled certificate stream: its SHA256 is
recorded rather than replacing literal equality by an isomorphism test.

## Reproduction and receipts

From the repository root:

```
python 05-knowledge/results/lca_ports_reverse_20261025/finite_head/check.py --census --output 05-knowledge/results/lca_ports_reverse_20261025/finite_head/normal.json
python -O 05-knowledge/results/lca_ports_reverse_20261025/finite_head/check.py --census --output 05-knowledge/results/lca_ports_reverse_20261025/finite_head/optimized.json
python 05-knowledge/results/lca_ports_reverse_20261025/finite_head/check.py --decoder-only --output 05-knowledge/results/lca_ports_reverse_20261025/finite_head/decoder_control.json
```

The runner invokes GCC 15.2.0 with
`g++ -O2 -std=c++17 -Wall -Wextra corrected.cpp -o <TEMP>/finite-head.exe`,
then invokes the executable as `8 4`, `9 4`, and `10 4`: no subcase,
range, random-sample, or stride filter. Temporary executables and stderr
captures stay outside the repository and are removed after execution.
Put `g++` on PATH, or pass `--compiler` followed by its absolute executable
path. Compiler discovery has no dependency on this machine's directory
layout. Without `--census` or `--decoder-only`, `check.py` performs only
the fast comparison of stored receipts and does not require a compiler.
The decoder-only invocation is `5 4 0 1`; it exercises startup controls
without claiming a census at order five.

The sole compiler diagnostic is the inherited compact-format
`-Wmisleading-indentation` warning at line 56; the branch semantics there
are unambiguous. The exact compiler output is retained in each receipt.
The checker uses explicit exceptions, not Python `assert`, for all
validity gates. The certificate-stream hashes and mathematical results
agree between ordinary and optimized runs. Elapsed times and timing
fields in the producer's original reporting are incidental and may differ.

The raw streams are processed completely but not committed; they can be
regenerated by compiling the preserved corrected source. The receipts
retain counts, unique-source cardinalities, stream hashes, original
non-certificate output and progress output. This is a fresh exhaustive
check of the specified three finite universes, not a new general search.

The [normal receipt](normal.json) and [optimized receipt](optimized.json)
give these identical LF-normalized certificate-stream SHA256 values:

```
n=8:  904be5eefec31502cda7f0d3a84fc08a7aacff758be48b19c7ca5ac93a6c61c4
n=9:  f3aaeb580fe0ba8237d709389f0af0d691fa5c94af6a5f56adc518d5c67dd706
n=10: 645f6553e90e5d00a06f9789a678529869612e16a157fbb447653e1a4850aece
```

Normal elapsed seconds were 0.211, 2.147 and 142.491; optimized elapsed
seconds were 0.237, 2.426 and 141.922. These timings describe this replay
only and are not a general complexity comparison.

[compare_receipts.py](compare_receipts.py) extracts a deterministic
[paired semantic summary](paired_semantics.json). It checks the corrected
source hash, PASS statuses, exact complete-universe counts, unfiltered
arguments, decoder-control presence, and matching certificate-stream
digests. It retains compiler version and flags, but strips elapsed times,
progress timing, absolute paths in diagnostics, and incidental timed
reporting. The raw observational receipts remain unchanged. The quick
normal/optimized commands

```
python 05-knowledge/results/lca_ports_reverse_20261025/finite_head/check.py
python -O 05-knowledge/results/lca_ports_reverse_20261025/finite_head/check.py
```

produce byte-identical summaries. They explicitly identify themselves as
comparisons of completed receipts; the `--census` commands above are the
full independent replay. The quick-mode switch was added after the two
exhaustive runs; it does not alter their source, decoder, certificate
validation, or source-universe enumeration.

## Scope and portfolio

The anchor is the complete finite head needed by the assembled closed
endpoint-four theorem. The niche is decoder correctness on disconnected
degree-valid candidates. The wildcard is transport from marked labelled
sources to arbitrary marked words without losing literal clone data.
The live concept board is: source completeness; actual-source blackness;
literal target; disconnected decoder hostile; small-order pure-blue target.
None of these can be replaced by path counts or quotient isomorphism.

This packet does not certify all other layers of the legacy census,
minimal detour lengths beyond the immediate one/two-move obstruction,
open endpoint-four paths, or all-order black-core connectivity.
