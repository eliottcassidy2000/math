# Tournament / automorphism-even graph work, merged with named responses

**Status:** INHERITED PROVED / CITED equinumerosity and specified bijections;
PROVED elementary Fourier specialization; FINITE-EXACT order-five controls.
This is a recovery and integration packet, not a claim to solve the classical
natural-bijection problem or the Hamilton spectrum.

The important connection is a change of basis for observables. A tournament
type is a point; an automorphism-even graph indexes a signed Fourier basis
function on those points. A bijection between the two index sets requires
additional choices. Prior work supplies several such choices, each preserving
different structure. Named two-port queries require the marked version of
this basis; forgetting the marks loses actual response information.

- [MERGE.md](MERGE.md): the recovered maps, their inverses, preserved data,
  failures, and the precise connection to cubic/quartic response coefficients.
- [SOURCE_LEDGER.md](SOURCE_LEDGER.md): authoritative prior files, correction
  lineage, the pasted sphere/MRB/projective proposal, and incoming frontier.
- [REVIEW.md](REVIEW.md): validity gate, independent checks and concept board.
- `receipt.json`: exact graph/tournament census through five and test receipt.
- `full_response.json.gz`: lossless complete order-five named response tables
  and their old-edge Fourier coefficients. No cubic or quartic is discarded.

At order five there are12 tournament types and12 sign-even graph types.
The latter have edge-parity split10+2. Converse has trace10-2=8, recovering
the eight self-converse tournament types. This is the inherited signed-trace
identity, not a volume or MRB calculation. The Euclidean unit-ball volume
does peak in integer dimension five, but the supplied common-decay claim
was already refuted in the November2 audit.

The new complete replay retains all1024 labelled five-tournaments,256
coefficients each,265,216 actual extension counts,102,400 independent cover
quartic checks, and twelve independent Hamilton-word full-polynomial checks.
It also checks51,200 forced zero quartic coefficients on odd old-edge Fourier
supports. Normal and optimized producers agree byte for byte.

Reproduce without changing the packet:

    python 05-knowledge/results/even_graph_response_merge_20261104/verify.py --replay

Incoming work now narrows a possible original-module XOR collision to highest
changed prime quotients of order at least nine, with a separate reviewed
proof. That improves the previous chat's nested boundary; unrestricted
full-response injectivity and completeness outside7/21 remain OPEN.
