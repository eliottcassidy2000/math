# Independent audit of all 715 nine-core/four-tail families

**Status: PASS. Complete independent full-comb reconstruction and all715
compact tree bindings passed in normal and optimized Python. The resulting
receipts are byte-identical, as are the control receipts.**

## Mathematical scope checked

For each nine-element set C contained in {1,...,13}, every four positive
integer speeds, including repeated speeds, admit a common time at radius
1/14. This is a bounded-core theorem, not unrestricted LRC(14).

The recursive proof uses two independent ingredients:

1. The partial-period conductor in the producer's P1 has constants
   c1=1/7, c2=1/3, c3=3/5, c4=1. For k>=2, c_k<=1 permits an interval
   of length c_k/u. Its smallest-speed bad measure is at most 2r/u;
   every other speed's bad measure is strictly less than
   2r*c_k/u+2r/u. The sum is strictly below c_k/u. The single-tail
   boundary is inclusive but may retain only isolated closed endpoints.
2. A branch is inherited only when a tail is an element of {1,...,13}
   absent from the original nine-core. This gives ten DISTINCT bounded
   speeds. The other three attachments include any already fixed prefix
   and all future tails; Session121 T3 applies without extra hypotheses.
   Neither a repeated core speed nor a new speed greater than13 is
   accepted by this sealed certificate's specific inheritance rule.
   A separate padding argument can strengthen the repeated-speed branch;
   it is not used by this four-tail certificate.

Every recursive internal cover must have a strictly positive largest
interval width. Empty or point-only internal covers cannot be silently
accepted. Every integer v in the sorted finite range lo<=v<cutoff is
accounted for. All larger sorted tails are covered analytically. Thus this
is a complete proof reduction to the finite tree, not a bounded sample.

## Independent implementation

`geometry.py` builds every speed's entire closed comb on [0,1] and uses a
two-pointer interval intersection. It does not import producer code or
its copied helper. In contrast, the producer clips tooth indices to the
current intervals and computes last-tail arrivals by phase. The auditor
constructs the whole intersection and takes its first left endpoint.
Every finite leaf witness is tested against all original speeds.

`check.py` independently reconstructs every complete recursive tree and
compares all per-core fields, canonical tree bytes and SHA256 commitment.
The receipt supplies no search cutoffs or widths to the auditor. Its
715 core lists must exactly equal the lexicographic combination universe.
Missing, repeated, reordered or substituted roots are rejected.

`controls.py` tests partial-period conductors at their radius boundaries,
closed endpoint equality, an undersized interval hostile, individually
successful but jointly failing attachments, and invalid inheritance
conditions. It mutates core coverage, coefficients, geometry, cutoffs,
counters, complete-tree branch coverage and inherited markers.

## Reproduction

From this directory, using the published producer receipt:

    python3 -B check.py ../lrc/four_tail_receipt.json > result.json
    python3 -B -O check.py ../lrc/four_tail_receipt.json > result.optimized.json
    cmp result.json result.optimized.json
    python3 -B controls.py ../lrc/four_tail_receipt.json > controls.json
    python3 -B -O controls.py ../lrc/four_tail_receipt.json > controls.optimized.json
    cmp controls.json controls.optimized.json

The preliminary independent traversal `replay.py` additionally records all
per-root original safe endpoints and per-level counts without reading the
producer receipt. It is redundant evidence, not a hidden input to check.py.

## Completed preliminary independent totals

- 715 exact nine-core roots.
- 23,281 three-tail width nodes.
- 218,081 two-tail width nodes.
- 407,207 one-tail width nodes.
- 44,133 literal finite leaves, each with original-speed witness.
- 52,347 correctly scoped inherited branches.
- Maximum cutoffs by remaining tail count1,2,3,4:72,94,104,92.

The preceding Session121 ten-core theorem is an explicit proved dependency;
its full 286-core certificate is not silently reclassified as a new result.

## Final binding and hostile controls

Both independent normal/optimized runs reproduced every canonical full tree,
including all 44,133 finite witnesses and52,347 inherited branches.
All715 compact commitments agree with independently reconstructed geometry.
The controls pass4,428 conductor boundary cases, seven hostile/endpoint checks,
and reject all eight certificate mutations. Exact results and receipt SHA256
are recorded in `result.json` and `controls.json`, with identical `.optimized`
companions. No external prover, Lean check, or human review is claimed.

## Lossless receipt transport

The full compact JSON is distributed as a gzip/base64 text archive to keep
repository API transfers small. From this directory, first run
`python3 -B ../materialize_receipts.py` to recover the exact named JSON files.
The decoder checks byte counts and SHA256. The proof still requires full
regeneration and independent replay; compression and hashes are not proofs.
The packet-level `verify.py` decodes into memory automatically and compares
the regenerated receipts byte for byte.
