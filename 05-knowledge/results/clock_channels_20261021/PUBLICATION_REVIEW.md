# Publication review: finite static observers and exact integer clocks

**Independent scoped PASS.** The [producer proof](FINITE_OBSERVERS_INTEGER_CLOCKS.md),
all original receipts and manifests, and the [independent audit](audit/AUDIT.md)
are preserved byte-for-byte. The audit authenticates six producer manifest entries;
its independent normal and optimized receipts each retain 46,083 checks. The
producer normal and optimized receipts each retain 263 checks.

The source family is explicitly grounded. Its primitive braid, shortcut ROOT,
and ordinary ROOT clocks are different quantities; finite static observations
cannot reconstruct the specified full integer clock on the infinite family.
This is not a Collatz convergence proof, a density statement, or an impossibility
claim about finite alphabets with unbounded words, counters, or rule schemes.
The inherited finite frontier is not exhaustively enumerated by this audit.

The companion hash in the frozen provenance identifies the
[original idoneal companion](../idoneal_information_channels_20261021/provenance/original/IDONEAL_INFORMATION_CHANNELS.md),
which the earlier publication preserved as historical provenance. It does not
identify the separately corrected main companion. Only its literal 63/273 and
3/7 examples are used and independently replayed here; the grounded family and
its unbounded-clock proof do not depend on global idoneal completeness or the
companion's analytic background. No frozen provenance hash is silently updated.

## Portable reproduction

From this directory, run both of the following with normal Python and with -O:

    python3 -B verify_clock_channels.py
    python3 -B audit/independent_clock_channels.py --packet .

The first output must equal verification.json and verification.optimized.json.
The second must equal audit/independent.normal.json and
 audit/independent.optimized.json. Audit producer.replayed.normal.json and
producer.replayed.optimized.json also exactly match the producer output.
Authenticate the unchanged packet and audit with:

    sha256sum -c SHA256SUMS
    (cd audit && sha256sum -c SHA256SUMS)

The explicit --packet option is necessary with this publication layout; the
frozen auditor's original sibling-directory default is deliberately unchanged.
