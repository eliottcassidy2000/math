# 3-adic past guards close a long-match cycle alternative

19 October 2026. **PROVED, scoped; independent review PASS.**

Start with [IMAGE_GUARD_DESCENT.md](IMAGE_GUARD_DESCENT.md). It proves:

1. Exact bank image types v_3(y-1): K gives 3, V gives 4, H gives 1,
   and F increments the input type, with the x=1 infinite-type exception.
2. A two-sided guard criterion: equal long future parity blocks and
   incompatible past-image residues force strict first-crossing endpoint
   descent. An optional CRT formula strengthens the source-height gap.
3. An unbounded family of greedy interleaved words, indexed by all upper
   convergents of a fixed logarithmic rotation, has strict endpoint descent.
   Its episode counts are unbounded and its large words avoid the previous
   macroscopic pure-expansion repetition test.

The arithmetic exclusion is concrete: the matched checkpoints are V and VF
outputs, respectively 82 and 1 modulo 243. Their exact equality is impossible.

[PERIODIC_TAIL_MISMATCH_HOSTILE.md](PERIODIC_TAIL_MISMATCH_HOSTILE.md) proves
why purely formal refinements miss a subfamily. The huge exact control has
1,833,237,826,965,784,108,483,587,883,633,758,044,204,176 Terras bits. At its
one-percent matching scale, every matching separator contracts and every
future continuation agrees. The new past-image guard nevertheless proves
endpoint descent, using a compressed certificate rather than constructing
this word or an orbit.

## Dependencies and reproduction

The height argument depends on the separately audited
[matched-block extension](../extensions/MATCHED_BLOCK_EXTENSION.md), the
[greedy construction](../extensions/REPETITION_POOR_GREEDY_BOUNDARY.md), and
the [primary interleaved theorem](../INTERLEAVED_REPETITION_BARRIER.md).
The primary result's explicit Matveev specialization remains a cited input.

Run:

    python3 verify_periodic_tail_hostile.py
    python3 -O verify_periodic_tail_hostile.py

Both outputs match the supplied receipts: 30,113 explicit checks, including
512-term rational logarithm enclosures, exact bank image types, the modulo
243 separation, the CRT factor 81, two compressed periodic-tail examples,
and 39 interval-certified upper-convergent controls. Universal statements
are proved in the notes, not inferred from this finite control list.

General Collatz convergence,
arbitrary bank entry and the full realized S_d boundary remain open.

## Independent reproduction

Read the [independent review](review/IMAGE_GUARD_DESCENT_REVIEW.md). From this package directory, using Python 3 and its standard library:

    python3 -B review/independent_image_guard_controls.py
    python3 -B -O review/independent_image_guard_controls.py

Compare stdout with `review/independent-receipt.json` and `review/independent-optimized.json`. Both runs perform 238,861 explicit checks. The independent logarithm enclosure uses a different basis and truncation from the submitted verifier. All checks remain effective under optimization. The exact primary and hostile notes are retained byte-for-byte; `SOURCE_PROVENANCE.json` records the review-status and package-path presentation overlays.
