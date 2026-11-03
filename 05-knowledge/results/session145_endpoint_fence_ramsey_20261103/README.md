# Session 145: sharp attachment carriers and scoped bridges

**Two verified checkpoints.** The main new result is the smaller pooled-parity
attachment carrier and its exact six-vertex minimum. The geometry, Ramsey,
knot, and nine-label lanes add separate proved classifications and explicit
boundaries; their shared vocabulary is not used to claim a grand equivalence.

[Research synthesis](SYNTHESIS.md) · [Verification](VERIFICATION.md)

## Attachment carrier: the strongest advance

[Proof](carrier/POOLED_PARITY_CARRIER.md) ·
[Independent audit](carrier/independent/AUDIT.md) ·
[Frozen manifest](carrier/MANIFEST.json)

Pool named internal vertices across old path components, but retain isolated
labels, endpoint pairing, and each component's parity. This gives an all-order
closed rational carrier for means under all fixed source-oblivious attachments.
Its coordinate counts begin **1,1,2,7,28,126,571,2668**. At order six, 90 pairs of
the old 661 coordinates merge, and an independently reconstructed matrix has
rank 571 modulo 101 and 103. Hence the exact rational linear minimum is **571**,
or **570 varying coordinates** for normalized laws. The inherited minima
7,28,126 and original all-order upper bound remain correct.

All-order linear minimality remains open. The all-order classification is the
coarsest pairwise signed-monomial quotient, which is not itself a rank theorem.
Conditioning, source-adaptive attachments, and shared-source higher moments
require additional information.

## Three separate bridge lanes

- **[Fences and product packing](portfolio/fences/PROOF.md).** In the supported, nonoverlapping fence model with each field area at most one,
  split-graph pseudoforests obey total field area at most n/4. Endpoint-only product cacti
  have an exact balanced-cycle optimum. The known twelve-edge theta value
  `(89+26 sqrt(13))/16` beats the cactus value `6+3 sqrt(3)`.
  A [separate theorem](portfolio/fences/theta_addendum/PROOF.md) proves that
  theta value optimal when its two-edge path separates the fields and the
  terminal chord lies inside the outer polygon. The global optimum is open.
- **[Sparse Ramsey and knot interfaces](portfolio/ramsey_knots/PROOFS.md).**
  The maximum avoiding orders under the four stated forest/pseudoforest
  promises are exactly `2,3,4,5` times `(k-1)`. Matching closures have loop
  polynomial `product_(j=0)^(r-1)(z+2j)`, but matching does not determine knot
  type. Pre-existing closed components require an additional count m.
  The Paley-11 number 1729 is an endpoint-channel sum and an affine-orbit
  count; a specified directed arc occurs in 17290 Hamilton words.
- **[Nine labels and magic maps](portfolio/magic_frequency/PROOF.md).**
  The actual nine decimal labels cannot be fully magic under any permutation
  or nonconstant affine relabelling: their required center 555 is absent.
  In the prescribed Lo Shu placement, every real-valued magic map is exactly
  `c+pU+qV`, and strict increase is exactly `p>0, q>2p`. The raw layout's one
  defective diagonal is explained by its nonzero V-squared term.

The [independent bridge audit](portfolio/independent_bridge/AUDIT.md) checks
these scopes with separate exact algorithms. Exact unrestricted Ramsey
numbers, general twelve-edge product optimality, knot classification, and
higher-order attachment minimality are not claimed solved.

Both science packets have immutable manifests. This navigation overview and
the synthesis may change without rewriting either frozen packet. All earlier
research packets and concurrently published work are preserved unchanged.
