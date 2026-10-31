# Independent bandwidth audit

**Pass.** Exact transfer and all-order width-two obstruction are sound within the stated presentation model. The finite width-three certificate has exactly the claimed coverage.

The carrier preserves endpoint pairing. An interior frontier vertex cannot accept any new selected edge; endpoint-only access is therefore complete. Distinct component occurrences contribute separate multiplicities. A forgotten end can always point to a future vertex, whereas no future vertex can point to a forgotten start. Discarding two forgotten starts is consequently permanent and safe. No isolated directed cycle is retained. These facts justify the state aggregation and the bijection with Hamiltonian paths.

The independent checker uses two different representations beyond the original transfer:

* It materializes abstract endpoint components as actual short directed paths, inserts a new vertex and candidate incident edges, reconstructs components, detects cycles by unreached vertices, and projects endpoints. It independently closes the carriers for widths0–4 and compares every transition to the producer.
* It exhausts actual edge subsets on small physical tournaments, filters indegree/outdegree violations and cycles, reconstructs full path components, then projects their endpoints. All281 tested full state distributions agree, not merely their accepted counts.

It also reconstructs physical adjacency independently, exhausts all10,923 width-two presentations through order8, and verifies all825 width-three positive witnesses. The47 witness is checked by all5,040 permutations of seven vertices. Hostile controls reject an altered count, invalid cut mask, and a single-edge mutation.

The all-order step is justified separately: SCC path counts multiply; every nontrivial strong component has count≥3; primality leaves exactly one count47 strong factor. The primary Busch paper was read during audit: Theorem1, not an inferred empirical bound, gives the sufficient uniform estimate. Induced width does not increase under compression of indices. Thus the n≤8 census exhausts every possible strong factor. The argument does not assert that the original tournament itself has at most eight vertices.

The positive witness file has one entry per target, no duplicate targets, valid masks, and maximum order9. It proves only the stated finite prefix. The scripts and receipts make no inference from an unused discovery search's completeness or statistics.

All published original source files are copied byte-for-byte. The curated README omits nonessential exploratory claims rather than silently treating them as verified dependencies. All verification checks are optimization-safe. The run-check harness verifies frozen output identities and normal/optimized equality. A copied package is replayed from a working directory outside the package before publication preparation.
