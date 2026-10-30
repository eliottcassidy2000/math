# Sources, provenance and claim boundaries

**CITED primary background / pinned repository inheritance.**
The startup routes, current guardrails, RESEARCH-PROTOCOL and CORE-PAPERS
were read before making literature claims. No historical reflection or
session log is used as the current truth source.

## Primary literature checked in this session

1. Joshua Brakensiek, Venkatesan Guruswami, Marcin Wrochna and Stanislav
   Živný, *The Power of the Combined Basic LP and Affine Relaxation for
   Promise CSPs*, SIAM Journal on Computing49(6),2020,1232–1248.
   [Version3 abstract and metadata](https://arxiv.org/abs/1907.04383v3);
   [publisher DOI](https://doi.org/10.1137/20M1312745).
   The checked abstract states the polynomial algorithm for infinitely
   many symmetric polymorphisms and its block-symmetric extension.
   This establishes prior context for P3–P4. The packet does not reprove
   or claim the full algorithmic characterization in that paper, and no
   claim of full-paper review is made. The arguments here stand on their
   displayed proofs; ordinary BLP alone is insufficient in our triangle
   example.

2. John Irving and Mohamed Omar, *Revisiting The Rédei-Berge Symmetric
   Functions via Matrix Algebra*, arXiv2412.10572v3,
   [section3.4, Corollaries19–20](https://arxiv.org/html/2412.10572v3).
   The primary formula credits Grinberg–Stanley and counts odd cycle
   collections with weight two per nontrivial cycle. Our use concerns
   Hamiltonian counts and their parity; it gives neither tournament
   isomorphism classification nor a new enumeration algorithm.
   The complement includes loops in the paper's convention, whereas a
   Hamiltonian path uses distinct vertices. For tournaments, reversal
   identifies its Hamiltonian paths with those of the loop-free converse.

These are background attributions, not claims that P1–P6 are new to the
literature. The contribution of this session is the scoped synthesis,
the exact cover-denominator formulation, its proof, the reused native
hostiles, and explicit obligations for the repo's target problems.

3. Pascal O. Vontobel and Ralf Koetter, *Graph-Cover Decoding and Finite-Length
   Analysis of Message-Passing Iterative Decoding of LDPC Codes*,2005,
   [arXiv:cs/0512078v1](https://arxiv.org/abs/cs/0512078v1).
   The checked primary abstract connects graph-cover decoding to LP decoding
   through the fundamental polytope. This is prior context for the
   cover/marginal correspondence; P4b gives its own proof for the stated
   finite constraint model. No claim of full-paper review is made.

## Repository provenance

The principal inherited sources are recorded with exact Git revisions
and SHA256 blob digests in manifest.json. The starting mathematical
surface is cd6d98ba58a2418fb81054659bc37962676c8dd0. The later Session118
README/frontier is pinned separately at its incoming revision.
The late Session119 README, signed-graph classification and weighted-cycle
classification are pinned at4b423e24b8. They were read during push integration;
P2's concurrent overlap is credited and the ordered-carry hostile is replayed.

Packet-local labels P1–P6 are not reservations of global theorem IDs.
Whenever a global theorem is used, the slug/path accompanies its ID.

| Local result or boundary | Prior source and role |
|---|---|
| P1, P6 | odd_prime_function_cycle_20261018d.md: existing odd-shade cover theorem, anti-invariant count, odd-degree square-class descent |
| P2 | attachment_lattices_square_sum_20261023/square_sum/SHARED_CLOCKS.md: existing tree root/pin arithmetic and collision boundaries |
| Bipartite sign gauge | aliquot_66_five_fifteen_20261028/PROOFS.md, P6 |
| Odd-function dictionary | LEM-004-tournaments-are-odd-functions.md |
| Odd-cycle formula | THM-002-ocf.md, inherited through the reviewed1–2–4 packet and the checked primary source |
| Corrected digit tower | THM-466-the-2-adic-digits-of-H-are-odd-cycle-collection-counts.md, with the October26 correction; liouville_mobius_11_13_20261026/RECOVERY.md |
| Named-port reconstruction | session27_repunits_lattices_20261022/session27-response/ODD_CYCLE_FAITHFULNESS.md |
| Corrected action mismatch | THM-1380-free-Z2-but-no-odd-map.md |
| Fixed residue fibres | collatz_mod6_20260917_floor_sums_odd_functions.md |
| Reflection parity and Lucas | odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md; lucas_mirror_gap_spaces_20261028/SYNTHESIS.md |
|81/27 native executable | regular_parent_symmetry_20261029/verify.py::hostile_checks |
| Current Collatz/LRC frontier | terminal_fees_phase_covers_20261030/README.md and ADAPTIVE.md; their large certificates are cited, not rerun here |
| Incoming circuit frontier | session118_shared_gate_deletions_20261030/README.md and frontier/FRONTIER.md; conditional external-application scope retained |
| Odd-modulus counting distinction | rectangular_counting_compiler_20261030/README.md, prior scaled-integer reconstruction |
| Affine XOR instead of median | scaled_count_fft_carriers_20261030/THEORY.md,S1–S2: the actual Fano-line Hamming syndrome carrier |
| Concurrent signed-graph overlap | session119_odd_transport_classifications_20261030/cycles/SIGNED_GRAPH_INHOMOGENEOUS.md: complete inhomogeneous integral and modular classification |
| Odd-prime/ordered-carry boundary | session119_odd_transport_classifications_20261030/arithmetic/WEIGHTED_CYCLE_OBSTRUCTIONS.md: cyclic obstruction and the two D=-139 clock words |

The native checker reruns four recovered odd-cycle examples and the81/27
helper. It does not rerun the older million-tournament census, all contextual
XOR queries, the large E11 certificate or all28 LRC phase-cover certificates.
No Sha, Collatz, LRC, Fermat-prime or circuit lower-bound theorem is cited
as solved.
