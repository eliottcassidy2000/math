# Verification scope and reproduction

All mathematical claims in this packet are scoped in their proof files. The exact finite computations are controls and explicit witness certificates, not proofs of universal odd-value coverage. Chan–Pak remains open.

## Canonical replay command

From this packet's directory, run:

    python3 -I -B verify_all.py > replay.json
    cmp replay.json REPLAY.json
    sha256sum -c SHA256SUMS

The aggregate checker executes every included producer and independent checker in both normal and optimized Python, then requires byte-for-byte equality with the committed deterministic receipt. A failed subprocess or mismatch raises an explicit exception. It uses only Python's standard library and does not require the rest of the repository, network access, or scratch probes.

The individual proof and audit reports specify their own finite universes and checks. Some earlier audit reports retain the original research-directory names in their execution history. The aggregate command above is the canonical portable replay for the assembled layout.

## Independent computational paths

- Graft images: all 5,405 rooted labelled presentations through order five, literal grafts, independently counted cyclic families and lattice probes. Producer graph recognition separately compares all 33,792 rooted-at-zero tournaments of orders five and six with a bank constructed by literal grafting and relabelling. The independent audit does not import producer code.
- Multi-root constructions: exact permutation-slot coefficients versus literal subset Hamilton DP, sharp degree and SCC controls, physical Vandermonde rank witnesses, and full-state collision graphs. The independent audit performs 170,515 explicit checks, including 12,078 literal grafts and ten collision graphs through order fifteen.
- Boundary escape: producer root-zero profiles through order five and every pair of the 75 labelled root-zero source/gadget presentations through order four. The independent path uses all 5,405 rooted labelled presentations, 7,878 pairwise defect bounds, and 5,625 graft pairs. Singleton and sink cases are included rather than filtered out.
- Finite-menu boundary: separately generated gadget statistics, literal grafts, mixed-menu word controls, exact chamber congruences, stationary sink-root words, large-U seed tails, and the no-common-residue-gap hostile. Each receipt states its precise finite menu and depth.
- Adaptive arithmetic: exact operator identities, literal physical switch witnesses, and the common-context sensitivity law. Short searches, if reported, are finite experiments only; they do not establish a universal absence of radix constructions.

The optimization runs retain all load-bearing checks. Individual manifests and audit input hashes preserve source bindings; the root manifest additionally binds the entire final deliverable. Status-only changes made after an audit are recorded separately where applicable.

## Repository checks

The [final close-out record](SCOPED_CHECKS.json) reports the exact commands and outcomes for whitespace, local links, source syntax, integrity, and the repository's bounded documentation checker. The latter requires five unmodified documentation files excluded by the sparse checkout; they are materialized from the same verified Git revision before running it. This is not a claim that every unrelated repository test or external CI job ran.

All publication checkpoints use nonforce updates with an expected parent. The exact changed tree is compared with the local committed tree, and the resulting remote head is fetched before completion is reported. Other checkout changes are not staged or modified.
