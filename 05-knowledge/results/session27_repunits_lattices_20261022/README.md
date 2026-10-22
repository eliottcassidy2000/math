# Session27: repunits, lattice observables, and reconstruction from one bit

Start with the [final synthesis](SYNTHESIS.md), [exploration record](EXPLORATION.md), and [final source-bound audit](session27-final-audit/AUDIT.md).

This independently audited research packet recovers prior repunit work, distinguishes labelled from unlabelled tournament counts, and states precise modular/lattice bridges and obstructions. No theorem IDs are reserved. The final checkpoint applies six narrowly reviewed historical repairs, recorded as MISTAKE-658; all other concurrent canon is preserved.

- [Repunit recovery and measure-aware rarity](session27-repunits/REPORT.md): provenance routes, full local-search catalogues, and the automorphism-weighted SCC counting law.
- [Typed modular and lattice synthesis](session27-modular/report.md): theta/eta normalizations, E8 and rank-24 characters, the rational-versus-positive-weight modular obstruction, and carefully separated 1729 mechanisms.
- [Independent audit](session27-repunit-modular-audit/REVIEW.md) and [source-bound manifest](session27-repunit-modular-audit/manifest.json).

The complete plain-text JSON catalogues are included: the initial and expanded passes record their own commits, available worktree files, omissions, line previews and hashes. They are discovery indexes, not exhaustive recovery of every historical Git blob or every use of 31, 63 and 127. Initial fetched sources retain their connector provenance. The bounded remote-search addendum and its raw results are now included, with explicit caps and missing-coverage limits. The E8-holes/1729 addendum is included in the second checkpoint with its separate, bounded recurrence audit; it was not part of the initial audit.

## Reproduction

From session27-repunits run `python3 verify_weighted.py`; from session27-modular run `python3 check.py` and `python3 check_camion.py`. Saved JSON outputs and independent replays are included. To repeat discovery against a checkout, run `python3 discover.py /path/to/math /tmp/catalogue.json`; checkout availability and later commits can change the result.

Run `python3 session27-repunit-modular-audit/independent_checks.py /path/to/math` from this packet directory. The checker needs the referenced Session24 order-eight CSV and catalogue source files in that checkout. It independently enumerates small tournaments, checks the two order-eight H=63 witnesses, root coordinates, oscillator coefficients and the Camion control. Catalogue hash differences from subsequent repository changes are reported separately. Without an argument it locates the containing Git repository.

Only the independent checker's two workspace paths were adapted for distribution: repository root is an explicit argument or Git-discovered, and catalogue is a sibling file. Mathematical code is unchanged. The original audit manifest retains the pre-portability script hash and original absolute source locations; these are historical provenance, not distributed-path requirements. The packaged replay and publication plan record the adapted script bytes. `bind_manifest.py` is the original audit-time provenance utility, not a portable reproduction command.

## Correction proposals and boundaries

The initial `session27-repunits/corrections.patch` and `session27-modular/historical_1729_correction_proposal.md` remain historical proposals, not the applied patch. The final authoritative six-file repair is `session27-corrections/corrections.patch`, with exact preimage/proposed hashes and rationale in that directory, applied at the original repository paths and logged as MISTAKE-658. Its scope differs from the initial proposals; do not apply the older proposal mechanically. Read the audit's scope and limits before relying on a claim. Imported identities remain cited results; finite checks corroborate but do not prove all-parameter statements. No result here solves LRC(14) or supplies a universal repunit-rarity law.


## Second checkpoint: manuscript controls and new mechanisms

- [Four-manuscript review](session27-papers-review/REVIEW.md), [independent lattice-controls audit](session27-lattice-controls-audit/AUDIT.md), and [independent Niemeier restriction audit](session27-niemeier-audit/independent-restriction-audit.md).
- [Rank-weighted orthogonal avoidance](session27-orthogonal/ORTHOGONAL_AVOIDANCE.md) and [independent audit](session27-orthogonal-audit/AUDIT.md): exact incidence and avoidance formulas, with a counterexample showing that pairwise compatibility alone loses span rank.
- [Reciprocal-power recurrence theorem](session27-recurrence/RECIPROCAL_POWER_THEOREM.md) and [independent audit](session27-recurrence-audit/AUDIT.md): the universal order-2d+1 annihilator, generic minimality, and the exceptional repeated-root boundary.
- [Recovered 1729 mechanisms and E8 holes](session27-modular/PRIOR_1729_AND_E8_HOLES.md): bounded historical recovery and typed null-quotient/hole statements. The recurrence audit checks the specified claims, not the whole inherited directory or complete 27-dimensional historical script.

The supplied attachment titles were:
1. An Extremal Even Unimodular Lattice in Dimension 96.
2. Dual-Horizon Spectral Rigidity in Dimension 96: Deterministic Cyclotomic-Symplectic Architecture, Rational Certification, and the Extremal Frontier.
3. On the Restriction of Genus-5 Siegel Theta Differences of Niemeier Lattices to the Jacobian and Hyperelliptic Loci.
4. Finite Algebraic Compression in ADE Root-System Moment Problems: An Exact Four-Moment Rigidity Threshold.

Only our reviews, computations, outputs and numeric matrix data are distributed. No supplied PDFs, full extracted manuscript texts, third-party paper copies, page images or compiled binaries are included. The old 96-dimensional matrix is inert numeric data obtained from https://srfp311t1.com/G96_reduced.txt, linked by the first manuscript; its exact SHA256 is 7a6a3b8f761a9f02e1a687be8040aaf403c7ac75c33d92f010f2dd131f17bc24. This matrix is not an extremality certificate or an identification with a later construction.

The ADE census is exact within the stated formal-multiset universe. Sparse signed-coordinate checks have their stated support bound only. Rounded-eigenvector and numerical Gram-LLL scripts are heuristics, not exhaustive minimum proofs; their saved integer witnesses are checked separately with Python exact arithmetic. None of these computations certifies a norm-10 minimum in dimension 96.

E8 source caveat: the cited Cambridge treatment contains an intermediate cross-polytope group ambiguity. The full abstract 16-vertex cross-polytope has symmetry W(B8); its E8 affine lattice-isometry stabilizer is the index-two subgroup W(D8). The recurrence audit supplies the independent justification. Do not conflate these groups.

### Second-checkpoint reproduction and provenance

In session27-papers-review, run `python3 controls.py`, `python3 root_counts.py`, `python3 collision_mechanism.py`, and `python3 check_matrix.py`; NumPy and SymPy are needed for the relevant controls. Compile `g++ -O3 audit_ade.cpp -o /tmp/session27_ade` and run `/tmp/session27_ade`. The heuristic scripts lll_control.py and short_vectors.py write their witness files beside themselves; run them on a copy if retaining saved provenance matters. sparse_hostile.py performs only its declared bounded search.

In session27-lattice-controls-audit, run `python3 independent_controls.py`. It writes PACKAGED_INPUT_SHA256SUMS.txt for distributed inputs only, leaving the historical input manifest intact. Compile `g++ -O3 sparse_independent.cpp -o /tmp/session27_sparse` and run `/tmp/session27_sparse ../session27-papers-review/downloaded_G96_reduced.txt`. In session27-niemeier-audit run `python3 count_a3a1.py`.

In session27-orthogonal run `python3 verify.py` and `python3 verify_graph_counterexample.py`; in session27-orthogonal-audit run `python3 independent_control.py`. In session27-recurrence run `python3 verify_reciprocal_powers.py`; in session27-recurrence-audit run `python3 independent_controls.py`. The inherited modular script is retained with its original replay and historical scope, not upgraded to a newly independent audit.

Seven scripts receive packaging-only path/output-manifest adaptations, recorded in SECOND_PORTABILITY.diff. Their mathematical operations are unchanged and reruns are recorded in SECOND_PORTABLE_REPLAYS.json. Original INPUT/OUTPUT checksum files bind the original workspace scripts and omitted source PDFs/texts; those original-source audit records do not assert that every referenced input is distributed or that adapted script hashes are unchanged. SECOND_SOURCE_PROVENANCE.json records source hashes before packaging. The original audit outputs remain historical evidence; portable reruns do not reproduce the omitted manuscript-reading portion of an audit.


## Final checkpoint: faithful marked responses and repaired history

The headline theorem reconstructs every arc of a labelled k-vertex tournament using attachment experiments with k-1 ordered base vertices. The second least significant bit ((H-1)/2) mod 2 of each Hamiltonian-path response suffices. Separately, the signed full-H reconstruction implies faithfulness for every raw modulus M>=3; raw H mod 2 is constant and fails. This distinguishes a transformed bit from dividing by 2 in an arbitrary residue ring. Two ports have exact tomography formulas and finite injectivity through k=6, not a proved all-order two-port theorem.

- [Integer reconstruction](session27-response/faithfulness-proof.md), [raw-modulus refinement](session27-response/MODULAR_FAITHFULNESS.md), and [odd-cycle cancellation proof](ODD_CYCLE_PARITY_CHARACTER.md), with the corresponding source-bound reviews in session27-response-audit.
- [Identical E8 root, different Hamiltonian counts](SCORE_ROOT_COLLISION.md), [Lorentzian score bridge and its domain boundary](LORENTZIAN_SCORE_BRIDGE.md), and [topology scope](TOPOLOGY_SCOPE.md). No numerical analogy is promoted to a target-preserving equivalence.
- [Bounded remote discovery](session27-repunits/REMOTE_SEARCH.md) adds five raw query records and three recovered full repository sources. Four queries hit their 100-result caps; this does not claim all remote history was recovered.
- [Final historical repairs](session27-corrections/RATIONALE.md): six source-bound corrections preserve valid numerical and proof survivors. MISTAKE-658 records the failed implications and missing coordinates. Original/proposed duplicate trees are omitted; the patch, source hashes and pinned preimages retain reconstruction provenance.

### Final reproduction and packaging

The response directory includes standard-library Python controls and census.cpp. From that directory run audit.py, structure.py, two_port.py, faithfulness.py, modular_faithfulness.py and odd_cycle_faithfulness.py; the saved normal and optimized-mode outputs are included. Compile census.cpp with a C++ compiler and follow its documented scope. The unfinished order-seven census output and all binaries/caches are omitted. The finite complete census ends at core order six.

From session27-response-audit run independent.py, modular_independent.py, reconstruction_independent.py and root_addendum_independent.py. The two historical workspace paths in independent.py and root_addendum_independent.py are adapted to packet-relative inputs; FINAL_PORTABILITY.diff records those changes and FINAL_PORTABLE_REPLAYS.json records byte-identical reruns. These adaptations do not alter mathematics. Original audit hash manifests retain the original script hashes and original source locations.

At packet root, run score_root_hostile.py, score_first_failure.py and check_lorentzian_score.py for the explicit score examples and Lorentzian identities. The added Lorentzian independent audit and checker are in session27-recurrence-audit. The final audit may refer to its original workspace and full-repository preimages: its source-bound findings remain distinct from portable computational replay, and omitted manuscript PDFs are not silently supplied by these checks.

The initial catalogue/manifest source hashes predate the deliberate six-file repairs. A later catalogue replay can therefore report those changed paths while retaining its own original snapshot and exact line previews. Historical snapshots and the newly corrected source surface are intentionally both preserved.

Run `python3 session27-final-audit/check_portable.py /path/to/math` from the packet directory to repeat the final 74 integration checks. This wrapper reconstructs the omitted correction preimages from their pinned Git blob IDs and derives proposed snapshots by applying the published patch in a temporary directory. It does not edit the checkout. A checkout with those original blobs is required; the original source-bound checker is unchanged. The portable replay is byte-identical to check_results.json.

## Late primary-source check: Smale and repunit factors

[Smale, binary repunits, and modular forms](session27-smale/ADDENDUM.md), with [source/arithmetic check](session27-smale/SOURCE_CHECK.md), closes the named topology thread. This CITED classical addendum distinguishes the Kervaire–Milnor subgroup order from the full homotopy-sphere group, verifies the factors 7, 31 and 127 in three exact examples, and records the Bernoulli normalization without inferring tournament rarity. It was checked after the frozen final integration audit and is outside that earlier audit scope. No earlier audited synthesis or theorem was changed.
