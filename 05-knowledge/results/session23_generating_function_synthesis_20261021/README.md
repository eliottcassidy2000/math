# Session23 generating-function synthesis

Research date: 2026-10-21. This packet contains the frozen, independently reviewed Session23 overview. Mathematical scope and publication provenance are separate: the enclosing verified repository checkpoint, not this wrapper, establishes publication status.

The [final synthesis](SESSION23_SYNTHESIS.md) has SHA-256 `d8815510d96f84aa15de8c52e25139a61de14de6a8c3ca37571039f3d097ec06`. Its [independent fidelity and scope review](session23-synthesis-review/REVIEW.md) passes these exact bytes. The [review receipt](session23-synthesis-review/receipt.json), [recorded controls](session23-synthesis-review/normal.json), and [review hashes](session23-synthesis-review/SHA256SUMS) retain the evidence boundary. The earlier draft and initial-binding snapshot are deliberately excluded.

## Six research packets, seven reviewed lanes

1. [Refined metagraph generating functions](../refined_metagraph_generating_functions_20261021/README.md): multivariate rationality, exact minimal eventual recurrences for every degree `3^r`, and the raw-automorphism-marker non-D-finiteness boundary. The first two recurrence orders were inherited; the all-r extension is the new scope.
2. [Fixed-presentation renewal limit laws](../renewal_tournament_limit_laws_20261021/README.md): fixed-h gap geometry, Gaussian laws, and an asymptotically uniform degree-`3^r` sampler. This is not unrestricted tournament sampling.
3. [Marked-path tournament series](../marked_path_tournament_series_20261021/README.md): rational normalized count channels and the narrow mixed-cycle centralizer correction. Correct tables and unrelated interpretations retain their scope.
4. [Fixed-symmetry-slack rationality](../fixed_symmetry_slack_generating_functions_20261021/README.md): the canonical finite compiler proves rationality at every fixed slack. The separately reviewed [slack-renewal transfer](../fixed_symmetry_slack_generating_functions_20261021/session23-root/SLACK_RENEWAL_TRANSFER.md) is bundled with this compiler packet, not a seventh standalone packet.
5. [Sharp symmetry-slack recurrences](../sharp_symmetry_slack_recurrences_20261021/README.md): exact ordinary, self-converse and converse-merged recurrence orders in every nonempty slack layer, with sufficient onset bounds. Numerators and initial coefficients are not supplied by a denominator alone.
6. [Split-word metagraph spectra](../split_word_metagraph_spectra_20261021/README.md): explicit slack-three spectra, two distinct Cartesian-square constructions, a Fibonacci-product bijection, and a literal wiring obstruction. Degree-support completeness is only within slack three; matching scalar invariants do not determine actual destinations.

## Evidence and replay boundary

All copied synthesis and allowlisted review files are unchanged. The review records 59 fidelity, binding, arithmetic and saved-output checks; its normal and optimized outputs are byte-identical. It compares seven current proof bindings with their independent reviews. It does not rerun all seven proof audits, constitute formal verification, or establish literature priority. Actual-tournament, formal-word and synthetic controls remain distinct universes.

The included `session23-synthesis-review/check_synthesis.py` is preserved as audit provenance. It expects the original session workspace layout, including the seven source/review lanes and the two recorded publication receipts. It is not a standalone checker runnable from this published overview directory. In the original workspace, the commands are:

    python3 -B session23-synthesis-review/check_synthesis.py
    python3 -B -O session23-synthesis-review/check_synthesis.py

Each research packet documents its own reproduction commands and dependencies. Source pins identify original source bytes; historical receipts describe their recorded observations rather than a fresh remote read in this packaging step.

## Release conditions and non-consequences

The six relative packet links and the bundled transfer must resolve to the reviewed release layout, with exact final proof and review bytes verified before this overview is released. These are release conditions, not an assertion that any unverified publication has happened or remains pending. Publication completion requires a verified remote commit and exact path readback; a local file, manifest, attempted push or platform approval request is not that evidence.

The optional HYP-9352 navigation patch is separate from this packet. Apply it only after checking its fresh preimage hash and the actual release layout. It preserves the existing scoped update, historical proposals and center frontier; it does not touch CURRENT-FRONTIER, the mistakes ledger, or any other canonical file.

Complete metagraph wiring, black connectivity, bridges, strong-atom/center coefficient inventories, and unrestricted tournament-enumeration complexity remain open. The 153N/all-horizon source and clock ledger is unchanged. No general Lonely Runner, Hadwiger–Nelson or Collatz theorem follows from these generating functions.

No Git operation, remote write or repository mutation was performed while preparing this packet.
