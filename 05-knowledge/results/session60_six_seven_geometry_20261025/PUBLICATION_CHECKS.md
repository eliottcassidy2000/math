# First checkpoint verification

Status: all listed local checks pass; theorem scope is in the proofs and independent audits. Publication is established by the containing Git commit and independent Git readback, not this file alone.

The first checkpoint contains the root-threshold packet, deck/K3/Hodge/Gale packet, hourglass-symmetry/CM packet, the independent threshold audit and independent deck/Kummer/Hodge audit. Fano/Campedelli classification, generic curve sections and Richelot extensions are excluded until separately audited.

Six current primary producer/audit scripts were executed normally and with `python3 -O` in two separate frozen packet copies. All twelve executions passed; stdout agrees between modes, every script has zero bare assertion sites, and every generated file matches the frozen publication bytes. `REPLAY_RECEIPT.json` records source hashes and actual runs. Standard Python is sufficient except for the symmetry producer and independent deck auditor, which also use SymPy.

The deck checker originally contained 28 bare assertion sites. Its preserved original and initial `deck.normal.json` / `deck.optimized.json` in `audits/deck-kummer` are historical evidence, not current verification. The corrected script executes 91,064 always-active gates in both modes; the independent auditor executes 825. The intentional wrong-genus mutant fails in both modes, as recorded in `NEGATIVE_CONTROL_RECEIPT.json`. The mutant is a negative control, not a passing checker. Producer freeze receipts retain their original pending-audit timestamps; current gates are the independent audit verdicts.

`SOURCE_VERIFICATION.json` records original source bytes and corresponding publication bytes. Portable audit/certificate/checksum paths were normalized; the synthesis states its explicit first-checkpoint boundary. No theorem or expected mathematical value was changed by packaging. Captured `/tmp` paths in provenance are historical locations, not required reproduction paths. Native checksum lists are relative to their containing directories, including the two-level audit-source references.

Run `python3 verify_packet.py` from any working directory to check the complete file set, byte sizes and SHA-256 hashes. The manifest excludes itself; Python cache files are ignored, not published. Normal and optimized runs of this verifier pass. Native checksum closures, local Markdown routes, UTF-8 decoding and the repository's unchanged bounded-documentation gate are checked separately. No external CI result is implied by local replay.
