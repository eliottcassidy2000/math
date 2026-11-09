# Reproduce Session233

Run from this packet directory:

    python3 -B replay.py

The standard-library-only driver runs eleven producer/audit programs both
normally and under Python -O. It compares stdout exactly and requires every
program's explicit deliberate-failure control to exit nonzero in both modes.
It writes REPLAY.json. The bridge auditor also rewrites its deterministic
receipt.json and values.json. No network, class database or compiled binary
is required.

Individual checks:

- graft/verify.py: inherited37/45 word certificates and original pair theorem.
- graft/verify_weighted.py: new13/29 word identities, weighted pairs, nested
  composition and direct CRT selector, including the L1000 endpoint hash.
- graft_audit/audit.py: separate permutation/split counts and exact coin boxes.
- graft_audit/audit_supplement.py: separate new-core identities, physical
  rectangles,377 CRT residues, all selected targets in the documented finite
  checks, nested-state checks and independent L1000 endpoint hash.
- endpoint/verify_endpoint.py and endpoint/response.py: binary selector and
  exact rooted response identity.
- endpoint_audit/audit.py: independent literal permutation/cut covers,
  recursive endpoints, actual stabilizers and arbitrary extensions.
- transducer/bridge_blocks.py and transducer/check_polynomial.py: transfer
  evaluation versus direct paths and direct odd-cycle packings.
- transducer_audit/audit.py: all three coefficient-level states through order8,
  direct paths through9, unpruned 262143 compositions through18, and original
  run-window adjacency controls.
- check_legacy_flip.py: all 120 single-arc flips of transitive tournaments at
  orders2..9, including the repaired adjacent-order argument and corollary.

The programs retain explicit exceptions under optimization. A failure-control
argument is not an ordinary mathematical input; see replay.py for each exact
option. Graft audits accept --source PATH if run outside the sibling layout.

## Scope

These checks validate literal small cores, finite word identities, exact
composition formulas and specified finite universes. The proof supplies all
repetition counts and the interval selector. No literal enumeration of paths
on the thousands-of-vertices constructed tournaments is claimed. The large
endpoints are exact integers; hashes encode them as unsigned big-endian bytes
without a leading zero byte.

Discovery generated small tournament/cut banks and explored finite carry
states. Generated banks, executables and unsuccessful scan logs are not proof
dependencies and are not included. Historical audit source hashes identify
pre-promotion files; final SHA256SUMS and final review identify the assembled
published versions, including portability and live-check hardening.
