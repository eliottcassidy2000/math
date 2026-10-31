# Verification and publication scope

`python3 -B verify_session.py` replays the six mathematical checkers in fresh
relocated copies, in both normal and optimized Python. All twelve executions
passed, expected outputs matched byte for byte, and six component integrity
manifests verified40 entries. The exact receipt is `SESSION_REPLAY.json`.
Python3 is required; the volume checker also requires SymPy, whose version
is recorded in its own output. No network or source checkout is needed.

The independent peer review covers the final all-order proof, the GSMA gates,
and the synthesis. The geometry review found no substantive gap; its three
wording/endpoint clarifications were applied before publication. Root review
retains the explicit sqrt(2)-weighted Minkowski metric in the >180 claim.

Local Markdown routes were checked across the session navigation and every
packet document. The final route receipt is retained with the publication
journal; no external-link availability audit is implied by this local check.
The owned whitespace check passes.

The repository documentation checker was run with its bounded sparse fixture:
34 theorem inputs (33 frontier files and one router-smoke candidate), not a
full theorem-graph audit. It reports one unchanged inherited failure:
`05-knowledge/hypotheses/INDEX.md` has121 maintained lines against budget120.
No unrelated index or ledger edit was made to suppress that failure.

The exact computational certificates check finite domains and algebraic
identities. They do not computationally reprove the cited Delone–Nagell,
Thue–Mahler, Legendre, or Minkowski theorems. Quartic/quintic complete solution
lists, all-height sparse-power exclusion, and a local Lean build are not
claimed. No CI success is asserted.
