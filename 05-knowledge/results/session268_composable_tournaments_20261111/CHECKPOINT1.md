# Checkpoint1: interaction and exact response compression

The starting main commit was ff37a34648c16c726e07ac1e36fdd7d9be1da350.
This checkpoint publishes complete proofs and exact controls for the new
interacting-cell family and regular-tower parity dynamics, together with the
response-inventory implementation available at this checkpoint. The research
continues toward contextual tomography, arithmetic-family classification and
an integrated synthesis.

The old shared checkout is not edited. An isolated sparse worktree had30
line-ending-only inbox artifacts immediately upon checkout; these are not
part of the research changes and are never staged. Publication uses an
explicit owned-path list on the latest main tree, preserving concurrent work.

Verification is deliberately local to new results: explicit seam and class
counts, small labelled controls for the response inventory, and direct binary
matrix controls. No inherited Lean rebuild or broad historical audit is
claimed. Python checks use explicit exceptions, so optimized mode does not
silently remove them.
