# Final synthesis audit

PASS for the reviewed Session47 root synthesis and README, after explicit
baseline, grounding, n>=2 and ordinary-phase scope clarifications.

[AUDIT.md](AUDIT.md) records the findings and remaining boundaries.
check.py independently rechecks compact certificate arithmetic and finite
projections:376,511 controls, identical in normal and optimized Python.
check_manifests.py verifies226 files across14 frozen packets.

Use --root-dir to select a workspace or published bundle root; omitting it
uses the parent directory. See AUDIT.md for exact nonmutating commands.
REVIEWED_INPUTS.json pins the inspected documents. inputs/*.txt contains
verbatim archived root prose, not relocated live Markdown.

The full native classifications and all-state28-slot connectivity remain
results of their separately frozen independent audits. No new all-order
native wiring claim, Git operation, or full-repository CI run is asserted.

check_portability.py tests an actual manifest-only copied publication bundle
with isolated Python, not a symlink view of the source workspace.
