# Publication review: two-unit corridor and first owner correction

**Independent scoped PASS.** All [producer proof files](two_unit_cells.md),
[first-owner proof](first_owner.md), primary receipts, manifests, and the
[independent audit](audit/AUDIT.md) are preserved byte-for-byte. Historical
"audit pending" wording in the producer files is retained as provenance; the
completed audit is authoritative for this exact frozen extension.

The all-order corridor is 0 <= delta <= 2k-1 for k >= 3. The first-owner
correction at delta = 2k is proved for k >= 4. The only nontrivial owner fibre
has one two-empty-gap matching presentation and two P3/one-unit-gap presentations
with the actual external port retained. The correction is K_k - 2H_k,
or -3(k-1)(2k-3)!! without the strong restriction.

The k=2 and k=3 first-owner layers are separately finite-exact. The k=4
C/D census is only those owner families, not the full mass layer. No result
for delta > 2k, unrestricted tournament recurrence, new palette-forcing theorem,
or revised approval of the unchanged base packet is included.

The incoming [deficit-k result](../zero_heavy_next_20261020h/README.md) is
credited rather than claimed as new. The independently reviewed
[base packet](../zero_heavy_robustness_20261021a/PUBLICATION_REVIEW.md) remains
a separate prerequisite. Its historical source pin differs from the corrected
current note. The exact original source is bundled as
[a historical capture](audit/frozen_base_zero_heavy.md), authenticated against
commit c4a7875652a5ea516e93a31c3c6410d2204773b2. Both lineages remain explicit;
no producer, audit, dependency pin, or mathematical statement was refreshed.

## Portable reproduction

From this directory, the complete producer verification is:

    python3 -B verify_addendum.py

Its output must match verify-addendum.out exactly. It authenticates the frozen
producer files and reruns all four producer checkers in normal and optimized
Python. The complete independent verification is:

    python3 -B audit/verify.py --source . --base ../zero_heavy_robustness_20261021a --repo ../../..

Its output must match audit/verification.json exactly. It authenticates both
frozen source packets, reads the exact current repository dependencies except
for its explicit historical base capture, and reruns the independent checker
in normal and optimized Python. Use these explicit paths with this publication
layout; the frozen default sibling-directory paths are deliberately unchanged.
Both verification entry points can themselves also be run under -O.

Preparation reconstructs the dependency inputs from one pinned Git snapshot
outside the repository before replay, so concurrent publication cannot change
the inputs mid-check. The original base packet is an explicit exact-byte prior
publication requirement; no substitute copy is added to this extension.

## Historical relative links

The unchanged historical capture contains two relative links that would point
inside this extension after relocation. Two explicitly labeled navigation-only
locators resolve them to the canonical [sparse theorem](../defect_boundary_20261020/sparse_theory.md)
and [boundary note](../defect_boundary_20261020/boundary.md). Both targets are
byte-identical at the historical source commit and the preparation snapshot.
No mathematical content is copied or modified by these locators.
