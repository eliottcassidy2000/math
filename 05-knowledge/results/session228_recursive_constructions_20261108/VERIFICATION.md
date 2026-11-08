# Final verification scope

## Reproducible exact runs

From this directory run:

    python3 -I -B verify_all.py > replay.local.json
    cmp REPLAY.json replay.local.json

The canonical aggregate runs 19 producer/independent-auditor entry points in both normal and optimized Python, for 38 executions. Every execution must exit successfully and reproduce its committed receipt byte for byte. The outer verifier and canonical correctness gates use explicit failures, so optimized Python does not remove them. The final aggregate passed after all packets were copied into their published locations.

Canonical scripts use the standard library and packet-local imports. The frozen exploratory coprime producer scripts under `coprime_audit/inputs/` are provenance only: their SymPy use, assertions and scratch imports are not dependencies and are not covered by the optimized-verification claim. The independent coprime verifier is the canonical check. The optional diagnostic `__main__` block in `mobile_triangle/transfer.py` also uses SymPy, but the canonical runs import its standard-library transfer functions and do not execute that diagnostic.

All-order claims rest on the supplied combinatorial proofs or exact matrix identities plus induction. A census through order six, checks through a fixed parameter, or 577 values through 4095 is not presented as proof beyond its finite scope.

## Package integrity and source checks

Run:

    python3 -I -B verify_package.py > package.local.json
    cmp PACKAGE_CHECK.json package.local.json

The package checker verifies every included SHA256SUMS manifest, parses every Python source file, and checks local Markdown destinations in authored documents. Frozen historical input documents retain their original relative references and are deliberately excluded from the Markdown-destination check. Their bytes are still hash-checked. Root SHA256SUMS covers every committed packet file except itself and PACKAGE_CHECK.json; the latter deterministically records the manifest hash and check totals.

The repository's `agents/check_docs.py` passed after the final navigation updates. Scoped `git diff --check` passed. Navigation links were checked to exist. No unrelated working-tree file was included in the commits.

## Audit status and portability

Producer proofs are immutable audit snapshots. A historical status such as “audit requested” or “pending independent audit” inside one of those snapshots is superseded by its linked independent PASS report and this final session status. The initial mobile proof, exact membership/strong extensions, cyclic inverse, recursive grammar, physical-state compression, coprime law, general divisibility, and unsaturated construction each have an independent report.

The cyclic-audit manifest had sibling source paths normalized for the packaged layout; the expected content hashes were unchanged. The recursive auditor was made packet-relative before its canonical replay. The final checks require no exploratory scratch directory.

## Publication versus CI

Each checkpoint was committed on main through a nonforce expected-head update, then checked by fetching origin/main and matching its entire tree to the local checkpoint tree. This verifies publication, not a CI run. GitHub reported no check runs on checkpoint 5 at inspection time. No GitHub Actions success is claimed; the mathematical and package checks above were run explicitly.
