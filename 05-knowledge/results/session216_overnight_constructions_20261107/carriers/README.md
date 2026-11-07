# Common right factors and universal pair restrictions

**PROVED, independently audited in the scopes below. Chan–Pak completeness remains OPEN.** No eventual spectrum, unrestricted coverage extension, or order-ten universal-pair exclusion is included in this release.

## Baseline common-factor and marked-module results

`baseline/PROOF.md` classifies a common rooted right factor on an already synchronized two-residue pair: adjusted-base synchronization survives exactly when the right deletion count is1. These factors are physically O/Z padding. This does not exclude asymmetric factors or source replacement.

The complete single-nonroot-module probe on six inherited pairs has four successful position pairs, all equivalent to existing O/Z families. The baseline finite asymmetric searches reach3,467,984 candidates after prefixes through length2 and find no persistent match involving a nontrivial right factor. These are supplied-bank, finite-menu statements only.

## New all-order universal-pair necessary conditions

`reciprocity/UNIVERSAL_PAIR_RECIPROCITY.md` proves rooted integer-polynomial reciprocity and applies it to the inherited all-O/Z/J synchronization criterion. Any universal width-eight pair must have opposite order parity, even starts in one of(4,2),(10,8),(16,14),(22,20), all eight common nonstart coordinates divisible by3, and the sharper cover congruences stated there.

These are necessary conditions, not proof of nonexistence of every pair. The separate809-source, depth-seven search contains2,653,520 presentations and finds no universal pair inside its declared menu. It does not establish any all-order or order-ten absence claim. Historical same-graph searches are finite controls now superseded, for their universal-pair goal, by the parity theorem.

## Reproduction and audit

Use Python3 and a C++17 compiler. The producer reciprocity polynomial verifier additionally uses SymPy; both independent Python auditors use only the standard library. Run commands inside the named directory because some baseline scripts recreate receipt files.

From `baseline/`:

    python3 probe.py
    python3 verify.py
    python3 -O verify.py
    g++ -O2 -std=c++17 amalgamation_prefix_reset.cpp -o /tmp/carrier-prefix-search
    /tmp/carrier-prefix-search BANK.txt

The two other amalgamation searches use the same build and BANK.txt invocation.

From `reciprocity/`:

    python3 verify_reciprocity.py
    python3 -O verify_reciprocity.py
    g++ -O2 -std=c++17 root_modular_candidates.cpp -o /tmp/carrier-root-search
    /tmp/carrier-root-search BANK.txt
    g++ -O2 -std=c++17 modular_deep_search.cpp -o /tmp/carrier-modular-search
    /tmp/carrier-modular-search ROOT_MODULAR_CANDIDATES.txt 7

From each `audit/carriers/` and `audit/reciprocity/` directory, run `python3 independent.py`. Additional literal and finite-search checks are documented in their REVIEW.md files. These audit folders preserve original reviewed inputs and receipts without status edits. Incidental macro-search/ALL_PAIRS files occur only in the baseline audit freeze to preserve its original exact manifest; they are explicitly outside the approved claims and absent from the baseline release payload. They are retained as historical audit inputs, not promoted results.

PROMOTION.diff records the three status-only changes. SHA256SUMS pins the complete release payload; PUBLICATION_FILES.txt lists files excluding those two manifest/list files. No binaries or caches are included. No publishing was performed by the release preparer.
