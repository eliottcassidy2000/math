# Fixed-defect endpoint spectra at a sole exit

**Status: PROVED elementary statements; independent proof and computational audit PASS. Chan–Pak remains OPEN. No literature-priority claim.**

This packet strengthens the *sole-exit part* of Session213's endpoint-equality obstruction from exact excess one to every fixed positive odd excess.

For a sole-exit mark a->r, put U=T-{a,r} and O=N+_U(r). Let K count ordered spanning two-path covers of U whose first endpoint and second start both lie in O. The exact identity is

    H(T)-2e_a(T)=H(U)+K.

For strong T at every fixed defect d, deleting suitable singleton strong components of U leaves a bounded physical skeleton of order at most

    d+floor(log2(max(1,d-1)))+4.

The removed singletons provide the only unbounded scalar coordinate. Consequently the endpoint and total-path spectra at fixed d are finite unions of affine dyadic rays, with only O_d(log X) attained values up to X. SCC factorization extends the finite affine-dyadic description to arbitrary T. This is an all-order theorem and an effective, potentially expensive, finite inverse.

At defect three, the total-path values are exactly

    {5} union {6*2^q+3:q>=0}
        union {6*2^q-1:q>=1} union {10*2^q-1:q>=1}.

The explicit constructions include 23,47,95,... and 19,39,79,...; `WITNESSES.json` supplies six checked adjacency-row examples. These families are not asserted to be previously unknown. The new content claimed here is the stated classification and bounded-skeleton proof, not literature priority.

## Scope

The marked endpoint in the **output** must itself have one exit. Deleting a sole-exit vertex from a cycle tournament does not automatically give that property at the new marked endpoint. This packet does not exclude arbitrary endpoint defects or general cycle-to-path conversions. It does not increase the inherited finite coverage bank or prove an eventual odd tail.

## Files

- `PROOF.md`: identities, bounded skeleton, exact inverse, nonstrong reduction, bounded-defect sparsity and complete defect-three classification.
- `verify.py`: self-contained finite checks with explicit exceptions, no repository imports.
- `RESULT.json`, `RESULT.optimized.json`: identical default runs.
- `RESULT.classes.json`: default tests plus supplied order-seven/eight class files.
- `WITNESSES.json`: concrete total-path realizers on six to eight vertices.
- `SOURCES.json`: source URLs, commit and verified blob identities at main `4c63fd58b566f3e006eb87eaedcf48824741cc84`.

## Reproduction

From this directory:

    python3 verify.py > /tmp/endpoint-defect-normal.json
    cmp RESULT.json /tmp/endpoint-defect-normal.json
    python3 -O verify.py > /tmp/endpoint-defect-optimized.json
    cmp RESULT.optimized.json /tmp/endpoint-defect-optimized.json

The default run checks all 33,867 labeled tournaments of orders one through six: 32,110 sole-exit marks and 22,654 strong compression instances. It also checks 1,098 literal path-cover configurations through order four, 165 transitive pair weights, 42 explicit ray constructions, 210 deterministic random attachments, and an actual unmarked H=13 hostile. Counts for these separate universes are not added as if all objects were distinct.

From the repository root after placement in `05-knowledge/results/session214_chan_pak_combined_proof_20261107/path-techniques`:

    python3 05-knowledge/results/session214_chan_pak_combined_proof_20261107/path-techniques/verify.py --class-files session212/cycle-levels/classes7.rows session212/cycle-levels/classes8.rows > /tmp/endpoint-defect-classes.json
    cmp 05-knowledge/results/session214_chan_pak_combined_proof_20261107/path-techniques/RESULT.classes.json /tmp/endpoint-defect-classes.json

The additional files contain 456 and 6,880 supplied isomorphism-class representatives. Their byte hashes are in the receipt. The verifier validates every tournament and theorem instance but does not independently prove the completeness or isomorphism-uniqueness of those inherited class files. The independent audit provides a separate literal path and path-cover counting route. It checked 32,108 nonempty-U sole-exit identities, 29,998 positive endpoint roots, 22,654 strong compression cases, 4,824 nontrivial compressions, and 6,024 defect-three roots through order six. It approved the frozen proof with SHA-256 `8dc84768893ecb8f0fcfcce4ed3ef89e16e98a37d0b38b6f22a4e5ea5b46c420`; this final audit supersedes the proof draft's pending-audit label.

## Inherited barriers and current-main check

Session213's two proofs that the binary cycle-comb has only O(X^(log_3 2)) total-path outputs count positive three-way allocations. Padding is included; the finite-menu induced-extension barrier follows from induced monotonicity. Its product closure uses a convergent Euler-product majorant for every exponent strictly above log_3 2. These arguments do not exclude arbitrary endpoint constructions, and the present theorem does not use a comb premise. The endpoint equality, capacity and atlas sources byte-match their later versions at the pinned main snapshot; the additional THM-4900 gate-tower asymptotics concern a different observable and do not supply an exact spectrum inverse. Source and inherited class-file identities are recorded in SOURCES.json.

The next unresolved step is outside the theorem: understand bounded endpoint defects when the output mark has multiple exits, or construct exact targets with an unbounded defect. No exploratory assertion in that direction is promoted here.
