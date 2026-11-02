# Verification and limits

## Checkpoint1

From this packet, run:

```sh
python3 -B replay_checkpoint1.py
sha256sum -c SHA256SUMS
```

`REPLAY_CHECKPOINT1.json` records a fresh relocation and execution from an
unrelated working directory. Both ordinary Python and `-O` reproduce the
frozen unit-distance, aliquot, and independent outputs byte-for-byte. Closed
aliquot/audit manifests are separately checked, including their own nested
normal/optimized replays. Root source-integrity checks also run in both modes.

Deliberate false claims must be rejected in both modes:57 changed to56 for
the lattice21 edge count; s(18)=21 changed to20;704-694=10 changed to20;
and a same-size proof-byte corruption. All failures use explicit exceptions,
not assertions disabled by optimization. Restoration is checked afterwards.
These are fail-closed controls, not a claim that every possible bug is caught.

The source manifest contains exact frozen selections. The final source report
checks original Git blob and SHA256 identities, snapshot bytes and current
publication links. The packet-wide manifest seals every published file except
itself. Historical source snapshots preserve their original link context.

Documentation checks are deliberately scoped: startup byte/line limits,
owned-file conflict markers, current publication links and unchanged unrelated
frontier boundaries. The full documentation checker was **not run**, because
it reads the opaque central mistakes ledger, which was neither read nor
hydrated. This is not a full-repository documentation pass.

## Mathematical scope

The all-size statements depend on the written proofs and independent review;
finite enumeration is corroboration. The unit-distance producer checks all
1100 labelled simple graphs of orders0..5 and1094 edge subdivisions, exact
algebraic products through four factors, and all210 pairs of the lattice21
witness. The independent graph checker instead enumerates all orientations
and all spanning edge subsets at orders1..5. Aliquot producer sieves and
multiplicative products agree through200000; independent direct-divisor
checks run through10000, with proven finite bounds for inverse fibres.

AMP's dated extremality theorem remains CITED, not independently re-enumerated.
No Lean build, exact HN value, Collatz convergence, or global aliquot
classification is certified. New source selections and proof refinements
need a later checkpoint and their own audit; no frozen first-checkpoint
scientific file may be silently rewritten.
